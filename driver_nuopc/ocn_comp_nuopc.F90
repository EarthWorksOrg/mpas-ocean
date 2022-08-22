module ocn_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for MPAS-Ocean
  !----------------------------------------------------------------------------
  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompFilterPhaseMap, NUOPC_IsUpdated, NUOPC_IsAtTime
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_CompSetClock
  use NUOPC                 , only : NUOPC_SetAttribute, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model           , only : model_routine_SS           => SetServices
  use NUOPC_Model           , only : SetVM
  use NUOPC_Model           , only : model_label_Advance        => label_Advance
  use NUOPC_Model           , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model           , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model           , only : model_label_CheckImport    => label_CheckImport
  use NUOPC_Model           , only : model_label_SetClock       => label_SetClock
  use NUOPC_Model           , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  !use perf_mod              , only : t_startf, t_stopf
  use ocn_import_export     , only : ocn_advertise_fields, ocn_realize_fields
  use ocn_import_export     , only : ocn_import, ocn_export, ocn_cpl_dt !, tlast_coupled
  use nuopc_shr_methods     , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use nuopc_shr_methods     , only : set_component_logging, get_component_instance, log_clock_advance
  use shr_kind_mod          , only : cl=>shr_kind_cl, cs=>shr_kind_cs, SHR_KIND_CX
  use shr_file_mod 
  use shr_sys_mod 
  use mpas_framework
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_stream_manager
  use mpas_abort
  use ocn_core_interface
  use shr_pio_mod
  use pio
  use perf_mod
  use ocn_config
  use ocn_gm
  use ocn_diagnostics
  !use mpas_ocn_constants, only : coupleAlarmID
  use ocn_tracer_ecosys
  use ocn_tracer_CFC
  use ocn_tracer_surface_restoring
  use ocn_tracer_short_wave_absorption_variable
  use ocn_analysis_driver
  use ocn_time_integration
  use ocn_frazil_forcing 
  use ocn_surface_land_ice_fluxes
  use ocn_forcing
  use ocn_time_average_coupled
  
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private                              ! By default make data private

  public  :: SetServices
  public  :: SetVM
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  integer, parameter  :: dbug = 1
  character(*), parameter :: u_FILE_u = &
       __FILE__

  integer           :: lmpicom
  integer           :: stdout 
  integer           :: nThreads        ! number of threads per mpi task for this component
  character(len=CL)   :: flds_scalar_name = ''
  integer             :: flds_scalar_num = 0
  integer             :: flds_scalar_index_nx = 0
  integer             :: flds_scalar_index_ny = 0

  type (core_type), pointer :: corelist => null()
  type (dm_info), pointer :: dminfo
  type (domain_type), pointer :: domain_ptr
  type (iosystem_desc_t), pointer :: io_system 
  integer :: itimestep   ! time step number for MPAS

  integer :: ocnLogUnit ! unit number for ocn log
  character (len=*), parameter :: coupleAlarmID = 'coupling'
   character(len=StrKIND) :: coupleTimeStamp

! !PRIVATE MODULE VARIABLES

  integer, private ::   &
      my_task
   character(len=StrKIND) :: runtype

!=======================================================================
contains
!=======================================================================

  subroutine SetServices(gcomp, rc)

    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    character(len=*),parameter  :: subname='ocn_comp_nuopc:(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
         specRoutine=DataInitialize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_CheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
         specRoutine=ModelCheckImport, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !--------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    use NUOPC, only : NUOPC_isConnected

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)     :: vm
    integer           :: iam
    integer           :: shrlogunit
    character(len=CL) :: logmsg
    character(len=CS) :: cvalue
    logical           :: isPresent, isSet
    character(len=*), parameter :: subname='ocn_comp_nuopc:(InitializeAdvertise) '
    !--------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! reset shr logging to my log file
    call set_component_logging(gcomp, iam==0, stdout, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    flds_scalar_name = trim(cvalue)
    call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue, *) flds_scalar_num
    write(logmsg,*) flds_scalar_num
    call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_nx
    write(logmsg,*) flds_scalar_index_nx
    call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_scalar_index_ny
    write(logmsg,*) flds_scalar_index_ny
    call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Advertise fields
    call ocn_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    !-----------------------------------------------------------------------
    !  first initializaiton phase of pop2
    !  initialize the timers, communication routines, global reductions,
    !  domain decomposition, grid, and overflows
    !-----------------------------------------------------------------------
    use ESMF               , only: ESMF_VMGet
      use mpas_stream_manager, only : MPAS_stream_mgr_init, MPAS_build_stream_filename, MPAS_stream_mgr_validate_streams
      use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
      use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string
      use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_time
      use mpas_bootstrapping, only : mpas_bootstrap_framework_phase1, mpas_bootstrap_framework_phase2
      use mpas_log


    ! Initialize MPASO

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    !  local variables
    integer                 :: iblock
    type(ESMF_VM)           :: vm
    type(ESMF_DistGrid)     :: distGrid
    type(ESMF_Mesh)         :: Emesh
    integer , allocatable   :: gindex_ocn(:)
    integer , allocatable   :: gindex_elim(:)
    integer , allocatable   :: gindex(:)
    integer                 :: globalID
    character(CL)           :: cvalue
    integer                 :: num_elim_global
    integer                 :: num_elim_local
    integer                 :: num_elim
    integer                 :: num_ocn
    integer                 :: num_elim_gcells ! local number of eliminated gridcells
    integer                 :: num_elim_blocks ! local number of eliminated blocks
    integer                 :: num_total_blocks
    integer                 :: my_elim_start
    integer                 :: my_elim_end
    integer                 :: lsize
    integer                 :: shrlogunit      ! old values
    integer                 :: npes
    integer                 :: iam, ierr_local
    logical                 :: mastertask
    character(len=32)       :: starttype
    integer                 :: n,i,j,iblk,jblk,ig,jg,ierr
    integer                 :: lbnum
    integer                 :: ocnid
    integer                 :: errorCode       ! error code
    type(ESMF_Time)         :: EcurrTime
    type(ESMF_Calendar)     :: Ecalendar
    logical                 :: associated_e
    integer(I8Kind)         :: s_e, sn_e, sd_e
    integer                 :: yy_e
    type(ESMF_CalKind_Flag)   :: type_e
    type(ESMF_TimeInterval)   :: timeStep        ! Model timestep
    type (MPAS_TimeInterval_type) :: denInterval, remInterval, zeroInterval
    type (MPAS_TimeInterval_Type) :: alarmTimeStep
    type (MPAS_Time_Type) :: alarmStartTime
    integer (kind=I8KIND) :: numDivs
    type (MPAS_Time_Type) :: currTime
    type (MPAS_timeInterval_type) :: mpastimeStep

    character(len=*), parameter  :: subname = "ocn_comp_nuopc:(InitializeRealize)"
    integer, dimension(:), pointer :: indexToCellID, indextocellid_0halo
    logical :: readNamelistArg, readStreamsArg
    character(len=SHR_KIND_CX) :: argument, namelistFile, streamsFile
    type (block_type), pointer :: block
    integer :: iCell, nCells, ncells_0halo
    integer, dimension(:), pointer :: nCellsArray
    type (mpas_pool_type), pointer :: meshPool, statePool, &
                                      forcingPool, &
                                      averagePool, scratchPool

      character(len=StrKIND) :: iotype
      logical :: streamsExists
      integer :: mesh_iotype
      type (c_ptr) :: mgr_p

      character(len=StrKIND) :: mesh_stream
      character(len=StrKIND) :: mesh_filename
      character(len=StrKIND) :: mesh_filename_temp
      character(len=StrKIND) :: ref_time_temp
      character(len=StrKIND) :: filename_interval_temp
      character(len=StrKIND) :: caseid
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
      character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_ref_time_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_filename_interval_temp
      character(kind=c_char), dimension(StrKIND+1) :: c_iotype

      integer :: blockID

      character(len=16) :: inst_suffix
      integer :: inst_index

      character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
      integer(kind=c_int) :: c_comm
      integer(kind=c_int) :: c_ierr

      type (MPAS_Time_type) :: start_time
      type (MPAS_Time_type) :: ref_time
      type (MPAS_TimeInterval_type) :: filename_interval
      character(len=StrKIND) :: timeStamp
      integer :: ocnPossibleErrUnit !< unit number to reserve for if a err log needs to be opened
      logical :: exists
      integer :: pio_iotype
      integer :: shrloglev !< shr log level and log unit
      logical, pointer :: tempLogicalConfig
      character(len=StrKIND), pointer :: tempCharConfig
      real (kind=RKIND) :: dt

     ! Added for coupling interval initialization
      integer, pointer :: index_avgZonalSSHGradient, index_avgMeridionalSSHGradient
      real (kind=RKIND), dimension(:),   pointer :: filteredSSHGradientZonal, filteredSSHGradientMeridional
      real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient

      interface
         subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
            use iso_c_binding, only : c_char, c_ptr, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            type (c_ptr), intent(inout) :: mgr_p
            integer(kind=c_int), intent(inout) :: comm
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_parser

         subroutine xml_stream_get_attributes(xmlname, streamname, comm, filename, ref_time, filename_interval, io_type, ierr) bind(c)
            use iso_c_binding, only : c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: xmlname
            character(kind=c_char), dimension(*), intent(in) :: streamname
            integer(kind=c_int), intent(inout) :: comm
            character(kind=c_char), dimension(*), intent(out) :: filename
            character(kind=c_char), dimension(*), intent(out) :: ref_time
            character(kind=c_char), dimension(*), intent(out) :: filename_interval
            character(kind=c_char), dimension(*), intent(out) :: io_type
            integer(kind=c_int), intent(out) :: ierr
         end subroutine xml_stream_get_attributes
      end interface

    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    errorCode = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)


    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, PetCount=npes, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, pet=iam, peCount=nthreads, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(nthreads==1) then
       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) nthreads
    endif

!$  call omp_set_num_threads(nThreads)

#if (defined _MEMTRACE)
    if (iam == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','InitializeRealize:start::',lbnum)
    endif
#endif

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) caseid

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    !-----------------------------------------------------------------------
    !  first initializaiton phase of mpas-ocean
    !  initialize the timers, communication routines, global reductions,
    !  domain decomposition, grid, and overflows
    !-----------------------------------------------------------------------


      readNamelistArg = .false.
      readStreamsArg = .false.
      !?!? do I need what follow in mpas_subdriver/mpas_init?
!---the rest of this routine was copied from mpas_subdriver.F, subroutine mpas_init
    allocate(corelist)
    nullify(corelist % next)

    allocate(corelist % domainlist)
    nullify(corelist % domainlist % next)

    domain_ptr => corelist % domainlist
    domain_ptr % core => corelist

    call mpas_allocate_domain(domain_ptr)


!-----------------------------------------------------------------------
!
!   first initializaiton phase of mpaso
!   call mpaso initialization routines
!
!-----------------------------------------------------------------------

    call t_startf('mpaso_init1')

    ! ----------------
    ! Set up log file information
    ! ----------------
    !?!?inst_suffix = seq_comm_suffix(ocnID) ! a suffix to append to log file name
    ocnLogUnit = shr_file_getUnit() ! reserve unit number for log unit
    ocnPossibleErrUnit = shr_file_getUnit() ! reserve unit number for possible error log file

!    ! Note: In following code, a file ocn_modelio.nml is queried to determine the log file name to use.
!    ! This file is generated when a case is created.  The default file name in there is currently ocn.log.DATE-TIME.
!    if (iam==0) then
!       inquire(file='ocn_modelio.nml'//trim(inst_suffix),exist=exists)
!       call shr_file_setio('ocn_modelio.nml'//trim(inst_suffix), ocnLogUnit)
!    endif

    ! Store shr log unit and level so we can reassign them later when OCN control is complete
    call shr_file_getLogUnit(shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    ! Set the shr log unit to be the ocnLogUnit
    call shr_file_setLogUnit(ocnLogUnit)

    ! Write to ocnLogUnit here, because the log module is not open yet.
    if (iam==0) write(ocnLogUnit,'(a,i6)') '=== Beginning InitializeRealize for ocn: rank=',iam

    call mpas_framework_init_phase1(domain_ptr % dminfo, lmpicom)
    
      call ocn_setup_core(corelist)
      call ocn_setup_domain(domain_ptr)

    ! ===========
    ! Initialize log manager
    call mpas_log_init(domain_ptr % logInfo, domain_ptr, unitNumbers=(/ocnLogUnit, ocnPossibleErrUnit/), err=ierr)
    if ( ierr /= 0 ) then
       write(ocnLogUnit,*) 'ERROR: log init failed for core ' // trim(domain_ptr % core % coreName)
       call mpas_dmpar_abort(domain_ptr % dminfo)
    end if

    ! Set core specific options here
    ! Disable output from all but the master task for E3SM!
    ! (This overrides the default set by mpas_log_init based on MPAS_DEBUG setting.)
    if (iam /= 0) then
       domain_ptr % logInfo % outputLog % isActive = .false.
    endif

    ! After core has had a chance to modify log defaults, open the output log
    call mpas_log_open(err=ierr)
    if ( ierr /= 0 ) then
       write(ocnLogUnit,*) 'ERROR: log open failed for core ' // trim(domain_ptr % core % coreName)
       call mpas_dmpar_abort(domain_ptr % dminfo)
    end if
    ! ===========


    ! ----------
    ! Process namelist and streams files
    ! ----------
    ! Override the names of the stream and namelist files
    domain_ptr % namelist_filename = 'mpaso_in'
    domain_ptr % streams_filename = 'streams.ocean'

    ! Setup namelist variables, and read the namelist
    ierr = domain_ptr % core % setup_namelist(domain_ptr % configs, domain_ptr % namelist_filename, domain_ptr % dminfo)
    if ( ierr /= 0 ) then
       call mpas_log_write('Namelist setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_framework_init_phase2(domain_ptr) !, io_system)


    ! Define package variables
    ierr = domain_ptr % core % define_packages(domain_ptr % packages)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package definition failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup packages (i.e. determine if they should be on or off)
    ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages, domain_ptr % ioContext)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup decompositions available for dimensions
    ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
    if ( ierr /= 0 ) then
       call mpas_log_write('Decomposition setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! ----------
    ! Override namelist options based on start type
    ! ----------
    if (trim(starttype) == trim('startup')) then
       runtype = "initial"
    else if (trim(starttype) == trim('continue') ) then
       runtype = "continue"
    else if (trim(starttype) == trim('branch')) then
       runtype = "continue"
    else
       !write(stdout,*) 'ocn_comp_nuopc: ERROR: unknown starttype'
       call ESMF_LogWrite(trim(subname)//'ERROR: unknown starttype',ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    end if

    if (runtype == "initial") then ! Start up run
        ! Turn off restart
        call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .false.

        ! Setup start time. Will be over written later when clocks are synchronized
        call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
        tempCharConfig = trim(tempCharConfig) // "0:00:00"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"
    else if (runtype == "continue" .or. runtype == "branch") then ! Restart run or branch run
        ! Turn on restart
        call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .true.

        ! Set start time to be read from file
        call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
        tempCharConfig = "file"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"
    end if

    ! Setup MPASO simulation clock
    ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
    if ( ierr /= 0 ) then
       call mpas_log_write('Clock setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_log_write('Reading streams configuration from file '//trim(domain_ptr % streams_filename))
    inquire(file=trim(domain_ptr % streams_filename), exist=streamsExists)

    if ( .not. streamsExists ) then
       call mpas_log_write('Streams file '//trim(domain_ptr % streams_filename)//' does not exist.', MPAS_LOG_CRIT)
    end if

    !
    ! Using information from the namelist, a graph.info file, and a file containing
    !    mesh fields, build halos and allocate blocks in the domain
    !
    ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
    if ( ierr /= 0 ) then
       call mpas_log_write('Failed to find mesh stream for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if


    call mpas_f_to_c_string(domain_ptr % streams_filename, c_filename)
    call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
    c_comm = domain_ptr % dminfo % comm
    call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                   c_mesh_filename_temp, c_ref_time_temp, &
                                   c_filename_interval_temp, c_iotype, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_get_attributes failed.', MPAS_LOG_CRIT)
    end if
    call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
    call mpas_c_to_f_string(c_ref_time_temp, ref_time_temp)
    call mpas_c_to_f_string(c_filename_interval_temp, filename_interval_temp)
    call mpas_c_to_f_string(c_iotype, iotype)

    if (trim(iotype) == 'pnetcdf') then
       mesh_iotype = MPAS_IO_PNETCDF
    else if (trim(iotype) == 'pnetcdf,cdf5') then
       mesh_iotype = MPAS_IO_PNETCDF5
    else if (trim(iotype) == 'netcdf') then
       mesh_iotype = MPAS_IO_NETCDF
    else if (trim(iotype) == 'netcdf4') then
       mesh_iotype = MPAS_IO_NETCDF4
    else
       mesh_iotype = MPAS_IO_PNETCDF
    end if

    start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
    if ( trim(ref_time_temp) == 'initial_time' ) then
        call mpas_get_time(start_time, dateTimeString=ref_time_temp, ierr=ierr)
    end if

    if ( trim(filename_interval_temp) == 'none' ) then
        call mpas_expand_string(ref_time_temp, -1, mesh_filename_temp, mesh_filename)
    else
        call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
        call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
        call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, -1, mesh_filename, ierr)
    end if

    ! Bootstrap framework (1). Here data structures are setup, but dimensions and arrays are not finalized.
    call mpas_log_write(' ** Attempting to bootstrap MPAS framework using stream: ' // trim(mesh_stream))
    call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

    !
    ! Set up run-time streams
    !
    call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % ioContext, domain_ptr % clock, &
                              domain_ptr % blocklist % allFields, domain_ptr % packages, domain_ptr % blocklist % allStructs)

    call add_stream_attributes(domain_ptr)

    ! Setup all immutable streams for the core
    ierr = domain_ptr % core % setup_immutable_streams(domain_ptr % streamManager)
    if ( ierr /= 0 ) then
       call mpas_log_write('Immutable streams setup failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Parse / read all streams configuration
    mgr_p = c_loc(domain_ptr % streamManager)
    call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_parser failed.', MPAS_LOG_CRIT)
    end if

    my_task = domain_ptr % dminfo % my_proc_id

    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain_ptr)
    call t_stopf('mpaso_init1')

    call t_startf('mpaso_init2')
    ! Initialize the MPASO core
    ierr = domain_ptr % core % core_init(domain_ptr, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_log_write('Core init failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    end if
    call t_stopf('mpaso_init2')

!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!   check for consistency of mpaso and sync clock initial time
!
!-----------------------------------------------------------------------

    call ESMF_ClockGet( clock, currTime=EcurrTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    !DD there hase to be a better way to go from esmf type to MPAS_Time_Type
    call ESMF_TimeGet(Ecurrtime, s_i8=s_e, sn_i8=sn_e, sd_i8=sd_e, yy=yy_e, calendar = ecalendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    currtime%t%basetime%S  = s_e
!    currtime%t%basetime%Sn = sn_e
!    currtime%t%basetime%Sd = sd_e
!    currtime%t%yr = yy_e
    currtime%t = Ecurrtime

    call ESMF_CalendarGet(Ecalendar, calkindflag=type_e, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!    if(type_e == ESMF_CALKIND_NOLEAP) then
!       currTime%t%calendar => noleapCal
!    elseif(type_e == ESMF_CALKIND_GREGORIAN) then
!       currTime%t%calendar => gregorianCal
!    else
!       rc = 1
!       if (ChkErr(rc,__LINE__,u_FILE_u)) return
!    endif

    if (runtype == 'initial') then
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    else if (runtype == 'continue' .or. runtype == 'branch') then
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    end if
    ! Doublecheck that clocks are synced here.  
    ! (Note that if this is an initial run, a section below will advance the ocean clock
    ! by one coupling interval, at which point we expect the clocks to be OUT of sync.)
    !DD add the check later
    !?!?call check_clocks_sync(domain % clock, Eclock, ierr)

    !-----------------------------------------------------------------------
    ! initialize necessary coupling info
    !-----------------------------------------------------------------------

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeIntervalGet( timeStep, s=ocn_cpl_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    call mpas_set_timeInterval(alarmTimeStep, S=ocn_cpl_dt, ierr=ierr)
    call mpas_get_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)

    ! Verify the mpas time step fits into a coupling interval
    call mpas_pool_get_config(domain_ptr % configs, 'config_dt', tempCharConfig)
    call mpas_set_timeInterval(denInterval, timeString=tempCharConfig, ierr=ierr)
    call mpas_set_timeInterval(zeroInterval, S=0, ierr=ierr)
    call mpas_interval_division(start_time, alarmTimeStep, denInterval, numDivs, remInterval)

    ierr = 0

    if ( alarmTimeStep < denInterval ) then
       ierr = 1
    end if
    ierr_local = ierr
    call mpas_dmpar_max_int(domain_ptr % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        OCN Model time step is: ' // trim(tempCharConfig), MPAS_LOG_ERR)
       call mpas_log_write('        The model time step cannot be longer then the coupling interval', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    if ( remInterval > zeroInterval ) then
       ierr = 1
    end if

    ierr_local = ierr
    call mpas_dmpar_max_int(domain_ptr % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        OCN Model time step is: ' // trim(tempCharConfig), MPAS_LOG_ERR)
       call mpas_log_write('        These are not synchronized, so time steps will not match to coupling interval boundaries.', MPAS_LOG_ERR)
       call mpas_log_write('        Please reconfigure either the coupling interval or the time step.', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    ! set coupling alarm
    alarmStartTime = currTime
    call mpas_add_clock_alarm(domain_ptr % clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
    call mpas_print_alarm(domain_ptr % clock, coupleAlarmID, ierr)
    call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

    mpastimeStep = mpas_get_clock_timestep(domain_ptr % clock, ierr=ierr)
    call mpas_get_timeInterval(mpastimeStep, dt=dt)

    ! Build forcing arrays.
    block => domain_ptr % blocklist
    do while(associated(block))
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)
        call mpas_timer_start("land_ice_build_arrays", .false.)
        call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                      forcingPool, scratchPool, statePool, dt, ierr)
        call mpas_timer_stop("land_ice_build_arrays")

        call ocn_frazil_forcing_build_arrays(domain_ptr, meshPool, forcingPool, statePool, ierr)

        block => block % next
    end do

    currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
    call mpas_log_write( 'Initial time '//trim(timeStamp))

    ! read initial data required for variable shortwave
    call mpas_timer_start('io_shortwave',.false.)
    call ocn_get_shortWaveData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .true.)
    call mpas_timer_stop('io_shortwave')

    ! read initial data required for ecosys forcing
    call mpas_pool_get_config(domain_ptr % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
    if (config_use_ecosysTracers) then
        call mpas_timer_start('io_ecosys',.false.)
        call ocn_get_ecosysData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .true.)
        call mpas_timer_stop('io_ecosys')
    endif

    ! read initial data required for CFC forcing
    call mpas_pool_get_config(domain_ptr % configs, 'config_use_CFCTracers', config_use_CFCTracers)
    if (config_use_CFCTracers) then
        call mpas_timer_start('io_CFC',.false.)
        call ocn_get_CFCData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .true.)
        call mpas_timer_stop('io_CFC')
    endif

    ! read initial data required for monthly surface salinity restoring
    call mpas_pool_get_config(domain_ptr % configs, 'config_use_activeTracers_surface_restoring',  &
                                                 config_use_activeTracers_surface_restoring)
    call mpas_pool_get_config(domain_ptr % configs, 'config_use_surface_salinity_monthly_restoring',  &
                                                 config_use_surface_salinity_monthly_restoring)
    if (config_use_activeTracers_surface_restoring .and.  config_use_surface_salinity_monthly_restoring) then
        call mpas_timer_start('io_monthly_surface_salinity',.false.)
        call ocn_get_surfaceSalinityData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .true.)
        call mpas_timer_stop('io_monthly_surface_salinity')
    endif

    ! Setup clock for initial runs
    if ( runtype == 'initial' ) then
       ! Advance clock one coupling interval to be in sync with the coupler clock.
       !do while (.not. mpas_is_alarm_ringing(domain_ptr % clock, coupleAlarmID, ierr=ierr))
       !   itimestep = itimestep + 1
       !   call mpas_advance_clock(domain_ptr % clock)
       !end do

       block => domain_ptr % blocklist
       do while(associated(block))
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(statePool, forcingPool, 1)
          block => block % next
       end do
    end if

!-----------------------------------------------------------------------
!
!   initialize coupling variables
!   NOTE: could be moved to subroutine
!
!-----------------------------------------------------------------------

    if ( runtype == "initial") then
       block => domain_ptr % blocklist
       do while(associated(block))
          call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
          call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
          call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)
          call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
          call mpas_pool_get_array(forcingPool, 'filteredSSHGradientZonal', filteredSSHGradientZonal)
          call mpas_pool_get_array(forcingPool, 'filteredSSHGradientMeridional', filteredSSHGradientMeridional)
          filteredSSHGradientZonal      = avgSSHGradient(index_avgZonalSSHGradient,     :)
          filteredSSHGradientMeridional = avgSSHGradient(index_avgMeridionalSSHGradient, :)
          block => block % next
       end do
    endif

    !---------------------------------------------------------------------------
    ! Determine the global index space needed for the distgrid
    !---------------------------------------------------------------------------

    n = 0
    block => domain_ptr % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )
       n = n + ncells
       block => block % next
    end do
    lsize = n
    allocate(gindex_ocn(lsize))

    n = 0
    block => domain_ptr % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       nCells = nCellsArray( 1 )
       do iCell = 1, nCells
          gindex_ocn(n+iCell) =indexToCellID(iCell)
       enddo
       n = n + ncells
       block => block % next
    end do

!       ! No eliminated land blocks
       num_ocn = size(gindex_ocn)
       allocate(gindex(num_ocn))
       do n = 1,num_ocn
          gindex(n) = gindex_ocn(n)
       end do

    !---------------------------------------------------------------------------
    ! Create distGrid from global index array
    !---------------------------------------------------------------------------
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the MPAS-O mesh
    !---------------------------------------------------------------------------

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, &
         elementDistgrid=Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !mastertask = iam == domain_ptr % dminfo % my_proc_id
    mastertask = iam == 0
    if (mastertask) then
       write(stdout,*)'mesh file for mpaso domain is ',trim(cvalue)
    end if

    !-----------------------------------------------------------------
    ! Realize the actively coupled fields
    !-----------------------------------------------------------------

    call ocn_realize_fields(gcomp, mesh=Emesh, flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num,                      &
         mastertask = mastertask, lmpicom = lmpicom,           &
         domain = domain_ptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeRealize


  !===============================================================================

  subroutine DataInitialize(gcomp, rc)

    !-----------------------------------------------------------------------
    !  second initializaiton phase of mpaso
    !?!  - initialize ???
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)          :: clock
    type(ESMF_State)          :: importState
    type(ESMF_State)          :: exportState
    type(ESMF_StateItem_Flag) :: itemType
    type(ESMF_StateItem_Flag) :: itemType1
    type(ESMF_StateItem_Flag) :: itemType2
    type(ESMF_TimeInterval)   :: timeStep        ! Model timestep
    type(ESMF_Time)           :: starttime
    character(CL)             :: cvalue
    integer                   :: ocn_cpl_dt
    integer                   :: pop_cpl_dt
    integer                   :: start_ymd
    integer                   :: start_tod
    integer                   :: start_year
    integer                   :: start_day
    integer                   :: start_month
    integer                   :: start_hour
    integer                   :: errorCode       ! error code
    integer                   :: shrlogunit      ! old values
    integer                   :: ocnid
    character(len=*), parameter  :: subname = "ocn_comp_nuopc:(DataInitialize)"
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (stdout)

    !-----------------------------------------------------------------------
    ! register non-standard incoming fields
    !-----------------------------------------------------------------------

    ! query the Component for its importState, exportState and clock
    call ESMF_GridCompGet(gcomp, importState=importState, exportState=exportState, clock=clock, rc=rc)

    !?call ESMF_StateGet(importState, 'Sa_co2prog', itemType, rc=rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?ldriver_has_atm_co2_prog = (itemType /= ESMF_STATEITEM_NOTFOUND)

    !?call ESMF_StateGet(importState, 'Sa_co2diag', itemType, rc=rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?ldriver_has_atm_co2_diag = (itemType /= ESMF_STATEITEM_NOTFOUND)

    !?call ESMF_StateGet(importState, 'Faxa_nhx', itemType1, rc=rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?call ESMF_StateGet(importState, 'Faxa_noy', itemType2, rc=rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?ldriver_has_ndep = ((itemType1 /= ESMF_STATEITEM_NOTFOUND) .or. (itemType2 /= ESMF_STATEITEM_NOTFOUND))

    !?if (ldriver_has_atm_co2_prog) then
    !?   call named_field_register('ATM_CO2_PROG', ATM_CO2_PROG_nf_ind)
    !?endif
    !?if (ldriver_has_atm_co2_diag) then
    !?   call named_field_register('ATM_CO2_DIAG', ATM_CO2_DIAG_nf_ind)
    !?endif
    !?if (ldriver_has_ndep) then
    !?   if (mastertask) write(stdout,'(" using ATM_NHx and ATM_NOy from mediator")')
    !?   call named_field_register('ATM_NHx', ATM_NHx_nf_ind)
    !?   call named_field_register('ATM_NOy', ATM_NOy_nf_ind)
    !?endif

    !?call register_string('pop_init_coupled')
    !?call flushm (stdout)

    !-----------------------------------------------------------------------
    ! second initialization phase of mpaso
    !-----------------------------------------------------------------------

    !?call t_startf ('mpaso_init2')

    !?call pop_init_phase2(errorCode)
    if (rc /= ESMF_Success) then
       call ESMF_LogWrite(trim(subname)//'MPASO_Initialize2: error in mpas_init_phase2',ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    endif
    !?! initialize driver-level flags and timers
    !?call access_time_flag ('stop_now', stop_now)
    !?call access_time_flag ('coupled_ts', cpl_ts)
    !?call get_timer(timer_total,'TOTAL', 1, distrb_clinic%nprocs)

    !?call t_stopf ('mpaso_init2')

    !-----------------------------------------------------------------------
    !  initialize time-stamp information
    !-----------------------------------------------------------------------

    !?call ccsm_char_date_and_time

    !?!-----------------------------------------------------------------------
    !?! check for consistency of pop and sync clock initial time
    !?!-----------------------------------------------------------------------

    !?if (runtype == 'initial') then

    !?   call ESMF_ClockGet( clock, startTime=startTime, rc=rc )
    !?   if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?   call ESMF_TimeGet( startTime, yy=start_year, mm=start_month, dd=start_day, s=start_tod, rc=rc )
    !?   if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?   call shr_cal_ymd2date(start_year,start_month,start_day,start_ymd)

!?!$OMP MASTER
    !?   if (iyear0 /= start_year) then
    !?      call document ('DataInitialize', 'iyear0      ', iyear0)
    !?      call document ('DataInitialize', 'imonth0     ', imonth0)
    !?      call document ('DataInitialize', 'iday0       ', iday0)
    !?      call document ('DataInitialize', 'start_year  ', start_year)
    !?      call document ('DataInitialize', 'start_month ', start_month)
    !?      call document ('DataInitialize', 'start_day   ', start_day)

    !?      ! skip exit_POP if pop is 1 day ahead across a year boundary
    !?      if (.not. (iyear0 == start_year + 1 .and. &
    !?                (imonth0 == 1 .and. start_month == 12) .and. &
    !?                (iday0 == 1   .and. start_day == 31))) then
    !?         call exit_POP(sigAbort,' iyear0 does not match start_year')
    !?      endif
    !?   else if (imonth0 /= start_month) then
    !?      call document ('DataInitialize', 'imonth0     ', imonth0)
    !?      call document ('DataInitialize', 'iday0       ', iday0)
    !?      call document ('DataInitialize', 'start_month ', start_month)
    !?      call document ('DataInitialize', 'start_day   ', start_day)
    !?      ! skip exit_POP if pop is 1 day ahead across a month boundary
    !?      !   this conditional doesn't confirm that start_day is the last day of the month,
    !?      !   only that iday0 is the first day of the month
    !?      if (.not. (imonth0 == start_month + 1 .and. iday0 == 1)) then
    !?         call exit_POP(sigAbort,' imonth0 does not match start_month')
    !?      endif

    !?   else if (iday0 /= start_day) then
    !?      call document ('DataInitialize', 'iday0       ', iday0)
    !?      call document ('DataInitialize', 'start_day   ', start_day)

    !?      ! skip exit_POP if pop is 1 day ahead
    !?      if (.not. (iday0 == start_day + 1)) then
    !?         call exit_POP(sigAbort,' iday0 does not match start_day')
    !?      endif
    !?   end if
!?!$OMP END MASTER
    !?end if
    
    !?!-----------------------------------------------------------------
    !?! Initialize MCT gsmaps and domains
    !?!-----------------------------------------------------------------

    !?call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?read(cvalue,*) ocnid  ! convert from string to integer

    !?call pop_mct_init(ocnid, mpi_communicator_ocn)
    !?if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !?!-----------------------------------------------------------------------
    !?! Initialize flags and shortwave absorption profile
    !?! Note that these cpl_write_xxx flags have no freqency options
    !?! set; therefore, they will retain a default value of .false.
    !?! unless they are explicitly set .true.  at the appropriate times
    !?!-----------------------------------------------------------------------

    !?call init_time_flag('cpl_write_restart',cpl_write_restart, owner = 'DataInitialize')
    !?call init_time_flag('cpl_write_history',cpl_write_history, owner = 'DataInitialize')
    !?call init_time_flag('cpl_write_tavg'   ,cpl_write_tavg,    owner = 'DataInitialize')
    !?call init_time_flag('cpl_diag_global'  ,cpl_diag_global,   owner = 'DataInitialize')
    !?call init_time_flag('cpl_diag_transp'  ,cpl_diag_transp,   owner = 'DataInitialize')

    !?lsmft_avail = .true.
    !?tlast_coupled = c0

    !?!-----------------------------------------------------------------------
    !?! initialize necessary coupling info
    !?!-----------------------------------------------------------------------

    !?call ESMF_ClockGet(clock, timeStep=timeStep, rc=rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?call ESMF_TimeIntervalGet( timeStep, s=ocn_cpl_dt, rc=rc )
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?pop_cpl_dt = seconds_in_day / ncouple_per_day

    !?if (pop_cpl_dt /= ocn_cpl_dt) then
    !?   write(stdout,*)'pop_cpl_dt= ',pop_cpl_dt,' ocn_cpl_dt= ',ocn_cpl_dt
    !?   call exit_POP(sigAbort,'ERROR pop_cpl_dt and ocn_cpl_dt must be identical')
    !?end if

    !-----------------------------------------------------------------------
    ! send export state
    !-----------------------------------------------------------------------

    call ocn_export(exportState, flds_scalar_name, domain_ptr,     &
                    errorCode, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?if (errorCode /= POP_Success) then
    !?   call POP_ErrorPrint(errorCode)
    !?   call exit_POP(sigAbort, 'ERROR in ocn_export')
    !?endif

    !?call State_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, &
    !?     flds_scalar_name, flds_scalar_num, rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?call State_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, &
    !?     flds_scalar_name, flds_scalar_num, rc)
    !?if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?if ( lsend_precip_fact ) then
    !?   call State_SetScalar(precip_fact, flds_scalar_index_precip_factor, exportState, &
    !?        flds_scalar_name, flds_scalar_num, rc)
    !?   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?else
    !?   call State_SetScalar(1._r8, flds_scalar_index_precip_factor, exportState, &
    !?        flds_scalar_name, flds_scalar_num, rc)
    !?   if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !?end if

#if (defined _MEMTRACE)
    if (iam  == 0) then
       lbnum=1
       call memmon_dump_fort('memmon.out','DataInitialize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    !?!-----------------------------------------------------------------------
    !?! Document orbital parameters
    !?!-----------------------------------------------------------------------

    !?if (registry_match('qsw_distrb_iopt_cosz')) then

    !?   call pop_orbital_init(gcomp, stdout, mastertask, rc)
    !?   if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?   call pop_orbital_update(clock, stdout, mastertask, orb_eccen, orb_obliqr, orb_lambm0, orb_mvelpp, rc)
    !?   if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?   write(stdout,*) ' '
    !?   call document ('DataInitialize', 'orb_eccen   ',  orb_eccen)
    !?   call document ('DataInitialize', 'orb_mvelpp  ',  orb_mvelpp)
    !?   call document ('DataInitialize', 'orb_lambm0  ',  orb_lambm0)
    !?   call document ('DataInitialize', 'orb_obliqr  ',  orb_obliqr)
    !?endif

    !-----------------------------------------------------------------------
    ! check whether all Fields in the exportState are "Updated"
    !-----------------------------------------------------------------------

    if (NUOPC_IsUpdated(exportState)) then
      call NUOPC_CompAttributeSet(gcomp, name="InitializeDataComplete", value="true", rc=rc)

      call ESMF_LogWrite("MPASo - Initialize-Data-Dependency SATISFIED!!!", ESMF_LOGMSG_INFO)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    !?!-----------------------------------------------------------------------
    !?! Now document all time flags, last step of pop initialization
    !?!-----------------------------------------------------------------------

    !?call document_time_flags

    !?!-----------------------------------------------------------------------
    !?! Output delimiter to log file
    !?!-----------------------------------------------------------------------

    !?if (mastertask) then
    !?   write(stdout,blank_fmt)
    !?   write(stdout,'(" End of initialization")')
    !?   write(stdout,blank_fmt)
    !?   write(stdout,ndelim_fmt)
    !?   call POP_IOUnitsFlush(POP_stdout)
    !?   call POP_IOUnitsFlush(stdout)
    !?endif

    !?!----------------------------------------------------------------------------
    !?! Reset shr logging to original values
    !?!----------------------------------------------------------------------------

    !?call shr_file_setLogUnit (shrlogunit)

  end subroutine DataInitialize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Run POP for a coupling interval
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    !  local variables
    type(ESMF_State)             :: importState
    type(ESMF_State)             :: exportState
    type(ESMF_Clock)             :: clock
    type(ESMF_Alarm)             :: alarm
    type(ESMF_Time)              :: nextTime
    character(ESMF_MAXSTR)       :: cvalue
    integer                      :: errorCode  ! error flag
    integer                      :: ymd        ! POP2 current date (YYYYMMDD)
    integer                      :: tod        ! POP2 current time of day (sec)
    integer                      :: ymd_sync   ! Sync clock current date (YYYYMMDD)
    integer                      :: tod_sync   ! Sync clcok current time of day (sec)
    integer                      :: lbnum
    integer                      :: ierr
    integer                      :: yr_sync
    integer                      :: mon_sync
    integer                      :: day_sync
    integer                      :: shrlogunit ! old values
    character(CL)          :: message
    logical                      :: first_time = .true.
    logical                      :: debugon
    character(len=*), parameter  :: subname = "ocn_comp_nuopc: (ModelAdvance)"

    type (mpas_pool_type), pointer :: meshPool, statePool, &
                                      forcingPool, &
                                      averagePool, scratchPool
    type (block_type), pointer :: block_ptr
      logical, pointer :: config_write_output_on_startup
      logical, pointer :: config_use_ecosysTracers
      logical, pointer :: config_use_CFCTracers
      logical, pointer :: config_use_activeTracers_surface_restoring
      logical, pointer :: config_use_surface_salinity_monthly_restoring
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND), pointer :: config_sw_absorption_type
      integer iam, SHRLOGLEV
      ! Added for coupling interval initialization
      integer, pointer :: index_avgZonalSSHGradient, index_avgMeridionalSSHGradient
      real (kind=RKIND), dimension(:),   pointer :: filteredSSHGradientZonal, filteredSSHGradientMeridional
      real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
      real (kind=RKIND), pointer :: config_ssh_grad_relax_timescale
      real (kind=RKIND) :: timeFilterFactor
      type (MPAS_timeInterval_type) :: timeStep
      character(len=StrKIND) :: timeStamp, streamName, WCstring
      integer  :: streamDirection
      logical :: streamActive, rstwr
      type (MPAS_Time_Type) :: currTime
      real (kind=RKIND) :: dt, current_wallclock_time

    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    errorCode = ESMF_Success


    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Note this logic triggers off of the component clock rather than the internal pop time
    ! The component clock does not get advanced until the end of the loop - not at the beginning

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

      iam = domain_ptr % dminfo % my_proc_id

      debugOn = .false.
#ifdef MPAS_DEBUG
      debugOn = .true.
#endif
 
      ! Set MPAS Log module instance
      mpas_log_info => domain_ptr % logInfo
      if (debugOn) call mpas_log_write("=== Beginning ocn_run_nuopc ===")

      call mpas_pool_get_config(domain_ptr % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(domain_ptr % configs, 'config_sw_absorption_type', config_sw_absorption_type)

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (ocnLogUnit)

      timeStep = mpas_get_clock_timestep(domain_ptr % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)
      call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

      ! Import state from coupler
       call ocn_import(importState, flds_scalar_name, domain_ptr,  &
                       errorCode, rc)
      ! Ensures MPAS AM write/compute startup steps are performed
      call ocn_analysis_compute_startup(domain_ptr, ierr)

      ! Handle writing initial state
      call mpas_pool_get_config(domain_ptr % configs, 'config_write_output_on_startup', config_write_output_on_startup)
      if (config_write_output_on_startup) then
         call mpas_stream_mgr_write(domain_ptr % streamManager, 'output', forceWriteNow=.true., ierr=ierr)

         ! Reset config to false, so we don't write the state every coupling interval.
         config_write_output_on_startup = .false.
      end if

      ! Initialize time average fields
      block_ptr => domain_ptr % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call ocn_time_average_coupled_init(forcingPool)
         block_ptr => block_ptr % next
      end do

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      ! This integration loop continues for a single coupling interval.
      do while (.not. mpas_is_alarm_ringing(domain_ptr % clock, coupleAlarmID, ierr=ierr))
         call mpas_stream_mgr_read(domain_ptr % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

         ! Build forcing arrays.
         if (debugOn) call mpas_log_write('   Building forcing arrays')
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

            call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)

            call mpas_timer_start("land_ice_build_arrays", .false.)
            call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                          forcingPool, scratchPool, statePool, dt, ierr)
            call mpas_timer_stop("land_ice_build_arrays")

            call ocn_frazil_forcing_build_arrays(domain_ptr, meshPool, forcingPool, statePool, ierr)

            ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
            if (config_use_Redi.or.config_use_GM) then
              call ocn_gm_compute_Bolus_velocity(statePool, meshPool, scratchPool, timeLevelIn=1)
            endif

            block_ptr => block_ptr % next
         end do
         if (debugOn) call mpas_log_write('   Finished building forcing arrays')

         ! Read forcing streams for next timestep
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Reading forcing streams')
         endif

         ! read next time level data required for variable shortwave
         call mpas_timer_start('io_shortwave',.false.)
         call ocn_get_shortWaveData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .false.)
         call mpas_timer_stop('io_shortwave')

         ! read next time level data required for ecosys forcing
         call mpas_pool_get_config(domain_ptr % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
         if (config_use_ecosysTracers) then
           call mpas_timer_start('io_ecosys',.false.)
           call ocn_get_ecosysData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .false.)
           call mpas_timer_stop('io_ecosys')
         endif

         ! read next time level data required for CFC forcing
         call mpas_pool_get_config(domain_ptr % configs, 'config_use_CFCTracers', config_use_CFCTracers)
         if (config_use_CFCTracers) then
           call mpas_timer_start('io_CFC',.false.)
           call ocn_get_CFCData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .false.)
           call mpas_timer_stop('io_CFC')
         endif

         ! read next time level data required for monthly surface salinity restoring
         call mpas_pool_get_config(domain_ptr % configs, 'config_use_activeTracers_surface_restoring',  &
                                                      config_use_activeTracers_surface_restoring)
         call mpas_pool_get_config(domain_ptr % configs, 'config_use_surface_salinity_monthly_restoring',  &
                                                      config_use_surface_salinity_monthly_restoring)
         if (config_use_activeTracers_surface_restoring .and.  config_use_surface_salinity_monthly_restoring) then
             if ( mpas_is_alarm_ringing(domain_ptr % clock, 'salinityDataReadAlarm', ierr=ierr) ) then
                 call mpas_reset_clock_alarm(domain_ptr % clock, 'salinityDataReadAlarm', ierr=ierr)
                 call mpas_timer_start('io_monthly_surface_salinity',.false.)
                 call ocn_get_surfaceSalinityData(domain_ptr % streamManager, domain_ptr, domain_ptr % clock, .false.)
                 call mpas_timer_stop('io_monthly_surface_salinity')
             endif
         endif
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Finished reading forcing streams')
         endif

         itimestep = itimestep + 1
         call mpas_advance_clock(domain_ptr % clock)

         currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         ! write time stamp at every step
         call mpas_dmpar_get_time(current_wallclock_time)
         write (WCstring,'(F18.3)') current_wallclock_time
         call mpas_log_write(trim(timeStamp) // '  WC time:' // WCstring)

         ! pre-timestep analysis computation
         if (debugOn) call mpas_log_write('   Starting analysis precompute', masterOnly=.true.)
         call ocn_analysis_precompute(domain_ptr, ierr)
         if (debugOn) call mpas_log_write('   Finished analysis precompute', masterOnly=.true.)

         if (debugOn) call mpas_log_write('   Performing forward update')
         call mpas_timer_start("time integration", .false.)
         call ocn_timestep(domain_ptr, dt, timeStamp)
         call mpas_timer_stop("time integration")
         if (debugOn) call mpas_log_write('   Finished performing forward update')

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do

         if (debugOn) call mpas_log_write('   Computing analysis members')
         call ocn_analysis_compute(domain_ptr, ierr) 
         if (iam==0.and.debugOn) then
            call mpas_log_write( '   Finished computing analysis members')
            call mpas_log_write('   Preparing analysis members for restart')
         endif
         call ocn_analysis_restart(domain_ptr, ierr) 
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Finished preparing analysis members for restart')
            call mpas_log_write('   Writing analysis member output')
         endif
         call ocn_analysis_write(domain_ptr, ierr)
         if (debugOn) call mpas_log_write('   Finished writing analysis member output')

         ! Reset any restart alarms to prevent restart files being written without the coupler requesting it.
         if (debugOn) call mpas_log_write('   Resetting restart alarms')
         call mpas_stream_mgr_begin_iteration(domain_ptr % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain_ptr % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID=streamName, ierr=ierr)
            end if
         end do
         if (debugOn) call mpas_log_write('   Finished resetting restart alarms')

         if (debugOn) call mpas_log_write('   Writing output streams')
         call mpas_stream_mgr_write(domain_ptr % streamManager, streamID='output', ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, streamID='output', ierr=ierr)

         call mpas_stream_mgr_write(domain_ptr % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Finished writing output streams')
            call mpas_log_write('   Validating ocean state')
         endif

         call ocn_validate_state(domain_ptr, timeLevel=1)
         if (debugOn) call mpas_log_write('   Completed validating ocean state')

         if (debugOn) call mpas_log_write('Completed timestep '//trim(timeStamp))
      end do

      ! update coupled variables that get calculated on coupling intervals
      ! NOTE: could be moved to subroutine
      ! time filter ssh gradient
      block_ptr => domain_ptr % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
         call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)
         call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
         call mpas_pool_get_array(forcingPool, 'filteredSSHGradientZonal', filteredSSHGradientZonal)
         call mpas_pool_get_array(forcingPool, 'filteredSSHGradientMeridional', filteredSSHGradientMeridional)
         call mpas_pool_get_config(domain_ptr % configs, 'config_ssh_grad_relax_timescale', &
                                                      config_ssh_grad_relax_timescale)
         if (config_ssh_grad_relax_timescale < real(ocn_cpl_dt,RKIND)) then
            filteredSSHGradientZonal = avgSSHGradient(index_avgZonalSSHGradient, :)
            filteredSSHGradientMeridional = avgSSHGradient(index_avgMeridionalSSHGradient, :)
         else
            timeFilterFactor = real(ocn_cpl_dt,RKIND) / config_ssh_grad_relax_timescale
            filteredSSHGradientZonal = filteredSSHGradientZonal * (1.0_RKIND - timeFilterFactor) + &
                 avgSSHGradient(index_avgZonalSSHGradient, :) * timeFilterFactor
            filteredSSHGradientMeridional = filteredSSHGradientMeridional * (1.0_RKIND - timeFilterFactor) + &
                 avgSSHGradient(index_avgMeridionalSSHGradient, :) * timeFilterFactor
         endif
         block_ptr => block_ptr % next
      end do

      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval

      if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
         if (ChkErr(rc,__LINE__,u_FILE_u)) return
         call ESMF_AlarmRingerOff( alarm, rc=rc )
         if (ChkErr(rc,__LINE__,u_FILE_u)) return

         if(trim(config_sw_absorption_type)=='ohlmann00') call ocn_shortwave_forcing_write_restart(domain_ptr)

         if (config_use_ecosysTracers) call ocn_ecosys_forcing_write_restart(domain_ptr)

         if (config_use_CFCTracers) call ocn_CFC_forcing_write_restart(domain_ptr)

         if (config_use_activeTracers_surface_restoring .and. config_use_surface_salinity_monthly_restoring)  &
             call ocn_salinity_restoring_forcing_write_restart(domain_ptr)
         ! Write a restart file, because the coupler asked for it.
         if (debugOn) call mpas_log_write('Writing restart streams')
         call mpas_stream_mgr_begin_iteration(domain_ptr % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain_ptr % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               if (debugOn) call mpas_log_write('   Writing stream ' // trim(streamName))
               call mpas_stream_mgr_write(domain_ptr % streamManager, forceWriteNow=.true., streamID=streamName, ierr=ierr)
               if (debugOn) call mpas_log_write('   Finished writing stream ' // trim(streamName))
            end if
         end do

         if ( iam == 0 ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if

         if (debugOn) call mpas_log_write('Finished writing restart streams')
      endif

      ! Export state to coupler
      if (debugOn) call mpas_log_write('Exporting ocean state')
      call ocn_export(exportState, flds_scalar_name, domain_ptr,     &
                      errorCode, rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      if (debugOn) call mpas_log_write('Finished exporting ocean state')

!DD adapt this statement to the esmf clock
!?!?      call check_clocks_sync(domain_ptr % clock, Eclock, ierr)

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

      call mpas_log_write('=== Completed coupling interval in ocn_run_nuopc. ===', flushNow=.true.)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelSetRunClock)'
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelCheckImport(model, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)  :: clock
    type(ESMF_Time)   :: currTime
    integer           :: yy  ! current date (YYYYMMDD)
    integer           :: mon ! current month
    integer           :: day ! current day
    integer           :: tod ! current time of day (sec)
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currTime, yy=yy, mm=mon, dd=day, s=tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !?if (mastertask) then
    !?   write(stdout,*)' CheckImport pop2 year = ',yy
    !?   write(stdout,*)' CheckImport pop2 mon  = ',mon
    !?   write(stdout,*)' CheckImport pop2 day  = ',day
    !?   write(stdout,*)' CheckImport pop2 tod  = ',tod
    !?end if

  end subroutine ModelCheckImport


  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)

    !--------------------------------
    ! MPASO finalization that shuts down MPASO gracefully (we hope).
    ! Exits the message environment and checks for successful execution.
    ! --------------------------------


    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(ocn_comp_nuopc) ',40a)"
    character(*), parameter :: F91   = "('(ocn_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname = 'ocn_comp_nuopc:(ModelFinalize) '
    !--------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (my_task == 0) then
       write(stdout,F91)
       write(stdout,F00) 'MPASO: end of main integration loop'
       write(stdout,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !! local variables
    !integer  :: ierr              ! error code
    !character(len=*),parameter :: subname='ocn_comp_nuopc:(ModelFinalize)'
    !integer  :: iam, SHRLOGLEV, SHRLOGUNIT
    !!--------------------------------

    !rc = ESMF_SUCCESS
    !if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

!-----------------------------------------------------------------------
    !iam = domain_ptr % dminfo % my_proc_id

    !! Set MPAS Log module instance
    !mpas_log_info => domain_ptr % logInfo

    !! Setup I/O logs
    !call shr_file_getLogUnit (shrlogunit)
    !call shr_file_getLogLevel(shrloglev)
    !call shr_file_setLogUnit (ocnLogUnit)

    !! Finalize MPASO
    !iErr = domain_ptr % core % core_finalize(domain_ptr)
    !if ( iErr /= 0 ) then
    !   call mpas_log_write('Core finalize failed for core ' // trim(domain_ptr % core % coreName), MPAS_LOG_CRIT)
    !end if

    !call mpas_timer_write()

    !call MPAS_stream_mgr_finalize(domain_ptr % streamManager)

    !call mpas_log_finalize(iErr)
    !if ( iErr /= 0 ) then
    !   write(ocnLogUnit,*) 'ERROR: log finalize failed for core ' // trim(domain_ptr % core % coreName)
    !   call mpas_dmpar_abort(domain_ptr % dminfo)
    !end if

    !call mpas_framework_finalize(domain_ptr % dminfo, domain_ptr, io_system)

    !! Reset I/O logs
    !call shr_file_setLogUnit (shrlogunit)
    !call shr_file_setLogLevel(shrloglev)

    !if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

   subroutine add_stream_attributes(domain)

      use mpas_stream_manager, only : MPAS_stream_mgr_add_att

      implicit none

      type (domain_type), intent(inout) :: domain

      type (MPAS_Pool_iterator_type) :: itr
      integer, pointer :: intAtt
      logical, pointer :: logAtt
      character (len=StrKIND), pointer :: charAtt
      real (kind=RKIND), pointer :: realAtt
      character (len=StrKIND) :: histAtt

      integer :: local_ierr

      if (domain % dminfo % nProcs < 10) then
          write(histAtt, '(A,I1,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100) then
          write(histAtt, '(A,I2,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 1000) then
          write(histAtt, '(A,I3,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 10000) then
          write(histAtt, '(A,I4,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100000) then
          write(histAtt, '(A,I5,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else
          write(histAtt, '(A,I6,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      end if
     
      call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'source', domain % core % source)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', domain % core % Conventions)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', domain % core % git_version)

      call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'is_periodic', domain % is_periodic)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'x_period', domain % x_period)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'y_period', domain % y_period)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
!     call MPAS_stream_mgr_add_att(domain % streamManager, 'history', domain % history)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', histAtt)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain %  parent_id)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)

      call mpas_pool_begin_iteration(domain % configs)

      do while (mpas_pool_get_next_member(domain % configs, itr))

         if ( itr % memberType == MPAS_POOL_CONFIG) then

            if ( itr % dataType == MPAS_POOL_REAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, realAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_INTEGER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, intAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, charAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
               if (logAtt) then
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'YES', ierr=local_ierr)
               else
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'NO', ierr=local_ierr)
               end if
            end if

          end if
      end do

   end subroutine add_stream_attributes


end module ocn_comp_nuopc
