module ocn_import_export

  use ESMF
  use NUOPC
  use NUOPC_Model
  use mpas_kind_types, only : r8 => r8kind
  use mpas_derived_types
  use mpas_field_routines
  use mpas_pool_routines
  use mpas_dmpar, only : mpas_dmpar_exch_halo_field, mpas_dmpar_sum_real
  use ocn_constants, only : RHO_SW, CP_SW, T0_KELVIN
  use ocn_config, only : CONFIG_FRAZIL_HEAT_OF_FUSION
  use ocn_equation_of_state, only : ocn_freezing_temperature
  
  use shr_kind_mod,          only: cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod,         only: shr_const_spval, radius=>SHR_CONST_REARTH, pi=>shr_const_pi
  use shr_string_mod,        only: shr_string_listGetNum, shr_string_listGetName
  use shr_mpi_mod,           only: shr_mpi_min, shr_mpi_max
  use shr_sys_mod 
  use nuopc_shr_methods,     only: chkerr

  implicit none
  public

  public  :: ocn_advertise_fields
  public  :: ocn_realize_fields
  public  :: ocn_import
  public  :: ocn_export

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_FldChk

  integer, public :: ocn_cpl_dt    ! length of coupling interval in seconds - set by coupler/ESMF
  
  ! Private module data

  type fld_list_type
    character(len=128) :: stdname
    integer :: ungridded_lbound = 0
    integer :: ungridded_ubound = 0
  end type fld_list_type

  integer, parameter       :: fldsMax = 100
  integer                  :: fldsToOcn_num = 0
  integer                  :: fldsFrOcn_num = 0
  type (fld_list_type)     :: fldsToOcn(fldsMax)
  type (fld_list_type)     :: fldsFrOcn(fldsMax)

  ! area correction factors for fluxes send and received from mediator
  real(r8), allocatable :: mod2med_areacor(:) ! ratios of model areas to input mesh areas
  real(r8), allocatable :: med2mod_areacor(:) ! ratios of input mesh areas to model areas

  interface state_getfldptr
     module procedure state_getfldptr_1d
     module procedure state_getfldptr_2d
  end interface state_getfldptr

  integer     , parameter :: dbug = 0        ! i/o debug messages
  integer, parameter  :: stdout = 6
  character(*), parameter :: u_FILE_u = &
       __FILE__

!==============================================================================
contains
!==============================================================================

  subroutine ocn_advertise_fields(gcomp, importState, exportState, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    integer       :: n
    character(CS) :: stdname
    character(CS) :: cvalue
    character(CS) :: cname
    integer       :: ice_ncat
    logical       :: flds_i2o_per_cat  ! .true. => select per ocn thickness category
    character(len=*), parameter :: subname='(ocn_advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)


    !-----------------
    ! advertise import fields
    !-----------------

    call fldlist_add(fldsToOcn_num, fldsToOcn, trim(flds_scalar_name))

    ! from ice
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_ifrac')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Si_bpress')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_melth')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_meltw')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Fioi_salt')

    ! from river
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofl')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_rofi')

    ! from mediator
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'So_duu10n')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_tauy')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_taux')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lat')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_sen')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_lwup')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_evap')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Foxx_swnet')

    ! from atmosphere
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Sa_pslv')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_lwdn')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_snow')
    call fldlist_add(fldsToOcn_num, fldsToOcn, 'Faxa_rain')

    do n = 1,fldsToOcn_num
       call NUOPC_Advertise(importState, standardName=fldsToOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    !-----------------
    ! advertise export fields
    !-----------------

    call fldlist_add(fldsFrOcn_num, fldsFrOcn, trim(flds_scalar_name))
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_omask')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_t')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_u')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_v')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_s')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdx')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_dhdy')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'So_bldepth')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Fioo_q')
    call fldlist_add(fldsFrOcn_num, fldsFrOcn, 'Fioo_frazil')

    do n = 1,fldsFrOcn_num
       call NUOPC_Advertise(exportState, standardName=fldsFrOcn(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_advertise_fields


  !==============================================================================
  subroutine ocn_realize_fields(gcomp, mesh, flds_scalar_name, flds_scalar_num, &
                                mastertask, lmpicom, domain, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    type(ESMF_Mesh)  , intent(in)  :: mesh
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: flds_scalar_num
    integer          , intent(in)  :: lmpicom  ! the ocean mpi communicator
    logical          , intent(in)  :: mastertask           
    type (domain_type), pointer, intent(in) :: domain
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)      :: importState
    type(ESMF_State)      :: exportState
    type(ESMF_Field)      :: lfield
    integer               :: spatialDim
    integer               :: numOwnedElements
    integer               :: i,j,iblock,n
    real(r8), allocatable :: mesh_areas(:)
    real(r8), allocatable :: model_areas(:)
    real(r8), pointer     :: dataptr(:)
    integer               :: num_ocn
    real(r8)              :: max_mod2med_areacor
    real(r8)              :: max_med2mod_areacor
    real(r8)              :: min_mod2med_areacor
    real(r8)              :: min_med2mod_areacor
    real(r8)              :: max_mod2med_areacor_glob
    real(r8)              :: max_med2mod_areacor_glob
    real(r8)              :: min_mod2med_areacor_glob
    real(r8)              :: min_med2mod_areacor_glob
    real(r8), pointer     :: ownedElemCoords(:)
    real(r8), pointer     :: latModel(:), latMesh(:)
    real(r8), pointer     :: lonModel(:), lonMesh(:)
    real(r8)              :: diff_lon
    real(r8)              :: diff_lat

    integer :: iCell, nCells
    integer, dimension(:), pointer :: nCellsArray
    real (kind=RKIND), dimension(:), pointer :: lonCell,                       &
                                                latCell,                       &
                                                areaCell 

    type (mpas_pool_type), pointer :: meshPool
    type (block_type), pointer :: block
    integer gcell, cell_offset
    real(r8) :: rad2deg = 180._r8/pi

    character(len=*), parameter :: subname='(ocn_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    ! Get import and export states
    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    ! Realize import and export states
    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrOcn, &
         numflds=fldsFrOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MPASO_Export',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToOcn, &
         numflds=fldsToOcn_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MPASO_Import',&
         mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine mesh lats and lons
    call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numownedelements))
    allocate(lonMesh(numOwnedElements))
    allocate(latMesh(numOwnedElements))
    allocate(lonModel(numOwnedElements))
    allocate(latModel(numOwnedElements))
    call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       lonMesh(n) = ownedElemCoords(2*n-1)
       latMesh(n) = ownedElemCoords(2*n)
    end do
    deallocate(ownedElemCoords)

    ! Compare mesh lats/lons to model generated lats/lons
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )
       call mpas_pool_get_array(meshPool, 'lonCell', lonModel)  
       call mpas_pool_get_array(meshPool, 'latCell', latModel)  
       
       do iCell = 1, nCells
          gcell = iCell + cell_offset

          diff_lon = abs(lonMesh(gcell) - lonModel(iCell)*rad2deg)
          if ( (diff_lon > 1.e2  .and. abs(diff_lon - 360.) > 1.e-1) .or.&
               (diff_lon > 1.e-3 .and. diff_lon < 1._r8) ) then
             write(stdout,'(a,i6,2(f21.13,3x),d21.5)') &
                 'ERROR: MPASO n, lonMesh, lonModel, diff_lon = ',&
                 gcell,lonMesh(gcell),lonModel(icell)*rad2deg, diff_lon
             call shr_sys_abort()
          end if
          if (abs(latMesh(gcell) - latModel(icell)*rad2deg) > 1.e-1) then
             write(stdout,'(a,i6,2(f21.13,3x),d21.5)') &
                  'ERROR: MPASO n, latMesh, latModel, diff_lat = ', &
                  gcell,latMesh(gcell),latModel(icell)*rad2deg, abs(latMesh(gcell)-latModel(iCell)*rad2deg)
             call shr_sys_abort()
          end if
       end do

       cell_offset = cell_offset + nCells
       
       block => block % next
    end do

    ! Determine mesh areas used in regridding
    lfield = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_areas(numOwnedElements))
    mesh_areas(:) = dataptr(:)
    call ESMF_FieldDestroy(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine flux correction factors (module variables)
    allocate(model_areas(numOwnedElements))
    allocate(mod2med_areacor(numOwnedElements))
    allocate(med2mod_areacor(numOwnedElements))
    mod2med_areacor(:) = 1._r8
    med2mod_areacor(:) = 1._r8
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)  
       
       do iCell = 1, nCells
          gcell = iCell + cell_offset

          model_areas(gcell) = areaCell(icell)/(radius*radius)   
          mod2med_areacor(gcell) = model_areas(gcell) / mesh_areas(gcell)
          med2mod_areacor(gcell) = mesh_areas(gcell) / model_areas(gcell)
       end do

       cell_offset = cell_offset + nCells
       
       block => block % next
    end do
    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, lmpicom)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, lmpicom)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, lmpicom)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, lmpicom)

    if (mastertask) then
       write(stdout,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'MPASO'
       write(stdout,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'MPASO'
    end if

    deallocate(model_areas)
    deallocate(mesh_areas)

  end subroutine ocn_realize_fields


  !==============================================================================
  subroutine ocn_import( importState, flds_scalar_name, domain,    &
                         errorCode, rc )

    !-----------------------------------------------------------------------
    ! swnet  -- net short-wave heat flux                 (W/m2   )
    ! lwup   -- longwave radiation (up)                  (W/m2   )
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)   , intent(in)  :: importState
    character(len=*)   , intent(in)  :: flds_scalar_name
    type (domain_type), pointer, intent(in) :: domain
    integer            , intent(out) :: errorCode
    integer            , intent(out) :: rc

    ! local variables
    character (cl) :: label,  message
    integer              :: i,j,k,n,ncol,iblock,nfld,nf
    real (r8)            :: m2percm2, gsum
    real (r8), pointer   :: dataptr1d(:)
    real (r8), pointer   :: dataptr2d(:,:)
    ! from mediator (virtual ocn)
    real (r8), pointer   :: foxx_swnet(:)
    real (r8), pointer   :: foxx_taux(:)
    real (r8), pointer   :: foxx_tauy(:)
    real (r8), pointer   :: foxx_lwup(:)
    real (r8), pointer   :: foxx_sen(:)
    real (r8), pointer   :: foxx_lat(:)
    real (r8), pointer   :: foxx_evap(:)
    real (r8), pointer   :: foxx_rofl(:)
    real (r8), pointer   :: foxx_rofi(:)
    real (r8), pointer   :: so_duu10n(:)
    ! from atm
    real (r8), pointer   :: sa_pslv(:)
    real (r8), pointer   :: faxa_rain(:)
    real (r8), pointer   :: faxa_snow(:)
    real (r8), pointer   :: faxa_lwdn(:)
    ! from ice
    real (r8), pointer   :: fioi_meltw(:)
    real (r8), pointer   :: fioi_melth(:)
    real (r8), pointer   :: fioi_salt(:)
    real (r8), pointer   :: Si_ifrac(:)
    real (r8), pointer   :: Si_bpress(:)
    !
    integer              :: fieldCount
    character (cl) :: fldname
    type(ESMF_StateItem_Flag) :: itemflag
!?#ifdef _HIRES
!?    real (r8)            :: qsw_eps = -1.e-3_r8
!?#else
    real (r8)            :: qsw_eps = 0._r8
!?#endif

    ! local variables - mpas names
    integer :: iCell, nCells
    integer, dimension(:), pointer :: nCellsArray
    type (mpas_pool_type), pointer :: meshPool, forcingPool
    type (block_type), pointer :: block
    integer gcell, cell_offset
    real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional, &
                                                evaporationFlux, rainFlux,             &
                                                snowFlux, longWaveHeatFluxDown,        &
                                                latentHeatFlux, sensibleHeatFlux,      &
                                                riverRunoffFlux, iceRunoffFlux,        &
                                                shortWaveHeatFlux, longWaveHeatFluxUp, &
                                                seaIceFreshWaterFlux, seaIceHeatFlux,  &
                                                seaIceSalinityFlux, atmosphericPressure, &
                                                iceFraction, removedRiverRunoffFlux,  &
                                                removedIceRunoffFlux, seaIcePressure, latCell

    type (field1DReal),         pointer :: windStressZonalField, windStressMeridionalField, &
                                           evaporationFluxField, rainFluxField,             &
                                           snowFluxField, longWaveHeatFluxDownField,        &
                                           latentHeatFluxField, sensibleHeatFluxField,      &
                                           riverRunoffFluxField, iceRunoffFluxField,        &
                                           shortWaveHeatFluxField, longWaveHeatFluxUpField, &
                                           seaIceFreshWaterFluxField, seaIceHeatFluxField,  &
                                           seaIceSalinityFluxField, atmosphericPressureField, &
                                           iceFractionField, seaIcePressureField

    character (cl), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_import)'
    logical, pointer :: config_remove_AIS_coupler_runoff
    real (kind=RKIND), pointer :: totalRemovedRiverRunoffFlux, totalRemovedIceRunoffFlux
    real (kind=RKIND) :: removedRiverRunoffFluxThisProc, removedIceRunoffFluxThisProc
    real (kind=RKIND) :: removedRiverRunoffFluxReduced, removedIceRunoffFluxReduced
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !?! NOTE: if there are code changes associated with changing the names or
    !?!       the number of fluxes received from the coupler, then subroutine
    !?!       update_ghost_cells_coupler_fluxes will need to be modified also

    !?! NOTE : RCALCT acts as a KMT mask, its 1 if KMT>=1 and 0 otherwise

    !-----------------------------------------------------------------------
    ! from mediator (virtual ocean)
    !-----------------------------------------------------------------------

    ! zonal wind and meridional wind stress  (W/m2)
    call state_getfldptr(importState, 'Foxx_taux', foxx_taux, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, 'Foxx_tauy', foxx_tauy, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! evaporation flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_evap', foxx_evap, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! sensible heat flux (W/m2)
    call state_getfldptr(importState, 'Foxx_lat', foxx_lat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! sensible heat flux (W/m2)
    call state_getfldptr(importState, 'Foxx_sen', foxx_sen, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! long wave up flux  (W/m2)
    call state_getfldptr(importState, 'Foxx_lwup', foxx_lwup, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! shortwave net flux  (W/m2)
    call state_getfldptr(importState, 'Foxx_swnet', foxx_swnet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! 10m wind speed squared (m^2/s^2)
    call state_getfldptr(importState, 'So_duu10n', so_duu10n, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
 
    !-----------------------------------------------------------------------
    ! from atmosphere
    !-----------------------------------------------------------------------

    ! sea-level pressure (Pa)
    call state_getfldptr(importState, 'Sa_pslv', sa_pslv, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! snow
    call state_getfldptr(importState, 'Faxa_snow', faxa_snow, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! rain
    call state_getfldptr(importState, 'Faxa_rain', faxa_rain, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! longwave radiation (down) (W/m2)
    call state_getfldptr(importState, 'Faxa_lwdn', faxa_lwdn, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------------
    ! from sea-ice
    !-----------------------------------------------------------------------

    ! ice fraction
    call state_getfldptr(importState, 'Si_ifrac', si_ifrac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ice base pressure
    call state_getfldptr(importState, 'Si_bpress', si_bpress, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! heat flux from sea ice snow & ice melt (W/m2)
    call state_getfldptr(importState, 'Fioi_melth', fioi_melth, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! salt from sea ice (kg(salt)/m2/s)
    call state_getfldptr(importState, 'Fioi_salt', fioi_salt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! snow melt flux from sea ice (kg/m2/s)
    call state_getfldptr(importState, 'Fioi_meltw', fioi_meltw, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------------
    ! from wave
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! from river
    !-----------------------------------------------------------------------

    ! liquid runoff flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_rofl', foxx_rofl, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! ice runoff flux (kg/m2/s)
    call state_getfldptr(importState, 'Foxx_rofi', foxx_rofi, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------------------------------------------------
    ! copy into component variables
    !-----------------------------------------------------------------------
    removedRiverRunoffFluxThisProc = 0.0_RKIND
    removedIceRunoffFluxThisProc = 0.0_RKIND
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )

       ! from mediator
       call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
       call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)
       call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
       call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
       call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
       call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
       call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)

       ! from atmosphere
       call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
       call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
       call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
       call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)

       ! from seaice
       call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
       call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
       call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
       call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
       call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)

       ! from river
       call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
       call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)

       call mpas_pool_get_config(domain % configs, 'config_remove_AIS_coupler_runoff', config_remove_AIS_coupler_runoff)
       if (config_remove_AIS_coupler_runoff) then
          call mpas_pool_get_array(forcingPool, 'removedRiverRunoffFlux', removedRiverRunoffFlux)
          call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFlux)
          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          ! Initialize these fields
          removedRiverRunoffFlux(:) = 0.0_RKIND
          removedIceRunoffFlux(:) = 0.0_RKIND
       endif

       do iCell = 1, nCells
          gcell = iCell + cell_offset
             
          windStressZonal(iCell)      = foxx_taux (gcell) * med2mod_areacor(gcell)
          windStressMeridional(iCell) = foxx_tauy (gcell) * med2mod_areacor(gcell)
          evaporationFlux(iCell)      = foxx_evap (gcell) * med2mod_areacor(gcell)
          latentHeatFlux(iCell)       = foxx_lat  (gcell) * med2mod_areacor(gcell)
          sensibleHeatFlux(iCell)     = foxx_sen  (gcell) * med2mod_areacor(gcell)
          riverRunoffFlux(iCell)      = foxx_rofl (gcell) * med2mod_areacor(gcell)
           if (config_remove_AIS_coupler_runoff) then
              if (latCell(i) < -1.04719666667_RKIND) then ! 60S in radians
                 removedRiverRunoffFlux(i) = riverRunoffFlux(i)
                 riverRunoffFlux(i) = 0.0_RKIND
                 removedRiverRunoffFluxThisProc = removedRiverRunoffFluxThisProc + removedRiverRunoffFlux(i)
               endif
           endif
          iceRunoffFlux(iCell)        = foxx_rofi (gcell) * med2mod_areacor(gcell)
           if(iceRunoffFlux(iCell) < 0.0_RKIND) then
               call shr_sys_abort ('Error: incoming rofi_F is negative')
           end if
           if (config_remove_AIS_coupler_runoff) then
              if (latCell(iCell) < -1.04719666667_RKIND) then ! 60S in radians
                 removedIceRunoffFlux(iCell) = iceRunoffFlux(iCell)
                 iceRunoffFlux(iCell) = 0.0_RKIND
                 removedIceRunoffFluxThisProc = removedIceRunoffFluxThisProc + removedIceRunoffFlux(iCell)
              endif
           endif
          shortWaveHeatFlux(iCell)    = foxx_swnet(gcell) * med2mod_areacor(gcell)
          longWaveHeatFluxUp(iCell)   = foxx_lwup (gcell) * med2mod_areacor(gcell)
          atmosphericPressure(iCell)  = Sa_pslv   (gcell) * med2mod_areacor(gcell)
          rainFlux(iCell)             = faxa_rain (gcell) * med2mod_areacor(gcell)
          snowFlux(iCell)             = faxa_snow (gcell) * med2mod_areacor(gcell)
          longWaveHeatFluxDown(iCell) = faxa_lwdn (gcell) * med2mod_areacor(gcell)
          seaIceFreshWaterFlux(iCell) = fioi_meltw(gcell) * med2mod_areacor(gcell)
          seaIceHeatFlux(iCell)       = fioi_melth(gcell) * med2mod_areacor(gcell)
          seaIceSalinityFlux(iCell)   = fioi_salt (gcell) * med2mod_areacor(gcell)
          iceFraction(iCell)          = Si_ifrac  (gcell) * med2mod_areacor(gcell)
          seaIcePressure(iCell)       = Si_bpress (gcell) * med2mod_areacor(gcell)
       end do

       if (ANY(shortWaveHeatFlux < qsw_eps)) then
         do iCell = 1, nCells
           gcell = iCell + cell_offset
           write(stdout,*)'ERROR: gcell,shortWaveHeatFlux = ',&
                      gcell,shortWaveHeatFlux(icell)
         enddo
         call shr_sys_abort('(set_surface_forcing) ERROR: SHF_QSW < qsw_eps in set_surface_forcing')
       endif

       cell_offset = cell_offset + nCells
       
       block => block % next
    end do

    !-----------------------------------------------------------------------
    ! update ghost cells for fluxes received from the mediator
    !-----------------------------------------------------------------------
      call mpas_pool_get_field(forcingPool, 'windStressZonal', windStressZonalField)
      call mpas_pool_get_field(forcingPool, 'windStressMeridional', windStressMeridionalField)
      call mpas_pool_get_field(forcingPool, 'latentHeatFlux', latentHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'sensibleHeatFlux', sensibleHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUpField)
      call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDownField)
      call mpas_pool_get_field(forcingPool, 'evaporationFlux', evaporationFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceHeatFlux', seaIceHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
      call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
      call mpas_pool_get_field(forcingPool, 'atmosphericPressure', atmosphericPressureField)
      call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIcePressure', seaIcePressureField)
      if ( windStressZonalField % isActive ) &
         call mpas_dmpar_exch_halo_field(windStressZonalField)
      if ( windStressMeridionalField % isActive ) &
         call mpas_dmpar_exch_halo_field(windStressMeridionalField)
      if ( latentHeatFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(latentHeatFluxField)
      if ( sensibleHeatFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(sensibleHeatFluxField)
      if ( longWaveHeatFluxUpField % isActive ) &
         call mpas_dmpar_exch_halo_field(longWaveHeatFluxUpField)
      if ( longWaveHeatFluxDownField % isActive ) &
         call mpas_dmpar_exch_halo_field(longWaveHeatFluxDownField)
      if ( evaporationFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(evaporationFluxField)
      if ( seaIceHeatFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(seaIceHeatFluxField)
      if ( snowFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(snowFluxField)
      if ( seaIceFreshWaterFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(seaIceFreshWaterFluxField)
      if ( seaIceSalinityFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(seaIceSalinityFluxField)
      if ( riverRunoffFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(riverRunoffFluxField)
      if ( iceRunoffFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(iceRunoffFluxField)
      if ( shortWaveHeatFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(shortWaveHeatFluxField)
      if ( rainFluxField % isActive ) &
         call mpas_dmpar_exch_halo_field(rainFluxField)
      if ( atmosphericPressureField % isActive ) &
         call mpas_dmpar_exch_halo_field(atmosphericPressureField)
      if ( iceFractionField % isActive ) &
         call mpas_dmpar_exch_halo_field(iceFractionField)
      if ( seaIcePressureField % isActive ) &
         call mpas_dmpar_exch_halo_field(seaIcePressureField)

   ! global sum of removed runoff
   if (config_remove_AIS_coupler_runoff) then
      call MPAS_dmpar_sum_real(domain % dminfo, removedRiverRunoffFluxThisProc, removedRiverRunoffFluxReduced)
      call MPAS_dmpar_sum_real(domain % dminfo, removedIceRunoffFluxThisProc, removedIceRunoffFluxReduced)
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

         call mpas_pool_get_array(forcingPool, 'totalRemovedRiverRunoffFlux', totalRemovedRiverRunoffFlux)
         call mpas_pool_get_array(forcingPool, 'totalRemovedIceRunoffFlux', totalRemovedIceRunoffFlux)
         totalRemovedRiverRunoffFlux = removedRiverRunoffFluxReduced
         totalRemovedIceRunoffFlux = removedIceRunoffFluxReduced

         block => block % next
      end do
   endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_import

  !==============================================================================

  subroutine ocn_export(exportState, flds_scalar_name, domain,    &
                        errorCode, rc)

    !-----------------------------------------------------------------------
    ! Create export state
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)                 :: exportState
    character(len=*)   , intent(in)  :: flds_scalar_name
    type (domain_type), intent(in), pointer :: domain
    integer            , intent(out) :: errorCode  ! pop error code
    integer            , intent(out) :: rc         ! returned error code

    ! local variables
    integer              :: n,i,j,k,iblock,nfld,lev
    character (cl) :: label
    real (r8)            :: m2percm2
    real (r8)            :: gsum
    integer              :: fieldCount
    character (cl), allocatable :: fieldNameList(:)
    character(len=*), parameter :: subname='(ocn_import_export:ocn_export)'

    ! local variables - mpas names
    integer :: iCell, nCells
    integer, dimension(:), pointer :: nCellsArray
    real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepth,     &
                          seaIceEnergy, accumulatedFrazilIceMass, frazilSurfacePressure
    integer, pointer :: index_avgZonalSSHGradient, index_avgMeridionalSSHGradient
    real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient,avgSurfaceVelocity
    real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, layerThickness
    integer, pointer :: index_temperatureSurfaceValue, index_salinitySurfaceValue, &
                        index_avgZonalSurfaceVelocity, index_avgMeridionalSurfaceVelocity

    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (block_type), pointer :: block
    integer gcell, cell_offset

    ! local variables - coupler names
    real (r8), pointer   :: So_t(:)
    real (r8), pointer   :: So_s(:)
    real (r8), pointer   :: So_u(:)
    real (r8), pointer   :: So_v(:)
    real (r8), pointer   :: So_bldepth(:)
    real (r8), pointer   :: So_dhdx(:)
    real (r8), pointer   :: So_dhdy(:)
    real (r8), pointer   :: So_omask(:)
    real (r8), pointer   :: Fioo_q(:)
    real (r8), pointer   :: Fioo_frazil(:)
    logical, pointer :: frazilIceActive
    logical :: keepFrazil
    real (kind=RKIND) :: surfaceFreezingTemp

    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !-----------------------------------------------------------------------
    ! ocean mask
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_omask', So_omask, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_omask(:) = 1.0  ! because grid exists only on ocean points

    !-----------------------------------------------------------------------
    ! surface velocities
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_u', So_u, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_u(:) = shr_const_spval

    call state_getfldptr(exportState, 'So_v', So_v, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_v(:) = shr_const_spval

    !-----------------------------------------------------------------------
    ! surface temperature
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_t', So_t, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_t(:) = shr_const_spval
    n = 0

    !-----------------------------------------------------------------------
    ! pack salinity
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_s', So_s, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_s(:) = shr_const_spval

    !-----------------------------------------------------------------------
    !  boundary layer depth
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_bldepth', So_bldepth, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_bldepth(:) = shr_const_spval

    !-----------------------------------------------------------------------
    !  Frazil mass flux
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'Fioo_frazil', Fioo_frazil, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    Fioo_frazil(:) = shr_const_spval

    !-----------------------------------------------------------------------
    !  freezing_melting_potential
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'Fioo_q', Fioo_q, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    Fioo_q(:) = shr_const_spval

    !-----------------------------------------------------------------------
    ! ssh gradients
    !-----------------------------------------------------------------------

    call state_getfldptr(exportState, 'So_dhdx', So_dhdx, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_dhdx(:) = shr_const_spval
    call state_getfldptr(exportState, 'So_dhdy', So_dhdy, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    So_dhdy(:) = shr_const_spval

    !-----------------------------------------------------------------------
    ! copy into mediator variables
    !-----------------------------------------------------------------------
    call mpas_pool_get_package(domain % packages, 'frazilIceActive', frazilIceActive)

    Fioo_q     (:) = 0.0_RKIND
    Fioo_frazil(:) = 0.0_RKIND
    
    cell_offset = 0
    block => domain % blocklist
    do while (associated(block))
       call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
       call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       nCells = nCellsArray( 1 )

       ! from forcing
       call mpas_pool_get_dimension(forcingPool, 'index_avgTemperatureSurfaceValue', index_temperatureSurfaceValue)
       call mpas_pool_get_dimension(forcingPool, 'index_avgSalinitySurfaceValue', index_salinitySurfaceValue)
       call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
       call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityZonal', index_avgZonalSurfaceVelocity)
       call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityMeridional', index_avgMeridionalSurfaceVelocity)
       call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
       call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
       call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)
       call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
       if ( frazilIceActive ) then
          call mpas_pool_get_subpool(block % structs, 'state', statePool)
          call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
          call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
          call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMass, 1)
          call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
       end if

       ! from forcing
       call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)

       do iCell = 1, nCells
          gcell = iCell + cell_offset
          So_t      (gcell) = avgTracersSurfaceValue(index_temperatureSurfaceValue,icell)
          So_s      (gcell) = avgTracersSurfaceValue(index_salinitySurfaceValue,icell)
          So_u      (gcell) = avgSurfaceVelocity(index_avgZonalSurfaceVelocity,icell)
          So_v      (gcell) = avgSurfaceVelocity(index_avgMeridionalSurfaceVelocity,icell)
          So_dhdx   (gcell) = avgSSHGradient(index_avgZonalSSHGradient,iCell)
          So_dhdy   (gcell) = avgSSHGradient(index_avgMeridionalSSHGradient,iCell)
          So_bldepth(gcell) = boundaryLayerDepth(iCell)
       end do

       if ( frazilIceActive ) then
          ! negative when frazil ice can be melted

          do iCell = 1, nCells
             gcell = iCell + cell_offset

             if ( accumulatedFrazilIceMass(iCell) > 0.0_RKIND ) then

              seaIceEnergy(iCell) = accumulatedFrazilIceMass(iCell) * config_frazil_heat_of_fusion 

             ! Otherwise calculate the melt potential where avgTracersSurfaceValue represents only the
             ! top layer of the ocean
             else

              surfaceFreezingTemp = ocn_freezing_temperature(salinity=avgTracersSurfaceValue(index_salinitySurfaceValue, iCell), &
                 pressure=0.0_RKIND,  inLandIceCavity=.false.) 

              seaIceEnergy(iCell) = min(rho_sw*cp_sw*layerThickness(1, iCell)*( surfaceFreezingTemp + T0_Kelvin &
                              - avgTracersSurfaceValue(index_temperatureSurfaceValue, iCell) ), 0.0_RKIND )

             end if

             Fioo_q     (gcell) = seaIceEnergy(iCell) / ocn_cpl_dt
             Fioo_frazil(gcell) = accumulatedFrazilIceMass(iCell) / ocn_cpl_dt

             ! Reset SeaIce Energy and Accumulated Frazil Ice
             seaIceEnergy(iCell) = 0.0_RKIND
             accumulatedFrazilIceMass(iCell) = 0.0_RKIND
             frazilSurfacePressure(iCell) = 0.0_RKIND
          end do
       endif

       cell_offset = cell_offset + nCells
       
       block => block % next
    end do
 

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ocn_export

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)

    ! input/output variables
    integer             , intent(inout) :: num
    type(fld_list_type) , intent(inout) :: fldlist(:)
    character(len=*)    , intent(in)    :: stdname
    integer, optional   , intent(in)    :: ungridded_lbound
    integer, optional   , intent(in)    :: ungridded_ubound

    ! local variables
    character(len=*), parameter :: subname='(fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call shr_sys_abort(trim(subname)//": ERROR num > fldsMax "//trim(stdname))
    endif
    fldlist(num)%stdname = trim(stdname)

    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC, only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(ocn_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO)

             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                                         ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                                         ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                                         gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if ! if not scalar field

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(fldlist_realize:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc) ! num of scalar values
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get 1d pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)     :: State
    character(len=*)  , intent(in)     :: fldname
    real(r8), pointer , intent(inout)  :: fldptr(:)
    integer, optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ocn_import_export:State_GetFldPtr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_1d

  !===============================================================================
  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    ! ----------------------------------------------
    ! Get 2d pointer to a state field
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)  , intent(in)     :: State
    character(len=*)  , intent(in)     :: fldname
    real(r8), pointer , intent(inout)  :: fldptr(:,:)
    integer, optional , intent(out)    :: rc

    ! local variables
    type(ESMF_Field) :: lfield
    character(len=*),parameter :: subname='(ocn_import_export:State_GetFldPtr_2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine State_GetFldPtr_2d

  !===============================================================================
  logical function State_FldChk(State, fldname)
    ! ----------------------------------------------
    ! Determine if field is in state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State) , intent(in)  :: State
    character(len=*) , intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag) :: itemType
    ! ----------------------------------------------

    call ESMF_StateGet(State, trim(fldname), itemType)
    State_FldChk = (itemType /= ESMF_STATEITEM_NOTFOUND)

  end function State_FldChk

end module ocn_import_export
