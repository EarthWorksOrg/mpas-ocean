! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for the split explicit
!>  time integration scheme
!
!-----------------------------------------------------------------------

module ocn_time_integration_split

   use mpaso_derived_types
   use mpaso_pool_routines
   use mpaso_constants
   use mpaso_dmpar
   use mpaso_vector_reconstruction
   use mpaso_spline_interpolation
   use mpaso_timer
   use mpaso_threading
   use mpaso_timekeeping
   use mpaso_log

   use ocn_config
   use ocn_mesh
   use ocn_tendency
   use ocn_diagnostics_variables
   use ocn_diagnostics
   use ocn_gm
   use ocn_submesoscale_eddies

   use ocn_compute_nonhydrostatic_pressure
   use ocn_equation_of_state
   use ocn_vmix
   use ocn_vertical_advection
   use ocn_vertical_remap
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice
   use ocn_surface_land_ice_fluxes
   use ocn_transport_tests

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_time_integrator_split, &
             ocn_time_integration_split_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   character (len=*), parameter :: &
      subcycleGroupName = 'subcycleFields', &! name for cached halo exch
      finalBtrGroupName = 'finalBtrFields'   ! name for cached halo exch

   logical :: &
      unsplit           ! flag for unsplit_explicit time option

   integer, dimension(:), allocatable :: &
      numClinicIterations  ! number of baroclinic iterations for each
                           ! outer timestep iteration

   integer :: &
      neededHalos,     &! number of halo levels needed
      haloDecrement,   &! number of halos to invalidate each cycle
      numTSIterations, &! number of outer timestep iterations
      nBtrSubcycles     ! number of barotropic subcycles

   integer :: ssh_sal_on

   real (kind=RKIND) :: &
      useVelocityCorrection,&! mask for velocity correction
      splitFact, &! mask for terms in split or unsplit cases
      self_attraction_and_loading_beta

!***********************************************************************

   contains

!***********************************************************************
!
!  ocn_time_integration_split
!
!> \brief MPAS ocean split explicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine integrates the ocean model forward one master time
!>  step (dt) using a split-explicit method in which the fast
!>  barotropic mode is sub-cycled using a shorter explicit time step.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_split(domain, dt)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step (sec) to move forward

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (dm_info) :: dminfo

      type (block_type), pointer :: &
         block ! structure with subdomain data

      type (mpas_pool_type), pointer :: &
         statePool,         &! structure holding state variables
         tracersPool,       &! structure holding tracers
         meshPool,          &! structure holding mesh variables
         verticalMeshPool,  &! structure holding vertical mesh variables
         tendPool,          &! structure holding tendencies
         tracersTendPool,   &! structure holding tracer tendencies
         forcingPool,       &! structure holding forcing variables
         scratchPool,       &! structure holding temporary variables
         swForcingPool,     &! structure holding short-wave forcing vars
         tracersIdealAgeFieldsPool     ! structure holding idealAge-related fields

      logical :: &
         activeTracersOnly   ! only compute tendencies for active tracers

      logical, pointer :: &
         config_use_tracerGroup  ! flag for using each tracer group

      integer :: &
         iCell, iEdge, k, &! loop iterators for cell, edge, vert loops
         nCells, nEdges,  &! number of cells or edges (incl halos)
         i,j,             &! generic loop iterators
         cell1, cell2,    &! neighbor cell addresses
         eoe,             &! index for edge on edge
         err,             &! local error flag
         splitExplicitStep, &! loop index for outer timestep loop
         oldBtrSubcycleTime,&! time index for old barotropic value
         newBtrSubcycleTime,&! time index for new barotropic value
         uPerpTime,       &! time index to use for uPerp calculation
         BtrCorIter,      &! loop index for barotropic coriolis cycle
         stage1_tend_time,&! time index for stage 1 tendencies
         edgeHaloComputeCounter, &! halo counters to reduce
         cellHaloComputeCounter   ! halo updates during cycling

      real (kind=RKIND) :: &
         normalThicknessFluxSum, &! sum of thickness flux in column
         thicknessSum,     &! sum of thicknesses in column
         flux,             &! temp for computing flux for barotropic
         sshEdge,          &! sea surface height at edge
         CoriolisTerm,     &! temp for computing coriolis term (fuperp)
         normalVelocityCorrection, &! velocity correction
         temp,             &! temp for holding vars at new time
         temp_h,           &! temporary for phi
         temp_mask,        &! temporary mask (edgeMask)
         lat,              &! cell latitude
         sshCell1,         &! sea sfc height in neighboring cells
         sshCell2

      real (kind=RKIND), dimension(:,:), allocatable:: &
         uTemp

      real (kind=RKIND), dimension(:), allocatable :: &
         btrvel_temp, &
         bottomDepthEdge
!$acc declare create(uTemp,bottomDepthEdge,btrvel_temp)

      ! State Array Pointers
      real (kind=RKIND), dimension(:), pointer :: &
         sshSubcycleCur,              &! subcycl ssh    at current time
         sshSubcycleNew,              &! subcycl ssh    at new     time
         sshSubcycleCurWithTides,     &! subcycl ssh    at current time
         sshSubcycleNewWithTides,     &! subcycl ssh    at new     time
         normalBarotropicVelocitySubcycleCur,&! barotropic vel subcyc
         normalBarotropicVelocitySubcycleNew,&! at current, new times
         sshCur,                      &! sea sfc height at current time
         sshNew,                      &! sea sfc height at new     time
         normalBarotropicVelocityCur, &! barotropic vel at current time
         normalBarotropicVelocityNew   ! barotropic vel at new     time

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalBaroclinicVelocityCur, &! baroclinic vel at current time
         normalBaroclinicVelocityNew, &! baroclinic vel at new     time
         normalVelocityCur,           &! full velocity  at current time
         normalVelocityNew,           &! full velocity  at new     time
         layerThicknessCur,           &! layer thick    at current time
         layerThicknessNew,           &! layer thick    at new     time
         highFreqThicknessCur,        &! high frq thick at current time
         highFreqThicknessNew,        &! high frq thick at new     time
         lowFreqDivergenceCur,        &! low frq div    at current time
         lowFreqDivergenceNew,        &! low frq div    at new     time
         tracerGroupIdealAgeMask,     & ! mask for resetting surface ideal age
         vertVelocityNonhydroCur,         &! vertical velocity at current time
         vertVelocityNonhydroNew,         &! vertical velocity at new time
         nonhydrostaticPressureNew,      &! nonhydro pressure at new time
         nonhydrostaticPressureCur     ! nonhdyro pressure at current time

      type (field1DReal), pointer :: &
         effectiveDensityField         ! field pointer for halo update

      ! State/Tracer arrays, info
      real (kind=RKIND), dimension(:,:,:), pointer ::      &!
         tracersGroupCur,      &! old, new tracer arrays
         tracersGroupNew

      integer :: &
         startIndex, endIndex, &! start, end index for tracer groups
         indexSalinity          ! tracer index for salinity

      integer :: iTracer, numTracers, numActiveTracers

      type (mpas_pool_iterator_type) :: &
         groupItr               ! iterator for tracer groups

      character (len=StrKIND) :: &
         modifiedGroupName,    &! constructed tracer group names
         configName             ! constructed config names for tracers

      ! Tendency Array Pointers

      real (kind=RKIND), dimension(:), pointer :: &
         sshTend                ! sea-surface height tendency

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocityTend,    &! normal velocity tendency
         highFreqThicknessTend, &! thickness tendency for high-freq iter
         lowFreqDivergenceTend, &! thickness tendency for low freq
         layerThicknessTend,    &! main thickness tendency
         vertVelocityNonhydroTend    ! vertical velocity tendency

      real (kind=RKIND), dimension(:,:,:), pointer :: &
         tracersGroupTend,      &! tendencies for tracer groups
         activeTracersTend       ! active tracer tendencies

      ! Forcing pool
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersNew

      ! Remap variables
      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThicknessLagNew

      ! These are public variables used from the diagnostics module
      ! indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      ! indexSSHGradientZonal, indexSSHGradientMeridional
      ! barotropicForcing, barotropicThicknessFlux
      ! layerThickEdgeFlux, layerThickEdgeMean,
      ! normalTransportVelocity, normalGMBolusVelocity
      ! vertAleTransportTop, normalMLEvelocity
      ! velocityX, velocityY, velocityZ
      ! velocityZonal, velocityMeridional
      ! gradSSH, gradSSHX, gradSSHY, gradSSHZ
      ! gradSSHZonal, gradSSHMeridional
      ! surfaceVelocity, SSHGradient
      ! temperatureShortWaveTendency
      ! activeTracerHorizontalAdvectionTendency
      ! activeTracerVerticalAdvectionTendency
      ! activeTracerSurfaceFluxTendency
      ! activeTracerNonLocalTendency
      ! activeTracerHorMixTendency
      ! activeTracerHorizontalAdvectionEdgeFlux
      ! vertVelocityTop

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Start main time step timer
      call mpas_timer_start("se timestep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Prep variables before first iteration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se prep")

      !*** Retrieve model state, pools
      !*** No longer support sub-blocks so this retrieves the only block

      block => domain%blocklist
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block%structs, 'verticalMesh', &
                                                 verticalMeshPool)
      call mpas_pool_get_subpool(block%structs, 'state', &
                                                 statePool)
      call mpas_pool_get_subpool(block%structs, 'forcing', &
                                                 forcingPool)
      call mpas_pool_get_subpool(block%structs, 'shortwave', &
                                                 swForcingPool)
      call mpas_pool_get_subpool(block%structs, 'tend', &
                                                 tendPool)
      call mpas_pool_get_subpool(block%structs, 'scratch', &
                                                 scratchPool)
      call mpas_pool_get_subpool(statePool, 'tracers', &
                                             tracersPool)
      call mpas_pool_get_subpool(tendPool,  'tracersTend', &
                                             tracersTendPool)

      !*** Retrieve state variables at two time levels
      if(config_enable_nonhydrostatic_mode) then
         write(*,*) "config_enable_nonhydrostatic_mode NOT EXECUTING"
         call mpas_pool_get_array(statePool, 'vertVelocityNonhydro', vertVelocityNonhydroCur, 1)
         call mpas_pool_get_array(statePool, 'vertVelocityNonhydro', vertVelocityNonhydroNew, 2)
         call mpas_pool_get_array(statePool, 'nonhydrostaticPressure', nonhydrostaticPressureNew, 2)
         call mpas_pool_get_array(statePool, 'nonhydrostaticPressure', nonhydrostaticPressureCur, 1)
      end if

      call mpas_pool_get_array_gpu(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocityCur, 1)
      call mpas_pool_get_array_gpu(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocityCur, 1)
      call mpas_pool_get_array_gpu(statePool, 'normalVelocity', &
                                           normalVelocityCur, 1)

      call mpas_pool_get_array_gpu(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocityNew, 2)
      call mpas_pool_get_array_gpu(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocityNew, 2)
      call mpas_pool_get_array_gpu(statePool, 'normalVelocity', &
                                           normalVelocityNew, 2)

      call mpas_pool_get_array_gpu(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array_gpu(statePool, 'ssh', sshNew, 2)

      call mpas_pool_get_array_gpu(statePool, 'layerThickness', &
                                           layerThicknessCur, 1)
      call mpas_pool_get_array_gpu(statePool, 'layerThickness', &
                                           layerThicknessNew, 2)

      call mpas_pool_get_array_gpu(statePool, 'highFreqThickness', &
                                           highFreqThicknessCur, 1)
      call mpas_pool_get_array_gpu(statePool, 'highFreqThickness', &
                                           highFreqThicknessNew, 2)

      call mpas_pool_get_array_gpu(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergenceCur, 1)
      call mpas_pool_get_array_gpu(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergenceNew, 2)

      call mpas_pool_get_dimension_scalar(tracersPool, 'index_salinity', &
                                                        indexSalinity)

      !*** Retrieve tendency variables

      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                          highFreqThicknessTend)
      call mpas_pool_get_array_gpu(tendPool, 'normalVelocity', &
                                          normalVelocityTend)
      call mpas_pool_get_array_gpu(tendPool, 'ssh', &
                                          sshTend)
      call mpas_pool_get_array_gpu(tendPool, 'layerThickness', &
                                          layerThicknessTend)
      call mpas_pool_get_array_gpu(tendPool, 'normalVelocity', &
                                          normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                          highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', &
                                          lowFreqDivergenceTend)

      call mpas_pool_get_array_gpu(tracersPool, 'activeTracers', activeTracersNew, 2)

      allocate(bottomDepthEdge(nEdgesAll+1))

      if (config_transport_tests_flow_id > 0) then
        write(*,*) "config_transport_tests_flow_id NOT EXECUTING"
        ! This is a transport test. Write advection velocity from prescribed
        ! flow field.
        call ocn_transport_test_velocity(meshPool, daysSinceStartOfSim, &
          0.5*dt, normalVelocityCur)
      endif

      ! Initialize * variables that are used to compute baroclinic
      ! tendencies below.

      ! The baroclinic velocity needs be recomputed at the beginning
      ! of a timestep because the implicit vertical mixing is
      ! conducted on the total u.  We keep normalBarotropicVelocity
      ! from the previous timestep.
      ! Note that normalBaroclinicVelocity may now include a
      ! barotropic component, because the weights layerThickness
      ! have changed.  That is OK, because the barotropicForcing
      ! variable subtracts out the barotropic component from the
      ! baroclinic.


!$acc data present(minLevelCell,maxLevelCell,&
!$acc              normalBaroclinicVelocityCur,normalVelocityCur,normalBarotropicVelocityCur,&
!$acc              normalVelocityNew,normalBaroclinicVelocityNew,sshNew,sshCur,layerThicknessNew,&
!$acc              layerThicknessCur,highFreqThicknessNew,highFreqThicknessCur,lowFreqDivergenceNew,&
!$acc              lowFreqDivergenceCur)
      !$omp parallel
      !$omp do schedule(runtime) private(k)
!$acc parallel vector_length(1024)
!$acc loop gang vector collapse(2)
      do iEdge = 1,nEdgesAll
      do k = 1,nVertLevels

         normalBaroclinicVelocityCur(k,iEdge) = &
                                 normalVelocityCur(k,iEdge) - &
                         normalBarotropicVelocityCur(iEdge)

         normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)

         normalBaroclinicVelocityNew(k,iEdge) = &
         normalBaroclinicVelocityCur(k,iEdge)
      end do
      end do
!$acc end parallel
      !$omp end do

      !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker 
      do iCell = 1, nCellsAll
         sshNew(iCell) = sshCur(iCell)
!$acc loop vector
         do k = 1,nVertLevels
            layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
         end do
      end do
!$acc end parallel
      !$omp end do
      !$omp end parallel

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr))
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array_gpu(tracersPool, groupItr%memberName,&
                                     tracersGroupCur, 1)
            call mpas_pool_get_array_gpu(tracersPool, groupItr%memberName,&
                                     tracersGroupNew, 2)

            if ( associated(tracersGroupCur) .and. &
                 associated(tracersGroupNew) ) then

               !$omp parallel
               !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
               do iCell = 1, nCellsAll
!$acc loop vector
               do k = minLevelCell(iCell), maxLevelCell(iCell)
                  do iTracer = 1,size(tracersGroupCur,1)
                     tracersGroupNew(iTracer,k,iCell) = &
                     tracersGroupCur(iTracer,k,iCell)
                  end do
               end do
               end do
!$acc end parallel
               !$omp end do
               !$omp end parallel
            end if
         end if
      end do

      if (associated(highFreqThicknessNew)) then
         !$omp parallel
         !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang vector collapse(2)
         do iCell = 1, nCellsAll
         do k = 1,nVertLevels
            highFreqThicknessNew(k,iCell) = &
            highFreqThicknessCur(k,iCell)
         end do
         end do
!$acc end parallel
         !$omp end do
         !$omp end parallel
      end if

      if (associated(lowFreqDivergenceNew)) then
         !$omp parallel
         !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang vector collapse(2)
         do iCell = 1, nCellsAll
         do k = 1,nVertLevels
            lowFreqDivergenceNew(k,iCell) = &
            lowFreqDivergenceCur(k,iCell)
         end do
         end do
!$acc end parallel
         !$omp end do
         !$omp end parallel
      endif

      if(config_enable_nonhydrostatic_mode) then
         !$omp parallel
         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCellsAll
            do k = 2, maxLevelCell(iCell)
               vertVelocityNonhydroNew(k,iCell) = vertVelocityNonhydroCur(k,iCell)
            end do
         end do
         !$omp end do
         !$omp end parallel
      end if

!$acc end data

      call mpas_timer_stop("se prep")

      call mpas_pool_get_array_gpu(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array_gpu(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array_gpu(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

      if (landIcePressureOn) then 
!NV!     !!! WORKING !!!
         call mpas_pool_get_array_gpu(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array_gpu(forcingPool, 'landIceDraft', landIceDraft)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large outer timestep iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do splitExplicitStep = 1, numTSIterations

         if (config_disable_thick_all_tend .and. &
             config_disable_vel_all_tend .and. &
             config_disable_tr_all_tend) then
            exit ! don't compute in loop meant to update velocity,
                 ! thickness, and tracers
         end if

         call mpas_timer_start('se loop')

         stage1_tend_time = min(splitExplicitStep,2)

         call mpas_pool_get_array_gpu(statePool, 'normalVelocity', &
                           normalVelocityCur, stage1_tend_time)

         ! ---  update halos for diagnostic ocean boundary layer depth
         if (config_use_cvmix_kpp) then
!NV!        !!! WORKING !!!
            call mpas_timer_start("se halo diag obd")
            call mpas_dmpar_field_halo_exch(domain,'boundaryLayerDepth')
            call mpas_timer_stop("se halo diag obd")
         end if

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("se halo diag")

         call mpas_dmpar_field_halo_exch(domain, &
                                 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
!NV!        !!! WORKING !!!
            call mpas_dmpar_field_halo_exch(domain, 'divergence')
            call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("se halo diag")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with
         !           long timestep
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then
!NV!        !!! NOT EXECUTING !!!
            write(*,*) "config_use_freq_filtered_thickness NOT EXECUTING"

            call mpas_timer_start("se freq-filtered-thick computations")

            call ocn_tend_freq_filtered_thickness(tendPool, statePool, &
                                                  stage1_tend_time)

            call mpas_timer_stop("se freq-filtered-thick computations")

            call mpas_timer_start("se freq-filtered-thick halo update")

            call mpas_dmpar_field_halo_exch(domain,'tendHighFreqThickness')
            call mpas_dmpar_field_halo_exch(domain,'tendLowFreqDivergence')

            call mpas_timer_stop("se freq-filtered-thick halo update")

            !$omp parallel
            !$omp do schedule(runtime) private(k)
            do iCell = 1, nCellsAll
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! this is h^{hf}_{n+1}
               highFreqThicknessNew(k,iCell) = &
               highFreqThicknessCur(k,iCell) + dt* &
               highFreqThicknessTend(k,iCell)
            end do
            end do
            !$omp end do
            !$omp end parallel

         endif ! freq filtered thickness

         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("se bcl vel")
         call mpas_timer_start('se bcl vel tend')

!NV!         call mpas_pool_get_array_gpu(statePool, 'normalVelocity', &
!NV!                           normalVelocityCur, stage1_tend_time)
         call mpas_pool_get_array(statePool, 'vertVelocityNonhydro', &
                           vertVelocityNonhydroCur, stage1_tend_time)

         ! compute vertAleTransportTop.  Use u (rather than &
         ! normalTransportVelocity) for momentum advection.
         ! Use the most recent time level available.

         if (associated(highFreqThicknessNew)) then
            call ocn_vert_transport_velocity_top_gpu( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalVelocityCur, &
                 sshCur, dt, vertAleTransportTop, &
                 err, highFreqThicknessNew)
         else
            call ocn_vert_transport_velocity_top_gpu( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalVelocityCur, &
                 sshCur, dt, vertAleTransportTop, err)
         endif

         call ocn_tend_vel_gpu(domain, tendPool, statePool, forcingPool, &
                           stage1_tend_time, domain % dminfo, dt)

         call mpas_timer_stop('se bcl vel tend')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,numClinicIterations(splitExplicitStep)

            call mpas_timer_start('bcl iters on linear Coriolis')

            ! Put f*normalBaroclinicVelocity^{perp} in
            ! normalVelocityNew as a work variable
            call ocn_fuperp_gpu(statePool, meshPool, 2)

            !TODO: look at loop optimizations here
            ! Only need to loop over owned cells, since there is a halo
            ! exchange immediately after this computation.

            allocate(uTemp(nVertLevels,nEdgesOwned))

!$acc data present(cellsOnEdge,minLevelEdgeBot,maxLevelEdgeTop,dcEdge,nEdgesHalo,&
!$acc              normalBaroclinicVelocityCur,normalVelocityTend,normalVelocityNew,sshNew,&
!$acc              layerThickEdgeFlux,barotropicForcing,normalBaroclinicVelocityNew,uTemp,bottomDepthEdge)
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k, cell1, cell2, uTemp, &
            !$omp         normalThicknessFluxSum, thicknessSum)
!$acc parallel
!$acc loop gang worker private(cell1,cell2,normalThicknessFluxSum,thicknessSum)
            do iEdge = 1, nEdgesOwned
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
!$acc loop vector
               do k = 1, nVertLevels
                  uTemp(k,iEdge) = 0.0_RKIND
               enddo
!$acc loop vector
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  ! normalBaroclinicVelocityNew =
                  ! normalBaroclinicVelocityOld +
                  !                dt*(-f*normalBaroclinicVelocityPerp
                  !                    + T(u*,w*,p*) + g*grad(SSH*) )
                  ! Here uNew is a work variable containing
                  !  -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                  uTemp(k,iEdge) = normalBaroclinicVelocityCur(k,iEdge) &
                           + dt * (normalVelocityTend(k,iEdge) &
                           + normalVelocityNew(k,iEdge) &
                           + splitFact*gravity * &
                             (sshNew(cell2) - sshNew(cell1)) &
                             /dcEdge(iEdge) )
               enddo ! vertical

               ! thicknessSum is initialized outside the loop because
               ! on land boundaries maxLevelEdgeTop=0, but we want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               normalThicknessFluxSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)* &
                                        uTemp(minLevelEdgeBot(iEdge),iEdge)
               thicknessSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)
!$acc loop vector reduction(+:normalThicknessFluxSum,thicknessSum)
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                 layerThickEdgeFlux(k,iEdge)*uTemp(k,iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo
               barotropicForcing(iEdge) = splitFact* &
                              normalThicknessFluxSum/thicknessSum/dt
!$acc loop vector
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  ! These two steps are together here:
                  !{\bf u}'_{k,n+1} =
                  !    {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                  !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left(
                  !        {\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right)
                  ! so that normalBaroclinicVelocityNew is at time n+1/2
                  normalBaroclinicVelocityNew(k,iEdge) = 0.5_RKIND*( &
                  normalBaroclinicVelocityCur(k,iEdge) + uTemp(k,iEdge) - &
                         dt * barotropicForcing(iEdge))

               enddo
            enddo ! iEdge
!$acc end parallel
            !$omp end do
            !$omp end parallel

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(cell1, cell2, k, thicknessSum)
!$acc parallel
!$acc loop gang vector private(cell1,cell2,thicknessSum)
            do iEdge = 1, nEdgesHalo(config_num_halos+1)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               thicknessSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)
!$acc loop seq
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo
               bottomDepthEdge(iEdge) = thicknessSum &
                  - 0.5_RKIND*(sshNew(cell1) + sshNew(cell2))
            enddo ! iEdge
!$acc end parallel
            !$omp end do
            !$omp end parallel
!$acc end data

            deallocate(uTemp)

            call mpas_timer_start("se halo normalBaroclinicVelocity")
            call mpas_dmpar_field_halo_exch(domain, &
                              'normalBaroclinicVelocity', timeLevel=2)
            call mpas_timer_stop("se halo normalBaroclinicVelocity")

            call mpas_timer_stop('bcl iters on linear Coriolis')

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_start('se halo barotropicForcing')
         call mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         call mpas_timer_stop('se halo barotropicForcing')

         call mpas_timer_stop("se bcl vel")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 2: Barotropic velocity prediction, explicitly subcycled
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call mpas_timer_start("se btr vel")

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2

         ! unsplit_explicit option
         if (unsplit) then
!NV!        !!! NOT EXECUTING !!!
            write(*,*) "unsplit NOT EXECUTING"

            call mpas_timer_start('btr vel ue')

            ! For Split_Explicit unsplit, simply set
            !    normalBarotropicVelocityNew=0,
            !    normalBarotropicVelocitySubcycle=0, and
            !    uNew=normalBaroclinicVelocityNew

            !$omp parallel
            !$omp do schedule(runtime) private(k)
            do iEdge = 1, nEdgesAll
               normalBarotropicVelocityNew(iEdge) = 0.0_RKIND
               do k = 1, nVertLevels
                  normalVelocityNew(k, iEdge) = &
                     normalBaroclinicVelocityNew(k, iEdge)

                  ! This is u used in advective terms for layerThickness and tracers
                  ! in tendency calls in stage 3.
                  normalTransportVelocity(k,iEdge) = normalBaroclinicVelocityNew(k,iEdge)


               enddo ! vertical
            end do  ! iEdge
            !$omp end do
            !$omp end parallel

            !add GM and MLE (submesoscale) velocities if requested
            call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdgesAll, &
                                             nVertLevels)
            call ocn_MLE_add_to_transport_vel(normalTransportvelocity, nEdgesAll)

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k)
            do iEdge=1, nEdgesAll
               do k = 1, nVertLevels
                  normalTransportVelocity(k,iEdge) = edgeMask(k,iEdge)* &
                                  normalTransportVelocity(k,iEdge)
               end do
            end do
            !$omp end do
            !$omp end parallel

            call mpas_timer_stop('btr vel ue')

         else ! split explicit option

            ! Initialize variables for barotropic subcycling
            call mpas_timer_start('btr vel se init')

            call mpas_pool_get_array_gpu(statePool, 'sshSubcycle', &
                                     sshSubcycleCur, oldBtrSubcycleTime)
            call mpas_pool_get_array_gpu(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleCur, &
                                  oldBtrSubcycleTime)

            if (config_filter_btr_mode) then
!NV!           !!! NOT EXECUTING !!!
               write(*,*) "config_filter_btr_mode NOT EXECUTING"
               !$omp parallel
               !$omp do schedule(runtime)
               do iEdge = 1, nEdgesAll
                  barotropicForcing(iEdge) = 0.0_RKIND
               end do
               !$omp end do
               !$omp end parallel
            endif

            !$omp parallel
            !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector 
            do iCell = 1, nCellsAll
               ! sshSubcycleOld = sshOld
               sshSubcycleCur(iCell) = sshCur(iCell)
            end do
!$acc end parallel
            !$omp end do

            !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector
            do iEdge = 1, nEdgesAll

               ! normalBarotropicVelocitySubcycleOld =
               ! normalBarotropicVelocityOld
               normalBarotropicVelocitySubcycleCur(iEdge) = &
               normalBarotropicVelocityCur(iEdge)

               ! normalBarotropicVelocityNew = BtrOld
               !  This is the first for the summation
               normalBarotropicVelocityNew(iEdge) = &
               normalBarotropicVelocityCur(iEdge)

               ! barotropicThicknessFlux = 0
               barotropicThicknessFlux(iEdge) = 0.0_RKIND
            end do
!$acc end parallel
            !$omp end do
            !$omp end parallel

            call mpas_timer_stop('btr vel se init')

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BEGIN Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Allocate subcycled scratch fields before starting
            ! subcycle loop
            allocate(btrvel_temp(nEdgesAll+1))
!$acc parallel
!$acc loop gang vector
            do j = 1, nEdgesAll+1
               btrvel_temp(j) = 0.0_RKIND
            enddo
!$acc end parallel

            cellHaloComputeCounter = 0
            edgeHaloComputeCounter = 0

            call mpas_timer_start('btr se subcycle loop')
            do j = 1, nBtrSubcycles * config_btr_subcycle_loop_factor

               ! Update halos if needed
               if (cellHaloComputeCounter < neededHalos) then

                  call mpas_timer_start('se halo subcycle')
                  call mpas_dmpar_exch_group_reuse_halo_exch(domain, &
                       subcycleGroupName, timeLevel=oldBtrSubcycleTime)
                  call mpas_timer_stop('se halo subcycle')

                  ! Reset halo counters after halo update
                  cellHaloComputeCounter = config_num_halos - &
                                           haloDecrement
                  edgeHaloComputeCounter = config_num_halos + 1 - &
                                           haloDecrement
               end if ! halo update

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! only do this part if it is needed in next SSH solve
               if (config_btr_gam1_velWt1 > 1.0e-12_RKIND) then
                  uPerpTime = oldBtrSubcycleTime

                  call mpas_pool_get_array_gpu(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleCur, &
                                  uPerpTime)
                  call mpas_pool_get_array_gpu(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleNew, &
                                  newBtrSubcycleTime)
                  call mpas_pool_get_array_gpu(statePool, 'sshSubcycle', &
                                  sshSubcycleCur, oldBtrSubcycleTime)

                  ! TODO: Modify this to eliminate pointer reassigning
                  ! Subtract tidal potential from ssh, if needed
                  ! Subtract the tidal potential from the current
                  !    subcycle ssh and store and a work array.
                  ! Then point sshSubcycleCur to the work array so the
                  ! tidal potential terms are included in the grad
                  ! operator inside the edge loop.
                  if (config_use_tidal_potential_forcing) then
!NV!                 !!! NOT EXECUTING !!!
                     write(*,*) "config_use_tidal_potential_forcing NOT EXECUTING"
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleCurWithTides', &
                                            sshSubcycleCurWithTides)
                     call mpas_pool_get_array(forcingPool, &
                                             'tidalPotentialEta', &
                                              tidalPotentialEta)

                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        sshSubcycleCurWithTides(iCell) = &
                                 sshSubcycleCur(iCell) - &
                              tidalPotentialEta(iCell) - &
                           ssh_sal_on * ssh_sal(iCell) - &
                           (1 - ssh_sal_on) * self_attraction_and_loading_beta* &
                                 sshSubcycleCur(iCell)
                     end do
                     !$omp end do
                     !$omp end parallel

                     call mpas_pool_get_array_gpu(forcingPool, &
                                          'sshSubcycleCurWithTides', &
                                           sshSubcycleCur)

                  end if ! tidal potential forcing

                  if (edgeHaloComputeCounter <= 1) then
                     nEdges = nEdgesOwned
                  else
                     nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
                  endif

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(i, cell1, cell2, eoe, CoriolisTerm, &
                  !$omp         temp_mask)
!$acc parallel
!$acc loop gang vector private(temp_mask,cell1,cell2,CoriolisTerm,eoe)
                  do iEdge = 1, nEdges

                     temp_mask = edgeMask(1, iEdge)

                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0_RKIND
!$acc loop seq
                     do i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm + &
                                       weightsOnEdge(i,iEdge)* &
                             normalBarotropicVelocitySubcycleCur(eoe)* &
                                       fEdge(eoe)
                     end do

                     normalBarotropicVelocitySubcycleNew(iEdge) = &
                           temp_mask* &
                           (normalBarotropicVelocitySubcycleCur(iEdge) &
                            + dt/nBtrSubcycles*(CoriolisTerm - gravity &
                             *(sshSubcycleCur(cell2) - &
                               sshSubcycleCur(cell1))/dcEdge(iEdge) + &
                               barotropicForcing(iEdge)))

                  end do ! edges
!$acc end parallel
                  !$omp end do
                  !$omp end parallel

               endif ! config_btr_gam1_velWt1>1.0e-12

               ! Halo from edges is corrupted here, so reduce the edge halo
               ! counter by 1
               edgeHaloComputeCounter = edgeHaloComputeCounter - 1

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: SSH PREDICTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               call mpas_pool_get_array_gpu(statePool, 'sshSubcycle', &
                                                    sshSubcycleCur, &
                                                    oldBtrSubcycleTime)
               call mpas_pool_get_array_gpu(statePool, 'sshSubcycle', &
                                                    sshSubcycleNew, &
                                                    newBtrSubcycleTime)
               call mpas_pool_get_array_gpu(statePool, &
                             'normalBarotropicVelocitySubcycle', &
                              normalBarotropicVelocitySubcycleCur, &
                              oldBtrSubcycleTime)
               call mpas_pool_get_array_gpu(statePool, &
                             'normalBarotropicVelocitySubcycle', &
                              normalBarotropicVelocitySubcycleNew, &
                              newBtrSubcycleTime)

               if (cellHaloComputeCounter <= 1) then
                  nCells = nCellsOwned
               else
                  nCells = nCellsHalo(cellHaloComputeCounter-1)
               endif
               if (edgeHaloComputeCounter <= 1) then
                  nEdges = nEdgesOwned
               else
                  nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
               endif

               ! config_btr_gam1_velWt1 sets the forward weighting of
               !    velocity in the SSH computation
               ! config_btr_gam1_velWt1=  1
               !    flux = normalBarotropicVelocityNew*H
               ! config_btr_gam1_velWt1=0.5
               !    flux = 1/2*(normalBarotropicVelocityNew +
               !                normalBarotropicVelocityOld)*H
               ! config_btr_gam1_velWt1=  0
               !    flux = normalBarotropicVelocityOld*H

               !$omp parallel
               !$omp do schedule(runtime) &
               !$omp private(i, iEdge, cell1, cell2, sshEdge, &
               !$omp         thicknessSum, flux)
!$acc parallel
!$acc loop gang vector private(iEdge,cell1,cell2,sshEdge,thicknessSum,flux)
               do iCell = 1, nCells
                  sshTend(iCell) = 0.0_RKIND
!$acc loop seq
                  do i = 1, nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(i, iCell)
!NV!                     if (maxLevelEdgeTop(iEdge).eq.0) cycle
                     if (maxLevelEdgeTop(iEdge).ne.0) then

                     cell1 = cellsOnEdge(1, iEdge)
                     cell2 = cellsOnEdge(2, iEdge)

                     sshEdge = 0.5_RKIND*(sshSubcycleCur(cell1) + &
                                          sshSubcycleCur(cell2))

                     ! Compute barotropic thickness at the edge. Note this
                     ! matches the sum of baroclinic thicknesses at this edge.
                     thicknessSum = sshEdge + bottomDepthEdge(iEdge)

                     flux = ((1.0-config_btr_gam1_velWt1)* &
                            normalBarotropicVelocitySubcycleCur(iEdge) &
                          +       config_btr_gam1_velWt1 * &
                            normalBarotropicVelocitySubcycleNew(iEdge))&
                           *thicknessSum

                     sshTend(iCell) = sshTend(iCell) + &
                                      edgeSignOncell(i, iCell)*flux* &
                                      dvEdge(iEdge)
                     endif
                  end do ! edges on cell

                  ! SSHnew = SSHold + dt/J*(-div(Flux))
                  sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                                        + dt/nBtrSubcycles* &
                                          sshTend(iCell)/areaCell(iCell)
               end do ! cells
!$acc end parallel
               !$omp end do
               !$omp end parallel

               ! avoid redundant computations when
               ! config_btr_solve_SSH2 is true
               if (config_btr_solve_SSH2) then
!NV!              !!! NOT EXECUTING !!!
                  write(*,*) "config_btr_solve_SSH2 NOT EXECUTING"
                  ! If config_btr_solve_SSH2=.true.,
                  ! then do NOT accumulate barotropicThicknessFlux
                  ! in this SSH predictor section, because it will be
                  ! accumulated in the SSH corrector section.

               else

                  ! otherwise, DO accumulate barotropicThicknessFlux
                  ! in this SSH predictor section
                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(cell1, cell2, sshEdge,thicknessSum,flux)
!$acc parallel
!$acc loop gang vector private(cell1,cell2,sshEdge,thicknessSum,flux)
                  do iEdge = 1, nEdges
!NV!                     if (maxLevelEdgeTop(iEdge).eq.0) cycle
                     if (maxLevelEdgeTop(iEdge).ne.0) then
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     sshEdge = 0.5_RKIND*(sshSubcycleCur(cell1) + &
                                          sshSubcycleCur(cell2))

                     thicknessSum = sshEdge + bottomDepthEdge(iEdge)

                     flux = ((1.0-config_btr_gam1_velWt1)* &
                            normalBarotropicVelocitySubcycleCur(iEdge) &
                            +     config_btr_gam1_velWt1 * &
                            normalBarotropicVelocitySubcycleNew(iEdge))&
                            *thicknessSum

                     barotropicThicknessFlux(iEdge) = &
                     barotropicThicknessFlux(iEdge) + flux
                     endif
                  end do ! edges
!$acc end parallel
                  !$omp end do
                  !$omp end parallel

               endif

               ! a cell halo layer is now corrupted - decrement counter
               cellHaloComputeCounter = cellHaloComputeCounter - 1

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: VELOCITY CORRECTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               do BtrCorIter = 1, config_n_btr_cor_iter
                  uPerpTime = newBtrSubcycleTime

                  call mpas_pool_get_array_gpu(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleCur, &
                                  oldBtrSubcycleTime)
                  call mpas_pool_get_array_gpu(statePool, &
                                 'normalBarotropicVelocitySubcycle', &
                                  normalBarotropicVelocitySubcycleNew, &
                                  newBtrSubcycleTime)
                  call mpas_pool_get_array_gpu(statePool, 'sshSubcycle', &
                                                       sshSubcycleCur,&
                                                    oldBtrSubcycleTime)
                  call mpas_pool_get_array_gpu(statePool, 'sshSubcycle', &
                                                       sshSubcycleNew,&
                                                    newBtrSubcycleTime)

                  ! TODO: change to avoid pointer reassignment
                  ! Subtract tidal potential from ssh, if needed
                  ! Subtract the tidal potential from the current and
                  !    new subcycle ssh and store in work arrays.
                  ! Then point sshSubcycleCur and ssh SubcycleNew to
                  ! the work arrays so the tidal potential terms are
                  ! included in the grad operator inside the edge loop.
                  if (config_use_tidal_potential_forcing) then
!NV!                 !!! NOT EXECUTING !!!
                     write(*,*) "config_use_tidal_potential_forcing NOT EXECUTING"
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleCurWithTides',  &
                                            sshSubcycleCurWithTides)
                     call mpas_pool_get_array(forcingPool, &
                                           'sshSubcycleNewWithTides',  &
                                            sshSubcycleNewWithTides)
                     call mpas_pool_get_array(forcingPool,  &
                                             'tidalPotentialEta',  &
                                              tidalPotentialEta)

                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        sshSubcycleCurWithTides(iCell) = &
                                 sshSubcycleCur(iCell) - &
                              tidalPotentialEta(iCell) - &
                           ssh_sal_on * ssh_sal(iCell) - &
                           (1 - ssh_sal_on) * self_attraction_and_loading_beta* &
                                 sshSubcycleCur(iCell)

                        sshSubcycleNewWithTides(iCell) = &
                                 sshSubcycleNew(iCell) - &
                              tidalPotentialEta(iCell) - &
                           ssh_sal_on * ssh_sal(iCell) - &
                           (1 - ssh_sal_on) * self_attraction_and_loading_beta* &
                                 sshSubcycleNew(iCell)
                     end do
                     !$omp end do
                     !$omp end parallel

                     call mpas_pool_get_array_gpu(forcingPool, &
                                    'sshSubcycleCurWithTides',  &
                                     sshSubcycleCur)
                     call mpas_pool_get_array(forcingPool, &
                                    'sshSubcycleNewWithTides',  &
                                     sshSubcycleNew)
                  end if ! tidal potential forcing

                  ! Need to initialize btr_vel_temp over one more halo
                  ! than we are computing over
                  nEdges = nEdgesHalo(min(edgeHaloComputeCounter, &
                                          config_num_halos))

                  !$omp parallel
                  !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector
                  do iEdge = 1, nEdges+1
                     btrvel_temp(iEdge) = &
                         normalBarotropicVelocitySubcycleNew(iEdge)
                  end do
!$acc end parallel
                  !$omp end do
                  !$omp end parallel

                  if (edgeHaloComputeCounter <= 1) then
                     nEdges = nEdgesOwned
                  else
                     nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
                  endif

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(i, cell1, cell2, eoe, coriolisTerm, &
                  !$omp         sshCell1, sshCell2, temp_mask)
!$acc parallel
!$acc loop gang vector private(temp_mask,cell1,cell2,CoriolisTerm,eoe,sshCell1,sshCell2)
                  do iEdge = 1, nEdges

                     temp_mask = edgeMask(1,iEdge)

                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! Compute the barotropic Coriolis term, -f*uPerp
                     CoriolisTerm = 0.0_RKIND
!$acc loop seq
                     do i = 1, nEdgesOnEdge(iEdge)
                        eoe = edgesOnEdge(i,iEdge)
                        CoriolisTerm = CoriolisTerm &
                                     + weightsOnEdge(i,iEdge)* &
                                       btrvel_temp(eoe)*fEdge(eoe)
                     end do

                     ! In this final solve for velocity, SSH is a
                     ! linear combination of SSHold and SSHnew.
                     sshCell1 = (1 - config_btr_gam2_SSHWt1) * &
                                         sshSubcycleCur(cell1) &
                              +      config_btr_gam2_SSHWt1  * &
                                         sshSubcycleNew(cell1)
                     sshCell2 = (1 - config_btr_gam2_SSHWt1) * &
                                         sshSubcycleCur(cell2) &
                              +      config_btr_gam2_SSHWt1  * &
                                         sshSubcycleNew(cell2)

                     ! normalBarotropicVelocityNew =
                     ! normalBarotropicVelocityOld + dt/J*
                     !     (-f*normalBarotropicVelocityoldPerp
                     !      - g*grad(SSH) + G)
                     normalBarotropicVelocitySubcycleNew(iEdge) = &
                          temp_mask * &
                          (normalBarotropicVelocitySubcycleCur(iEdge) &
                           + dt/nBtrSubcycles &
                           *(CoriolisTerm - gravity* &
                             (sshCell2 - sshCell1)/dcEdge(iEdge) &
                           + barotropicForcing(iEdge)))

                  end do
!$acc end parallel
                  !$omp end do
                  !$omp end parallel

                  edgeHaloComputeCounter = edgeHaloComputeCounter - 1
                  if (BtrCorIter >= 1 .or. &
                      .not. config_btr_solve_SSH2) &
                     cellHaloComputeCounter = cellHaloComputeCounter-1

               end do !do BtrCorIter=1,config_n_btr_cor_iter

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: SSH CORRECTOR STEP
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               if (config_btr_solve_SSH2) then
!NV!              !!! NOT EXECUTING !!!
                  write(*,*) "config_btr_solve_SSH2 NOT EXECUTING"
                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                       sshSubcycleCur,&
                                                    oldBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, 'sshSubcycle', &
                                                       sshSubcycleNew,&
                                                    newBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, &
                                  'normalBarotropicVelocitySubcycle', &
                                   normalBarotropicVelocitySubcycleCur, &
                                   oldBtrSubcycleTime)
                  call mpas_pool_get_array(statePool, &
                                  'normalBarotropicVelocitySubcycle', &
                                   normalBarotropicVelocitySubcycleNew, &
                                   newBtrSubcycleTime)


                  if (cellHaloComputeCounter <= 1) then
                     nCells = nCellsOwned
                  else
                     nCells = nCellsHalo(cellHaloComputeCounter-1)
                  endif
                  if (edgeHaloComputeCounter <= 1) then
                     nEdges = nEdgesOwned
                  else
                     nEdges = nEdgesHalo(edgeHaloComputeCounter-1)
                  endif

                  ! config_btr_gam3_velWt2 sets the forward weighting
                  ! of velocity in the SSH computation
                  ! config_btr_gam3_velWt2=  1
                  !    flux = normalBarotropicVelocityNew*H
                  ! config_btr_gam3_velWt2=0.5
                  !    flux = 1/2*(normalBarotropicVelocityNew +
                  !                normalBarotropicVelocityOld)*H
                  ! config_btr_gam3_velWt2=  0
                  !    flux = normalBarotropicVelocityOld*H

                  !$omp parallel
                  !$omp do schedule(runtime) &
                  !$omp private(i, iEdge, cell1, cell2, thicknessSum, &
                  !$omp         sshCell1, sshCell2, sshEdge, flux)
                  do iCell = 1, nCells
                     sshTend(iCell) = 0.0_RKIND
                     do i = 1, nEdgesOnCell(iCell)
                        iEdge = edgesOnCell(i, iCell)
                        if (maxLevelEdgeTop(iEdge).eq.0) cycle

                        cell1 = cellsOnEdge(1,iEdge)
                        cell2 = cellsOnEdge(2,iEdge)

                        ! SSH is a linear combination of SSHold
                        ! and SSHnew.
                        sshCell1 = (1-config_btr_gam2_SSHWt1)* &
                                         sshSubcycleCur(cell1) &
                                 +    config_btr_gam2_SSHWt1 * &
                                         sshSubcycleNew(cell1)
                        sshCell2 = (1-config_btr_gam2_SSHWt1)* &
                                         sshSubcycleCur(cell2) &
                                 +    config_btr_gam2_SSHWt1 * &
                                         sshSubcycleNew(cell2)

                        sshEdge = 0.5_RKIND*(sshCell1 + sshCell2)

                        thicknessSum = sshEdge + bottomDepthEdge(iEdge)

                        flux = ((1.0-config_btr_gam3_velWt2)* &
                           normalBarotropicVelocitySubcycleCur(iEdge) &
                             +       config_btr_gam3_velWt2 * &
                           normalBarotropicVelocitySubcycleNew(iEdge))&
                             *thicknessSum

                        sshTend(iCell) = sshTend(iCell) + &
                                         edgeSignOnCell(i, iCell)* &
                                         flux*dvEdge(iEdge)

                     end do ! edges on cell

                     ! SSHnew = SSHold + dt/J*(-div(Flux))
                     sshSubcycleNew(iCell) = sshSubcycleCur(iCell) &
                                        + dt/nBtrSubcycles * &
                                          sshTend(iCell)/areaCell(iCell)
                  end do ! cell loop for ssh
                  !$omp end do

                  ! compute barotropic thickness flux on edges
                  !$omp do schedule(runtime) &
                  !$omp private(cell1, cell2, thicknessSum, &
                  !$omp         sshCell1, sshCell2, sshEdge, flux)
                  do iEdge = 1, nEdges
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)

                     ! SSH is linear combination of SSHold and SSHnew
                     sshCell1 = (1-config_btr_gam2_SSHWt1)* &
                                      sshSubcycleCur(cell1) &
                              +    config_btr_gam2_SSHWt1 * &
                                      sshSubcycleNew(cell1)
                     sshCell2 = (1-config_btr_gam2_SSHWt1)* &
                                      sshSubcycleCur(cell2) &
                              +    config_btr_gam2_SSHWt1 * &
                                      sshSubcycleNew(cell2)
                     sshEdge = 0.5_RKIND * (sshCell1 + sshCell2)

                     thicknessSum = sshEdge + bottomDepthEdge(iEdge)

                     flux = ((1.0-config_btr_gam3_velWt2) * &
                         normalBarotropicVelocitySubcycleCur(iEdge) &
                          +       config_btr_gam3_velWt2  * &
                         normalBarotropicVelocitySubcycleNew(iEdge)) &
                            * thicknessSum

                     barotropicThicknessFlux(iEdge) = &
                     barotropicThicknessFlux(iEdge) + flux

                  end do ! edge barotropic thickness flux
                  !$omp end do
                  !$omp end parallel

                  edgeHaloComputeCounter = config_num_halos + 1
               endif ! config_btr_solve_SSH2

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! Barotropic subcycle: Accumulate running sums, advance
               !                      timestep pointers
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               call mpas_pool_get_array_gpu(statePool, &
                                  'normalBarotropicVelocitySubcycle', &
                                   normalBarotropicVelocitySubcycleNew, &
                                   newBtrSubcycleTime)

               ! normalBarotropicVelocityNew =
               ! normalBarotropicVelocityNew +
               ! normalBarotropicVelocitySubcycleNEW
               ! This accumulates the sum. If the Barotropic Coriolis
               ! iteration is limited to one, this could
               ! be merged with the above code.

               !$omp parallel
               !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector 
               do iEdge = 1, nEdgesAll
                  normalBarotropicVelocityNew(iEdge) = &
                  normalBarotropicVelocityNew(iEdge) + &
                  normalBarotropicVelocitySubcycleNew(iEdge)
               end do  ! iEdge
!$acc end parallel
               !$omp end do
               !$omp end parallel

               ! advance time pointers
               oldBtrSubcycleTime = mod(oldBtrSubcycleTime,2)+1
               newBtrSubcycleTime = mod(newBtrSubcycleTime,2)+1

            end do ! j=1,nBtrSubcycles
            call mpas_timer_stop('btr se subcycle loop')

            deallocate(btrvel_temp)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END Barotropic subcycle loop
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Normalize Barotropic subcycle sums: ssh,
            ! normalBarotropicVelocity, and F
            call mpas_timer_start('btr se norm')

            !$omp parallel
            !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector
            do iEdge = 1, nEdgesOwned
               barotropicThicknessFlux(iEdge) = &
               barotropicThicknessFlux(iEdge) &
                      /(nBtrSubcycles*config_btr_subcycle_loop_factor)

               normalBarotropicVelocityNew(iEdge) = &
               normalBarotropicVelocityNew(iEdge) &
                     / (nBtrSubcycles*config_btr_subcycle_loop_factor + 1)
            end do
!$acc end parallel
            !$omp end do
            !$omp end parallel

            call mpas_timer_stop('btr se norm')

            ! boundary update on F
            call mpas_timer_start("se halo F and btr vel")
            call mpas_dmpar_exch_group_create(domain, finalBtrGroupName)

            call mpas_dmpar_exch_group_add_field(domain, &
                                              finalBtrGroupName, &
                                             'barotropicThicknessFlux')
            call mpas_dmpar_exch_group_add_field(domain, &
                                              finalBtrGroupName, &
                                             'normalBarotropicVelocity', &
                                              timeLevel=2)

            call mpas_dmpar_exch_group_full_halo_exch(domain, &
                                                      finalBtrGroupName)

            call mpas_dmpar_exch_group_destroy(domain,finalBtrGroupName)
            call mpas_timer_stop("se halo F and btr vel")

            ! Check that you can compute SSH using the total sum or the
            ! individual increments over the barotropic subcycles.
            ! efficiency: This next block of code is really a check for
            ! debugging, and can be removed later.
            call mpas_timer_start('btr se ssh verif')

            ! Correction velocity
            !   normalVelocityCorrection = (Flux - Sum(h u*))/H
            ! or, for the full latex version:
            !{\bf u}^{corr} = \left( {\overline {\bf F}}
            !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}
            ! {\bf u}_k^{avg} \right)
            ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right.

            nEdges = nEdgesHalo(config_num_halos-1 )
            allocate(uTemp(nVertLevels,nEdges))

            !Compute uTemp
            ! velocity for normalVelocityCorrection is
            ! normalBarotropicVelocity +
            ! normalBaroclinicVelocity + uBolus
            !$omp parallel
            !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector collapse(2)
            do iEdge = 1, nEdges
               do k = 1, nVertLevels
                  uTemp(k,iEdge) = normalBarotropicVelocityNew(iEdge) &
                               + normalBaroclinicVelocityNew(k,iEdge)
               enddo
            end do
!$acc end parallel
            !$omp end do
            !$omp end parallel

            call ocn_GM_add_to_transport_vel_gpu(uTemp, nEdges, &
                                             nVertLevels)
!NV!        !!! NOT EXECUTING !!!
            call ocn_MLE_add_to_transport_vel_gpu(uTemp, nEdges)

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k)
!$acc parallel
!$acc loop gang vector collapse(2)
            do iEdge = 1, nEdges
               do k = 1, nVertLevels
                  normalTransportVelocity(k,iEdge) =  &
                             normalBarotropicVelocityNew(iEdge)   + &
                             normalBaroclinicVelocityNew(k,iEdge)
               end do
            end do
!$acc end parallel
            !$omp end do
            !$omp end parallel

            call ocn_GM_add_to_transport_vel_gpu(normalTransportVelocity, nEdges, &
                                             nVertLevels)
!NV!        !!! NOT EXECUTING !!!
            call ocn_MLE_add_to_transport_vel_gpu(normalTransportVelocity, nEdges)

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k, normalThicknessFluxSum, &
            !$omp         thicknessSum, normalVelocityCorrection)
!$acc parallel
!$acc loop gang worker private(normalThicknessFluxSum,thicknessSum,normalVelocityCorrection)
            do iEdge = 1, nEdges

               ! thicknessSum is initialized outside the loop because
               ! on land boundaries maxLevelEdgeTop=0, but I want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               normalThicknessFluxSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)* &
                                        uTemp(minLevelEdgeBot(iEdge),iEdge)
               thicknessSum  = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)

!$acc loop vector reduction(+:normalThicknessFluxSum,thicknessSum)
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                           layerThickEdgeFlux(k,iEdge)*&
                                           uTemp(k,iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo

               normalVelocityCorrection = useVelocityCorrection* &
                             ((barotropicThicknessFlux(iEdge) - &
                               normalThicknessFluxSum)/thicknessSum)

!$acc loop vector
               do k = 1, nVertLevels

                  ! normalTransportVelocity = normalBarotropicVelocity
                  !                         + normalBaroclinicVelocity
                  !                         + normalGMBolusVelocity
                  !                         + normalMLEvelocity 
                  !                         + normalVelocityCorrection
                  ! This is the velocity used in advective terms for layerThickness
                  ! and tracers in tendency calls in stage 3.
                  normalTransportVelocity(k,iEdge) = edgeMask(k,iEdge)*( &
                             normalTransportVelocity(k,iEdge) + &
                             normalVelocityCorrection)
               enddo

            end do ! iEdge
!$acc end parallel
            !$omp end do
            !$omp end parallel

            deallocate(uTemp)

            call mpas_timer_stop('btr se ssh verif')

         endif ! split_explicit

         call mpas_timer_stop("se btr vel")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 3: Tracer, density, pressure, vert velocity prediction
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! only compute tendencies for active tracers on last large iteration
         if (splitExplicitStep < numTSIterations) then
            activeTracersOnly = .true.
         else
            activeTracersOnly = .false.
         endif

         ! Thickness tendency computations and thickness halo updates
         ! are completed before tracer tendency computations to allow
         ! monotonic advection.

         call mpas_timer_start('se thick tend')

         ! compute vertAleTransportTop.  Use normalTransportVelocity
         ! for advection of layerThickness and tracers.
         ! Use time level 1 values of layerThickness and
         ! layerThickEdgeFlux because layerThickness has not yet
         ! been computed for time level 2.
         call mpas_timer_start('thick vert trans vel top')
         if (associated(highFreqThicknessNew)) then
            call ocn_vert_transport_velocity_top_gpu( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalTransportVelocity, &
                 sshCur, dt, vertAleTransportTop, &
                 err, highFreqThicknessNew)
         else
            call ocn_vert_transport_velocity_top_gpu( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalTransportVelocity, &
                 sshCur, dt, vertAleTransportTop, err)
         endif
         call mpas_timer_stop('thick vert trans vel top')

         call ocn_tend_thick_gpu(tendPool, forcingPool)

         call mpas_timer_stop('se thick tend')

         ! update halo for thickness tendencies
         call mpas_timer_start("se halo thickness")

         call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

         call mpas_timer_stop("se halo thickness")

         call mpas_timer_start('se tracer tend', .false.)

         call ocn_tend_tracer_gpu(tendPool, statePool, forcingPool, &
                              meshPool, swForcingPool, &
                              dt, activeTracersOnly, 2)

         call mpas_pool_get_array_gpu(tracersPool, 'activeTracers', &
                                                tracersGroupCur, 1)
         call mpas_pool_get_array_gpu(tracersPool, 'activeTracers', &
                                                tracersGroupNew, 2)
         call mpas_pool_get_array_gpu(statePool,   'normalVelocity', &
                                                normalVelocityCur, 1)
         call mpas_pool_get_array_gpu(tracersTendPool,'activeTracersTend', &
                                                   activeTracersTend)


         call mpas_timer_stop('se tracer tend')

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers")

         call mpas_pool_begin_iteration(tracersTendPool)
         do while (mpas_pool_get_next_member(tracersTendPool, &
                                             groupItr) )
            if (groupItr%memberType == MPAS_POOL_FIELD ) then
               ! Only compute tendencies for active tracers if
               ! activeTracersOnly flag is true.
               if (.not. activeTracersOnly .or. &
                   trim(groupItr%memberName)=='activeTracersTend') then
                  call mpas_dmpar_field_halo_exch(domain, &
                                                  groupItr%memberName)
               end if
            end if
         end do
         call mpas_timer_stop("se halo tracers")
         if(config_enable_nonhydrostatic_mode) call ocn_tend_vvel(tendPool, statePool, &
                          forcingPool, meshPool, 2, dt)

         call mpas_timer_start('se loop fini')

         call mpas_pool_get_array_gpu(tracersPool, 'activeTracers', &
                                                tracersGroupCur, 1)
         call mpas_pool_get_array_gpu(tracersPool, 'activeTracers', &
                                                tracersGroupNew, 2)
         if(config_enable_nonhydrostatic_mode) then
            write(*,*) "config_enable_nonhydrostatic_mode NOT EXECUTING"
            call mpas_pool_get_array(statePool, 'vertVelocityNonhydro', &
                                                vertVelocityNonhydroCur, 1)
            call mpas_pool_get_array(statePool, 'vertVelocityNonhydro', &
                                                vertVelocityNonhydroNew, 2)
            call mpas_pool_get_array(tendPool, 'vertVelocityNonhydroTend', &
                                                vertVelocityNonhydroTend)
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  If iterating, reset variables for next iteration
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (splitExplicitStep < numTSIterations) then

            ! Get indices for dynamic tracers (Includes T&S).
            call mpas_pool_get_dimension_scalar(tracersPool,'activeGRP_start',&
                                                             startIndex)
            call mpas_pool_get_dimension_scalar(tracersPool,'activeGRP_end', &
                                                             endIndex)

            ! Only need T & S for earlier iterations,
            ! then all the tracers needed the last time through.

            !$omp parallel
            !$omp do schedule(runtime) private(i, k, temp_h, temp)
!$acc parallel
!$acc loop gang worker 
            do iCell = 1, nCellsAll
!$acc loop vector private(temp_h,temp)
            do k = minLevelCell(iCell), maxLevelCell(iCell)

               ! this is h_{n+1}
               temp_h = layerThicknessCur(k,iCell) + dt* &
               layerThicknessTend(k,iCell)

               ! this is h_{n+1/2}
               layerThicknessNew(k,iCell) = 0.5* &
               (layerThicknessCur(k,iCell) + temp_h)
               do i = startIndex, endIndex
                  ! This is Phi at n+1
                  temp = (tracersGroupCur(i,k,iCell)* &
                          layerThicknessCur(k,iCell) + dt* &
                          activeTracersTend(i,k,iCell))/temp_h

                  ! This is Phi at n+1/2
                  tracersGroupNew(i,k,iCell) = 0.5_RKIND* &
                            (tracersGroupCur(i,k,iCell) + temp)
               end do ! tracer index
            end do ! vertical
            end do ! iCell
!$acc end parallel
            !$omp end do
            !$omp end parallel

            if(config_enable_nonhydrostatic_mode) then
!NV!          !!! NOT EXECUTING !!!
              write(*,*) "config_enable_nonhydrostatic_mode NOT EXECUTING"
              !$omp parallel
              !$omp do schedule(runtime) private(k, temp_h)
              do iCell = 1, nCells
                 do k = 2, maxLevelCell(iCell)
                    temp_h = vertVelocityNonhydroCur(k,iCell) + dt * vertVelocityNonhydroTend(k,iCell)
                    vertVelocityNonhydroNew(k,iCell) = 0.5_RKIND*(vertVelocityNonhydroCur(k,iCell) + temp_h)
                 end do
              end do
              !$omp end do
              !$omp end parallel
            end if

            if (config_use_freq_filtered_thickness) then
!NV!          !!! NOT EXECUTING !!!
              write(*,*) "config_use_freq_filtered_thickness NOT EXECUTING"
               !$omp parallel
               !$omp do schedule(runtime) private(k, temp)
               do iCell = 1, nCellsAll
               do k = minLevelCell(iCell), maxLevelCell(iCell)

                  ! h^{hf}_{n+1} was computed in Stage 1

                  ! this is h^{hf}_{n+1/2}
                  highFreqThicknessNew(k,iCell) = 0.5_RKIND* &
                             (highFreqThicknessCur(k,iCell) + &
                              highFreqThicknessNew(k,iCell))

                  ! this is D^{lf}_{n+1}
                  temp = lowFreqDivergenceCur(k,iCell) + dt* &
                         lowFreqDivergenceTend(k,iCell)

                  ! this is D^{lf}_{n+1/2}
                  lowFreqDivergenceNew(k,iCell) = 0.5_RKIND* &
                           (lowFreqDivergenceCur(k,iCell) + temp)
               end do
               end do
               !$omp end do
               !$omp end parallel
            end if

            !$omp parallel
            !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang vector collapse(2)
            do iEdge = 1, nEdgesAll
            do k = 1, nVertLevels

               ! u = normalBarotropicVelocity + normalBaroclinicVelocity
               ! here normalBaroclinicVelocity is at time n+1/2
               ! This is u used in next iteration or step
               normalVelocityNew(k,iEdge) = edgeMask(k,iEdge)* &
                         ( normalBarotropicVelocityNew(iEdge) + &
                           normalBaroclinicVelocityNew(k,iEdge) )

            enddo
            end do ! iEdge
!$acc end parallel
            !$omp end do
            !$omp end parallel

            ! Efficiency note: We really only need this to compute
            ! layerThickEdgeFlux, density, pressure, and SSH
            ! in this diagnostics solve.
            call ocn_diagnostic_solve_gpu(dt, statePool, forcingPool, &
                                      meshPool, scratchPool, &
                                      tracersPool, 2, full=.false.)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If large iteration complete, compute all variables at
         ! time n+1
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         elseif (splitExplicitStep == numTSIterations) then

            !$omp parallel
            !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
            do iCell = 1, nCellsAll
!$acc loop vector
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! this is h_{n+1}
               layerThicknessNew(k,iCell) = &
                                     layerThicknessCur(k,iCell) + &
                                  dt*layerThicknessTend(k,iCell)
            end do
            end do
!$acc end parallel
            !$omp end do
            !$omp end parallel

            if(config_enable_nonhydrostatic_mode) then
!NV!           !!! NOT EXECUTING !!!
               write(*,*) "config_enable_nonhydrostatic_mode NOT EXECUTING"
               !$omp parallel
               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCells
                  do k = 1, maxLevelCell(iCell)
                     vertVelocityNonhydroNew(k,iCell) = vertVelocityNonhydroCur(k,iCell) + dt * vertVelocityNonhydroTend(k,iCell)
                  end do
               end do
               !$omp end do
               !$omp end parallel
            end if

            if (config_compute_active_tracer_budgets) then

               !$omp parallel
               !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
               do iEdge = 1, nEdgesAll
!$acc loop vector
               do k= minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  do iTracer = 1, size(activeTracerHorizontalAdvectionEdgeFlux,1)
                     activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) = &
                     activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) / &
                          layerThickEdgeFlux(k,iEdge)
                  enddo
               enddo
               enddo
!$acc end parallel
               !$omp end do

               !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
               do iCell = 1, nCellsAll
!$acc loop vector
               do k= minLevelCell(iCell), maxLevelCell(iCell)
               do iTracer = 1, size(activeTracerHorizontalAdvectionTendency,1)
                  activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) = &
                  activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) / &
                            layerThicknessNew(k,iCell)

                  activeTracerVerticalAdvectionTendency(iTracer,k,iCell) = &
                  activeTracerVerticalAdvectionTendency(iTracer,k,iCell) / &
                            layerThicknessNew(k,iCell)

                  activeTracerHorMixTendency(iTracer,k,iCell) = &
                  activeTracerHorMixTendency(iTracer,k,iCell) / &
                             layerThicknessNew(k,iCell)

                  activeTracerSurfaceFluxTendency(iTracer,k,iCell) = &
                  activeTracerSurfaceFluxTendency(iTracer,k,iCell) / &
                             layerThicknessNew(k,iCell)

                  activeTracerNonLocalTendency(iTracer,k,iCell) = &
                  activeTracerNonLocalTendency(iTracer,k,iCell) / &
                             layerThicknessNew(k,iCell)
               end do

                  temperatureShortWaveTendency(k,iCell) = &
                  temperatureShortWaveTendency(k,iCell) / &
                             layerThicknessNew(k,iCell)
               end do
               end do
!$acc end parallel
               !$omp end do
               !$omp end parallel
            endif


            call mpas_pool_begin_iteration(tracersPool)
            do while (mpas_pool_get_next_member(tracersPool, &
                                                groupItr) )
               if (groupItr%memberType == MPAS_POOL_FIELD) then
                  configName = 'config_use_'//trim(groupItr%memberName)
                  call mpas_pool_get_config(domain%configs, &
                                configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array_gpu(tracersPool, &
                                              groupItr%memberName, &
                                              tracersGroupCur, 1)
                     call mpas_pool_get_array_gpu(tracersPool, &
                                              groupItr%memberName, &
                                              tracersGroupNew, 2)

                     modifiedGroupName = &
                             trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array_gpu(tracersTendPool, &
                                              modifiedGroupName, &
                                              tracersGroupTend)

                     !$omp parallel
                     !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
                     do iCell = 1, nCellsAll
!$acc loop vector
                     do k = minLevelCell(iCell), maxLevelCell(iCell)
                     do iTracer = 1, size(tracersGroupNew,1)
                        tracersGroupNew(iTracer,k,iCell) = &
                       (tracersGroupCur(iTracer,k,iCell) * &
                        layerThicknessCur(k,iCell) + dt* &
                       tracersGroupTend(iTracer,k,iCell))/ &
                        layerThicknessNew(k,iCell)
                     end do
                     end do
                     end do
!$acc end parallel
                     !$omp end do
                     !$omp end parallel

                     ! limit salinity in separate loop
                     if (trim(groupItr%memberName) == &
                         'activeTracers' ) then
                        !$omp parallel
                        !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
                        do iCell = 1, nCellsAll
!$acc loop vector
                        do k = minLevelCell(iCell), maxLevelCell(iCell)
                           tracersGroupNew(indexSalinity,k,iCell) = &
                           max(0.001_RKIND,  &
                           tracersGroupNew(indexSalinity,k,iCell))
                        end do
                        end do
!$acc end parallel
                        !$omp end do
                        !$omp end parallel
                     end if

                     ! Reset debugTracers to fixed value at the surface
                     if (trim(groupItr % memberName) == &
                         'debugTracers' .and. &
                         config_reset_debugTracers_near_surface) then
!NV!                    !!! NOT EXECUTING !!!
                        write(*,*) "config_reset_debugTracers_near_surface NOT EXECUTING"

                        !$omp parallel
                        !$omp do schedule(runtime) private(k, lat)
                        do iCell = 1, nCellsAll

                           ! Reset tracer1 to 2 in top n layers
                           do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                tracersGroupNew(1,k,iCell) = 2.0_RKIND
                           end do

                           ! Reset tracer2 to 2 in top n layers
                           ! in zonal bands, and 1 outside
                           lat = latCell(iCell)*180./3.1415
                           if (     lat>-60.0.and.lat<-55.0 &
                                .or.lat>-40.0.and.lat<-35.0 &
                                .or.lat>- 2.5.and.lat<  2.5 &
                                .or.lat> 35.0.and.lat< 40.0 &
                                .or.lat> 55.0.and.lat< 60.0 ) then
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(2,k,iCell) = 2.0_RKIND
                              end do
                           else
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(2,k,iCell) = 1.0_RKIND
                              end do
                           end if

                           ! Reset tracer3 to 2 in top n layers
                           ! in zonal bands, and 1 outside
                           lat = latCell(iCell)*180./3.1415
                           if (     lat>-55.0.and.lat<-50.0 &
                                .or.lat>-35.0.and.lat<-30.0 &
                                .or.lat>-15.0.and.lat<-10.0 &
                                .or.lat> 10.0.and.lat< 15.0 &
                                .or.lat> 30.0.and.lat< 35.0 &
                                .or.lat> 50.0.and.lat< 55.0 ) then
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(3,k,iCell) = 2.0_RKIND
                              end do
                           else
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(3,k,iCell) = 1.0_RKIND
                              end do
                           end if
                        end do ! cells
                        !$omp end do
                        !$omp end parallel
                     end if ! debug tracers
                  end if ! use tracer group
               end if ! tracer
            end do ! tracer group

            if (config_use_freq_filtered_thickness) then
!NV!           !!! NOT EXECUTING !!!
               write(*,*) "config_use_freq_filtered_thickness NOT EXECUTING"
               !$omp parallel
               !$omp do schedule(runtime) private(k)
               do iCell = 1, nCellsAll
               do k = minLevelCell(iCell), maxLevelCell(iCell)

                  ! h^{hf}_{n+1} was computed in Stage 1

                  ! this is D^{lf}_{n+1}
                  lowFreqDivergenceNew(k,iCell) = &
                  lowFreqDivergenceCur(k,iCell) + dt* &
                  lowFreqDivergenceTend(k,iCell)
               end do
               end do
               !$omp end do
               !$omp end parallel
            end if

            ! Recompute final u to go on to next step.
            ! u_{n+1} = normalBarotropicVelocity_{n+1} +
            !           normalBaroclinicVelocity_{n+1}
            ! Right now normalBaroclinicVelocityNew is at time n+1/2,
            ! so back compute to get normalBaroclinicVelocity at
            ! time n+1 using normalBaroclinicVelocity_{n+1/2} =
            !            1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
            ! so the following lines are
            ! u_{n+1} = normalBarotropicVelocity_{n+1} +
            !         2*normalBaroclinicVelocity_{n+1/2} -
            !         normalBaroclinicVelocity_n
            ! note that normalBaroclinicVelocity is recomputed at the
            ! beginning of the next timestep due to Imp Vert mixing,
            ! so normalBaroclinicVelocity does not have to be
            ! recomputed here.

            !$omp parallel
            !$omp do schedule(runtime) private(k)
!$acc parallel
!$acc loop gang worker
            do iEdge = 1, nEdgesAll
!$acc loop vector
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               normalVelocityNew(k,iEdge) = &
                               normalBarotropicVelocityNew(iEdge) + &
                             2*normalBaroclinicVelocityNew(k,iEdge) - &
                               normalBaroclinicVelocityCur(k,iEdge)
            end do
            end do ! iEdges
!$acc end parallel
            !$omp end do
            !$omp end parallel

         endif ! splitExplicitStep


         call mpas_timer_stop('se loop fini')
         call mpas_timer_stop('se loop')

      end do  ! outer timestep loop (splitExplicitStep)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se implicit vert mix")

      ! Call ocean diagnostic solve in preparation for vertical mixing.
      ! Note it is called again after vertical mixing, because u and
      ! tracers change. For Richardson vertical mixing, only density,
      ! layerThickEdgeFlux, and kineticEnergyCell need to be computed.
      ! For kpp, more variables may be needed.  Either way, this could
      ! be made more efficient by only computing what is needed for the
      ! implicit vmix routine that follows.
      call ocn_diagnostic_solve_gpu(dt, statePool, forcingPool, meshPool, &
                                scratchPool, tracersPool, 2)

      call mpas_dmpar_field_halo_exch(domain, 'surfaceFrictionVelocity')

      call mpas_dmpar_field_halo_exch(domain, 'surfaceFrictionVelocity')

      call ocn_vmix_implicit_gpu(dt, meshPool, statePool, forcingPool, &
                             scratchPool, err, 2)

      ! Update halo on u and tracers, which were just updated for
      ! implicit vertical mixing.  If not done, this leads to lack of
      ! volume conservation.  It is required because halo updates in
      ! stage 3 are only conducted on tendencies, not on the velocity
      ! and tracer fields.  So this update is required to communicate
      ! the change due to implicit vertical mixing across the boundary.

      if(config_enable_nonhydrostatic_mode) then
!NV!    !!! NOT EXECUTING !!!
        write(*,*) "config_enable_nonhydrostatic_mode NOT EXECUTING"
        dminfo = domain % dminfo
         call ocn_nonhydrostatic_pressure_tend(normalVelocityNew, &
              vertVelocityNonhydroNew, layerThicknessNew, layerThickEdgeFlux, &
              nhPressureCorrection, dminfo % comm, dt)

         call mpas_timer_start('se NH pressure tend halo pressure')
         call mpas_dmpar_field_halo_exch(domain, 'nonhydrostaticPressure', timeLevel=2)
         call mpas_dmpar_field_halo_exch(domain, 'nonhydrostaticPressure', timeLevel=1)
         call mpas_dmpar_field_halo_exch(domain, 'nhPressureCorrection')
         call mpas_timer_stop('se NH pressure tend halo pressure')

         call ocn_nonhydrostatic_pressure_update_velocity(normalVelocityNew, &
              vertVelocityNonhydroNew, layerThicknessNew, layerThickEdgeFlux, &
              nonhydrostaticPressureNew, nonhydrostaticPressureCur, &
              nhPressureCorrection, dt)
      end if


      call mpas_timer_start('se vmix halos')

      call mpas_timer_start('se vmix halos normalVelFld')
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', &
                                      timeLevel=2)
      call mpas_timer_stop('se vmix halos normalVelFld')

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if (groupItr%memberType == MPAS_POOL_FIELD) then
            call mpas_dmpar_field_halo_exch(domain, &
                         groupItr%memberName, timeLevel=2)
            ! Reset iAge to zero where mask == 0
            if (config_use_idealAgeTracers.and.trim(groupItr%memberName) == 'idealAgeTracers') then
!NV!            !!! NOT EXECUTING !!!
                write(*,*) "config_use_idealAgeTracers NOT EXECUTING"
                call mpas_pool_get_array(tracersPool, &
                                         groupItr%memberName, &
                                         tracersGroupNew, 2)
                call mpas_pool_get_subpool(forcingPool, &
                                           'tracersIdealAgeFields', &
                                            tracersIdealAgeFieldsPool)
                call mpas_pool_get_array(tracersIdealAgeFieldsPool, &
                                         'idealAgeTracersIdealAgeMask', &
                                          tracerGroupIdealAgeMask)

                !$omp parallel
                !$omp do schedule(runtime)
                do iCell = 1, nCellsOwned
                   tracersGroupNew(1,1,iCell) = tracerGroupIdealAgeMask(1,iCell)*tracersGroupNew(1,1,iCell)
                end do ! cells
                !$omp end do
                !$omp end parallel
             endif

            ! Halo update on all tracer groups
            call mpas_dmpar_field_halo_exch(domain, groupItr%memberName, timeLevel=2)

         end if
      end do

      call mpas_timer_stop('se vmix halos')

      call mpas_timer_stop("se implicit vert mix")

      if ( configVertAdvMethod == vertAdvRemap ) then
!NV!     !!! NOT EXECUTING !!!
         write(*,*) "configVertAdvMethod == vertAdvRemap NOT EXECUTING"
         
         call mpas_timer_start("vertical remap")

         call mpas_pool_get_array(statePool, 'layerThicknessLag', layerThicknessLagNew, 2)

         ! Store the current, Lagrangian location for computation of the
         ! Eulerian vertical velocity induced by remapping
         layerThicknessLagNew = layerThicknessNew

         ! Perform the remapping
         call ocn_remap_vert_state(block, err)

         call mpas_timer_stop("vertical remap")

      end if

      call mpas_timer_start('se fini')


      if (config_prescribe_velocity) then
!NV!     !!! NOT EXECUTING !!!
         write(*,*) "config_prescribe_velocity NOT EXECUTING"
         !$omp parallel
         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdgesAll
         do k=1,nVertLevels
            normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)
         end do
         end do
         !$omp end do
         !$omp end parallel
      end if

      if (config_prescribe_thickness) then
!NV!     !!! NOT EXECUTING !!!
         write(*,*) "config_prescribe_thickness NOT EXECUTING"
         !$omp parallel
         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCellsAll
         do k=1,nVertLevels
            layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
         end do
         end do
         !$omp end do
         !$omp end parallel
      end if

      call ocn_diagnostic_solve_gpu(dt, statePool, forcingPool, meshPool, &
                                scratchPool, tracersPool, 2)

      ! Update the effective desnity in land ice if we're coupling to
      ! land ice
!NV!  !!! NOT EXECUTING !!!
      call ocn_effective_density_in_land_ice_update_gpu(forcingPool, &
                                                    statePool, err)

      call mpas_timer_start('se final mpas reconstruct gpu', .false.)

      call mpas_reconstruct_gpu(meshPool, normalVelocityNew,       &
                                velocityX, velocityY, velocityZ,   &
                                velocityZonal, velocityMeridional, &
                                includeHalos = .true.)

      call mpas_reconstruct_gpu(meshPool, gradSSH,               &
                                gradSSHX, gradSSHY, gradSSHZ,    &
                                gradSSHZonal, gradSSHMeridional, &
                                includeHalos = .true.)

      call mpas_timer_stop('se final mpas reconstruct gpu')

      !$omp parallel
      !$omp do schedule(runtime)
!$acc parallel
!$acc loop gang vector
      do iCell = 1, nCellsAll
         surfaceVelocity(indexSurfaceVelocityZonal,iCell) = &
                                   velocityZonal(minLevelCell(iCell),iCell)
         surfaceVelocity(indexSurfaceVelocityMeridional,iCell) = &
                                   velocityMeridional(minLevelCell(iCell),iCell)

         SSHGradient(indexSSHGradientZonal,iCell) = gradSSHZonal(iCell)
         SSHGradient(indexSSHGradientMeridional,iCell) = &
                                               gradSSHMeridional(iCell)
      end do
!$acc end parallel
      !$omp end do
      !$omp end parallel

      call ocn_time_average_coupled_accumulate_gpu(statePool,forcingPool,2)

      if (config_use_GM .or. config_submesoscale_enable) then
         call ocn_reconstruct_eddy_vectors_gpu(meshPool)
      end if

      if (trim(config_land_ice_flux_mode) == 'coupled') then
!NV!     !!! NOT EXECUTING !!!
         write(*,*) "config_land_ice_flux_mode NOT EXECUTING"
         call mpas_timer_start("se effective density halo")
         call mpas_pool_get_field(statePool, &
                                 'effectiveDensityInLandIce', &
                                  effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("se effective density halo")
      end if
      deallocate(bottomDepthEdge)

      call mpas_timer_stop('se fini')
      call mpas_timer_stop("se timestep")

!!!$acc update host(pressure,thermExpCoeff,salineContractCoeff,zMid,density,potentialDensity)

   end subroutine ocn_time_integrator_split!}}}

!***********************************************************************
!
!  routine ocn_time_integration_split_init
!
!> \brief   Initialize split-explicit time stepping within ocean
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes variables required for the split-explicit
!>  method of integrating the ocean model forward in time.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_split_init(domain)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain      !< [inout] data structure containing most variables

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: &
         statePool,      &! structure containing model state
         meshPool,       &! structure containing mesh fields
         tracersPool      ! structure containing tracer fields

      integer ::         &
         iCell, iEdge, k,&! loop indices for cell, edge and vertical
         kmax,           &! index of deepest active edge
         ierr,           &! local error flag
         cell1, cell2     ! neighbor cell indices across edge

      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) ::       &
         normalThicknessFluxSum, &! vertical sum of thick flux
         layerThicknessSum,      &! vertical sum of layer thickness
         layerThicknessEdge1      ! layer thickness on edge

      real (kind=RKIND), dimension(:), pointer :: &
         refBottomDepth,         &! reference bottom depth
         normalBarotropicVelocity ! normal barotropic velocity

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness,          &! layer thickness cell center
         normalBaroclinicVelocity,&! normal baroclinic velocity
         normalVelocity            ! normal velocity (total)

      type (mpas_time_type) :: &
         nowTime          ! current model time

      type (mpas_timeInterval_type) :: &! various time step calculations
         fullTimeStep,       &! full model time step
         barotropicTimeStep, &! barotropic time step
         remainder,          &! remaining time after interval division
         zeroInterval         ! zero timestep for comparing remainder

      integer, pointer :: nVertLevels, nCells, nEdges

      integer (kind=I8KIND) :: nBtrSubcyclesI8

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Determine the number of barotropic subcycles based on the
      !*** ratio of time steps

      nowTime = mpas_get_clock_time(domain%clock, MPAS_NOW, ierr)

      call mpas_set_timeInterval( zeroInterval, S=0 )
      call mpas_set_timeInterval( fullTimeStep , timeString=config_dt )
      call mpas_set_timeInterval( barotropicTimeStep, &
                                             timeString=config_btr_dt )

      !*** transfer to I8 for division step
      nBtrSubcyclesI8 = nBtrSubcycles
      call mpas_interval_division( nowTime, fullTimeStep, &
                       barotropicTimeStep, nBtrSubcyclesI8, remainder )
      nBtrSubcycles = nBtrSubcyclesI8

      if ( remainder > zeroInterval ) then
         nBtrSubcycles = nBtrSubcycles + 1
      end if

      !*** Set mask for using velocity correction
      if (config_vel_correction) then
         useVelocityCorrection = 1.0_RKIND
      else
         useVelocityCorrection = 0.0_RKIND
      endif

      !*** Compute some halo needs to reduce the number
      !*** of halo updates during subcycling
      neededHalos = 1 + config_n_btr_cor_iter
      haloDecrement = mod(config_num_halos, neededHalos)

      !*** Determine the time integration type and set associated masks

      select case (trim(config_time_integrator))

      case ('unsplit_explicit')
         unsplit   = .true.
         splitFact = 0.0_RKIND

      case ('split_explicit')
         unsplit   = .false.
         splitFact = 1.0_RKIND
         call mpas_log_write( '*******************************************************************************')
         call mpas_log_write( 'The split explicit time integration is configured to use: $i barotropic subcycles', &
            intArgs=(/ nBtrSubcycles /) )
         call mpas_log_write( '*******************************************************************************')

      case default
         call mpas_log_write('Incorrect choice config_time_integrator',&
                             MPAS_LOG_CRIT)

      end select

      !*** set number of baroclinic iterations on each outer
      !*** time step iteration (number can be different on the
      !*** first and last time step iteration)

      numTSIterations     = config_n_ts_iter
      allocate(numClinicIterations(numTSIterations))

      numClinicIterations    = config_n_bcl_iter_mid ! most iterations
      numClinicIterations(1) = config_n_bcl_iter_beg ! first iteration
      numClinicIterations(numTSIterations)=config_n_bcl_iter_end !last

      !*** Initialize some variables if not a restart

      if (.not. config_do_restart) then

         !*** Initialize z-level mesh variables from h, read in from
         !*** input file.

         block => domain % blocklist
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_array_gpu(statePool, 'layerThickness', &
                                              layerThickness, 1)
         call mpas_pool_get_array_gpu(statePool, 'normalVelocity', &
                                              normalVelocity, 1)
         call mpas_pool_get_array_gpu(statePool, 'normalBarotropicVelocity', &
                                              normalBarotropicVelocity, 1)
         call mpas_pool_get_array_gpu(statePool, 'normalBaroclinicVelocity', &
                                              normalBaroclinicVelocity, 1)

         call mpas_pool_get_array_gpu(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array_gpu(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array_gpu(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array_gpu(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
         call mpas_pool_get_array_gpu(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
!$acc update host(layerThickness,normalVelocity,normalBarotropicVelocity,normalBaroclinicVelocity,&
!$acc             refBottomDepth,cellsOnEdge,minLevelCell,minLevelEdgeBot,maxLevelEdgeTop)

         !*** Compute barotropic velocity at first timestep
         !*** This is only done upon start-up.
         if (unsplit) then

            do iEdge = 1, nEdgesAll
               normalBarotropicVelocity(  iEdge) = 0.0_RKIND
               normalBaroclinicVelocity(:,iEdge) = &
                         normalVelocity(:,iEdge)
            end do

         else ! split explicit

            if (config_filter_btr_mode) then
               do iCell = 1, nCells
                  layerThickness(minLevelCell(iCell),iCell) = refBottomDepth(minLevelCell(iCell))
               enddo
            endif

            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               kmax  = maxLevelEdgeTop(iEdge)

               ! normalBarotropicVelocity = sum(h*u)/sum(h) on each edge
               ! ocn_diagnostic_solve has not yet been called, so
               ! compute hEdge just for this edge.

               ! thicknessSum is initialized outside the loop because on
               ! land boundaries maxLevelEdgeTop=0, but we want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               layerThicknessEdge1 = 0.5_RKIND* &
                                     (layerThickness(minLevelCell(cell1),cell1) + &
                                      layerThickness(minLevelCell(cell2),cell2) )
               normalThicknessFluxSum = layerThicknessEdge1* &
                                        normalVelocity(minLevelEdgeBot(iEdge),iEdge)
               layerThicknessSum = layerThicknessEdge1

               do k=minLevelEdgeBot(iEdge)+1, kmax
                  layerThicknessEdge1 = 0.5_RKIND* &
                                       (layerThickness(k,cell1) + &
                                        layerThickness(k,cell2))

                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                           layerThicknessEdge1* &
                                           normalVelocity(k,iEdge)
                  layerThicknessSum = layerThicknessSum + &
                                      layerThicknessEdge1

               enddo
               normalBarotropicVelocity(iEdge) = &
                     normalThicknessFluxSum/layerThicknessSum

               ! normalBaroclinicVelocity = normalVelocity -
               !     normalBarotropicVelocity
               do k = minLevelEdgeBot(iEdge), kmax
                  normalBaroclinicVelocity(k,iEdge) = &
                            normalVelocity(k,iEdge) - &
                    normalBarotropicVelocity(iEdge)
               enddo

               ! normalBaroclinicVelocity=0,
               ! normalVelocity=0 on land cells
               do k = kmax+1, nVertLevels
                  normalBaroclinicVelocity(k,iEdge) = 0.0_RKIND
                  normalVelocity(k,iEdge) = 0.0_RKIND
               enddo
            enddo ! edge loop

            if (config_filter_btr_mode) then
               ! filter normalBarotropicVelocity out of initial condition

               normalVelocity(:,:) = normalBaroclinicVelocity(:,:)
               normalBarotropicVelocity(:) = 0.0_RKIND

            endif

         endif ! split explicit

      end if ! not restart

      ! Create reusable halo exchange for subcycling
      call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
      call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                           'sshSubcycle')
      call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                  'normalBarotropicVelocitySubcycle')
      call mpas_dmpar_exch_group_build_reusable_buffers(domain, &
                                                     subcycleGroupName)

      !-----------------------------------------------------------------
      if (config_use_self_attraction_loading) then
         ! set ssh_sal_on to 1
         ssh_sal_on = 1
      else
         ssh_sal_on = 0
      endif

      self_attraction_and_loading_beta = config_self_attraction_and_loading_beta

!$acc update device(normalBarotropicVelocity,normalBaroclinicVelocity,normalVelocity)
   end subroutine ocn_time_integration_split_init!}}}

!***********************************************************************

end module ocn_time_integration_split

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
