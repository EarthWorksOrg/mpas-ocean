! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_diagnostics_variables
!
!> \brief MPAS ocean diagnostics variable definitions
!> \author Matt Turner
!> \date   13 May 2020
!> \details
!>  This module contains the definitions of all variables contianed in the
!>  diagnosticsPool, and also gets the variables from the Pool in init.
!
!-----------------------------------------------------------------------

module ocn_diagnostics_variables

   use mpaso_derived_types
   use mpaso_pool_routines
   use mpaso_constants
   use mpaso_stream_manager
   use mpaso_io_units

   use ocn_config

   implicit none
   public
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   real (kind=RKIND), dimension(:,:), pointer :: zMid
   real (kind=RKIND), dimension(:,:), pointer :: zTop
   real (kind=RKIND), dimension(:,:), pointer :: divergence
   real (kind=RKIND), dimension(:,:), pointer :: circulation
   real (kind=RKIND), dimension(:,:), pointer :: relativeVorticity
   real (kind=RKIND), dimension(:,:), pointer :: relativeVorticityCell
   real (kind=RKIND), dimension(:,:), pointer :: normPlanetVortEdge
   real (kind=RKIND), dimension(:,:), pointer :: normRelVortEdge
   real (kind=RKIND), dimension(:,:), pointer :: normalizedRelativeVorticityCell
   real (kind=RKIND), dimension(:,:), pointer :: density
   real (kind=RKIND), dimension(:,:), pointer :: displaceddensity
   real (kind=RKIND), dimension(:,:), pointer :: potentialdensity
   real (kind=RKIND), dimension(:,:), pointer :: montgomeryPotential
   real (kind=RKIND), dimension(:,:), pointer :: pressure
   real (kind=RKIND), dimension(:,:), pointer :: thermExpCoeff
   real (kind=RKIND), dimension(:,:), pointer :: salineContractCoeff
   real (kind=RKIND), dimension(:,:), pointer :: BruntVaisalaFreqTop
   real (kind=RKIND), dimension(:,:), pointer :: tangentialVelocity
   real (kind=RKIND), dimension(:,:), pointer :: layerThickEdgeFlux
   real (kind=RKIND), dimension(:,:), pointer :: layerThickEdgeMean
   real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell
   real (kind=RKIND), dimension(:,:), pointer :: vertVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: vertTransportVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: vertGMBolusVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: vertMLEBolusVelocityTop
   real (kind=RKIND), dimension(:,:), pointer :: normalGMBolusVelocity
   real (kind=RKIND), dimension(:,:), pointer :: normalMLEvelocity
   real (kind=RKIND), dimension(:,:), pointer :: normalTransportVelocity
   real (kind=RKIND), dimension(:,:), pointer :: RiTopOfCell
   real (kind=RKIND), dimension(:,:), pointer :: RiTopOfEdge

   real (kind=RKIND), dimension(:), pointer :: gradSSH
   real (kind=RKIND), dimension(:), pointer :: pressureAdjustedSSH
   real (kind=RKIND), dimension(:), pointer :: ssh_sal

   real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue
   real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceLayerValue
   real (kind=RKIND), dimension(:),   pointer :: normalVelocitySurfaceLayer
   real (kind=RKIND), dimension(:),   pointer :: indexSurfaceLayerDepth

   real (kind=RKIND), dimension(:), pointer :: sfcFlxAttCoeff
   real (kind=RKIND), dimension(:), pointer :: surfaceFluxAttenuationCoefficientRunoff

   real (kind=RKIND), dimension(:), pointer :: landIceFrictionVelocity
   real (kind=RKIND), dimension(:), pointer :: topDrag
   real (kind=RKIND), dimension(:), pointer :: topDragMag
   real (kind=RKIND), dimension(:), pointer :: topographic_wave_drag
   real (kind=RKIND), dimension(:,:), pointer :: landIceBoundaryLayerTracers
   real (kind=RKIND), dimension(:,:), pointer :: landIceTracerTransferVelocities

   integer :: indexBLT, indexBLS, &
              indexHeatTrans, indexSaltTrans

   real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFluxOBL
   real (kind=RKIND), dimension(:), pointer :: bottomLayerShortwaveTemperatureFlux

   real (kind=RKIND), dimension(:), pointer :: surfaceBuoyancyForcing
   real (kind=RKIND), dimension(:), pointer :: surfaceFrictionVelocity
   real (kind=RKIND), dimension(:), pointer :: surfacePressure

   real (kind=RKIND), dimension(:,:), pointer :: &
      transportVelocityX, transportVelocityY, transportVelocityZ, transportVelocityZonal, &
      transportVelocityMeridional, eddyVelocityX, eddyVelocityY, &
      eddyVelocityZ, GMBolusVelocityZonal, GMBolusVelocityMeridional, gmStreamFuncTopOfEdge, &
      GMStreamFuncX, GMStreamFuncY, GMStreamFuncZ, GMStreamFuncZonal, GMStreamFuncMeridional, &
      mleVelocityZonal, mleVelocityMeridional

   real (kind=RKIND), dimension(:,:), pointer :: rx1Edge
   real (kind=RKIND), dimension(:,:), pointer :: rx1Cell
   real (kind=RKIND), dimension(:), pointer :: rx1MaxEdge
   real (kind=RKIND), dimension(:), pointer :: rx1MaxCell
   real (kind=RKIND), pointer :: globalRx1Max

   integer, dimension(:), pointer :: smoothingMask
   real (kind=RKIND), dimension(:,:), pointer :: verticalStretch

   integer, dimension(:), pointer :: modifyLandIcePressureMask

   real (kind=RKIND), dimension(:,:), pointer :: velocityMeridional
   real (kind=RKIND), dimension(:,:), pointer :: velocityZonal

   real(kind=RKIND), dimension(:,:), pointer :: gmKappaScaling
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaScaling
   real(kind=RKIND), dimension(:,:), pointer :: RediKappaSfcTaper
   real(kind=RKIND), dimension(:,:), pointer :: k33
   real(kind=RKIND), dimension(:,:), pointer :: gradDensityEddy
   real(kind=RKIND), dimension(:), pointer   :: gmBolusKappa
   real(kind=RKIND), dimension(:), pointer   :: cGMphaseSpeed
   real(kind=RKIND), dimension(:,:,:), pointer :: slopeTriadUp, slopeTriadDown

   integer, dimension(:), pointer :: indMLD
   real(kind=RKIND), dimension(:), pointer :: dThreshMLD
   real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ

   real (kind=RKIND), dimension(:,:), pointer :: vertAleTransportTop

   real (kind=RKIND), dimension(:,:), pointer :: viscosity
   real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfEdge
   real (kind=RKIND), dimension(:,:), pointer :: vertViscTopOfCell
   real (kind=RKIND), dimension(:,:), pointer :: wettingVelocity
   type (field2DReal), pointer :: wettingVelocityField

   real (kind=RKIND), dimension(:,:), pointer :: vertDiffTopOfCell
   real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
   integer, dimension(:,:), pointer :: rediLimiterCount
   real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerSurfaceFluxTendency
   real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerNonLocalTendency
   real (kind=RKIND), dimension(:,:,:), pointer :: activeTracerHorMixTendency
   real (kind=RKIND), dimension(:,:), pointer :: temperatureShortWaveTendency
   real (kind=RKIND), dimension(:), pointer :: RediKappa

   ! pointers for tendencies used in diagnostic budget computation
   real (kind=RKIND), dimension(:,:,:), pointer, contiguous :: &
      activeTracerHorizontalAdvectionTendency, &
      activeTracerVerticalAdvectionTendency,   &
      activeTracerHorizontalAdvectionEdgeFlux, &
      activeTracerVerticalAdvectionTopFlux,    &
      activeTracerVertMixTendency

   real (kind=RKIND), pointer :: daysSinceStartOfSim
   real (kind=RKIND), pointer :: Time
   real (kind=RKIND), dimension(:), pointer :: Time_bnds
   character (len=StrKIND), pointer :: xtime, simulationStartTime

   real (kind=RKIND), dimension(:,:), pointer :: bulkRichardsonNumber
   real (kind=RKIND), dimension(:,:), pointer :: bulkRichardsonNumberBuoy
   real (kind=RKIND), dimension(:,:), pointer :: bulkRichardsonNumberShear

   real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepth
   real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepthSmooth
   real (kind=RKIND), dimension(:), pointer :: indexBoundaryLayerDepth
   integer :: index_vertNonLocalFluxTemp

   real (kind=RKIND), dimension(:), pointer :: salinitySurfaceRestoringTendency

   real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity
   real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
   real (kind=RKIND), dimension(:), pointer :: gradSSHX
   real (kind=RKIND), dimension(:), pointer :: gradSSHY
   real (kind=RKIND), dimension(:), pointer :: gradSSHZ

   real (kind=RKIND), dimension(:), pointer :: barotropicForcing
   real (kind=RKIND), dimension(:), pointer :: barotropicCoriolisTerm
   real (kind=RKIND), dimension(:), pointer :: barotropicThicknessFlux

   real (kind=RKIND), dimension(:, :), pointer :: SSHGradient
   integer :: indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
   integer :: indexSSHGradientZonal, indexSSHGradientMeridional

   type (field1DInteger), pointer :: smoothingMaskField
   type (field2DReal), pointer :: verticalStretchField

   real (kind=RKIND), pointer :: globalVerticalStretchMax, globalVerticalStretchMin
   real (kind=RKIND), dimension(:), pointer :: betaEdge
   real (kind=RKIND), dimension(:), pointer :: RediHorizontalTaper, gmHorizontalTaper, RossbyRadius

   ! Semi-implicit Array Pointers
   real (kind=RKIND), dimension(:), pointer :: SIvec_r0,SIvec_r00,SIvec_r1 ,SIvec_rh0,SIvec_rh1,SIvec_ph0,SIvec_ph1
   real (kind=RKIND), dimension(:), pointer :: SIvec_v0,SIvec_s0 ,SIvec_s1 ,SIvec_sh0,SIvec_sh1
   real (kind=RKIND), dimension(:), pointer :: SIvec_t0,SIvec_t1 ,SIvec_q0 ,SIvec_q1 ,SIvec_qh0
   real (kind=RKIND), dimension(:), pointer :: SIvec_w0,SIvec_w1 ,SIvec_wh0,SIvec_wh1,SIvec_y0
   real (kind=RKIND), dimension(:), pointer :: SIvec_z0,SIvec_z1 ,SIvec_zh0,SIvec_zh1

   real (kind=RKIND), dimension(:,:), pointer :: nhPressureCorrection

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_diagnostics_variables_init
   public :: ocn_diagnostics_variables_destroy

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_diagnostics_variables_init
!
!> \brief   Initializes variables in the diagnosticsPool
!> \author  Matt Turner
!> \date    13 May 2020
!> \details
!>  This routine initializes all of the pointers associated with
!>  variables contained in the diagnosticsPool.
!
!-----------------------------------------------------------------------
   subroutine ocn_diagnostics_variables_init(domain, jenkinsOn, hollandJenkinsOn, err)!{{{
      type (domain_type), intent(in) :: domain

      logical, intent(in) :: jenkinsOn
      logical, intent(in) :: hollandjenkinsOn

      integer, intent(out) :: err !< Output: Error flag

      type (mpas_pool_type), pointer :: diagnosticsPool

      ! Local pointers for pool retrievals
      integer, pointer :: &
         indexBLTPtr, indexBLSPtr, &
         indexHeatTransPtr, indexSaltTransPtr, &
         index_vertNonLocalFluxTempPtr, &
         indexSurfaceVelocityZonalPtr, &
         indexSurfaceVelocityMeridionalPtr, &
         indexSSHGradientZonalPtr, &
         indexSSHGradientMeridionalPtr

      err = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)

      if (config_use_topographic_wave_drag) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'topographic_wave_drag', topographic_wave_drag)
      end if

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'wettingVelocity', &
                   wettingVelocity)
         call mpas_pool_get_field(diagnosticsPool, 'wettingVelocity', &
                   wettingVelocityField)
      end if

      if (config_use_GM.or.config_use_Redi) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'RediKappa', &
                   RediKappa)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'RediKappaScaling', &
                   RediKappaScaling)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'RediKappaSfcTaper', &
                   RediKappaSfcTaper)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'rediLimiterCount', &
                   rediLimiterCount)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'gmBolusKappa',  &
                   gmBolusKappa)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gmKappaScaling', &
                   gmKappaScaling)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'RediHorizontalTaper', &
                   RediHorizontalTaper)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gmHorizontalTaper', &
                   gmHorizontalTaper)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'RossbyRadius', &
                   RossbyRadius)
      end if

      if (config_use_GM .or. config_submesoscale_enable) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradDensityEddy', &
                   gradDensityEddy)
      end if

      if ( config_compute_active_tracer_budgets ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerSurfaceFluxTendency', &
                   activeTracerSurfaceFluxTendency)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'temperatureShortWaveTendency', &
                   temperatureShortWaveTendency)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerNonLocalTendency', &
                   activeTracerNonLocalTendency)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerHorMixTendency', &
                   activeTracerHorMixTendency)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerVerticalAdvectionTopFlux', &
                   activeTracerVerticalAdvectionTopFlux)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerHorizontalAdvectionEdgeFlux', &
                   activeTracerHorizontalAdvectionEdgeFlux)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerHorizontalAdvectionTendency', &
                   activeTracerHorizontalAdvectionTendency)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerVerticalAdvectionTendency',  &
                   activeTracerVerticalAdvectionTendency)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'activeTracerVertMixTendency', &
                   activeTracerVertMixTendency)
      end if

      if ( config_use_activeTracers_surface_restoring ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'salinitySurfaceRestoringTendency', &
                   salinitySurfaceRestoringTendency)
      end if

      ! Retrieve pointers for configurations where land_ice_flux_mode is enabled
      if ( trim(config_land_ice_flux_mode) /= 'off' ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'topDrag', topDrag)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'topDragMagnitude', topDragMag)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceBoundaryLayerTemperature', indexBLTPtr)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceBoundaryLayerSalinity', indexBLSPtr)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceHeatTransferVelocity', indexHeatTransPtr)
         call mpas_pool_get_dimension(diagnosticsPool, &
                   'index_landIceSaltTransferVelocity', indexSaltTransPtr)
         indexBLT = indexBLTPtr
         indexBLS = indexBLSPtr
         indexHeatTrans = indexHeatTransPtr
         indexSaltTrans = indexSaltTransPtr

         call mpas_pool_get_array_gpu(diagnosticsPool, &
                   'landIceBoundaryLayerTracers', landIceBoundaryLayerTracers)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                   'landIceTracerTransferVelocities', &
                    landIceTracerTransferVelocities)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                   'landIceFrictionVelocity', landIceFrictionVelocity)
      end if

      if ( trim(config_ocean_run_mode) == 'init' ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'rx1Edge', rx1Edge, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'rx1Cell', rx1Cell, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'rx1MaxEdge', rx1MaxEdge, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'rx1InitSmoothingMask', &
                   smoothingMask, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'verticalStretch', &
                   verticalStretch, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'globalVerticalStretchMax', &
                   globalVerticalStretchMax, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'globalVerticalStretchMin', &
                   globalVerticalStretchMin, 1)
      end if

      if (config_use_topographic_wave_drag) then
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'topographic_wave_drag', &
                   topographic_wave_drag)
      end if

      if (config_use_tidal_potential_forcing) then
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'ssh_sal', &
                   ssh_sal)
      end if

      if ( trim(config_ocean_run_mode) == 'init' .or. &
                config_AM_debugDiagnostics_enable ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'rx1MaxCell', &
                   rx1MaxCell, 1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'globalRx1Max', &
                   globalRx1Max, 1)
      end if

      if ( trim(config_ocean_run_mode) == 'forward' .or. &
           trim(config_ocean_run_mode) == 'analysis' ) then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'k33', k33)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradSSH', gradSSH)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradSSHX', gradSSHX)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradSSHY', gradSSHY)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradSSHZ', gradSSHZ)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradSSHZonal', &
                   gradSSHZonal)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gradSSHMeridional', &
                   gradSSHMeridional)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SSHGradient', SSHGradient)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'viscosity', viscosity)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'divergence', divergence)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'circulation', circulation)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'velocityX', velocityX)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'velocityY', velocityY)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'velocityZ', velocityZ)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'velocityMeridional', &
                   velocityMeridional)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'velocityZonal', &
                   velocityZonal)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'RiTopOfEdge', RiTopOfEdge)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'slopeTriadUp', &
                   slopeTriadUp)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'slopeTriadDown', &
                   slopeTriadDown)

         ! GM-related fields
         call mpas_pool_get_array_gpu(diagnosticsPool, 'cGMphaseSpeed', &
                   cGMphaseSpeed)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'gmStreamFuncTopOfEdge', &
                   gmStreamFuncTopOfEdge)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMStreamFuncX', &
                   GMStreamFuncX)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMStreamFuncY', &
                   GMStreamFuncY)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMStreamFuncZ', &
                   GMStreamFuncZ)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMStreamFuncZonal', &
                   GMStreamFuncZonal)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMStreamFuncMeridional', &
                   GMStreamFuncMeridional)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'mleVelocityZonal', &
                   mleVelocityZonal)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'mleVelocityMeridional', &
                   mleVelocityMeridional)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'eddyVelocityX', &
                   eddyVelocityX)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'eddyVelocityY', &
                   eddyVelocityY)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'eddyVelocityZ', &
                   eddyVelocityZ)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMBolusVelocityZonal', &
                   GMBolusVelocityZonal)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'GMBolusVelocityMeridional', &
                   GMBolusVelocityMeridional)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'normalGMBolusVelocity', &
                   normalGMBolusVelocity)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'normalMLEvelocity', &
                   normalMLEvelocity)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertGMBolusVelocityTop', &
                   vertGMBolusVelocityTop)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertMLEBolusVelocityTop', &
                   vertMLEBolusVelocityTop)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'layerThicknessEdgeFlux', &
                   layerThickEdgeFlux)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'layerThicknessEdgeMean', &
                   layerThickEdgeMean)

         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'normalizedRelativeVorticityEdge', &
                   normRelVortEdge)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'normalizedPlanetaryVorticityEdge', &
                   normPlanetVortEdge)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'relativeVorticity', &
                   relativeVorticity)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'relativeVorticityCell', &
                   relativeVorticityCell)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'surfaceVelocity', &
                   surfaceVelocity)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertVelocityTop', &
                   vertVelocityTop)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertNonLocalFlux', &
                   vertNonLocalFlux)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertViscTopOfEdge', &
                   vertViscTopOfEdge)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertViscTopOfCell', &
                   vertViscTopOfCell)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertDiffTopOfCell', &
                   vertDiffTopOfCell)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertAleTransportTop', &
                   vertAleTransportTop)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'vertTransportVelocityTop', &
                   vertTransportVelocityTop)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'barotropicForcing', &
                   barotropicForcing)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'kineticEnergyCell', &
                   kineticEnergyCell)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'transportVelocityX', &
                   transportVelocityX)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'transportVelocityY', &
                   transportVelocityY)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'transportVelocityZ', &
                   transportVelocityZ)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'transportVelocityZonal', &
                   transportVelocityZonal)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'transportVelocityMeridional', &
                   transportVelocityMeridional)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'tangentialVelocity', &
                   tangentialVelocity)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'montgomeryPotential', &
                   montgomeryPotential)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'BruntVaisalaFreqTop', &
                   BruntVaisalaFreqTop)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'bulkRichardsonNumber', &
                   bulkRichardsonNumber)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'bulkRichardsonNumberBuoy',&
                   bulkRichardsonNumberBuoy)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'bulkRichardsonNumberShear',&
                   bulkRichardsonNumberShear)

         ! set pointers for fields related forcing at ocean surface
         call mpas_pool_get_array_gpu(diagnosticsPool, 'surfaceBuoyancyForcing', &
                   surfaceBuoyancyForcing)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'surfaceFrictionVelocity', &
                   surfaceFrictionVelocity)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'barotropicThicknessFlux', &
                   barotropicThicknessFlux)

         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'normalTransportVelocity', &
                   normalTransportVelocity)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'normalVelocitySurfaceLayer', &
                   normalVelocitySurfaceLayer)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'penetrativeTemperatureFluxOBL', &
                   penetrativeTemperatureFluxOBL)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'normalizedRelativeVorticityCell', &
                   normalizedRelativeVorticityCell)

         call mpas_pool_get_array_gpu(diagnosticsPool, 'tracersSurfaceLayerValue', &
                   tracersSurfaceLayerValue)

         ! Set pointer for index variables
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_surfaceVelocityZonal', &
                   indexSurfaceVelocityZonalPtr)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_surfaceVelocityMeridional', &
                   indexSurfaceVelocityMeridionalPtr)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_vertNonLocalFluxTemp', &
                   index_vertNonLocalFluxTempPtr)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_SSHGradientMeridional', &
                   indexSSHGradientMeridionalPtr)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'indexSurfaceLayerDepth', &
                   indexSurfaceLayerDepth)
         call mpas_pool_get_array_gpu(diagnosticsPool, &
                  'indexBoundaryLayerDepth', &
                   indexBoundaryLayerDepth)
         call mpas_pool_get_dimension(diagnosticsPool, &
                  'index_SSHGradientZonal', &
                   indexSSHGradientZonalPtr)
         indexSurfaceVelocityZonal = indexSurfaceVelocityZonalPtr
         indexSurfaceVelocityMeridional = indexSurfaceVelocityMeridionalPtr
         index_vertNonLocalFluxTemp = index_vertNonLocalFluxTempPtr
         indexSSHGradientMeridional = indexSSHGradientMeridionalPtr
         indexSSHGradientZonal = indexSSHGradientZonalPtr

      end if

      call mpas_pool_get_array_gpu(diagnosticsPool, &
               "daysSinceStartOfSim", &
                daysSinceStartOfSim)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'Time', Time)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'Time_bnds', Time_bnds)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'xtime', xtime)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'simulationStartTime', &
                simulationStartTime)

      call mpas_pool_get_array_gpu(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'density', density)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'pressure', pressure)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'displacedDensity', &
                displacedDensity)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'potentialDensity', &
                potentialDensity)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'inSituThermalExpansionCoeff', &
                thermExpCoeff)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'inSituSalineContractionCoeff', &
                salineContractCoeff)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'pressureAdjustedSSH', &
                pressureAdjustedSSH)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'tracersSurfaceValue', &
                tracersSurfaceValue)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'surfaceFluxAttenuationCoefficient', &
                sfcFlxAttCoeff)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'surfaceFluxAttenuationCoefficientRunoff', &
                surfaceFluxAttenuationCoefficientRunoff)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'boundaryLayerDepth', &
                boundaryLayerDepth)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'boundaryLayerDepthSmooth', &
                boundaryLayerDepthSmooth)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'barotropicCoriolisTerm', &
                barotropicCoriolisTerm)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'betaEdge', betaEdge)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'bottomLayerShortwaveTemperatureFlux', &
                bottomLayerShortwaveTemperatureFlux)
      call mpas_pool_get_field(diagnosticsPool, 'rx1InitSmoothingMask', &
                smoothingMaskField, 1)
      call mpas_pool_get_field(diagnosticsPool, 'verticalStretch', &
                verticalStretchField, 1)
      call mpas_pool_get_array_gpu(diagnosticsPool, &
               'modifyLandIcePressureMask', &
                modifyLandIcePressureMask)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'indMLD', indMLD)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'dThreshMLD', dThreshMLD)
      call mpas_pool_get_array_gpu(diagnosticsPool, 'surfacePressure', &
                surfacePressure)

      call mpas_pool_get_array_gpu(diagnosticsPool, 'nhPressureCorrection', nhPressureCorrection)

      ! Semi-implicit Array Pointer retrievals
      if (trim(config_time_integrator) == 'split_implicit') then
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_r0', SIvec_r0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_r1', SIvec_r1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_v0', SIvec_v0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_w0', SIvec_w0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_w1', SIvec_w1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_t0', SIvec_t0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_q0', SIvec_q0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_s0', SIvec_s0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_s1', SIvec_s1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_t0', SIvec_t0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_t1', SIvec_t1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_y0', SIvec_y0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_z0', SIvec_z0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_z1', SIvec_z1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_r00', SIvec_r00)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_rh0', SIvec_rh0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_rh1', SIvec_rh1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_wh0', SIvec_wh0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_wh1', SIvec_wh1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_ph0', SIvec_ph0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_ph1', SIvec_ph1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_qh0', SIvec_qh0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_sh0', SIvec_sh0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_sh1', SIvec_sh1)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_zh0', SIvec_zh0)
         call mpas_pool_get_array_gpu(diagnosticsPool, 'SIvec_zh1', SIvec_zh1)
      end if

    end subroutine ocn_diagnostics_variables_init!}}}

!*******************************************************************************
!
!  ocn_diagnostics_variables_destroy
!
!> \brief Remove diagnostics variables from device
!> \author Rob Aulwes and Phil Jones
!> \date   14 Jan 2020
!> \details
!> This module removes the diagnostics variables from the device
!
!-------------------------------------------------------------------------------

   subroutine ocn_diagnostics_variables_destroy(err) !{{{

      ! Input variables

      ! Output variables

      integer, intent(out) :: &
         err                   ! returned error flag

      ! Local variables

      !***
      !*** end of preamble, begin code
      !***

      err = 0

      ! Nullify pointers
      if ( trim(config_land_ice_flux_mode) /= 'off' ) then
         nullify(topDrag,                         &
                 topDragMag,                      &
                 landIceBoundaryLayerTracers,     &
                 landIceTracerTransferVelocities, &
                 landIceFrictionVelocity)
      end if
      if ( trim(config_ocean_run_mode) == 'init' ) then
         nullify(rx1Edge,                  &
                 rx1Cell,                  &
                 rx1MaxEdge,               &
                 smoothingMask,            &
                 verticalStretch,          &
                 globalVerticalStretchMax, &
                 globalVerticalStretchMin)
      end if
      if ( trim(config_ocean_run_mode) == 'init' .or. &
                config_AM_debugDiagnostics_enable ) then
         nullify(rx1MaxCell,               &
                 globalRx1Max)
      end if
      if ( trim(config_ocean_run_mode) == 'forward' .or. &
           trim(config_ocean_run_mode) == 'analysis' ) then
         nullify(k33,                             &
                 gradSSH,                         &
                 gradSSHX,                        &
                 gradSSHY,                        &
                 gradSSHZ,                        &
                 viscosity,                       &
                 velocityX,                       &
                 velocityY,                       &
                 velocityZ,                       &
                 divergence,                      &
                 RiTopOfCell,                     &
                 RiTopOfEdge,                     &
                 SSHGradient,                     &
                 circulation,                     &
                 slopeTriadUp,                    &
                 gradSSHZonal,                    &
                 cGMphaseSpeed,                   &
                 velocityZonal,                   &
                 GMStreamFuncX,                   &
                 GMStreamFuncY,                   &
                 GMStreamFuncZ,                   &
                 layerThickEdgeFlux,              &
                 layerThickEdgeMean,              &
                 slopeTriadDown,                  &
                 normRelVortEdge,                 &
                 surfaceVelocity,                 &
                 vertVelocityTop,                 &
                 eddyVelocityX,                   &
                 eddyVelocityY,                   &
                 eddyVelocityZ,                   &
                 vertNonLocalFlux,                &
                 vertViscTopOfEdge,               &
                 gradSSHMeridional,               &
                 barotropicForcing,               &
                 vertViscTopOfCell,               &
                 vertDiffTopOfCell,               &
                 GMStreamFuncZonal,               &
                 relativeVorticity,               &
                 kineticEnergyCell,               &
                 normPlanetVortEdge,              &
                 velocityMeridional,              &
                 transportVelocityX,              &
                 transportVelocityY,              &
                 transportVelocityZ,              &
                 tangentialVelocity,              &
                 montgomeryPotential,             &
                 BruntVaisalaFreqTop,             &
                 vertAleTransportTop,             &
                 GMBolusVelocityZonal,            &
                 mleVelocityZonal,                &
                 bulkRichardsonNumber,            &
                 gmStreamFuncTopOfEdge,           &
                 relativeVorticityCell,           &
                 normalGMBolusVelocity,           &
                 normalMLEvelocity,               &
                 vertGMBolusVelocityTop,          &
                 vertMLEBolusVelocityTop,         &
                 transportVelocityZonal,          &
                 GMStreamFuncMeridional,          &
                 indexSurfaceLayerDepth,          &
                 surfaceBuoyancyForcing,          &
                 surfaceFrictionVelocity,         &
                 indexBoundaryLayerDepth,         &
                 barotropicThicknessFlux,         &
                 normalTransportVelocity,         &
                 vertTransportVelocityTop,        &
                 bulkRichardsonNumberBuoy,        &
                 tracersSurfaceLayerValue,        &
                 mleVelocityMeridional,           &
                 GMBolusVelocityMeridional,       &
                 bulkRichardsonNumberShear,       &
                 normalVelocitySurfaceLayer,      &
                 transportVelocityMeridional,     &
                 penetrativeTemperatureFluxOBL,   &
                 normalizedRelativeVorticityCell, &
                 nhPressureCorrection)
      end if
      if ( trim(config_ocean_run_mode) == 'forward' ) then
         nullify(wettingVelocity)
      end if
      if (config_use_GM.or.config_use_Redi) then
         nullify(RediKappa,         &
                 gmBolusKappa,      &
                 gmKappaScaling,    &
                 RediKappaScaling,  &
                 rediLimiterCount,  &
                 RediKappaSfcTaper, &
                 RediHorizontalTaper, &
                 gmHorizontalTaper, &
                 RossbyRadius)
      end if
      if (config_use_GM.or.config_submesoscale_enable) then
         nullify(gradDensityEddy)
      end if
      if ( config_compute_active_tracer_budgets ) then
         nullify(activeTracerHorMixTendency,              &
                 activeTracerVertMixTendency,             &
                 temperatureShortWaveTendency,            &
                 activeTracerNonLocalTendency,            &
                 activeTracerSurfaceFluxTendency,         &
                 activeTracerVerticalAdvectionTopFlux,    &
                 activeTracerVerticalAdvectionTendency,   &
                 activeTracerHorizontalAdvectionTendency, &
                 activeTracerHorizontalAdvectionEdgeFlux)
      end if
      if ( config_use_activeTracers_surface_restoring ) then
         nullify(salinitySurfaceRestoringTendency)
      end if
      if ( config_use_self_attraction_loading ) then
         nullify(ssh_sal)
      end if
      if ( trim(config_time_integrator) == 'split_implicit' ) then
         nullify(SIvec_r0,              &
                 SIvec_r00,             &
                 SIvec_r1,              &
                 SIvec_rh0,             &
                 SIvec_rh1,             &
                 SIvec_ph0,             &
                 SIvec_ph1,             &
                 SIvec_v0,              &
                 SIvec_v0,              &
                 SIvec_s0,              &
                 SIvec_s1,              &
                 SIvec_sh0,             &
                 SIvec_sh1,             &
                 SIvec_t0,              &
                 SIvec_t1,              &
                 SIvec_q0,              &
                 SIvec_q1,              &
                 SIvec_qh0,             &
                 SIvec_w0,              &
                 SIvec_w1,              &
                 SIvec_wh0,             &
                 SIvec_wh1,             &
                 SIvec_y0,              &
                 SIvec_z0,              &
                 SIvec_z1,              &
                 SIvec_zh0,             &
                 SIvec_zh1,             &
                 barotropicCoriolisTerm)
      end if

      nullify(zMid,                                    &
              zTop,                                    &
              xtime,                                   &
              indMLD,                                  &
              dThreshMLD,                              &
              density,                                 &
              betaEdge,                                &
              pressure,                                &
              thermExpCoeff,                           &
              sfcFlxAttCoeff,                          &
              surfacePressure,                         &
              potentialDensity,                        &
              displacedDensity,                        &
              boundaryLayerDepth,                      &
              salineContractCoeff,                     &
              pressureAdjustedSSH,                     &
              tracersSurfaceValue,                     &
              daysSinceStartOfSim,                     &
              Time,                                    &
              Time_bnds,                               &
              simulationStartTime,                     &
              boundaryLayerDepthSmooth,                &
              surfaceFluxAttenuationCoefficientRunoff)

    end subroutine ocn_diagnostics_variables_destroy!}}}

end module ocn_diagnostics_variables

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
