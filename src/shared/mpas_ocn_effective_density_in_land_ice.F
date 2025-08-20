! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_effective_density_in_land_ice
!
!> \brief MPAS ocean effective density in land ice
!> \author Xylar Asay-Davis
!> \date   10/03/2015
!> \details
!>  This module contains routines for computing the effective seawater
!>  density in land ice using Arhimedes' principle.
!
!-----------------------------------------------------------------------

module ocn_effective_density_in_land_ice

   use mpaso_constants
   use mpaso_kind_types
   use mpaso_derived_types
   use mpaso_pool_routines

   use ocn_constants
   use ocn_config
   use ocn_mesh

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

   public :: ocn_effective_density_in_land_ice_update, &
             ocn_effective_density_in_land_ice_update_gpu

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_effective_density_in_land_ice_update
!
!> \brief updates effective density in land ice
!> \author Xylar Asay-Davis
!> \date   10/03/2015
!> \details
!>  This routine updates the value of the effective seawater density
!>  displaced by land ice, based on Archimedes' principle.  The
!>  effective density is smoothed and extrapolated by averaging with
!>  nearest neighbors (cellsOnCell).
!
!-----------------------------------------------------------------------

   subroutine ocn_effective_density_in_land_ice_update(forcingPool, &
                                                    statePool, ierr)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         forcingPool !< [in] Forcing information

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool !< [inout] state information

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: ierr !< [out] Error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, iCell,  &! loop indices for neighbor loop and cell loop
         cell2       ! neighbor cell address

      real (kind=RKIND) :: &
         weightSum   ! local sum of weights for masked points

      real (kind=RKIND), dimension(:), pointer :: &
         ssh,                 &! sea surface height
         landIcePressure,     &! land ice pressure
         effectiveDensityCur, &! effective density at current time
         effectiveDensityNew   ! effective density at new time level
                                                  
      integer, dimension(:), pointer :: &
         landIceMask           ! mask for land ice presence on cell

      ! Scratch Arrays
      ! effectiveDensityScratch: effective seawater density in land
      ! ice before horizontal averaging units: kg m^{-3}
      real (kind=RKIND), dimension(:), allocatable :: &
         effectiveDensityScratch

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Set default return code and exit if not needed
      ierr = 0
      if ( (trim(config_land_ice_flux_mode) /= 'coupled') ) return

      ! Retrieve forcing and state variables
      call mpas_pool_get_array(forcingPool, 'landIceMask', &
                                             landIceMask)
      call mpas_pool_get_array(forcingPool, 'landIcePressure', &
                                             landIcePressure)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 2)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', &
                                           effectiveDensityCur, 1)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', &
                                           effectiveDensityNew, 2)

      allocate(effectiveDensityScratch(nCellsAll))

      ! Compute effective density in each cell
      ! TODO: should only apply to floating land ice,
      !       once wetting/drying is supported
      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCellsAll
         if (landIceMask(iCell) == 1) then
            ! land ice is present to update the effective density
            effectiveDensityScratch(iCell) = -landIcePressure(iCell)/ &
                                              (ssh(iCell)*gravity)
         else
            ! we copy the previous effective density
            effectiveDensityScratch(iCell) = effectiveDensityCur(iCell)
         end if
      end do
      !$omp end do

      ! smooth/extrapolate density by averaging with nearest neighbors
      !$omp do schedule(runtime) private(weightSum, i, cell2)
      do iCell = 1, nCellsAll
         weightSum = 1.0_RKIND
         effectiveDensityNew(iCell) = effectiveDensityScratch(iCell)
         do i = 1, nEdgesOnCell(iCell)
            cell2 = cellsOnCell(i,iCell)
            effectiveDensityNew(iCell) = effectiveDensityNew(iCell) &
                                       + cellMask(1,cell2)* &
                                         effectiveDensityScratch(cell2)
            weightSum = weightSum + cellMask(1,cell2)
         end do
         effectiveDensityNew(iCell) = effectiveDensityNew(iCell)/ &
                                      weightSum
      end do
      !$omp end do
      !$omp end parallel

      deallocate(effectiveDensityScratch)

   !--------------------------------------------------------------------

   end subroutine ocn_effective_density_in_land_ice_update !}}}

   subroutine ocn_effective_density_in_land_ice_update_gpu(forcingPool, &
                                                    statePool, ierr)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         forcingPool !< [in] Forcing information

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool !< [inout] state information

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: ierr !< [out] Error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, iCell,  &! loop indices for neighbor loop and cell loop
         cell2       ! neighbor cell address

      real (kind=RKIND) :: &
         weightSum   ! local sum of weights for masked points

      real (kind=RKIND), dimension(:), pointer :: &
         ssh,                 &! sea surface height
         landIcePressure,     &! land ice pressure
         effectiveDensityCur, &! effective density at current time
         effectiveDensityNew   ! effective density at new time level
                                                  
      integer, dimension(:), pointer :: &
         landIceMask           ! mask for land ice presence on cell

      ! Scratch Arrays
      ! effectiveDensityScratch: effective seawater density in land
      ! ice before horizontal averaging units: kg m^{-3}
      real (kind=RKIND), dimension(:), allocatable :: &
         effectiveDensityScratch

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Set default return code and exit if not needed
      ierr = 0
      if ( (trim(config_land_ice_flux_mode) /= 'coupled') ) return
!NV!  !!! NOT EXECUTING !!!
      write(*,*) "config_land_ice_flux_mode NOT EXECUTING"

      ! Retrieve forcing and state variables
      call mpas_pool_get_array(forcingPool, 'landIceMask', &
                                             landIceMask)
      call mpas_pool_get_array(forcingPool, 'landIcePressure', &
                                             landIcePressure)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 2)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', &
                                           effectiveDensityCur, 1)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', &
                                           effectiveDensityNew, 2)

      allocate(effectiveDensityScratch(nCellsAll))

      ! Compute effective density in each cell
      ! TODO: should only apply to floating land ice,
      !       once wetting/drying is supported
      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCellsAll
         if (landIceMask(iCell) == 1) then
            ! land ice is present to update the effective density
            effectiveDensityScratch(iCell) = -landIcePressure(iCell)/ &
                                              (ssh(iCell)*gravity)
         else
            ! we copy the previous effective density
            effectiveDensityScratch(iCell) = effectiveDensityCur(iCell)
         end if
      end do
      !$omp end do

      ! smooth/extrapolate density by averaging with nearest neighbors
      !$omp do schedule(runtime) private(weightSum, i, cell2)
      do iCell = 1, nCellsAll
         weightSum = 1.0_RKIND
         effectiveDensityNew(iCell) = effectiveDensityScratch(iCell)
         do i = 1, nEdgesOnCell(iCell)
            cell2 = cellsOnCell(i,iCell)
            effectiveDensityNew(iCell) = effectiveDensityNew(iCell) &
                                       + cellMask(1,cell2)* &
                                         effectiveDensityScratch(cell2)
            weightSum = weightSum + cellMask(1,cell2)
         end do
         effectiveDensityNew(iCell) = effectiveDensityNew(iCell)/ &
                                      weightSum
      end do
      !$omp end do
      !$omp end parallel

      deallocate(effectiveDensityScratch)

   !--------------------------------------------------------------------

   end subroutine ocn_effective_density_in_land_ice_update_gpu !}}}


!***********************************************************************

end module ocn_effective_density_in_land_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
