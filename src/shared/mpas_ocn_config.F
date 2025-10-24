! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_config
!
!> \brief MPAS ocean specific config
!> \details
!>  This module contains config specific to the ocean model.
!
!-----------------------------------------------------------------------

module ocn_config

   use mpaso_derived_types
   use mpaso_pool_routines
   use mpaso_kind_types

   implicit none
   public
   save

   #include "../inc/config_declare.inc"
   real (kind=RKIND) :: mod_lic_coeff_0, mod_lic_coeff_S, mod_lic_coeff_p, mod_lic_coeff_pS
   real (kind=RKIND) :: mod_oo_coeff_0, mod_oo_coeff_S, mod_oo_coeff_p, mod_oo_coeff_pS,mod_oo_az1_liq
   integer :: mod_jerlov_water_type

!$acc declare create(mod_lic_coeff_0,mod_lic_coeff_S,mod_lic_coeff_p,mod_lic_coeff_pS,&
!$acc                mod_oo_coeff_0,mod_oo_coeff_S,mod_oo_coeff_p,mod_oo_coeff_pS,mod_oo_az1_liq,&
!$acc                mod_jerlov_water_type)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_config_init
!
!> \brief   Initializes the ocean config
!> \details
!>  This routine sets up config for use in the ocean model.
!
!-----------------------------------------------------------------------
   subroutine ocn_config_init(configPool)!{{{
       type (mpas_pool_type), pointer :: configPool

       #include "../inc/config_get.inc"
       mod_lic_coeff_0 = config_land_ice_cavity_freezing_temperature_coeff_0
       mod_lic_coeff_S = config_land_ice_cavity_freezing_temperature_coeff_S
       mod_lic_coeff_p = config_land_ice_cavity_freezing_temperature_coeff_p
       mod_lic_coeff_pS = config_land_ice_cavity_freezing_temperature_coeff_pS
       mod_oo_coeff_0 = config_open_ocean_freezing_temperature_coeff_0
       mod_oo_coeff_S = config_open_ocean_freezing_temperature_coeff_S
       mod_oo_coeff_p = config_open_ocean_freezing_temperature_coeff_p
       mod_oo_coeff_pS = config_open_ocean_freezing_temperature_coeff_pS
       mod_oo_az1_liq = config_open_ocean_freezing_temperature_coeff_mushy_az1_liq
       mod_jerlov_water_type = config_jerlov_water_type

!$acc update device(mod_lic_coeff_0,mod_lic_coeff_S,mod_lic_coeff_p,mod_lic_coeff_pS,&
!$acc               mod_oo_coeff_0,mod_oo_coeff_S,mod_oo_coeff_p,mod_oo_coeff_pS,mod_oo_az1_liq,&
!$acc               mod_jerlov_water_type)

   end subroutine ocn_config_init!}}}

!***********************************************************************

end module ocn_config

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
