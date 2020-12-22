! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Account for effects of refraction and geometry related to the spherical shell Earth model
! Igor Polonsky
! ipolonsk@aer.com
! 
! service functions: 
! 1. cmpalt        computational altitude grid corresponding to 
!                    pressure, temperature and water vapor   
! 2  indexRefraction - computed air index of refraction
!
! Input:
!  
! -------------------------------------------------------------------------------------------------
module mo_rte_spherical_correction
  use mo_rte_kind, only: wp, wl
  implicit none
  private


  public :: spherical_angles
  public :: indexRefraction

  real(wp), parameter, public :: BADANGLE=-9999.0

  ! public :: lw_spherical_angles

  real(wp), parameter :: pi = acos(-1._wp)
  real(wp), parameter :: PZERO  = 1013.25 ! [hPa]
  real(wp), parameter :: TZERO = 273.15D0 ![K]
  real(wp), parameter :: AIRMWT = 28.964D0 ! air molecular weight (grams/mole)
  real(wp), parameter :: H2OMWT = 18.015D0 ! air molecular weight (grams/mole)
  real(wp), parameter :: XMASS_RATIO = H2OMWT/AIRMWT

  real(wp), PARAMETER :: GASCON = 8.314472D+07
  real(wp), parameter :: RE = 6371.23D0 ! km
!                    RE        radius of earth (km)
!                     a)  MODEL 0,2,3,6    RE = 6371.23 km
!                     b)        1          RE = 6378.39 km
!                     c)        4,5        RE = 6356.91 km

  interface spherical_angles
    module procedure spherical_angles_simple
    module procedure spherical_angles_refr
  end interface spherical_angles


contains

  ! -------------------------------------------------------------------------------------------------
  !
  ! Top-level kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! pseudo-spherical angles: convert solar zenith angle at surface to layer average angle
  !  
  !
  ! ---------------------------------------------------------------
  subroutine spherical_angles_simple(nlev, refAlt, ErthRad, alt, refMu, mu)
    integer,                               intent(in   ) :: nlev         ! Number of levels
    real(wp),                              intent(in   ) :: refAlt       ! reference altitude at which solar zenith angle is provided [km]
    real(wp), dimension(nlev),             intent(in   ) :: alt          ! level altitude grid [km]
    real(wp),                              intent(in   ) :: ErthRad      ! The Earth curvature radius [km]
    real(wp),                              intent(in   ) :: refMu        ! cosine of the solar zenith angle at reference altitude
    real(wp), dimension(nlev),             intent(inout) :: mu           ! cosine of the solar zenith angle at level grid

    integer             :: ilev
    real(wp)            :: impact            
    ! ------------------------------------
    impact = (ErthRad + refAlt)*sqrt(1_wp - refMu**2)
    ! ------------------------------------
    do ilev=1, nlev
      mu(ilev) = sqrt(1_wp - ( impact/(ErthRad + alt(ilev)) )**2 )
    end do
  end subroutine spherical_angles_simple

  ! -------------------------------------------------------------------------------------------------
  !
  ! pseudo-spherical angles: convert solar zenith angle at reference altitude to level solar zenith angle cosine 
  ! note: if sine at any altitude > 1 (super refraction) BADANGLE returned
  ! presence of BADANGLE at any level means that the solar radiance cannot penetrate this altitiute and deeper
  !
  ! ---------------------------------------------------------------
  subroutine spherical_angles_refr(nlev, refAlt, ErthRad, alt, refMu, mu, refRefrIndex, refrIndex) 
    integer,                               intent(in   ) :: nlev         ! Number of levels
    real(wp),                              intent(in   ) :: refAlt       ! reference altitude at which solar zenith angle is provided [km]
    real(wp), dimension(nlev),             intent(in   ) :: alt          ! level altitude grid [km]
    real(wp),                              intent(in   ) :: ErthRad      ! The Earth curvature radius [km]
    real(wp),                              intent(in   ) :: refMu        ! cosine of the solar zenith angle at reference altitude
    real(wp),                              intent(in   ) :: refRefrIndex ! air index of refraction at reference altitude
    real(wp), dimension(nlev),             intent(in   ) :: refrIndex    ! air index of refraction at altitude grid
    real(wp), dimension(nlev),             intent(inout) :: mu           ! cosine of the solar zenith angle at altitude grid

    integer             :: ilev
    real(wp)            :: impact, sn            
    ! ------------------------------------
    impact = (ErthRad + refAlt)*refRefrIndex*sqrt(1_wp - refMu**2)
    do ilev=1, nlev
      sn = impact/(ErthRad + alt(ilev))/refrIndex(ilev)
      if (sn > 1) then
         mu(ilev) = BADANGLE
      else   
         mu(ilev) = sqrt(1_wp - ( impact/(ErthRad + alt(ilev))/refrIndex(ilev) )**2 )
      endif  
    end do
  end subroutine spherical_angles_refr

  
  
    subroutine indexRefraction(nlev, TM, PM, q, wavenumber, RefIndex)
      !  compute air refraction index for range >0.23 mkm and IR
      !  fails for MV and radio wave
      ! 
      !  Bengt Edlén (1966): The Refractive Index of Air, Metrologia, 2, 71-80
      ! 
       integer,                               intent(in   ) :: nlev         ! Number of levels
       real(wp), dimension(nlev),             intent(in   ) :: PM           ! Pressure profile [mbar]
       real(wp), dimension(nlev),             intent(in   ) :: TM           ! Temperature profile [K]
       real(wp), dimension(nlev),             intent(in   ) :: q            !  Water vapor profile [g/g]
       real(wp),                              intent(in   ) :: wavenumber   !  wavenumber  [cm^{-1}]

       real(wp), dimension(nlev),             intent(inout) :: RefIndex     ! air refraction index

       integer  :: ilev
       real(wp) :: c1, c2
       real(wp), parameter  :: T15  = TZERO + 15_wp![K]

      do ilev=1, nlev
   !
   !   Approximation to refraction index (from LOWTRAN6)
   !   Kneizys, F. X. Chetwynd, J. H., Jr. Clough, S. A. Shettle, E. P. Abreu, L. W. (1983)
   !  Atmospheric Transmittance/Radiance: Computer Code LOWTRAN 6
   !  Environmental research papers
   !  AIR FORCE GEOPHYSICS LAB HANSCOM AFB MA
   !  https://apps.dtic.mil/sti/citations/ADA137786

       c1 = (83.42_wp+(185.08_wp/(1_wp-(wavenumber/1.14E+5)**2))+&
            (4.11_wp/(1_wp-(wavenumber/6.24E+4)**2)))*(T15/TM(ilev)) 

       c2 = (43.49_wp -(wavenumber/1.7E+4)**2)*(q(ilev)/(XMASS_RATIO + q(ilev)))
       RefIndex(ilev) = 1_wp + 1d-6*(c1 - c2)*(PM(ilev)/PZERO)
      enddo
    end subroutine indexRefraction

end module mo_rte_spherical_correction
