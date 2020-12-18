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


  public :: spherical_angles, cmpalt
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

  
    !**************************************************************         
    !     AUTHOR: TONY CLOUGH, JENNIFER DELAMERE, JOHN WARDEN               
    !             JANUARY 2001                                              
    !     PROGRAM TO CALCULATE ALTITUDE LEVEL (alt) GIVEN                  
    !     PRESSURE (p), TEMPERATURE (T) AND THE NUMBER DENSITY            
    !     OF WATER VAPOR (q) USING THE HYDROSTATIC EQUATION              
    !                                                                       
    !     INPUT:                                                            
    !      A) PRESSURE (MBAR)                                               
    !      B) TEMPERATURE (KELVIN)                                          
    !      C) NUMBER DENSITY OF WATER VAPOR                                 
    !                                                                       
    !     OUTPUT:                                                           
    !      A) ALTITUDE (KM)                                                 
    !      IDEAL GAS LAW: P.E.CIDDOR (1996), Refractive index of 
    !      air: New equations for the visible and near infrared,
    !      Applied Optics, 35(9), 1566-1573.
    !    
    !     changed layer numbering convention from LBLRT to OSSTRAN
    !**************************************************************         
  subroutine cmpalt(nlev, p, T, q, surfAlt, alt, top_at_1, lat)
    !-----------------------------------------------------------------------                                                                        
    integer,                               intent(in   ) :: nlev         ! Number of levels
    real(wp),                              intent(in   ) :: surfAlt      ! surface altitude  [km]
    logical(wl),                           intent(in   ) :: top_at_1     ! TOA atnosphere at index 1
    real(wp), dimension(nlev),             intent(in   ) :: p            ! Pressure profile [mbar]
    real(wp), dimension(nlev),             intent(in   ) :: T            ! Temperature profile [K]
    real(wp), dimension(nlev),             intent(in   ) :: q            !  Water vapor profile [g/g]

    real(wp), optional,                    intent(in   ) :: lat          ! latitude
    real(wp), dimension(nlev),             intent(inout) :: alt          ! level altitude [km]
    
    real(wp), PARAMETER :: CA0 = 1.58123D-6, CA1 = -2.9331D-8, CA2 = 1.1043D-10
    real(wp), PARAMETER :: CB0 = 5.707D-6,   CB1 = -2.051D-8
    real(wp), PARAMETER :: CC0 = 1.9898D-4,  CC1 = -2.376D-6
    real(wp), PARAMETER :: CD  = 1.83D-11,   CE  = -0.0765D-8

    integer :: I, J
    real(wp) :: COMP_FACTOR(nlev)
    real(wp) :: G0, GAVE
    real(wp) :: Y
    real(wp) :: CHI0, DCHI
    real(wp) :: T0, DT
    real(wp) :: CHIM
    real(wp) :: C1, C2, C3
    real(wp) :: A, B, ALPHA
    real(wp) :: XINT_TOT


    ! CALCULATE GRAVITY AT REFERENCE LATITUDE AT SURFACE                    

    G0 = GRAV_CONST(lat)

    ! CALCULATE THE NUMBER DENSITY OF TOTAL AIR MOLECULES [MOLEC/CM^3]      
    ! CALCULATE THE COMPRESSIBILITY FACTOR (COMP_FAC) FOR THE               
    ! IDEAL GAS LAW                                                         
    DO J = 1, nlev
        DT = T(J) - TZERO
        CHIM = q(J)
        COMP_FACTOR(J) = 1_wp - &
                        (p(J) * 100_wp/T(J)) * &
                        (CA0 + CA1 * DT + CA2 * DT**2 + (CB0 + CB1 * DT) * CHIM + (CC0 + CC1 * DT) * CHIM**2) + &
                        (CD + CE * CHIM**2)*(p(J) * 100_wp/T(J))**2
    ENDDO


    ! CONVERT REFERENCE ALTITUDE TO METERS                                  
    if ( top_at_1 ) then
        alt(nlev) = surfAlt
        DO I = nlev, 2, -1
            GAVE = G0 * (RE/(RE + alt(I)))**2
            Y = LOG(p(I - 1)/p(I))
            IF (Y .NE. 0.0_wp) THEN
                B = 1_wp + q(I)
                A = (q(I - 1) - q(I))
                ALPHA = A/B
                ! TODO should be added a p-log interpolation to 
                ! make it more accurate
                !
                ! IF (ABS(ALPHA) .GE. 0.01d0) THEN
                !     PRINT*, 'LAYER TOO THICK'
                !     STOP
                ! END IF
                CHI0 = q(I)/XMASS_RATIO
                DCHI = A/XMASS_RATIO

                T0 = T(I)
                DT = (T(I - 1) - T(I))

                C1 = T0 + T0 * CHI0
                C2 = T0 * DCHI + DT * CHI0 + DT - C1 * ALPHA
                C3 = DT * DCHI - C2 * ALPHA

                XINT_TOT = C1  + 0.5_wp * C2  + 0.3333_wp * C3
                XINT_TOT = -Y*XINT_TOT * (GASCON * 1.0D-7)/(AIRMWT * GAVE * B)

                alt(I - 1) = alt(I) + XINT_TOT * COMP_FACTOR(I)
            ELSE
                alt(I - 1) = alt(I)
            END IF
        END DO
      else
        alt(1) = surfAlt
        DO I = 1, nlev -1
            GAVE = G0 * (RE/(RE + alt(I)))**2
            Y = LOG(p(I + 1)/p(I))
            IF (Y .NE. 0.0_wp) THEN
                B = 1_wp + q(I)
                A = (q(I + 1) - q(I))
                ALPHA = A/B
                ! TODO should be added a p-log interpolation to 
                ! make it more accurate
                ! IF (ABS(ALPHA) .GE. 0.01d0) THEN
                !     PRINT*, 'LAYER TOO THICK'
                !     STOP
                ! END IF

                CHI0 = q(I)/XMASS_RATIO
                DCHI = A/XMASS_RATIO

                T0 = T(I)
                DT = (T(I + 1) - T(I))

                C1 = T0 + T0 * CHI0
                C2 = T0 * DCHI + DT * CHI0 + DT - C1 * ALPHA
                C3 = DT * DCHI - C2 * ALPHA

                XINT_TOT = C1  + 0.5_wp * C2 + 0.3333_wp * C3
                XINT_TOT = -Y*XINT_TOT * (GASCON * 1.0D-7)/(AIRMWT * GAVE * B)
                alt(I + 1) = alt(I) + XINT_TOT * COMP_FACTOR(I)
            ELSE
                alt(I + 1) = alt(I)
            END IF
        END DO
      endif
    RETURN

end subroutine cmpalt
    
!
! calculate gravitational acceleration constant
!        

function GRAV_CONST(LATITUDE)
        REAL(wp), INTENT(IN), OPTIONAL :: LATITUDE ! in degrees
        REAL(wp) :: GRAV_CONST ! in meters/s^2
        ! local
        REAL(wp) lat
        !         Latitude for which gravitational constant desired
        IF (PRESENT(LATITUDE)) THEN
            lat = LATITUDE
        ELSE
            lat = 45.d0
        END IF
        !         Gravitational constant for Earth in meters/s^2
        GRAV_CONST = 9.80665d0 - 0.02586d0 * COS(2.0d0 * PI * lat/180.0d0)

  end function GRAV_CONST

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
