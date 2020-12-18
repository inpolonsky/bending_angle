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
! test for mo_rte_spherical_correction module
! Igor Polonsky
! ipolonsk@aer.com
! 
! reads test data file ref_data.dat
!  

program test_bending
   use mo_rte_spherical_correction 
   use mo_rte_kind, only: wp, wl
	implicit none


integer, parameter      :: nlev = 97 ! number of vertical levels
real(wp),parameter      :: RE = 6371.23 ! in [km]

real(wp), dimension(nlev)   :: temp, hum, RefIndex
!  no refraction
real(wp), dimension(nlev)   :: mu      ! solar zenith angle cosine
real(wp), dimension(nlev)   :: SZ      ! solar zenith angle
real(wp), dimension(nlev)   :: refSZ   ! reference solar zenith angle
!  with refraction
real(wp), dimension(nlev)   :: muRef      ! solar zenith angle cosine
real(wp), dimension(nlev)   :: SZref      ! solar zenith angle
real(wp), dimension(nlev)   :: refSZref   ! reference solar zenith angle
real(wp), dimension(nlev)   :: alt, p

real(wp), dimension(nlev)   :: refMu, refMuRef

logical(wl)             :: top_at_1     ! TOA atnosphere at index 1
logical(wl)             :: isError

integer                 :: jj
real(wp)                :: g0, zsurf, mu0, wavenumber, refRefIndex
!
! read meterological profile 
open(100, file='ref_data.dat', action='read', status='old')
!  skip comments
do jj=1,5
	read(100,*)
enddo

do jj=1,nlev
	read(100,*) alt(jj), p(jj), temp(jj), hum(jj), szRef(jj), refSZref(jj)
enddo
!  convert to [g/g]
hum=hum*1e-3
close(100)

! initial condition
!  solar zenith anle is 90 at surface
!  wavelength is 500 nm

mu0 = cos(acos(-1d0)/180_wp*90_wp)
wavenumber = 1e4/0.5

top_at_1 = p(1) < p(2)
!  compute index of refraction

call indexRefraction(nlev, temp, p, hum, wavenumber, RefIndex)

if (top_at_1)  then
   refRefIndex = RefIndex(nlev)
   zsurf = alt(nlev)
else
   refRefIndex = RefIndex(1)
   zsurf = alt(1)
endif

! if only pressure, temperature and water vapor 
! the altitude can be computed using cmpalt
! uncomment the following line

! call cmpalt(nlev, p, temp, hum, zsurf, alt, top_at_1, 45.0_wp)

! compute solar zenith cosine without effect of refraction
call spherical_angles(nlev, zsurf, RE, alt, mu0, mu) 


! compute solar zenith cosine including the effect of refraction
call spherical_angles(nlev, zsurf, RE, alt, mu0, muRef, refRefIndex, RefIndex) 

! print out results
print *, ' altitude    SZA cosine    SZA cosine        SZA        SZA      Refractivity'
print *, '             no refrac.    with refrac.     no ref.   with ref.   (n-1)*1e6'
print *, '     [km]                                   [deg]      [deg]     [ppm]'
 
do jj=1, nlev
	print '(4x, f6.3, 2(2x,E12.4), 3(4x, F7.2))', alt(jj), mu(jj), muRef(jj), 180_wp/acos(-1.)*acos(mu(jj)), &
	     & 180_wp/acos(-1.)*acos(muRef(jj)), (RefIndex(jj) - 1.)*1e6
enddo
print *,''
print *, 'Compare with reference:'

!  compare with reference vaules
SZ = 180_wp/acos(-1.)*acos(mu) 
isError = any(abs(SZ-szRef) > 1e-6)

if (isError) then
	print *, 'bias in solar zenith angle (no refration)'
   print *, "level     angle       ref.angle        bias"
	do jj =1, nlev
      if (abs(SZ(jj)-szRef(jj)) > 1e-6) &
           print "(I6, 3(F12.6, 2x)  )", jj, SZ(jj), szRef(jj), SZ(jj) - szRef(jj)
	enddo
else
	print *, '  solar zenith angle (no refration)  : PASSED'
endif

SZref = 180_wp/acos(-1.)*acos(muRef)

if (any(abs(SZref-refSZref) > 1e-6)) then
	print *, 'bias in solar zenith angle (with refration)'
   print *, "level     angle       ref.angle        bias"
	do jj =1, nlev
		if (abs(SZref(jj)-refSZref(jj)) > 1e-6) print "(I6, 3(F12.6, 2x)  )", jj, SZref(jj), &
                      refSZref(jj), SZref(jj) - refSZref(jj)
	enddo
else
	print *, '  solar zenith angle (with refration): PASSED'
endif
end program test_bending