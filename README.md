bending_angle repository is provided support for RRTMGP package 
(https://github.com/RobertPincus/rte-rrtmgp)
that provides solar zenith angle computed in spherical shell model
for the plane-parallel codes
the algorithm are discussed in 
Gallery W., F. X. Kneizys, S.A. Clough (1983): 
Air Mass Computer Program for Atmospheric Transmittance/Radiance Calculation: FSCATM 
AFGL-TR-83-0065, ENVIRONMENTAL RESEARCH PAPERS, NO. 828, 
AIR FORCE GEOPHYSICS LABORATORY
(https://apps.dtic.mil/sti/pdfs/ADA132108.pdf)


The package consists of: 
   mo_rte_spherical_correction.F90 - module that provides all functionality
   test_spherical_correction.f90 - example demonstrating the usage
   ref_data.dat - ASCII that provides data for the above test


mo_rte_spherical_correction.F90 provides two ways to compute spherical correction:
with and without use of refraction, both are rely on the call of subroutine 
spherical_angles with a different set of parameters

1. spherical_angles(nlev, refAlt, ErthRad, alt, refMu, mu)
2. spherical_angles(nlev, refAlt, ErthRad, alt, refMu, mu, refRefrIndex, refrIndex)

    integer,                  :: nlev         - Number of levels
    real(wp),                 :: refAlt       - reference altitude at which solar zenith angle is provided [km]
    real(wp), dimension(nlev) :: alt          - level altitude grid [km]
    real(wp),                 :: ErthRad      - The Earth curvature radius [km]
    real(wp),                 :: refMu        - cosine of the solar zenith angle at reference altitude
    real(wp),                 :: refRefrIndex - air index of refraction at reference altitude
    real(wp), dimension(nlev) :: refrIndex    - air index of refraction at altitude grid
    real(wp), dimension(nlev) :: mu           - cosine of the solar zenith angle at altitude grid


The user is expected to provide the solar zenith angle (and the air refractive index if effect of 
refraction is included) at any fixed altitude, the subroutine would provide the solar zenith cosine 
at the user specified altitude grid.

The module also provides two auxiliary subroutine:

1. subroutine cmpalt(nlev, p, T, q, surfAlt, alt, top_at_1, lat) 
    integer,                  :: nlev        - Number of levels
    real(wp),                 :: surfAlt     - surface altitude  [km]
    logical(wl),              :: top_at_1    - TOA atnosphere at index 1
    real(wp), dimension(nlev) :: p           - Pressure profile [mbar]
    real(wp), dimension(nlev) :: T           - Temperature profile [K]
    real(wp), dimension(nlev) :: q           -  Water vapor profile [g/g]

    real(wp), optional,       :: lat         - latitude
    real(wp), dimension(nlev) :: alt         - level altitude [km]
    
    computes altitude based on the user provided pressure, temperature and water vapor mixing ratio profiles.
    The altitude corresponding the highest pressure as well as geographical latitude are expected.

2.  subroutine indexRefraction(nlev, TM, PM, q, wavenumber, RefIndex)
       integer,                  :: nlev         - Number of levels
       real(wp), dimension(nlev) :: PM           - Pressure profile [mbar]
       real(wp), dimension(nlev) :: TM           - Temperature profile [K]
       real(wp), dimension(nlev) :: q            -  Water vapor profile [g/g]
       real(wp),                 :: wavenumber   -  wavenumber  [cm^{-1}]

       real(wp), dimension(nlev) :: RefIndex     - air refraction index

   computes the air refraction index for range >0.23 mkm and IR  fails for MV and radio wave implemented based on LOWTRAN 6 approach

 Kneizys, F. X. Chetwynd, J. H., Jr. Clough, S. A. Shettle, E. P. Abreu, L. W. (1983)
 Atmospheric Transmittance/Radiance: Computer Code LOWTRAN 6
 Environmental research papers
 AIR FORCE GEOPHYSICS LAB HANSCOM AFB MA
 https://apps.dtic.mil/sti/citations/ADA137786