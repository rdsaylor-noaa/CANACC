! ======================================================================================================================!
!                                                                                                                      !
!     Program:      CANACC                                                                                             !
!                   Canopy physics model for ACCESS                                                                    !
!                                                                                                                      !
!     Version:      3.0.0                                                                                              !
!                   SEB; Flux upper BCs; Coweeta Kv; Measured LWdn_sky; Enhanced diagnostics                           !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Contact:      Rick D. Saylor, PhD                                                                                !
!                   Physical Scientist                                                                                 !
!                   U. S. Department of Commerce                                                                       !
!                   National Oceanic and Atmospheric Administration                                                    !
!                   Air Resources Laboratory                                                                           !
!                   Atmospheric Turbulence and Diffusion Division                                                      !
!                   456 S. Illinois Ave                                                                                !
!                   Oak Ridge, TN 37830                                                                                !
!                   email: Rick.Saylor@noaa.gov                                                                        !
!                                                                                                                      !
!**********************************************************************************************************************!
!                   NOTE: See Legal Notice in Main.f90                                                                 !
!**********************************************************************************************************************!
!                                                                                                                      ! 
!     Module:       GlobalData                                                                                         !
!                                                                                                                      ! 
!     Description:  Contains all global variables used by the model                                                    !
!                                                                                                                      !
!======================================================================================================================!
module GlobalData
  implicit none
  public
  save

  ! double precision kind parameter
  integer, parameter :: dp=kind(0.d0)
  ! 4 bytes integer
  integer, parameter :: i4=kind(1)

  ! tstart - time of simulation start
  real(kind=dp) :: tstart
  ! tend = time of simulation end
  real(kind=dp) :: tend
  ! dt - integration time step (seconds) 
  real(kind=dp) :: dt
  ! dthalf - half of integration time step (seconds)
  real(kind=dp) :: dthalf
  ! dtout - time step for output
  !       MUST be even multiple of dt
  real(kind=dp) :: dtout
  ! ntout - total number of output time steps for simulation
  integer(kind=i4) :: ntout
  ! zadeg - zenith angle at current meteorological data time step (degrees)
  real(kind=dp) :: zadeg
  ! zarad - zenith angle at current meteorological data time step (radians)
  real(kind=dp) :: zarad

  ! t - current simulation time
  real(kind=dp) :: t
  ! nt - current output time step number
  integer(kind=i4) :: nt
  ! sdatetime - current simulation time as a datetime string
  character(len=19) :: sdatetime
  ! zcurrent - current z value (only valid within chemistry integration!) (m)
  real(kind=dp) :: zcurrent
  ! zpt - current vertical grid point (only valid within chemistry integration!)
  integer(kind=i4) :: zpt

  ! number of grid nodes total
  integer(kind=i4), parameter :: npts = 87
  ! number of grid nodes in canopy
  integer(kind=i4)            :: ncnpy
  ! within canopy grid resolution (cm)
  real(kind=dp)               :: dzhc

  ! column bottom and top (cm)
  real(kind=dp)               :: z0
  real(kind=dp)               :: zi

  ! canopy height (m)
  real(kind=dp)  :: hc
  ! canopy height (cm)
  real(kind=dp)  :: hccm
  ! stretched grid parameter ()
  real(kind=dp)  :: alfa

  ! roughness length for momentum (cm)
  real(kind=dp)  :: z0m
  ! roughness length for heat (cm)
  real(kind=dp)  :: z0h
  ! zero-plane displacement height (cm)
  real(kind=dp)  :: d
  ! measurement reference height (cm)
  real(kind=dp)  :: zref

  ! vertical grid definition (cm)
  real(kind=dp), dimension(npts) :: z

  ! layer leaf area density of canopy (cm2 leaf/cm3)
  real(kind=dp), dimension(npts) :: lad
  ! cumulative leaf area index of canopy (cm2 leaf/cm2)
  real(kind=dp), dimension(npts) :: clai
  ! layer leaf area index of canopy (cm2 leaf/cm2)
  real(kind=dp), dimension(npts) :: lai
  ! total leaf area index of canopy (cm2 leaf/cm2)
  real(kind=dp)                  :: laitot


  ! ninteg, number of integrated species
  integer(kind=i4), parameter    :: ninteg=1

  ! ntotal, total number of integrated species
  integer(kind=i4), parameter    :: ntotal=ninteg

  ! q - emissions rate of all species at the current time (molecules/cm3-s)
  real(kind=dp), dimension(npts, ninteg) :: q
  ! gp - canopy compensation point :: (molecules/cm3)
  real(kind=dp), dimension(npts, ninteg) :: gp
  ! gpst - stomatal compensation point :: (molecules/cm3)
  real(kind=dp), dimension(npts, ninteg) :: gpst
  ! gpg - soil compensation point :: (molecules/cm3)
  real(kind=dp), dimension(ninteg) :: gpg

  ! vd - dry deposition exchange coefficient (cm/s)
  real(kind=dp), dimension(npts, ninteg) :: vd
  ! rb - leaf boundary layer resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rb
  ! rs - stomatal resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rs
  ! rw - cuticular resistance (s/cm)
  real(kind=dp), dimension(npts, ninteg) :: rw

  ! rakv - aerodynamic resistance calculated from Kv profile (s/cm)
  real(kind=dp), allocatable :: rakv(:)
  ! ra - aerodynamic resistance calculated from Viney et al. (s/cm)
  real(kind=dp), allocatable :: ra(:)
  ! Monin-Obukhov stability parameter, (zref-d)/L ()
  real(kind=dp), allocatable :: zolout(:)
  ! gaero - aerodynamic conductance (mol/m2-s)
  real(kind=dp), allocatable :: gaero(:)

  ! Soil exchange data
  ! vs - soil exchange coefficients (cm/s)
  real(kind=dp), dimension(ninteg) :: vs
  ! vsh2o - soil exchange coefficient for water vapor (cm/s)
  real(kind=dp)                    :: vsh2o
  ! qsoil - effective soil molar humidity (mol/cm3)
  real(kind=dp)                    :: qsoil
  ! ksoil - thermal conductivity of soil (W/m-K)
  real(kind=dp)                    :: ksoil
  ! tsoilk - soil/litter temperature (K)
  real(kind=dp)                    :: tsoilk
  ! tk0 - lower boundary condition for air temperature
  !        derived from the surface energy balance (K)
  real(kind=dp)                    :: tk0
  ! dtdzsoil - temperature gradient in soil surface (K/m)
  real(kind=dp)                    :: dtdzsoil
  ! gbg - ground boundary layer conductance (mol/m2-s)
  real(kind=dp)                    :: gbg
  ! effrhsoil - effective soil fractional relative humidity
  real(kind=dp)                    :: effrhsoil
  ! rbg - ground boundary layer resistance (s/cm)
  real(kind=dp)                    :: rbg
  ! rsoil - resistance to diffusion thru soil pore space (s/cm)
  real(kind=dp)                    :: rsoil
  ! nsoilep - number of species with non-zero soil emission potentials
  integer(kind=i4) :: nsoilep
  ! gammaso - soil emission potential (molar ratio) = [NH4+]/[H+] for ammonia
  real(kind=dp), dimension(ninteg) :: gammaso
  ! isoiltype - soil type
  integer(kind=i4) :: isoiltype
  ! stheta - volumetric soil water content (m3/m3)
  real(kind=dp)                    :: stheta
  ! sattheta - saturation volumetric soil water content (m3/m3)
  real(kind=dp)                    :: sattheta
  ! rtheta - residual volumetric soil water content (m3/m3)
  real(kind=dp)                    :: rtheta
  ! sbcoef - Clapp and Hornberger exponent
  real(kind=dp)                    :: sbcoef
  ! satphi - saturation matric potential (suction) (m)
  real(kind=dp)                    :: satphi
  ! dsoil - depth of topsoil (cm)
  real(kind=dp)                    :: dsoil
  ! csoil - soil compensation points (molecules/cm3)
  real(kind=dp), dimension(ninteg) :: csoil

  ! saved soil physics data
  ! vsh2oout - soil exchange coefficient for water vapor (cm/s)
  real(kind=dp), allocatable     :: vsh2oout(:)
  ! qsoilout - effective soil molar humidity (mol/cm3)
  real(kind=dp), allocatable     :: qsoilout(:)
  ! effrhsoilout - effective soil fractional relative humidity ()
  real(kind=dp), allocatable     :: effrhsoilout(:)
  ! rbgout - ground boundary layer resistance (s/cm)
  real(kind=dp), allocatable     :: rbgout(:)
  ! gbgout - ground aerodynamic conductance (mol/m2-s)
  real(kind=dp), allocatable     :: gbgout(:)
  ! rsoilout - resistance to diffusion through soil pore space (s/cm)
  real(kind=dp), allocatable     :: rsoilout(:) 
  ! tsoilkout - soil temperature over simulation (K)
  real(kind=dp), allocatable     :: tsoilkout(:)
  ! tk0out - lower boundary condition for air temperature over simulation (K)
  real(kind=dp), allocatable     :: tk0out(:)

  ! Soil type data
  ! Soil type strings
   character(len=16), parameter, dimension(11) ::  &
                     ssoiltype     = (/         &
                                   'sand            ', &
                                   'loamy sand      ', &
                                   'sandy loam      ', &
                                   'loam            ', &
                                   'silt loam       ', &
                                   'sandy clay loam ', &
                                   'clay loam       ', &
                                   'silty clay loam ', &
                                   'sandy clay      ', &
                                   'silty clay      ', &
                                   'clay            ' &
                                             /)
  ! xsattheta - volumetric soil water content at saturation (m3/m3)
  real(kind=dp), dimension(11)     :: xsattheta
  data xsattheta /0.437,0.437,0.453,0.463,0.501,0.398,0.464,0.471,0.430,0.479,0.475/
  ! xrtheta - residual volumetric soil water content (m3/m3)
  real(kind=dp), dimension(11)     :: xrtheta
  data xrtheta /0.020,0.035,0.041,0.027,0.015,0.068,0.075,0.040,0.109,0.056,0.090/
  ! xsbcoef - Clapp-Hornberger exponents
  real(kind=dp), dimension(11)     :: xsbcoef
  data xsbcoef /4.05,4.38,4.90,5.39,5.30,7.12,8.52,7.75,10.4,10.4,11.4/
  ! xsatphi - soil matric potential (suction) at saturation (cm)
  real(kind=dp), dimension(11)     :: xsatphi
  data xsatphi /12.1,9.0,21.8,78.6,47.8,29.9,35.6,63.0,15.3,49.0,40.5/

  ! Soil thermal conductivity data
  ! Estimated from:
  ! Kersten, M. S. (1949) Thermal Properties of Soils. University of Minnesota, Engineering
  !  Experiment Station Bulletin No. 28 (http://hdl.handle.net/11299/124271).
  !
  ! ksoildry - thermal conductivity for dry soil
  real(kind=dp), dimension(11)     :: ksoildry
  data ksoildry /0.30,0.28,0.26,0.20,0.16,0.20,0.16,0.16,0.20,0.16,0.16/
  ! ksoilwet - thermal conductivity for wet soil
  real(kind=dp), dimension(11)     :: ksoilwet
  data ksoilwet /2.20,2.10,2.04,1.66,1.28,1.70,1.28,1.28,1.74,1.28,1.28/

  ! ncnpep - number of species with non-zero canopy emission potentials
  integer(kind=i4) :: ncnpep
  ! gammast - stomatal emission potential = [NH4+]/[H+]stomata for ammonia
  real(kind=dp), dimension(ninteg)       :: gammast

  ! tkout - saved temperature profiles (K)
  real(kind=dp), allocatable :: tkout(:,:)
  ! pmbout - saved pressure profiles (mb)
  real(kind=dp), allocatable :: pmbout(:,:)
  ! qhout - saved specific humidity profiles (g/kg)
  real(kind=dp), allocatable :: qhout(:,:)
  ! rhout - saved relative humidity profiles (%)
  real(kind=dp), allocatable :: rhout(:,:)
  ! ubarout - saved mean wind speed profiles (cm/s)
  real(kind=dp), allocatable :: ubarout(:,:)
  ! kvout - saved eddy diffusivity profiles (cm2/s)
  real(kind=dp), allocatable :: kvout(:,:)
  ! ppfdout - saved PPFD profiles (umol/m2-s)
  real(kind=dp), allocatable :: ppfdout(:,:)
  ! cairout - saved cair profiles (molecules/cm3)
  real(kind=dp), allocatable :: cairout(:,:)
  ! h2oout - saved h2o profiles (molecules/cm3)
  real(kind=dp), allocatable :: h2oout(:,:) 

  ! timeout - saved t values corresponding to output saved in cout
  real(kind=dp), allocatable :: timeout(:)
  ! sdtout - saved sdatetime values corresponding to output saved in cout
  character(len=19), allocatable :: sdtout(:)

  ! Canopy physics data
  ! ppfd_direct - direct beam component of incoming PPFD (umol/m2-s)
  real(kind=dp)                  :: ppfd_direct
  ! ppfd_diffus - diffuse component of incoming PPFD (umol/m2-s)
  real(kind=dp)                  :: ppfd_diffus
  ! nir_direct - direct beam component of incoming NIR (W/m2)
  real(kind=dp)                  :: nir_direct
  ! nir_diffus - diffuse component of incoming NIR (W/m2)
  real(kind=dp)                  :: nir_diffus
  ! ppfd_sun - sunlit fraction PPFD (umol/m2-s)
  real(kind=dp), dimension(npts) :: ppfd_sun
  ! ppfd_shd - shaded fraction PPFD (umol/m2-s)
  real(kind=dp), dimension(npts) :: ppfd_shd
  ! ppfd_wgt - sun/shade weighted PPFD (umol/m2-s)
  real(kind=dp), dimension(npts) :: ppfd_wgt
  ! nir_sun - sunlit fraction NIR (W/m2)
  real(kind=dp), dimension(npts) :: nir_sun
  ! nir_shd - shaded fraction NIR (W/m2)
  real(kind=dp), dimension(npts) :: nir_shd
  ! nir_wgt - sun/shade weighted NIR (W/m2)
  real(kind=dp), dimension(npts) :: nir_wgt
  ! lw_up - upwelling long-wave radiation (W/m2)
  real(kind=dp), dimension(npts) :: lw_up
  ! lw_dn - downwelling long-wave radiation (W/m2)
  real(kind=dp), dimension(npts) :: lw_dn
  ! rt_sun - total absorbed PPFD and NIR radiation in sunlit fraction (W/m2)
  real(kind=dp), dimension(npts) :: rt_sun
  ! rt_shd - total absorbed PPFD and NIR radiation in shaded fraction (W/m2)
  real(kind=dp), dimension(npts) :: rt_shd
  ! rt_wgt - sun/shade weighted total absorbed PPFD and NIR radiation (W/m2)
  real(kind=dp), dimension(npts) :: rt_wgt
  ! rabs_sun - total absorbed radiation (PPFD+NIR+IR) in sunlit fraction (W/m2)
  real(kind=dp), dimension(npts) :: rabs_sun
  ! rabs_shd - total absorbed radiation (PPFD+NIR+IR) in shaded fraction (W/m2)
  real(kind=dp), dimension(npts) :: rabs_shd
  ! rabs_wgt - sun/shade weighted total absorbed radiation (PPFD+NIR+IR) (W/m2)
  real(kind=dp), dimension(npts) :: rabs_wgt
  ! fsun - sunlit fraction of canopy ()
  real(kind=dp), dimension(npts) :: fsun
  ! fshd - shaded fraction of canopy ()
  real(kind=dp), dimension(npts) :: fshd
  ! tl_sun - leaf temperature in sunlit fraction (K)
  real(kind=dp), dimension(npts) :: tl_sun
  ! tl_shd - leaf temperature in shaded fraction (K)
  real(kind=dp), dimension(npts) :: tl_shd
  ! tlp_sun - previous iteration leaf temperature in sunlit fraction (K)
  real(kind=dp), dimension(npts) :: tlp_sun
  ! tlp_shd - previous iteration leaf temperature in shaded fraction (K)
  real(kind=dp), dimension(npts) :: tlp_shd
  ! tl_wgt - sun/shade weighted leaf temperature (K)
  real(kind=dp), dimension(npts) :: tl_wgt
  ! gs_sun - leaf stomatal conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gs_sun
  ! gs_shd - leaf stomatal conductance in shaded fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gs_shd
  ! gs_wgt - sun/shade weighted leaf stomatal conductance (mol/m2-s)
  real(kind=dp), dimension(npts) :: gs_wgt

  ! gl_sun - leaf water vapor conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gl_sun  
  ! gl_shd - leaf water vapor conductance in shaded fraction (mol/m2-s)
  real(kind=dp), dimension(npts) :: gl_shd
  ! vl_sun - leaf water vapor exchange velocity in sunlit fraction (cm/s)
  real(kind=dp), dimension(npts) :: vl_sun  
  ! vl_shd - leaf water vapor exchange velocity in shaded fraction (cm/s)
  real(kind=dp), dimension(npts) :: vl_shd
  ! gb - leaf boundary layer conductance (mol/m2-s)
  real(kind=dp), dimension(npts) :: gb
  ! esun - evapotranspiration flux for the sunlit canopy fraction (mol/cm2-s)
  real(kind=dp), dimension(npts) :: esun
  ! eshd - evapotranspiration flux for the shaded canopy fraction (mol/cm2-s)
  real(kind=dp), dimension(npts) :: eshd
  ! etot - total evapotranspiration flux for the canopy (mol/cm2-s)
  real(kind=dp), dimension(npts) :: etot
  ! hsun - sensible heat flux for the sunlit canopy fraction (W/cm2)
  real(kind=dp), dimension(npts) :: hsun
  ! hshd - sensible heat flux for the shaded canopy fraction (W/cm2)
  real(kind=dp), dimension(npts) :: hshd
  ! htot - total sensible heat flux for the canopy (W/cm2)
  real(kind=dp), dimension(npts) :: htot

  ! rs_sun - leaf stomatal resistance in sunlit fraction (s/m)
  real(kind=dp), dimension(npts) :: rs_sun
  ! rs_shd - leaf stomatal resistance in shaded fraction (s/m)
  real(kind=dp), dimension(npts) :: rs_shd
  ! rs_wgt - sun/shade weighted leaf stomatal resistance (s/m)
  real(kind=dp), dimension(npts) :: rs_wgt
  ! anet_sun - net photosynthetic assimilation in sunlit fraction (umol/m2-s)
  real(kind=dp), dimension(npts) :: anet_sun
  ! anet_shd - net photosynthetic assimilation in shaded fraction (umol/m2-s)
  real(kind=dp), dimension(npts) :: anet_shd
  ! anet_wgt - sun/shade weighted net photosynthetic assimilation (umol/m2-s)
  real(kind=dp), dimension(npts) :: anet_wgt
  ! epsleaf - leaf emissivity ()
  real(kind=dp), parameter       :: epsleaf=0.97
  ! epsgrnd - soil/litter surface emissivity ()
  real(kind=dp), parameter       :: epsgrnd=0.95
  ! sbsig - Stefan-Boltzmann constant (W/m2-K4)
  real(kind=dp), parameter       :: sbsig=5.67D-08
  ! cpair - heat capacity of air, J/mol-K
  real(kind=dp), parameter       :: cpair=29.3
  ! cca - atmospheric CO2 mixing ratio (umol/mol)
  real(kind=dp), parameter       :: cca=405.0

  ! x - leaf angle distribution parameter in ellipsoidal kb function
  real(kind=dp)                  :: x
  ! kd - diffuse radiation extinction coefficient
  real(kind=dp)                  :: kd
  ! g0 - Medlyn et al. stomatal conductance parameter
  real(kind=dp)                  :: g0
  ! g1 - Medlyn et al. stomatal conductance parameter
  real(kind=dp)                  :: g1
  ! gsmax - maximum stomatal conductance (mol/m2-s)
  real(kind=dp), parameter       :: gsmax=10.0
  ! dleaf - characteristic leaf dimension = 0.75*(leaf width) (m)
  real(kind=dp)                  :: dleaf
  ! ltol - iteration tolerance used in calculating leaf temperature
  real(kind=dp), parameter       :: ltol=0.001
  ! rtol - iteration tolerance used in calculating overall canopy physics
  real(kind=dp), parameter       :: rtol=0.01
 
  ! saved canopy physics data
  ! ppfddirout - direct beam component of incoming PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfddirout(:)
  ! ppfddifout - diffuse component of incoming PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfddifout(:)
  ! nirdirout - direct beam component of incoming NIR (W/m2)
  real(kind=dp), allocatable     :: nirdirout(:)
  ! nirdifout - diffuse component of incoming NIR (W/m2)
  real(kind=dp), allocatable     :: nirdifout(:)
  ! ppfdsunout - sunlit fraction PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfdsunout(:,:)
  ! ppfdshdout - shaded fraction PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfdshdout(:,:)
  ! ppfdwgtout - sun/shade weighted PPFD (umol/m2-s)
  real(kind=dp), allocatable     :: ppfdwgtout(:,:)
  ! nirsunout - sunlit fraction NIR (W/m2)
  real(kind=dp), allocatable     :: nirsunout(:,:)
  ! nirshdout - shaded fraction NIR (W/m2)
  real(kind=dp), allocatable     :: nirshdout(:,:)
  ! nirwgtout - sun/shade weighted NIR (W/m2)
  real(kind=dp), allocatable     :: nirwgtout(:,:)
  ! lwupout - upwelling long-wave radiation (W/m2)
  real(kind=dp), allocatable     :: lwupout(:,:)
  ! lwdnout - downwelling long-wave radiation (W/m2)
  real(kind=dp), allocatable     :: lwdnout(:,:)
  ! rtsunout - total absorbed PPFD and NIR radiation in sunlit fraction (W/m2)
  real(kind=dp), allocatable     :: rtsunout(:,:)
  ! rtshdout - total absorbed PPFD and NIR radiation in shaded fraction (W/m2)
  real(kind=dp), allocatable     :: rtshdout(:,:)
  ! rtwgtout - sun/shade weighted total absorbed PPFD and NIR radiation (W/m2)
  real(kind=dp), allocatable     :: rtwgtout(:,:)
  ! rabssunout - total absorbed radiation (PPFD+NIR+IR) in sunlit fraction (W/m2)
  real(kind=dp), allocatable     :: rabssunout(:,:)
  ! rabsshdout - total absorbed radiation (PPFD+NIR+IR) in shaded fraction (W/m2)
  real(kind=dp), allocatable     :: rabsshdout(:,:)
  ! rabswgtout - sun/shade weighted total absorbed radiation (PPFD+NIR+IR) (W/m2)
  real(kind=dp), allocatable     :: rabswgtout(:,:)
  ! fsunout - sunlit fraction of canopy ()
  real(kind=dp), allocatable     :: fsunout(:,:)
  ! fshdout - shaded fraction of canopy ()
  real(kind=dp), allocatable     :: fshdout(:,:)
  ! tlsunout - leaf temperature in sunlit fraction (K)
  real(kind=dp), allocatable     :: tlsunout(:,:)
  ! tlshdout - leaf temperature in shaded fraction (K)
  real(kind=dp), allocatable     :: tlshdout(:,:)
  ! tlwgtout - sun/shade weighted leaf temperature (K)
  real(kind=dp), allocatable     :: tlwgtout(:,:)
  ! gssunout - leaf stomatal conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), allocatable     :: gssunout(:,:)
  ! gsshdout - leaf stomatal conductance in shaded fraction (mol/m2-s)
  real(kind=dp), allocatable     :: gsshdout(:,:)
  ! gswgtout - sun/shade weighted leaf stomatal conductance (mol/m2-s)
  real(kind=dp), allocatable     :: gswgtout(:,:)
  ! rssunout - leaf stomatal resistance in sunlit fraction (s/m)
  real(kind=dp), allocatable     :: rssunout(:,:)
  ! rsshdout - leaf stomatal resistance in shaded fraction (s/m)
  real(kind=dp), allocatable     :: rsshdout(:,:)
  ! rswgtout - sun/shade weighted leaf stomatal resistance (s/m)
  real(kind=dp), allocatable     :: rswgtout(:,:)

  ! hsunout - sensible heat flux for the sunlit fraction of the canopy (W/cm2)
  real(kind=dp), allocatable     :: hsunout(:,:)
  ! hshdout - sensible heat flux for the shaded fraction of the canopy (W/cm2)
  real(kind=dp), allocatable     :: hshdout(:,:)
  ! htotout - total sensible heat flux for the canopy (W/cm2)
  real(kind=dp), allocatable     :: htotout(:,:)

  ! esunout - evapotranspiration flux for the sunlit fraction of the canopy (W/cm2)
  real(kind=dp), allocatable     :: esunout(:,:)
  ! eshdout - evapotranspiration flux for the shaded fraction of the canopy (W/cm2)
  real(kind=dp), allocatable     :: eshdout(:,:)
  ! etotout - total evapotranspiration flux for the canopy (W/cm2)
  real(kind=dp), allocatable     :: etotout(:,:)

  ! gbout - leaf boundary layer conductance (mol/m2-s)
  real(kind=dp), allocatable     :: gbout(:,:)

  ! gl_sunout - leaf water vapor conductance in sunlit fraction (mol/m2-s)
  real(kind=dp), allocatable     :: gl_sunout(:,:)
  ! gl_shdout - leaf water vapor conductance in shaded fraction (mol/m2-s)
  real(kind=dp), allocatable     :: gl_shdout(:,:)

  ! vl_sunout - leaf water vapor exchange velocity in sunlit fraction (cm/s)
  real(kind=dp), allocatable     :: vl_sunout(:,:)
  ! vl_shdout - leaf water vapor exchange velocity in shaded fraction (cm/s)
  real(kind=dp), allocatable     :: vl_shdout(:,:)

  ! anetsunout - net photosynthetic assimilation in sunlit fraction (umol/m2-s)
  real(kind=dp), allocatable     :: anetsunout(:,:)
  ! anetshdout - net photosynthetic assimilation in shaded fraction (umol/m2-s)
  real(kind=dp), allocatable     :: anetshdout(:,:)
  ! anetwgtout - sun/shade weighted net photosynthetic assimilation (umol/m2-s)
  real(kind=dp), allocatable     :: anetwgtout(:,:)

  ! Surface Energy Balance components
  ! rnout - net radiation flux vertical profiles (W/m2)
  real(kind=dp), allocatable     :: rnout(:,:)
  ! hout - sensible heat flux vertical profiles (W/m2)
  real(kind=dp), allocatable     :: hout(:,:)
  ! leout - latent heat flux vertical profiles (W/m2)
  real(kind=dp), allocatable     :: leout(:,:)
  ! gout - ground heat flux (W/m2)
  real(kind=dp), allocatable     :: gout(:)
  ! hflux - measured sensible heat flux at current time (W/m2)
  real(kind=dp)                  :: hflux
  ! leflux - measured latent heat flux at current time (W/m2)
  real(kind=dp)                  :: leflux
  ! gflux - measured ground heat flux (W/m2)
  real(kind=dp)                  :: gflux

  ! Meteorological data
  ! Met data at reference height
  ! pmbzref - air pressure (mb)
  real(kind=dp)                  :: pmbzref
  ! umszref - wind speed (m/s)
  real(kind=dp)                  :: umszref
  ! ubzref - wind speed (cm/s)
  real(kind=dp)                  :: ubzref
  ! rhzref - relative humidity (%)
  real(kind=dp)                  :: rhzref
  ! qzref - molar water vapor concentration (moles/m3)
  real(kind=dp)                  :: qzref
  ! tkzref - air temperature (K)
  real(kind=dp)                  :: tkzref
  ! ppfdzref - PPFD (umol/m2-s) 
  real(kind=dp)                  :: ppfdzref
  ! sradzref - solar radiation (W/m2)
  real(kind=dp)                  :: sradzref
  ! lwdnzref - downwelling longwave radiation (W/m2)
  real(kind=dp)                  :: lwdnzref
  ! razref - aerodynamic resistance (s/cm)
  real(kind=dp)                  :: razref
  ! kvzref - measured Kv at z=H (cm2/s)
  real(kind=dp)                  :: kvzref

  ! eatm - water vapor pressure (kPa)
  real(kind=dp)                  :: eatm
  ! tk - vertical profile of temperature at current simulation time (K)
  real(kind=dp), dimension(npts) :: tk
  ! qh - vertical profile of specific humidity at current simulation time (g/kg)
  real(kind=dp), dimension(npts) :: qh 
  ! tkp - previous iteration vertical profile of temperature at current simulation time (K)
  real(kind=dp), dimension(npts) :: tkp
  ! qhp - previous iteration vertical profile of specific humidity at current simulation time (g/kg)
  real(kind=dp), dimension(npts) :: qhp 
  ! pmb - vertical profile of pressure at current simulation time  (mb)
  real(kind=dp), dimension(npts) :: pmb
  ! kv - eddy diffusivity at current simulation time (cm2/s)
  real(kind=dp), dimension(npts) :: kv
  ! ppfd - vertical profile of photosynthetic photon flux (umol/m2-s) at current simulation time
  real(kind=dp), dimension(npts) :: ppfd
  ! ubar - vertical profile of mean wind speed at current simulation time (cm/s)
  real(kind=dp), dimension(npts) :: ubar

  ! Calculated meteorological data
  ! cair - vertical air density at current simulation time (molec/cm3)
  real(kind=dp), dimension(npts) :: cair
  ! h2o - vertical water vapor concentration at current simulation time (molec/cm3)
  real(kind=dp), dimension(npts) :: h2o

  ! Physical-Chemical data
  ! mdiffstp - molecular diffusivities of species in air at 0degC and 1 atm (cm2/s)
  real(kind=dp), dimension(ntotal) :: mdiffstp
  ! mdiffstp_default - default value of mdiffstp (cm2/s) for species with no reliable data
  real(kind=dp), parameter         :: mdiffstp_default=0.100
  ! hstar - effective Henry's Law coefficients (M/atm)
  real(kind=dp), dimension(ntotal) :: hstar
  ! hstar_default - default value of hstar (M/atm)
  real(kind=dp), parameter         :: hstar_default=1.000
  ! f0 - reactivity parameter used in resistance calculation (dimensionless)
  real(kind=dp), dimension(ntotal) :: f0
  ! f0_default - default value of f0
  real(kind=dp), parameter         :: f0_default=0.0
  ! molecmass - molecular mass of simulated species
  real(kind=dp), dimension(ntotal) :: molecmass

  ! simulation control file unit number
  integer(kind=i4), parameter :: UCTRL=11
  ! simulation control file name
  character(len=14), parameter :: filectrl='canaccCTRL.dat'
  ! met data file unit number
  integer(kind=i4), parameter :: UENV=13
  ! met data file name
  character(len=16) :: envfile
  ! canopy morphology file unit number
  integer(kind=i4), parameter :: UCNPY=14
  ! canopy morphology data file name
  character(len=16) :: cnpyfile
  ! soil data file unit number
  integer(kind=i4), parameter :: USOIL=19
  ! soil data file name
  character(len=16) :: slfile

  ! simdescr - descriptive simulation identification
  character(len=100) :: simdescr
  ! simname - unique simulation name
  character(len=16) :: simname

  ! slat - latitude (deg)
  real(kind=dp) :: slat
  ! slon - longitude (deg)
  real(kind=dp) :: slon

  ! year - year of simulation
  integer(kind=i4) :: year
  ! month - month of simulation (at start)
  integer(kind=i4) :: month
  ! daymonth - day of the month (at start)
  integer(kind=i4) :: daymonth
  ! doy - day of year at beginning of simulation
  integer(kind=i4) :: doy
  ! hz, mz, sz - local time at beginning of simulation
  integer(kind=i4) :: hz, mz, sz
  ! tzone - time zone differential
  integer(kind=i4) :: tzone

  ! unit number for output CTRL simulation summary
  integer(kind=i4), parameter :: USUMM=20
  ! name of CTRL simulation summary file
  character(len=32) :: simsummary
  ! unit number for runtime output file
  integer(kind=i4), parameter :: URUN=21
  ! name of runtime output file
  character(len=32) :: simrunfile
  ! grid definition file unit number
  integer(kind=i4), parameter :: UGRID=18
  ! grid definition file name
  character(len=16) :: grdfile
  ! begining output file unit number for species files
  integer(kind=i4), parameter :: UOUT=30

  ! stdsp - indices of species printed to STDOUT
  integer(kind=i4) :: nstdsp
  integer(kind=i4), allocatable :: stdsp(:)

  ! options
  ! OPT_SIMTYPE - dynamic or steady-state simulation
  integer(kind=i4)            :: OPT_SIMTYPE
  integer(kind=i4), parameter :: DYNAMIC=1      ! solar time changes with t
  integer(kind=i4), parameter :: SSTATE=2       ! solar time fixed at tstart
  ! INTGWVAP - option for water vapor integration
  logical :: INTGWVAP
  ! INTGTAIR - option for air temperature integration
  logical :: INTGTAIR
  !
  ! KVSELECT - Kv parameterization selection
  integer(kind=i4)            :: KVSELECT
  integer(kind=i4), parameter :: KVRA=1         ! Aerodynamic resistance parameterization
  integer(kind=i4), parameter :: KVSTULL=2      ! Stull-based
  integer(kind=i4), parameter :: KVULKE=3       ! Ulke-based
  integer(kind=i4), parameter :: KVCOWEETA=4    ! Based on measured EC data at Coweeta

  ! BCTYPE - Flux or Dirichlet BC at domain top
  integer(kind=i4)            :: BCTYPE
  integer(kind=i4), parameter :: BCFLUX=1
  integer(kind=i4), parameter :: BCCONST=2

  ! SKYLW - how to determine downwelling longwave radiation
  integer(kind=i4)            :: SKYLW
  integer(kind=i4), parameter :: MEASURED=1
  integer(kind=i4), parameter :: ESTIMATED=2

  ! Sensitivity simulation factors
  ! Kv profile sensitivity
  real(kind=dp)            :: senskv

  ! Avogadro's number
  real(kind=dp), parameter :: navo=6.02D+023

  ! pi
  real(kind=dp), parameter :: pi=4.0*atan(1.0)

end module GlobalData
