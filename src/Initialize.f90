!======================================================================================================================!
!                                                                                                                      !
!     Program:      CANACC                                                                                             !
!                   Canopy physics model for ACCESS                                                                    !
!                                                                                                                      !
!     Version:      2.0.0                                                                                              !
!                                                                                                                      !
!======================================================================================================================!
!                                                                                                                      !
!     Initiated:    Mar 2017                                                                                           !
!     Last Update:  July 2018                                                                                          !
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
!     Module:       Initialize                                                                                         !
!                                                                                                                      ! 
!     Description:  Contains routines related to model initialization                                                  !
!                                                                                                                      !
!======================================================================================================================!
module Initialize
  use GlobalData
  use EnvironData
  use PhysChemData
  use Utils
  use Output
  implicit none

  private SetSimulationData, ReadVerticalGridData, GetCanopyData, GetSoilData
  public InitializeModel

contains

!**********************************************************************************************************************!
! InitializeModel - performs all the model initialization steps 
!                   Called from Main.f90
!**********************************************************************************************************************!
subroutine InitializeModel()
  real(kind=dp)    :: delta, leap, mjd, hrlocal, hrutc, time21
  integer(kind=i4) :: m
  integer(kind=i4), dimension(12) :: months

  data months /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

  ! read CTRL file
  call SetSimulationData()

  ! calculate doy (day-of-year)
  doy = 0.0
  do m=1,month-1
    doy = doy + dble(months(m))
  end do
  doy = doy + dble(daymonth)

  ! calculate decimal UTC hour
  hrlocal = dble(hz) + (dble(mz)/60.) + (dble(sz)/3600.)
  hrutc = hrlocal - tzone

  ! years from 1949
  delta = dble(year)-1949.

  ! leap days over this period (note correction below)
  leap = aint(delta/4.)

  ! modified Julian Day
  mjd=32916.5+delta*365.+leap+dble(doy)+hrutc/24.

  ! the last year of century is not a leap year unless divisible by 400
  if (dmod(dble(year),100.0D+00) .eq. 0.0 .and. dmod(dble(year),400.0D+00) .ne. 0.0) then
    mjd=mjd-1.
  end if

  ! time in days since Jan 1, 2000 Greenwich Noon
  time21 = mjd-51545.0 
 
  ! time in seconds since Jan 1, 2000 Greenwich Noon
  ! simulation time is kept in seconds
  tstart=time21*24.*3600.

  zadeg = SolarZenithAngle(time21)
  zarad = zadeg*pi/180.

  ! ending time
  tend=tstart+dtout*dble(ntout)

  dthalf = 0.5_dp*dt                   ! half time step
  t=tstart                             ! simulation time start
  nt=0                                 ! output number start

  ! allocate storage for final output arrays
  allocate(ra(0:ntout-1))
  allocate(rakv(0:ntout-1))
  allocate(zolout(0:ntout-1))
  allocate(gaero(0:ntout-1))

  allocate(tkout(npts,0:ntout-1))
  allocate(pmbout(npts,0:ntout-1))
  allocate(qhout(npts,0:ntout-1)) 
  allocate(ubarout(npts,0:ntout-1))
  allocate(kvout(npts,0:ntout-1))
  allocate(ppfdout(npts,0:ntout-1))
  allocate(cairout(npts,0:ntout-1))
  allocate(h2oout(npts,0:ntout-1))
  allocate(rhout(npts,0:ntout-1))

  allocate(ppfddirout(0:ntout-1))
  allocate(ppfddifout(0:ntout-1))
  allocate(nirdirout(0:ntout-1))
  allocate(nirdifout(0:ntout-1))
  allocate(ppfdsunout(npts,0:ntout-1))
  allocate(ppfdshdout(npts,0:ntout-1))
  allocate(ppfdwgtout(npts,0:ntout-1))
  allocate(nirsunout(npts,0:ntout-1))
  allocate(nirshdout(npts,0:ntout-1))
  allocate(nirwgtout(npts,0:ntout-1))
  allocate(lwupout(npts,0:ntout-1))
  allocate(lwdnout(npts,0:ntout-1))
  allocate(rtsunout(npts,0:ntout-1))
  allocate(rtshdout(npts,0:ntout-1))
  allocate(rtwgtout(npts,0:ntout-1))
  allocate(rabssunout(npts,0:ntout-1))
  allocate(rabsshdout(npts,0:ntout-1))
  allocate(rabswgtout(npts,0:ntout-1))
  allocate(fsunout(npts,0:ntout-1))
  allocate(fshdout(npts,0:ntout-1))
  allocate(tlsunout(npts,0:ntout-1))
  allocate(tlshdout(npts,0:ntout-1))
  allocate(tlwgtout(npts,0:ntout-1))
  allocate(gssunout(npts,0:ntout-1))
  allocate(gsshdout(npts,0:ntout-1))
  allocate(gswgtout(npts,0:ntout-1))
  allocate(rssunout(npts,0:ntout-1))
  allocate(rsshdout(npts,0:ntout-1))
  allocate(rswgtout(npts,0:ntout-1))
  allocate(anetsunout(npts,0:ntout-1))
  allocate(anetshdout(npts,0:ntout-1))
  allocate(anetwgtout(npts,0:ntout-1))

  allocate(timeout(0:ntout-1))
  allocate(sdtout(0:ntout-1))

  ! make sure arrays are zeroed initially
  ra=0.0_dp
  rakv=0.0_dp
  zolout=0.0_dp
  gaero=0.0_dp
  tkout=0.0_dp
  pmbout=0.0_dp
  qhout=0.0_dp
  ubarout=0.0_dp
  kvout=0.0_dp
  ppfdout=0.0_dp
  cairout=0.0_dp
  h2oout=0.0_dp
  rhout=0.0_dp
  ppfdsunout=0.0_dp
  ppfdshdout=0.0_dp
  nirsunout=0.0_dp
  nirshdout=0.0_dp
  lwupout=0.0_dp
  lwdnout=0.0_dp
  rtsunout=0.0_dp
  rtshdout=0.0_dp
  rabssunout=0.0_dp
  rabsshdout=0.0_dp
  fsunout=0.0_dp
  fshdout=0.0_dp
  tlsunout=0.0_dp
  tlshdout=0.0_dp
  gssunout=0.0_dp
  gsshdout=0.0_dp
  rssunout=0.0_dp
  rsshdout=0.0_dp
  anetsunout=0.0_dp
  anetshdout=0.0_dp

  ! read z values of vertical grid
  call ReadVerticalGridData()

  ! read canopy morphology data
  call GetCanopyData()

  ! read soil data file
  call GetSoilData()

  return
end subroutine InitializeModel

!**********************************************************************************************************************!
! ReadVerticalGridData - read vertical grid definition from an input file
!
!**********************************************************************************************************************!
subroutine ReadVerticalGridData()
  integer(kind=i4) :: i, k, nptsdef
  real(kind=dp)    :: zm, dzm, dzm2

  ! read GridDef file
  open(unit=UGRID,file=('./data/' // grdfile))
  read(UGRID,*) nptsdef, hc, ncnpy, alfa

  ! check file consistency
  if (nptsdef .ne. npts) then
    write(*,101) grdfile
    write(*,102) nptsdef, npts
    close(UGRID)
    stop
  end if

  ! ok, so read data
  do i=1,npts
    read(UGRID,*) k, zm, dzm
    z(k) = zm*100.0
    if (k == 2) dzm2=dzm
  end do

  dzhc = dzm2*100.0

  ! domain boundaries
  z0 = z(1)
  zi = z(npts)

  ! canopy height in cm
  hccm = hc*100.0

  ! set aerodynamic parameters
  !  Based on the recommendations of ...
  !  Monteith and Unsworth (2013) Principles of Environmental Physics
  d   = 0.667*hccm
  z0m = 0.097*hccm
  z0h = z0m/7.3

  close(UGRID)

101 format('***GridDef file ', a, ' is inconsistent with current CANACC configuration!')
102 format('***GridDef npts = ', i4, /'***GlobalData npts = ', i4)
  return
end subroutine ReadVerticalGridData

!**********************************************************************************************************************!
! SetSimulationData - reads input data from CTRL file for simulation
!                     and writes to STDOUT and the simulation summary file
!**********************************************************************************************************************!
subroutine SetSimulationData()
  integer(kind=i4)  :: j, isp
  character(len=35) :: outdir
  character(len=35) :: strmkdir
  logical           :: eoutdir
  integer(kind=i4), parameter   :: nhlines=6

  ! read CTRL file
  open(unit=UCTRL,file=('./ctrl/' // filectrl))
  ! skip header lines
  do j=1,nhlines
    read(UCTRL,*)
  end do
  read(UCTRL,101) simdescr
  read(UCTRL,*) 
  read(UCTRL,102) simname 
  read(UCTRL,*) OPT_SIMTYPE
  read(UCTRL,*) slat
  read(UCTRL,*) slon
  read(UCTRL,*) year
  read(UCTRL,*) month
  read(UCTRL,*) daymonth
  read(UCTRL,*) hz, mz, sz
  read(UCTRL,*) tzone
  read(UCTRL,*) dtout
  read(UCTRL,*) ntout
  read(UCTRL,*) INTGWVAP
  read(UCTRL,*) INTGTAIR
  read(UCTRL,*)
  read(UCTRL,*) grdfile
  read(UCTRL,*)
  read(UCTRL,*) envfile
  read(UCTRL,*)
  read(UCTRL,*) cnpyfile
  read(UCTRL,*)
  read(UCTRL,*) slfile
  close(UCTRL)

101 format(a100)
102 format(a16)

  ! write CTRL data to STDOUT
  write(6,300)
  write(6,2000)
  write(6,*) simdescr
  write(6,2010) simname 
  if (OPT_SIMTYPE .eq. DYNAMIC) then
    write(6,2015)
  else
    write(6,2016)
  end if
  write(6,2020) slat
  write(6,2030) slon
  write(6,2031) year
  write(6,2040) month
  write(6,2041) daymonth
  write(6,2050) hz, mz, sz
  write(6,2051) tzone
  write(6,2070) dtout
  write(6,2071) ntout
  write(6,2200) 
  write(6,2172) INTGWVAP
  write(6,2173) INTGTAIR
  write(6,2210) 
  write(6,2179) grdfile
  write(6,2190) envfile
  write(6,2300) cnpyfile
  write(6,2304) slfile
  write(6,300)

  ! if output directory in 'out' does not exist, create it
  ! and all necessary subdirectories
  outdir='./out/' // trim(simname) // '/.'
  inquire(file=outdir, exist=eoutdir)
  if(eoutdir .eqv. .FALSE.) then
    strmkdir = 'mkdir ./out/' // trim(simname)
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/grid'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/met'
    call system(strmkdir)
    strmkdir = 'mkdir ./out/' // trim(simname) // '/canopy'
    call system(strmkdir)
  end if

  ! open simulation runtime file to capture STDOUT
  ! runtime file is normally closed by Utils:CleanUp
  ! at the end of the simulation
  simrunfile='./out/' // trim(simname) // '/CANACC_STD.out'
  simrunfile=trim(simrunfile)
  open(unit=URUN,file=simrunfile)

  ! write header for runtime file
  write(URUN,300)
  write(URUN,*) simdescr
  write(URUN,2010) simname 
  write(URUN,300)

  ! write CTRL data to USUMM 
  simsummary='./out/' // trim(simname) // '/CANACC_SUMM.out'
  simsummary=trim(simsummary)
  open(unit=USUMM,file=simsummary)
  write(USUMM,300)
  write(USUMM,2000)
  write(USUMM,*) simdescr
  write(USUMM,2010) simname 
  if (OPT_SIMTYPE .eq. DYNAMIC) then
    write(USUMM,2015)
  else
    write(USUMM,2016)
  end if
  write(USUMM,2020) slat
  write(USUMM,2030) slon
  write(USUMM,2031) year
  write(USUMM,2040) month
  write(USUMM,2041) daymonth
  write(USUMM,2050) hz, mz, sz
  write(USUMM,2051) tzone
  write(USUMM,2070) dtout
  write(USUMM,2071) ntout
  write(USUMM,2200) 
  write(USUMM,2172) INTGWVAP
  write(USUMM,2173) INTGTAIR
  write(USUMM,2210) 
  write(USUMM,2179) grdfile
  write(USUMM,2190) envfile
  write(USUMM,2300) cnpyfile
  write(USUMM,2304) slfile
  write(USUMM,300)
  close(USUMM)

2000 format(/'Starting CANACC simulation...'/)
2010 format(' Short sim name = ', a)
2015 format(' OPT_SIMTYPE    = DYNAMIC')
2016 format(' OPT_SIMTYPE    = SSTATE')
2020 format(' Latitude       = ', f7.2, ' deg')
2030 format(' Longitude      = ', f7.2, ' deg')
2031 format(' Year           = ', i4)
2040 format(' Month          = ', i4)
2041 format(' Day            = ', i4)
2050 format(' Start time     = ', i2.2, ':', i2.2, ':', i2.2, ' LT')
2051 format(' Time zone diff = ', i3)
2070 format(' Input dt       = ', f7.1, ' s')
2071 format(' # of input dt  = ', i4)
2090 format(100(1x,a))
2172 format(' INTGWVAP       = ', l2)
2173 format(' INTGTAIR       = ', l2)
2179 format(' GRD file name  = ', a)
2190 format(' ENV file name  = ', a)
2200 format(/' Model Options:')
2210 format(/' Input Files:')
2300 format(' CNPY file name  = ', a)
2304 format(' SOIL  file name = ', a)
300 format(80('='))

  return
end subroutine SetSimulationData

!**********************************************************************************************************************!
! GetSoilData - reads soil data file
!**********************************************************************************************************************!
subroutine GetSoilData()
  integer(kind=i4)  :: j
  integer(kind=i4), parameter  :: nhlines=21

  ! read soil file
  open(unit=USOIL,file=('./data/' // slfile))
  ! skip header lines
  do j=1,nhlines
    read(USOIL,*)
  end do

  ! read soil type
  read(USOIL,*) isoiltype

  ! read volumetric soil water content
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*) stheta

  ! read depth of topsoil
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*)
  read(USOIL,*) dsoil

  close(USOIL)

  ! Based on soil type, set soil data
  sattheta = xsattheta(isoiltype)
  rtheta   = xrtheta(isoiltype)
  sbcoef   = xsbcoef(isoiltype)

  return
end subroutine GetSoilData

!**********************************************************************************************************************!
! GetCanopyData - reads canopy morphology data from input canopy file
!**********************************************************************************************************************!
subroutine GetCanopyData()
  integer(kind=i4) :: i, j
  integer(kind=i4), parameter :: nchlines=15
  integer(kind=i4) :: cm_npts
  real(kind=dp)    :: cm_z0, cm_zi, cm_hc, cm_alfa

  open(unit=UCNPY,file=('./data/' // cnpyfile))
  do j=1,nchlines
    read(UCNPY,*)
  end do
  ! read grid consistency data and check
  read(UCNPY,*)   cm_npts
  read(UCNPY,*)   cm_z0
  read(UCNPY,*)   cm_zi 
  read(UCNPY,*)   cm_hc
  read(UCNPY,*)   cm_alfa
  if (cm_npts .ne. npts) then
    write(*,201) cnpyfile
    write(*,203) cm_npts, npts
    close(UCNPY)
    stop
  else if (cm_z0 .ne. z0) then
    write(*,201) cnpyfile
    write(*,204) cm_z0, z0
    close(UCNPY)
    stop
  else if (cm_zi .ne. zi) then
    write(*,201) cnpyfile
    write(*,205) cm_zi, zi
    close(UCNPY)
    stop
  else if (cm_hc .ne. hc) then
    write(*,201) cnpyfile
    write(*,206) cm_hc, hc
    close(UCNPY)
    stop
  else if (cm_alfa .ne. alfa) then
    write(*,201) cnpyfile
    write(*,207) cm_alfa, alfa
    close(UCNPY)
    stop
  end if 

  ! passed checks, now read canopy data

  ! leaf angle distribution parameter
  read(UCNPY,*)
  read(UCNPY,*)
  read(UCNPY,*)   x
  read(UCNPY,*)
  read(UCNPY,*)

  ! diffuse radiation extinction coefficient
  read(UCNPY,*)   kd 
  read(UCNPY,*)
  read(UCNPY,*)

  ! g0, Medlyn et al. stomatal conductance parameter
  read(UCNPY,*)   g0 
  read(UCNPY,*)
  read(UCNPY,*)
 
  ! g1, Medlyn et al. stomatal conductance parameter
  read(UCNPY,*)   g1 
  read(UCNPY,*)
  read(UCNPY,*)
 
  ! dleaf, characteristic leaf dimension
  read(UCNPY,*)   dleaf
  read(UCNPY,*)
  read(UCNPY,*)
  read(UCNPY,*)

  ! zero arrays
  do i=1,npts
    lad(i) =0.0
    lai(i) =0.0
    clai(i)=0.0
  end do

  ! read canopy lad profile
  laitot=0.0
  do i=ncnpy,1,-1
    read(UCNPY,*) lad(i)
    lai(i)=lad(i)*dzhc          ! within canopy grid resolution is always constant!
    laitot=laitot+lad(i)*dzhc
    clai(i)=laitot
  end do
  close(UCNPY)

201 format('***Canopy data file ', a, ' is inconsistent with current ACCESS configuration!')
203 format('***Canopy npts = ', i4 /'***ACCESS npts = ', i4)
204 format('***Canopy z0 = ', e12.4 /'***ACCESS z0 = ', e12.4)
205 format('***Canopy zi = ', e12.4 /'***ACCESS zi = ', e12.4)
206 format('***Canopy hc = ', e12.4 /'***ACCESS hc = ', e12.4)
207 format('***Canopy alfa = ', e12.4 /'***ACCESS alfa = ', e12.4)
  return

end subroutine GetCanopyData

end module Initialize
