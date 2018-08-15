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
    lai(i)=lad(i)*dzhc             ! within canopy grid resolution is always constant!
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

