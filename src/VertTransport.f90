!======================================================================================================================!
!                                                                                                                      !
!     Program:      ACCESS                                                                                             !
!                   Atmospheric Chemistry and Canopy Exchange Simulation System                                        !
!                   Full BVOC chemistry version                                                                        !
!                                                                                                                      !
!     Version:      3.0.0                                                                                              !
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
!     Module:       VertTransport                                                                                      !
!                                                                                                                      !
!     Description:  Contains routines for integrating the diffusion equation in the vertical domain                    !
!                                                                                                                      !
!======================================================================================================================!
module VertTransport
  use GlobalData
  use Utils
  implicit none

  public SubIntegScalarConstBC, SubIntegScalarFluxBC, SubIntegScalarMixedBCa, SubIntegScalarMixedBCb, &
          SubIntegTempConstBC, SubIntegTempMixedBCb

contains

!**********************************************************************************************************************!
! subroutine  SubIntegScalarConstBC - integrates the tracer diffusion equation over time step "dttot"
!                                     using nvxtot timesteps of dtstep timestep
!                                     Constant BCs at surface and domain top
!
!     Follows the recommendation of Venkatram (1993) to use the mass conservative
!           form of the vertical transport equation for tracers.
!
!     Venkatram, A. (1993) The parameterization of the vertical dispersion of a 
!      scalar in the atmospheric boundary layer, Atmos. Environ., 27A, 1963-1966. 
!
!     Discretization follows the Crank-Nicolson treatment of the 1-D diffusion (conduction) equation of
!      Patankar (1980).
!
!     Patankar, S. V. (1980) Numerical Heat Transfer and Fluid Flow, Hemisphere Publishing Co., Washington, D.C.
!**********************************************************************************************************************!
subroutine SubIntegScalarConstBC(phi, phip, kv, xa, xs, rho, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip, rho
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: xa, xs
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: rho2m1, rho2p1, rhoim1, rhoip1, rhohm1, rhohp1
  real(kind=dp)                    :: k2m1, k2p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dz2m1, dz2p1, dzim1, dzip1, dzhm1, dzhp1
  real(kind=dp)                    :: dz2, dzi, dznm1
  real(kind=dp), dimension(npts-2) :: a, b, c, r
  real(kind=dp), dimension(npts-2) :: x
  integer(kind=i4) :: nvx, nvxtot, i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !
    !  Surface
    rho2m1 = 0.5*(rho(2)+rho(1))
    rho2p1 = 0.5*(rho(3)+rho(2))
    k2m1   = 0.5*(kv(2)+kv(1))
    k2p1   = 0.5*(kv(3)+kv(2))
    dz2m1  = z(2)-z(1)   
    dz2p1  = z(3)-z(2)
    dz2    = 0.5*(z(3)-z(1))
    a(1) = 0.0
    b(1) = (0.5*rho2m1*k2m1/(rho(2)*dz2m1)) + (dz2/dtstep) + ( 0.5*rho2p1*k2p1/(rho(2)*dz2p1))
    c(1) = -0.5*rho2p1*k2p1/(rho(3)*dz2p1)
    r(1) = rho2m1*k2m1*xs/(rho(1)*dz2m1) +  &
             (-0.5*rho2m1*k2m1/(rho(2)*dz2m1) + (dz2/dtstep) -0.5*rho2p1*k2p1/(rho(2)*dz2p1))*phip(2) &
           + 0.5*rho2p1*k2p1*phip(3)/(rho(3)*dz2p1)

    !  Interior points
    !   "k" is the tridiagonal matrix row
    !   "i" is the domain vertical level
    do k=2,npts-3
      i = k+1
      rhoim1 = 0.5*(rho(i)+rho(i-1))
      rhoip1 = 0.5*(rho(i+1)+rho(i))
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*rhoim1*kim1/(rho(i-1)*dzim1)
      b(k) = 0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) + 0.5*rhoip1*kip1/(rho(i)*dzip1)
      c(k) = -0.5*rhoip1*kip1/(rho(i+1)*dzip1)
      r(k) = 0.5*rhoim1*kim1*phip(i-1)/(rho(i-1)*dzim1) +  &
             (-0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) - 0.5*rhoip1*kip1/(rho(i)*dzip1))*phip(i)  &
             + 0.5*rhoip1*kip1*phip(i+1)/(rho(i+1)*dzip1)
    end do 

    !  Domain top
    rhohm1 = 0.5*(rho(npts-1)+rho(npts-2))
    rhohp1 = 0.5*(rho(npts)+rho(npts-1))
    khm1   = 0.5*(kv(npts-1)+kv(npts-2))
    khp1   = 0.5*(kv(npts)+kv(npts-1))
    dzhm1  = z(npts-1)-z(npts-2)
    dzhp1  = z(npts)-z(npts-1)
    dznm1  = 0.5*(z(npts)-z(npts-2))
    a(npts-2) = -0.5*rhohm1*khm1/(rho(npts-2)*dzhm1)
    b(npts-2) = 0.5*rhohm1*khm1/(rho(npts-1)*dzhm1) + (dznm1/dtstep) + 0.5*rhohp1*khp1/(rho(npts-1)*dzhp1)
    c(npts-2) = 0.0
    r(npts-2) = 0.5*rhohm1*khm1*phip(npts-2)/(rho(npts-2)*dzhm1) +  &
              (-0.5*rhohm1*khm1/(rho(npts-1)*dzhm1)+(dznm1/dtstep)-0.5*rhohp1*khp1/(rho(npts-1)*dzhp1))*phip(npts-1) &
              + rhohp1*khp1*xa/(rho(npts)*dzhp1)

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-2)

    do k=1,npts-2
      i = k+1
      phi(i) = x(k)
    end do

    phi(1)    = xs
    phi(npts) = xa 

    phip = phi

  end do

  return
end subroutine SubIntegScalarConstBC


!**********************************************************************************************************************!
! subroutine  SubIntegScalarMixedBCa - integrates the diffusion equation for dttot
!                                           using nvxtot timesteps of dtstep timestep
!                                           Constant BC at domain top, Flux BC at surface
!
!     Follows the recommendation of Venkatram (1993) to use the mass conservative
!           form of the vertical transport equation.
!
!     Venkatram, A. (1993) The parameterization of the vertical dispersion of a 
!      scalar in the atmospheric boundary layer, Atmos. Environ., 27A, 1963-1966. 
!
!     Discretization follows the Crank-Nicolson treatment of the 1-D diffusion (conduction) equation of
!      Patankar (1980).
!
!     Patankar, S. V. (1980) Numerical Heat Transfer and Fluid Flow, Hemisphere Publishing Co., Washington, D.C.
!**********************************************************************************************************************!
subroutine SubIntegScalarMixedBCa(phi, phip, kv, xa, xs, vs, rho, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip, rho
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: xa, xs, vs, alf
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: rho1p1, rhoim1, rhoip1, rhohm1, rhohp1
  real(kind=dp)                    :: k1p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dzim1, dzip1, dzhm1, dzhp1, dznm1
  real(kind=dp)                    :: dz0, dz02, dz2, dzi
  real(kind=dp), dimension(npts-1) :: a, b, c, r
  real(kind=dp), dimension(npts-1) :: x
  integer(kind=i4)                 :: nvx, nvxtot
  integer(kind=i4)                 :: i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !
    !  Surface
    rho1p1 = 0.5*(rho(2)+rho(1))
    k1p1   = 0.5*(kv(2)+kv(1))
    dz0    = z(2)-z(1)
    dz02   = 2*dz0
    alf    = vs*dz02/(rho(1)*kv(1))
    a(1) = 0.0
    b(1) = (0.5*kv(1)/dz0)*(1.0+alf*rho(1)) + (2.0*dz0/dtstep) + (0.5*rho1p1*k1p1/(rho(1)*dz0))
    c(1) = -0.5*rho(1)*kv(1)/(rho(2)*dz0) - 0.5*rho1p1*k1p1/(rho(2)*dz0)
    r(1) = rho(1)*kv(1)*alf*xs/dz0 +  &
             ((-0.5*kv(1)/dz0)*(1.0+alf*rho(1)) + (2.0*dz0/dtstep) - 0.5*rho1p1*k1p1/(rho(1)*dz0))*phip(1) &
           + (0.5*rho(1)*kv(1)/(rho(2)*dz0) + 0.5*rho1p1*k1p1/(rho(2)*dz0))*phip(2)

    !  Interior points
    !   "k" is the tridiagonal matrix row
    !   "i" is the domain vertical level
    do k=2,npts-2
      i = k
      rhoim1 = 0.5*(rho(i)+rho(i-1))
      rhoip1 = 0.5*(rho(i+1)+rho(i))
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*rhoim1*kim1/(rho(i-1)*dzim1)
      b(k) = 0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) + 0.5*rhoip1*kip1/(rho(i)*dzip1)
      c(k) = -0.5*rhoip1*kip1/(rho(i+1)*dzip1)
      r(k) = 0.5*rhoim1*kim1*phip(i-1)/(rho(i-1)*dzim1) +  &
             (-0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) - 0.5*rhoip1*kip1/(rho(i)*dzip1))*phip(i)  &
             + 0.5*rhoip1*kip1*phip(i+1)/(rho(i+1)*dzip1)
    end do 

    !  Domain top
    rhohm1 = 0.5*(rho(npts-1)+rho(npts-2))
    rhohp1 = 0.5*(rho(npts)+rho(npts-1))
    khm1   = 0.5*(kv(npts-1)+kv(npts-2))
    khp1   = 0.5*(kv(npts)+kv(npts-1))
    dzhm1  = z(npts-1)-z(npts-2)
    dzhp1  = z(npts)-z(npts-1)
    dznm1  = 0.5*(z(npts)-z(npts-2))
    a(npts-1) = -0.5*rhohm1*khm1/(rho(npts-2)*dzhm1)
    b(npts-1) = 0.5*rhohm1*khm1/(rho(npts-1)*dzhm1) + (dznm1/dtstep) + 0.5*rhohp1*khp1/(rho(npts-1)*dzhp1)
    c(npts-1) = 0.0
    r(npts-1) = 0.5*rhohm1*khm1*phip(npts-2)/(rho(npts-2)*dzhm1) +  &
              (-0.5*rhohm1*khm1/(rho(npts-1)*dzhm1)+(dznm1/dtstep)-0.5*rhohp1*khp1/(rho(npts-1)*dzhp1))*phip(npts-1) &
              + rhohp1*khp1*xa/(rho(npts)*dzhp1)

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-1)

    do k=1,npts-1
      i = k
      phi(i) = x(k)
    end do

    phi(npts) = xa 

    phip = phi

  end do

  return
end subroutine SubIntegScalarMixedBCa

!**********************************************************************************************************************!
! subroutine  SubIntegScalarFluxBC  - integrates the tracer diffusion equation over a timestep "dttot"
!                                     using nvxtot timesteps of dtstep timestep
!                                     Flux BC at domain top, Flux BC at surface
!
!     Follows the recommendation of Venkatram (1993) to use the mass conservative
!           form of the vertical transport equation.
!
!     Venkatram, A. (1993) The parameterization of the vertical dispersion of a 
!      scalar in the atmospheric boundary layer, Atmos. Environ., 27A, 1963-1966. 
!
!     Discretization follows the Crank-Nicolson treatment of the 1-D diffusion (conduction) equation of
!      Patankar (1980).
!
!     Patankar, S. V. (1980) Numerical Heat Transfer and Fluid Flow, Hemisphere Publishing Co., Washington, D.C.
!**********************************************************************************************************************!
subroutine SubIntegScalarFluxBC(phi, phip, kv, xs, vs, vflux, rho, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip, rho
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: xs, vs, vflux, alf
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: rho1p1, rhoim1, rhoip1, rhohm1, rhohp1
  real(kind=dp)                    :: k1p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dzim1, dzip1, dzh, dznm1
  real(kind=dp)                    :: dz0, dz02, dz2, dzi
  real(kind=dp), dimension(npts)   :: a, b, c, r
  real(kind=dp), dimension(npts)   :: x
  integer(kind=i4)                 :: nvx, nvxtot
  integer(kind=i4)                 :: i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !
    !  Surface
    rho1p1 = 0.5*(rho(2)+rho(1))
    k1p1   = 0.5*(kv(2)+kv(1))
    dz0    = z(2)-z(1)
    dz02   = 2*dz0
    alf    = vs*dz02/(rho(1)*kv(1))
    a(1) = 0.0
    b(1) = (0.5*kv(1)/dz0)*(1.0+alf*rho(1)) + (2.0*dz0/dtstep) + (0.5*rho1p1*k1p1/(rho(1)*dz0))
    c(1) = -0.5*rho(1)*kv(1)/(rho(2)*dz0) - 0.5*rho1p1*k1p1/(rho(2)*dz0)
    r(1) = rho(1)*kv(1)*alf*xs/dz0 +  &
             ((-0.5*kv(1)/dz0)*(1.0+alf*rho(1)) + (2.0*dz0/dtstep) - 0.5*rho1p1*k1p1/(rho(1)*dz0))*phip(1) &
           + (0.5*rho(1)*kv(1)/(rho(2)*dz0) + 0.5*rho1p1*k1p1/(rho(2)*dz0))*phip(2)

    !  Interior points
    !   "k" is the tridiagonal matrix row
    !   "i" is the domain vertical level
    do k=2,npts-1
      i = k
      rhoim1 = 0.5*(rho(i)+rho(i-1))
      rhoip1 = 0.5*(rho(i+1)+rho(i))
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*rhoim1*kim1/(rho(i-1)*dzim1)
      b(k) = 0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) + 0.5*rhoip1*kip1/(rho(i)*dzip1)
      c(k) = -0.5*rhoip1*kip1/(rho(i+1)*dzip1)
      r(k) = 0.5*rhoim1*kim1*phip(i-1)/(rho(i-1)*dzim1) +  &
             (-0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) - 0.5*rhoip1*kip1/(rho(i)*dzip1))*phip(i)  &
             + 0.5*rhoip1*kip1*phip(i+1)/(rho(i+1)*dzip1)
    end do 

    !  Domain top
    rhohm1 = 0.5*(rho(npts)+rho(npts-1))
    rhohp1 = rho(npts)
    khm1   = 0.5*(kv(npts)+kv(npts-1))
    khp1   = kv(npts)
    dzh  = z(npts)-z(npts-1)
    a(npts) = -(0.5*rhohm1*khm1/(rho(npts-1)*dzh)) - (0.5*rho(npts)*kv(npts)/(rho(npts)*dzh))
    b(npts) = 0.5*rhohm1*khm1/(rho(npts)*dzh) + (dzh/dtstep) + 0.5*rhohp1*khp1/(rho(npts)*dzh)
    c(npts) = 0.0
    r(npts) = (0.5*rhohm1*khm1/(rho(npts-1)*dzh) + 0.5*rhohp1*khp1/(rho(npts)*dzh))*phip(npts-1) +  &
              (-0.5*rhohm1*khm1/(rho(npts)*dzh)+(dzh/dtstep)-0.5*rhohp1*khp1/(rho(npts)*dzh))*phip(npts) &
              - 2.0*khp1*vflux

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts)

    do k=1,npts
      i = k
      phi(i) = x(k)
    end do

    phip = phi

  end do

  return
end subroutine SubIntegScalarFluxBC

!**********************************************************************************************************************!
! subroutine  SubIntegScalarMixedBCb - integrates the tracer diffusion equation over time step "dttot"
!                                      using nvxtot timesteps of dtstep timestep
!                                      Constant BC at the surface and flux BC at the domain top
!
!     Follows the recommendation of Venkatram (1993) to use the mass conservative
!           form of the vertical transport equation.
!
!     Venkatram, A. (1993) The parameterization of the vertical dispersion of a 
!      scalar in the atmospheric boundary layer, Atmos. Environ., 27A, 1963-1966. 
!
!     Discretization follows the Crank-Nicolson treatment of the 1-D diffusion (conduction) equation of
!      Patankar (1980).
!
!     Patankar, S. V. (1980) Numerical Heat Transfer and Fluid Flow, Hemisphere Publishing Co., Washington, D.C.
!**********************************************************************************************************************!
subroutine SubIntegScalarMixedBCb(phi, phip, kv, vflux, xs, rho, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip, rho
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: vflux, xs
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: rho2m1, rho2p1, rhoim1, rhoip1, rhohm1, rhohp1
  real(kind=dp)                    :: k2m1, k2p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dz2m1, dz2p1, dzim1, dzip1, dzh
  real(kind=dp)                    :: dz2, dzi, dznm1
  real(kind=dp), dimension(npts-1) :: a, b, c, r
  real(kind=dp), dimension(npts-1) :: x
  integer(kind=i4) :: nvx, nvxtot, i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !
    !  Surface
    rho2m1 = 0.5*(rho(2)+rho(1))
    rho2p1 = 0.5*(rho(3)+rho(2))
    k2m1   = 0.5*(kv(2)+kv(1))
    k2p1   = 0.5*(kv(3)+kv(2))
    dz2m1  = z(2)-z(1)   
    dz2p1  = z(3)-z(2)
    dz2    = 0.5*(z(3)-z(1))
    a(1) = 0.0
    b(1) = (0.5*rho2m1*k2m1/(rho(2)*dz2m1)) + (dz2/dtstep) + ( 0.5*rho2p1*k2p1/(rho(2)*dz2p1))
    c(1) = -0.5*rho2p1*k2p1/(rho(3)*dz2p1)
    r(1) = rho2m1*k2m1*xs/(rho(1)*dz2m1) +  &
             (-0.5*rho2m1*k2m1/(rho(2)*dz2m1) + (dz2/dtstep) -0.5*rho2p1*k2p1/(rho(2)*dz2p1))*phip(2) &
           + 0.5*rho2p1*k2p1*phip(3)/(rho(3)*dz2p1)

    !  Interior points
    !   "k" is the tridiagonal matrix row
    !   "i" is the domain vertical level 
    do k=2,npts-2
      i = k+1
      rhoim1 = 0.5*(rho(i)+rho(i-1))
      rhoip1 = 0.5*(rho(i+1)+rho(i))
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*rhoim1*kim1/(rho(i-1)*dzim1)
      b(k) = 0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) + 0.5*rhoip1*kip1/(rho(i)*dzip1)
      c(k) = -0.5*rhoip1*kip1/(rho(i+1)*dzip1)
      r(k) = 0.5*rhoim1*kim1*phip(i-1)/(rho(i-1)*dzim1) +  &
             (-0.5*rhoim1*kim1/(rho(i)*dzim1) + (dzi/dtstep) - 0.5*rhoip1*kip1/(rho(i)*dzip1))*phip(i)  &
             + 0.5*rhoip1*kip1*phip(i+1)/(rho(i+1)*dzip1)
    end do 

    !  Domain top
    rhohm1 = 0.5*(rho(npts)+rho(npts-1))
    rhohp1 = rho(npts)
    khm1   = 0.5*(kv(npts)+kv(npts-1))
    khp1   = kv(npts)
    dzh  = z(npts)-z(npts-1)
    a(npts-1) = -(0.5*rhohm1*khm1/(rho(npts-1)*dzh)) - (0.5*rho(npts)*kv(npts)/(rho(npts)*dzh))
    b(npts-1) = 0.5*rhohm1*khm1/(rho(npts)*dzh) + (dzh/dtstep) + 0.5*rhohp1*khp1/(rho(npts)*dzh)
    c(npts-1) = 0.0
    r(npts-1) = (0.5*rhohm1*khm1/(rho(npts-1)*dzh) + 0.5*rhohp1*khp1/(rho(npts)*dzh))*phip(npts-1) +  &
              (-0.5*rhohm1*khm1/(rho(npts)*dzh)+(dzh/dtstep)-0.5*rhohp1*khp1/(rho(npts)*dzh))*phip(npts) &
              - 2.0*khp1*vflux

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-1)

    do k=1,npts-1
      i = k+1
      phi(k+1) = x(k)
    end do

    phi(1)    = xs

    phip = phi

  end do

  return
end subroutine SubIntegScalarMixedBCb

!**********************************************************************************************************************!
! subroutine  SubIntegTempConstBC - integrates the heat equation over time step "dttot"
!                                    using nvxtot timesteps of dtstep timestep
!                                    Constant BCs at surface and domain top
!
!     Discretization follows the Crank-Nicolson treatment of the 1-D diffusion (conduction) equation of
!      Patankar (1980).
!
!     Patankar, S. V. (1980) Numerical Heat Transfer and Fluid Flow, Hemisphere Publishing Co., Washington, D.C.
!
!**********************************************************************************************************************!
subroutine SubIntegTempConstBC(phi, phip, kv, xa, xs, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: xa, xs
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: k2m1, k2p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dz2m1, dz2p1, dzim1, dzip1, dzhm1, dzhp1
  real(kind=dp)                    :: dz2, dzi, dznm1
  real(kind=dp), dimension(npts-2) :: a, b, c, r
  real(kind=dp), dimension(npts-2) :: x
  integer(kind=i4) :: nvx, nvxtot, i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !  Surface
    k2m1   = 0.5*(kv(2)+kv(1))
    k2p1   = 0.5*(kv(3)+kv(2))
    dz2m1  = z(2)-z(1)   
    dz2p1  = z(3)-z(2)
    dz2    = 0.5*(z(3)-z(1))
    a(1) = 0.0
    b(1) = (0.5*k2m1/dz2m1) + (dz2/dtstep) + (0.5*k2p1/dz2p1)
    c(1) = -0.5*k2p1/dz2p1
    r(1) = (k2m1*xs/dz2m1) +  &
             ((-0.5*k2m1/dz2m1) + (dz2/dtstep) - (0.5*k2p1/dz2p1))*phip(2) + (0.5*k2p1/dz2p1)*phip(3)

    !  Interior points
    !    "k" is the tridiagonal matrix row
    !    "i" is the domain vertical level
    do k=2,npts-3
      i = k+1
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*kim1/dzim1
      b(k) = (0.5*kim1/dzim1) + (dzi/dtstep) + (0.5*kip1/dzip1)
      c(k) = -0.5*kip1/dzip1
      r(k) = (0.5*kim1/dzim1)*phip(i-1) +  &
             ( (-0.5*kim1/dzim1) + (dzi/dtstep) - (0.5*kip1/dzip1) )*phip(i) + (0.5*kip1/dzip1)*phip(i+1)
    end do 

    !  Domain top
    khm1   = 0.5*(kv(npts-1)+kv(npts-2))
    khp1   = 0.5*(kv(npts)+kv(npts-1))
    dzhm1  = z(npts-1)-z(npts-2)
    dzhp1  = z(npts)-z(npts-1)
    dznm1  = 0.5*(z(npts)-z(npts-2))
    a(npts-2) = -0.5*khm1/dzhm1
    b(npts-2) = (0.5*khm1/dzhm1) + (dznm1/dtstep) + (0.5*khp1/dzhp1)
    c(npts-2) = 0.0
    r(npts-2) = (0.5*khm1/dzhm1)*phip(npts-2) +  &
              ( (-0.5*khm1/dzhm1) + (dznm1/dtstep) - (0.5*khp1/dzhp1) )*phip(npts-1) + (khp1*xa/dzhp1)

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-2)

    do k=1,npts-2
      i = k+1
      phi(i) = x(k)
    end do

    phi(1)    = xs
    phi(npts) = xa 

    phip = phi

  end do

  return
end subroutine SubIntegTempConstBC


!**********************************************************************************************************************!
! subroutine  SubIntegTempMixedBCb - integrates the heat equation over "dttot"
!                                    using nvxtot timesteps of dtstep timestep
!                                    Constant BC at the surface and and flux BC at the domain top
!
!     Discretization follows the Crank-Nicolson treatment of the 1-D diffusion (conduction) equation of
!      Patankar (1980).
!
!     Patankar, S. V. (1980) Numerical Heat Transfer and Fluid Flow, Hemisphere Publishing Co., Washington, D.C.
!
!**********************************************************************************************************************!
subroutine SubIntegTempMixedBCb(phi, phip, kv, vflux, xs, dttot)
  real(kind=dp), dimension(npts)   :: phi, phip
  real(kind=dp), dimension(npts)   :: kv
  real(kind=dp)                    :: vflux, xs
  real(kind=dp)                    :: dttot, dtstep
  real(kind=dp)                    :: k2m1, k2p1, kim1, kip1, khm1, khp1
  real(kind=dp)                    :: dz2m1, dz2p1, dzim1, dzip1
  real(kind=dp)                    :: dz2, dzi, dzh
  real(kind=dp), dimension(npts-1) :: a, b, c, r
  real(kind=dp), dimension(npts-1) :: x
  integer(kind=i4) :: nvx, nvxtot, i, k

  dtstep = 0.25      ! vertical transport time step (s)

  ! number of dtstep time steps to take to integrate over dttot
  nvxtot = int(dttot/dtstep)

  do nvx=1,nvxtot
 
    ! Fill tridiagonal matrix and rhs
    !
    !  Surface
    k2m1   = 0.5*(kv(2)+kv(1))
    k2p1   = 0.5*(kv(3)+kv(2))
    dz2m1  = z(2)-z(1)   
    dz2p1  = z(3)-z(2)
    dz2    = 0.5*(z(3)-z(1))
    a(1) = 0.0
    b(1) = (0.5*k2m1/dz2m1) + (dz2/dtstep) + (0.5*k2p1/dz2p1)
    c(1) = -0.5*k2p1/dz2p1
    r(1) = (k2m1*xs/dz2m1) +  &
             ( (-0.5*k2m1/dz2m1) + (dz2/dtstep) - (0.5*k2p1/dz2p1) )*phip(2) + (0.5*k2p1/dz2p1)*phip(3)

    !  Interior points
    !  "k" is the tridiagonal matrix row
    !  "i" is the domain vertical level
    do k=2,npts-2
      i = k+1
      kim1   = 0.5*(kv(i)+kv(i-1))
      kip1   = 0.5*(kv(i+1)+kv(i))
      dzim1  = z(i)-z(i-1)
      dzip1  = z(i+1)-z(i)
      dzi    = 0.5*(z(i+1)-z(i-1))
      a(k) = -0.5*kim1/dzim1
      b(k) = (0.5*kim1/dzim1) + (dzi/dtstep) + (0.5*kip1/dzip1)
      c(k) = -0.5*kip1/dzip1
      r(k) = (0.5*kim1/dzim1)*phip(i-1) +  &
             ( (-0.5*kim1/dzim1) + (dzi/dtstep) - (0.5*kip1/dzip1) )*phip(i) + (0.5*kip1/dzip1)*phip(i+1)
    end do 

    !  Domain top
    khm1   = 0.5*(kv(npts)+kv(npts-1))
    khp1   = kv(npts)
    dzh  = z(npts)-z(npts-1)
    a(npts-1) = -(0.5*khm1/dzh) - (0.5*khp1/dzh)
    b(npts-1) = (0.5*khm1/dzh) + (dzh/dtstep) + (0.5*khp1/dzh)
    c(npts-1) = 0.0
    r(npts-1) = ( (0.5*khm1/dzh) + (0.5*khp1/dzh) )*phip(npts-1) +  &
              ( (-0.5*khm1/dzh) + (dzh/dtstep) - (0.5*khp1/dzh) )*phip(npts) - 2.0*khp1*vflux

    ! Call tridiagonal matrix solver
    call SolveTridiag(a, b, c, r, x, npts-1)

    do k=1,npts-1
      i = k+1
      phi(i) = x(k)
    end do

    phi(1)    = xs

    phip = phi

  end do

  return
end subroutine SubIntegTempMixedBCb

end module VertTransport
!======================================================================================================================!
