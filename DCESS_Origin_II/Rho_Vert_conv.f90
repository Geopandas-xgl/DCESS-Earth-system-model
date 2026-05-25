subroutine rho_vert_conv(psi,pres,rho,kv)
! Calculates vertical diffusion (stability-dependent)
! Inputs: 
! psi: ocen tracers (see OdeExp for data structure)
! pres: ocean pressure
! Outputs:
! rho: ocean density (kg/m3)
! kv: vertical diffusion profile (m2/s)

use parameters, only: wp,nto,n,d,g,rho0
implicit none
real(wp), intent(in) :: psi(nto,n), pres(n)
real(wp), intent(out) :: rho(n), kv(n-1)
real(wp) :: kv0,bvf0,alp,kmax
real(wp), dimension(n-1) :: drdz,bvf

call poly_teos10_bsq_z(psi(2,:),psi(1,:),pres,n,rho)
rho = rho - 1.e3_wp


! Scale values
kv0  = 3.e-5_wp
bvf0 = 1.e-2_wp
alp  = .5_wp
kmax = 0.8e-3_wp         ! Max vertical diffusivity

drdz =  1/d*(rho(1:n-1) - rho(2:n))

bvf = merge(0._wp, (-g / rho0 * drdz)**0.5, drdz >= 0._wp)

kv = merge(kmax,kv0*bvf0**alp/bvf**alp,drdz>=0._wp)
kv = min(kv,kmax)


end subroutine rho_vert_conv

include 'Poly_TEOS10_bsq_z.f90'