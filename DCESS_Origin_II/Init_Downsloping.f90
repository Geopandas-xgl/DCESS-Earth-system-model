subroutine init_downsloping(soc,slso,shl_dat,z,alfi,buo0,irq0,ctoi,saoi,rsoi,r_ds0,p_ds)
! Calculates initial conditions for the downsloping flow
! Inputs:
! soc: Southern Ocean tracers
! slso: Shelf tracers
! shl_dat: Antarctic slope data
! Outputs:
! alfi: interpolated Antarctic slope data
! buo0: buoyancy at the shelf break
! irq0: Initial influx rate
! ctoi: Interpolated Southern Ocean temperature
! saoi: Interpolated Southern Ocean salinity
! rsoi: Interpolated Southern Ocean density
! r_ds0: density at the shelf break
! p_ds: pressure at shelf

use parameters, only: wp,nto,n,nz_all,nz_bot,nz_upp,nz_mid,d,shl,zmid,shl_dz,rho0,hgt0,to_f,g
use :: lib_array, only: interp1d
implicit none
real(wp), intent(in) :: soc(nto,n), slso(nto),shl_dat(nz_all,3)
real(wp), dimension(nz_bot), intent(out) :: z,ctoi,saoi,alfi,p_ds
real(wp), intent(out) :: rsoi(nz_all),buo0,irq0,r_ds0(nz_upp)

real(wp) :: shlz,zini,zend,rso_shl,alfix(nz_all),pshi(nz_all),psoi(nz_all)
real(wp), dimension(nz_all) :: zi,ctoix,saoix,rsoix
real(wp) :: zic(nz_mid),p(nz_upp),s_shv(nz_upp),t_shv(nz_upp)
integer :: i,ii,in
logical :: msk_iz(nz_all),msk_rso(nz_all),msk_shl(nz_all)

shlz  = d*shl                       ! shelf level (meters)    

zi  = (/(real(i,wp), i=0,int(n*d),int(shl_dz))/)
zic = (/(real(i,wp), i=int(zmid(1)),int(zmid(n)),int(shl_dz))/)

ii = findloc(zi .EQ. zmid(1), .TRUE.,1)
in = findloc(zi .EQ. zmid(n), .TRUE.,1)

ctoix(:) = 0._wp
saoix(:) = 0._wp

ctoix(ii:in)        = interp1d(zmid,soc(1,:),zic)
saoix(ii:in)        = interp1d(zmid,soc(2,:),zic)
ctoix(1      :ii-1) = soc(1,1)
ctoix(in+1:nz_all ) = soc(1,n)
saoix(1      :ii-1) = soc(2,1)
saoix(in+1:nz_all ) = soc(2,n)

alfix = shl_dat(:,1)
pshi  = shl_dat(:,2)
psoi  = shl_dat(:,3)

call poly_teos10_bsq_z(saoix,ctoix,psoi,nz_all,rsoix)

zini   = shlz
zend   = n*d
msk_iz = ( zi>=zini .AND. zi<=zend )
z      = pack(zi   ,msk_iz)
alfi   = pack(alfix,msk_iz)
ctoi   = pack(ctoix,msk_iz)
saoi   = pack(saoix,msk_iz)

msk_rso = ( zi<=zend )
rsoi = pack(rsoix,msk_rso)

! SO rho at at shelf break
ii = findloc(zi .EQ. zini       , .TRUE.,1)
rso_shl = rsoi(ii)

msk_shl = ( zi <= shlz )

p = pack(pshi,msk_shl)

s_shv = slso(2)
t_shv = to_f

call poly_teos10_bsq_z(s_shv,t_shv,p,nz_upp,r_ds0)

p_ds = pack(pshi,msk_iz)

buo0 = g*(r_ds0(size(r_ds0))-rso_shl)/rho0      ! Buoyancy (G)
if (buo0 < 0._wp) then
        print*, 'rho_shelf < rho_SO'
end if
irq0 = hgt0*(buo0*hgt0)**0.5                    ! influx rate (Q)

end subroutine init_downsloping