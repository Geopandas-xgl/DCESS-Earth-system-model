subroutine entrainment(sa_shl,buo0,irq0,rds0,alfi,tsoi,ssoi,pds,z,rsoi,fdsup,fe,ent_zlimits)
! Calculates entrainment flux for the downsloping flow from Southern Ocean shelf

! Inputs:
! sa_shl: shelf salinity
! buo0: buoyancy at the shelf break
! irq0: Initial influx rate
! rds0: density at the shelf break
! alfi: Antarctic slope
! tsoi,ssoi,rsoi: temperature, salinity and pressure of Southern Ocean
! Outputs:
! fdsup: flux leaving the shelf (Sv)
! fe: Entrainment flux (Sv)
! ent_zlimits: enter level of plume and plume height in z-levels
 
use parameters, only: wp,nz_upp,nz_bot,nz_all,n,hgt0,to_f,shl_dz,dr,rho0,g,d,shl
implicit none
real(wp), intent(in) :: sa_shl,buo0,irq0,rds0(nz_upp),alfi(nz_bot),tsoi(nz_bot),ssoi(nz_bot),pds(nz_bot) ! shelf salinity
real(wp), intent(in) :: z(nz_bot),rsoi(nz_all)
real(wp), intent(out) :: fdsup,fe(n)
integer, intent(out) :: ent_zlimits(2)
real(wp), dimension(nz_bot) :: irq
real(wp) :: rds(nz_all),buo,hgt,tds,sds,eps,ds,ri,ec,dirq,dhgt,rdsaux,i0,limz,deltaq(n)
integer :: i,ip,ii,zlim,i1,i2,nzin
real(wp), parameter :: cd    = 0.0020_wp,&          ! Drag coefficient Baines 2008, Table 1
                       s2    = 0.80_wp,&
                       f     = 5.4_wp/5.6766_wp*5.4_wp/6.569_wp*1.7_wp*15._wp/25._wp,&
                       ls    = f*0.85e5_wp          ! Length scale to convert to Sverdrups (10^6 m3/s)

irq = 0._wp
rds = 0._wp

buo           = buo0
hgt           = hgt0
irq(1)        = irq0
tds           = to_f
sds           = sa_shl
rds(1:nz_upp) = rds0
ip            = nz_upp-1
eps           = 0.00004_wp

do i=2,nz_bot
        ds = shl_dz/tan(alfi(i)*dr)                                     ! scaled horizontal step 
        ri = buo*hgt**3*cos(alfi(i)*dr)/irq(i-1)**2                     ! Richardson number
        ec  = merge(0.2_wp*(1._wp-ri/0.25_wp) ,0._wp,ri <= 0.25_wp)     ! entrainment coefficient
        
        dirq = (ec*irq(i-1)/hgt)*ds
        dhgt = (2._wp*ec+cd-s2*ri*tan(alfi(i)*dr))*ds

        tds = ( irq(i-1)*tds + dirq*tsoi(i-1) )/(irq(i-1)+dirq)
        sds = ( irq(i-1)*sds + dirq*ssoi(i-1) )/(irq(i-1)+dirq)
        
        call poly_teos10_bsq(sds,tds,pds(i),rdsaux)
        rds(i+ip) = rdsaux
        irq(i)    = irq(i-1) + dirq
        hgt       = hgt      + dhgt

        i0  = ( z(i) - hgt/2._wp )/shl_dz
        ii  = floor(i0)

        buo = g*(rds(ii)-rsoi(ii))/rho0 ! buoyancy

        if (buo <= eps) then
        limz=z(i)
        exit
        else
        limz=z(i)
        end if

end do
i=i-1

irq(i+1:nz_bot) = irq(i)
zlim            = ceiling(limz/d)

do i=shl+1,zlim
        i1 = findloc(z .EQ.  i   *d,.TRUE.,1)
        i2 = findloc(z .EQ. (i-1)*d,.TRUE.,1)
        deltaq(i) = irq(i1) - irq(i2)
end do

fdsup        = irq0*ls                          ! Upper layer flux (m3/s)
fe(1:n)      = 0._wp                            ! Entrainment flux (m3/s)
fe(shl+1:zlim) = deltaq(shl+1:zlim)*ls


nzin        = ceiling(hgt/d)
ent_zlimits = (/zlim,nzin/)

end subroutine entrainment

include 'Poly_TEOS10_bsq.f90'