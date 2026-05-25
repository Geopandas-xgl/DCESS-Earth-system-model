subroutine land_new_prod(t,np_t0,np_g0,np_b0,fac_t,fac_g,fac_b,pco2,np_t,np_g,np_b,npp)
! Calculates net primary production on land as function of co2 and light

use parameters, only: wp,co2fer,pco2int
implicit none
real(wp), intent(in) :: t,np_t0,np_g0,np_b0,fac_t,fac_g,fac_b,pco2
real(wp), intent(out) :: np_t,np_g,np_b,npp(3)
real(wp) :: light_t,light_g,light_b

call light_time(t,light_t,light_g,light_b)

! New production for the 3 vegetation types.
np_t = light_t*np_t0*fac_t*(1._wp + co2fer*log(pco2/pco2int))       ! Primary production tropical forest
np_g = light_g*np_g0*fac_g*(1._wp + co2fer*log(pco2/pco2int))       ! Primary production grassland sd
np_b = light_b*np_b0*fac_b*(1._wp + co2fer*log(pco2/pco2int))       ! Primary production boreal forest

npp = (/np_t,np_g,np_b/)

contains
subroutine light_time(t,funtf,fungs,funbf)
implicit none
real(wp), intent(in) :: t
real(wp), intent(out) :: funtf,fungs,funbf
real(wp) :: a,b,c
! TF
a     = 0.4864_wp
b     = 1.672e7_wp
c     = 16590000._wp
funtf = 0.62_wp + a*exp(-((t-b)/c)**2)

! GSD
a     = 0.456_wp
b     = 1.672e7_wp
c     = 8295000._wp
fungs = 0.78_wp + a*exp(-((t-b)/c)**2)

! BF
a     = 2.62746_wp
b     = 1.672e7_wp
c     = 5574240._wp
funbf = 0.12_wp + a*exp(-((t-b)/c)**2)


end subroutine light_time
end subroutine land_new_prod