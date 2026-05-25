subroutine water_vapor_18o(hna,mna,ena,esa,msa,mnp,enp,esp,msp,ptan,ptas,taen,&
                           o18a,o18p,o18_z)
! Calculates seawater d18O ratio for all ocean boxes
! Inputs:
! Ocean tarcers
! Factors for meridional temperature atmospheric profile
! Outputs:
! o18a,o18p,o18_z: seawater d18O ratio

use parameters, only: wp,nto,n,fbr,fdh,fdl,nob

implicit none
real(wp), dimension(nto,n), intent(in) :: hna,mna,ena,esa,msa,mnp,enp,esp,msp
real(wp), intent(in) :: ptan(3),ptas(3),taen
real(wp), intent(out) :: o18a(nob),o18p(nob-1),o18_z
real(wp), external :: ta_at_phi
! Atlantic
o18a(1) = wavao18(hna(1,1),hna(10,1),ptan,fbr)
o18a(2) = wavao18(mna(1,1),mna(10,1),ptan,fdh)
o18a(3) = wavao18(ena(1,1),ena(10,1),ptan,fdl)
o18a(4) = wavao18(esa(1,1),esa(10,1),ptas,fdl)
o18a(5) = wavao18(msa(1,1),msa(10,1),ptas,fdh)

! Pacific
o18p(1) = wavao18(mnp(1,1),mnp(10,1),ptan,fdh)
o18p(2) = wavao18(enp(1,1),enp(10,1),ptan,fdl)
o18p(3) = wavao18(esp(1,1),esp(10,1),ptas,fdl)
o18p(4) = wavao18(msp(1,1),msp(10,1),ptas,fdh)

! Zonal transport ena-enp
o18_z = wavao18_zon(ena(1,1),ena(10,1),taen)


contains
real(wp) function wavao18(tw,ow,pta,phi) result(ao18c)
implicit none
real(wp), intent(in) :: tw,ow,pta(3),phi
real(wp) :: td,h,k18,alp18w,alp18d,alp18c,eww,ewd,fv,ao18w
! Water Vapour 18O
! tw: SST
! ow: d18O water
! ptaa: values to make Ta(phi) 
! tw and ow from evaorative zone

td     = ta_at_phi(pta,phi)
h      = .75_wp
k18    = .006_wp
alp18w = (.9884_wp + 1.025e-4_wp*(tw) - 3.57e-7_wp*(tw)**2)**(-1)
alp18d = (.9884_wp + 1.025e-4_wp*(td) - 3.57e-7_wp*(td)**2)**(-1)
alp18c = .5_wp*(alp18w+alp18d)
eww    = 10**(9.4_wp - 2.35e3_wp/(tw+273._wp))
ewd    = 10**(9.4_wp - 2.35e3_wp/(td+273._wp))
fv     = ewd/eww
ao18w  = ( 1._wp/alp18w * (ow+1._wp)*(1._wp-k18)/(1._wp-h*k18)-1._wp )
ao18c  = fv**(alp18c-1._wp)*(1._wp+AO18w)-1._wp

end function wavao18



real(wp) function wavao18_zon(tw,ow,td) result(ao18c)
implicit none
real(wp), intent(in) :: tw,ow,td
real(wp) :: h,k18,alp18w,alp18d,alp18c,eww,ewd,fv,ao18w

h      = .75_wp
k18    = .006_wp
alp18w = (.9884_wp + 1.025e-4_wp*(tw) - 3.57e-7_wp*(tw)**2)**(-1)
alp18d = (.9884_wp + 1.025e-4_wp*(td) - 3.57e-7_wp*(td)**2)**(-1)
alp18c = .5_wp*(alp18w+alp18d)
eww    = 10**(9.4_wp - 2.35e3_wp/(tw+273._wp))
ewd    = 10**(9.4_wp - 2.35e3_wp/(td+273._wp))
fv     = ewd/eww
ao18w  = ( 1._wp/alp18w * (ow+1._wp)*(1._wp-k18)/(1._wp-h*k18)-1._wp )
ao18c  = fv**(alp18c-1._wp)*(1._wp+AO18w)-1._wp

end function wavao18_zon

end subroutine water_vapor_18o
