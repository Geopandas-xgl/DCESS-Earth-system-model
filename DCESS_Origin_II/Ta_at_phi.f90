real(wp) function ta_at_phi(pta,phi) result(ta_phi)
! Calculates atmospheric temperature at specific latitude phi
use parameters, only: wp
implicit none
real(wp), intent(in) :: pta(3), phi
ta_phi = pta(1) + pta(2)*(1._wp/2._wp)*(3._wp*sin(phi)**2-1._wp)+&
         pta(3)*(1._wp/8._wp)*(35._wp*sin(phi)**4-30._wp*sin(phi)**2+3._wp)
end function