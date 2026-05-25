real(wp) function dtdy_at_phi(pta,phi) result(ty)
use parameters, only: wp
implicit none
real(wp), intent(in) :: pta(3),phi
ty = sin(phi)*cos(phi)*( 3._wp*pta(2) + (35._wp/2._wp)*pta(3)*sin(phi)**2-(15._wp/2._wp)*pta(3) )
end function dtdy_at_phi