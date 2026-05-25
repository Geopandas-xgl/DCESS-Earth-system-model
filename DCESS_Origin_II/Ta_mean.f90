real(wp) function ta_mean(pta,fdhl,fdll) result(ta)
! Calculates mean atmospheric temperature between two latitudes
use parameters, only: wp
implicit none
real(wp), intent(in) :: pta(3),fdhl,fdll
real(wp) :: sl,sc,sq
! fdhl: fd high latitude
! fdll: fd low latitude

sl    = sin(fdhl)**1 - sin(fdll)**1
sc    = sin(fdhl)**3 - sin(fdll)**3
sq    = sin(fdhl)**5 - sin(fdll)**5

ta = 1._wp/(sin(fdhl)-sin(fdll))*&
        ( pta(1)*sl + 1._wp/2._wp*pta(2)*(sc - sl) + 1._wp/8._wp*pta(3)*( 7*sq - 10*sc + 3*sl ) )

end function ta_mean