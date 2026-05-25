function fun_snow(x,c,d,e) result(y)
use parameters, only: wp,tsnow
implicit none
real(wp), intent(in) :: x,c,d,e
real(wp) :: y

y = 35._wp/8._wp*e*sin(x)**4 + (3._wp/2._wp*d - 30._wp/8._wp*e)*sin(x)**2 + 3._wp/8._wp*e -d/2._wp + c - tsnow

end function fun_snow

