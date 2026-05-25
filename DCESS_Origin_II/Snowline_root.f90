subroutine snowline_root(ptanh,fsn_n,fsn_s)
! Calculates position of snowline based on tsnow

use parameters, only: wp,pi,bsh
implicit none
real(wp), intent(in) :: ptanh(3)
real(wp), intent(out) :: fsn_n,fsn_s
real(wp) :: zbrent_snow
real(wp), external :: fun_snow
real(wp) :: a,b,tol
real(wp) :: phi(101),taan(101)
integer :: i,ipos

phi  = (/ ( (real(i,wp)-1._wp) * (pi/2._wp) / 100._wp, i = 1, 101 ) /)
taan = ptanh(1) + ptanh(2)*(1._wp/2._wp)*(3._wp*sin(phi)**2-1._wp) + &
       ptanh(3)*(1._wp/8._wp)*(35._wp*sin(phi)**4 - 30._wp*sin(phi)**2 + 3._wp)

ipos = count(taan > 0._wp)

a   = 0._wp
b   = pi/2._wp
tol = 1.0e-6_wp

if (ipos==size(phi)) then
    fsn_n = pi/2._wp - 1.e-3
else
    fsn_n = zbrent_snow(fun_snow, a, b, tol,ptanh(1),ptanh(2),ptanh(3))
    fsn_n = min(fsn_n,pi/2._wp)
end if

fsn_s = bsh ! fixed position

end subroutine snowline_root

include 'zbrent_snow.f90'
include 'fun_snow.f90'