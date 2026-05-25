function cumsum(x,ist,iend) result(cumx)
use parameters, only: wp,n
implicit none
real(wp), dimension(n), intent(in) :: x
integer, intent(in) :: ist,iend
real(wp), dimension(n) :: cumx
integer :: i,k

cumx(:) =0._wp
if (ist < iend) then
        cumx = (/ (sum(x(ist:i)), i=ist,iend) /)
elseif (ist > iend) then
        k=0
        do i=ist,iend,-1    
        k=k+1
        cumx(k) = sum(x(ist:i:-1))
        end do
end if

end function cumsum