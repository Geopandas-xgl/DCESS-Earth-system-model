subroutine med_flux(rho_n,rho_s,dy,width,phi,vy)
! Calculates ocean meridional velocity crossing an ocean box boundary
! Inputs:
! rho_n,rho_s: north and south ocean density at the box boundary
! dy: distance between two adjacent boxes
! witdh: ocean box width
! phi: boundary latitude
! Outputs:
! vy: meridional flux (m3/s)

use parameters, only: wp,rho0,rq,d,e_r,pi,n,g
! Meridional flux
implicit none
real(wp), dimension(n), intent(in) :: rho_n,rho_s
real(wp), intent(in) :: dy,width,phi
real(wp), dimension(n), intent(out) :: vy
real(wp), dimension(n) :: cum,py

cum     = cumsum((rho_n - rho_s)*d,1,n)
py      = g*( cum - (rho_n - rho_s)*.5_wp*d)/( dy*e_r )
vy      = -1._wp/(rho0*rq)*py
vy(1)   = 0._wp
vy(2:n) = vy(2:n) - sum(vy(2:n))/size(vy(2:n))
vy = vy*d*2._wp*pi*width*e_r*cos(phi)

contains
function cumsum(x,ist,iend) result(cumx)
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
end subroutine med_flux