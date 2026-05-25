subroutine ekman(rho_f,rho_t,evf,evfx,edl)
! Calculates vertical level (equivalent density level) where a flux (ekman) es injected
! between two adjacent boxes
! Inputs:
! rho_f: rho from
! rho_t: rho to
! evf: volume flux
! Outputs: 
! edl  : equivalent density level of equatorward ekman flow
! evfx : injection vector
use parameters, only: wp,n

implicit none
real(wp), intent(in) :: rho_f(n),rho_t(n),evf
real(wp), intent(out) :: evfx(n)
integer, intent(out) :: edl
logical, dimension(n) :: mask
integer :: i,cont

mask = ( (rho_f(1) - rho_t) < 0._wp)
cont = count(mask)
if (cont == 0) then ! Set edl to bottom if rho_f(1)>rho_to
        edl=n
else
        edl = minval(pack([(i, i = 1, n)], mask))
end if

evfx(1:n) = 0._wp
if (edl==n) then 
        evfx(edl-1) = 0.25_wp*evf
        evfx(edl  ) = 0.75_wp*evf
elseif (edl == 1) then
        evfx(edl  ) = 0.75_wp*evf
        evfx(edl+1) = 0.25_wp*evf
else
        evfx(edl-1) = 0.25_wp*evf
        evfx(edl  ) = 0.50_wp*evf
        evfx(edl+1) = 0.25_wp*evf
end if


end subroutine ekman