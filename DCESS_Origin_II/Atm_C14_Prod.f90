subroutine atm_c14_prod(at,p14c)
use parameters,only: wp,nta,r14oas,r13pdb,nab
! Calculate for all atmopsheric boxes 14C production
! Inputs:
! at: atmospheric tracers
! Outputs:
! p14c: atmopsheric 14c production [hn mn en es ms hs]

implicit none
real(wp), dimension(nta,nab), intent(in) :: at
real(wp), dimension(nab), intent(out) :: p14c
integer :: i

do i=1,nab
    p14c(i) = c14_atm(at(5,i),at(6,i),at(7,i))
end do

contains

real(wp) function c14_atm(c12,c13,c14) result(p14c)
implicit none
real(wp), intent(in) :: c12,c13,c14
real(wp) :: d14a

! Calculate delta^13C and Delta^14C values in permil for the atmosphere 
! and set ^14C production to maintain Delta^14C approx = 0

d14a = (c14/c12/r14oas*(r13pdb*.975_wp/(c13/c12))**2-1._wp)*1.e3_wp
p14c = max(0._wp,1.e5_wp*(0._wp-d14a))                  ! Production [atoms/(sm^2)]

end function c14_atm

end subroutine atm_c14_prod