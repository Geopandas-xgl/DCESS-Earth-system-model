subroutine atm_met(at,mdrn,mdrmn,mdren,mdres,mdrms,mdrs)
use parameters, only: wp,nta,rvahp,rvahm,rvahe,pch4int,mdts,sy,nab
! Calculates methane oxidation in the atmosphere
! Inputs:
! at: atmospheric tracers
! Outputs: ch4 oxidation for each atmospheric box  and for each isotope (12-13)
! mdrn,mdrmn,mdren,mdres,mdrms,mdrs


implicit none
real(wp), intent(in) :: at(nta,nab)
real(wp), dimension(2), intent(out) :: mdrn,mdrmn,mdren,mdres,mdrms,mdrs
integer :: i12,i13

i12 = 3         ! 12-CH4
i13 = 8         ! 13-CH4
mdrn  = methane(at(i12,1),at(i13,1),rvahp)
mdrmn = methane(at(i12,2),at(i13,2),rvahm)
mdren = methane(at(i12,3),at(i13,3),rvahe)
mdres = methane(at(i12,4),at(i13,4),rvahe)
mdrms = methane(at(i12,5),at(i13,5),rvahm)
mdrs  = methane(at(i12,6),at(i13,6),rvahp)

contains
function methane(p12ch4,p13ch4,rva) result(mdr)
implicit none
real(wp), intent(in) :: p12ch4,p13ch4,rva
real(wp), dimension(2) :: mdr
real(wp) :: m,a,b,lamch4

m = (p12ch4-pch4int)/pch4int
a = 0.78_wp
b = 11._wp
lamch4 = rva/(mdts*sy)*(1._wp-a*m/(m+b))
mdr(1) =  lamch4*p12ch4                         ! atmospheric methane oxidation to CO2 - C12
mdr(2) =  lamch4*p12ch4*p13ch4/p12ch4           ! atmospheric methane oxidation to CO2 - C13

end function methane
end subroutine atm_met