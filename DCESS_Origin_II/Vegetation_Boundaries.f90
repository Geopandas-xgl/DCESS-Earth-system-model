subroutine vegetation_boundaries(ta,pc,ftrop,fgrass)
! Calculates vegetation boundaries
! Inputs:
! ta: Annual mean atm temperature
! pc: polynomial coeffs
! Outputs:
! ftrop: Latitude line of Tropical Forest
! fgrass: Latitude line of grassland, desert, savanna

use parameters, only:wp,dr

implicit none
real(wp), intent(in) :: ta,pc(12)
real(wp), intent(out) :: ftrop,fgrass
real(wp) :: del_ta,y1,y2

del_ta  = ta - 15._wp ! deviation of global temperature from PI

y1 = pc(1)*del_ta**5 + pc(2)*del_ta**4 + pc(3)*del_ta**3 + pc(4) *del_ta**2 + pc(5) *del_ta + pc(6)
y2 = pc(7)*del_ta**5 + pc(8)*del_ta**4 + pc(9)*del_ta**3 + pc(10)*del_ta**2 + pc(11)*del_ta + pc(12)

ftrop  = abs(y1)*dr         ! Latitude line of Tropical Forest
fgrass = abs(y2)*dr         ! Latitude line of grassland, desert, savanna


end subroutine vegetation_boundaries