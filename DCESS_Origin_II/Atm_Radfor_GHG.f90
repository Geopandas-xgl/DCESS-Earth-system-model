subroutine atm_radfor_ghg(pco2,pch4,pn2o,rco2,rch4,rn2o,aol)
! Calculates radiatve forcing for co2,ch4 and n2o based on Byrne and Goldblatt (2014)
! Inputs:
! Partial pressure of co2, ch4 and n2o (pco2,pch4,pn2o) (atm)
! Outputs:
! Radiative forcing for co2, ch4 and n2o (rco2,rch4,rn2o,aol) (W/m2)

use parameters, only: wp,pco2int,pch4int,pn2oint,nab

implicit none
real(wp), dimension(nab), intent(in) :: pco2,pch4,pn2o
real(wp), dimension(nab), intent(out) :: rco2,rch4,rn2o,aol
real(wp) :: pch4int2,a,b,c,d,m1
integer :: i

rco2 = 5.32_wp*log(pco2/pco2int) + 0.39_wp*log(pco2/pco2int)**2

! Constants for new CH4 forcing fit at high CH4
a = -4.054_wp
b = -6482._wp
c = 9.542_wp
d = -60.99_wp

pch4int2=pch4int-0.005e-6_wp
m1=2.5e-6_wp

do i=1,nab
        if (pch4(i) < 2.5e-6_wp) then
        rch4(i) = 1173._wp*(pch4(i)**(.5)-pch4int2**(.5))-71636._wp*(pch4(i)**(.5)-pch4int2**(.5))**2
        elseif (pch4(i) > 100e-6_wp) then
        rch4(i) = a*exp(b*pch4(i)) + c*exp(d*pch4(i))-0.867_wp
        else
        rch4(i) = 0.824_wp+0.8_wp*log(pch4(i)/m1)+0.2_wp*log(pch4(i)/m1)**2
        end if
end do
rn2o = 3899._wp*(pn2o**(.5) - pn2oint**(.5))+38256._wp*(pn2o**(.5) - pn2oint**(.5))**2

! Overlap
! aol  = -16.16_wp*exp(-0.036_wp*(log(pco2-pco2int)-0.0024_wp)**2&
!        -0.05_wp *(log(pn2o-pn2oint)+6.5_wp)**2)-24._wp*exp(-0.02_wp*(log(pch4-pch4int2)-0.01_wp)**2&
!        -0.044_wp*(log(pn2o-pn2oint)+7.73_wp)**2)
aol = 0._wp

end subroutine atm_radfor_ghg