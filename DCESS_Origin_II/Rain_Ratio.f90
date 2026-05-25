real(wp) function rain_ratio(ome,sst) result(rp)
! Calculates calcite rain ratio
! Inputs:
! ome: calcite saturation state
! sst: surface ocean temperature
! Outputs: calcite rain ratio

use parameters, only:wp,rpm
implicit none
real(wp), intent(in) :: ome,sst
real(wp) :: tr,p1,p2

tr  = 10._wp                                           ! [oC] Reference temperatur, st. val 10
p1  = 1._wp                                            ! [-] O.Marchal 98, Maier-Reimer 93
p2  = 0.1565_wp                                        ! [oC^-1] Maier-Reimer 93

if (ome < 1._wp) then
    rp = 0._wp
else
    rp = rpm* (p1*exp(p2*(sst-tr))/(1._wp+p1*exp(p2*(sst-tr))))*((ome-1._wp)/(1._wp+(ome-1._wp)))
end if

end function rain_ratio