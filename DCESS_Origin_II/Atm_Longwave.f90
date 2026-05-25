subroutine atm_longwave(pch4,pn2o,pco2,t_atm,tot_lw)
! Calculates longwave radiation for each atmospheric box
! Inputs:
! pch4, pn2o, pco2: partial pressures of ch4, n2o and co2 (mu atm)
! t_atm: atmospheric temperature
! Outputs:
! tot_lw: total longwave radiation for each atmospheric box (W)

use parameters, only : wp,aa_h,aa_m,aa_e,nab

implicit none
real(wp), dimension(nab),intent(in) :: pch4,pn2o,pco2,t_atm
real(wp), dimension(nab), intent(out) :: tot_lw
real(wp), dimension(nab) :: rco2,rch4,rn2o,aol,a_tot
real(wp), parameter :: ao = 206.56_wp,&         ! Zero degree outgoing radiation  [W/m2]             
                       b  = 2.20_wp             ! LW radiation sensitivity  [W/m2/K]

! Radiative forcing for co2, ch4 and n2o
call atm_radfor_ghg(pco2,pch4,pn2o,rco2,rch4,rn2o,aol)

a_tot = ao - ( rco2 + rch4 + rn2o + aol )

tot_lw(1) =  aa_h*( (a_tot(1)  + b*t_atm(1) ) )
tot_lw(2) =  aa_m*( (a_tot(2)  + b*t_atm(2) ) )
tot_lw(3) =  aa_e*( (a_tot(3)  + b*t_atm(3) ) )
tot_lw(4) =  aa_e*( (a_tot(4)  + b*t_atm(4) ) )
tot_lw(5) =  aa_m*( (a_tot(5)  + b*t_atm(5) ) )
tot_lw(6) =  aa_h*( (a_tot(6)  + b*t_atm(6) ) )


end subroutine atm_longwave

include 'Atm_Radfor_GHG.f90'