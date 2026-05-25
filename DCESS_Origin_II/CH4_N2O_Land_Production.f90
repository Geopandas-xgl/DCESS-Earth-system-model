
subroutine ch4_n2o_land_production(mp_t0,mp_b0,np_t0,np_g0,np_b0,ms_t,ms_g,ms_b,&
                                   ms_t0,ms_g0,ms_b0,ta_t,ta_g,ta_b,ta_t0,ta_g0,ta_b0,&
                                   mp_t,mp_g,mp_b,np_t,np_g,np_b)
! Calculates ch4 and n2o land production

use parameters, only: wp,q10_met
implicit none
real(wp), intent(in) :: mp_t0,mp_b0,np_t0,np_g0,np_b0,ms_t,ms_g,ms_b,ms_t0,ms_g0,ms_b0,ta_t,ta_g,ta_b,ta_t0,ta_g0,ta_b0
real(wp), intent(out) :: mp_t,mp_g,mp_b,np_t,np_g,np_b

! Assumption of methane emissions only from wet areas: T and B
! Assumption that T and B areas both have about the same PI CH4 emissions
mp_t = mp_t0*ms_t/(ms_t0)*q10_met**((ta_t - ta_t0)/10._wp)    ! tropical forest methane production with Q10 T dependence
mp_g = 0._wp
mp_b = mp_b0*ms_b/(ms_b0)*q10_met**((ta_b - ta_b0)/10._wp)    ! Extratropical forest methane production with Q10 T dependence

! N2O production for the 3 different vegetation types
np_t = np_t0*ms_t/(ms_t0)*q10_met**((ta_t - ta_t0)/ 10._wp)     ! N2O production tropical forest
np_g = np_g0*ms_g/(ms_g0)*q10_met**((ta_g - ta_g0)/ 10._wp)     ! N2O production grassland sd
np_b = np_b0*ms_b/(ms_b0)*q10_met**((ta_b - ta_b0)/ 10._wp)     ! N2O production boreal forest

end subroutine ch4_n2o_land_production