subroutine land_decay_rates(np_t0,np_g0,np_b0,&
                            gr_t0,wo_t0,li_t0,sl_t0,mp_t0,&
                            gr_g0,wo_g0,li_g0,sl_g0,mp_g0,&
                            gr_b0,wo_b0,li_b0,sl_b0,mp_b0,&
                            ta_t,ta_t0,ta_g,ta_g0,ta_b,ta_b0,&
                            agr_t,awo_t,ali_t,asl_t,&
                            agr_g,awo_g,ali_g,asl_g,&
                            agr_b,awo_b,ali_b,asl_b)
! Calculates decay rates for every carbon reservoir and vegetatin zone
use parameters, only: wp,q10
implicit none
real(wp), intent(in) :: np_t0,np_g0,np_b0,&
                        gr_t0,wo_t0,li_t0,sl_t0,mp_t0,&
                        gr_g0,wo_g0,li_g0,sl_g0,mp_g0,&
                        gr_b0,wo_b0,li_b0,sl_b0,mp_b0,&
                        ta_t,ta_t0,ta_g,ta_g0,ta_b,ta_b0
real(wp), intent(out) :: agr_t,awo_t,ali_t,asl_t,&
                         agr_g,awo_g,ali_g,asl_g,&
                         agr_b,awo_b,ali_b,asl_b

!-----------------------------------------------------
! Tropical forest
!-----------------------------------------------------
agr_t = (35._wp/60._wp )*np_t0/gr_t0                                        ! decay rate for leafy biomass, GtC/yr
awo_t = (25._wp/60._wp )*np_t0/wo_t0                                        ! decay rate for woody biomass, GtC/yr 
ali_t = ((55._wp/60._wp)*np_t0/li_t0)      *q10**((ta_t-ta_t0)/10._wp)      ! decay rate for litter biomass with Q10 T dependence, GtC/yr  
asl_t = ((15._wp/60._wp)*np_t0-mp_t0)/sl_t0*q10**((ta_t-ta_t0)/10._wp)      ! decay rate for soil biomass with Q10 T dependence, GtC/yr

!-----------------------------------------------------
! Grassland, desert, savanna
!-----------------------------------------------------
agr_g = (35._wp/60._wp )*np_g0/gr_g0                                        ! decay rate for leafy biomass, GtC/yr
awo_g = (25._wp/60._wp )*np_g0/wo_g0                                        ! decay rate for woody biomass, GtC/yr 
ali_g = ((55._wp/60._wp)*np_g0/li_g0)      *q10**((ta_g-ta_g0)/10._wp)      ! decay rate for litter biomass with Q10 T dependence, GtC/yr  
asl_g = ((15._wp/60._wp)*np_g0-mp_g0)/sl_g0*q10**((ta_g-ta_g0)/10._wp)      ! decay rate for soil biomass with Q10 T dependence, GtC/yr  

!-----------------------------------------------------
! Extratropical forest
!-----------------------------------------------------
agr_b = (35._wp/60._wp )*np_b0/gr_b0                                        ! decay rate for leafy biomass, GtC/yr
awo_b = (25._wp/60._wp )*np_b0/wo_b0                                        ! decay rate for woody biomass, GtC/yr 
ali_b = ((55._wp/60._wp)*np_b0/li_b0)      *q10**((ta_b-ta_b0)/10._wp)      ! decay rate for litter biomass with Q10 T dependence, GtC/yr  
asl_b = ((15._wp/60._wp)*np_b0-mp_b0)/sl_b0*q10**((ta_b-ta_b0)/10._wp)      ! decay rate for soil biomass with Q10 T dependence, GtC/yr  


end subroutine land_decay_rates