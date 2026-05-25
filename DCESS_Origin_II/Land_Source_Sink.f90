subroutine land_source_sink(mp_b,mp_g,mp_t,np_b,np_g,np_t,npp_b,npp_g,npp_t,&
                            ali_b,ali_g,ali_t,asl_b,asl_g,asl_t,lb,pco2,&
                            alss)
! Calculates source/sinks of land biosphere for every atmopsheric gas

use parameters, only: wp,mgt,sy,nlb
implicit none
real(wp), intent(in) :: mp_b,mp_g,mp_t,np_b,np_g,np_t,npp_b,npp_g,npp_t,&
                        ali_b,ali_g,ali_t,asl_b,asl_g,asl_t,lb(nlb),pco2(3)
real(wp), intent(out) :: alss(8,3)
real(wp) :: f13al,f13lm,f14al,f14lm,p12co2,p13co2,p14co2,f

! redefine
p12co2 = pco2(1)
p13co2 = pco2(2)
p14co2 = pco2(3)

f13al = 0.9819_wp                   ! Carbon ion fractionation factors for land air exchange
f13lm = 0.9700_wp                   ! Carbon ion fractionation factor during land biosphere methane production
f14al = 1._wp-2._wp*(1._wp-f13al)   ! Carbon ion fractionation factors for land air exchange
f14lm = 1._wp-2._wp*(1._wp-f13lm)   ! Carbon ion fractionation factor during land biosphere methane production

f = 45._wp/55._wp

! CH4-12
alss(3,:) = (/mp_b,mp_g,mp_t/)*mgt/sy
! N2O
alss(4,:) = (/np_b,np_g,np_t/)
! CO2-12
alss(5,:) = (/-npp_b + f*ali_b*lb(11) + asl_b*lb(12),&
              -npp_g + f*ali_g*lb(7 ) + asl_g*lb(8 ),&
              -npp_t + f*ali_t*lb(3 ) + asl_t*lb(4 ) /)*mgt/sy
! CO2-13
alss(6,:) = (/-npp_b*p13co2/p12co2*f13al + f*ali_b*lb(23) + asl_b*lb(24),&
              -npp_g*p13co2/p12co2*f13al + f*ali_g*lb(19) + asl_g*lb(20),&
              -npp_t*p13co2/p12co2*f13al + f*ali_t*lb(15) + asl_t*lb(16) /)*mgt/sy
! CO2-14
alss(7,:) = (/-npp_b*p14co2/p12co2*f14al + f*ali_b*lb(35) + asl_b*lb(36) + mp_b*lb(36)/lb(12)*f14lm,&
              -npp_g*p14co2/p12co2*f14al + f*ali_g*lb(31) + asl_g*lb(32) + mp_g*lb(32)/lb(8 )*f14lm,&
              -npp_t*p14co2/p12co2*f14al + f*ali_t*lb(27) + asl_t*lb(28) + mp_t*lb(28)/lb(4 )*f14lm/)*mgt/sy
! CH4-13
alss(8,:) = (/ mp_b*lb(24)/lb(12)*f13lm,&
               mp_g*lb(20)/lb(8 )*f13lm,&
               mp_t*lb(16)/lb(4 )*f13lm /)*mgt/sy

end subroutine land_source_sink