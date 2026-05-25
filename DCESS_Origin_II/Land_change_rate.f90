subroutine land_change_rate(np_t ,np_g,np_b,&
                            agr_t,agr_g,agr_b,&
                            awo_t,awo_g,awo_b,&
                            ali_t,ali_g,ali_b,&
                            asl_t,asl_g,asl_b,&
                            mp_t ,mp_g,mp_b,lb,pco2,al)
! Calculates land reservoirs rate change

use parameters, only:wp,sy,dc14,nlb
implicit none
real(wp), intent(in) :: np_t ,np_g,np_b,agr_t,agr_g,agr_b,awo_t,awo_g,awo_b,&
                        ali_t,ali_g,ali_b,asl_t,asl_g,asl_b,mp_t,mp_g,mp_b,lb(nlb),pco2(3)
real(wp), intent(out) :: al(nlb)
real(wp) :: p12co2,p13co2,p14co2,f1,f2,f3,f4,f5,f13al,f13lm,f14al,f14lm

! redefine
p12co2 = pco2(1) ! (atm)
p13co2 = pco2(2)
p14co2 = pco2(3)

f1 = 35._wp/60._wp
f2 = 25._wp/60._wp
f3 = 20._wp/25._wp
f4 = 10._wp/55._wp
f5 = 5._wp/25._wp

!-----------------------------------------------------
! 12C 
!-----------------------------------------------------
! Tropical forest
!-----------------------------------------------------
al(1) = (f1*np_t - agr_t*lb(1))/sy                                      ! 12C rate of change for leafy biomass, GtC/s
al(2) = (f2*np_t - awo_t*lb(2))/sy                                      ! 12C rate of change for wood, Gt/s
al(3) = (agr_t*lb(1) + f3*awo_t*lb(2) - ali_t*lb(3))/sy                 ! 12C rate of change for litter, GtC/s
al(4) = (f4*ali_t*lb(3) + f5*awo_t*lb(2) - asl_t*lb(4)-mp_t)/sy         ! 12C rate of change for soil, GtC/s

!-----------------------------------------------------
! Grasslands, savanna, deserts
!-----------------------------------------------------
al(5) = (f1*np_g - agr_g*lb(5))/sy                                      ! 12C rate of change for leafy biomass, GtC/s
al(6) = (f2*np_g - awo_g*lb(6))/sy                                      ! 12C rate of change for wood, Gt/s
al(7) = (agr_g*lb(5) + f3*awo_g*lb(6) - ali_g*lb(7))/sy                 ! 12C rate of change for litter, GtC/s
al(8) = (f4*ali_g*lb(7) + f5*awo_g*lb(6) - asl_g*lb(8)-mp_g)/sy         ! 12C rate of change for soil, GtC/s

!-----------------------------------------------------
! Extratropical forest
!-----------------------------------------------------
al(9)  = ( f1*np_b - agr_b*lb(9)  )/sy                                  ! 12C rate of change for leafy biomass, GtC/s
al(10) = ( f2*np_b - awo_b*lb(10) )/sy                                  ! 12C rate of change for wood, Gt/s
al(11) = ( agr_b*lb(9) + f3*awo_b*lb(10) - ali_b*lb(11) )/sy            ! 12C rate of change for litter, GtC/s
al(12) = ( f4*ali_b*lb(11) + f5*awo_b*lb(10) - asl_b*lb(12)-mp_b )/sy   ! 12C rate of change for soil, GtC/s


!-----------------------------------------------------
! 13C 
f13al = 0.9819_wp       ! Carbon ion fractionation factors for land air exchange
f13lm = 0.9700_wp       ! Carbon ion fractionation factor during land biosphere methane production
!-----------------------------------------------------
! Tropical forest
!-----------------------------------------------------
al(13) = (f1*np_t*p13co2/p12co2*f13al - agr_t*lb(13))/sy    ! 13C rate of change for leafy biomass, GtC/s   
al(14) = (f2*np_t*p13co2/p12co2*f13al - awo_t*lb(14))/sy    ! 13C rate of change for wood, GtC/s
al(15) = (agr_t*lb(13) + f3*awo_t*lb(14) - ali_t*lb(15))/sy ! 13C rate of change for litter, GtC/s
al(16) = (f4*ali_t*lb(15) + f5*awo_t*lb(14) - &
             asl_t*lb(16) - mp_t*lb(16)/lb(4)*f13lm)/sy           ! 13C rate of change for soil biomass, GtC/s 

!-----------------------------------------------------
! Grasslands, savanna, deserts
!-----------------------------------------------------
al(17) = (f1*np_g*p13co2/p12co2*f13al - agr_g*lb(17))/sy    ! 13C rate of change for leafy biomass, GtC/s   
al(18) = (f2*np_g*p13co2/p12co2*f13al - awo_g*lb(18))/sy    ! 13C rate of change for wood, GtC/s
al(19) = (agr_g*lb(17) + f3*awo_g*lb(18) - ali_g*lb(19))/sy ! 13C rate of change for litter, GtC/s
al(20) = (f4*ali_g*lb(19) + f5*awo_g*lb(18) - &
             asl_g*lb(20)-mp_g*lb(20)/lb(8)*f13lm)/sy          ! 13C rate of change for soil biomass, GtC/s 

!-----------------------------------------------------
! Extratropical forest
!-----------------------------------------------------
al(21) = (f1*np_b*p13co2/p12co2*f13al - agr_b*lb(21))/sy    ! 13C rate of change for leafy biomass, GtC/s   
al(22) = (f2*np_b*p13co2/p12co2*f13al - awo_b*lb(22))/sy    ! 13C rate of change for wood, GtC/s
al(23) = (agr_b*lb(21) + f3*awo_b*lb(22) - ali_b*lb(23))/sy ! 13C rate of change for litter, GtC/s
al(24) = (f4*ali_b*lb(23) + f5*awo_b*lb(22) - &
             asl_b*lb(24)-mp_b*lb(24)/lb(12)*f13lm)/sy          ! 13C rate of change for soil biomass, GtC/s 

!-----------------------------------------------------
! 14C
f14al = 1._wp-2._wp*(1._wp-f13al)   ! Carbon ion fractionation factors for land air exchange
f14lm = 1._wp-2._wp*(1._wp-f13lm)   ! Carbon ion fractionation factor during land biosphere methane production
!-----------------------------------------------------
! Tropical forest
!-----------------------------------------------------
al(25) = (f1*np_t*p14co2/p12co2*f14al - (agr_t+dc14*sy)*lb(25))/sy    ! 14C rate of change for leafy biomass, GtC/s   
al(26) = (f2*np_t*p14co2/p12co2*f14al - (awo_t+dc14*sy)*lb(26))/sy    ! 14C rate of change for wood, GtC/s
al(27) = (agr_t*lb(25) + f3*awo_t*lb(26) - (ali_t+dc14*sy)*lb(27))/sy ! 14C rate of change for litter, GtC/s
al(28) = (f4*ali_t*lb(27) + f5*awo_t*lb(26) - &
           (asl_t+dc14*sy)*lb(28)-mp_t*lb(28)/lb(4)*f14lm)/sy           ! 14C rate of change for soil biomass, GtC/s 

!-----------------------------------------------------
! Grasslands, savanna, deserts
!-----------------------------------------------------
al(29) = (f1*np_g*p14co2/p12co2*f14al - (agr_g+dc14*sy)*lb(29))/sy    ! 14C rate of change for leafy biomass, GtC/s   
al(30) = (f2*np_g*p14co2/p12co2*f14al - (awo_g+dc14*sy)*lb(30))/sy    ! 14C rate of change for wood, GtC/s
al(31) = (agr_g*lb(29) + f3*awo_g*lb(30) - (ali_g+dc14*sy)*lb(31))/sy ! 14C rate of change for litter, GtC/s
al(32) = (f4*ali_g*lb(31) + f5*awo_g*lb(30) - &
            (asl_g+dc14*sy)*lb(32)-mp_g*lb(32)/lb(8)*f14lm)/sy           ! 14C rate of change for soil biomass, GtC/s 

!-----------------------------------------------------
! Extratropical forest
!-----------------------------------------------------
al(33) = (f1*np_b*p14co2/p12co2*f14al - (agr_b+dc14*sy)*lb(33))/sy    ! 14C rate of change for leafy biomass, GtC/s   
al(34) = (f2*np_b*p14co2/p12co2*f14al - (awo_b+dc14*sy)*lb(34))/sy    ! 14C rate of change for wood, GtC/s
al(35) = (agr_b*lb(33) + f3*awo_b*lb(34) - (ali_b+dc14*sy)*lb(35))/sy ! 14C rate of change for litter, GtC/s
al(36) = (f4*ali_b*lb(35) + f5*awo_b*lb(34) - &
            (asl_b+dc14*sy)*lb(36)-mp_b*lb(36)/lb(12)*f14lm)/sy          ! 14C rate of change for soil biomass, GtC/s 
end subroutine land_change_rate