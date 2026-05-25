subroutine land_bios(t,taam,ptanh,ptash,at,lbn,lbs,olf_bios,&
                     nppn,npps,aln,als,alhn,almn,alen,ales,alms,alhs)
! Calculates Land Biomass changes and land-sea gas exchanges of CO2, CH4 and N2O
! Inputs: 
! t: time
! taam: annual mean atmospheric temperature
! ptanh,ptash: factors for atmospheric temperature meridional profile
! at: atmopsheric tracers
! lbn, lbs: Land biosphere tracers for NH and SH
! olf_bios: ocean land distribution as function of latitude

! Outputs:
! aln, als: rate changes. Structure is:
! ----------------------------------------------------------------------------------
!          al(1)  - increment change in leafy biomass   12C in tropical forests (GtC/s)
!          al(2)  - increment change in woody biomass   12C in tropical forests (GtC/s)
!          al(3)  - increment change in litter biomass  12C in tropical forests (GtC/s)
!          al(4)  - increment change in soil biomass    12C in tropical forests (GtC/s)
! ----------------------------------
!          al(5)  - increment change in leafy biomass   12C in grassland, desert, savanna (GtC/s)
!          al(6)  - increment change in woody biomass   12C in grassland, desert, savanna (GtC/s)
!          al(7)  - increment change in litter biomass  12C in grassland, desert, savanna (GtC/s)
!          al(8)  - increment change in soil biomass    12C in grassland, desert, savanna (GtC/s)
! ----------------------------------
!          al(9)  - increment change in leafy biomass  12C in extratropical forests (GtC/s)
!          al(10) - increment change in woody biomass  12C in extratropical forests (GtC/s)
!          al(11) - increment change in litter biomass 12C in extratropical forests (GtC/s)
!          al(12) - increment change in soil biomass   12C in extratropical forests (GtC/s)
! ----------------------------------------------------------------------------------
!          al(13) - increment change in leafy biomass  13C in tropical forests (GtC/s)
!          al(14) - increment change in woody biomass  13C in tropical forests (GtC/s)
!          al(15) - increment change in litter biomass 13C in tropical forests (GtC/s)
!          al(16) - increment change in soil biomass   13C in tropical forests (GtC/s)
! ----------------------------------
!          al(17) - increment change in leafy biomass  13C in grassland, desert, savanna (GtC/s)
!          al(18) - increment change in woody biomass  13C in grassland, desert, savanna (GtC/s)
!          al(19) - increment change in litter biomass 13C in grassland, desert, savanna (GtC/s)
!          al(20) - increment change in soil biomass   13C in grassland, desert, savanna (GtC/s)
! ----------------------------------
!          al(21) - increment change in leafy biomass  13C in extratropical forests (GtC/s)
!          al(22) - increment change in woody biomass  13C in extratropical forests (GtC/s)
!          al(23) - increment change in litter biomass 13C in extratropical forests (GtC/s)
!          al(24) - increment change in soil biomass   13C in extratropical forests (GtC/s)
! ----------------------------------------------------------------------------------
!          al(25) - increment change in leafy biomass  14C in tropical forests (GtC/s)
!          al(26) - increment change in woody biomass  14C in tropical forests (GtC/s)
!          al(27) - increment change in litter biomass 14C in tropical forests (GtC/s)
!          al(28) - increment change in soil biomass   14C in tropical forests (GtC/s)
! ----------------------------------
!          al(29) - increment change in leafy biomass  14C in grassland, desert, savanna (GtC/s)
!          al(30) - increment change in woody biomass  14C in grassland, desert, savanna (GtC/s)
!          al(31) - increment change in litter biomass 14C in grassland, desert, savanna (GtC/s)
!          al(32) - increment change in soil biomass   14C in grassland, desert, savanna (GtC/s)
! ----------------------------------
!          al(33) - increment change in leafy biomass  14C in extratropical forests (GtC/s)
!          al(34) - increment change in woody biomass  14C in extratropical forests (GtC/s)
!          al(35) - increment change in litter biomass 14C in extratropical forests (GtC/s)
!          al(36) - increment change in soil biomass   14C in extratropical forests (GtC/s)
! ----------------------------------------------------------------------------------
! alhn,almn,alen,ales,alms,alhs: Atmospheric fluxes of CO2, CH4 and N2O
! nppn,npps: net primary production for NH and SH

use parameters, only: wp,nta,nlb,np_tn0,np_gn0,np_bn0,np_ts0,np_gs0,np_bs0,&
     lbmp_tn0,lbmp_bn0,lbnp_tn0,lbnp_gn0,lbnp_bn0,sl_tn0,sl_gn0,sl_bn0,&
     lbmp_ts0,lbmp_bs0,lbnp_ts0,lbnp_gs0,lbnp_bs0,sl_ts0,sl_gs0,sl_bs0,&
     ta_tn0,ta_gn0,ta_bn0,ta_ts0,ta_gs0,ta_bs0,&
     gr_tn0,wo_tn0,li_tn0,gr_gn0,wo_gn0,li_gn0,lbmp_gn0,gr_bn0,wo_bn0,li_bn0,&
     gr_ts0,wo_ts0,li_ts0,gr_gs0,wo_gs0,li_gs0,lbmp_gs0,gr_bs0,wo_bs0,li_bs0,&
     fdh,fdl,fde,sy,nab
implicit none

real(wp),intent(in) :: t,taam(nab),ptanh(3),ptash(3),at(nta,nab),lbn(nlb),lbs(nlb),olf_bios(1801,2)
real(wp), intent(out) :: nppn(3),npps(3),aln(nlb),als(nlb),alen(8),almn(8),alhn(8),ales(8),alms(8),alhs(8)
real(wp) :: tahn,tamn,taen,taes,tams,tahs,ta_mean_n,ta_mean_s,&
p12co2n,p13co2n,p14co2n,p12co2s,p13co2s,p14co2s,p_nh(12),p_sh(12),&
ftropn,fgrassn,ftrops,fgrasss,pta_amnh(3),pta_amsh(3),fsn_n,fsn_s,&
ta_tn,ta_gn,ta_bn,ta_ts,ta_gs,ta_bs,factfn,facgsn,facbfn,factfs,facgss,facbfs,&
tevaln,tevals,nptfn,npgsn,npbfn,nptfs,npgss,npbfs,&
mp_tn,mp_gn,mp_bn,mp_ts,mp_gs,mp_bs,np_tn,np_gn,np_bn,np_ts,np_gs,np_bs,&
agr_tn,awo_tn,ali_tn,asl_tn,agr_gn,awo_gn,ali_gn,asl_gn,agr_bn,awo_bn,ali_bn,asl_bn,&
agr_ts,awo_ts,ali_ts,asl_ts,agr_gs,awo_gs,ali_gs,asl_gs,agr_bs,awo_bs,ali_bs,asl_bs,pco2(3),&
alssn(8,3),alsss(8,3),alenx,almnx,alhnx,alesx,almsx,alhsx
integer :: i
real(wp), external :: ta_mean

tahn = taam(1)
tamn = taam(2)
taen = taam(3)
taes = taam(4)
tams = taam(5)
tahs = taam(6)

ta_mean_n = tahn*(1._wp-sin(fdh)) + tamn*(sin(fdh)-sin(fdl)) + taen*sin(fdl)
ta_mean_s = tahs*(1._wp-sin(fdh)) + tams*(sin(fdh)-sin(fdl)) + taes*sin(fdl)

p12co2n = at(5,1)*(1._wp-sin(fdh)) + at(5,2)*(sin(fdh)-sin(fdl)) + at(5,3)*sin(fdl)
p13co2n = at(6,1)*(1._wp-sin(fdh)) + at(6,2)*(sin(fdh)-sin(fdl)) + at(6,3)*sin(fdl)
p14co2n = at(7,1)*(1._wp-sin(fdh)) + at(7,2)*(sin(fdh)-sin(fdl)) + at(7,3)*sin(fdl)

p12co2s = at(5,6)*(1._wp-sin(fdh)) + at(5,5)*(sin(fdh)-sin(fdl)) + at(5,4)*sin(fdl)
p13co2s = at(6,6)*(1._wp-sin(fdh)) + at(6,5)*(sin(fdh)-sin(fdl)) + at(6,4)*sin(fdl)
p14co2s = at(7,6)*(1._wp-sin(fdh)) + at(7,5)*(sin(fdh)-sin(fdl)) + at(7,4)*sin(fdl)

!-----------------------------------------------------------------
! Calculate vegetation boundaries using annual mean temperature
!-----------------------------------------------------------------
p_nh = (/ -1.803e-5_wp, -0.0005809_wp, -0.005168_wp, 0.0497_wp , 1.092_wp, 11.28_wp,&
1.152e-5_wp, -0.0001785_wp, -0.004557_wp, 0.04156_wp, 1.017_wp, 37.77_wp /)

p_sh = (/ 8.413e-5_wp, 0.0007339_wp, -0.008333_wp, -0.08764_wp, -0.6965_wp, -19.66_wp,&
-6.51e-6_wp , 0.0002131_wp,  0.001857_wp, -0.041_wp  , -0.6615_wp, -35.63_wp /)

call vegetation_boundaries(ta_mean_n,p_nh,ftropn,fgrassn)
call vegetation_boundaries(ta_mean_s,p_sh,ftrops,fgrasss)

call get_poly_pta(tahn,tamn,taen,taes,tams,tahs,pta_amnh,pta_amsh)
call snowline_root(pta_amnh,fsn_n,fsn_s)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Mean temperature for each zone
!-----------------------------------------------------------------
ta_tn = ta_mean(ptanh,ftropn ,fde    )
ta_gn = ta_mean(ptanh,fgrassn,ftropn )
ta_bn = ta_mean(ptanh,fsn_n  ,fgrassn)

ta_ts = ta_mean(ptash,ftrops ,fde    )
ta_gs = ta_mean(ptash,fgrasss,ftrops )
ta_bs = ta_mean(ptash,fdh    ,fgrasss)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Area change in the different zones
!-----------------------------------------------------------------
call olf_lat_nh(ftropn ,fgrassn ,fsn_n,olf_bios,factfn,facgsn,facbfn)
call olf_lat_sh(-ftrops,-fgrasss,-fdh ,olf_bios,factfs,facgss,facbfs)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! New production for the 3 vegetation types
!-----------------------------------------------------------------
tevaln = mod(t     , sy)
if (tevaln .eq. 0._wp) then
tevaln = sy
end if
tevals = mod(t+sy/2, sy)
if (tevals .eq. 0._wp) then
tevals = sy
end if

call land_new_prod(tevaln,np_tn0,np_gn0,np_bn0,factfn,facgsn,facbfn,p12co2n,nptfn,npgsn,npbfn,nppn)
call land_new_prod(tevals,np_ts0,np_gs0,np_bs0,factfs,facgss,facbfs,p12co2s,nptfs,npgss,npbfs,npps)
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! CH4 and N2O production for the 3 different vegetation types
!-----------------------------------------------------------------
call ch4_n2o_land_production(lbmp_tn0,lbmp_bn0,lbnp_tn0,lbnp_gn0,lbnp_bn0,&
            lbn(4),lbn(8),lbn(12),sl_tn0,sl_gn0,sl_bn0,&
            ta_tn,ta_gn,ta_bn,ta_tn0,ta_gn0,ta_bn0,&
            mp_tn,mp_gn,mp_bn,np_tn,np_gn,np_bn)

call ch4_n2o_land_production(lbmp_ts0,lbmp_bs0,lbnp_ts0,lbnp_gs0,lbnp_bs0,&
            lbs(4),lbs(8),lbs(12),sl_ts0,sl_gs0,sl_bs0,&
            ta_ts,ta_gs,ta_bs,ta_ts0,ta_gs0,ta_bs0,&
            mp_ts,mp_gs,mp_bs,np_ts,np_gs,np_bs)
!-----------------------------------------------------------------



!-----------------------------------------------------------------
! Decay rates for all different vegetation types
!-----------------------------------------------------------------
call land_decay_rates(np_tn0,np_gn0,np_bn0,&
     gr_tn0,wo_tn0,li_tn0,sl_tn0,lbmp_tn0,&
     gr_gn0,wo_gn0,li_gn0,sl_gn0,lbmp_gn0,&
     gr_bn0,wo_bn0,li_bn0,sl_bn0,lbmp_bn0,&
     ta_tn,ta_tn0,ta_gn,ta_gn0,ta_bn,ta_bn0,&
     agr_tn,awo_tn,ali_tn,asl_tn,&
     agr_gn,awo_gn,ali_gn,asl_gn,&
     agr_bn,awo_bn,ali_bn,asl_bn)

call land_decay_rates(np_ts0,np_gs0,np_bs0,&
     gr_ts0,wo_ts0,li_ts0,sl_ts0,lbmp_ts0,&
     gr_gs0,wo_gs0,li_gs0,sl_gs0,lbmp_gs0,&
     gr_bs0,wo_bs0,li_bs0,sl_bs0,lbmp_bs0,&
     ta_ts,ta_ts0,ta_gs,ta_gs0,ta_bs,ta_bs0,&
     agr_ts,awo_ts,ali_ts,asl_ts,&
     agr_gs,awo_gs,ali_gs,asl_gs,&
     agr_bs,awo_bs,ali_bs,asl_bs)
!-----------------------------------------------------------------
pco2 = (/p12co2n,p13co2n,p14co2n/)
call land_change_rate(nptfn,npgsn,npbfn,&
     agr_tn,agr_gn,agr_bn,&
     awo_tn,awo_gn,awo_bn,&
     ali_tn,ali_gn,ali_bn,&
     asl_tn,asl_gn,asl_bn,&
     mp_tn ,mp_gn,mp_bn,lbn,pco2,aln)

pco2 = (/p12co2s,p13co2s,p14co2s/)
call land_change_rate(nptfs,npgss,npbfs,&
     agr_ts,agr_gs,agr_bs,&
     awo_ts,awo_gs,awo_bs,&
     ali_ts,ali_gs,ali_bs,&
     asl_ts,asl_gs,asl_bs,&
     mp_ts ,mp_gs,mp_bs,lbs,pco2,als)
!-----------------------------------------------------------------


!-----------------------------------------------------------------
! Calculate sink/source to atmosphere from changes in LB (in mol/s)
! for CO2 and CH4 for each isotope
!-----------------------------------------------------------------
pco2 = (/p12co2n,p13co2n,p14co2n/)
call land_source_sink(mp_bn,mp_gn,mp_tn,np_bn,np_gn,np_tn,npbfn,npgsn,nptfn,&
     ali_bn,ali_gn,ali_tn,asl_bn,asl_gn,asl_tn,lbn,pco2,&
     alssn)

pco2 = (/p12co2s,p13co2s,p14co2s/)
call land_source_sink(mp_bs,mp_gs,mp_ts,np_bs,np_gs,np_ts,npbfs,npgss,nptfs,&
     ali_bs,ali_gs,ali_ts,asl_bs,asl_gs,asl_ts,lbs,pco2,&
     alsss)
!-----------------------------------------------------------------


!-----------------------------------------------------------------
! Fluxes for each atmospheric box
!-----------------------------------------------------------------
do i=3,8
call land_flx_nh(alssn(i,:),ftropn,fgrassn,fsn_n,alenx,almnx,alhnx)
alen(i) = alenx; almn(i) = almnx; alhn(i) = alhnx
call land_flx_sh(alsss(i,:),ftrops,fgrasss,alesx,almsx,alhsx)
ales(i) = alesx; alms(i) = almsx; alhs(i) = alhsx
end do

!-----------------------------------------------------------------
end subroutine land_bios

include 'Vegetation_Boundaries.f90'
include 'olf_lat_NH.f90'
include 'olf_lat_SH.f90'
include 'Land_New_Prod.f90'
include 'CH4_N2O_Land_Production.f90'
include 'Land_decay_rates.f90'
include 'Land_change_rate.f90'
include 'Land_Source_Sink.f90'
include 'Land_flx_NH.f90'
include 'Land_flx_SH.f90'