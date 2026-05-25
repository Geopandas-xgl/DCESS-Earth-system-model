!==========================================================================
! Main Atm
!==========================================================================
subroutine atm(t,at,to_a,to_p,to_ar,to_so,inso,&
               hsen,hlat,htot,h_fde,fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s,&
               tot_sw,tot_lw,ase_ar,ase_at,ase_pa,ase_so,&
               dfi_na,dfi_np,dfis,aoif_at,aoif_pa,aoif_ar,aoif_so,pta_nh,pta_sh,fsn_n,fsn_s)
! Calculates terms for the atmosphere energy balance
! Inputs:
! t: time (s)
! at: atmospheric tracers
! Surface Ocean temperatures: to_a,to_p,to_ar,to_so (°C)
! Outputs:
! Meridional heat (W) and freshwater (Sv) transport: hsen,hlat,htot,h_fde,fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s
! Total short and longwave radiation: tot_sw,tot_lw (W)
! Air-sea heat exchange: ase_ar,ase_at,ase_pa,ase_so (W)
! rate change of sea-ice lines: dfi_na,dfi_np,dfis (1/s)
! Ice-free ocean areas: aoif_at,aoif_pa,aoif_ar,aoif_so (m2)
! Terms for meridional Ta profile
! Snow lines: fsn_n,fsn_s (radians)

use parameters, only: wp,nta,fdh,fdl,fbr,b_north,bsh,nab,nob
implicit none

real(wp), intent(in) :: t,at(nta,nab),to_a(nob),to_p(nob),to_ar,to_so,inso(19,14)
real(wp), dimension(4), intent(out) :: hsen,hlat,htot
real(wp), intent(out) :: fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s,h_fde
real(wp), dimension(nab), intent(out) :: tot_sw,tot_lw
real(wp), intent(out) :: ase_ar,ase_at(nob),ase_pa(nob),ase_so,fsn_n,fsn_s
real(wp), intent(out) :: dfi_na,dfi_np,dfis,aoif_at(nob),aoif_pa(nob),aoif_ar,aoif_so
real(wp), dimension(3), intent(out) :: pta_nh,pta_sh
real(wp), dimension(nab) :: pch4,pn2o,pco2
real(wp) :: ta_n,ta_mn,ta_en,ta_es,ta_ms,ta_s
real(wp) :: to_na,to_mna,to_ena,to_esa,to_msa
real(wp) :: to_np,to_mnp,to_enp,to_esp,to_msp
real(wp) :: fsi_na,fsi_np,fsi_so
real(wp) :: ta_nh,ta_sh,dta_nh,dta_sh
real(wp) :: r1,r2,r3,r4,r5,r6,r7,r8,vegrat_nh,vegrat_sh
real(wp) :: to_si_np,to_si_na,to_si_so

! Redefine
pch4 = at(3,:)
pn2o = at(4,:)
pco2 = at(5,:)

ta_n  = at(1,1)
ta_mn = at(1,2)
ta_en = at(1,3)
ta_es = at(1,4)
ta_ms = at(1,5)
ta_s  = at(1,6)

to_na  = to_a(1)
to_mna = to_a(2)
to_enA = to_a(3)
to_esA = to_a(4)
to_msA = to_a(5)

to_np  = to_p(1)
to_mnp = to_p(2)
to_enp = to_p(3)
to_esp = to_p(4)
to_msp = to_p(5)

fsi_na = at(10,1)
fsi_np = at(10,2)
fsi_so = at(10,6)

ta_nh = ta_n*(1._wp-sin(fdh)) + ta_mn*(sin(fdh)-sin(fdl)) + ta_en*sin(fdl)
ta_sh = ta_s*(1._wp-sin(fdh)) + ta_ms*(sin(fdh)-sin(fdl)) + ta_es*sin(fdl)

dta_nh  = ta_nh - 15._wp
dta_sh  = ta_sh - 15._wp
r1 = 3.294e-08_wp; r2 = -3.220e-07_wp; r3 = -1.2790e-06_wp; r4 = 3.287e-05_wp;
r5 = 9.811e-05_wp; r6 = -0.0001164_wp; r7 = -0.03321_wp;    r8 = 1._wp;

vegrat_nh = r1*dta_nh**7 + r2*dta_nh**6 + r3*dta_nh**5 + r4*dta_nh**4 +&
     r5*dta_nh**3 + r6*dta_nh**2 + r7*dta_nh   + r8
vegrat_sh = r1*dta_sh**7 + r2*dta_sh**6 + r3*dta_sh**5 + r4*dta_sh**4 +&
     r5*dta_sh**3 + r6*dta_sh**2 + r7*dta_sh   + r8

! Establish the atmospheric temperature profile; Decoupled N/S
call get_poly_pta(ta_n,ta_mn,ta_en,ta_es,ta_ms,ta_s,pta_nh,pta_sh)

! Sea-ice model
to_si_np = seaice_temp(to_mnp,to_np,fsi_np,fdh)
to_si_na = seaice_temp(to_na ,to_ar,fsi_na,fbr)
to_si_so = to_so

call sea_ice(pta_nh,pta_sh,fsi_na,fsi_np,fsi_so,to_si_na,to_si_np,to_si_so,dfi_na,dfi_np,dfis)

! Snow line
call snowline_root(pta_nh,fsn_n,fsn_s)


! Shortwave radiation
call atm_shortwave(fsi_na,fsi_np,fsi_so,fsn_n,fsn_s,vegrat_nh,vegrat_sh,t,inso,tot_sw)

! Longwave radiation
call atm_longwave(pch4,pn2o,pco2,at(1,:),tot_lw)

! Heat and water vapor transport
call atm_heat_fw_transport(pta_nh,pta_sh,ta_en,ta_es,&
                    fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s,&
                    hsen,hlat,htot,h_fde)

! Air-sea heat exchange
call atm_airsea_heat_exc(fsi_na,fsi_np,fsi_so,pta_nh,pta_sh,&
                         ta_mn,ta_en,ta_es,ta_ms,&
                         to_ar,to_na,to_mna,to_ena,to_esa,to_msa,to_np,to_mnp,to_enp,to_esp,to_msp,to_so,&
                         ase_ar,ase_at,ase_pa,ase_so,aoif_at,aoif_pa,aoif_ar,aoif_so)

contains
real(wp) function seaice_temp(tw,tc,si,phi) result(to_si)
use parameters, only: to_f
real(wp),intent(in) :: tw,tc,si,phi
real(wp) :: amp,b

amp = tw-tc ! amplitude
b   = 1.e-2_wp
to_si = tc + amp/2._wp*(1._wp - tanh((si-phi)/b))
to_si = max(to_si,to_f)
end function
end subroutine atm
!==========================================================================

include 'Get_Poly_PTa.f90'
include 'Sea_Ice.f90'
include 'Snowline_root.f90'
include 'Atm_Shortwave.f90'
include 'Atm_Longwave.f90'
include 'Atm_Heat_Fw_Transport.f90'
include 'Atm_AirSea_heat_exc.f90'
