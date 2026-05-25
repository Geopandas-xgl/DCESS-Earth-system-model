subroutine sea_ice(pta_nh,pta_sh,fsi_na,fsi_np,fsi_so,to_na,to_np,to_so,dfi_na,dfi_np,dfis)
! Calculates rate change of sea-ice lines
! Inputs:
! factors for meridional profile of atmospheric temperature: pta_nh, pta_sh
! sea ice lines: fsi_na,fsi_np,fsi_so
! ocen surface temperatures: to_na,to_np,to_so
! Outputs:
! change rates of sea ice lines: dfi_na (North Atlantic/Arctic), dfi_np (North Pacific), dfis (Southern Ocean)

use parameters, only: wp,rc,d,to_f,ko,pi,b_north,bsh,fdl,fbr

implicit none
real(wp), intent(in) :: pta_nh(3),pta_sh(3),fsi_na,fsi_np,fsi_so,to_na,to_np,to_so
real(wp), intent(out) :: dfi_na,dfi_np,dfis

real(wp) :: ta_i_na,ta_i_np,ta_i_s
real(wp) :: rho_ice,l_ice,d_ice,k_ice,k_io
real(wp) :: iqasna,iqasnp,iqass,iqina,iqinp,iqis
real(wp) :: logi1,logi2,logi3,logi4
real(wp), external :: ta_at_phi

rho_ice = 917._wp           ! Density of ice
l_ice   = 3.34e5_wp         ! Melting heat 
d_ice   = 2._wp             ! Scale thickness of ice
k_ice   = 2._wp             ! Thermal diffusivity of ice 
k_io    = 150._wp           ! Ice-Ocean energy transfer

! Atmospheric temperature at ice edge
ta_i_na = ta_at_phi(pta_nh,fsi_na)
ta_i_np = ta_at_phi(pta_nh,fsi_np)
ta_i_s  = ta_at_phi(pta_sh,fsi_so)

iqasna = 1._wp/(rc*d*(to_na - to_f) + rho_ice*l_ice*d_ice)*ko*(ta_i_na-to_f)
iqasnp = 1._wp/(rc*d*(to_np - to_f) + rho_ice*l_ice*d_ice)*ko*(ta_i_np-to_f)
iqass  = 1._wp/(rc*d*(to_so - to_f) + rho_ice*l_ice*d_ice)*ko*(ta_i_s -to_f)


iqina = -k_ice/(d_ice**2*l_ice*rho_ice)*(to_f-ta_i_na)+k_io/(rho_ice*l_ice*d_ice)*(to_na - to_f)
iqinp = -k_ice/(d_ice**2*l_ice*rho_ice)*(to_f-ta_i_np)+k_io/(rho_ice*l_ice*d_ice)*(to_np - to_f)
iqis  = -k_ice/(d_ice**2*l_ice*rho_ice)*(to_f-ta_i_s )+k_io/(rho_ice*l_ice*d_ice)*(to_so - to_f)


if (iqasna <  0._wp) then;logi1=1._wp;else;logi1=0._wp;end if
if (fsi_na >= fdl  ) then;logi2=1._wp;else;logi2=0._wp;end if
if (iqina  >  0._wp) then;logi3=1._wp;else;logi3=0._wp;end if
if (fsi_na <= b_north   ) then;logi4=1._wp;else;logi4=0._wp;end if
dfi_na =  iqasna*logi1*logi2 + iqina*logi3*logi4

if (iqasnp <  0._wp) then;logi1=1._wp;else;logi1=0._wp;end if
if (fsi_np >= fdl  ) then;logi2=1._wp;else;logi2=0._wp;end if
if (iqinp  >  0._wp) then;logi3=1._wp;else;logi3=0._wp;end if
if (fsi_np <= fbr  ) then;logi4=1._wp;else;logi4=0._wp;end if
dfi_np =  iqasnp*logi1*logi2 + iqinp*logi3*logi4

if (iqass  <  0._wp) then;logi1=1._wp;else;logi1=0._wp;end if
if (fsi_so >= fdl  ) then;logi2=1._wp;else;logi2=0._wp;end if
if (iqis   >  0._wp) then;logi3=1._wp;else;logi3=0._wp;end if
if (fsi_so <= bsh  ) then;logi4=1._wp;else;logi4=0._wp;end if
dfis  =  iqass *logi1*logi2 + iqis *logi3*logi4




end subroutine sea_ice

include 'Ta_at_phi.f90'