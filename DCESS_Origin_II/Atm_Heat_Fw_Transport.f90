subroutine atm_heat_fw_transport(pta_n,pta_s,ta_en,ta_es,&
        fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s,&
        hsen,hlat,htot,h_fde)
! Calculates heat (sensible and latent) and freshwater transport across atmopsheric boundaries
! Inputs:
! atmospheric terms for meridional profile of Ta (pta_n,pta_s)
! atmospheric temperatures of tropical zones
! Outputs:
! Heat (W) and freshwater (Sv) transports at every atmospheric boundary


use parameters, only: wp,fde,fbr,fdh,fdl

implicit none
real(wp), intent(in) :: pta_n(3),pta_s(3),ta_en,ta_es
real(wp), intent(out) :: fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s,h_fde
real(wp), intent(out), dimension(4) :: hsen,hlat,htot
real(wp) :: hs_fdh_n,hs_fdl_n,hs_fdl_s,hs_fdh_s
real(wp) :: hl_fdh_n,hl_fdl_n,hl_fdl_s,hl_fdh_s
real(wp) :: ty_fbr,ty_fdh_n,ty_fdl_n,ty_fdl_s,ty_fdh_s
real(wp) :: ta_fbr,ta_fdh_n,ta_fdl_n,ta_fde_n,ta_fde_s,ta_fdl_s,ta_fdh_s
real(wp) :: nl,nh,cfwl,cfwh,cfweq,csl,csh,cseq,lv,teq_pi
real(wp), external :: dtdy_at_phi,ta_at_phi
nl     = 2.5_wp                   ! Gradient exponent
nh     = 1.7_wp
cfwl   = 1.56e10_wp             ! Freshwater transport
cfwh   = 14.66_wp*cfwl
csl    = 4.21e11_wp       ! Sensible heat transport 
csh    = 7.85_wp*csl
lv     =  2.25e9_wp         ! Latent heat of vaporization [J/1000kg] = [J/m^3]   
cfweq  = -1.6e15_wp/lv
cseq   = 2.4e15_wp
teq_pi = 29._wp

ty_fbr   = dtdy_at_phi(pta_n,fbr)
ty_fdh_n = dtdy_at_phi(pta_n,fdh)
ty_fdl_n = dtdy_at_phi(pta_n,fdl)
ty_fdl_s = dtdy_at_phi(pta_s,fdl)
ty_fdh_s = dtdy_at_phi(pta_s,fdh)

ta_fbr   = ta_at_phi(pta_n,fbr)
ta_fdh_n = ta_at_phi(pta_n,fdh)
ta_fdl_n = ta_at_phi(pta_n,fdl)
ta_fde_n = ta_at_phi(pta_n,fde)
ta_fde_s = ta_at_phi(pta_s,fde)
ta_fdl_s = ta_at_phi(pta_s,fdl)
ta_fdh_s = ta_at_phi(pta_s,fdh)

! Freshwater transport [m3/s]
fw_fbr   = fw_flux(ta_fbr  ,ty_fbr  ,nh,0.70_wp*cfwh)
fw_fdh_n = fw_flux(ta_fdh_n,ty_fdh_n,nh,cfwh     )
fw_fdl_n = fw_flux(ta_fdl_n,ty_fdl_n,nl,cfwl     )
fw_fdl_s = fw_flux(ta_fdl_s,ty_fdl_s,nl,cfwl     )
fw_fdh_s = fw_flux(ta_fdh_s,ty_fdh_s,nh,cfwh     )

! Heat transport - Sensible [W]
hs_fdh_n = sens_heat(ty_fdh_n,nh,csh)
hs_fdl_n = sens_heat(ty_fdl_n,nl,csl)
hs_fdl_s = sens_heat(ty_fdl_s,nl,csl)
hs_fdh_s = sens_heat(ty_fdh_s,nh,csh)

! Heat transport - Latent [W]
hl_fdh_n = lv*fw_fdh_n
hl_fdl_n = lv*fw_fdl_n
hl_fdl_s = lv*fw_fdl_s
hl_fdh_s = lv*fw_fdh_s


hsen = (/hs_fdh_s,hs_fdl_s,hs_fdl_n,hs_fdh_n/)
hlat = (/hl_fdh_s,hl_fdl_s,hl_fdl_n,hl_fdh_n/)
htot = hsen + hlat

h_fde = cseq*(ta_es-ta_en) + 1.1_wp +&
(lv*cfweq*(ta_es-ta_en)-1.1_wp)*exp(-5420._wp/((ta_fde_n+ta_fde_s)/2._wp+273._wp))/exp(-5420._wp/(teq_pi+273._wp))

! write(*,'(4e25.15)') hsen
! write(*,'(4e25.15)') hlat
! write(*,'(4e25.15)') htot
! write(*,'(1e25.15)') h_fde
! write(*,'(5e25.15)') fw_fbr,fw_fdh_n,fw_fdl_n,fw_fdl_s,fw_fdh_s

! print*, 'heat_transport'

!==========================================================================
contains

real(wp) function fw_flux(ta,ty,n,cfw) result(fw)
implicit none
real(wp), intent(in) :: ta,ty,n,cfw
fw = cfw*exp(-5420._wp/(ta+273._wp))*abs(ty)**n
end function

real(wp) function sens_heat(ty,n,cs) result(hs)
implicit none
real(wp), intent(in) :: ty,n,cs
hs = cs * abs(ty)**n
end function

end subroutine atm_heat_fw_transport

include 'dtdy_at_phi.f90'
! include 'Ta_at_phi.f90'