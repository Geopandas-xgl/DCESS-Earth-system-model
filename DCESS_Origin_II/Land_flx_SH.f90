subroutine land_flx_sh(alss,fitf,figs,&
                      all,alm,alh)
! Calculates land fluxes of GHG for SH

use parameters, only: wp,fdl,fdh
implicit none
real(wp), intent(in) :: alss(3),fitf,figs
real(wp), intent(out) :: all,alm,alh
real(wp) :: fl_bf,fl_gs,fl_tf,wl,wm

fl_bf = alss(1)
fl_gs = alss(2)
fl_tf = alss(3)

! Case 1
if (fiTF < fdl .AND. figs <= fdl) then
    wl = ( sin(fdl)-sin(figs) )/( sin(fdh)-sin(figs) )
    wm = ( sin(fdh)-sin(fdl)  )/( sin(fdh)-sin(figs) )

    all = fl_tf + fl_gs + fl_bf*wl
    alm = fl_bf*wm
    alh = 0._wp

! Case 2
elseif (fitf < fdl .AND. figs > fdl .AND.  figs < fdh) then
    wl = ( sin(fdl)-sin(fitf) )/( sin(figs)-sin(fitf) )
    wm = ( sin(figs)-sin(fdl) )/( sin(figs)-sin(fitf) )

    all = fl_tf + fl_gs*wl
    alm = fl_gs*wm + fl_bf
    alh = 0._wp

end if
end subroutine land_flx_sh