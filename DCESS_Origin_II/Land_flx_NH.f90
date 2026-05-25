subroutine land_flx_nh(alss,fitf,figs,fsn,&
                      all,alm,alh)
! Calculates land fluxes of GHG for NH

use parameters, only: wp,fdh,fdl
implicit none
real(wp), intent(in) :: alss(3),fitf,figs,fsn
real(wp), intent(out) :: all,alm,alh
real(wp) :: fl_bf,fl_gs,fl_tf,wl,wm,wh,wlgs,wmgs,wmbf,whbf
! alss: [BF GSD TF]

fl_bf = alss(1)
fl_gs = alss(2)
fl_tf = alss(3)

! Case 1
if (fitf < fdl .AND. figs <= fdl .AND. fsn > fdh) then
    wl = ( sin(fdl)-sin(figs) )/( sin(fsn)-sin(figs) )
    wm = ( sin(fdh)-sin(fdl)  )/( sin(fsn)-sin(figs) )
    wh = ( sin(fsn)-sin(fdh)  )/( sin(fsn)-sin(figs) )

    all = fl_tf + fl_gs + fl_bf*wl
    alm = fl_bf*wm
    alh = fl_bf*wh

! Case 2
elseif (fitf < fdl .AND. figs <= fdl .AND. fsn > fdl .AND. fsn <= fdh) then
    wl = ( sin(fdl)-sin(figs) )/( sin(fsn)-sin(figs) )
    wm = ( sin(fsn)-sin(fdl)  )/( sin(fsn)-sin(figs) )

    all = fl_tf + fl_gs + fl_bf*wl
    alm = fl_bf*wm
    alh = 0._wp
! Case 3
elseif (fitf < fdl .AND. figs > fdl .AND. figs < fdh .AND. fsn <= fdh) then
    wl = ( sin(fdl)-sin(fitf) )/( sin(figs)-sin(fitf) )
    wm = ( sin(figs)-sin(fdl) )/( sin(figs)-sin(fitf) )

    all = fl_tf + fl_gs*wl
    alm = fl_gs*wm + fl_bf
    alh = 0._wp
! Case 4
elseif (fitf < fdl .AND. figs > fdl .AND. figs < fdh .AND. fsn > fdh) then
    wlgs = ( sin(fdl)-sin(fitf) )/( sin(figs)-sin(fitf) )
    wmgs = ( sin(figs)-sin(fdl) )/( sin(figs)-sin(fitf) )
    wmbf = ( sin(fdh)-sin(figs) )/( sin(fsn) -sin(figs) )
    whbf = ( sin(fsn)-sin(fdh)  )/( sin(fsn) -sin(figs) )

    all = fl_tf + fl_gs*wlgs
    alm = fl_gs*wmgs + fl_bf*wmbf
    alh = fl_bf*whbf
end if
end subroutine land_flx_nh