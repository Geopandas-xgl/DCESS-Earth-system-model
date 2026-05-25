subroutine atm_shortwave(fsi_na,fsi_np,fsi_so,fsn_n,fsn_s,vegrat_nh,vegrat_sh,t,inso,tot_sw)
! Calculates shortwave radiation for each atmospheric box
! Inputs:
! Sea ice lines: fsi_na,fsi_np,fsi_so (radians)
! Snow lines: fsn_n,fsn_s (radians)
! vegetation ratios: vegrat_nh,vegrat_sh
! t: time 
! inso: insolation data
! Outputs:
! tot_sw: total shortwave radiation (W)

use parameters, only : wp,h_area,bsh,dr,pi,fbr,fdh,fdl,&
                        olf_na,olf_mna,olf_ena,olf_esa,olf_msa,&
                        olf_np,olf_mnp,olf_enp,olf_esp,olf_msp,olf_ar,nab


implicit none
real(wp), intent(in) :: fsi_na,fsi_np,fsi_so,fsn_n,fsn_s,vegrat_nh,vegrat_sh,t,inso(19,14)
real(wp), dimension(nab), intent(out) :: tot_sw
real(wp) :: gamma,al0_n,al0_s,al2,ali,fg,sg,phii(180)
real(wp) :: olf_n,olf_mn,olf_en,olf_es,olf_ms
real(wp), dimension(90) :: phi,al_n,al_s,alsi_na,alsi_np,alsn_n,width,alsi_s,alsn_s,qf_n,qf_s,a1phi
real(wp) :: lsw_ar,lsw_hn,lsw_mn,lsw_en,lsw_es,lsw_ms,lsw_hs
real(wp) :: osw_ar,osw_hna,osw_mna,osw_ena,osw_hnp,osw_mnp,osw_enp,osw_es,osw_ms,osw_hs
integer :: i

gamma = 0.02_wp
al0_n  =  0.7_wp + gamma*(1-vegrat_nh)
al0_s  =  0.7_wp + gamma*(1-vegrat_sh)

al2     =  -0.175_wp
ali     =  .62_wp                ! Albedo of ice covered area

! Equatorward Extension of Glaciation (fg) 
fg       =  70.0_wp*dr                  ! [rad] 

! Calculate the poleward glaciation increase (sg), long per lat (rad/rad), 
! assuming total land glaciation  at the pole.

fg = fsn_n*merge(1._wp,0._wp,fg>fsn_n) + fg*merge(1,0,fg<=fsn_n)    ! [rad] 
sg = 2._wp*pi*(1._wp-olf_ar)/(pi/2._wp-fg)        ! [rad/rad]                                        

!  One degree latitude array [0.5:89.5], used for numerical integration of insolation.
phi = (/(real(i,wp)+0.5_wp, i=0,89,1 )/)*dr

! Calculate the latutudinal zonal mean albedo values
al_n   = 1._wp - al0_n - al2/2._wp*(3._wp*sin(phi)**2-1._wp)    ! ice free albedo
al_s   = 1._wp - al0_s - al2/2._wp*(3._wp*sin(phi)**2-1._wp)    ! ice free albedo

alsi_na = al_n*merge(1._wp,0._wp,phi<fsi_na) + ali*merge(1._wp,0._wp,phi>=fsi_na)    ! albedo north, including ice line
alsi_np = al_n*merge(1._wp,0._wp,phi<fsi_np) + ali*merge(1._wp,0._wp,phi>=fsi_np)    ! albedo north, including ice line
alsn_n  = al_n*merge(1._wp,0._wp,phi<fsn_n ) + ali*merge(1._wp,0._wp,phi>=fsn_n )    ! albedo north, including ice line

width  = sg*(phi-fg)                                        ! width of glaciation as function of latitude
where (width < 0._wp) width=0._wp
where (phi   < fg   ) width=0._wp
alsn_n = (alsn_n*(2._wp*pi-width) + ali*width)/(2._wp*pi)   ! albedo north, including ice line and glaciation
alsi_s  = al_s*merge(1._wp,0._wp,phi<fsi_so) + ali*merge(1._wp,0._wp,phi>=fsi_so)       ! albedo south, including ice line
alsn_s  = al_s*merge(1._wp,0._wp,phi<fsn_s ) + ali*merge(1._wp,0._wp,phi>=fsn_s )        ! albedo south, including ice line

call get_insolation(t,inso,qf_n,qf_s,phii)

! Area of one degree latitude bands
a1phi   = h_area*( sin(phi+.5_wp*dr)-sin(phi-.5_wp*dr) )

olf_n  = olf_na  + olf_np
olf_mn = olf_mna + olf_mnp
olf_en = olf_ena + olf_enp
olf_es = olf_esa + olf_esp
olf_ms = olf_msa + olf_msp

! Integrate insolation over area for each land sector
lsw_ar = sum( (a1phi*qf_n*(1._wp-alsn_n))*merge(1._wp,0._wp,phi>fBr  )                            )
lsw_hn = sum( (a1phi*qf_n*(1._wp-alsn_n))*merge(1._wp,0._wp,phi>fdh  )*merge(1._wp,0._wp,phi<fBr) )
lsw_mn = sum( (a1phi*qf_n*(1._wp-alsn_n))*merge(1._wp,0._wp,phi>fdl  )*merge(1._wp,0._wp,phi<fdh) )
lsw_en = sum( (a1phi*qf_n*(1._wp-alsn_n))*merge(1._wp,0._wp,phi>0._wp)*merge(1._wp,0._wp,phi<fdl) )
lsw_es = sum( (a1phi*qf_s*(1._wp-alsn_s))*merge(1._wp,0._wp,phi>0._wp)*merge(1._wp,0._wp,phi<fdl) )
lsw_ms = sum( (a1phi*qf_s*(1._wp-alsn_s))*merge(1._wp,0._wp,phi>fdl  )*merge(1._wp,0._wp,phi<fdh) )
lsw_hs = sum( (a1phi*qf_s*(1._wp-alsn_s))*merge(1._wp,0._wp,phi>bsh  )                            )

! Integrate insolation over area for each ocean sector
osw_ar  = sum( (a1phi*qf_n*(1._wp-alsi_na) )*merge(1._wp,0._wp,phi>fBr  )                            )
osw_hna = sum( (a1phi*qf_n*(1._wp-alsi_na) )*merge(1._wp,0._wp,phi>fdh  )*merge(1._wp,0._wp,phi<fBr) )
osw_mna = sum( (a1phi*qf_n*(1._wp-alsi_na) )*merge(1._wp,0._wp,phi>fdl  )*merge(1._wp,0._wp,phi<fdh) )
osw_ena = sum( (a1phi*qf_n*(1._wp-alsi_na) )*merge(1._wp,0._wp,phi>0._wp)*merge(1._wp,0._wp,phi<fdl) )

osw_hnp = sum( (a1phi*qf_n*(1._wp-alsi_np))*merge(1._wp,0._wp,phi>fdh  )*merge(1._wp,0._wp,phi<fBr) )
osw_mnp = sum( (a1phi*qf_n*(1._wp-alsi_np))*merge(1._wp,0._wp,phi>fdl  )*merge(1._wp,0._wp,phi<fdh) )
osw_enp = sum( (a1phi*qf_n*(1._wp-alsi_np))*merge(1._wp,0._wp,phi>0._wp)*merge(1._wp,0._wp,phi<fdl) )

osw_es  = sum( (a1phi*qf_s*(1._wp-alsi_s ))*merge(1._wp,0._wp,phi>0._wp)*merge(1._wp,0._wp,phi<fdl) )
osw_ms  = sum( (a1phi*qf_s*(1._wp-alsi_s ))*merge(1._wp,0._wp,phi>fdl  )*merge(1._wp,0._wp,phi<fdh) )
osw_hs  = sum( (a1phi*qf_s*(1._wp-alsi_s ))*merge(1._wp,0._wp,phi>fdh  )*merge(1._wp,0._wp,phi<bsh) )


tot_sw(1) = olf_na *osw_hna + olf_np *osw_hnp  + olf_ar*osw_ar + (1._wp-olf_ar)*lsw_ar + (1._wp-olf_n)*lsw_hn
tot_sw(2) = olf_mna*osw_mna + olf_mnp*osw_mnp                  + (1._wp-olf_mn)*lsw_mn
tot_sw(3) = olf_ena*osw_ena + olf_enp*osw_enp                  + (1._wp-olf_en)*lsw_en
tot_sw(4) = olf_es *osw_es                                     + (1._wp-olf_es)*lsw_es
tot_sw(5) = olf_ms *osw_ms                                     + (1._wp-olf_ms)*lsw_ms
tot_sw(6) =         osw_hs                                     +                lsw_hs


end subroutine atm_shortwave

include 'Get_Insolation.f90'