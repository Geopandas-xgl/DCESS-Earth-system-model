subroutine oce(arc,soc,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,slso,&
        fw_fbr,fw_fdhn,fw_fdln,fw_fdls,fw_fdhs,&
        oc_pres,shl_dat,&
        kv_hna,kv_mna,kv_ena,kv_esa,kv_msa,kv_hnp,kv_mnp,kv_enp,kv_esp,kv_msp,kv_ar,kv_so,&
        kh_fbra,kh_fdhna,kh_fdlna,kh_fdea,kh_fdlsa,kh_fdhsa,kh_fdhnp,kh_fdlnp,kh_fdep,kh_fdlsp,kh_fdhsp,kh_msap,&
        q_fbr,q_fdhna,q_fdhnp,q_fdlna,q_fdlnp,q_fdea,q_fdep,q_fdlsa,q_fdlsp,q_fdhsa,q_fdhsp,&
        ievfsa,idlsa,ievfsp,idlsp,ievfbr,idlar,fdsup,feso,ent_lim,rfds_v,&
        w_hna,w_mna,w_ena,w_esa,w_msa,w_hnp,w_mnp,w_enp,w_esp,w_msp,w_so,w_ar,qshl)

! Calculates ocean dynamics
! Inputs:
! Ocean tracers: arc,soc,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,slso
! freshwater fluxes: fw_fbr,fw_fdhn,fw_fdln,fw_fdls,fw_fdhs (Sv)
! Outputs:
! kv: vertical diffusion
! kh: horizontal diffusion
! q : meridional advection
! w : vertical advection
! iev: northward volume fluxes crossing 55 °S and 65°N (Bering Strait)
! fdsup: volume flux leaving the ocean shelf
! feso: entrainment flux
! rfds_v: return flux to the ocean shelf
! qshl: heat release to the atmosphere from ocean shelf

use parameters, only: wp,nto,n,g,zmid,kh_s,kh_d,z_kh,d,pi,fbr,fdh,fdl,fde,n_ds,n_dp,e_r,rho0,rq,&
               nz_bot,nz_all,nz_upp,&
               ly_Ar,ly_n,ly_mn,ly_en,ly_es,ly_ms,ly_s,&
               w_65na,w_55na,w_35na,w_eqa,w_35sa,w_55sa,&
               w_55np,w_35np,w_eqp,w_35sp,w_55sp,&
               c_ar,c_na,c_mna,c_msa,c_sa,c_np,c_mnp,c_msp,c_sp,c_arx,&
               brvf,evfn,fw_z,evfsa,evfsp,fice,fgla_sh,fgla_so,rc,shl,to_f

implicit none
real(wp), dimension(nto,n), intent(in) :: arc,soc,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp
real(wp), intent(in) :: fw_fbr,fw_fdhn,fw_fdln,fw_fdls,fw_fdhs,oc_pres(n,7),shl_dat(nz_all,3),slso(nto)
real(wp), dimension(n) :: p_ar,p_hn,p_mn,p_en,p_es,p_ms,p_so!,p_sh
real(wp), dimension(n) :: rho_hna,rho_mna,rho_ena,rho_esa,rho_msa,&
                   rho_hnp,rho_mnp,rho_enp,rho_esp,rho_msp,rho_ar,rho_so
real(wp), dimension(n-1), intent(out) :: kv_hna,kv_mna,kv_ena,kv_esa,kv_msa,&
                                         kv_hnp,kv_mnp,kv_enp,kv_esp,kv_msp,kv_ar,kv_so
real(wp), dimension(n), intent(out) :: kh_fbra,kh_fdhna,kh_fdlna,kh_fdea,kh_fdlsa,kh_fdhsa,&
                                kh_fdhnp,kh_fdlnp,kh_fdep,kh_fdlsp,kh_fdhsp,kh_msap
real(wp), dimension(n), intent(out) :: q_fbr,q_fdhna,q_fdhnp,q_fdlna,q_fdlnp,q_fdea,q_fdep,q_fdlsa,q_fdlsp,q_fdhsa,q_fdhsp
real(wp), dimension(n), intent(out) :: ievfsa,ievfsp,ievfbr
integer, intent(out) :: idlsa,idlsp,idlar
real(wp), intent(out) :: fdsup,feso(n),qshl
integer, intent(out) :: ent_lim(2)

real(wp), dimension(n), intent(out) :: rfds_v,w_hna,w_mna,w_ena,w_esa,w_msa,w_hnp,w_mnp,w_enp,w_esp,w_msp,w_so,w_ar
real(wp), dimension(n) :: w_so_a,w_so_b,cs1,cs2,cs3

real(wp), dimension(n) :: kh, kh_lows,kh_lown,cum,fds_v,zerop
real(wp), dimension(n) :: py_fbr,py_fdhsa,py_fdhsp
real(wp), dimension(n) :: vy_fbr,vy_fdhna,vy_fdhnp,vy_fdlna,vy_fdlnp,vy_fdea,vy_fdep,vy_fdlsa,vy_fdlsp,vy_fdhsa,vy_fdhsp
real(wp) :: dy_fbr,dy_fdhn,dy_fdln,dy_fde,dy_fdls,dy_fdhs,vsur,vbot,vel(n),dy,fa,fp,n_1,n_nb,mtso
integer, dimension(:), allocatable :: izl
integer :: i
real(wp) :: z_shl(nz_bot),alfi(nz_bot),ctoi(nz_bot),saoi(nz_bot),pds(nz_bot),buo0,irq0,rsoi(nz_all),rds0(nz_upp)


p_ar = oc_pres(:,1)
p_hn = oc_pres(:,2)
p_mn = oc_pres(:,3)
p_en = oc_pres(:,4)
p_es = oc_pres(:,5)
p_ms = oc_pres(:,6)
p_so = oc_pres(:,7)


! Density and vertical diffusion profiles
call rho_vert_conv(hlna,p_hn,rho_hna,kv_hna)
call rho_vert_conv(mlna,p_mn,rho_mna,kv_mna)
call rho_vert_conv(llna,p_en,rho_ena,kv_ena)
call rho_vert_conv(llsa,p_es,rho_esa,kv_esa)
call rho_vert_conv(mlsa,p_ms,rho_msa,kv_msa)

call rho_vert_conv(hlnp,p_hn,rho_hnp,kv_hnp)
call rho_vert_conv(mlnp,p_mn,rho_mnp,kv_mnp)
call rho_vert_conv(llnp,p_en,rho_enp,kv_enp)
call rho_vert_conv(llsp,p_es,rho_esp,kv_esp)
call rho_vert_conv(mlsp,p_ms,rho_msp,kv_msp)

call rho_vert_conv(arc,p_ar,rho_ar,kv_ar)
call rho_vert_conv(soc,p_so,rho_so,kv_so)

! dy
dy_fbr  = ly_Ar - ly_n
dy_fdhn = ly_n  - ly_mn
dy_fdln = ly_mn - ly_en
dy_fde  = ly_en - ly_es
dy_fdls = ly_es - ly_ms
dy_fdhs = ly_ms - ly_s 

! Horizontal diffusion
kh      = kh_d +        kh_s*exp(-zmid/z_kh)      ! Combined horizontal diffusion [m2/s]
kh_lows = kh_d + 0.1_wp*kh_s*exp(-zmid/z_kh)
kh_lown = kh_d + 0.5_wp*kh_s*exp(-zmid/z_kh)

kh_fbra  = kh      * d*2._wp*pi*w_65na*cos(fbr)/dy_fbr
kh_fdhna = kh      * d*2._wp*pi*w_55na*cos(fdh)/dy_fdhn
kh_fdlna = kh      * d*2._wp*pi*w_35na*cos(fdl)/dy_fdln
kh_fdea  = kh      * d*2._wp*pi*w_eqa *cos(fde)/dy_fde
kh_fdlsa = kh_lows * d*2._wp*pi*w_35sa*cos(fdl)/dy_fdls
kh_fdhsa = kh_lows * d*2._wp*pi*w_55sa*cos(fdh)/dy_fdhs

kh_fdhnp = kh_lown * d*2._wp*pi*w_55np*cos(fdh)/dy_fdhn
kh_fdlnp = kh_lown * d*2._wp*pi*w_35np*cos(fdl)/dy_fdln
kh_fdep  = kh      * d*2._wp*pi*w_eqp *cos(fde)/dy_fde
kh_fdlsp = kh_lows * d*2._wp*pi*w_35sp*cos(fdl)/dy_fdls
kh_fdhsp = kh_lows * d*2._wp*pi*w_55sp*cos(fdh)/dy_fdhs

kh_fbra(n_ds+1:n) = 0._wp ! There is no flux below "Denmark Strait"


! Zonal mixing between MLSA/MLSP
vsur    = 23.e-2_wp ! m/s
vbot    = 10.e-2_wp ! m/s
vel     = vbot + (vsur - vbot)*exp(-zmid/z_kh)   ! m/s
dy      = 2222.e3_wp ! distance m
kh_msap = 0.08_wp*vel*dy*d    ! m3/s (sum(KhMLSAP) ~ 100 Sv)


! Meridional flow, Rayleigh friction
fa=1._wp
fp=1._wp-fa
n_1 = real(n,wp)-1._wp

cum = cumsum((rho_ar - rho_hna)*d,1,n)
py_fbr         = g*( cum - (rho_ar - rho_hna)*.5_wp*d)/( dy_fbr*e_r )
vy_fbr         = -1._wp/(rho0*rq)*py_fbr
vy_fbr(1)      = 0._wp
vy_fbr(2:n_ds) = vy_fbr(2:n_ds) - sum(vy_fbr(2:n_ds))/size(vy_fbr(2:n_ds))
q_fbr(1:n)     = 0._wp
q_fbr(2:n_ds)  = vy_fbr(2:n_ds)*d*2._wp*pi*w_65na*e_r*cos(fbr) + evfn/(real(n_ds,wp)-1._wp) - &
          c_ar*fw_fbr/(real(n_ds,wp)-1._wp) - fa*c_arx*fw_fdhn/(real(n_ds,wp)-1._wp) - &
          brvf/(real(n_ds,wp)-1._wp) - fp*c_arx*fw_fdhn/(real(n_ds,wp)-1._wp)
q_fbr(1)       = -evfn

call med_flux(rho_hna,rho_mna,dy_fdhn,w_55na,fdh,vy_fdhna)
q_fdhna(2:n) = vy_fdhna(2:n) - c_na*fw_fdhn/n_1 - brvf/n_1 - fp*c_arx*fw_fdhn/n_1
q_fdhna(1)   = 0._wp

call med_flux(rho_hnp,rho_mnp,dy_fdhn,w_55np,fdh,vy_fdhnp)
q_fdhnp(2:n) = vy_fdhnp(2:n) - c_np*fw_fdhn/n_1 + brvf/n_1 + fp*c_arx*fw_fdhn/n_1
q_fdhnp(1)   = 0._wp

call med_flux(rho_mna,rho_ena,dy_fdln,w_35na,fdl,vy_fdlna)
q_fdlna(2:n) = vy_fdlna(2:n) - c_mna*fw_fdln/n_1  - brvf/n_1 - fp*c_arx*fw_fdhn/n_1
q_fdlna(1)   = 0._wp

call med_flux(rho_mnp,rho_enp,dy_fdln,w_35np,fdl,vy_fdlnp)
q_fdlnp(2:n)  = vy_fdlnp(2:n) - c_mnp*fw_fdln/n_1 + brvf/n_1 + fp*c_arx*fw_fdhn/n_1
q_fdlnp(1)   = 0._wp

call med_flux(rho_ena,rho_esa,dy_fde,w_eqa,fde,vy_fdea)
q_fdea = vy_fdea + fw_z/n_1 - brvf/n_1 - fp*c_arx*fw_fdhn/n_1

call med_flux(rho_enp,rho_esp,dy_fde,w_eqp,fde,vy_fdep)
q_fdep = vy_fdep - fw_z/n_1 + brvf/n_1 + fp*c_arx*fw_fdhn/n_1

call med_flux(rho_esa,rho_msa,dy_fdls,w_35sa,fdl,vy_fdlsa)
q_fdlsa(2:n) = vy_fdlsa(2:n) + c_msa*fw_fdls/n_1 + fw_z/n_1 - brvf/n_1 - fp*c_arx*fw_fdhn/n_1
q_fdlsa(1)    = 0._wp

call med_flux(rho_esp,rho_msp,dy_fdls,w_35sp,fdl,vy_fdlsp)
q_fdlsp(2:n)  = vy_fdlsp(2:n) + c_msp*fw_fdls/n_1 - fw_z/n_1 + brvf/n_1 + fp*c_arx*fw_fdhn/n_1
q_fdlsp(1)    = 0._wp

n_nb = real(n,wp)-real(n_dp,wp)
cum = cumsum((rho_msa - rho_so)*d,1,n)
py_fdhsa           = g*( cum - (rho_msa - rho_so)*.5_wp*d)/( dy_fdhs*e_r )
vy_fdhsa           = -1._wp/(rho0*rq)*py_fdhsa
vy_fdhsa(1)        = 0._wp
vy_fdhsa(n_dp+1:n) = vy_fdhsa(n_dp+1:n) - sum(vy_fdhsa(n_dp+1:n))/size(vy_fdhsa(n_dp+1:n))
q_fdhsa(n_dp+1:n)  = vy_fdhsa(n_dp+1:n)*d*2._wp*pi*w_55sa*e_r*cos(fdh) - evfsa/n_nb + &
              c_sa*fw_fdhs/n_nb + fw_z/n_nb - brvf/n_nb - fp*c_arx*fw_fdhn/n_nb

cum = cumsum((rho_msp - rho_so)*d,1,n)
py_fdhsp           = g*( cum - (rho_msp - rho_so)*.5_wp*d)/( dy_fdhs*e_r )
vy_fdhsp           = -1._wp/(rho0*rq)*py_fdhsp
vy_fdhsp(1)        = 0._wp
vy_fdhsp(n_dp+1:n) = vy_fdhsp(n_dp+1:n) - sum(vy_fdhsp(n_dp+1:n))/size(vy_fdhsp(n_dp+1:n))
q_fdhsp(n_dp+1:n)  = vy_fdhsp(n_dp+1:n)*d*2._wp*pi*w_55sp*e_r*cos(fdh) - evfsp/n_nb + &
              c_sp*fw_fdhs/n_nb - fw_z/n_nb + brvf/n_nb + fp*c_arx*fw_fdhn/n_nb

! Southern Hemisphere Ekman transport
call ekman(rho_so ,rho_msa,evfsa,ievfsa,idlsa)
call ekman(rho_so ,rho_msp,evfsp,ievfsp,idlsp)
call ekman(rho_hnp,rho_ar ,brvf ,ievfbr,idlar)

! Downsloping flow
call init_downsloping(soc,slso,shl_dat,z_shl,alfi,buo0,irq0,ctoi,saoi,rsoi,rds0,pds)
call entrainment(slso(2),buo0,irq0,rds0,alfi,ctoi,saoi,pds,z_shl,rsoi,fdsup,feso,ent_lim)
! print*, sum(feso)


allocate(izl(ent_lim(2)))
izl = (/ (i, i=ent_lim(1)-ent_lim(2)+1,ent_lim(1)) /)
fds_v=0._wp
fds_v(izl) = (fdsup+sum(feso))/ent_lim(2) ! m3/s
deallocate(izl)
! print*, fdsup
rfds_v = 0._wp
rfds_v(1:shl) = (fdsup+fice-fgla_sh+fgla_so)/shl ! Shelf return flux

mtso = sum(soc(1,1:shl))/shl
qshl = rc*( fice*(-to_f) + (fdsup+fice)* (mtso - to_f) )

! Vertical advection
zerop = 0._wp

if (idlar .eq. 1) then
w_ar  = vert_adv(q_fbr          ,zerop )
else
w_ar  = vert_adv(q_fbr + ievfbr ,zerop )
end if

w_hna  = vert_adv(q_fdhna        ,q_fbr  )
w_mna  = vert_adv(q_fdlna        ,q_fdhna)
w_ena  = vert_adv(q_fdea         ,q_fdlna)
w_esa  = vert_adv(q_fdlsa        ,q_fdea )
w_msa  = vert_adv(q_fdhsa+ievfsa ,q_fdlsa)

w_hnp  = vert_adv(q_fdhnp        ,zerop  )
w_mnp  = vert_adv(q_fdlnp        ,q_fdhnp)
w_enp  = vert_adv(q_fdep         ,q_fdlnp)
w_esp  = vert_adv(q_fdlsp        ,q_fdep )
w_msp  = vert_adv(q_fdhsp+ievfsp ,q_fdlsp)


w_so_a = vert_adv(zerop,q_fdhsa)
w_so_b = vert_adv(zerop,q_fdhsp)
w_so   = w_so_a + w_so_b

cs1 = cumsum(fds_v ,n,1)
cs2 = cumsum(rfds_v,n,1)
cs3 = cumsum(feso  ,n,1)

w_so(n:1:-1) = w_so(n:1:-1) + (/0._wp,cs1(1:n-1)-cs2(1:n-1)-cs3(1:n-1)/)



contains
function cumsum(x,ist,iend) result(cumx)
implicit none
real(wp), dimension(n), intent(in) :: x
integer, intent(in) :: ist,iend
real(wp), dimension(n) :: cumx
integer :: i,k

cumx(:) =0._wp
if (ist < iend) then
        cumx = (/ (sum(x(ist:i)), i=ist,iend) /)
elseif (ist > iend) then
        k=0
        do i=ist,iend,-1    
        k=k+1
        cumx(k) = sum(x(ist:i:-1))
        end do
end if

end function cumsum

function vert_adv(qs,qn) result(w)
implicit none
real(wp), dimension(n), intent(in) :: qs,qn
real(wp), dimension(n) :: csq,w

csq    = cumsum(qs-qn,n,1)
w(n:1:-1) = (/0._wp,csq(1:n-1)/)

end function vert_adv
end subroutine oce

include 'Rho_Vert_conv.f90'
include 'Med_Flux.f90'
include 'Ekman.f90'
include 'Init_Downsloping.f90'
include 'Entrainment.f90'