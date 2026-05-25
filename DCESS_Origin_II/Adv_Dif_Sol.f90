subroutine adv_dif_sol(psii,psin,psis,kv,khn,khs,w,qn,qs,src,ao,ga,kpsi)
! Calculates advection-diffusion equation + internal source/sinks for each ocean box

! Inputs:
! psii: matrix containing ocean tracers of box "i"
! psin: matrix containing ocean tracers of box at north of psii
! psis: matrix containing ocean tracers of box at south of psii
! kv: vertical diffusion
! khn: horizontal diffusion at north boundary of psii
! khs: horizontal diffusion at south boundary of psii
! w: vertical advection
! qn: horizontal advection at north boundary of psii
! qs: horizontal advection at south boundary of psii
! src: internal/external source/sink terms of box "i"
! ao: Ocean area
! ga: vertical profile of area fraction
! Outputs:
! kpsi: time change rate of all tracers in psii ([]/s)

use parameters, only: wp,nto,n,dm,d

implicit none
real(wp), dimension(nto,n), intent(in) :: psii,psin,psis,src
real(wp), dimension(n), intent(in) :: khn,khs,w,qn,qs,ga
real(wp), dimension(n-1), intent(in) :: kv
real(wp), intent(in) :: ao
real(wp), dimension(nto,n), intent(out) :: kpsi
real(wp), dimension(nto,n) :: v_diff,h_diff,v_adv,h_adv,vol
real(wp), dimension(nto) :: zeronto

zeronto = 0._wp
v_diff = reshape( (/ -kv(1)*ao/(3._wp*d)*(8._wp*psii(:,1) - 9._wp*psii(:,2) + psii(:,3) ) ,&
                      spread(kv(2:n-1),1,nto)*ao/d*(psii(:,3:n) - psii(:,2:n-1))          , zeronto  /), (/nto,n/) ) - &
         reshape( (/ zeronto , -kv(1)*ao/(3._wp*d)*(8._wp*psii(:,1) - 9._wp*psii(:,2) + psii(:,3) ),&
                      spread(kv(2:n-1),1,nto)*ao/d*(psii(:,3:n) - psii(:,2:n-1)) /), (/nto,n/) )

h_diff = - spread(khn,1,nto)*(psii-psin) + spread(khs,1,nto)*(psis-psii)

v_adv = reshape(  (/ spread(w(1:n-1),1,nto)*(psii(:,1:n-1) + psii(:,2:n))*0.5_wp,zeronto /),(/nto,n/)   ) - &
        reshape(  (/ zeronto, spread(w(1:n-1),1,nto)*(psii(:,1:n-1) + psii(:,2:n))*0.5_wp /), (/nto,n/) )

h_adv = - spread( merge(1._wp,0._wp,qn .GT. 0._wp)*qn,1,nto )*psii - spread( merge(1._wp,0._wp,qn .LT. 0._wp)*qn,1,nto )*psin + &
          spread( merge(1._wp,0._wp,qs .GT. 0._wp)*qs,1,nto )*psis + spread( merge(1._wp,0._wp,qs .LT. 0._wp)*qs,1,nto )*psii

vol = spread( ao*(/dm,d*ga(2:n)/),1,nto )

kpsi = 1._wp/vol*( v_diff + h_diff + v_adv + h_adv + src )


end subroutine adv_dif_sol