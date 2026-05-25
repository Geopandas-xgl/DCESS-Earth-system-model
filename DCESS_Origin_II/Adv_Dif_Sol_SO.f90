subroutine adv_dif_sol_so(psii,psinp,psina,psis,kv,khnp,khna,khs,w,qnp,qna,qs,src,ao,ga,kpsi)
! Calculates advection-diffusion equation + internal source/sinks for a tracer Psi
!
! Inputs:
! psii: matrix containing ocean tracers of box "i"
! psinp/a: north sector tracers Atlantic/Pacific (msp/msa)
! kv: vertical diffusion
! khnp/a: horizontal diffusion Atlantic/Pacific (msp/msa)
! khs: horizontal diffusion south (zero for Southern Ocean)
! w: vertical advection
! qnp/a: horizontal advection Atlantic/Pacific (msp/msa)
! qs: horizontal advection south (zero for Southern Ocean)
! src: internal/external source/sink terms of box "i"
! ao: Ocean area
! ga: vertical profile of area fraction
! Outputs:
! kpsi: time change rate of all tracers in psi ([]/s)

use parameters, only: wp,nto,n,dm,d

implicit none
real(wp), dimension(nto,n), intent(in) :: psii,psinp,psina,psis,src
real(wp), dimension(n), intent(in) :: khnp,khna,khs,w,qnp,qna,qs,ga
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

h_diff = - spread(khna,1,nto)*(psii-psina) + spread(khs,1,nto)*(psis-psii) + &
         - spread(khnp,1,nto)*(psii-psinp) + spread(khs,1,nto)*(psis-psii)

v_adv = reshape(  (/ spread(w(1:n-1),1,nto)*(psii(:,1:n-1) + psii(:,2:n))*0.5_wp,zeronto /),(/nto,n/)   ) - &
        reshape(  (/ zeronto, spread(w(1:n-1),1,nto)*(psii(:,1:n-1) + psii(:,2:n))*0.5_wp /), (/nto,n/) )
        
h_adv = - spread( merge(1._wp,0._wp,qna .GT. 0._wp)*qna,1,nto )*psii  &
        - spread( merge(1._wp,0._wp,qna .LT. 0._wp)*qna,1,nto )*psina & 
        + spread( merge(1._wp,0._wp,qs  .GT. 0._wp)*qs ,1,nto )*psis  &
        + spread( merge(1._wp,0._wp,qs  .LT. 0._wp)*qs ,1,nto )*psii  & 
        - spread( merge(1._wp,0._wp,qnp .GT. 0._wp)*qnp,1,nto )*psii  &
        - spread( merge(1._wp,0._wp,qnp .LT. 0._wp)*qnp,1,nto )*psinp &
        + spread( merge(1._wp,0._wp,qs  .GT. 0._wp)*qs ,1,nto )*psis  &
        + spread( merge(1._wp,0._wp,qs  .LT. 0._wp)*qs ,1,nto )*psii
        
vol = spread( ao*(/dm,d*ga(2:n)/),1,nto )
kpsi = 1._wp/vol*( v_diff + h_diff + v_adv + h_adv + src )

end subroutine adv_dif_sol_so