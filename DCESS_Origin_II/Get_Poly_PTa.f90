subroutine get_poly_pta(ta_n,ta_mn,ta_en,ta_es,ta_ms,ta_s,pta_nh,pta_sh)
! Calculates factors to build atmospheric temperature meridional distribution
! Inputs: atmospheric temperature
! Outputs: factors pta_nh,pta_sh

use parameters, only: wp,fdh,fdl

real(wp), intent(in) :: ta_n,ta_mn,ta_en,ta_es,ta_ms,ta_s
real(wp), dimension(3), intent(out) :: pta_nh,pta_sh
real(wp) :: c0,c12,c13,c22,c23,c32,c33,auxden

c0 = 1._wp
c12 = (1._wp/2_wp)*(sin(fdl)**2-1._wp)
c13 = (1._wp/4._wp)*(3._wp/2._wp - 5._wp*sin(fdl)**2 + (7._wp/2._wp)*sin(fdl)**4)
c22 = (1._wp/2._wp)*(sin(fdh)**3 - sin(fdl)**3 - sin(fdh) + sin(fdl))/(sin(fdh)-sin(fdl))
c23 = (1._wp/4._wp)*( (3._wp/2._wp)*(sin(fdh)-sin(fdl))-5._wp*(sin(fdh)**3-sin(fdl)**3)+&
        (7._wp/2._wp)*(sin(fdh)**5-sin(fdl)**5) )/(sin(fdh)-sin(fdl))
c32 = (1._wp/2._wp)*(sin(fdh)-sin(fdh)**3)/(1._wp-sin(fdh))
c33 = -(1._wp/4._wp)*( (3._wp/2._wp)*sin(fdh)-5._wp*sin(fdh)**3+(7._wp/2._wp)*sin(fdh)**5 )/(1._wp-sin(fdh))
auxden=1._wp-( (c12-c32)*(c13-c23*c0)+c13*(c22*c0-c12) )/(c33*(c22*c0-c12))

pta_nh(3)=1._wp/auxden*1._wp/c33*( ta_n-ta_en+ (c12-c32)*(ta_mn*c0-ta_en)/(c22*c0-c12) )
pta_nh(2)=1._wp/(c22-c12)*( ta_mn-ta_en+(c13-c23)*pta_nh(3) )
pta_nh(1)=1._wp/c0*( ta_en-c12*pta_nh(2)-c13*pta_nh(3) )

pta_sh(3)=1._wp/auxden*1._wp/c33*( ta_s-ta_es+ (c12-c32)*(ta_ms*c0-ta_es)/(c22*c0-c12) )
pta_sh(2)=1._wp/(c22-c12)*( ta_ms-ta_es+(c13-c23)*pta_sh(3) )
pta_sh(1)=1._wp/c0*( ta_es-c12*pta_sh(2)-c13*pta_sh(3) )


end subroutine get_poly_pta