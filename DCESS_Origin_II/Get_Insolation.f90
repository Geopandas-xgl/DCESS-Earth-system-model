subroutine get_insolation(t,inso,qf_n,qf_s,phii)
! Get meridional distribution of insolation from observed data inso

use parameters, only: wp,sy,dr,pi
use :: lib_array, only: interp1d

implicit none
real(wp), intent(in) :: t,inso(19,14)
real(wp), intent(out) :: qf_n(90),qf_s(90),phii(180)
real(wp) :: mth,upm,dupm,lom,dlom,q(19,14),qf10(19),phi10(19),qf(180)
real(wp), dimension(90) :: phi
integer :: i

phi = (/(real(i,wp)+0.5_wp, i=0,89,1 )/)*dr

! Set decimal month from input time (t)
mth  = mod((t/sy)/(1._wp/12._wp),12._wp)+1._wp

! Upper and Lower month, used for interpolation
upm  = real(ceiling(mth),wp)
dupm = 1._wp-(upm-mth)
if (upm>12) then
        upm  = 1._wp
        dupm = 1._wp - (upm-mth+12._wp)
end if
lom  = real(floor(mth),wp)
dlom = 1._wp-(mth-lom)

q = inso
qf10   = .4843_wp*( q(:,2+int(upm))*dupm + q(:,2+int(lom))*dlom )/(dupm+dlom) ! (19,1)
phi10  = q(:,2)*dr

! Interpolate Qf to the one degree latitude resolution, (f) 
! interp1d(x,y,xi)
phii = (/-pi/2._wp+phi,phi/)
qf = interp1d(phi10,qf10,phii)
where (isnan(qf)) qf=0._wp

qf_s = qf(90:1:-1) 
qf_n = qf(91:180:1)

! Remove unphysical negative interpolation values of qf at the ends, occuring when 
! insolation is zero for more than one ten-degree latitude band (winter extreme)
qf_n  = qf_n*merge(1._wp,0._wp,qf_n>0._wp) + 0._wp*merge(1._wp,0._wp,qf_n<=0._wp)
qf_s  = qf_s*merge(1._wp,0._wp,qf_s>0._wp) + 0._wp*merge(1._wp,0._wp,qf_s<=0._wp)


end subroutine get_insolation
