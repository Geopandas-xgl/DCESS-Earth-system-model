subroutine river_input(at,srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp,srar,srso,&
                       rorgat,rorgpa,rorgar,rorgso,rcarat,rcarpa,rcarar,rcarso)
! Call river input for all ocean boxes

use parameters, only: wp,nto,n,nta,nab,nob

real(wp), dimension(nto,n), intent(inout) :: srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp,srar,srso
real(wp), intent(in) :: at(nta,nab),rorgat(nob),rorgpa(nob),rorgar,rorgso,rcarat(nob),rcarpa(nob),rcarar,rcarso

call riv_in(srhna,rorgat(1),rcarat(1),at(5:7,1))
call riv_in(srmna,rorgat(2),rcarat(2),at(5:7,2))
call riv_in(srena,rorgat(3),rcarat(3),at(5:7,3))
call riv_in(sresa,rorgat(4),rcarat(4),at(5:7,4))
call riv_in(srmsa,rorgat(5),rcarat(5),at(5:7,5))

call riv_in(srhnp,rorgpa(1),rcarpa(1),at(5:7,1))
call riv_in(srmnp,rorgpa(2),rcarpa(2),at(5:7,2))
call riv_in(srenp,rorgpa(3),rcarpa(3),at(5:7,3))
call riv_in(sresp,rorgpa(4),rcarpa(4),at(5:7,4))
call riv_in(srmsp,rorgpa(5),rcarpa(5),at(5:7,5))

call riv_in(srar,rorgar,rcarar,at(5:7,1))
call riv_in(srso,rorgso,rcarso,at(5:7,6))

end subroutine river_input
