subroutine ocean_rad(t,inso,qsom,qarm,qm)
! Calculates solar radiation used for ocean new production as function of time


use parameters, only: wp,b_south,fdh,fdl,fde,fbr,nob
implicit none
real(wp), intent(in) :: t,inso(19,14)
real(wp), intent(out) :: qsom,qarm,qm(nob)
real(wp) :: qfn(90),qfs(90),iq(180),phi(180),qhnm,qmnm,qenm,qesm,qmsm
logical, dimension(180) :: iso,ims,ies,ien,imn,ihn,iar

call get_insolation(t,inso,qfn,qfs,phi)
qfs = qfs(size(qfs):1:-1)

iq  = (/qfs,qfn/)

iso = ( phi >= -b_south .AND. phi < -fdh )
ims = ( phi >= -fdh     .AND. phi < -fdl )
ies = ( phi >= -fdl     .AND. phi <  fde )
ien = ( phi >=  fde     .AND. phi <  fdl )
imn = ( phi >=  fdl     .AND. phi <  fdh )
ihn = ( phi >=  fdh     .AND. phi <  fbr )
iar = ( phi >=  fbr)

qsom = sum( pack(iq,iso)*cos(pack(phi,iso)) )/sum(cos(pack(phi,iso)))
qmsm = sum( pack(iq,ims)*cos(pack(phi,ims)) )/sum(cos(pack(phi,ims)))
qesm = sum( pack(iq,ies)*cos(pack(phi,ies)) )/sum(cos(pack(phi,ies)))
qenm = sum( pack(iq,ien)*cos(pack(phi,ien)) )/sum(cos(pack(phi,ien)))
qmnm = sum( pack(iq,imn)*cos(pack(phi,imn)) )/sum(cos(pack(phi,imn)))
qhnm = sum( pack(iq,ihn)*cos(pack(phi,ihn)) )/sum(cos(pack(phi,ihn)))
qarm = sum( pack(iq,iar)*cos(pack(phi,iar)) )/sum(cos(pack(phi,iar)))

qm = (/qhnm,qmnm,qenm,qesm,qmsm/)

end subroutine ocean_rad
