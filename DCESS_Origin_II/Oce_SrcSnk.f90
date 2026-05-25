subroutine oce_srcsnk(hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc,&
                      aoifa,aoifp,aoifar,aoifso,gaa,gap,gaar,gaso,&
                      omea,omep,omear,omeso,co2a,co2p,co2ar,co2so,&
                      co3a,co3p,co3ar,co3so,porga,porgp,porgar,porgso,&
                      pcala,pcalp,pcalar,pcalso,qsom,qarm,qm,&
                      srar,srso,srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp)
! Calculates source/sink terms for all ocean boxes
! Inputs:
! Ocean tracers:hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc
! Ice-free ocean area: aoifa,aoifp,aoifar,aoifso
! Vertical distribution of ocean area: gaa,gap,gaar,gaso
! Calcite saturation state: omea,omep,omear,omeso
! Ocean carbon species: co2,co3
! percentage organic carbon rem: porg
! percentage organic calcite dissolution: pcal
! sun radiation: qsom,qarm,qm
! Outputs:
! sr: internal source/sink vertical distributions

use parameters, only: wp,nto,n,lfar,lfso,lf,b_south,fdh,fdl,fde,fbr,sy,dm,rcp,rcd,rda,rca,rno,rcop,mino2,&
rpm,knut,kcar,q10,zmid,d,lmdcar,nob
implicit none
real(wp), intent(in) :: aoifar,aoifso,omear,omeso,&
                        co2ar,co2so,co3ar,co3so
real(wp), dimension(nob), intent(in) :: aoifa,aoifp,omea,omep,co2a,co2p,co3a,co3p
real(wp), dimension(nob,n), intent(in) :: gaa,gap,porga,porgp,pcala,pcalp
real(wp), dimension(n), intent(in) :: gaar,gaso,porgar,porgso,pcalar,pcalso
real(wp), dimension(nto,n), intent(in) :: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc
real(wp), intent(in) :: qsom,qarm,qm(nob)
real(wp), dimension(nto,n), intent(out) :: srar,srso,srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp
integer :: i

call orgflx(arc,aoifar,lfar,qarm,gaar,omear,co2ar,co3ar,porgar,pcalar,srar)

call orgflx(soc,aoifso,lfso,qsom,gaso,omeso,co2so,co3so,porgso,pcalso,srso)

i=1
call orgflx(hna,aoifa(i),lf(i),qm(i),gaa(i,:),omea(i),co2a(i),co3a(i),porga(i,:),pcala(i,:),srhna)
i=2
call orgflx(mna,aoifa(i),lf(i),qm(i),gaa(i,:),omea(i),co2a(i),co3a(i),porga(i,:),pcala(i,:),srmna)
i=3
call orgflx(ena,aoifa(i),lf(i),qm(i),gaa(i,:),omea(i),co2a(i),co3a(i),porga(i,:),pcala(i,:),srena)
i=4
call orgflx(esa,aoifa(i),lf(i),qm(i),gaa(i,:),omea(i),co2a(i),co3a(i),porga(i,:),pcala(i,:),sresa)
i=5
call orgflx(msa,aoifa(i),lf(i),qm(i),gaa(i,:),omea(i),co2a(i),co3a(i),porga(i,:),pcala(i,:),srmsa)

i=1
call orgflx(hnp,aoifp(i),lf(i),qm(i),gap(i,:),omep(i),co2p(i),co3p(i),porgp(i,:),pcalp(i,:),srhnp)
i=2
call orgflx(mnp,aoifp(i),lf(i),qm(i),gap(i,:),omep(i),co2p(i),co3p(i),porgp(i,:),pcalp(i,:),srmnp)
i=3
call orgflx(enp,aoifp(i),lf(i),qm(i),gap(i,:),omep(i),co2p(i),co3p(i),porgp(i,:),pcalp(i,:),srenp)
i=4
call orgflx(esp,aoifp(i),lf(i),qm(i),gap(i,:),omep(i),co2p(i),co3p(i),porgp(i,:),pcalp(i,:),sresp)
i=5
call orgflx(msp,aoifp(i),lf(i),qm(i),gap(i,:),omep(i),co2p(i),co3p(i),porgp(i,:),pcalp(i,:),srmsp)



end subroutine oce_srcsnk
include 'Org_Flx.f90'
