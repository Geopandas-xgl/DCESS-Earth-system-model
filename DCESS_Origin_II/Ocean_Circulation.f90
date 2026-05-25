subroutine ocean_circulation(at,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,soc,arc,slso,&
                             inso,oc_pres,shl_dat,&
                             q_atl,q_pac,ekman)
use parameters, only: wp,nto,n,nob,nta,nab,nz_all
implicit none
real(wp), intent(in), dimension(nto,n) :: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,soc,arc
real(wp), intent(in) , dimension(nta,nab) :: at
real(wp), intent(in) :: inso(19,14),oc_pres(n,7),shl_dat(nz_all,3)
real(wp), intent(in) :: slso(nto)

real(wp), intent(out) :: q_atl(n,6),q_pac(n,5),ekman(n,2)

real(wp) :: to_a(nob),to_p(nob),to_ar,to_so
real(wp), dimension(4) :: hsen,hlat,htot
real(wp) :: fwfbr,fwfdhn,fwfdln,fwfdls,fwfdhs,hfde
real(wp), dimension(nab) :: totsw,totlw
real(wp) :: asear,aseat(nob),asepa(nob),aseso,fsnn,fsns
real(wp) :: dfina,dfinp,dfis,aoifa(nob),aoifp(nob),aoifar,aoifso
real(wp), dimension(3) :: ptanh,ptash

real(wp), dimension(n-1):: kvhna,kvmna,kvena,kvesa,kvmsa,&
                           kvhnp,kvmnp,kvenp,kvesp,kvmsp,kvar,kvso
real(wp), dimension(n) ::  khfbra,khfdhna,khfdlna,khfdea,khfdlsa,khfdhsa,&
                           khfdhnp,khfdlnp,khfdep,khfdlsp,khfdhsp,khmsap,&
                           qfbra,qfdhna,qfdhnp,qfdlna,qfdlnp,qfdea,qfdep,&
                           qfdlsa,qfdlsp,qfdhsa,qfdhsp,ievfsa,ievfsp,ievfbr,&
                           whna,wmna,wena,wesa,wmsa,whnp,wmnp,wenp,wesp,wmsp,wso,war
real(wp) :: fdsup,feso(n),rfds_v(n),qshl
integer :: ent_lim(2),idlsa,idlsp,idlar

to_a  = (/hna(1,1),mna(1,1),ena(1,1),esa(1,1),msa(1,1)/)
to_p  = (/hnp(1,1),mnp(1,1),enp(1,1),esp(1,1),msp(1,1)/)
to_ar = arc(1,1)
to_so = soc(1,1)

call atm(0._wp,at,to_a,to_p,to_ar,to_so,inso,&
         hsen,hlat,htot,hfde,fwfbr,fwfdhn,fwfdln,fwfdls,fwfdhs,&
         totsw,totlw,asear,aseat,asepa,aseso,&
         dfina,dfinp,dfis,aoifa,aoifp,aoifar,aoifso,ptanh,ptash,fsnn,fsns)

call oce(arc,soc,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,slso,&
         fwfbr,fwfdhn,fwfdln,fwfdls,fwfdhs,oc_pres,shl_dat,&
         kvhna,kvmna,kvena,kvesa,kvmsa,kvhnp,kvmnp,kvenp,kvesp,kvmsp,kvar,kvso,&
         khfbra,khfdhna,khfdlna,khfdea,khfdlsa,khfdhsa,khfdhnp,khfdlnp,khfdep,khfdlsp,khfdhsp,khmsap,&
         qfbra,qfdhna,qfdhnp,qfdlna,qfdlnp,qfdea,qfdep,qfdlsa,qfdlsp,qfdhsa,qfdhsp,&
         ievfsa,idlsa,ievfsp,idlsp,ievfbr,idlar,fdsup,feso,ent_lim,rfds_v,&
         whna,wmna,wena,wesa,wmsa,whnp,wmnp,wenp,wesp,wmsp,wso,war,qshl)

q_atl(:,1) = qfdhsa
q_atl(:,2) = qfdlsa
q_atl(:,3) = qfdea
q_atl(:,4) = qfdlna
q_atl(:,5) = qfdhna
q_atl(:,6) = qfbra

q_pac(:,1) = qfdhsp
q_pac(:,2) = qfdlsp
q_pac(:,3) = qfdep
q_pac(:,4) = qfdlnp
q_pac(:,5) = qfdhnp

ekman(:,1) = ievfsa
ekman(:,2) = ievfsp
end subroutine ocean_circulation