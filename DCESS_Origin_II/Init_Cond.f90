subroutine init_cond(at,atam,lbn,lbs,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,soc,arc,slso,&
                     pcala,pcalp,pcalar,pcalso,porga,porgp,porgar,porgso,&
                     dwpcala,dwpcalssa,dwporga,dwporgssa,fima,wseda,dwpcalp,dwpcalssp,dwporgp,dwporgssp,fimp,wsedp,&
                     dwpcalar,dwpcalssar,dwporgar,dwporgsar,fimar,wsedar,dwpcalso,dwpcalssso,dwporgso,dwporgssso,fimso,wsedso)

! Load initial conditions to run the model

use parameters, only: wp,nta,nto,nlb,n,nab,nob

implicit none
real(wp), dimension(nto,n), intent(out) :: hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,soc,arc
real(wp), intent(out) :: slso(nto)
real(wp), dimension(nta,nab), intent(out) :: at,atam
real(wp), dimension(nlb), intent(out) :: lbn,lbs
real(wp), dimension(nob,n), intent(out) :: pcala,pcalp,porga,porgp
real(wp), dimension(n), intent(out) :: pcalar,pcalso,porgar,porgso
real(wp), dimension(nob,n), intent(out) :: dwpcala,dwpcalssa,dwporga,dwporgssa,fima,wseda,&
                                         dwpcalp,dwpcalssp,dwporgp,dwporgssp,fimp,wsedp
real(wp), dimension(n), intent(out) :: dwpcalar,dwpcalssar,dwporgar,dwporgsar,fimar,wsedar,&
                                       dwpcalso,dwpcalssso,dwporgso,dwporgssso,fimso,wsedso


integer :: i

open(unit=67,file='Input_Data/AT.txt',status='old')
read(67,*) (at(i,:), i=1,nta); close(67)

open(unit=68,file='Input_Data/ATam.txt',status='old')
read(68,*) (atam(i,:), i=1,nta); close(68)

open(unit=69,file='Input_Data/LBN.txt',status='old')
read(69,*) (lbn(i), i=1,nlb); close(69)

open(unit=70,file='Input_Data/LBS.txt',status='old')
read(70,*) (lbs(i), i=1,nlb); close(70)

open(unit=71,file='Input_Data/HLNA.txt',status='old')
read(71,*) (hlna(i,:), i=1,nto); close(71)

open(unit=72,file='Input_Data/MLNA.txt',status='old')
read(72,*) (mlna(i,:), i=1,nto); close(72)

open(unit=73,file='Input_Data/LLNA.txt',status='old')
read(73,*) (llna(i,:), i=1,nto); close(73)

open(unit=74,file='Input_Data/LLSA.txt',status='old')
read(74,*) (llsa(i,:), i=1,nto); close(74)

open(unit=75,file='Input_Data/MLSA.txt',status='old')
read(75,*) (mlsa(i,:), i=1,nto); close(75)

open(unit=76,file='Input_Data/HLNP.txt',status='old')
read(76,*) (hlnp(i,:), i=1,nto); close(76)

open(unit=77,file='Input_Data/MLNP.txt',status='old')
read(77,*) (mlnp(i,:), i=1,nto); close(77)

open(unit=78,file='Input_Data/LLNP.txt',status='old')
read(78,*) (llnp(i,:), i=1,nto); close(78)

open(unit=79,file='Input_Data/LLSP.txt',status='old')
read(79,*) (llsp(i,:), i=1,nto); close(79)

open(unit=80,file='Input_Data/MLSP.txt',status='old')
read(80,*) (mlsp(i,:), i=1,nto); close(80)

open(unit=81,file='Input_Data/SOc.txt',status='old')
read(81,*) (soc(i,:), i=1,nto); close(81)

open(unit=82,file='Input_Data/Arc.txt',status='old')
read(82,*) (arc(i,:), i=1,nto); close(82)

open(unit=83,file='Input_Data/SlSO.txt',status='old')
read(83,*) (slso(i), i=1,nto); close(83)

! Sediment
open(unit=84,file='Input_Data/pCalA.txt',status='old')
read(84,*) (pcala(i,:), i=1,nob); close(84)

open(unit=85,file='Input_Data/pCalP.txt',status='old')
read(85,*) (pcalp(i,:), i=1,nob); close(85)

open(unit=86,file='Input_Data/pCalAr.txt',status='old')
read(86,*) (pcalar(i), i=1,n); close(86)

open(unit=87,file='Input_Data/pCalSO.txt',status='old')
read(87,*) (pcalso(i), i=1,n); close(87)

open(unit=88,file='Input_Data/pOrgA.txt',status='old')
read(88,*) (porga(i,:), i=1,nob); close(88)

open(unit=89,file='Input_Data/pOrgP.txt',status='old')
read(89,*) (porgp(i,:), i=1,nob); close(89)

open(unit=90,file='Input_Data/pOrgAr.txt',status='old')
read(90,*) (porgar(i), i=1,n); close(90)

open(unit=91,file='Input_Data/pOrgSO.txt',status='old')
read(91,*) (porgso(i), i=1,n); close(91)

open(unit=73,file='Input_Data/dwpCalA.txt',status='old')
read(73,*) (dwpcala(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpCalssA.txt',status='old')
read(73,*) (dwpcalssa(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpOrgA.txt',status='old')
read(73,*) (dwporga(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpOrgssA.txt',status='old')
read(73,*) (dwporgssa(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/fimA.txt',status='old')
read(73,*) (fima(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/wsedA.txt',status='old')
read(73,*) (wseda(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpCalAr.txt',status='old')
read(73,*) (dwpcalar(i), i=1,n); close(73)

open(unit=64,file='Input_Data/dwpCalssAr.txt',status='old')
read(64,*) (dwpcalssar(i), i=1,n); close(64)

open(unit=73,file='Input_Data/dwpOrgAr.txt',status='old')
read(73,*) (dwporgar(i), i=1,n); close(73)

open(unit=73,file='Input_Data/dwpOrgssAr.txt',status='old')
read(73,*) (dwporgsar(i), i=1,n); close(73)

open(unit=73,file='Input_Data/fimAr.txt',status='old')
read(73,*) (fimar(i), i=1,n); close(73)

open(unit=73,file='Input_Data/wsedAr.txt',status='old')
read(73,*) (wsedar(i), i=1,n); close(73)

open(unit=73,file='Input_Data/dwpCalP.txt',status='old')
read(73,*) (dwpcalp(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpCalssP.txt',status='old')
read(73,*) (dwpcalssp(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpOrgP.txt',status='old')
read(73,*) (dwporgp(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpOrgssP.txt',status='old')
read(73,*) (dwporgssp(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/fimP.txt',status='old')
read(73,*) (fimp(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/wsedP.txt',status='old')
read(73,*) (wsedp(i,:), i=1,nob); close(73)

open(unit=73,file='Input_Data/dwpCalSO.txt',status='old')
read(73,*) (dwpcalso(i), i=1,n); close(73)

open(unit=64,file='Input_Data/dwpCalssSO.txt',status='old')
read(64,*) (dwpcalssso(i), i=1,n); close(64)

open(unit=73,file='Input_Data/dwpOrgSO.txt',status='old')
read(73,*) (dwporgso(i), i=1,n); close(73)

open(unit=73,file='Input_Data/dwpOrgssSO.txt',status='old')
read(73,*) (dwporgssso(i), i=1,n); close(73)

open(unit=73,file='Input_Data/fimSO.txt',status='old')
read(73,*) (fimso(i), i=1,n); close(73)

open(unit=73,file='Input_Data/wsedSO.txt',status='old')
read(73,*) (wsedso(i), i=1,n); close(73)
end subroutine init_cond