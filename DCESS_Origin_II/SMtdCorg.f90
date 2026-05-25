subroutine smtdcorg(psii,nnpp,dwpccalss,dwpoorgss,&
                    dwpccal,dwpoorg,fim,wsed,pccal,poorg,ko,h)

! "Complete" time dependent sediment model in a 10 cm thick, bioturbated
! sediment layer, driven by water column concentrations from the model
! and carbonate flxes, simpliest version with dissolution not depending
! on organic carbon remineralization, approximate carbonate solution
! (ignoring [CO2]),all 5 sediment layers, Carbonate and Non-carbonate, fractions only


use corg_data
use parameters, only: wp,d,dm,sy,n,nto,rcp,mino2,lmdcar,kcalcite,mmca,rhom,ncf,caf


implicit none
real(wp), dimension(nto,n), intent(in) :: psii
real(wp), intent(in) :: nnpp,ko,h
real(wp), dimension(n), intent(in) :: dwpccalss,dwpoorgss
real(wp), dimension(n), intent(inout) :: dwpccal,dwpoorg,fim,wsed
real(wp), dimension(n), intent(out) :: pccal, poorg

real(wp) ::llhh(nto,n),ddbo,rrmmorgc
real(wp), dimension(n) :: ttllhh,ccoo3llhh,ccoo3llhhs
real(wp) :: kc,nnccffo,mmoc,rhoc,compc,jo2o,rr,difzb(7),dtc,fzo,rpm,ttr,p1,p2,rrmmcar,&
            nncc,ffcal,mmy,ddo2,ddccoo3,pperccaccoo3o,alf,fimax,fimeano,fiboto,pperoorgo,&
            rhmeano,pperccaccoo3,pperoorg,rhmean,fimeanb,fimean,fibot,&
            ttcaldis,ttcaldis1,ttcaldis2,ttcaldis3,ttcaldis4,ttcaldis5,ttcaldis6,ttcaldis7,s3,&
            alf2,alf3,alf4,alf5,alf6,alf7,co,x(14),zo,&
            c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,&
            c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,cc(14,14),b(14),ipvt(14),&
            kcal,ttremin,dpp,doorg,rres1,rres2,rres3,rres4,rres5,rres6,rres7
real(wp) :: zbrent
integer :: i,k,pt
real(wp), external :: corg1a7w_eq,corg2a7w_eq,corg3a7w_eq,corg4a7w_eq,corg5a7w_eq,&
                      corg6a7w_eq,corg7a7w_eq
real(wp) :: pon(n),poc(n),pca(n),rp
real(wp), dimension(n) :: co3_p,co3s_p,co2_p,hco3_p,hp_p,llhhccoo3,llhhccoo3s
real(wp), external :: rain_ratio
real(wp) :: pco2,k0,co2,co3,hco3,ome,nnccff,ccaaff

llhh = psii
nnccff = ncf
ccaaff = caf

call carsys_pres(llhh,co3_p,co3s_p,co2_p,hco3_p,hp_p)
llhhccoo3  = co3_p
llhhccoo3s = co3s_p

call carsys_sur(llhh(1,1),llhh(2,1),llhh(5,1),llhh(8,1),pco2,k0,co2,co3,hco3,ome)
rp = rain_ratio(ome,llhh(1,1))

where (llhh(9,:) <= mino2) llhh(9,:)=mino2



ttllhh     = llhh(1,:)+273.16_wp
ccoo3llhh  = llhhccoo3*1.e-6_wp
ccoo3llhhs = llhhccoo3s*1.e-6_wp
kc         = kcalcite/(3600._wp*24._wp)     ! Calcite dissolution rate constant [s-1], st. val, 0.0003
nnccffo    = nnccff/1000._wp/sy             ! Non-carbonate rain [g/cm^2/s]

! Carbon burial parameters
mmoc  = 12._wp
rhoc  = 1.1_wp                      ! organic matter density, g per cm3
compc = 2.7_wp                      ! g organic matter per g C 
ddbo  = 0.0232_wp*18.774_wp/sy      !  Archer et al 2002, surface value ! 18.774 factor is 31.5^^0.85
jo2o  = 1.e-9_wp                               
rr    = 1.4_wp                      ! O to C ratio


zb    = (/0._wp,0.2_wp,0.5_wp,1._wp,1.8_wp,3.2_wp,6._wp,10._wp/)  ! Chosen layer boundaries in sediment layer, cm 
difzb =(/zb(2)-zb(1),zb(3)-zb(2),zb(4)-zb(3),zb(5)-zb(4),zb(6)-zb(5),zb(7)-zb(6),zb(8)-zb(7)/)

dtc  = ko*h*sy              ! time step, must be consistent with time step in Thilda!
fzo  = 1._wp                ! porosity at the sediment surface
! print*, rp
! Surface flux of CaCO3 and organic carbon 
rrmmcar = -nnpp*rcp*rp/(1.e4_wp)  ! [mol/cm2/s]
rrmmorgc= -nnpp*rcp/(1.e4_wp)
call parmat_depth(llhh,pon,poc,pca)

! Loop layers
do i =1,n

! Non carbonate flux, rapidly decreasing with depth to
! mimic high levels in coastal areas
nncc   = ccaaff*nnccffo*exp((1._wp-real(i,wp))/2._wp)+ nnccffo    ! [g/cm2/s]

! Flux of calcite, exponential decay
ffcal  = rrmmcar*exp(-lmdcar*(real(i,wp)*d-dm))      ! [mol/cm2/s]

! Organic carbon flux, exponential decay
ffcc = rrmmorgc*poc(i)    !Org C flux in mole C cm-2 s-1  

iinioo2 = llhh(9,i)*1.e-6_wp
ddb     = ddbo/18.774_wp*((ffcc*1.e6_wp*sy)**0.85)*iinioo2/(iinioo2+20._wp*1.e-9_wp)
jo2     = jo2o*ddb/ddbo
beta    = 0.1_wp*(ffcc*1.e6_wp*sy/31.5_wp)**(-0.3)

! Temperature dependent diffusion of O2 and CO3
mmy     = 2.31_wp/(1._wp + 0.036_wp*(ttllhh(i)-273.16_wp) + 1.85e-4_wp*(ttllhh(i)-273.16_wp)**2)
ddo2    = 1.25e-5_wp*ttllhh(i)/278.16_wp*1.95_wp/mmy
ddccoo3 = 7.9e-6_wp*(ttllhh(i)/291.16_wp)*(1.3525_wp/mmy)

pperccaccoo3o=dwpccalss(i)
alf     = 0.25_wp*pperccaccoo3o+3._wp*(1._wp-pperccaccoo3o)
fimax   = 1._wp-(0.483_wp+0.45_wp*pperccaccoo3o)/2.5_wp
fimeano = fimax-alf*(1._wp-fimax)*(exp(-zb(8)/alf)-1._wp)/zb(8)  ! mean porosity in 10 cm bioturbated layer
fiboto  = fimax-alf*(fzo-fimax)*(exp(-zb(8)/alf)-exp(-zb(7)/alf))/(zb(8)-zb(7))         ! mean porosity of bottom layer

pperoorgo = dwpoorgss(i)
rhmeano =(1._wp-pperoorgo)*rhom+pperoorgo*rhoc

pperccaccoo3 = dwpccal(i)
pperoorg     = dwpoorg(i)
rhmean       = (1._wp-pperoorg)*rhom+pperoorg*rhoc
fimeanb      = fim(i)
w            = wsed(i)

if (ccoo3llhh(i) >= ccoo3llhhs(i)) then
    ! Case: super saturated, no calcite dissolution:  
    alf    = 0.25_wp*pperccaccoo3+3._wp*(1._wp-pperccaccoo3)
    fimax  = 1._wp-(0.483_wp+0.45_wp*pperccaccoo3)/2.5_wp
    fimean = fimax-alf*(fzo-fimax)*(exp(-zb(8)/alf)-1._wp)/(zb(8))
    fibot  = fimax+(fzo-fimax)*exp(-zb(8)/alf)
    ttcaldis  = 0._wp
    ttcaldis1 = 0._wp
    ttcaldis2 = 0._wp
    ttcaldis3 = 0._wp
    ttcaldis4 = 0._wp
    ttcaldis5 = 0._wp
    ttcaldis6 = 0._wp
    ttcaldis7 = 0._wp
    s3        = 0._wp

    ! Calculate mean sediment formation factor F for each sediment model layer
    do k=1,7
    ff(k) = (zb(k+1)-zb(k))/((fimax**3)*(zb(k+1)-zb(k))-3._wp*fimax**2*(fzo-fimax)*alf*&
            (exp(-zb(k+1)/alf)-exp(-zb(k)/alf))-(3._wp/2._wp)*fimax*(fzo-fimax)**2*alf*&
            (exp(-2._wp*zb(k+1)/alf)-exp(-2._wp*zb(k)/alf))-(1._wp/3._wp)*(fzo-fimax)**3*alf*&
            (exp(-3._wp*zb(k+1)/alf)-exp(-3._wp*zb(k)/alf)))
    fi(k) = fimax-alf*(fzo-fimax)*(exp(-zb(k+1)/alf)-exp(-zb(k)/alf))/(zb(k+1)-zb(k))
    end do

    fimean = sum(fi*difzb)/zb(8)

else
    ! Case: under saturated, calcite dissolution:
    alf    = 0.25_wp*pperccaccoo3+3._wp*(1._wp-pperccaccoo3)
    fimax  = 1._wp-(0.483_wp+0.45_wp*pperccaccoo3)/2.5_wp
    fimean = fimax-alf*(fzo-fimax)*(exp(-zb(8)/alf)-1._wp)/(zb(8))
    fibot  = fimax+(fzo-fimax)*exp(-zb(8)/alf)

    do k=1,7
    ff(k) = (zb(k+1)-zb(k))/((fimax**3)*(zb(k+1)-zb(k))-3._wp*fimax**2*(fzo-fimax)*alf*&
            (exp(-zb(k+1)/alf)-exp(-zb(k)/alf))-(3._wp/2._wp)*fimax*(fzo-fimax)**2*alf*&
            (exp(-2._wp*zb(k+1)/alf)-exp(-2._wp*zb(k)/alf))-(1._wp/3._wp)*(fzo-fimax)**3*alf*&
            (exp(-3._wp*zb(k+1)/alf)-exp(-3._wp*zb(k)/alf)))
    fi(k) = fimax-alf*(fzo-fimax)*(exp(-zb(k+1)/alf)-exp(-zb(k)/alf))/(zb(k+1)-zb(k))
    end do

    fimean = sum(fi*difzb)/zb(8)

    ! Apply analytical solutions for CO3 concentration,
    s3 = ((kc*pperccaccoo3*rhmean*(1._wp-fi(1))*ff(1))/(ccoo3llhhs(i)*mmca*ddccoo3))**(0.5)

    alf2 = ((1._wp-fi(2))*ff(2)/((1._wp-fi(1))*ff(1)))**0.5
    alf3 = ((1._wp-fi(3))*ff(3)/((1._wp-fi(1))*ff(1)))**0.5
    alf4 = ((1._wp-fi(4))*ff(4)/((1._wp-fi(1))*ff(1)))**0.5
    alf5 = ((1._wp-fi(5))*ff(5)/((1._wp-fi(1))*ff(1)))**0.5
    alf6 = ((1._wp-fi(6))*ff(6)/((1._wp-fi(1))*ff(1)))**0.5
    alf7 = ((1._wp-fi(7))*ff(7)/((1._wp-fi(1))*ff(1)))**0.5

    co=(ccoo3llhh(i)- ccoo3llhhs(i))

    c1=exp(s3*zb(2))
    c2=exp(-s3*zb(2))
    c3=exp(alf2*s3*zb(2))
    c4=exp(-alf2*s3*zb(2))
    c5=ff(1)/ff(2)*alf2*c3
    c6=ff(1)/ff(2)*alf2*c4
    c7=exp(alf2*s3*zb(3))
    c8=exp(-alf2*s3*zb(3))
    c9=exp(alf3*s3*zb(3))
    c10=exp(-alf3*s3*zb(3))
    c11=ff(2)/ff(3)*alf3*c9
    c12=ff(2)/ff(3)*alf3*c10
    c13=exp(alf3*s3*zb(4))
    c14=exp(-alf3*s3*zb(4))
    c15=exp(alf4*s3*zb(4)) 
    c16=exp(-alf4*s3*zb(4))
    c17=ff(3)/ff(4)*alf4*c15
    c18=ff(3)/ff(4)*alf4*c16
    c19=exp(alf4*s3*zb(5))
    c20=exp(-alf4*s3*zb(5))
    c21=exp(alf5*s3*zb(5))
    c22=exp(-alf5*s3*zb(5))
    c23=ff(4)/ff(5)*alf5*c21 
    c24=ff(4)/ff(5)*alf5*c22
    c25=exp(alf5*s3*zb(6))
    c26=exp(-alf5*s3*zb(6))
    c27=exp(alf6*s3*zb(6)) 
    c28=exp(-alf6*s3*zb(6))
    c29=ff(5)/ff(6)*alf6*c27 
    c30=ff(5)/ff(6)*alf6*c28
    c31=exp(alf6*s3*zb(7))
    c32=exp(-alf6*s3*zb(7))
    c33=exp(alf7*s3*zb(7))
    c34=exp(-alf7*s3*zb(7))
    c35=ff(6)/ff(7)*alf7*c33
    c36=ff(6)/ff(7)*alf7*c34
    c37=exp(alf7*s3*zb(8))
    c38=exp(-alf7*s3*zb(8))

    cc=transpose(reshape((/1._wp,1._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            c1,c2,-c3,-c4,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            c1,-c2,-c5,c6,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,c7,c8,-c9,-c10,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,alf2*c7,-alf2*c8,-c11,c12,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,c13,c14,-c15,-c16,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,alf3*c13,-alf3*c14,-c17,c18,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,c19,c20,-c21,-c22,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,alf4*c19,-alf4*c20,-c23,c24,0._wp,0._wp,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,c25,c26,-c27,-c28,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,alf5*c25,-alf5*c26,-c29,c30,0._wp,0._wp,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,c31,c32,-c33,-c34,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,alf6*c31,-alf6*c32,-c35,c36,&
                            0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,c37,c38/),(/14,14/)))


    x = (/co,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp/)                            
    CALL DGESV(14,1,cc,14, ipiv,x,14,info )

    kcal= kc*pperccaccoo3*rhmean/(mmca*ccoo3llhhs(i))

    ttcaldis1 = kcal*((1._wp-fi(1))*(-x(1)/s3*(exp(s3*zb(2))-1._wp)-x(2)/(-s3)*(exp(-s3*zb(2))-1._wp)))
    ttcaldis2 = kcal*((1._wp-fi(2))*(-x(3)/(alf2*s3)*(exp(alf2*s3*zb(3))-exp(alf2*s3*zb(2)))-&
                x(4)/(-alf2*s3)*(exp(-alf2*s3*zb(3))-exp(-alf2*s3*zb(2)))))
    ttcaldis3 = kcal*((1._wp-fi(3))*(-x(5)/(alf3*s3)*(exp(alf3*s3*zb(4))-exp(alf3*s3*zb(3)))-&
                x(6)/(-alf3*s3)*(exp(-alf3*s3*zb(4))-exp(-alf3*s3*zb(3)))))            
    ttcaldis4 = kcal*((1._wp-fi(4))*(-x(7)/(alf4*s3)*(exp(alf4*s3*zb(5))-exp(alf4*s3*zb(4)))-&
                x(8)/(-alf4*s3)*(exp(-alf4*s3*zb(5))-exp(-alf4*s3*zb(4)))))
    ttcaldis5 = kcal*((1._wp-fi(5))*(-x(9)/(alf5*s3)*(exp(alf5*s3*zb(6))-exp(alf5*s3*zb(5)))-&
                x(10)/(-alf5*s3)*(exp(-alf5*s3*zb(6))-exp(-alf5*s3*zb(5)))))
    ttcaldis6 = kcal*((1._wp-fi(6))*(-x(11)/(alf6*s3)*(exp(alf6*s3*zb(7))-exp(alf6*s3*zb(6)))-&
                x(12)/(-alf6*s3)*(exp(-alf6*s3*zb(7))-exp(-alf6*s3*zb(6)))))
    ttcaldis7 = kcal*((1._wp-fi(7))*(-x(13)/(alf7*s3)*(exp(alf7*s3*zb(8))-exp(alf7*s3*zb(7)))-&
                x(14)/(-alf7*s3)*(exp(-alf7*s3*zb(8))-exp(-alf7*s3*zb(7)))))

    ! Total calcite dissolution rate in the model sediment.
    ttcaldis = ttcaldis1+ttcaldis2+ttcaldis3+ttcaldis4+ttcaldis5+ttcaldis6+ttcaldis7

end if

del=rr*jo2*(1._wp-fi)*ff/ddo2

call evaluate_s_corg()
rres1 = corg1a7w_eq(zb(2))
rres2 = corg2a7w_eq(zb(3))
rres3 = corg3a7w_eq(zb(4))
rres4 = corg4a7w_eq(zb(5))
rres5 = corg5a7w_eq(zb(6))
rres6 = corg6a7w_eq(zb(7))
rres7 = corg7a7w_eq(zb(8))

if (rres1 .GT. 0.0_wp) then
    zo    = zbrent(corg1a7w_eq,0.001_wp  ,zb(2),1.0D-05)
    rres1 = corg1a7w_eq(zo)
    pt    = 1
elseif (rres2 .GT. 0.0_wp) then
    zo    = zbrent(corg2a7w_eq,zb(2),zb(3),1.0D-05)
    rres2 = corg2a7w_eq(zo)
    pt    = 2
elseif (rres3 .GT. 0.0_wp) then
    zo    = zbrent(corg3a7w_eq,zb(3),zb(4),1.0D-05)
    rres3 = corg3a7w_eq(zo)
    pt    = 3
elseif (rres4 .GT. 0.0_wp) then
    zo    = zbrent(corg4a7w_eq,zb(4),zb(5),1.0D-05)
    rres4 = corg4a7w_eq(zo)
    pt    = 4
elseif (rres5 .GT. 0.0_wp) then
    zo    = zbrent(corg5a7w_eq,zb(5),zb(6),1.0D-05)
    rres5 = corg5a7w_eq(zo)
    pt    = 5
elseif (rres6 .GT. 0.0_wp) then
    zo    = zbrent(corg6a7w_eq,zb(6),zb(7),1.0D-05)
    rres6 = corg6a7w_eq(zo)
    pt    = 6
elseif (rres7 .GT. 0.0_wp) then
    zo    = zbrent(corg7a7w_eq,zb(6),zb(8),1.0D-05)
    rres7 = corg7a7w_eq(zo)
    pt    = 7
else
    zo=zb(8)
    pt=8
end if

ttremin = rrm(1)+rrm(2)+rrm(3)+rrm(4)+rrm(5)+rrm(6)+rrm(7)
! Total organic C remineralization in the sediment. This is equal
! to total source rate to the bottom of the water column of CO2 ion,
! i.e. one times the source for DIC 


if ((ffcal-ttcaldis)*mmca +(ffcc-ttremin)*mmoc*compc+nncc-&
    (fimeanb-fimean)*zb(8)*rhmean/dtc >= 0.) then

    w=((ffcal-ttcaldis)*mmca +(ffcc-ttremin)*mmoc*compc+nncc-&
       (fimeanb-fimean)*zb(8)*rhmean/dtc)/((1._wp-fi(7))*rhmean)

    dpp=((ffcal-ttcaldis)*mmca*dtc-(w*(1._wp-fi(7))+zb(8)*&
         (fimeanb-fimean)/dtc)*rhmean*pperccaccoo3*dtc)/&
        ((1._wp-fimean)*zb(8)*rhmean)

    doorg=((ffcc-ttremin)*mmoc*compc*dtc-(w*(1._wp-fi(7))+&
          zb(8)*(fimeanb-fimean)/dtc)*rhmean*pperoorg*dtc)/&
          ((1._wp-fimean)*zb(8)*rhmean)
else

    w=((ffcal-ttcaldis)*mmca +(ffcc-ttremin)*mmoc*compc+nncc-&
       (fimeanb-fimean)*zb(8)*rhmean/dtc)/((1._wp-fiboto)*rhmeano)

    dpp=((ffcal-ttcaldis)*mmca*dtc-(w*(1._wp-fiboto)*pperccaccoo3o*&
          rhmeano+zb(8)*(fimeanb-fimean)/dtc*pperccaccoo3*&
          rhmean)*dtc)/((1._wp-fimean)*zb(8)*rhmean)

    doorg=((ffcc-ttremin)*mmoc*compc*dtc-(w*(1._wp-fiboto)*&
            pperoorgo*rhmeano+zb(8)*(fimeanb-fimean)/dtc*&
            pperoorg*rhmean)*dtc)/((1._wp-fimean)*zb(8)*rhmean)
end if

pperccaccoo3 = pperccaccoo3+dpp
pperoorg     = pperoorg+doorg
fimeanb      = fimean
! print*, ffcc
! Output parameters
! - percentage of calcite dissolution
pccal(i)   = ttcaldis/ffcal
! - percentage organic carbon dissolution
poorg(i)   = ttremin/ffcc
! - dry weight fraction of calcite in sediment
dwpccal(i) = pperccaccoo3
dwpoorg(i) = pperoorg
wsed(i)    = w
fim(i)     = fimeanb

end do
end subroutine smtdcorg

include 'CarSys_pres.f90'
include "Corg_Data.f90"
include "corg1a7w_eq.f90"
include "corg2a7w_eq.f90"
include "corg3a7w_eq.f90"
include "corg4a7w_eq.f90"
include "corg5a7w_eq.f90"
include "corg6a7w_eq.f90"
include "corg7a7w_eq.f90"
include "evaluate_s_corg.f90"
include "zbrent.f90"