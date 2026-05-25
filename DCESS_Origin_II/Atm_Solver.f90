subroutine atm_solver(at,totsw,totlw,aseat,asepa,aseso,asear,htot,hfde,qsh,&
                      alhn,almn,alen,ales,alms,alhs,mdrhn,mdrmn,mdren,mdres,mdrms,mdrhs,&
                      ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,asar,asso,assh,&
                      wcar,wsil,worg,vol,p14c,dfina,dfinp,dfis,kat)
! Calculates time evolution of atmospheric tracers
! Inputs:
! at: atmospheric tracers
! total short and longwave radiation: totsw, totlw (W)
! air-sea heat exchange: aseat,asepa,aseso,asear,qsh (W)
! meridional heat transport: htot,hfde
! land biosphere greenhouse gases fluxes: alhn,almn,alen,ales,alms,alhs
! CH4-CO2 oxidation terms: mdrhn,mdrmn,mdren,mdres,mdrms,mdrhs
! air-sea gas exchange: ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,asar,asso,assh
! Weathering rates and volvanic inputs: wcar,wsil,worg,vol
! p14c: 14c atmospheric production
! sea-ice line rate changes: dfina,dfinp,dfis
! Outputs:
! kat: time evolution of atmospheric tracers

use parameters, only: wp,nta,rc,aa_h,aa_m,aa_e,tn2o,r13pdb,d13corgcpi,d13volpi,dc14,avg,rvahp,rvahm,rvahe,sy,nab,nob
implicit none
real(wp), dimension(nta,nab), intent(in) :: at
real(wp), dimension(nab), intent(in) :: totsw,totlw,wcar,wsil,worg,vol,p14c
real(wp), intent(in) :: aseat(nob),asepa(nob),aseso,asear,htot(4),hfde,qsh,dfina,dfinp,dfis
real(wp), dimension(8), intent(in) :: alhn,almn,alen,ales,alms,alhs
real(wp), dimension(2), intent(in) :: mdrhn,mdrmn,mdren,mdres,mdrms,mdrhs
real(wp), dimension(9), intent(in) :: ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,asar,asso,assh
real(wp), dimension(nta,nab), intent(out) :: kat
real(wp) :: ase(nab),gam(nab),hdiv(nab),atm_a(nab)
real(wp), parameter :: datmp  =  5._wp,&                ! Equivalent water depth of atmosphere, high, mid and low lat
                       datmm  =  10._wp,&
                       datme  =  20._wp,&
                       tice   =  3.e-2_wp                   ! Inverse timescale for sea-ice changes  
integer :: it
real(wp) :: tgas(nab-1),dgas(nab),ss(nab),n2o_dec(nab),f13c_vol,f13c_org,f1314f12r(nab)

kat =0._wp

!----------------------------------------------------------------------
! Temp
ase =(/ asear+aseat(1)+asepa(1), aseat(2)+asepa(2), aseat(3)+asepa(3),&
              aseat(4)+asepa(4), aseat(5)+asepa(5), aseso /)
atm_a = (/ aa_h,aa_m,aa_e,aa_e,aa_m,aa_h/)

gam = rc*atm_a*(/ datmp,datmm,datme,datme,datmm,datmp /)
hdiv = (/ htot(4),htot(3)-htot(4),-htot(3)+hfde,-htot(2)-hfde,htot(2)-htot(1), htot(1) /)

kat(1,:) = 1._wp/gam*( totsw - totlw - ase + hdiv ) 
kat(1,6) = kat(1,6) + 1._wp/gam(6)*qsh

!----------------------------------------------------------------------
! 12-CH4
it=3
tgas = transport(at(it,:))
dgas = gas_div(tgas)
ss(1) = -alhn(it) + mdrhn(1) + dgas(1)
ss(2) = -almn(it) + mdrmn(1) + dgas(2)
ss(3) = -alen(it) + mdren(1) + dgas(3)
ss(4) = -ales(it) + mdres(1) + dgas(4)
ss(5) = -alms(it) + mdrms(1) + dgas(5)
ss(6) = -alhs(it) + mdrhs(1) + dgas(6)

gam = (/rvahp,rvahm,rvahe,rvahe,rvahm,rvahp/)

kat(it,:) = -1._wp/gam*ss

!----------------------------------------------------------------------
! N2O
it=4
tgas = transport(at(it,:))
dgas = gas_div(tgas)
n2o_dec = at(it,:)*gam/(tn2o*sy)

ss(1) = -alhn(it) + n2o_dec(1) + dgas(1)
ss(2) = -almn(it) + n2o_dec(2) + dgas(2)
ss(3) = -alen(it) + n2o_dec(3) + dgas(3)
ss(4) = -ales(it) + n2o_dec(4) + dgas(4)
ss(5) = -alms(it) + n2o_dec(5) + dgas(5)
ss(6) = -alhs(it) + n2o_dec(6) + dgas(6)

kat(it,:) = -1._wp/gam*ss

!----------------------------------------------------------------------
! 12-CO2
it=5
tgas = transport(at(it,:))
dgas = gas_div(tgas)
ase = (/ asar(it)+ashna(it)+ashnp(it),asmna(it)+asmnp(it),asena(it)+asenp(it),&
                  asesa(it)+asesp(it),asmsa(it)+asmsp(it),asso(it)+assh(it) /)

ss(1) = ase(1) - alhn(it) - mdrhn(1) + wcar(1) + 2._wp*wsil(1) - worg(1) - vol(1) + dgas(1)
ss(2) = ase(2) - almn(it) - mdrmn(1) + wcar(2) + 2._wp*wsil(2) - worg(2) - vol(2) + dgas(2)
ss(3) = ase(3) - alen(it) - mdren(1) + wcar(3) + 2._wp*wsil(3) - worg(3) - vol(3) + dgas(3)
ss(4) = ase(4) - ales(it) - mdres(1) + wcar(4) + 2._wp*wsil(4) - worg(4) - vol(4) + dgas(4)
ss(5) = ase(5) - alms(it) - mdrms(1) + wcar(5) + 2._wp*wsil(5) - worg(5) - vol(5) + dgas(5)
ss(6) = ase(6) - alhs(it) - mdrhs(1) + wcar(6) + 2._wp*wsil(6) - worg(6) - vol(6) + dgas(6)

kat(it,:) = -1._wp/gam*ss

!----------------------------------------------------------------------
! 13-CO2
it = 6
tgas = transport(at(it,:))
dgas = gas_div(tgas)
ase = (/ asar(it)+ashna(it)+ashnp(it),asmna(it)+asmnp(it),asena(it)+asenp(it),&
                  asesa(it)+asesp(it),asmsa(it)+asmsp(it),asso(it)+assh(it) /)
f13c_vol = r13pdb*(d13volpi   + 1._wp)
f13c_org = r13pdb*(d13corgcpi + 1._wp)
f1314f12r = at(6,:)/at(5,:)
ss(1) = ase(1) - alhn(it) - mdrhn(2) + (wcar(1)+2._wp*wsil(1))*f1314f12r(1) - worg(1)*f13c_org - vol(1)*f13c_vol + dgas(1)
ss(2) = ase(2) - almn(it) - mdrmn(2) + (wcar(2)+2._wp*wsil(2))*f1314f12r(2) - worg(2)*f13c_org - vol(2)*f13c_vol + dgas(2)
ss(3) = ase(3) - alen(it) - mdren(2) + (wcar(3)+2._wp*wsil(3))*f1314f12r(3) - worg(3)*f13c_org - vol(3)*f13c_vol + dgas(3)
ss(4) = ase(4) - ales(it) - mdres(2) + (wcar(4)+2._wp*wsil(4))*f1314f12r(4) - worg(4)*f13c_org - vol(4)*f13c_vol + dgas(4)
ss(5) = ase(5) - alms(it) - mdrms(2) + (wcar(5)+2._wp*wsil(5))*f1314f12r(5) - worg(5)*f13c_org - vol(5)*f13c_vol + dgas(5)
ss(6) = ase(6) - alhs(it) - mdrhs(2) + (wcar(6)+2._wp*wsil(6))*f1314f12r(6) - worg(6)*f13c_org - vol(6)*f13c_vol + dgas(6)

kat(it,:) = -1._wp/gam*ss

!----------------------------------------------------------------------
! 14-CO2
it = 7
tgas = transport(at(it,:))
dgas = gas_div(tgas)
ase = (/ asar(it)+ashna(it)+ashnp(it),asmna(it)+asmnp(it),asena(it)+asenp(it),&
                  asesa(it)+asesp(it),asmsa(it)+asmsp(it),asso(it)+assh(it) /)
f1314f12r = at(7,:)/at(5,:)

ss(1) = ase(1) - alhn(it) + (wcar(1) + 2._wp*wsil(1))*f1314f12r(1) + dgas(1)
ss(2) = ase(2) - almn(it) + (wcar(2) + 2._wp*wsil(2))*f1314f12r(2) + dgas(2)
ss(3) = ase(3) - alen(it) + (wcar(3) + 2._wp*wsil(3))*f1314f12r(3) + dgas(3)
ss(4) = ase(4) - ales(it) + (wcar(4) + 2._wp*wsil(4))*f1314f12r(4) + dgas(4)
ss(5) = ase(5) - alms(it) + (wcar(5) + 2._wp*wsil(5))*f1314f12r(5) + dgas(5)
ss(6) = ase(6) - alhs(it) + (wcar(6) + 2._wp*wsil(6))*f1314f12r(6) + dgas(6)

kat(it,:) = -1._wp/gam*ss - dc14*at(it,:) + atm_a*p14c/gam*1._wp/avg

!----------------------------------------------------------------------
! 13-CH4
it = 8
tgas = transport(at(it,:))
dgas = gas_div(tgas)

ss(1) = - alhn(it) + mdrhn(2) + dgas(1)
ss(2) = - almn(it) + mdrmn(2) + dgas(2)
ss(3) = - alen(it) + mdren(2) + dgas(3)
ss(4) = - ales(it) + mdres(2) + dgas(4)
ss(5) = - alms(it) + mdrms(2) + dgas(5)
ss(6) = - alhs(it) + mdrhs(2) + dgas(6)

kat(it,:) = -1._wp/gam*ss

!----------------------------------------------------------------------
! Oxygen

!----------------------------------------------------------------------
! Sea-ice
it = 10
kat(it,1) = dfina*tice ! Arctic
kat(it,2) = dfinp*tice ! North Pacific
kat(it,6) = dfis *tice ! Southern Ocean


contains
function transport(gas) result(tgas)
implicit none
real(wp), dimension(6), intent(in) :: gas(6)
real(wp) :: gn,gmn,gen,ges,gms,gs,k,tgas(5)

gn  = gas(1)
gmn = gas(2)
gen = gas(3)
ges = gas(4)
gms = gas(5)
gs  = gas(6)

k = 6.e-8_wp ! transfer coefficient (1/s)

tgas(1) = -k*(gmn - gn ) ! atm/s
tgas(2) = -k*(gen - gmn)
tgas(3) =  k*(ges - gen)
tgas(4) = -k*(ges - gms)
tgas(5) = -k*(gms - gs )

end function transport


function gas_div(tgas) result(gdiv)
implicit none
real(wp), intent(in) :: tgas(5)
real(wp) :: gdiv(6)

gdiv = (/                  tgas(1)*rvahm,&
            tgas(2)*rvahe - tgas(1)*rvahm,&
            -tgas(2)*rvahe - tgas(3)*rvahe,&
            -tgas(4)*rvahe + tgas(3)*rvahe,&
            tgas(4)*rvahe - tgas(5)*rvahm,&
                            tgas(5)*rvahm  /)

end function gas_div

end subroutine atm_solver