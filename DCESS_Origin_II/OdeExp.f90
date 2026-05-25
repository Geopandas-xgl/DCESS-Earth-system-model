subroutine odeexp(t,at,taam,lbn,lbs,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc,slso,&
                  porga,porgp,porgar,porgso,pcala,pcalp,pcalar,pcalso,&
                  inso,oc_pres,shl_dat,olf_bios,gaa,gap,gaar,gaso,&
                  kat,klbn,klbs,&
                  karc,khna,kmna,kena,kesa,kmsa,khnp,kmnp,kenp,kesp,kmsp,ksoc,kslso)

! Solve differential equations for DCESS II
! Inputs:
! t: time
! at: atmospheric tracers
! taam: annual mean atmospheric temperature
! lbn, lbs: land biosphere tracers NH and SH
! Ocean tracers: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc
! Ocean shelf tracers: slso
! porg: percentage organic carbon rem (from sediment)
! pcal: percentage organic calcite dissolution
! inso,oc_pres,shl_dat,olf_bios: solar radiation, ocean pressure, shelf data and ocean land fraction meridional profile
! gaa,gap,gaar,gaso: vertical profile of ocean area fraction
! Outputs: 
! kat: Rate of change (per sec) of atmospheric tracers
! klbn,klbs: Rate of change (per sec) of Land Biosphere tracers
! Rate of change (per sec) of ocean tracers: karc,khna,kmna,kena,kesa,kmsa,khnp,kmnp,kenp,kesp,kmsp,ksoc
! kslso: rate of change (per sec) pf ocean shelf tracers
!--------------------------------------------
! Data ordered as:
!--------------------------------------------
! Atmospheric tracers: at
!                       1 Temperature (°C)
!                       3 pCH4 (atm)
!                       4 pN2O (atm)
!                       5 p12CO2 (atm)
!                       6 p13CO2 (atm)
!                       7 p14CO2 (atm)
!                       8 p13CH4 (atm)
!                       9 Oxygen (atm)
!                       10 Sea ice lines (rad)
!--------------------------------------------
! Ocean tracers: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc and slso (look at Fig1 in manuscript)
!                       1 Temperature (°C)
!                       2 Salinity (g/kg)
!                       3 Phosphate (mol/m3)
!                       5 DIC-12 (mol/m3)
!                       6 DIC-13 (mol/m3)
!                       7 DIC-14 (mol/m3)
!                       8 ALK (mol/m3)
!                       9 Dissolved Oxygen (mol/m3)
!                       10 Isotope ration d18O (-))
!--------------------------------------------
! Land tracers
!                       1  - 12: ^12C
!                       1  - 4 : Tropical forest
!                       5  - 8 : Grassland, savanna, desert
!                       9  - 12: Extratropical forest
!                       13 - 24: ^13C
!                       13 - 16: Tropical forest
!                       17 - 20: Grassland, savanna, desert
!                       21 - 24: Extratropical forest
!                       25 - 36: ^14C
!                       25 - 28: Tropical forest
!                       29 - 32: Grassland, savanna, desert
!                       33 - 36: Extratropical forest
!--------------------------------------------
use parameters, only: wp,nta,nto,nlb,n,nz_all,rc,sice,sgla,o18fg,fice,fgla_so,fgla_sh,evfs,brvf,d,dm,to_f,shl,&
ao_ar,ao_na,ao_mna,ao_ena,ao_esa,ao_msa,ao_np,ao_mnp,ao_enp,ao_esp,ao_msp,ao_so,ao_sh,&
rvahp,rvahm,rvahe,dc14,c_ar,c_na,c_mna,c_msa,c_sa,c_np,c_mnp,c_msp,c_sp,fw_z,nob,nab

implicit none
real(wp), intent(in) :: t
real(wp), intent(in) , dimension(nta,nab) :: at
real(wp), intent(in) , dimension(nto,n) :: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc
real(wp), intent(in) , dimension(nto) :: slso
real(wp), intent(in) :: inso(19,14),oc_pres(n,7),shl_dat(nz_all,3),olf_bios(1801,2)
real(wp), intent(in) :: lbn(nlb),lbs(nlb),taam(nab)
real(wp), intent(in) :: gaa(nob,n),gap(nob,n),gaar(n),gaso(n)
real(wp), intent(in), dimension(nob,n) :: pcala,pcalp,porga,porgp
real(wp), intent(in), dimension(n) :: pcalar,pcalso,porgar,porgso
! Atm
real(wp) :: to_a(nob),to_p(nob),to_ar,to_so
real(wp), dimension(4) :: hsen,hlat,htot
real(wp) :: fwfbr,fwfdhn,fwfdln,fwfdls,fwfdhs,hfde
real(wp), dimension(nab) :: totsw,totlw
real(wp) :: asear,aseat(nob),asepa(nob),aseso,fsnn,fsns
real(wp) :: dfina,dfinp,dfis,aoifa(nob),aoifp(nob),aoifar,aoifso
real(wp), dimension(3) :: ptanh,ptash
! Water vapor
real(wp) :: o18a(nob),o18p(nob-1),o18_z
! Oce
real(wp), dimension(n-1):: kvhna,kvmna,kvena,kvesa,kvmsa,&
                           kvhnp,kvmnp,kvenp,kvesp,kvmsp,kvar,kvso
real(wp), dimension(n) ::  khfbra,khfdhna,khfdlna,khfdea,khfdlsa,khfdhsa,&
                           khfdhnp,khfdlnp,khfdep,khfdlsp,khfdhsp,khmsap,&
                           qfbra,qfdhna,qfdhnp,qfdlna,qfdlnp,qfdea,qfdep,&
                           qfdlsa,qfdlsp,qfdhsa,qfdhsp,ievfsa,ievfsp,ievfbr,&
                           whna,wmna,wena,wesa,wmsa,whnp,wmnp,wenp,wesp,wmsp,wso,war
real(wp) :: fdsup,feso(n),rfds_v(n),qshl
! Land biosphere
real(wp) :: nppn(3),npps(3),aln(nlb),als(nlb),alhn(8),almn(8),alen(8),ales(8),alms(8),alhs(8)
! External Forcing
real(wp) :: rcarat(nob),rcarpa(nob),rcarar,rcarso,rorgat(nob),rorgpa(nob),rorgar,rorgso,wcar(nab),wsil(nab),worg(nab),vol(nab)
! Atm_Met, C14 atm
real(wp), dimension(2) :: mdrn,mdrmn,mdren,mdres,mdrms,mdrs
real(wp), dimension(nab) :: p14c
! CarSys
real(wp),dimension(nob) :: pco2wa,k0a,co2a,co3a,hco3a,omea,pco2wp,k0p,co2p,co3p,hco3p,omep
real(wp) :: pco2war,k0ar,co2ar,co3ar,hco3ar,omear,pco2wso,k0so,co2so,co3so,hco3so,omeso,&
            pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh
! Gas Exc
real(wp), dimension(9) :: ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,asar,asso,assh
real(wp) :: co2_atm(3),o2_atm
! OrgFlx
real(wp), dimension(nto,n) :: srar,srso,srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp
real(wp), dimension(nto) :: shv,glv
real(wp) :: zerot(nto,n),zerop(n)
real(wp), dimension(nto,n), intent(out) :: karc,khna,kmna,kena,kesa,kmsa,khnp,kmnp,kenp,kesp,kmsp,ksoc
real(wp), intent(out) :: kslso(nto),kat(nta,nab),klbn(nlb),klbs(nlb)
real(wp) :: varc(n),vmsa(n),vmsp(n),vsoc(n),vshl,asshl(nto)
integer :: i,j,ent_lim(2),idlsa,idlsp,idlar,inx(4),zlim,nzin
integer, dimension(:), allocatable :: zin
real(wp) :: qsom,qarm,qm(nob)

to_a  = (/hna(1,1),mna(1,1),ena(1,1),esa(1,1),msa(1,1)/)
to_p  = (/hnp(1,1),mnp(1,1),enp(1,1),esp(1,1),msp(1,1)/)
to_ar = arc(1,1)
to_so = soc(1,1)

! Atmospheric energy balance
call atm(t,at,to_a,to_p,to_ar,to_so,inso,&
         hsen,hlat,htot,hfde,fwfbr,fwfdhn,fwfdln,fwfdls,fwfdhs,&
         totsw,totlw,asear,aseat,asepa,aseso,&
         dfina,dfinp,dfis,aoifa,aoifp,aoifar,aoifso,ptanh,ptash,fsnn,fsns)

! Water 18O
call water_vapor_18o(hna,mna,ena,esa,msa,mnp,enp,esp,msp,ptanh,ptash,at(1,3),&
                     o18a,o18p,o18_z)

! Ocean dynamics
call oce(arc,soc,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,slso,&
         fwfbr,fwfdhn,fwfdln,fwfdls,fwfdhs,oc_pres,shl_dat,&
         kvhna,kvmna,kvena,kvesa,kvmsa,kvhnp,kvmnp,kvenp,kvesp,kvmsp,kvar,kvso,&
         khfbra,khfdhna,khfdlna,khfdea,khfdlsa,khfdhsa,khfdhnp,khfdlnp,khfdep,khfdlsp,khfdhsp,khmsap,&
         qfbra,qfdhna,qfdhnp,qfdlna,qfdlnp,qfdea,qfdep,qfdlsa,qfdlsp,qfdhsa,qfdhsp,&
         ievfsa,idlsa,ievfsp,idlsp,ievfbr,idlar,fdsup,feso,ent_lim,rfds_v,&
         whna,wmna,wena,wesa,wmsa,whnp,wmnp,wenp,wesp,wmsp,wso,war,qshl)

! Land biosphere
call land_bios(t,taam,ptanh,ptash,at,lbn,lbs,olf_bios,&
               nppn,npps,aln,als,alhn,almn,alen,ales,alms,alhs)

! External Forcing
call external_forcing(taam,rcarat,rcarpa,rcarar,rcarso,rorgat,rorgpa,rorgar,rorgso,wcar,wsil,worg,vol)

! Atmospheric methane sinks
call atm_met(at,mdrn,mdrmn,mdren,mdres,mdrms,mdrs)

! Atmospheric 14C production
call atm_c14_prod(at,p14c)

! Ocean carbon system
! Iterate Carbonate system for pCO2 - partial pressure of surface CO2, and K0 - solubility of CO2
call oc_carbon_sys_surface(hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,soc,arc,slso,&
                           pco2wa ,k0a ,co2a ,co3a ,hco3a ,omea ,&
                           pco2wp ,k0p ,co2p ,co3p ,hco3p ,omep ,&
                           pco2war,k0ar,co2ar,co3ar,hco3ar,omear,&
                           pco2wso,k0so,co2so,co3so,hco3so,omeso,&
                           pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh)


! Gas Exchange
call gas_exchange(at,arc,soc,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,slso,&
                  pco2wa,pco2wp,pco2war,pco2wso,pco2wsh,k0a,k0p,k0ar,k0so,k0sh,&
                  co3a,co3p,co3ar,co3so,co3sh,aoifa,aoifp,aoifar,aoifso,&
                  ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,&
                  asar,asso,assh,co2_atm,o2_atm)

! Biogeochemical source/sinks
call ocean_rad(t,inso,qsom,qarm,qm)

call oce_srcsnk(hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc,&
                aoifa,aoifp,aoifar,aoifso,gaa,gap,gaar,gaso,&
                omea,omep,omear,omeso,co2a,co2p,co2ar,co2so,&
                co3a,co3p,co3ar,co3so,porga,porgp,porgar,porgso,&
                pcala,pcalp,pcalar,pcalso,qsom,qarm,qm,&
                srar,srso,srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp)


! Add air-sea heat flux
srar(1,1)  = asear/rc

srhna(1,1) = aseat(1)/rc
srmna(1,1) = aseat(2)/rc
srena(1,1) = aseat(3)/rc
sresa(1,1) = aseat(4)/rc
srmsa(1,1) = aseat(5)/rc

srhnp(1,1) = asepa(1)/rc
srmnp(1,1) = asepa(2)/rc
srenp(1,1) = asepa(3)/rc
sresp(1,1) = asepa(4)/rc
srmsp(1,1) = asepa(5)/rc

srso(1,1)  = aseso/rc

! Southern Ocean and North Pacific surface volume fluxes balance
shv = (/0._wp,sice,slso(3),0._wp,slso(5:nto)/)
glv = (/0._wp,sgla,0._wp  ,0._wp,co2_atm,0._wp,o2_atm,o18fg/)
do i=1,nto
    srso(i,1)  = srso(i,1)  + fice*shv(i) + fgla_so*glv(i) - fgla_sh*glv(i) - evfs*soc(i,1)
    srhnp(i,1) = srhnp(i,1) - brvf*hnp(i,1)
end do

! Interior sources for DI^14C
! Radiocarbon decay is given via the decay constant (dc14).
srhna(7,:)  = srhna(7,:) - dc14*hna(7,:)*ao_na *gaa(1,:)*d
srmna(7,:)  = srmna(7,:) - dc14*mna(7,:)*ao_mna*gaa(2,:)*d
srena(7,:)  = srena(7,:) - dc14*ena(7,:)*ao_ena*gaa(3,:)*d
sresa(7,:)  = sresa(7,:) - dc14*esa(7,:)*ao_esa*gaa(4,:)*d
srmsa(7,:)  = srmsa(7,:) - dc14*msa(7,:)*ao_msa*gaa(5,:)*d

srhnp(7,:)  = srhnp(7,:) - dc14*hnp(7,:)*ao_np *gap(1,:)*d
srmnp(7,:)  = srmnp(7,:) - dc14*mnp(7,:)*ao_mnp*gap(2,:)*d
srenp(7,:)  = srenp(7,:) - dc14*enp(7,:)*ao_enp*gap(3,:)*d
sresp(7,:)  = sresp(7,:) - dc14*esp(7,:)*ao_esp*gap(4,:)*d
srmsp(7,:)  = srmsp(7,:) - dc14*msp(7,:)*ao_msp*gap(5,:)*d

srar(7,:)  = srar(7,:)   - dc14*arc(7,:) *ao_ar *gaar*d
srso(7,:)  = srso(7,:)   - dc14*soc(7,:) *ao_so *gaso*d

! Combine surface fluxes and interior sources for DIC-12,13,14 and O2

inx = (/5,6,7,9/)
do j=1,4
    i=inx(j)
    srar(i,1)  = srar(i,1)  + asar(i)
    srhna(i,1) = srhna(i,1) + ashna(i)
    srmna(i,1) = srmna(i,1) + asmna(i)
    srena(i,1) = srena(i,1) + asena(i)
    sresa(i,1) = sresa(i,1) + asesa(i)
    srmsa(i,1) = srmsa(i,1) + asmsa(i)

    srhnP(i,1) = srhnP(i,1) + ashnp(i)
    srmnP(i,1) = srmnP(i,1) + asmnp(i)
    srenP(i,1) = srenP(i,1) + asenp(i)
    sresP(i,1) = sresP(i,1) + asesp(i)
    srmsP(i,1) = srmsP(i,1) + asmsp(i)

    srso(i,1)  = srso(i,1)  + asso(i)
end do

i=10
srar(i,1)  = srar(i,1)                                                + c_ar*fwfbr*o18a(1)
srhna(i,1) = srhna(i,1)                        + c_na *fwfdhn*o18a(2) - c_ar*fwfbr*o18a(1)

srmna(i,1) = srmna(i,1) + c_mna*fwfdln*o18a(3) - c_na *fwfdhn*o18a(2)
srena(i,1) = srena(i,1) - c_mna*fwfdln*o18a(3) - fw_z*o18_z

sresa(i,1) = sresa(i,1) - c_msa*fwfdls*o18a(4)
srmsa(i,1) = srmsa(i,1) + c_msa*fwfdls*o18a(4) - c_sa*fwfdhs*o18a(5)

srhnp(i,1) = srhnp(i,1)                        + c_np*fwfdhn*o18p(1)
srmnp(i,1) = srmnp(i,1) + c_mnp*fwfdln*o18p(2) - c_np*fwfdhn*o18p(1)
srenp(i,1) = srenp(i,1) - c_mnp*fwfdln*o18p(2) + fw_z*o18_z

sresp(i,1) = sresp(i,1) - c_msp*fwfdls*o18p(3)
srmsp(i,1) = srmsp(i,1) + c_msp*fwfdls*o18p(3) - c_sp*fwfdhs*o18p(4)

srso(i,1) = srso(i,1)   + c_sa*fwfdhs*o18a(5) + c_sp*fwfdhs*o18p(4)



! River input
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


! Bathymetry
kvar = kvar*gaar(1:n-1)
kvso = kvso*gaso(1:n-1)

kvhna = kvhna*gaa(1,1:n-1)
kvmna = kvmna*gaa(2,1:n-1)
kvena = kvena*gaa(3,1:n-1)
kvesa = kvesa*gaa(4,1:n-1)
kvmsa = kvmsa*gaa(5,1:n-1)

kvhnp = kvhnp*gap(1,1:n-1)
kvmnp = kvmnp*gap(2,1:n-1)
kvenp = kvenp*gap(3,1:n-1)
kvesp = kvesp*gap(4,1:n-1)
kvmsp = kvmsp*gap(5,1:n-1)

khfbra  = khfbra *gaa(1,1:n)
khfdhna = khfdhna*gaa(2,1:n)
khfdlna = khfdlna*gaa(3,1:n)
khfdea  = khfdea *gaa(4,1:n)
khfdlsa = khfdlsa*gaa(5,1:n)
khfdhsa = khfdhsa*gaa(5,1:n)

khfdhnp = khfdhnp*gap(2,1:n)
khfdlnp = khfdlnp*gap(3,1:n)
khfdep  = khfdep *gap(4,1:n)
khfdlsp = khfdlsp*gap(5,1:n)
khfdhsp = khfdhsp*gap(5,1:n)

khmsap = khmsap*gaa(5,1:n)


! Solve ocean diff-adv equation
zerot = 0._wp
zerop = 0._wp

call adv_dif_sol(arc,zerot,hna,kvar ,zerop  ,khfbra ,war ,zerop ,qfbra ,srar ,ao_ar ,gaar    ,karc)

call adv_dif_sol(hna,arc  ,mna,kvhna,khfbra ,khfdhna,whna,qfbra ,qfdhna,srhna,ao_na ,gaa(1,:),khna)
call adv_dif_sol(mna,hna  ,ena,kvmna,khfdhna,khfdlna,wmna,qfdhna,qfdlna,srmna,ao_mna,gaa(2,:),kmna)
call adv_dif_sol(ena,mna  ,esa,kvena,khfdlna,khfdea ,wena,qfdlna,qfdea ,srena,ao_ena,gaa(3,:),kena)
call adv_dif_sol(esa,ena  ,msa,kvesa,khfdea ,khfdlsa,wesa,qfdea ,qfdlsa,sresa,ao_esa,gaa(4,:),kesa)
call adv_dif_sol(msa,esa  ,soc,kvmsa,khfdlsa,khfdhsa,wmsa,qfdlsa,qfdhsa,srmsa,ao_msa,gaa(5,:),kmsa)

call adv_dif_sol(hnp,zerot,mnp,kvhnp,zerop  ,khfdhnp,whnp,zerop ,qfdhnp,srhnp,ao_np ,gap(1,:),khnp)
call adv_dif_sol(mnp,hnp  ,enp,kvmnp,khfdhnp,khfdlnp,wmnp,qfdhnp,qfdlnp,srmnp,ao_mnp,gap(2,:),kmnp)
call adv_dif_sol(enp,mnp  ,esp,kvenp,khfdlnp,khfdep ,wenp,qfdlnp,qfdep ,srenp,ao_enp,gap(3,:),kenp)
call adv_dif_sol(esp,enp  ,msp,kvesp,khfdep ,khfdlsp,wesp,qfdep ,qfdlsp,sresp,ao_esp,gap(4,:),kesp)
call adv_dif_sol(msp,esp  ,soc,kvmsp,khfdlsp,khfdhsp,wmsp,qfdlsp,qfdhsp,srmsp,ao_msp,gap(5,:),kmsp)

call adv_dif_sol_so(soc,msp,msa,zerot,kvso,khfdhsp,khfdhsa,zerop,wso,qfdhsp,qfdhsa,zerop,srso,ao_so,gaso,ksoc)

! Add extra volume fluxes
varc = ao_ar *(/dm, d*gaar(2:n) /)
vmsa = ao_msa*(/dm, d*gaa(5,2:n)/)
vmsp = ao_msp*(/dm, d*gap(5,2:n)/)
vsoc = ao_so *(/dm, d*gaso(2:n) /)

do i=1,nto
    do j=1,n
        karc(i,j) = karc(i,j) + 1._wp/varc(j)*( ievfbr(j)*hnp(i,1) )
        kmsa(i,j) = kmsa(i,j) + 1._wp/vmsa(j)*( ievfsa(j)*soc(i,1) + khmsap(j)*(msp(i,j) - msa(i,j)) )
        kmsp(i,j) = kmsp(i,j) + 1._wp/vmsp(j)*( ievfsp(j)*soc(i,1) - khmsap(j)*(msp(i,j) - msa(i,j)) )
        ksoc(i,j) = ksoc(i,j) + 1._wp/vsoc(j)*( -feso(j)*soc(i,j)  - rfds_v(j)*soc(i,j) )
    end do
end do

! Input flux to SO at z-end of gravity current
zlim = ent_lim(1) ! reached depth of gravity current
nzin = ent_lim(2) ! # of layers from zlim
allocate(zin(nzin))
zin = (/(i , i=zlim-nzin+1,zlim)/)
shv = (/to_f,slso(2:nto)/)

do i=1,nto
    ksoc(i,zin) = ksoc(i,zin) + 1._wp/vsoc(zin)*(  ( fdsup*shv(i) + sum(feso*soc(i,:)) )/nzin )
end do
deallocate(zin)

! Shelf mass balance
kslso=0._wp
shv = (/0._wp,sice,slso(3),0._wp,slso(5:nto)/)

vshl = d*ao_sh*sum(gaso(1:shl))
asshl = (/0._wp,0._wp,0._wp,0._wp,assh(5:7),0._wp,assh(9),0._wp/)
do i=2,nto
    kslso(i) = 1._wp/vshl*( fgla_sh*glv(i) - fice*shv(i) - fgla_so*glv(i) + sum(rfds_v*soc(i,:)) - fdsup*slso(i) + asshl(i) )
end do


! Solve Atmospheric tracers
call atm_solver(at,totsw,totlw,aseat,asepa,aseso,asear,htot,hfde,qshl,&
                alhn,almn,alen,ales,alms,alhs,mdrn,mdrmn,mdren,mdres,mdrms,mdrs,&
                ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,asar,asso,assh,&
                wcar,wsil,worg,vol,p14c,dfina,dfinp,dfis,kat)

! Solve land biosphere tracers
klbn = aln
klbs = als


end subroutine odeexp

include 'Atm.f90'
include 'Water_Vapor_18O.f90'
include 'Oce.f90'
include 'Land_Bios.f90'
include 'External_Forcing.f90'
include 'Atm_Met.f90'
include 'Atm_C14_Prod.f90'
include 'Oc_Carbon_Sys_Surface.f90'
include 'Gas_Exchange.f90'
include 'Ocean_Rad.f90'
include 'Oce_SrcSnk.f90'
include 'Riv_In.f90'
include 'Adv_Dif_Sol.f90'
include 'Adv_Dif_Sol_SO.f90'
include 'Atm_Solver.f90'
