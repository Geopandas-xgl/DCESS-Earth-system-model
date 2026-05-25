subroutine oc_carbon_sys_surface(hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,soc,arc,slso,&
                                    pco2wa,k0a,co2a,co3a,hco3a,omea,&
                                    pco2wp,k0p,co2p,co3p,hco3p,omep,&
                                    pco2war,k0ar,co2ar,co3ar,hco3ar,omear,&
                                    pco2wso,k0so,co2so,co3so,hco3so,omeso,&
                                    pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh)
! Calculates ocean carbon system at the ocean surface for all ocena boxes
! Inputs:
! Ocean boxes: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,soc,arc,slso
! Outputs:
! pco2w: partial pressure of water co2
! k0: Solubility of CO2
! co2: water co2
! co3: carbonate ion
! hco3: bicarbonate
! ome: calcite saturation state

use parameters, only: wp,nto,n,to_f,nob
implicit none
real(wp), dimension(nto,n), intent(in) :: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,soc,arc
real(wp), dimension(nto), intent(in) :: slso
real(wp),intent(out) :: pco2war,k0ar,co2ar,co3ar,hco3ar,omear,pco2wso,k0so,co2so,co3so,hco3so,omeso,&
                        pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh
real(wp), dimension(nob),intent(out) :: pco2wa,k0a,co2a,co3a,hco3a,omea,&
                                      pco2wp,k0p,co2p,co3p,hco3p,omep
real(wp) :: pco2w,k0,co2,co3,hco3,ome
real(wp), dimension(nto,n) :: psi

psi = arc
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2war,k0ar,co2ar,co3ar,hco3ar,omear)
psi = soc
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2wso,k0so,co2so,co3so,hco3so,omeso)

psi = hna
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wa(1) = pco2w; k0a(1) = k0; co2a(1) = co2; co3a(1) = co3; hco3a(1) = hco3; omea(1)  = ome
psi = mna
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wa(2) = pco2w; k0a(2) = k0; co2a(2) = co2; co3a(2) = co3; hco3a(2) = hco3; omea(2)  = ome
psi = ena
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wa(3) = pco2w; k0a(3) = k0; co2a(3) = co2; co3a(3) = co3; hco3a(3) = hco3; omea(3)  = ome
psi = esa
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wa(4) = pco2w; k0a(4) = k0; co2a(4) = co2; co3a(4) = co3; hco3a(4) = hco3; omea(4)  = ome
psi = msa
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wa(5) = pco2w; k0a(5) = k0; co2a(5) = co2; co3a(5) = co3; hco3a(5) = hco3; omea(5)  = ome

psi = hnp
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wp(1) = pco2w; k0p(1) = k0; co2p(1) = co2; co3p(1) = co3; hco3p(1) = hco3; omep(1)  = ome
psi = mnp
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wp(2) = pco2w; k0p(2) = k0; co2p(2) = co2; co3p(2) = co3; hco3p(2) = hco3; omep(2)  = ome
psi = enp
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wp(3) = pco2w; k0p(3) = k0; co2p(3) = co2; co3p(3) = co3; hco3p(3) = hco3; omep(3)  = ome
psi = esp
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wp(4) = pco2w; k0p(4) = k0; co2p(4) = co2; co3p(4) = co3; hco3p(4) = hco3; omep(4)  = ome
psi = msp
call carsys_sur(psi(1,1),psi(2,1),psi(5,1),psi(8,1),pco2w,k0,co2,co3,hco3,ome)
pco2wp(5) = pco2w; k0p(5) = k0; co2p(5) = co2; co3p(5) = co3; hco3p(5) = hco3; omep(5)  = ome

call carsys_sur(to_f,slso(2),slso(5),slso(8),pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh)


end subroutine oc_carbon_sys_surface

include 'CarSys_sur.f90'