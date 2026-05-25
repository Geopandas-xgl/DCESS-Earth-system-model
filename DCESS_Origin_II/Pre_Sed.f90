subroutine pre_sed(at,taam,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,slso,&
                   aoif_at,aoif_pa,aoif_ar,aoif_so,rorgat,rorgpa,rorgar,rorgso,&
                   co2a,co3a,omea,co2p,co3p,omep,co2ar,co3ar,omear,co2so,co3so,omeso)
! Calculates preliminar variables to be used in the sediment module

use parameters, only: wp,nto,n,nta,lf,lfar,lfso,&
ao_ar,ao_na,ao_mna,ao_ena,ao_esa,ao_msa,ao_np,ao_mnp,ao_enp,ao_esp,ao_msp,ao_so,nab,nob

implicit none
real(wp), intent(in) :: at(nta,nab),taam(nab),slso(nto)
real(wp), intent(in), dimension(nto,n) :: hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc
real(wp),intent(out) :: aoif_at(nob),aoif_pa(nob),aoif_ar,aoif_so
real(wp),intent(out) :: co2a(nob),co3a(nob),omea(nob),co2p(nob),co3p(nob),omep(nob),co2ar,co3ar,omear,co2so,co3so,omeso


real(wp) :: to_a(nob),to_p(nob),to_ar,to_so,pta_nh(3),pta_sh(3)
real(wp) :: ase_ar,ase_at(nob),ase_pa(nob),ase_so
real(wp),dimension(nob) :: pco2wa,k0a,hco3a,pco2wp,k0p,hco3p
real(wp) :: pco2war,k0ar,hco3ar,pco2wso,k0so,hco3so,&
            pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh
real(wp) :: rcarat(nob),rcarpa(nob),rcarar,rcarso,wcar(nab),wsil(nab),worg(nab),vol(nab)
real(wp), intent(out) :: rorgat(nob),rorgpa(nob),rorgar,rorgso


to_a  = (/hlna(1,1),mlna(1,1),llna(1,1),llsa(1,1),mlsa(1,1)/)
to_p  = (/hlnp(1,1),mlnp(1,1),llnp(1,1),llsp(1,1),mlsp(1,1)/)
to_ar = arc(1,1)
to_so = soc(1,1)

call get_poly_pta(at(1,1),at(1,2),at(1,3),at(1,4),at(1,5),at(1,6),pta_nh,pta_sh)
call atm_airsea_heat_exc(at(10,1),at(10,2),at(10,6),pta_nh,pta_sh,&
                            at(1,2),at(1,3),at(1,4),at(1,5),&
                            to_ar,to_a(1),to_a(2),to_a(3),to_a(4),to_a(5),to_p(1),to_p(2),to_p(3),to_p(4),to_p(5),to_so,&
                            ase_ar,ase_at,ase_pa,ase_so,aoif_at,aoif_pa,aoif_ar,aoif_so)

call oc_carbon_sys_surface(hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,soc,arc,slso,&
                           pco2wa ,k0a ,co2a ,co3a ,hco3a ,omea ,&
                           pco2wp ,k0p ,co2p ,co3p ,hco3p ,omep ,&
                           pco2war,k0ar,co2ar,co3ar,hco3ar,omear,&
                           pco2wso,k0so,co2so,co3so,hco3so,omeso,&
                           pco2wsh,k0sh,co2sh,co3sh,hco3sh,omesh)


call external_forcing(taam,rcarat,rcarpa,rcarar,rcarso,rorgat,rorgpa,rorgar,rorgso,wcar,wsil,worg,vol)


end subroutine pre_sed
