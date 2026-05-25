subroutine atm_airsea_heat_exc(fsi_na,fsi_np,fsi_so,pta_n,pta_s,&
        ta_mn,ta_en,ta_es,ta_ms,&
        to_ar,to_na,to_mna,to_ena,to_esa,to_msa,to_np,to_mnp,to_enp,to_esp,to_msp,to_so,&
        ase_ar,ase_at,ase_pa,ase_so,aoif_at,aoif_pa,aoif_ar,aoif_so)
! Calculates air-sea heat exchange for each atmospheric/ocean box
! Inputs:
! Ice lines for North Atlantic/Arctic, North Pacific and Southern Ocean (fsi_na,fsi_np,fsi_so)
! Terms for atmopsheric temperature profile (pta_n,pta_s)
! Atmospheric temperatures: ta_mn,ta_en,ta_es,ta_ms
! Ocean temperatures for each ocean box: to_ar,to_na,to_mna,to_ena,to_esa,to_msa,to_np,to_mnp,to_enp,to_esp,to_msp,to_so
! Outputs:
! ase_ar: air-sea heat exchange Arctic Ocean
! ase_so: air-sea heat exchange Southern Ocean
! ase_at,ase_pa: air-sea heat exchange for Atlantic and Pacific Ocean [hn mn en es ms]
! Ice-free ocean areas
! aoif_ar,aoif_so,aoif_at,aoif_pa

use parameters, only: wp,fbr,fdh,fdl,b_north,ko,h_area,olf_ar,olf_na,olf_np,olf_mnp,&
ao_na,ao_mna,ao_ena,ao_esa,ao_msa,ao_mnp,ao_enp,ao_esp,ao_msp,dr,nob

implicit none
real(wp), intent(in) :: fsi_na,fsi_np,fsi_so,pta_n(3),pta_s(3)
real(wp), intent(in) :: to_ar,to_na,to_mna,to_ena,to_esa,to_msa,to_np,to_mnp,to_enp,to_esp,to_msp,to_so
real(wp), intent(in) :: ta_mn,ta_en,ta_es,ta_ms
real(wp), intent(out) :: ase_ar,ase_at(nob),ase_pa(nob),ase_so,aoif_at(nob),aoif_pa(nob),aoif_ar,aoif_so
real(wp) :: aoif_na,aoif_mna,aoif_ena,aoif_esa,aoif_msa
real(wp) :: aoif_np,aoif_mnp,aoif_enp,aoif_esp,aoif_msp
real(wp) :: taif_ar,taif_na,taif_mna,taif_ena,taif_esa,taif_msa,taif_np,taif_mnp,taif_enp,taif_esp,taif_msp,taif_so
real(wp) :: d_heat(3)
real(wp), external :: ta_mean

if (fsi_na > fbr .and. fsi_na<= b_north) then
aoif_ar = h_area*(sin(fsi_na)-sin(fbr))*olf_ar
taif_ar = ta_mean(pta_n,fsi_na,fbr)

aoif_na  = ao_na
taif_na  = ta_mean(pta_n,fBr,fdh)
elseif (fsi_na <= fbr .and. fsi_na > fdh) then
aoif_ar = 0._wp
taif_ar = 0._wp

aoif_na = h_area*(sin(fsi_na)-sin(fdh))*olf_na
taif_na = ta_mean(pta_n,fsi_na,fdh)
else !if (fsi_na < fdh) then
print*, 'sea-ice hna out of limits'
end if
aoif_mna = ao_mna
aoif_ena = ao_ena
aoif_esa = ao_esa
aoif_msa = ao_msa

! if (fsi_np > fdh .and. fsi_np <= fbr) then
! aoif_np = h_area*(sin(fsi_np)-sin(fdh))*olf_np
! taif_np = ta_mean(pta_n,fsi_np,fdh)

! aoif_mnp = ao_mnp
! taif_mnp = ta_mn
! elseif (fsi_np <= fdh .and. fsi_np>fdl) then
! aoif_np = 0._wp
! taif_np = 0._wp

! aoif_mnp = h_area*(sin(fsi_np)-sin(fdl))*olf_mnp
! taif_mnp = ta_mean(pta_n,fsi_np,fdl)
! else
!         ! write(*,'(f15.4)') fsi_np/dr
! ! print*, 'sea-ice hnp out of limits'
! end if

aoif_np = h_area*(sin(fsi_np)-sin(fdh))*olf_np
taif_np = ta_mean(pta_n,fsi_np,fdh)

aoif_mnp = ao_mnp
taif_mnp = ta_mean(pta_n,fdh,fdl)

aoif_so = h_area*(sin(fsi_so)-sin(fdh))
taif_so = ta_mean(pta_s,fsi_so,fdh)

aoif_enp = ao_enp
aoif_esp = ao_esp
aoif_msp = ao_msp

taif_mna = ta_mn
taif_ena = ta_en
taif_esa = ta_es
taif_msa = ta_ms

taif_enp = ta_en
taif_esp = ta_es
taif_msp = ta_ms

d_heat = (/40._wp,20._wp,0._wp/)

ase_ar  = aoif_ar *( d_heat(3) + ko*(taif_ar  - to_ar ) )

ase_at(1) = aoif_na *( d_heat(3) + ko*(taif_na  - to_na ) )
ase_at(2) = aoif_mna*( d_heat(2) + ko*(taif_mna - to_mna) )
ase_at(3) = aoif_ena*( d_heat(1) + ko*(taif_ena - to_ena) )
ase_at(4) = aoif_esa*( d_heat(1) + ko*(taif_esa - to_esa) )
ase_at(5) = aoif_msa*( d_heat(2) + ko*(taif_msa - to_msa) )

ase_pa(1) = aoif_np *( d_heat(3) + ko*(taif_np  - to_np ) )
ase_pa(2) = aoif_mnp*( d_heat(2) + ko*(taif_mnp - to_mnp) )
ase_pa(3) = aoif_enp*( d_heat(1) + ko*(taif_enp - to_enp) )
ase_pa(4) = aoif_esp*( d_heat(1) + ko*(taif_esp - to_esp) )
ase_pa(5) = aoif_msp*( d_heat(2) + ko*(taif_msp - to_msp) )

ase_so  = aoif_so *( d_heat(3) + ko*(taif_so - to_so ) )

aoif_at = (/aoif_na,aoif_mna,aoif_ena,aoif_esa,aoif_msa/)
aoif_pa = (/aoif_np,aoif_mnp,aoif_enp,aoif_esp,aoif_msp/)

! write(*,'(f10.4)') fsi_na/dr
! ! write(*,'(1e25.15)') aoif_na! bien
! write(*,'(1e25.15)') taif_na
! ! write(*,'(1e25.15)') to_na !bien

! write(*,'(5e25.15)') ase_at
! write(*,'(5e25.15)') ase_pa
! write(*,'(2e25.15)') ase_ar,ase_so

! print*, 'heat_exchange'

end subroutine atm_airsea_heat_exc

include 'Ta_mean.f90'