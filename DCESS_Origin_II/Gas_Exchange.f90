subroutine gas_exchange(at,arc,soc,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,slso,&
                        pco2wa,pco2wp,pco2war,pco2wso,pco2wsh,k0a,k0p,k0ar,k0so,k0sh,&
                        co3a,co3p,co3ar,co3so,co3sh,aoa,aop,aoar,aoso,&
                        ashna,asmna,asena,asesa,asmsa,ashnp,asmnp,asenp,asesp,asmsp,&
                        asar,asso,assh,co2_atm,o2_atm)
! Calculates air-sea gas exchange for 12,13,14co2 and oxygen
! Inputs: 
! at: atmospheric tracers
! Ocean tracers: arc,soc,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,slso
! pco2w: water co2
! k0: solubility of co2
! co3: carbonate ion
! ao: ocean areas
! Outputs:
! as: air-sea gas exchange array
! co2, o2 _atm: atmospheric co2 and o2 content in ocean model units

use parameters, only: wp,nto,n,nta,ao_sh,to_f,nab,nob
implicit none
real(wp), dimension(nto,n), intent(in) :: arc,soc,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp
real(wp), intent(in) :: at(nta,nab),slso(nto),aoar,aoso,pco2war,pco2wso,pco2wsh,k0ar,k0so,k0sh,co3ar,co3so,co3sh
real(wp), dimension(nob), intent(in) :: aoa,aop,pco2wa,pco2wp,k0a,k0p,co3a,co3p
real(wp), dimension(9), intent(out) :: ashna,asmna,asena,asesa,asmsa,&
                                       ashnp,asmnp,asenp,asesp,asmsp,&
                                       asar,asso,assh
real(wp), intent(out) :: co2_atm(3),o2_atm
real(wp) :: aoifsh
real(wp), parameter ::  wsar = 5.6627_wp,&          ! wind speed (m/s)
                        wshn = 7.5164_wp,&
                        wsmn = 7.7432_wp,&
                        wsen = 5.9326_wp,&
                        wses = 6.1087_wp,&
                        wsms = 9.1370_wp,&
                        wsso = 7.0378_wp,&
                        mvi  = 22.41361_wp         ! Molar volume for an ideal gas [liter/mol]


call gas_exc(arc(:,1),at(:,1),wsar,aoar,pco2war,k0ar,co3ar,asar)
call gas_exc(hna(:,1),at(:,1),wshn,aoa(1),pco2wa(1),k0a(1),co3a(1),ashna)
call gas_exc(mna(:,1),at(:,2),wsmn,aoa(2),pco2wa(2),k0a(2),co3a(2),asmna)
call gas_exc(ena(:,1),at(:,3),wsen,aoa(3),pco2wa(3),k0a(3),co3a(3),asena)
call gas_exc(esa(:,1),at(:,4),wses,aoa(4),pco2wa(4),k0a(4),co3a(4),asesa)
call gas_exc(msa(:,1),at(:,5),wsms,aoa(5),pco2wa(5),k0a(5),co3a(5),asmsa)

call gas_exc(hnp(:,1),at(:,1),wshn,aop(1),pco2wp(1),k0p(1),co3p(1),ashnp)
call gas_exc(mnp(:,1),at(:,2),wsmn,aop(2),pco2wp(2),k0p(2),co3p(2),asmnp)
call gas_exc(enp(:,1),at(:,3),wsen,aop(3),pco2wp(3),k0p(3),co3p(3),asenp)
call gas_exc(esp(:,1),at(:,4),wses,aop(4),pco2wp(4),k0p(4),co3p(4),asesp)
call gas_exc(msp(:,1),at(:,5),wsms,aop(5),pco2wp(5),k0p(5),co3p(5),asmsp)

call gas_exc(soc(:,1),at(:,6),wsso,aoso,pco2wso,k0so,co3so,asso)

aoifsh = 0.2_wp*ao_sh
call gas_exc_sh(slso,at(:,6),wsso,aoifsh,pco2wsh,k0sh,co3sh,assh,co2_atm,o2_atm)

contains


subroutine gas_exc(psi_s,at,ws,ao,pco2w,k0,co3,as)
! Calculates the air-sea gas exchange of 12,13,14CO2 and O2 for a single zone
implicit none
real(wp), intent(in) :: psi_s(nto),at(nta),ws,ao,pco2w,k0,co3
real(wp), intent(out) :: as(9)
real(wp) :: f13as,f13sa,f14as,f14sa,sc,kw,tk,beta,k0o2

! Carbon ion fractionation factors for air-sea exchange
call frac_factors(psi_s(1),psi_s(5),co3,f13as,f13sa,f14as,f14sa)

!-----------------------------------------------------------------------------
! Air sea 12,13,14C02 exchange
!-----------------------------------------------------------------------------
! Gas exchange coefficient kw (m/s), formulated in the long-term averaged 
! winds U (m/s) (Wanninkhof 1992) and temperature dependent Schmidt number 
! for CO2 (Wanninkhof 1992).
sc = (2073.1_wp - 125.62_wp*psi_s(1) + 3.6276_wp*psi_s(1)**2 - 0.043219_wp*psi_s(1)**3) ! T in oC
kw = gas_exc_coeff(sc,ws)

as(5) = ao*kw*k0*(                   at(5) -                         pco2w ) ! [mol/s]
as(6) = ao*kw*k0*( f13as*at(6)/at(5)*at(5) - f13sa*psi_s(6)/psi_s(5)*pco2w )
as(7) = ao*kw*k0*( f14as*at(7)/at(5)*at(5) - f14sa*psi_s(7)/psi_s(5)*pco2w )


!-----------------------------------------------------------------------------
! Air sea 02 exchange
! -----------------------------------------------------------------------------
! Gas exchange coefficient kw (m/s), formulated in the long-term averaged 
! winds U (m/s) (Wanninkhof 1992) and temperature dependent Schmidt number 
! for O2 (Keeling et al 1998).
sc = 1638._wp - 81.83_wp*psi_s(1) + 1.483_wp*psi_s(1)**2 - 0.008004_wp*psi_s(1)**3 ! T in oC
kw = gas_exc_coeff(sc,ws)

! Bunsen solubility coefficient for oxygen (Weiss 1970)
tk   = psi_s(1) + 273.15_wp             ! Kelvin temperature
beta = exp( -58.3877_wp + 85.8079_wp*(100._wp/tk) + 23.8439_wp*log(tk/100._wp) + &
            psi_s(2)*(-0.034892_wp + 0.015568_wp*(tk/100._wp) - 0.0019387_wp*(tk/100._wp)**2) )

! Conversion of Bunsen coefficient [atm-1] to solubility K0 [mol/m3/atm]. 
! The conversion ignores the the correction for non-ideality, ok for oxygen.

k0o2  = beta/mvi*1.e3_wp     ! [mol/m3/atm]
as(9) = ao*kw*(k0o2*at(9)-psi_s(9))  ! [mol/s]

end subroutine gas_exc

subroutine gas_exc_sh(slso,at,ws,ao,pco2w,k0,co3,&
                        as,co2_atm,o2_atm)
implicit none
real(wp), intent(in) :: slso(nto),at(nta),ws,ao,pco2w,k0,co3
real(wp), intent(out) :: as(9),co2_atm(3),o2_atm
real(wp) :: f13as,f13sa,f14as,f14sa,sc,kw,c12_atm,c13_atm,c14_atm,&
tk,beta,k0o2

! Carbon ion fractionation factors for air-sea exchange
call frac_factors(to_f,slso(5),co3,f13as,f13sa,f14as,f14sa)

!-----------------------------------------------------------------------------
! Air sea 12,13,14C02 exchange
!-----------------------------------------------------------------------------
sc = (2073.1_wp - 125.62_wp*to_f + 3.6276_wp*to_f**2 - 0.043219_wp*to_f**3)     ! T in oC
kw = gas_exc_coeff(sc,ws)

as(5) = ao*kw*k0*(       at(5)             -                       pco2w ) ! [mol/s]
as(6) = ao*kw*k0*( f13as*at(6)/at(5)*at(5) - f13sa*slso(6)/slso(5)*pco2w )
as(7) = ao*kw*k0*( f14as*at(7)/at(5)*at(5) - f14sa*slso(7)/slso(5)*pco2w )

c12_atm = k0*  at(5)
c13_atm = k0*( f13as*at(6)/at(5)*at(5) )
c14_atm = k0*( f14as*at(7)/at(5)*at(5) )
co2_atm = (/c12_atm,c13_atm,c14_atm/)

! Oxygen
sc   = 1638._wp - 81.83_wp*to_f + 1.483_wp*to_f**2 - 0.008004_wp*to_f**3
kw = gas_exc_coeff(sc,ws)

tk   = to_f+273.15_wp ! Kelvin temperature
beta = exp( -58.3877_wp + 85.8079_wp*(100._wp/tk) + 23.8439_wp*log(tk/100._wp) + &
                slso(2)*(-0.034892_wp + 0.015568_wp*(tk/100._wp) - 0.0019387_wp*(tk/100._wp)**2) )

k0o2  = beta/mvi*1.e3_wp     ! [mol/m3/atm]
as(9) = ao*kw*(k0o2*at(9)-slso(9))  ! [mol/s]

o2_atm = atm2oce_o2(0._wp,at(9))

end subroutine gas_exc_sh

subroutine frac_factors(t,dic,co3,f13as,f13sa,f14as,f14sa)
implicit none
real(wp), intent(in) :: t,dic,co3
real(wp), intent(out) :: f13as,f13sa,f14as,f14sa
real(wp) :: w13w,f13a,f13dg
! Carbon ion fractionation factors for air-sea exchange
! C13
w13w  = 0.99901_wp + 8.7e-6_wp*t
f13a  = 0.99869_wp + 4.9e-6_wp*t
f13dg = 1._wp + ((0.014_wp*co3/dic - 0.107_wp)*t+10.53_wp)/1000._wp     ! From Zhang et al, 1994
f13as = w13w*f13a
f13sa = f13as/f13dg

! C14
f14as = 1._wp - 2._wp*(1._wp-f13as)
f14sa = 1._wp - 2._wp*(1._wp-f13sa)

end subroutine frac_factors

real(wp) function atm2oce_o2(sa,o2atm) result(o2oce)
implicit none
real(wp), intent(in) :: sa,o2atm
real(wp) :: tf,tk,beta,k0

tf   = 0._wp
tk   = tf + 273.15_wp
beta = exp( -58.3877_wp + 85.8079_wp*(100._wp/tk) + 23.8439_wp*log(tk/100._wp) + &
        sa*(-0.034892_wp + 0.015568_wp*(tk/100._wp) - 0.0019387_wp*(tk/100._wp)**2) )
k0   = beta/mvi*1.e3_wp     ! [mol/m3/atm]

o2oce = o2atm*k0  ! mol/m3

end function atm2oce_o2

real(wp) function gas_exc_coeff(sc,ws) result(kw)
implicit none
real(wp), intent(in) :: sc,ws
! Calculate gas exchange coefficient according to Wanninkhof 1992
kw = 0.39_wp/(100._wp*60._wp**2)*ws**2*(sc/660._wp)**(-.5)

end function gas_exc_coeff
end subroutine gas_exchange