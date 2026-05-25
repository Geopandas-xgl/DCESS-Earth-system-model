subroutine carsys_sur(ti,si,dici,alki,pco2,k0,co2,co3,hco3,ome)
! Calculates solubilities and iterate the carbonate system
! Inputs:
! ti: input temp
! si: input salinity
! dici: input dic
! alki: input alk
! Output variables
!      pCO2       [atm]
!      K0         [mol/m3/atm]
!      CO2        [mol/m3]
!      CO3        [mol/m3]
!      HCO3       [mol/m3]
!      Ome

use parameters, only:wp,rho0,ca0,s0
implicit none
real(wp), intent(in) :: ti,si,dici,alki
real(wp), intent(out) :: pco2,k0,co2,co3,hco3,ome
real(wp) :: ro,dic,alk,t,s,k1,k2,kw,kb,kc
real(wp) :: co3s,br,kr,a,z,hp,ca
integer :: i,imx

! Convert input DIC and ALK to [mol/kg] from [mol/m3] using a fixed density
ro   = rho0    ! [kg/m3]  
dic  = dici/ro  ! [mol/kg]
alk  = alki/ro  ! [mol/kg]

! Convert temperature to kelvin scale
t = ti + 273.15_wp            ! [K]
s = si

call carbon_constants(t,s,k0,k1,k2,kw,kb,kc)
! print*, k0,k1,k2,kw,kb,kc

! Use a salinity corrected value for the concentration of Ca(2+) 
! [Millero 1982] to  estimate the saturation concentration of CO3(2-) 
co3s = kc/(ca0*s/s0)

! Total boric acid in sea-water, [mol/kg] Millero 1982, see Millero 1995 below eq 53
br   = 0.000416_wp*s/(35._wp)

! Iterate the carbon dioxide in solution CO_2 using the recursive formulation 
! of Antonie and Morel 1995. The values for the hydrogen ion activity, H+, 
! and the carbonate alkalinity are also computed, CA = HCO3- + 2 CO3--
kr   = k1/k2
a    = alk
imx  = 500

do i=1,imx
    z    = ( (dic*kr)**2+a*kr*(2._wp*dic-a)*(4._wp-kr) )**(.5)
    co2  = dic-a+(a*kr-dic*kr-4._wp*a+z)/(2._wp*(kr-4._wp))
    hp   = co2*k1/(2._wp*a) + (( ((co2*k1)**2+8._wp*a*co2*k1*k2) )**(.5))/(2._wp*a)
    ca   = alk-br*kb/(kb+hp)-kw/hp+hp
    if (abs(ca-a) < 1.e-8_wp) then
        exit
    end if
    a    = ca
end do
if ( i .GT. imx) then
    print*, '!! Carbonate system not converging, CarSys_sur !!'
end if

pco2 = dic/(k0*(1._wp+k1/hp+k1*k2/hp**2))   ! [atm] or alternatively give as CO2/K0
hco3 = k0*k1   *pco2/hp                     ! [mol/kg] 
co3  = k0*k1*k2*pco2/hp**2                  ! [mol/kg] 
ome  = co3/co3s

! Convert output values to model units
k0   = k0  *ro                       ! [mol/m3/atm] 
co2  = co2 *ro                       ! [mol/m3] 
co3  = co3 *ro                       ! [mol/m3] 
hco3 = hco3*ro                       ! [mol/m3] 

end subroutine carsys_sur

include 'Carbon_Constants.f90'