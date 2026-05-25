subroutine carsys_pres(psi,co3,co3s,co2,hco3,hp)
!  Calculates solubilities and iterate the pressure dependent carbonate system 
!  Input : psi  - vertical tracer distribution
!                see OdeExp.f90 for the data structure.
!                The sub-set used in CarSysPres are
!                      1:   - temperature [oC]
!                      2:   - salinity []
!                      5:   - DIC [mol/m3] 
!                      8:   - alkalinity [eq/m3]
!  Output variables
!       CO2        [mol/m3]
!       CO3        [mol/m3]
!       CO3s       [mol/m3]
!       HCO3       [mol/m3]
!       pH         []


use parameters, only: wp,nto,n,rho0,dm,d,ca0

implicit none
real(wp), dimension(nto,n), intent(in) :: psi
real(wp), dimension(n), intent(out) :: co3,co3s,co2,hco3,hp
real(wp) :: ro,dic(n),alk(n),tc(n),t(n),s(n)
real(wp) :: k0x,k1x,k2x,kwx,kbx,kcx,r
real(wp), dimension(n) :: k1,k2,kw,kb,kc,br,p,molarvol,disscompr,ratio,kw_p,k1_p,k2_p,kb_p,kc_p
real(wp), dimension(n) :: a,kr,z,ca
integer :: i,imx

! Convert input DIC and ALK to [mol/kg] from [mol/m3] using a fixed density
ro   = rho0             ! [kg/m3]  
dic  = psi(5,:)/ro      ! [mol/kg]
alk  = psi(8,:)/ro      ! [mol/kg]

! Convert temperature to kelvin scale
tc = psi(1,:) 
t = tc + 273.15_wp            ! [K]
s = psi(2,:)

do i=1,n
call carbon_constants(t(i),s(i),k0x,k1x,k2x,kwx,kbx,kcx)
k1(i) = k1x
k2(i) = k2x
kw(i) = kwx
kb(i) = kbx
kc(i) = kcx
end do

! Total boric acid in sea-water, [mol/kg] Millero 1982, see Millero 1995 below eq 53
br   = 0.000416_wp*s/(35._wp)

! Pressure dependence of system thermodynamics
r = 83.131_wp                         ! [bar cm3 mol-1 K-1] Gas constant
p = 1._wp + (/dm/2._wp, (real(i,wp)*d - d/2._wp + dm, i=1,n-1)/)/10._wp     ! [bar] Pressure

! Molar volume and compressibility changes of dissociation reactions
! Kw: H2O
molarvol    =  (-2.00E01_wp) + ( 1.120E-01_wp)*tc + (-1.410E-03_wp)*(tc**2)
disscompr   =  (-5.13E00_wp) + ( 7.940E-02_wp)*tc
ratio       = -(molarvol/(r*t))*p+(.5E-03_wp*disscompr/(r*t))*(p**2)
kw_p        =  kw*exp(ratio)
! print*, molarvol( (/1,15,27,39,45,55/))
! print*, disscompr( (/1,15,27,39,45,55/))
! print*, ratio( (/1,15,27,39,45,55/))


!K1: H2CO3
molarvol   =  (-2.550E+01_wp) + ( 1.271E-01_wp)*tc + (0.E-03_wp)*(tc**2)
disscompr  =  (-3.080E+00_wp) + ( 8.770E-02_wp)*tc
ratio      = -(molarvol/(r*t))*p+(.5E-03_wp*disscompr/(r*t))*(p**2)
k1_p       =  k1*exp(ratio)

!K2: HCO3(-)
molarvol  =  (-1.582E+01_wp) + (-2.190E-02_wp)*tc + (    0E-03_wp)*(tc**2) 
disscompr =  ( 1.130E+00_wp) + (-1.475E-01_wp)*tc
ratio     = -(molarvol/(r*t))*p+(.5E-03_wp*disscompr/(r*t))*(p**2)
k2_p      =  k2*exp(ratio)

!Kb: B(OH)3
molarvol   =  (-2.948E+01_wp) + ( 1.622E-01_wp)*tc + ( 2.608E-03_wp)*(tc**2)
disscompr  =  (-2.840E+00_wp) + (     0E-02_wp)*tc
ratio      = -(molarvol/(r*t))*p+(.5E-03_wp*disscompr/(r*t))*(p**2)
kb_p       =  kb*exp(ratio)

!Kc: CaCO3(cal)
molarvol  =  (-4.876E+01_wp) + ( 5.304E-01_wp)*tc + (     0E-03_wp)*(tc**2) 
disscompr =  (-1.176E+01_wp) + ( 3.692E-01_wp)*tc
ratio     = -(molarvol/(r*t))*p+(.5E-03_wp*disscompr/(r*t))*(p**2)
kc_p      = kc*exp(ratio)

! Use a salinity corrected value for the concentration of Ca(2+) 
! [Millero 1982] to  estimate the saturation concentration of CO3(2-) 
co3s = kc_p/(ca0*s/35._wp)


! Iterate the carbon dioxide in solution CO_2 using the recursive formulation 
! of Antonie and Morel 1995. The values for the hydrogen ion activity, H+, 
! and the carbonate alkalinity are also computed, CA = HCO3- + 2 CO3--
kr   = k1_p/k2_p
a    = alk
imx  = 500

do i=1,imx
    z    = ( (dic*kr)**2+a*kr*(2._wp*dic-a)*(4._wp-kr) )**(.5)
    co2  = dic-a+(a*kr-dic*kr-4._wp*a+z)/(2._wp*(kr-4._wp))
    hp   = co2*k1_p/(2._wp*a) + (( ((co2*k1_p)**2+8._wp*a*co2*k1_p*k2_p) )**(.5))/(2._wp*a)
    ca   = alk-br*kb_p/(kb_p+hp)-kw_p/hp+hp
    if ( maxval(abs(ca-a)) .LT. 1.e-8_wp ) then
        exit
    end if
    a    = ca
end do

if ( i .GT. imx) then
    print*, '!! Carbonate system not converging, CarSys_press !!'
end if


co3  = k1_p*k2_p*(  dic/( (1._wp+k1_p/hp+k1_p*k2_p/hp**2) )  )/hp**2 ! [mol/kg]

! Convert output values to model units

co3  = co3 *ro                          ! [mol/m3] 
co3s = co3s*ro                          ! [mol/m3] 
co2  = co2 *ro                          ! [mol/m3] 
hco3 = dic*ro-co3-co2                   ! [mol/m3] 


end subroutine carsys_pres