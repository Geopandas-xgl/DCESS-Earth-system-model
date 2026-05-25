subroutine carbon_constants(t,s,k0,k1,k2,kw,kb,kc)
! Calculates thermodynamics constants for ocean carbon system

use parameters, only: wp
implicit none
real(wp), intent(in) :: t,s
real(wp),intent(out) :: k0,k1,k2,kw,kb,kc

! Solubility of CO2 [mol/kg/atm] given by Weiss 1974, see Millero 1995 eq. 26
k0   =  exp(- 60.2409_wp + 93.4517_wp*(100._wp/t) + 23.3585_wp*log(t/100._wp) + &
            s*( 0.023517_wp  - 0.023656_wp*(t/100._wp) + 0.0047036_wp*(t/100._wp)**2))

! Thermodynamic constants for the dissociation of carbonic acid, 
! k1 and k2 [mol/kg] (Mehrbach et al. 1973), see Millero 95 eq 35/36
! The Mehrbach constants are recomended by Lee and Millero 1997 
k1   = 10**(-(3670.7_wp/t - 62.008_wp + 9.7944_wp*log(t) - 0.0118_wp*s+0.000116_wp*s**2))
k2   = 10**(-(1394.7_wp/t + 4.777_wp                     - 0.0184_wp*s+0.000118_wp*s**2))

! Dissociation constant for sea-water [mol/kg], Millero 95 (eq 63) 
kw   = exp(148.9802_wp - 13847.26_wp/t - 23.6521_wp*log(t) + &
           (-5.977_wp + 118.67_wp/t + 1.0495_wp*log(t))*s**(.5) - 0.01615_wp*s)

! Dissociation constant for boric acid [mol/kg], Dickson 1990, recommended by Millero 95 (eq 52) 
kb   = exp((-8966.90_wp - 2890.51_wp*s**(.5) - 77.942_wp*s  + 1.726_wp*s**(1.5) - 0.0993_wp*s**2)/t + &
            148.0248_wp + 137.194_wp*s**(.5) + 1.62247_wp*s + &
            (-24.4344_wp - 25.085_wp*s**(.5) - 0.2474_wp*s)*log(t) + 0.053105_wp*s**(.5)*t)

! Dissociation constant for calcite [mol/kg], Mucci 1983
kc   = 10**(-171.9065_wp - 0.077993_wp*t + 2839.319_wp/t + 71.595_wp*log10(t) + &
            (-0.77712_wp + 0.0028426_wp*t+178.34_wp/t)*s**(0.5) - 0.07711_wp*s + 0.0041249_wp*s**(1.5))

end subroutine carbon_constants