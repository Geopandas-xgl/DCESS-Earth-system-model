subroutine external_forcing(taam,&
                            rcarat,rcarpa,rcarar,rcarso,&
                            rorgat,rorgpa,rorgar,rorgso,&
                            wcar,wsil,worg,vol)
! Calculates external forcing like river inputs, weathering rates and volcanic emissions
! Inputs: 
! taam: annual mean atmospheric temperature 
! Outputs:
! River input of carbonate and organic carbon rcar, rorg (m/s)
! Weathering rates of carbonate, silicate and organic carbon wcar,wsil,worg (mol/s)
! Volcanic emissions: vol (mol/s)

use parameters, only: wp,tamnpi,taenpi,taespi,tamspi,q10,friva,frivp,frivar,frivso,&
bcarpi,borgpi,bcorgpi,gamma_sil,volpi,fvol,nab,nob

implicit none
real(wp), intent(in) :: taam(nab)
real(wp), intent(out) :: rcarat(nob),rcarpa(nob),rcarar,rcarso,rorgat(nob),rorgpa(nob),rorgar,rorgso,&
 wcar(nab),wsil(nab),worg(nab),vol(nab)
real(wp) :: tamn,taen,taes,tams,dtahn,dtamn,dtaen,dtaes,dtams,dtahs,&
lmdhn,lmdmn,lmden,lmdes,lmdms,lmdhs,lmd(nab),wcarg,wsilg,worgg,&
frhn,frmn,fren,fres,frms,frhs,friv(nab)

tamn = taam(2)
taen = taam(3)
taes = taam(4)
tams = taam(5)

! dta = 0 --> PI steady state
dtahn = 0._wp ! tamn - tamnpi
dtamn = 0._wp ! tamn - tamnpi
dtaen = 0._wp ! taen - taenpi
dtaes = 0._wp ! taes - taespi
dtams = 0._wp ! tams - tamspi
dtahs = 0._wp ! tams - tamspi

lmdhn = q10**(dtahn/10._wp)
lmdmn = q10**(dtamn/10._wp)
lmden = q10**(dtaen/10._wp)
lmdes = q10**(dtaes/10._wp)
lmdms = q10**(dtams/10._wp)
lmdhs = q10**(dtahs/10._wp)

lmd = (/lmdhn,lmdmn,lmden,lmdes,lmdms,lmdhs/)
!-----------------------------------------------------
! River inputs
!-----------------------------------------------------
! Carbonate
rcarat = friva*lmd(1:5)*bcarpi
rcarpa = frivp*lmd(1:5)*bcarpi
rcarar = frivar*bcarpi*lmdhn
rcarso = frivso*bcarpi*lmdhs

! Phosphorus
rorgat = friva*lmd(1:5)*borgpi
rorgpa = frivp*lmd(1:5)*borgpi
rorgar = frivar*borgpi*lmdhn
rorgso = frivso*borgpi*lmdhs
!-----------------------------------------------------

!-----------------------------------------------------
! Global atm source/sink rates
!-----------------------------------------------------
wcarg = bcarpi/(1+gamma_sil)
wsilg = gamma_sil*wcarg
worgg = (bcorgpi + gamma_sil*bcarpi/(1._wp+gamma_sil) - volpi)

frhn = friva(1) + frivp(1) + frivar
frmn = friva(2) + frivp(2)
fren = friva(3) + frivp(3)
fres = friva(4) + frivp(4)
frms = friva(5) + frivp(5)
frhs = frivso
friv = (/frhn,frmn,fren,fres,frms,frhs/)

wcar = wcarg*friv*lmd
wsil = wsilg*friv*lmd
worg = worgg*friv*lmd

! Volcanic flux
vol = volpi*fvol

end subroutine external_forcing