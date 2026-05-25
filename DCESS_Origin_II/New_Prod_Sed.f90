subroutine new_prod_sed(srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp,srar,srso,&
                        rorga,rorgp,rorgar,rorgso,npat,nppa,npar,npso)
! Calculates the ocean new production for the sediment model
! Inputs:
! sr: internal source/sink matrix
! rorg: river input of organic carbon
! Outputs:
! np: ocean new production

use parameters, only: wp,nto,n,&
ao_ar,ao_na,ao_mna,ao_ena,ao_esa,ao_msa,ao_np,ao_mnp,ao_enp,ao_esp,ao_msp,ao_so,nob
implicit none

real(wp), dimension(nto,n), intent(in) :: srhna,srmna,srena,sresa,srmsa,srhnp,srmnp,srenp,sresp,srmsp,srar,srso
real(wp), intent(in) :: rorga(nob),rorgp(nob),rorgar,rorgso
real(wp), intent(out) :: npat(nob),nppa(nob),npar,npso
real(wp), external :: np_sed
real(wp) :: nphna,npmna,npena,npesa,npmsa,nphnp,npmnp,npenp,npesp,npmsp

nphna = np_sed(srhna,rorga(1),ao_na )
npmna = np_sed(srmna,rorga(2),ao_mna)
npena = np_sed(srena,rorga(3),ao_ena)
npesa = np_sed(sresa,rorga(4),ao_esa)
npmsa = np_sed(srmsa,rorga(5),ao_msa)

nphnp = np_sed(srhnp,rorgp(1),ao_np )
npmnp = np_sed(srmnp,rorgp(2),ao_mnp)
npenp = np_sed(srenp,rorgp(3),ao_enp)
npesp = np_sed(sresp,rorgp(4),ao_esp)
npmsp = np_sed(srmsp,rorgp(5),ao_msp)

npat = (/nphna,npmna,npena,npesa,npmsa/)
nppa = (/nphnp,npmnp,npenp,npesp,npmsp/)

npar  = np_sed(srar,rorgar,ao_ar)
npso  = np_sed(srso,rorgso,ao_so)

end subroutine new_prod_sed
include 'NP_Sed.f90'