subroutine riv_in(sr,rorg,rcar,co2a)
! Add river input to the source/sink terms (sr)
! Inputs:
! sr: ocean source/sink terms
! rorg: river input of organic carbon 
! rcar: river input of carbonate
! co2a: atmospheric co2
! Outputs:
! sr: ocean source/sink terms plus river inputs

use parameters, only: wp,nto,n,gamma_sil,r13pdb,d13ccarpi,rda
implicit none
real(wp), dimension(nto,n), intent(inout) :: sr
real(wp), intent(in) :: rorg,rcar,co2a(3)


! PO4
sr(3,1) = sr(3,1) + rorg

! DIC-12
sr(5,1) = sr(5,1) + 2._wp*rcar

! DIC-13
sr(6,1) = sr(6,1) + rcar/(1._wp+gamma_sil)*((r13pdb*(d13ccarpi+1)) + (1._wp+2._wp*gamma_sil)*co2a(2)/co2a(1))

! DIC-14
sr(7,1) = sr(7,1) + rcar/(1._wp+gamma_sil)*(                         (1._wp+2._wp*gamma_sil)*co2a(3)/co2a(1))

! ALK
sr(8,1) = sr(8,1) +(2._wp*rcar + rorg*rda)
end subroutine riv_in
