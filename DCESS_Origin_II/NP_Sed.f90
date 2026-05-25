real(wp) function np_sed(sr,rin,ao) result(np)
! Calculates ocena new production for the sediment model
! Inputs:
! sr: ocean source/sinks matrix
! rin: river input or organic matter
! ao: ocean area
! Outputs:
! np: new production

use parameters, only: wp,nto,n
implicit none
real(wp), intent(in) :: sr(nto,n),rin,ao

np = (sr(3,1)-rin)/ao

end function np_sed