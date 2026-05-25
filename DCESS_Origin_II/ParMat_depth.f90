subroutine parmat_depth(psi,pon,poc,pca)
! Determine particulate material profiles as temperature-dependent (nutrient, carbon) and constant for calcite


use parameters, only: wp,nto,n,knut,kcar,q10,n,zmid,d,lmdcar
! PON: particulate organic nutrient
! POC: particulate organic carbon
! PCA: particulate calcium carbonate
implicit none
real(wp), intent(in) :: psi(nto,n)
real(wp), dimension(n), intent(out) :: pon,poc,pca
real(wp) :: z(n),pnut(n),pcar(n),meant
integer :: i

meant = 9.3195_wp ! global mean (500 m) ocean temp

pnut(1)=1._wp
pcar(1)=1._wp
do i=2,n
    pnut(i)=pnut(i-1)*(1._wp - knut*exp(0.1_wp*log(q10)*(psi(1,i)-meant)))
    pcar(i)=pcar(i-1)*(1._wp - kcar*exp(0.1_wp*log(q10)*(psi(1,i)-meant)))
end do

z   = zmid - d/2._wp
pon = pnut
poc = pcar
pca = exp(-lmdcar *z)


end subroutine parmat_depth