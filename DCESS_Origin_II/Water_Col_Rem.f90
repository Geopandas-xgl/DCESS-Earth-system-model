subroutine water_col_rem(pon,poc,pca,np,aoif,ga,porg,pcal,rp,rmorg,rmorgc,rmcar)
! Water column remineralization
! Set RM to the divergence of the flux profile for each layer
! Account for bathymetri (GA [0-1]) and apply 
! dissolution/remineralization percentages from the 
! sediment model
! Inputs:
! pon,poc,pca: particulate material profile for nutrient (pon), organic carbon (poc) and calcite (pca) 
! np: ocean new production
! aoif: ice-free ocean area
! ga: ocea area fraction profile
! porg: percentage organic carbon rem (from sediment)
! pcal: percentage organic calcite dissolution
! rp: rain ratio
! Outputs:
! rmorg,rmorgc,rmcar: remineralization profiles of nutrient (rmorg), organic carbon (rmorgc) and calcite (rmcar)

use parameters, only: wp,n,rcp

implicit none
real(wp), dimension(n), intent(in) :: pon,poc,pca,ga,porg,pcal
real(wp), intent(in) :: np,aoif,rp
real(wp), dimension(n), intent(out) :: rmorg,rmorgc,rmcar

rmorg = aoif*np*(/ -pon(1),&
                    (pon(1:n-2)-pon(2:n-1))*ga(2:n-1) + pon(1:n-2)*(ga(1:n-2)-ga(2:n-1))*porg(1:n-2),&
                    (pon(n-1)  -pon(n)    )*ga(n)     + pon(n-1)  *(ga(n-1)  -ga(n)    )*porg(n-1) + &
                    pon(n)*ga(n)*porg(n)/) ! mol P/s

rmorgc = aoif*np*(/ -poc(1),&
                        (poc(1:n-2)-poc(2:n-1))*ga(2:n-1) + poc(1:n-2)*(ga(1:n-2)-ga(2:n-1))*porg(1:n-2),&
                        (poc(n-1)  -poc(n)    )*ga(n)     + poc(n-1)  *(ga(n-1)  -ga(n)    )*porg(n-1) + &
                        poc(n)*ga(n)*porg(n)/) ! mol orgC/s

rmcar = aoif*np*rcp*rp*(/ -pca(1),&
                            (pca(1:n-2)-pca(2:n-1))*ga(2:n-1) + pca(1:n-2)*(ga(1:n-2)-ga(2:n-1))*pcal(1:n-2),&
                            (pca(n-1)  -pca(n))*ga(n)+pca(n-1)*(ga(n-1)-ga(n))*pcal(n-1)+ &
                            pca(n)*ga(n)*pcal(n)/) ! mol Carbonate/s

end subroutine water_col_rem