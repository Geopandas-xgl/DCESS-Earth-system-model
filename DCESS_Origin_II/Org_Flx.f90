subroutine orgflx(psi,aoif,lf,iq,ga,ome,co2,co3,porg,pcal,sr)
! Calculates biogeochemical internal source/sinks terms for each ocean tracers
! Inputs:
! psi: Ocean tracer (see OdeExp for data distribution)
! aoif: ice-free ocean area (m2)
! lf: limitation factor (-)
! iq: solar radiation (w/m2)
! ga: ocean area distribution
! ome: calcite saturation state
! co2: ocean co2
! co3: carbonate content
! porg: percentage organic carbon rem (from sediment)
! pcal: percentage organic calcite dissolution
! Outputs
! sr: internal biogeochemical source/sinks (tracers/s)

use parameters, only: wp,sy,nto,n,dm,rcp,rcd,rda,rca,rno,rcop,mino2

implicit none
real(wp), intent(in) :: psi(nto,n),lf,iq,ome,co2,co3,aoif,ga(n),porg(n),pcal(n)
real(wp), intent(out) :: sr(nto,n)
real(wp) :: bpe,phlf,ihlf,npp,np,rp,fcp,fnp,f13dg,fco2,f14np,f14co2,f14fcp
real(wp), dimension(n) :: pon,poc,pca,rmorg,rmorgc,rmcar
real(wp), external :: rain_ratio

! New production of organic matter (euphotic zone phosphate uptake)
bpe  = 1._wp/sy                                         ! Bio-Production Efficiency [s-1]
phlf = 1.e-6_wp                                         ! Half-saturation constant [molP/m3], st. val. 1e-6  
ihlf = 100._wp
npp  = bpe*lf*dm*iq/(iq+ihlf)*psi(3,1)**2/(phlf+psi(3,1))
np   = npp

rp = rain_ratio(ome,psi(1,1))

! Carbon ion fractionation for organic and carbonate new production
fcp   = 1._wp
fnp   = ((-(25._wp - (116.96_wp*1000._wp*lf*psi(3,1) + 81.42_wp)/(1000._wp*co2)))/1000._wp+1._wp)   ! Pagani et al, 1999
f13dg = 1._wp + ((0.014_wp*co3/psi(5,1)-0.107_wp)*psi(1,1)+10.53_wp)/1000._wp
fco2  = (0.99869_wp + 4.9e-6_wp *psi(1,1))/f13dg

! Particulate material profiles
call parmat_depth(psi,pon,poc,pca)

! Organic matter remineralization and calcite dissolution
call water_col_rem(pon,poc,pca,np,aoif,ga,porg,pcal,rp,rmorg,rmorgc,rmcar)

sr=0._wp

! Phosphate
sr(3,:) = rmorg

! DIC-12
sr(5,:) = (rmorgc*rcp + rmcar*rcd)

! DIC-13
sr(6,:) = psi(6,1)/psi(5,1)*(rmorgc*rcp*fnp*fco2 + rmcar*rcd*fcp)

! DIC-14
f14np  = 1._wp - 2._wp*(1._wp - fnp )
f14co2 = 1._wp - 2._wp*(1._wp - fco2)
f14fcp = 1._wp - 2._wp*(1._wp - fcp )
sr(7,:) = psi(7,1)/psi(5,1)*(rmorgc*f14np*f14co2*rcp + rmcar*f14fcp*rcd)

! ALK
sr(8,:) = (rmorg*rda + rmcar*rca)

! O2
sr(9,:) = (rmorg*rno + rmorgc*rcop)*merge(1._wp,0._wp,psi(9,:)>mino2)
sr(9,1) = rmorg(1)*(rno+rcop)


end subroutine orgflx

include 'Rain_Ratio.f90'
include 'ParMat_depth.f90'
include 'Water_Col_Rem.f90'