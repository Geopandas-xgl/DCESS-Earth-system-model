subroutine olf_lat_sh(ftrop,fgrass,fboreal,olf_bios,f_tf,f_gs,f_bf)
! Determines land fractions according to the time-dependent vegetation zone limits

use parameters, only: wp,dr,fde,fdh,olfts0,olfgs0,olfbs0,ftrops0,fgrasss0
implicit none
real(wp), intent(in) :: ftrop,fgrass,fboreal,olf_bios(1801,2)
real(wp), intent(out) :: f_tf,f_gs,f_bf
real(wp) :: ftropx,fgrassx,fborealx,lat(1801),olf(1801),olfmt,olfmg,olfmb
real(wp), dimension(:), allocatable :: lat_tf,lat_gs,lat_bf,olf_tf,olf_gs,olf_bf
logical, dimension(1801) :: msk_tf,msk_gs,msk_bf

lat = olf_bios(:,1)
olf = olf_bios(:,2)

ftropx   = ftrop/dr
fgrassx  = fgrass/dr
fborealx = fboreal/dr

msk_tf = ( lat < 0._wp   .AND. lat >= ftropx   )
msk_gs = ( lat < ftropx  .AND. lat >= fgrassx  )
msk_bf = ( lat < fgrassx .AND. lat >= fborealx )

allocate(lat_tf(count(msk_tf)));allocate(olf_tf(count(msk_tf)))
allocate(lat_gs(count(msk_gs)));allocate(olf_gs(count(msk_gs)))
allocate(lat_bf(count(msk_bf)));allocate(olf_bf(count(msk_bf)))

lat_tf = pack(lat,msk_tf); olf_tf = pack(olf,msk_tf)
lat_gs = pack(lat,msk_gs); olf_gs = pack(olf,msk_gs)
lat_bf = pack(lat,msk_bf); olf_bf = pack(olf,msk_bf)

olfmt = sum(cos(lat_tf*dr)*olf_tf)/sum(cos(lat_tf*dr))
olfmg = sum(cos(lat_gs*dr)*olf_gs)/sum(cos(lat_gs*dr))
olfmb = sum(cos(lat_bf*dr)*olf_bf)/sum(cos(lat_bf*dr))

deallocate(lat_tf);deallocate(olf_tf)
deallocate(lat_gs);deallocate(olf_gs)
deallocate(lat_bf);deallocate(olf_bf)

f_tf = (1._wp - olfmt)*(sin(abs(ftrop))  - sin(fde)        )/( (1._wp - olfts0)*(sin(ftrops0)  - sin(fde)     ) )     ! Tropical forest area factor of change
f_gs = (1._wp - olfmg)*(sin(abs(fgrass)) - sin(abs(ftrop)) )/( (1._wp - olfgs0)*(sin(fgrasss0) - sin(ftrops0) ) )     ! Grassland area factor of change
f_bf = (1._wp - olfmb)*(sin(fdh)         - sin(abs(fgrass)))/( (1._wp - olfbs0)*(sin(fdh)      - sin(fgrasss0)) )     ! Boreal forest area factor of change

end subroutine olf_lat_sh