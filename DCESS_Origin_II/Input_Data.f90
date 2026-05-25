subroutine input_data(inso,o_p,shl_dat,olf_bios,gaa,gap,gaar,gaso)
! Load input data used in the model:
! inso: solar insolation
! o_p: ocean pressure
! shl_dat: Antarctic slope data
! olf_bios: ocen to land fraction as function of latitude
! gaa,gap,gaar,gaso: vertical profile of ocean area fraction

use parameters, only: wp,n,nto,nta,nz_all,nob
implicit none
character(len=200) :: header,line
real(wp),intent(out) :: shl_dat(nz_all,3),o_p(n,7),inso(19,14),olf_bios(1801,2),&
                        gaa(nob,n),gap(nob,n),gaar(n),gaso(n)
integer :: i

! Insolation
open(unit=82,file='Input_Data/Insolation.dat',status='old', action='read')
read(82,*) (inso(i,:), i=1,19)
close(82)

! Ocean pressure
open(unit=10, file='Input_Data/Ocean_pressure.dat', status='old', action='read')
read(10, '(A)') header
do i = 1, n
    read(10, '(A)') line
    read(line, *) o_p(i, :)
end do
close(10)

! Shelf data
open(unit=41, file='Input_Data/Shelf_data.dat', status='old', action='read')
read(41, '(A)') header
do i = 1, nz_all
    read(41, '(A)') line
    read(line, *) shl_dat(i, :)
end do
close(41)

! Land biosphere data (olfBios)
open(unit=42, file='Input_Data/olf_Bios.dat', status='old', action='read')
read(42, '(A)') header
do i = 1, 1801
    read(42, '(A)') line
    read(line, *) olf_bios(i, :)
end do
close(42)

! Hypsometric profiles
open(unit=83,file='Input_Data/GAA.dat',status='old', action='read')
read(83,*) (gaa(i,:), i=1,nob)
close(83)

open(unit=83,file='Input_Data/GAP.dat',status='old', action='read')
read(83,*) (gap(i,:), i=1,nob)
close(83)

open(unit=83,file='Input_Data/GAAr.dat',status='old', action='read')
read(83,*) (gaar(i), i=1,n)
close(83)

open(unit=83,file='Input_Data/GASO.dat',status='old', action='read')
read(83,*) (gaso(i), i=1,n)
close(83)

end subroutine input_data