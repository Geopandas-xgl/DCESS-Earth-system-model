subroutine save_1d(name_file,name_var,var,l,dim)
use parameters, only: wp
use netcdf

implicit none
CHARACTER(LEN = 100), intent(in) :: name_file,name_var
integer, intent(in) :: l,dim
real(wp), dimension(l), intent(in) :: var
integer :: ncid,varid


CALL CHECK( NF90_OPEN(name_file,NF90_WRITE,ncid))
CALL CHECK( NF90_REDEF(ncid))
CALL CHECK( NF90_DEF_VAR(ncid, name_var,NF90_DOUBLE, (/dim/)        , varid)  )
CALL CHECK( NF90_ENDDEF(ncid) )
CALL CHECK( NF90_PUT_VAR(ncid, varid , var      ) )
CALL CHECK( NF90_CLOSE(ncid))

end subroutine save_1d