subroutine save_2d(name_file,name_var,var,l1,l2,dim1,dim2)
use parameters, only: wp
use netcdf

implicit none
CHARACTER(LEN = 100), intent(in) :: name_file,name_var
integer, intent(in) :: l1,l2,dim1,dim2
real(wp), dimension(l1,l2), intent(in) :: var
integer :: ncid,varid


CALL CHECK( NF90_OPEN(name_file,NF90_WRITE,ncid))
CALL CHECK( NF90_REDEF(ncid))
CALL CHECK( NF90_DEF_VAR(ncid, name_var,NF90_DOUBLE, (/dim1,dim2/)        , varid)  )
CALL CHECK( NF90_ENDDEF(ncid) )
CALL CHECK( NF90_PUT_VAR(ncid, varid , var      ) )
CALL CHECK( NF90_CLOSE(ncid))

end subroutine save_2d