subroutine gen_vn(Ngrid,Z,r,vn)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: Z,r(Ngrid)
real(8), intent(out) :: vn(Ngrid)
integer :: ir

vn=-Z/r

end subroutine
