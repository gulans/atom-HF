subroutine gen_vn(Ngrid,Z,r,vn)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: Z,r(Ngrid)
real(8), intent(out) :: vn(Ngrid)
integer :: ir
real(8) :: sigma

sigma=1d-5
vn=-Z/r
if(.false.)then
  vn=(-Z/r)*erf(r/(dsqrt(2d0)*sigma))
endif

do ir=1, 60
write(*,*)r(ir),vn(ir), -Z/r(ir)
enddo
end subroutine
