subroutine gen_vn(Ngrid,Z,r,sigma,vn)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: Z,r(Ngrid),sigma
real(8), intent(out) :: vn(Ngrid)
integer :: ir

  if (sigma.lt.1d-30)then
    vn=-Z/r 
  else
    vn=(-Z/r)*erf(r/(dsqrt(2d0)*sigma)) 
  endif
  
!do ir=1, 60
!  write(*,*)r(ir),vn(ir), -Z/r(ir)
!enddo
end subroutine
