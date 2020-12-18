subroutine gengrid(Ngrid,Rmax,r)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: Rmax 
real(8), intent(out) :: r(Ngrid)
integer :: ir

do ir=1,Ngrid
  r(ir)=dble(ir)/dble(Ngrid)*Rmax
enddo

end subroutine
