subroutine wfnorm(Ngrid,r,psi)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid) 
real(8), intent(out) :: psi(Ngrid)
integer :: ir
real(8) :: hh,norm

norm=0
do ir=1, Ngrid-1
  hh=r(ir+1)-r(ir)
  norm=norm+4d0*Pi*r(ir)**2d0*psi(ir)**2d0*hh
enddo
psi=psi/sqrt(norm)

end subroutine
