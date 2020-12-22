subroutine wfnorm(Ngrid,r,psi)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid) 
real(8), intent(out) :: psi(Ngrid)
integer :: ir
real(8) :: hh,norm

hh=r(2)-r(1) !RUPJI

norm=4d0*Pi*sum(r**2d0*psi**2d0*hh)
psi=psi/sqrt(norm)

end subroutine
