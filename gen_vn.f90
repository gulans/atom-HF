subroutine gen_vn(Ngrid,Z,r,sigma,vn)
use modinteg
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: Z,r(Ngrid),sigma
real(8), intent(out) :: vn(Ngrid)
integer :: ir

real(8), PARAMETER :: Pi = 3.1415926535897932384d0


real(8) ::ksi


  if (sigma.lt.1d-30)then
    vn=-Z/r 
  else
    ksi=3d0/(2d0*sigma**2)
    vn=(-Z/r)*erf(dsqrt(ksi)*r)
!old version
!    vn=(-Z/r)*erf(r/(dsqrt(2d0)*sigma)) 
  endif

  

end subroutine
