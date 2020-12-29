subroutine getvhtrapez(Ngrid,r,rho,vh)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),rho(Ngrid) 
real(8), intent(out) :: vh(Ngrid)
integer :: i
real(8) :: hh,f0,f1
real(8) :: integ1(Ngrid), integ2(Ngrid)




hh=(r(2)-r(1))/2 !RUPJI

i=1
f0=4d0*Pi*r(i)**2*rho(i)*hh
integ1(i)=f0

i=2
f1=4d0*Pi*r(i)**2*rho(i)*hh
integ1(i)=integ1(i-1)+f1

do i=3,Ngrid
  f0=f1
  f1=4d0*Pi*r(i)**2*rho(i)*hh
  integ1(i)=integ1(i-1)+f0+f1
enddo
integ1=integ1/r


i=Ngrid
f0=4d0*Pi*r(i)*rho(i)*hh
integ2(i)=f0

i=Ngrid-1
f1=4d0*Pi*r(i)*rho(i)*hh
integ2(i)=integ2(i+1)+f1


do i=Ngrid-2,1,-1
  f0=f1
  f1=4d0*Pi*r(i)*rho(i)*hh
  integ2(i)=integ2(i+1)+f0+f1
enddo

vh=integ1+integ2

end subroutine
