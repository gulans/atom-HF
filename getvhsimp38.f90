subroutine getvhsimp38(Ngrid,r,rho,vh)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),rho(Ngrid) 
real(8), intent(out) :: vh(Ngrid)
integer :: i,i_interp
real(8) :: hh,f0,f1,f2,f3,rho1,rho2
real(8) :: integ1(Ngrid), integ2(Ngrid)
real(8) :: x(4),y(4) 



i=1
f3=4d0*Pi*r(i)**2*rho(i) !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
integ1(i)=0

do i=1,Ngrid-1
  hh=(r(i+1)-r(i))/3
  if (i.lt.2) then
    i_interp=1
  else if(i.gt.Ngrid-3) then
    i_interp=Ngrid-3
  else
    i_interp=i-1
  endif

  x=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3)/)
  y=(/rho(i_interp),rho(i_interp+1),rho(i_interp+2),rho(i_interp+3)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call interp(4,x,y,r(i)+hh,rho1)
  call interp(4,x,y,r(i)+2d0*hh,rho2)

  f0=f3 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
  f1=4d0*Pi*(r(i)+hh)**2d0*rho1
  f2=4d0*Pi*(r(i)+2d0*hh)**2d0*rho2
  f3=4d0*Pi*r(i+1)**2d0*rho(i+1)

  integ1(i+1)=integ1(i)+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo
integ1=integ1/r


i=Ngrid
f0=4d0*Pi*r(i)*rho(i)
integ2(i)=0


do i=Ngrid-1,1,-1
  hh=(r(i+1)-r(i))/3
  if (i.lt.2) then
    i_interp=1
  else if(i.gt.Ngrid-3) then
    i_interp=Ngrid-3
  else
    i_interp=i-1
  endif

  x=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3)/)
  y=(/rho(i_interp),rho(i_interp+1),rho(i_interp+2),rho(i_interp+3)/)

  !interpolation is done twice for the same points, for saving time i could create an array to store the results from the first time  
  call interp(4,x,y,r(i)+hh,rho1)
  call interp(4,x,y,r(i)+2d0*hh,rho2)


  f3=f0 !every iteration f3 is previous iterations f0 (becouse of the oposit intergration direction)  
  f0=4d0*Pi*(r(i))*rho(i)
  f1=4d0*Pi*(r(i)+hh)*rho1
  f2=4d0*Pi*(r(i)+2d0*hh)*rho2


  integ2(i)=integ2(i+1)+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo

vh=integ1+integ2

end subroutine

