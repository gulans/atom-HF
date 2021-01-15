subroutine integ_sph_s38_value(Ngrid,r,f,rez)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),f(Ngrid) 
real(8), intent(out) :: rez
integer :: i,i_interp
real(8) :: hh
real(8) :: x(4),y(4),f0,f1,f2,f3,f1_interp,f2_interp


rez=0d0
i=1
f3=4d0*Pi*r(i)**2d0*f(i) !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  

do i=1,Ngrid-1
  hh=(r(i+1)-r(i))/3d0
  if (i.lt.2) then
    i_interp=1
  else if(i.gt.Ngrid-3) then
    i_interp=Ngrid-3
  else
    i_interp=i-1
  endif

  x=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3)/)
  y=(/f(i_interp),f(i_interp+1),f(i_interp+2),f(i_interp+3)/)
  call interp(4,x,y,r(i)+hh,f1_interp)
  call interp(4,x,y,r(i)+2d0*hh,f2_interp)

  f0=f3 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
  f1=4d0*Pi*(r(i)+1d0*hh)**2d0* f1_interp
  f2=4d0*Pi*(r(i)+2d0*hh)**2d0* f2_interp
  f3=4d0*Pi*r(i+1)**2d0*        f(i+1)

  rez=rez+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo

end subroutine
