subroutine get_energy_ext(Ngrid,r,rho,Z,e_ext)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0

integer, intent(in)  :: Ngrid
real(8), intent(in)  :: r(Ngrid), rho(Ngrid),Z
real(8), intent(out) :: e_ext

!  call integ_sph_s38_value(Ngrid,r,-rho*Z/r,e_ext)
integer :: i,i_interp
real(8) :: hh,rez
real(8) :: x(4),y(4),f0,f1,f2,f3,rho_interp1,rho_interp2


rez=0d0
i=1
f3=-r(i)*rho(i)*Z !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  

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
  y=(/rho(i_interp),rho(i_interp+1),rho(i_interp+2),rho(i_interp+3)/)
  call interp(4,x,y,r(i)+hh,rho_interp1)
  call interp(4,x,y,r(i)+2d0*hh,rho_interp2)

  f0=f3 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
  f1=-(r(i)+1d0*hh)*Z * rho_interp1
  f2=-(r(i)+2d0*hh)*Z * rho_interp2
  f3=-r(i+1)*       Z * rho(i+1)

  rez=rez+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo

e_ext=rez

end subroutine
