subroutine integ_sph_B6_value(Ngrid,r,finf,rezi)
        
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),finf(Ngrid) 
real(8), intent(out) :: rezi

real(8) :: rez(Ngrid),fin(Ngrid)
integer :: i,i_interp
real(8) :: hh,f0,f1,f2,f3,f4,fin1,fin2,fin3
real(8) :: x(6),y(6) 


fin=r**2*finf

i=1
f4=fin(i) !in numerical recepies (4.1.5) every steps f0 is previous steps f3  
rez(i)=0d0

do i=1,Ngrid-1
  hh=(r(i+1)-r(i))/4d0
  if (i.lt.3) then
    i_interp=1
  else if(i.gt.Ngrid-3) then
    i_interp=Ngrid-5
  else
    i_interp=i-2
  endif

  x=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3),r(i_interp+4),r(i_interp+5)/)
  y=(/fin(i_interp),fin(i_interp+1),fin(i_interp+2),fin(i_interp+3),fin(i_interp+4),fin(i_interp+5)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call interp(6,x,y,r(i)+hh,fin1)
  call interp(6,x,y,r(i)+2d0*hh,fin2)
  call interp(6,x,y,r(i)+3d0*hh,fin3)

  f0=f4 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
  f1=fin1
  f2=fin2
  f3=fin3
  f4=fin(i+1)

!  rez(i+1)=rez(i)+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
  rez(i+1)=rez(i)+hh*(14d0*f0+64d0*f1+24d0*f2+64d0*f3+14d0*f4)/45d0

  enddo

rezi=rez(Ngrid)
end subroutine

