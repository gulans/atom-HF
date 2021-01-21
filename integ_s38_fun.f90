subroutine integ_s38_fun(Ngrid,r,fin,dir,rez)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid,dir
real(8), intent(in) :: r(Ngrid),fin(Ngrid) 
real(8), intent(out) :: rez(Ngrid)
integer :: i,i_interp
real(8) :: hh,f0,f1,f2,f3,fin1,fin2
real(8) :: x(4),y(4) 


if (dir.eq.1) then
i=1
f3=fin(i) !in numerical recepies (4.1.5) every steps f0 is previous steps f3  
rez(i)=0

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
  y=(/fin(i_interp),fin(i_interp+1),fin(i_interp+2),fin(i_interp+3)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call interp(4,x,y,r(i)+hh,fin1)
  call interp(4,x,y,r(i)+2d0*hh,fin2)

  f0=f3 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
  f1=fin1
  f2=fin2
  f3=fin(i+1)

  rez(i+1)=rez(i)+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo

elseif (dir.eq.-1) then

i=Ngrid
f0=fin(i)
rez(i)=0

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
  y=(/fin(i_interp),fin(i_interp+1),fin(i_interp+2),fin(i_interp+3)/)

  !interpolation is done twice for the same points, for saving time i could create an array to store the results from the first time  
  call interp(4,x,y,r(i)+hh,fin1)
  call interp(4,x,y,r(i)+2d0*hh,fin2)


  f3=f0 !every iteration f3 is previous iterations f0 (becouse of the oposit intergration direction)  
  f0=fin(i)
  f1=fin1
  f2=fin2


  rez(i)=rez(i+1)+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo

endif

end subroutine

