subroutine rderivative_lagr3(Ngrid,r,fin,rez)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),fin(Ngrid) 
real(8), intent(out) :: rez(Ngrid)
integer :: i,i_interp
real(8) :: x(4),y(4) 



do i=1,Ngrid

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
  call lagr_interp_d(r(i),x,y,rez(i),4)
enddo
end subroutine

subroutine lagr_interp_d(xx,x,y,rez,k)
implicit none
integer, intent(in) :: k
real(8), intent(in) :: x(k),y(k),xx
real(8), intent(out):: rez

real(8) :: l(k),prod
integer :: i,j,m

!Get coeficients
do j=1,k
  l(j)=0d0
  do i=1,k
    if (i.ne.j) then
      prod=1d0
      do m=1,k
         if ((m.ne.i).and.(m.ne.j)) then 
         prod=prod*(xx-x(m))/(x(j)-x(m))
         endif
      enddo
      l(j)=l(j)+prod/(x(j)-x(i))
    endif 
  enddo
enddo
!Use coeficients
rez=0d0
do j=1,k
  rez=rez+y(j)*l(j)
enddo

end subroutine



