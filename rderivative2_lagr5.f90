subroutine rderivative2_lagr5(Ngrid,r,fin,rez)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),fin(Ngrid) 
real(8), intent(out) :: rez(Ngrid)
integer :: i,i_interp
real(8) :: x(6),y(6) 



do i=1,Ngrid

  if (i.lt.3) then
    i_interp=1
  else if(i.gt.Ngrid-4) then
    i_interp=Ngrid-5
  else
    i_interp=i-2
  endif

  x=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3),r(i_interp+4),r(i_interp+5)/)
  y=(/fin(i_interp),fin(i_interp+1),fin(i_interp+2),fin(i_interp+3),fin(i_interp+4),fin(i_interp+5)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call lagr_interp_d2(r(i),x,y,rez(i),6)
enddo
end subroutine




