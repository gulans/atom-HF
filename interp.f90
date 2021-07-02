subroutine interp(m,x,y,xp,yp)
 implicit none
 integer, intent(in) :: m
 real(8), intent(in) :: x(m),y(m),xp
 real(8), intent(out):: yp

 integer :: i, j
 real(8) :: sk,sauc
 yp=0d0
 do i=1,m
   sk=1.0d0
   sauc=1.0d0
   do j=1,m
     if (i.ne.j) then
       sk=sk*(xp-x(j))
       sauc=sauc*(x(i)-x(j))

     endif
 enddo
 yp=yp+sk*y(i)/sauc
 enddo


end subroutine
