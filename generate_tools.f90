subroutine generate_tools(Ngrid,r,tools,tools_info)
implicit none
integer, intent(in) :: Ngrid,tools_info(3)
real(8), intent(in) :: r(Ngrid)
real(8), intent(out) :: tools(Ngrid,tools_info(1))

integer :: Npoints,ir,i
real(8) :: x(tools_info(2)+1), lc(tools_info(2)+1),h
real(8) :: xx(tools_info(3)+1), lc1(tools_info(3)+1), lc2(tools_info(3)+1),lc3(tools_info(3)+1)



!!!!!!!!!!!!!!!!!!!!Interpolation coeficients for derivative calculation!!!!!!!!!!!!!!!!!!!!!
Npoints=tools_info(2)+1
do i=1, Npoints
  x(i)=r(i)
enddo

do ir=1,int(Npoints/2)
   call lagr_c_d(r(ir),Npoints,x,lc)
   do i=1, Npoints
     tools(ir,i)=lc(i)
   enddo
!    write(*,*)r(ir),x(:)

 enddo


do ir=int(Npoints/2)+1, Ngrid-int(Npoints/2)-1
  do i=1, Npoints
    if (MOD(Npoints,2) .eq. 0) then
     !Npoints even
     x(i)=r(ir-int(Npoints/2)+i)
    else
     !odd
      x(i)=r(ir-int(Npoints/2)-1+i)
    endif
    enddo
!    write(*,*)r(ir),x(:)
   call lagr_c_d(r(ir),Npoints,x,lc)
   do i=1, Npoints
     tools(ir,i)=lc(i)
   enddo
enddo


do i=1, Npoints
  x(i)=r(Ngrid-Npoints+i)
enddo

do ir=Ngrid-int(Npoints/2),Ngrid
   call lagr_c_d(r(ir),Npoints,x,lc)
   do i=1, Npoints
     tools(ir,i)=lc(i)
   enddo
!    write(*,*)r(ir),x(:)

enddo

!!!!!!!!!!!!!!!!!!!!Interpolation coeficients for integral calculation!!!!!!!!!!!!!!!!!!!!!


Npoints=tools_info(3)+1
do i=1, Npoints
  xx(i)=r(i)
enddo

do ir=1,int(Npoints/2)
   h=(r(ir+1)-r(ir))/4d0
   call lagr_c(r(ir)+ h,Npoints,xx,lc1)
   call lagr_c(r(ir)+ 2d0*h,Npoints,xx,lc2)
   call lagr_c(r(ir)+ 3d0*h,Npoints,xx,lc3)
   do i=1, Npoints
     tools(ir,10+i)=lc1(i)
     tools(ir,20+i)=lc2(i)
     tools(ir,30+i)=lc3(i)
   enddo
!    write(*,*)r(ir)+0.5d0,xx(:)

 enddo
!write(*,*)

do ir=int(Npoints/2)+1, Ngrid-int(Npoints/2)-1
   h=(r(ir+1)-r(ir))/4d0
   do i=1, Npoints
    if (MOD(Npoints,2) .eq. 0) then
     !Npoints even
     xx(i)=r(ir-int(Npoints/2)+i)
    else
     !odd
     xx(i)=r(ir-int(Npoints/2)-1+i)
    endif
    enddo
!   write(*,*)r(ir)+0.5d0,xx(:)
   call lagr_c(r(ir)+ h,Npoints,xx,lc1)
   call lagr_c(r(ir)+ 2d0*h,Npoints,xx,lc2)
   call lagr_c(r(ir)+ 3d0*h,Npoints,xx,lc3)

   do i=1, Npoints
     tools(ir,10+i)=lc1(i)
     tools(ir,20+i)=lc2(i)
     tools(ir,30+i)=lc3(i)
   enddo
enddo


do i=1, Npoints
  xx(i)=r(Ngrid-Npoints+i)
enddo

do ir=Ngrid-int(Npoints/2),Ngrid-1
   h=(r(ir+1)-r(ir))/4d0
   call lagr_c(r(ir)+ h,Npoints,xx,lc1)
   call lagr_c(r(ir)+ 2d0*h,Npoints,xx,lc2)
   call lagr_c(r(ir)+ 3d0*h,Npoints,xx,lc3)

  do i=1, Npoints
     tools(ir,10+i)=lc1(i)
     tools(ir,20+i)=lc2(i)
     tools(ir,30+i)=lc3(i)
   enddo
!    write(*,*)r(ir)+0.5d0,xx(:)

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine

subroutine lagr_c(xx,k,x,l)
implicit none
integer, intent(in) :: k
real(8), intent(in) :: x(k),xx
real(8), intent(out):: l(k)

real(8) :: prod
integer :: i,j,m

!Get coeficients
do j=1,k
 prod=1d0
  do m=1,k
    if (m.ne.j) then
      prod=prod*(xx-x(m))/(x(j)-x(m))
    endif
  enddo
  l(j)=prod
enddo
end subroutine

subroutine lagr_c_d(xx,k,x,l)
implicit none
integer, intent(in) :: k
real(8), intent(in) :: x(k),xx
real(8), intent(out):: l(k)

real(8) :: prod
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

end subroutine


