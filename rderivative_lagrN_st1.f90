subroutine rderivative_lagrN_st1(Ngrid,r,tools,tools_info,fin,rez)
implicit none
integer, intent(in) :: Ngrid,tools_info(3)
real(8), intent(in) :: r(Ngrid),fin(Ngrid),tools(Ngrid,tools_info(1)) 
real(8), intent(out) :: rez(Ngrid)
integer :: ord,ir,i,Npoints
real(8) :: y(tools_info(2)+1)

integer, parameter :: Nstart=2
real(8) :: xstart(Nstart),ystart(Nstart)
real(8) :: lcstart(Nstart)
ord=tools_info(2)
rez=0d0*r
Npoints=ord+1








do ir=1,int(Npoints/2)
   do i=1, Nstart
!     write(*,*)"ir=",ir, "input=",ir+i-2
     xstart(i)=r(ir+i-1)
     ystart(i)=fin(ir+i-1)
   enddo
   call lagr_c_d(r(ir),Nstart,xstart,lcstart)
   do i=1, Nstart
     rez(ir)=rez(ir)+ystart(i)*lcstart(i)
   enddo


enddo






if (.false.)then
do i=1, Npoints
  y(i)=fin(i)
enddo

do ir=1,int(Npoints/2)
   do i=1, Npoints
    rez(ir)=rez(ir)+y(i)*tools(ir,i)
   enddo
 enddo
endif !!end comment

do ir=int(Npoints/2)+1, Ngrid-int(Npoints/2)-1
  do i=1, Npoints
    if (MOD(Npoints,2) .eq. 0) then
     !Npoints even
     y(i)=fin(ir-int(Npoints/2)+i)
    else
     !odd
      y(i)=fin(ir-int(Npoints/2)-1+i)
    endif
    enddo
   do i=1, Npoints
    rez(ir)=rez(ir)+y(i)*tools(ir,i)

   enddo
enddo

do i=1, Npoints
  y(i)=fin(Ngrid-Npoints+i)
enddo

do ir=Ngrid-int(Npoints/2),Ngrid
   do i=1, Npoints
    rez(ir)=rez(ir)+y(i)*tools(ir,i)
   enddo
enddo



end subroutine




