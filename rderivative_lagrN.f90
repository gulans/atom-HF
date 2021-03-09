subroutine rderivative_lagrN(Ngrid,r,tools,tools_info,fin,rez)
implicit none
integer, intent(in) :: Ngrid,tools_info(3)
real(8), intent(in) :: r(Ngrid),fin(Ngrid),tools(Ngrid,tools_info(1)) 
real(8), intent(out) :: rez(Ngrid)
integer :: ord,ir,i,Npoints
real(8) :: y(tools_info(2)+1)

ord=tools_info(2)
rez=0d0*r

Npoints=ord+1
do i=1, Npoints
  y(i)=fin(i)
enddo

do ir=1,int(Npoints/2)
   do i=1, Npoints
    rez(ir)=rez(ir)+y(i)*tools(ir,i)
   enddo
 enddo

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




