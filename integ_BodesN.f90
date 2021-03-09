subroutine integ_BodesN_fun(Ngrid,r,tools,tools_info,dir,fin,rez)
implicit none
integer, intent(in) :: Ngrid,tools_info(3),dir
real(8), intent(in) :: r(Ngrid),fin(Ngrid),tools(Ngrid,tools_info(1)) 
real(8), intent(out) :: rez(Ngrid)
integer :: ord,ir,i,Npoints
real(8) :: y(tools_info(3)+1),interp(3),h
ord=tools_info(3)
rez=0d0*r
Npoints=ord+1

!!!!!!!!!!!!!!!!!!!!Interpolation coeficients for integral calculation!!!!!!!!!!!!!!!!!!!!!
if (dir.eq.1) then

rez(1)=0d0

Npoints=tools_info(3)+1
do i=1, Npoints
  y(i)=fin(i)
enddo

do ir=1,int(Npoints/2)
   h=(r(ir+1)-r(ir))/4d0
   interp(1)=0d0
   interp(2)=0d0
   interp(3)=0d0
   do i=1, Npoints
     interp(1)=interp(1)+y(i)*tools(ir,10+i)
     interp(2)=interp(2)+y(i)*tools(ir,20+i)
     interp(3)=interp(3)+y(i)*tools(ir,30+i)
     enddo
  rez(ir+1)=rez(ir)+h*(14d0*fin(ir)+64d0*interp(1)+24d0*interp(2)+64d0*interp(3)+14d0*fin(ir+1))/45d0
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

   h=(r(ir+1)-r(ir))/4d0
   interp(1)=0d0
   interp(2)=0d0
   interp(3)=0d0
   do i=1, Npoints
     interp(1)=interp(1)+y(i)*tools(ir,10+i)
     interp(2)=interp(2)+y(i)*tools(ir,20+i)
     interp(3)=interp(3)+y(i)*tools(ir,30+i)
   enddo
   rez(ir+1)=rez(ir)+h*(14d0*fin(ir)+64d0*interp(1)+24d0*interp(2)+64d0*interp(3)+14d0*fin(ir+1))/45d0
    
    
enddo

do i=1, Npoints
  y(i)=fin(Ngrid-Npoints+i)
enddo

do ir=Ngrid-int(Npoints/2),Ngrid-1

   h=(r(ir+1)-r(ir))/4d0
   interp(1)=0d0
   interp(2)=0d0
   interp(3)=0d0
   do i=1, Npoints
     interp(1)=interp(1)+y(i)*tools(ir,10+i)
     interp(2)=interp(2)+y(i)*tools(ir,20+i)
     interp(3)=interp(3)+y(i)*tools(ir,30+i)
   enddo
   rez(ir+1)=rez(ir)+h*(14d0*fin(ir)+64d0*interp(1)+24d0*interp(2)+64d0*interp(3)+14d0*fin(ir+1))/45d0

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

elseif (dir.eq.-1)then

!       do ir=1,Ngrid
!          rez(ir)=rez(Ngrid)-rez(ir)
!       enddo

!!BEIGAS
rez(Ngrid)=0d0

do i=1, Npoints
  y(i)=fin(Ngrid-Npoints+i)
enddo

do ir=Ngrid-1,Ngrid-int(Npoints/2),-1

   h=(r(ir+1)-r(ir))/4d0
   interp(1)=0d0
   interp(2)=0d0
   interp(3)=0d0
   do i=1, Npoints
     interp(1)=interp(1)+y(i)*tools(ir,10+i)
     interp(2)=interp(2)+y(i)*tools(ir,20+i)
     interp(3)=interp(3)+y(i)*tools(ir,30+i)
   enddo
   rez(ir)=rez(ir+1)+h*(14d0*fin(ir)+64d0*interp(1)+24d0*interp(2)+64d0*interp(3)+14d0*fin(ir+1))/45d0

enddo

!!!VIDUS
do ir=Ngrid-int(Npoints/2)-1,int(Npoints/2)+1,-1 
  do i=1, Npoints
    if (MOD(Npoints,2) .eq. 0) then
     !Npoints even
     y(i)=fin(ir-int(Npoints/2)+i)
    else
     !odd
     y(i)=fin(ir-int(Npoints/2)-1+i)
    endif
  enddo

   h=(r(ir+1)-r(ir))/4d0
   interp(1)=0d0
   interp(2)=0d0
   interp(3)=0d0
   do i=1, Npoints
     interp(1)=interp(1)+y(i)*tools(ir,10+i)
     interp(2)=interp(2)+y(i)*tools(ir,20+i)
     interp(3)=interp(3)+y(i)*tools(ir,30+i)
   enddo
   rez(ir)=rez(ir+1)+h*(14d0*fin(ir)+64d0*interp(1)+24d0*interp(2)+64d0*interp(3)+14d0*fin(ir+1))/45d0
    
    
enddo

!!!Sākums


Npoints=tools_info(3)+1
do i=1, Npoints
  y(i)=fin(i)
enddo

do ir=int(Npoints/2),1,-1
   h=(r(ir+1)-r(ir))/4d0
   interp(1)=0d0
   interp(2)=0d0
   interp(3)=0d0
   do i=1, Npoints
     interp(1)=interp(1)+y(i)*tools(ir,10+i)
     interp(2)=interp(2)+y(i)*tools(ir,20+i)
     interp(3)=interp(3)+y(i)*tools(ir,30+i)
     enddo
  rez(ir)=rez(ir+1)+h*(14d0*fin(ir)+64d0*interp(1)+24d0*interp(2)+64d0*interp(3)+14d0*fin(ir+1))/45d0
 enddo


endif

end subroutine


subroutine integ_BodesN_value(Ngrid,r,tools,tools_info,fin,rezi)

integer, intent(in) :: Ngrid,tools_info(3)
real(8), intent(in) :: r(Ngrid),fin(Ngrid),tools(Ngrid,tools_info(1))
real(8), intent(out) :: rezi

real(8) :: rez(Ngrid)

call integ_BodesN_fun(Ngrid,r,tools,tools_info,1,fin,rez)

rezi=rez(Ngrid)

end subroutine
