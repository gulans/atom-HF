module modinteg
implicit none 

real(8),allocatable :: icf(:,:)
real(8),allocatable :: dcf(:,:)
real(8),allocatable :: icv(:)
integer,allocatable :: istart(:)
integer,allocatable :: dstart(:)
logical :: integrate_0_r1
integer :: d_order
integer :: i_order

contains

subroutine gen_icoef(Ngrid,r)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid) 


integer :: ir,i,j
integer :: iNpoints
integer :: dNpoints
real(8) :: hh(Ngrid)

real(8), allocatable :: w(:)
real(8), allocatable :: lc(:,:,:)
real(8), allocatable :: xx1(:)
real(8), allocatable :: xx2(:)

iNpoints=i_order+1
dNpoints=d_order+1
!!!!!module variables!!!!!!
allocate(icf(Ngrid,iNpoints))
allocate(dcf(Ngrid,dNpoints))
allocate(icv(Ngrid))
allocate(istart(Ngrid))
allocate(dstart(Ngrid))
!!!!!subroutine variables!!!!!!!
allocate(w(iNpoints))
allocate(lc(Ngrid,iNpoints,iNpoints-2))
allocate(xx1(iNpoints)) !Grid section
allocate(xx2(dNpoints)) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (iNpoints.eq.2) then !Trapezoidal rule
        w(1)=1d0
        w(2)=1d0
        w=w/2d0
elseif (iNpoints.eq.3) then !Simpson's rule
        w(1)=1d0
        w(2)=4d0
        w(3)=1d0
        w=1d0*w/3d0
elseif (iNpoints.eq.4) then !Simpson's 3/8 rule
        w(1)=1d0
        w(2)=3d0
        w(3)=3d0
        w(4)=1d0
        w=3d0*w/8d0
elseif (iNpoints.eq.5) then ! Boole's rule
        w(1)=7d0
        w(2)=32d0
        w(3)=12d0
        w(4)=32d0
        w(5)=7d0
        w=2d0*w/45d0
elseif (iNpoints.eq.6) then
        w(1)=19d0
        w(2)=75d0
        w(3)=50d0
        w(4)=50d0
        w(5)=75d0
        w(6)=19d0
        w=5d0*w/288d0
elseif (iNpoints.eq.7) then
        w(1)=41d0
        w(2)=216d0
        w(3)=27d0
        w(4)=272d0
        w(5)=27d0
        w(6)=216d0
        w(7)=41d0
        w=w/140d0
elseif (iNpoints.eq.8) then
        w(1)=751d0
        w(2)=3577d0
        w(3)=1323d0
        w(4)=2989d0
        w(5)=w(4)
        w(6)=w(3)
        w(7)=w(2)
        w(8)=w(1)
        w=7d0*w/17280d0
elseif (iNpoints.eq.9) then
        w(1)=989d0
        w(2)=5888d0
        w(3)=-928d0
        w(4)=10496d0
        w(5)=-4540d0
        w(6)=w(4)
        w(7)=w(3)
        w(8)=w(2)
        w(9)=w(1)
        w=4d0*w/14175d0
elseif (iNpoints.eq.10) then
        w(1)=2857d0
        w(2)=15741d0
        w(3)=1080d0
        w(4)=19344d0
        w(5)=5778d0
        w(6)=w(5)
        w(7)=w(4)
        w(8)=w(3)
        w(9)=w(2)
        w(10)=w(1)
        w=9d0*w/89600d0
elseif (iNpoints.eq.11) then
        w(1)=16067d0
        w(2)=106300d0
        w(3)=-48525d0
        w(4)=272400d0
        w(5)=-260550d0
        w(6)=427368d0
        w(7)=w(5)
        w(8)=w(4)
        w(9)=w(3)
        w(10)=w(2)
        w(11)=w(1)
        w=5d0*w/299376d0
endif


do ir=1,Ngrid-1
  hh(ir)=(r(ir+1)-r(ir))/(iNpoints-1)
enddo

!!!!!!!!!!!!!!!Generate istart array !!!!!!!!!!!!!!!

do ir=1,int(iNpoints/2)
   istart(ir)=1
enddo
do ir=int(iNpoints/2)+1, Ngrid-int(iNpoints/2)-1
 if (MOD(iNpoints,2) .eq. 0) then
   istart(ir)=ir-int(iNpoints/2)+1
 else
   istart(ir)=ir-int(iNpoints/2)
 endif
enddo
do ir=Ngrid-int(iNpoints/2),Ngrid-1
   istart(ir)=Ngrid-iNpoints+1
enddo


!!!!!!!!!!!!!!!END Generate istart array !!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!Generate dstart array !!!!!!!!!!!!!!!

do ir=1,int(dNpoints/2)
   dstart(ir)=1
enddo
do ir=int(dNpoints/2)+1, Ngrid-int(dNpoints/2)-1
 if (MOD(dNpoints,2) .eq. 0) then
   dstart(ir)=ir-int(dNpoints/2)+1
 else
   dstart(ir)=ir-int(dNpoints/2)
 endif
enddo
do ir=Ngrid-int(dNpoints/2),Ngrid
   dstart(ir)=Ngrid-dNpoints+1
enddo


!!!!!!!!!!!!!!!END Generate dstart array !!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!Interpolation coeficients !!!!!!!!!!!!!!!!!!!!!

xx1(1)=0d0
xx1(2:iNpoints)=r(1:iNpoints-1)
do j=1,iNpoints-2
   call lagr_c(0d0+j*r(1)/(iNpoints-1),iNpoints,xx1,lc(1,:,j))
enddo

do ir=1,Ngrid-1
   do j=1,iNpoints-2  !How many points have to be interpolated in between grid points
     xx1=r(istart(ir):istart(ir)+iNpoints-1)
     call lagr_c(r(ir)+j*hh(ir),iNpoints,xx1,lc(ir+1,:,j))
   enddo
enddo

!!!!!!!!!!!!!!!!!!!! END Interpolation coeficients!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!! Build array of coficients to get function

do i=1, iNpoints
  icf(:,i)=0d0*r
enddo

if (integrate_0_r1) then
  icf(1,1)=icf(1,1)+w(1)
  icf(1,2)=icf(1,2)+w(iNpoints)
  do i=1, iNpoints
    do j=1, iNpoints-2
      icf(1,i)=icf(1,i)+w(j+1)*lc(1,i,j)
    enddo
  enddo
endif

do ir=1, Ngrid-1
  do i=1, iNpoints
    if (istart(ir)+i-1.eq.ir) then
      icf(ir+1,i)=icf(ir+1,i)+w(1)
    elseif (istart(ir)+i-1.eq.ir+1) then
      icf(ir+1,i)=icf(ir+1,i)+w(iNpoints)
    endif
    do j=1, iNpoints-2
      icf(ir+1,i)=icf(ir+1,i)+w(j+1)*lc(ir+1,i,j)
    enddo
  enddo
enddo


do i=1, iNpoints
  icf(1,i)=icf(1,i)*r(1)/(iNpoints-1)
  do ir=1,Ngrid-1
     icf(ir+1,i)=icf(ir+1,i)*hh(ir)
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!END Build array of coficients to get function

!!!!!!!!!!!!!!!!!!!! Build array of coficients to get integral value

icv=r*0d0

do i=2, iNpoints
    icv(i-1)=icv(i-1)+icf(1,i)
enddo

do ir=1, Ngrid-1
  do i=1, iNpoints
     icv(istart(ir)+i-1)=icv(istart(ir)+i-1)+icf(ir+1,i)
  enddo

enddo


!!!!!!!!!!!!!!!!!!!!END Build array of coficients to get integral value




!!!!!!!!!!!!!!!!!!!! Build array of coficients to get derivative function

do ir=1,Ngrid
  xx2=r(dstart(ir):dstart(ir)+dNpoints-1)
  call lagr_c_d(r(ir),dNpoints,xx2,dcf(ir,:))
enddo


!!!!!!!!!!!!!!!!!!!! Build array of coficients to get derivative function





deallocate(w)
deallocate(lc)
deallocate(xx1)
deallocate(xx2)

end subroutine

subroutine deriv_f(Ngrid,r,fin,fout)

integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: fout(Ngrid)

real(8) :: r1
integer :: i,ir
integer :: Npoints

Npoints=d_order+1


do ir=1, Ngrid
  r1=0d0
  do i=1, Npoints
    r1=r1+dcf(ir,i)*fin(dstart(ir)+i-1)
  enddo
  fout(ir)=r1
enddo

end subroutine

subroutine integ_f(Ngrid,r,fin,fout)

integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: fout(Ngrid)

real(8) :: r1
integer :: ir,i,Npoints
Npoints=i_order+1

if (integrate_0_r1) then
  r1=0d0
  r1=r1+icf(1,1)*0d0
    do i=2, Npoints
      r1=r1+icf(1,i)*fin(i-1)
    enddo
  fout(1)=r1
else
  fout(1)=0d0
endif

do ir=1, Ngrid-1
  r1=0d0
  do i=1, Npoints
    r1=r1+icf(ir+1,i)*fin(istart(ir)+i-1)
  enddo
  fout(ir+1)=fout(ir)+r1
enddo

end subroutine

subroutine integ_f_rev(Ngrid,r,fin,fout)

integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: fout(Ngrid)

real(8) :: r1
integer :: ir,i,Npoints
Npoints=i_order+1

fout(Ngrid)=0d0

do ir=Ngrid-1, Ngrid-int(Npoints/2),-1
  !interpolation gives incorect results if the fin functioon is small in the current point, but is a few orders higher around the
  !point ,trapezoidal rule is used insead:
  fout(ir)=fout(ir+1)+(r(ir+1)-r(ir))*(fin(ir)+fin(ir+1))/2d0
enddo


do ir=Ngrid-int(Npoints/2)-1, 1,-1
  r1=0d0
  do i=1, Npoints
    r1=r1+icf(ir+1,i)*fin(istart(ir)+i-1)
  enddo
  fout(ir)=fout(ir+1)+r1

enddo

end subroutine





subroutine integ_v(Ngrid,r,fin,vout)

integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
real(8), intent(in) :: fin(Ngrid)
real(8), intent(out) :: vout

real(8) :: r1
integer :: ir

r1=0d0
do ir=1, Ngrid
  r1=r1+icv(ir)*fin(ir)
enddo

vout=r1

end subroutine

subroutine integ_cf(Ngrid,r,fin,rez)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
complex(8), intent(in)  :: fin(Ngrid)
complex(8), intent(out) :: rez(Ngrid)
integer :: ir

real(8) :: finRe(Ngrid),finIm(Ngrid),rezRe(Ngrid),rezIm(Ngrid)

do ir=1, Ngrid
  finRe(ir)=realpart(fin(ir))
  finIm(ir)=imagpart(fin(ir))
enddo

call integ_f(Ngrid,r,finRe,rezRe)
call integ_f(Ngrid,r,finIm,rezIm)

do ir=1, Ngrid
  rez(ir)=cmplx(rezRe(ir),rezIm(ir),8)
enddo
end subroutine

subroutine dealoc_icoef()
implicit none

deallocate(icf)
deallocate(icv)
deallocate(istart)
deallocate(dstart)
end subroutine

end module
