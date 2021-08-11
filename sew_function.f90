subroutine sew_function(Ngrid,r,tools,tools_info,Z,l,r_sew,eig,psi_in,psi_out)

implicit none
integer(8), intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
real(8), intent(in) :: Z
integer, intent(in) :: l
real(8), intent(in) :: r_sew
real(8), intent(in) :: eig
real(8), intent(in) :: psi_in(Ngrid)
real(8), intent(out) :: psi_out(Ngrid)


real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)
real(8), PARAMETER :: lspeed=137.035999084d0
integer, PARAMETER :: Nbase=10
real(8) :: fitc(Nbase)
real(8),allocatable :: fitbase(:,:)
integer :: ir_sew,ir,i
real(8) :: a,b,psi_sew(Ngrid),norm_sew,ftemp1(Ngrid)


do ir=1,Ngrid
if (r(ir).ge.r_sew)then
    ir_sew=ir
    exit
  endif
enddo

allocate(fitbase(ir_sew,Nbase))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.true.)then !analytic result r**a+b*r**(a+1)
  a=(-lspeed**2+sqrt(lspeed**4+lspeed**4*dble(l)+lspeed**4*dble(l)**2-lspeed**2*Z**2))/lspeed**2
  b=(-2d0*lspeed**4 - eig*Z**2 + 2d0*lspeed**2*(-Z**2 + Sqrt(lspeed**4*(1d0 +dble(l)+dble(l)**2)&
        - lspeed**2*Z**2)))/(Z*(lspeed**2 + 2d0*Sqrt(lspeed**4*(1d0 +dble(l)+dble(l)**2) - lspeed**2*Z**2)))

  psi_sew=r**a+b*r**(a+1d0)
  norm_sew=psi_in(ir_sew)/psi_sew(ir_sew)
  psi_sew=psi_sew*norm_sew
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.false.)then

 do ir=1, ir_sew
 do i=0,Nbase-1
    fitbase(ir,i+1)=r(ir)**i
  enddo
  enddo
call lin_fit(ir_sew,Nbase,r(1:ir_sew),psi_in(1:ir_sew),fitbase,fitc)
psi_sew=0d0*r
do ir=1, ir_sew
 do i=0,Nbase-1
    psi_sew(ir)=psi_sew(ir)+fitc(i+1)*fitbase(ir,i+1)
 enddo
 enddo

! do i=1, Nbase
! write(*,*)fitc(i)
! enddo
  norm_sew=psi_in(ir_sew)/psi_sew(ir_sew)
  psi_sew=psi_sew*norm_sew


! open(11,file="poly_fit.dat",status='replace')
! do ir=1, Ngrid
!  write(11,*)r(ir),psi_in(ir),psi_sew(ir)
! enddo
! close(11)
!stop
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
psi_out=psi_in
 do ir=1, ir_sew
   psi_out(ir)=psi_sew(ir)
  enddo

end subroutine




