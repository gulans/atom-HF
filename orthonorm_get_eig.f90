subroutine orthonorm_get_eig(Ngrid,r,tools,tools_info,Z,l,nmax,relativity,v_rel,hybx_w,&
        vxc,vh,vx_psi,vx_psi_sr,&
        psi,eig, vx_psi_out, vx_psi_sr_out)

use xc_f03_lib_m
implicit none
integer(8), intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
real(8), intent(in) :: Z
integer, intent(in) :: l
integer, intent(in) :: nmax
logical, intent(in) :: relativity
real(8), intent(in) :: v_rel(Ngrid)
real(8), intent(in) :: hybx_w(5,2)
real(8), intent(in) :: vxc(Ngrid)
real(8), intent(in) :: vh(Ngrid)
real(8), intent(in) :: vx_psi(Ngrid,nmax)
real(8), intent(in) :: vx_psi_sr(Ngrid,nmax)
real(8), intent(inout) :: psi(Ngrid,nmax) 
real(8), intent(out) :: eig(nmax)
real(8), intent(out) :: vx_psi_out(Ngrid,nmax)
real(8), intent(out) :: vx_psi_sr_out(Ngrid,nmax)



real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)

integer :: inn,inp,i,j
real(8) :: Snn,Hnn
real(8) :: S(nmax,nmax),H(nmax,nmax)
real(8) :: Hevec(nmax,nmax),Hp(nmax,nmax), Sevec(nmax,nmax),Seval(nmax), s12(nmax,nmax)
real(8) :: lambda_test(nmax,nmax),lambda(nmax)
real(8) :: x(nmax,nmax),xp(nmax,nmax),W(nmax,nmax),Winv(nmax,nmax)
real(8) :: psi_new(Ngrid,nmax)
real(8) :: f(Ngrid),f1(Ngrid),f2(Ngrid),f3(Ngrid),f4(Ngrid),f5(Ngrid)


vx_psi_out=0d0*vx_psi
vx_psi_sr_out=0d0*vx_psi

do inn=1,nmax !matrix is symetric !lower triangular matrix 
do inp=inn,nmax
    call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn)*psi(:,inp),Snn)
    S(inn,inp)=Snn
    S(inp,inn)=Snn
if(.not.relativity)then
!    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inp),f2)
!    call rderivative_lagrN(Ngrid,r,tools,tools_info,f2*r**2,f3)
!    call integ_BodesN_value(Ngrid,r,tools,tools_info,-0.5d0*psi(:,inn)*f3, H(inn,inp))
    
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inp),f4)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inn),f1)
    call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*f4*f1*r**2,H(inn,inp))
!write(*,*) "H1", H(inn,inp)

    f=(0.5d0*dble(l)*dble(l+1)/r**2-Z/r+vh+vxc)*psi(:,inp)+&
            hybx_w(4,1)*vx_psi(:,inp)+hybx_w(5,1)*vx_psi_sr(:,inp)
    call integ_BodesN_value(Ngrid,r,tools,tools_info,psi(:,inn)*f*r**2,Hnn)

    H(inn,inp)=H(inn,inp)+Hnn
    H(inp,inn)=H(inn,inp)
else
    f1=1d0/(1d0-v_rel*alpha2)
 
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inp),f2)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inn),f3)
    call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*f2*f1*f3*r**2,H(inn,inp))



    f=(0.5d0*f1*dble(l*(l+1))/r**2-Z/r+vh+vxc)*psi(:,inp)
    call integ_BodesN_value(Ngrid,r,tools,tools_info,psi(:,inn)*f*r**2,Hnn)
    H(inn,inp)=H(inn,inp)+Hnn
    H(inp,inn)=H(inn,inp)
endif
    enddo
enddo

if (.false.)then
  write(*,*)"S:"
  do inn=1,nmax
    write(*,*)S(inn,:)
  enddo
  write(*,*)"H:"
  do inn=1,nmax
    write(*,*)H(inn,:)
  enddo
endif!print and stop

!!!!!!!!Caulculate W=S^-0.5 matrix
 Sevec=s
 call diasym_small(Sevec,Seval,nmax)
 
!!!!!!!! construct s12 matrix, (eigenvalue diagonal-matrix with elements in power -0.5)
 do i=1, nmax
   do j=1, nmax
     if (i==j) Then
       s12(i,j)=Seval(i)**(-0.5d0)
     else
       s12(i,j)=0d0
     endif
   enddo
 enddo
 W=matmul(matmul(Sevec,s12),transpose(Sevec))

!!!!!!!!! Make a substitution x=Wx' un H'=WHW
call inver(W,Winv,nmax)
 xp=matmul(Winv,x)
 Hp=matmul(matmul(W,H),W)
 Hevec=Hp
!!!!!!!! SOLVE H'x'=x'
 call diasym_small(Hevec,lambda,nmax)

! print *,"lambda:"
! print *,lambda(:)
!if(iscl.eq.2)then
!stop
!endif


 do i=1, nmax
   do j=1, nmax
     if (i.eq.j) then
       lambda_test(i,j)=lambda(i)
       else
       lambda_test(i,j)=0d0
     endif
   enddo
 enddo

 x=matmul(W,Hevec)

! print *," Test if the equation is solved.. "
! temp1=matmul(H,x)
!
! temp2=matmul(matmul(S,x),lambda_test)
! write(*,*)"One side:"
! do i=1, nmax
!   print *,temp1(i,:)
! enddo 
!write(*,*)"Other side"
! do i=1, nmax 
!   print *,temp2(i,:)
! enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!CONSTRUCT wave functions
do inn=1,nmax
  psi_new(:,inn)=0d0*r
  do inp=1,nmax
    psi_new(:,inn)=psi_new(:,inn)+x(inp,inn)*psi(:,inp)
  enddo 
enddo

!!CONSTRUCT non-local exhange


if (abs(hybx_w(4,1)).gt.1d-20) then
do inn=1,nmax
  do inp=1,nmax
    vx_psi_out(:,inn)=vx_psi_out(:,inn)+x(inp,inn)*vx_psi(:,inp)
  enddo
enddo
endif

if (abs(hybx_w(5,1)).gt.1d-20) then
do inn=1,nmax
  do inp=1,nmax
    vx_psi_sr_out(:,inn)=vx_psi_sr_out(:,inn)+x(inp,inn)*vx_psi_sr(:,inp)
  enddo
enddo
endif




!!STORE WF and eigenvalues



do inn=1,nmax
  psi(:,inn)=psi_new(:,inn)
  eig(inn)=lambda(inn)
 
!  call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn+shell0)**2,norm)
!  write(*,*)"l=",l," eig(",inn,")=",eig(inn+shell0),"norm=",norm," eigp(",inn,")=",eigp(inn+shell0)
  enddo



end subroutine




