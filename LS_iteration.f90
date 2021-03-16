subroutine LS_iteration(Ngrid, r,tools,tools_info,version, Z,l, shell_l, nmax,&
                shell0, Nshell, vxc, vh,vx_psi, psi_in,norm_arr,psi,eig)
        ! Ngrid
        ! r
        ! vfull - potential (v_n+v_h+v_xc)
        ! l - quantum number l
        ! num - Number Of eigvals and eigfun to solve
        ! nummax - size of eigval array
        ! eigval (OUT) - erigval array
        ! eigfun (OUT) - eigfun array


!--input and output variables--
implicit none
integer, intent(in) :: tools_info(3),version
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: Ngrid, nmax, shell0, Nshell
integer, intent(in) :: l, shell_l(Nshell) 
real(8), intent(in) :: r(Ngrid),Z,vxc(Ngrid),vh(Ngrid)
real(8), intent(in) :: norm_arr(Nshell)
real(8), intent(in) :: psi_in(Ngrid,Nshell),vx_psi(Ngrid,Nshell)

real(8), intent(inout) :: eig(Nshell)
real(8), intent(out) :: psi(Ngrid,Nshell)
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer :: inn,inp,ish,i,j 
real(8) :: vx_chi(Ngrid,Nshell)
real(8) :: S(nmax,nmax),H(nmax,nmax),Snn,Hnn !Overlap and Hamiltonian matrix
real(8) :: Sevec(nmax,nmax),Seval(nmax),s12(nmax,nmax),W(nmax,nmax),test(nmax,nmax)
real(8) :: Winv(nmax,nmax),x(nmax,nmax),xp(nmax,nmax),Hevec(nmax,nmax),Heval(nmax),Hp(nmax,nmax)
real(8) :: lambda(nmax),temp1(nmax,nmax),temp2(nmax,nmax)
real(8) :: f1(Ngrid),f2(Ngrid),f3(Ngrid),f4(Ngrid),lambda_test(nmax,nmax)
real(8) :: phi(Ngrid,nmax),norm,eigp(Nshell)
integer :: iscl,maxscl,ir
maxscl=20

do iscl=1,maxscl

write(*,*)iscl,". eig-eigp: ",eig-eigp
if((maxval(abs((eig-eigp)/(eig+1d0)))).lt.1d-12)then !1d-13 bija
        write(*,*)"konverģence absolūtā: ",maxval(abs(eig-eigp))," relatīvā: ",maxval(abs((eig-eigp)/(eig+1d0)))

        exit
endif
do inn=1,nmax
  ish=shell0+inn
  call scrPoisson(Ngrid, r,tools,tools_info,version, Z, l, vh, vxc,vx_psi(:,ish),psi_in(:,ish), eig(ish), psi(:,ish))
   call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,ish)**2,norm)
   psi(:,ish)=psi(:,ish)/dsqrt(norm)

enddo

!Overlap and Hamiltonian matrix
if (version.eq.4) then
   do inn=1,nmax
   ish=inn+shell0
   call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,psi(:,ish),psi_in,vx_chi(:,ish))
   enddo
endif

do inn=1,nmax
  do inp=1,nmax
    call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn+shell0)*psi(:,inp+shell0),Snn)


    S(inn,inp)=Snn
    call rderivative_lagrN(Ngrid,r,tools,tools_info,r*psi(:,inp+shell0),f4)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,f4,f1)

    f1=f1/r
    call integ_BodesN_value(Ngrid,r,tools,tools_info,psi(:,inn+shell0)*f3*r**2,Hnn)

if (version.eq.5) then    
!LDA
    f2=(0.5d0*dble(l)*dble(l+1)/r**2-Z/r+vh+vxc)*psi(:,inp+shell0)
else if (version.eq.4) then
!Fock exchange
    f2=(0.5d0*dble(l)*dble(l+1)/r**2-Z/r+vh)*psi(:,inp+shell0)+vx_chi(:,inp+shell0)
endif
!Hydrogen
!    f2=(0.5d0*dble(l)*dble(l+1)/r**2-1d0/r)*psi(:,inp+shell0)

    f3=-0.5d0*f1+f2
     call integ_BodesN_value(Ngrid,r,tools,tools_info,psi(:,inn+shell0)*f3*r**2,Hnn)
    H(inn,inp)=Hnn
  enddo
enddo
!write(*,*)"S:"
do inn=1,nmax
!  write(*,*)S(inn,:)
enddo
!write(*,*)"H:"
do inn=1,nmax
!  write(*,*)H(inn,:)
enddo

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
! print *,"Ja W matrica ir S^(-1/2), tad WSW=I ?:"

!!!!!!!!! Make a substitution x=Wx' un H'=WHW
call inver(W,Winv,nmax)
 xp=matmul(Winv,x)
 Hp=matmul(matmul(W,H),W)
 Hevec=Hp
!!!!!!!! SOLVE H'x'=x'
 call diasym_small(Hevec,lambda,nmax)
! print *,"lambda:"
! print *,lambda(:)
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
  phi(:,inn)=0d0*r
  do inp=1,nmax
    phi(:,inn)=phi(:,inn)+x(inp,inn)*psi(:,inp+shell0)
  enddo 
enddo
!!STORE WF and eigenvalues

  eigp=eig

do inn=1,nmax
  psi(:,inn+shell0)=phi(:,inn)
  eig(inn+shell0)=lambda(inn)
 
  call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn+shell0)**2,norm)

  write(*,*)"l=",l," eig(",inn,")=",eig(inn+shell0),"norm=",norm," eigp(",inn,")=",eigp(inn+shell0)
  enddo

!  open(11,file='psi_wf.dat',status='replace')
!  write(11,*)"r psi1 psi2 psi3"
!   do ir = 1,Ngrid
!      if ((Z.lt.2.5).and.(Z.gt.0.5)) then
!      write(11,*)r(ir),psi(ir,1)
!      elseif ((Z.lt.4.5).and.(Z.gt.3.5)) then
!      write(11,*)r(ir),psi(ir,1),psi(ir,2)
!      elseif ((Z.lt.10.5).and.(Z.gt.9.5)) then
!      write(11,*)r(ir),psi(ir,1),psi(ir,2),psi(ir,3)  
!      endif
!   end do
!   close(11)




enddo
end subroutine

subroutine scrPoisson(Ngrid, r,tools,tools_info,version, Z, l, vh, vxc,vx_phi,phi, e, psi)

integer, intent(in) :: tools_info(3),version
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: Ngrid
integer, intent(in) :: l
real(8), intent(in) :: r(Ngrid),Z,vxc(Ngrid),vh(Ngrid),phi(Ngrid),vx_phi(Ngrid)
real(8), intent(inout) :: e

real(8), intent(out) :: psi(Ngrid)
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8) :: f(Ngrid),lam
real(8) :: besrezi(0:50), besrezk(0:50)
integer :: ri,i
real(8) :: f11(Ngrid),f12(Ngrid),f21(Ngrid),f22(Ngrid),int1(Ngrid),int2(Ngrid)
!result=f11*(integral_(0->r)f12)+f21(integral_(r->Ngrid)f22)
real(8) :: besi,besk


!write(*,*)"e=",e
if (e.gt.0) then
        write(*,*)"scrPoisson Error: positive eigenvalue!"
        e=-1d-3
endif

lam=dsqrt(-2d0*e)

if (version.eq.5)then
!LDA
f=-2d0*(-Z/r+vh+vxc)*phi
else if (version.eq.4)then
!Fock
f=-2d0*(-Z/r+vh)*phi-2d0*vx_phi
endif
!!! Ūdeņradis 1s
!lam=1d0
!f=2d0/r*(2d0*dexp(-1d0*r))
!!!!!!!!!!!!!!
!!! Ūdeņradis 2p
!lam=0.5d0 !E=-0.125d0
!f=1d0/(dsqrt(6d0))*(dexp(-0.5d0*r))
!l=1
!!!!!!!!!!!!!!
!!! Ūdeņradis 3d
!lam=0.333333333333333d0 !E=-0.055555555555556d0
!f=(2d0/81d0)*(dsqrt(2d0/15d0))*(dexp(-0.5d0*r))*r
!l=2
!!!!!!!!!!!!!!

f11=0d0*r
f12=0d0*r
f21=0d0*r
f22=0d0*r
!write(*,*)"argument       besi         besk"

do ri=1, Ngrid
if((lam*r(ri)).gt.100d0) then
exit
endif
call msbesseli(l,lam*r(ri), besrezi)
call msbesselk(l,lam*r(ri), besrezk)
besi=besrezi(l)
besk=besrezk(l)
!call besi_simp(l,lam*r(ri),besi)
!call besk_simp(l,lam*r(ri),besk)

f11(ri)=lam*besk
f12(ri)=besi*f(ri)
f21(ri)=lam*besi
f22(ri)=besk*f(ri)
!write(*,*)lam*r(ri),besi, besk
enddo

!call integ_s38_fun(Ngrid,r,f12*r**2,1,int1)
!call integ_s38_fun(Ngrid,r,f22*r**2,-1,int2)
!call integ_Bodes6_fun(Ngrid,r,f12*r**2,1,int1)
!call integ_Bodes6_fun(Ngrid,r,f22*r**2,-1,int2)

call integ_BodesN_fun(Ngrid,r,tools,tools_info,1,f12*r**2,int1)
call integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,f22*r**2,int2)
psi=f11*int1+f21*int2


!write(*,*)"r*lam","f11","int1","f21","int2","psi"

!if (e.lt.-300d0) then
!if (l.eq.2) then
!write(*,*)"arguments    ","besk*lam    ","integr1    ","besi*lam    ","integr2    ","psi"
do ri=1, Ngrid
!write(*,*)lam*r(ri),f11(ri),int1(ri),f21(ri),int2(ri),psi(ri)
enddo
!endif

end subroutine


subroutine laplacian(Ngrid,r,f1,upp) 
  IMPLICIT NONE
  integer, intent(in)    :: Ngrid 
  real(8), intent(in)    :: r(Ngrid),f1(Ngrid)
  real(8), intent(out)   :: upp(Ngrid)


  real(8)                :: u(Ngrid),up(Ngrid),up1,up2 
  integer :: ir
  u=f1*r
  up2=0d0
  do ir=1, Ngrid-1
    up1=up2
    up2=(u(ir+1)-u(ir))/(r(ir+1)-r(ir))
    up(ir)=(up1+up2)/2d0
  enddo
  up(Ngrid)=up2/2d0 

  up2=0d0
  do ir=1, Ngrid-1
    up1=up2
    up2=(up(ir+1)-up(ir))/(r(ir+1)-r(ir))
    upp(ir)=(up1+up2)/2d0
  enddo
  upp(Ngrid)=up2/2d0

  upp=upp/r
  return
end subroutine

subroutine diasym_small(a,eig,n)
 implicit none

 integer n,l,inf
 real*8  a(n,n),eig(n),work(n*(3+n/2))

 l=n*(3+n/2)
 call dsyev('V','U',n,a,n,eig,work,l,inf)

end subroutine

subroutine inver(A,Ainv,si)
    implicit none
    integer            :: si
    real(8)            :: A(si,si)
    real(8)            :: Ainv(si,si)
    real(8)            :: work(si*si)            ! work array for LAPACK
    integer         :: n,info,ipiv(si)     ! pivot indices

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = si
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) then
      print *, 'Matrix is numerically singular!'
      stop
    endif
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call DGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) then
      print *,  'Matrix inversion failed!'
      stop
    endif

end subroutine
        

!subroutine sph_i(l,lmax,x,rez)
!    implicit none
!    integer            :: l,lmax
!    real(8)            :: x,rez, A1(0:250),A2(0:250)
!    call SPHI(l,x,lmax,A1,A2)
!    rez=A1(l)
!end subroutine
!
!subroutine sph_k(l,lmax,x,rez)
!    implicit none
!    integer            :: l,lmax
!    real(8)            :: x,rez, A1(0:250),A2(0:250)
!    call SPHK(l,x,lmax,A1,A2)
!    rez=A1(l)
!end subroutine

