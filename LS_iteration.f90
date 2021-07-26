subroutine LS_iteration(Ngrid, r,tools,tools_info,rsfunC,Nrsfun,hybx_w, Z,l,sp,shell_l,shell_occ, nmax,&
                shell0, Nshell,Nspin,lmax, vxc, vh,vx_psi,vx_psi_sr, psi_in,psi,eig,Bess_ik)
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
integer, intent(in) :: Nrsfun,lmax
complex(8), intent(in) :: rsfunC(Nrsfun+1,2)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1)),hybx_w(5,2)
integer, intent(in) :: Ngrid, nmax, shell0, Nshell,Nspin
integer, intent(in) :: l,sp, shell_l(Nshell) 
real(8), intent(in) :: r(Ngrid),Z,vxc(Ngrid,Nspin),vh(Ngrid),shell_occ(Nshell,Nspin)
real(8), intent(in) :: psi_in(Ngrid,Nshell,Nspin),vx_psi(Ngrid,Nshell,Nspin),vx_psi_sr(Ngrid,Nshell,Nspin)
complex(8), intent(in) ::Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
real(8), intent(inout) :: eig(Nshell,Nspin)
real(8), intent(out) :: psi(Ngrid,Nshell,Nspin)
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer :: inn,inp,ish,i,j 
real(8) :: vx_chi(Ngrid,Nshell,Nspin),vx_chi_sr(Ngrid,Nshell,Nspin)
real(8) :: S(nmax,nmax),H(nmax,nmax),Snn,Hnn !Overlap and Hamiltonian matrix
real(8) :: Sevec(nmax,nmax),Seval(nmax),s12(nmax,nmax),W(nmax,nmax),test(nmax,nmax)
real(8) :: Winv(nmax,nmax),x(nmax,nmax),xp(nmax,nmax),Hevec(nmax,nmax),Heval(nmax),Hp(nmax,nmax)
real(8) :: lambda(nmax),temp1(nmax,nmax),temp2(nmax,nmax)
real(8) :: f1(Ngrid),f2(Ngrid),f3(Ngrid),f4(Ngrid),lambda_test(nmax,nmax)
real(8) :: phi(Ngrid,nmax),norm,eigp(Nshell,Nspin)
integer :: iscl,maxscl,ir

maxscl=20
do iscl=1,maxscl

!write(*,*)iscl,". eig-eigp: ",eig-eigp
!convergence check
if((maxval(abs((eig-eigp)/(eig-1d0)))).lt.1d-13)then 
       !write(*,*)"Convergence of internal cycle reached, iteration ",iscl
       !write(*,*)"max(eig-eigp) absolute : ",maxval(abs(eig-eigp))," relative: ",maxval(abs(eig-eigp)/(abs(eig-1d0)))
        exit
endif

do inn=1,nmax
  ish=shell0+inn
  call scrPoisson(Ngrid, r,tools,tools_info,hybx_w, Z, l, vh, vxc(:,sp),&
          vx_psi(:,ish,sp),vx_psi_sr(:,ish,sp),psi_in(:,ish,sp), eig(ish,sp), psi(:,ish,sp))
   call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,ish,sp)**2,norm)
   psi(:,ish,sp)=psi(:,ish,sp)/dsqrt(norm)
enddo

!Overlap and Hamiltonian matrix
if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
   do inn=1,nmax
   ish=inn+shell0
!   call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ,lmax,&
!           psi(:,ish),psi_in,vx_chi(:,ish),vx_chi_sr(:,ish),rsfunC,Nrsfun,hybx_w,Bess_ik)
   enddo
endif

do inn=1,nmax
  do inp=1,nmax
    call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn+shell0,sp)*psi(:,inp+shell0,sp),Snn)
    S(inn,inp)=Snn
!!!!!
!    call rderivative_lagrN(Ngrid,r,tools,tools_info,r*psi(:,inp+shell0),f4)
!    call rderivative_lagrN(Ngrid,r,tools,tools_info,f4,f1)
!    f1=f1/r
!    call integ_BodesN_value(Ngrid,r,tools,tools_info,-0.5d0*psi(:,inn+shell0)*f1*r**2,Hnn)
!!!!! 
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inp+shell0,sp),f4)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,inn+shell0,sp),f1)
    call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*f4*f1*r**2,Hnn)
!!!!!


!Hybrid exchange
    f2=(0.5d0*dble(l)*dble(l+1)/r**2-Z/r+vh+vxc(:,sp))*psi(:,inp+shell0,sp)+&
            hybx_w(4,1)*vx_chi(:,inp+shell0,sp)+hybx_w(5,1)*vx_chi_sr(:,inp+shell0,sp)
    H(inn,inp)=Hnn
    call integ_BodesN_value(Ngrid,r,tools,tools_info,psi(:,inn+shell0,sp)*f2*r**2,Hnn)
    H(inn,inp)=H(inn,inp)+Hnn
  enddo
enddo
!write(*,*)"S:"
!do inn=1,nmax
!  write(*,*)S(inn,:)
!enddo
!write(*,*)"H:"
!do inn=1,nmax
!  write(*,*)H(inn,:)
!enddo

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
    phi(:,inn)=phi(:,inn)+x(inp,inn)*psi(:,inp+shell0,sp)
  enddo 
enddo
!!STORE WF and eigenvalues

  eigp=eig

do inn=1,nmax
  psi(:,inn+shell0,sp)=phi(:,inn)
  eig(inn+shell0,sp)=lambda(inn)
 
!  call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn+shell0)**2,norm)
!  write(*,*)"l=",l," eig(",inn,")=",eig(inn+shell0),"norm=",norm," eigp(",inn,")=",eigp(inn+shell0)
  enddo

enddo
end subroutine

subroutine scrPoisson(Ngrid, r,tools,tools_info,hybx_w,&
                Z, l, vh, vxc,vx_phi,vx_phi_sr,phi, e, psi)

integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1)),hybx_w(5,2)
integer, intent(in) :: Ngrid
integer, intent(in) :: l
real(8), intent(in) :: r(Ngrid),Z,vxc(Ngrid),vh(Ngrid),phi(Ngrid),vx_phi(Ngrid),&
        vx_phi_sr(Ngrid)
real(8), intent(inout) :: e
real(8) :: vx(Ngrid),vc(Ngrid)
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
        !write(*,*)"scrPoisson Error: positive eigenvalue!"
        e=-1d-3
endif

lam=dsqrt(-2d0*e)


f=-2d0*( (-Z/r+vh+vxc)*phi + hybx_w(4,1)*vx_phi + hybx_w(5,1)*vx_phi_sr )


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

f11(ri)=lam*besk
f12(ri)=besi*f(ri)
f21(ri)=lam*besi
f22(ri)=besk*f(ri)
enddo

call integ_BodesN_fun(Ngrid,r,tools,tools_info,1,f12*r**2,int1)
call integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,f22*r**2,int2)
psi=f11*int1+f21*int2

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
