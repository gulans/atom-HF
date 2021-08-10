subroutine LS_iteration(Ngrid, r,tools,tools_info,rsfunC,Nrsfun,hybx_w, Z,l,sp,shell_l,shell_occ, nmax,&
                shell0, Nshell,Nspin,relativity,lmax,vxc,v_rel, vh,vx_psi,vx_psi_sr, psi_in,psi,eig,Bess_ik,&
                iner_loop,xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func)
        ! Ngrid
        ! r
        ! vfull - potential (v_n+v_h+v_xc)
        ! l - quantum number l
        ! num - Number Of eigvals and eigfun to solve
        ! nummax - size of eigval array
        ! eigval (OUT) - erigval array
        ! eigfun (OUT) - eigfun array
use xc_f03_lib_m


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
logical, intent(in) :: relativity
real(8), intent(in) :: v_rel(Ngrid)
complex(8), intent(in) ::Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
real(8), intent(inout) :: eig(Nshell,Nspin)
integer, intent(out) :: iner_loop
real(8), intent(out) :: psi(Ngrid,Nshell,Nspin)
integer, intent(in) ::xc1_num,xc2_num,xc3_num
TYPE(xc_f03_func_t) :: xc1_func,xc2_func,xc3_func

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)
integer :: inn,inp,ish,i,j 
real(8) :: vx_chi(Ngrid,Nshell,Nspin),vx_chi_sr(Ngrid,Nshell,Nspin)
real(8) :: S(nmax,nmax),H(nmax,nmax),Snn,Hnn !Overlap and Hamiltonian matrix
real(8) :: Sevec(nmax,nmax),Seval(nmax),s12(nmax,nmax),W(nmax,nmax),test(nmax,nmax)
real(8) :: Winv(nmax,nmax),x(nmax,nmax),xp(nmax,nmax),Hevec(nmax,nmax),Heval(nmax),Hp(nmax,nmax)
real(8) :: lambda(nmax),temp1(nmax,nmax),temp2(nmax,nmax)
real(8) :: f1(Ngrid),f2(Ngrid),f3(Ngrid),f4(Ngrid),f5(Ngrid),lambda_test(nmax,nmax)
real(8) :: f6(Ngrid),f7(Ngrid)
real(8) :: phi(Ngrid,nmax),norm,eigp(Nshell,Nspin)
integer :: iscl,maxscl,ir,isp
real(8) :: f(Ngrid)
!real(8) :: vx_psi(Ngrid,Nshell,Nspin),vx_psi_sr(Ngrid,Nshell,Nspin)
!real(8) :: vh1(Ngrid),vxc1(Ngrid,Nspin)
!real(8) :: rho1(Ngrid),exc1(Ngrid,Nspin) !not used
logical :: spin

do ish=shell0+1,shell0+nmax
  psi(:,ish,sp)=psi_in(:,ish,sp)
enddo

!vxc=vxc_in
!vh=vh_in
!vx_psi_sr=vx_psi_sr_in
!vx_psi=vx_psi_in

do isp=1, Nspin
do ish=1, Nshell
vx_chi(:,ish,isp)=0d0*r
vx_chi_sr(:,ish,isp)=0d0*r
enddo
enddo

if (Nspin.eq.2)then
        spin=.true.
else
        spin=.false.
endif


maxscl=20
iner_loop=maxscl+1
do iscl=1,maxscl

!write(*,*)l,eig,iscl!,". eig-eigp: ",eig-eigp
!convergence check
if((maxval(abs((eig-eigp)/(eig-1d0)))).lt.1d-13)then
        iner_loop=iscl-1
       !write(*,*)"Convergence of internal cycle reached, iteration ",iscl
       !write(*,*)"max(eig-eigp) absolute : ",maxval(abs(eig-eigp))," relative: ",maxval(abs(eig-eigp)/(abs(eig-1d0)))
        exit
endif

do inn=1,nmax
  ish=shell0+inn

  if (.not.relativity)then 
  f=-2d0*( (-Z/r+vh+vxc(:,sp))*psi_in(:,ish,sp) + hybx_w(4,1)*vx_psi(:,ish,sp)&
          + hybx_w(5,1)*vx_psi_sr(:,ish,sp) )
  else

 ! call rderivative_lagrN_st3(Ngrid,r,tools,tools_info,psi_in(:,ish,sp),f1)
 ! call rderivative_lagrN_st3(Ngrid,r,tools,tools_info,f1,f2)
 ! f3=1+alpha2*Z/r

 ! f4 =alpha2*r**(-2)*( -Z/f3 -(alpha2*Z**2/r)/f3**2 )*f1
 ! f6 =alpha2*r**(-2)*(Z*r/f3)*f2 

 ! f5 = -Z*r**(-1)*alpha2/(1d0+Z*r**(-1)*alpha2)
 ! f7= 2d0*(-Z/r+vh+vxc(:,sp))*psi_in(:,ish,sp)
 ! f=(-f4-f6+f5*dble(l*(l+1))/r**2*psi_in(:,ish,sp))+f7
 ! open(11,file='psi_f1_f2_f_1.dat',status='replace')
 ! write(11,*)"r psi -f4 -f6 f7"
 ! do ir=1, Ngrid
 !   write(11,*)r(ir),psi(ir,ish,sp),-f4(ir),-f6(ir),f7(ir)
 ! enddo
 !
 ! close(11)
 ! f=-f


 f1=v_rel*alpha2/(1d0-v_rel*alpha2)
 call rderivative_lagrN_st3(Ngrid,r,tools,tools_info,psi_in(:,ish,sp),f2)
 call rderivative_lagrN_st3(Ngrid,r,tools,tools_info,f1*f2*r**2,f3)
 f3=f3/r**2
 f=-f3+f1*dble(l*(l+1))/r**2*psi_in(:,ish,sp) + 2d0*(-Z/r+vh+vxc(:,sp))*psi_in(:,ish,sp) 
 f=-f
  endif
  call scrPoisson(Ngrid, r,tools,tools_info,l, f, eig(ish,sp), psi(:,ish,sp))


  call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,ish,sp)**2,norm)
  psi(:,ish,sp)=psi(:,ish,sp)/dsqrt(norm)

enddo



!Get non-local exchange
if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
   do inn=1,nmax
   ish=inn+shell0
   call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ(:,sp),lmax,&
           psi(:,ish,sp),psi_in(:,:,sp),vx_chi(:,ish,sp),vx_chi_sr(:,ish,sp),rsfunC,Nrsfun,hybx_w,Bess_ik)
   vx_chi(:,ish,sp)=vx_chi(:,ish,sp)*dble(Nspin)
   vx_chi_sr(:,ish,sp)=vx_chi_sr(:,ish,sp)*dble(Nspin)
   enddo

endif


eigp=eig




call orthonorm_get_eig(Ngrid,r,tools,tools_info,Z,l,nmax,relativity,v_rel,hybx_w,&
        vxc(:,sp),vh,vx_chi(:,shell0+1:shell0+nmax,sp),vx_chi_sr(:,shell0+1:shell0+nmax,sp),&
        psi(:,shell0+1:shell0+nmax,sp),eig(shell0+1:shell0+nmax,sp))



!Place for convergence check!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!recalculate potentials

!call get_local_exc_vxc_vh_rho(Ngrid,r,tools,tools_info,Nshell,shell_occ,spin,Nspin,psi,&
!         xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func,hybx_w,exc,vxc,vh,rho)
!
! ! Get non-local exchange
!
!if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
! do ish=1,Nshell
! do isp=1,Nspin !only one spin chanel must be updated????
! call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ(:,isp),lmax,psi(:,ish,isp),psi(:,:,isp),&
!           vx_psi(:,ish,isp),vx_psi_sr(:,ish,isp),rsfunC,Nrsfun,hybx_w,Bess_ik)
!   enddo
!   enddo
!vx_psi=vx_psi*dble(Nspin)
!vx_psi_sr=vx_psi_sr*dble(Nspin)
!
!endif !end get Fock exchange



enddo !self consistent loop
end subroutine

subroutine scrPoisson(Ngrid, r,tools,tools_info,l,f, e, psi)

integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: Ngrid
integer, intent(in) :: l
real(8), intent(in) :: r(Ngrid)
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
!        write(*,*)"scrPoisson Error: positive eigenvalue!"
        e=-1d-3
endif

lam=dsqrt(-2d0*e)


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
