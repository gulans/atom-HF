subroutine LS_iteration(Ngrid, r,tools,tools_info,rsfunC,Nrsfun,hybx_w, vn,l,sp,shell_l,shell_occ, nmax,&
                shell0, Nshell,Nspin,relativity,lmax,vxc,v_rel, vh,vx_psi,vx_psi_sr, psi_in,psi_in_p,psi,eig,Bess_ik,&
                iner_loop,xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func,F_mix)
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
real(8), intent(in) :: r(Ngrid),vn(Ngrid),vxc(Ngrid,Nspin),vh(Ngrid),shell_occ(Nshell,Nspin)
real(8), intent(in) :: psi_in(Ngrid,Nshell,Nspin),psi_in_p(Ngrid,Nshell,Nspin),F_mix
real(8), intent(inout) :: vx_psi(Ngrid,nmax),vx_psi_sr(Ngrid,nmax)

logical, intent(in) :: relativity
real(8), intent(in) :: v_rel(Ngrid)
complex(8), intent(in) ::Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
real(8), intent(inout) :: eig(nmax)
integer, intent(out) :: iner_loop
real(8), intent(out) :: psi(Ngrid,nmax)
integer, intent(in) ::xc1_num,xc2_num,xc3_num
TYPE(xc_f03_func_t) :: xc1_func,xc2_func,xc3_func

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)
integer :: inn,inp,ish,i,j 
real(8) :: vx_chi(Ngrid,nmax),vx_chi_sr(Ngrid,nmax)
real(8) :: vx_chi2(Ngrid,nmax),vx_chi2_sr(Ngrid,nmax)
real(8) :: f1(Ngrid),f2(Ngrid),f3(Ngrid),f4(Ngrid),f5(Ngrid)
real(8) :: f6(Ngrid),f7(Ngrid)
real(8) :: phi(Ngrid,nmax),norm,eigp(nmax)
integer :: iscl,maxscl,ir,isp
real(8) :: f(Ngrid)
logical :: spin

real(8) :: elimit, eshift
real(8) :: vx_psi_in(Ngrid,nmax),vx_psi_sr_in(Ngrid,nmax)
logical :: new_algorithm
logical :: eig_limiter

new_algorithm=.false.
eig_limiter=.false.
elimit=-0.0002d0

vx_psi_in=vx_psi
vx_psi_sr_in=vx_psi_sr

vx_chi=psi*0d0
vx_chi_sr=psi*0d0

if (Nspin.eq.2)then
        spin=.true.
else
        spin=.false.
endif


maxscl=10
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

          
if (new_algorithm)then
  f=-2d0*( (vn+vh+vxc(:,sp))*psi(:,inn) + hybx_w(4,1)*vx_psi(:,inn)&
          + hybx_w(5,1)*vx_psi_sr(:,inn) )
else!old version
  f=-2d0*( (vn+vh+vxc(:,sp))*psi_in(:,ish,sp) + hybx_w(4,1)*vx_psi_in(:,inn)&
          + hybx_w(5,1)*vx_psi_sr_in(:,inn) )
endif




  else

if(.true.)then
    f1=1d0/(1d0-alpha2*v_rel)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,log(f1),f2)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi_in(:,ish,sp),f3)
    f4=-f2*f3
    
    f2=2d0*(vn+vh+vxc(:,sp))*psi_in(:,ish,sp)/f1
    f3=2d0*alpha2*v_rel*eig(inn)*psi_in(:,ish,sp)
    f=f4+f2+f3
    f=-f
else
    f1=v_rel*alpha2/(1d0-v_rel*alpha2)
    ! call rderivative_lagrN_st3(Ngrid,r,tools,tools_info,psi_in(:,ish,sp),f2)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi_in(:,ish,sp),f2)
    ! call rderivative_lagrN_st3(Ngrid,r,tools,tools_info,f1*f2*r**2,f3)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,f1*f2*r**2,f3)
    f3=f3/r**2
    f=-f3+f1*dble(l*(l+1))/r**2*psi_in(:,ish,sp) + 2d0*(vn+vh+vxc(:,sp))*psi_in(:,ish,sp) 
    f=-f
endif



  endif




elimit=-2d0
if ((eig_limiter).and.(eig(inn).gt.elimit))then
  eshift=-eig(inn)+elimit
  eig(inn)=elimit
  if(new_algorithm)then
    call scrPoisson(Ngrid, r,tools,tools_info,l, f-2*(eshift)*psi(:,inn), eig(inn), psi(:,inn))
  else
    call scrPoisson(Ngrid, r,tools,tools_info,l, f-2*(eshift)*psi_in(:,ish,sp), eig(inn), psi(:,inn))
  endif    
  eig(inn)=eig(inn)-eshift

else
  call scrPoisson(Ngrid, r,tools,tools_info,l, f, eig(inn), psi(:,inn))
endif






  call integ_BodesN_value(Ngrid,r,tools,tools_info,r**2*psi(:,inn)**2,norm)
  psi(:,inn)=psi(:,inn)/dsqrt(norm)

enddo



!Get non-local exchange
if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
   do inn=1,nmax
   ish=inn+shell0
   call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ(:,sp),lmax,&
           psi(:,inn),psi_in(:,:,sp),vx_chi(:,inn),vx_chi_sr(:,inn),rsfunC,Nrsfun,hybx_w,Bess_ik)

   if (F_mix.lt.0.999999999d0) then
   call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ(:,sp),lmax,&
           psi(:,inn),psi_in_p(:,:,sp),vx_chi2(:,inn),vx_chi2_sr(:,inn),rsfunC,Nrsfun,hybx_w,Bess_ik)
   vx_chi(:,inn)=F_mix*vx_chi(:,inn)+(1d0-F_mix)*vx_chi2(:,inn)
   vx_chi_sr(:,inn)=F_mix*vx_chi_sr(:,inn)+(1d0-F_mix)*vx_chi2_sr(:,inn)

   endif


   vx_chi(:,inn)=vx_chi(:,inn)*dble(Nspin)
   vx_chi_sr(:,inn)=vx_chi_sr(:,inn)*dble(Nspin)




   enddo

endif


eigp=eig




call orthonorm_get_eig(Ngrid,r,tools,tools_info,vn,l,nmax,relativity,v_rel,hybx_w,&
        vxc(:,sp),vh,vx_chi,vx_chi_sr,&
        psi,eig,&
        vx_psi,vx_psi_sr)


!if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
! do ish=shell0+1,shell0+nmax
! call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ(:,sp),lmax,psi(:,ish,sp),psi_in(:,:,sp),&
!           vx_psi(:,ish,sp),vx_psi_sr(:,ish,sp),rsfunC,Nrsfun,hybx_w,Bess_ik)
!   enddo
!vx_psi=vx_psi*dble(Nspin)
!vx_psi_sr=vx_psi_sr*dble(Nspin)
!endif !end get Fock exchange


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
integer :: ri,i,endpoint
real(8) :: f11(Ngrid),f12(Ngrid),f21(Ngrid),f22(Ngrid),int1(Ngrid),int2(Ngrid)
!result=f11*(integral_(0->r)f12)+f21(integral_(r->Ngrid)f22)
real(8) :: besi,besk


!write(*,*)"e=",e
if (e.gt.0) then
!        write(*,*)"scrPoisson Error: positive eigenvalue!"
        e=-2d-2
endif

lam=dsqrt(-2d0*e)


f11=0d0*r
f12=0d0*r
f21=0d0*r
f22=0d0*r
!write(*,*)"argument       besi         besk"
endpoint=Ngrid
do ri=1, Ngrid
if((lam*r(ri)).gt.200d0) then
endpoint=ri
exit
endif
call msbesseli(l,lam*r(ri), besrezi)
call msbesselk(l,lam*r(ri), besrezk)


besi=besrezi(l)
besk=besrezk(l)
!write(*,*)lam*r(ri),besi,besk

f11(ri)=lam*besk
f12(ri)=besi*f(ri)
f21(ri)=lam*besi
f22(ri)=besk*f(ri)
enddo

call integ_BodesN_fun(Ngrid,r,tools,tools_info,1,f12*r**2,int1)
!call integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,f22*r**2,int2)
call integ_BodesN_fun(endpoint,r(1:endpoint),tools(1:endpoint,:),tools_info,-1,f22(1:endpoint)*r(1:endpoint)**2,int2(1:endpoint))
int2(endpoint+1:Ngrid)=0d0*r(endpoint+1:Ngrid)
!do ri=1, Ngrid
!write(*,*)r(ri),int2(ri)
!enddo
!stop
psi=f11*int1+f21*int2

! open(11,file='scr_poisson_out1.dat',status='replace')
!
!  write(11,*)"r,lam*r,f,f11, int1,f11*int1 , f21, int2,f21*int2, psi'"
!  do ir=1,Ngrid
!  write(11,*) r(ir),",",lam*r(ir),",",f(ir),",",f11(ir),",",int1(ir),",",f11(ir)*int1(ir)&
!          ,",",f21(ir),",",int2(ir),",",f21(ir)*int2(ir),",",psi(ir)
!
!  enddo
!  close(11)




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
