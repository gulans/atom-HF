subroutine get_Fock_ex(Ngrid,r,tools,tools_info,shell,Nshell,shell_l,lmax,psi,psi_all,&
                vx_psi,vx_psi_sr,vx_psi_lr,rs,rsfunC,Nrsfun,hybx_w,Bess_ik)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: tools_info(3),lmax
real(8), intent(in) :: tools(Ngrid,tools_info(1)),hybx_w(4)
complex(8), intent(in) ::rsfunC(Nrsfun+1,2) !Nrsfun+1-th element real part is rmax and imaginary is parameter mu
integer, intent(in)  :: Nshell,Ngrid,shell,shell_l(Nshell),rs,Nrsfun
real(8), intent(in)  :: psi_all(Ngrid,Nshell),r(Ngrid),psi(Ngrid)
complex(8), intent(in)  :: Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
real(8), intent(out) :: vx_psi(Ngrid),vx_psi_sr(Ngrid),vx_psi_lr(Ngrid)

complex(8) ::rez(Ngrid)
complex(8) :: integc1(Ngrid),integc2(Ngrid)
integer :: k,i,ir,ish,l,lpri,lpripri
real(8) :: gc,u_all(Ngrid,Nshell),u(Ngrid)
real(8) :: integ1(Ngrid),integ2(Ngrid),integ(Ngrid),integ_sr(Ngrid)
real(8) :: rsmu,rsrmax




u=psi*r

do ish=1, Nshell
  u_all(:,ish)=psi_all(:,ish)*r
enddo
l=shell_l(shell)
integ=0d0*r
integ_sr=0d0*r
do ish=1, Nshell
  lpri=shell_l(ish)

!  write(*,*)"l'' sum range",abs(l-lpri), l+lpri, "step 2"

  do lpripri=abs(l-lpri),l+lpri,2
    
  
    call wigner3j(l,lpri,lpripri,gc)
    gc=dble(2*lpri+1)*gc**2  !(2*lpri+1) should be replaced with shell_occ 
!    write(*,*)"(l,l',l'') (",l,",",lpri,",",lpripri,")", " Gaunt_coef=",gc
    if (gc.ne.0d0) then
!!without range seperation (Coulumb)
    call integ_BodesN_fun(Ngrid,r, tools, tools_info,1,  u_all(:,ish)*u*r**lpripri    ,integ1)

      integ1=integ1/r**(lpripri+1)
!      call integ_s38_fun(Ngrid,r,   u_all(:,ish)*u/r**(lpripri+1),-1,integ2)
!      call integ_Bodes6_fun(Ngrid,r,   u_all(:,ish)*u/r**(lpripri+1),-1,integ2)
      call integ_BodesN_fun(Ngrid,r,tools,tools_info,-1, u_all(:,ish)*u/r**(lpripri+1),integ2)


      integ2=integ2*r**lpripri
      integ=integ+gc*u_all(:,ish)*(-integ1-integ2)

if (rs.eq.1) then !! range seperation short-range component (erfc)/r
!write(*,*)"Foka funkcija ",Nrsfun

rez=cmplx(0d0,0d0,8)*r
do k=1, Nrsfun

!call integC_BodesN_fun(Ngrid,r, tools, tools_info,1,  besi(:,-k)*u_all(:,ish)*u    ,integc1)
!integc1=integc1*conjg(rsfunC(k,2))*besk(:,-k)
!call integC_BodesN_fun(Ngrid,r, tools, tools_info,-1,  besk(:,-k)*u_all(:,ish)*u    ,integc2)
!integc2=integc2*conjg(rsfunC(k,2))*besi(:,-k)
!rez=rez+conjg(rsfunC(k,1))*(integc1+integc2)

!varbūt var saīināt ja integrāļi arī ir kompleksi saistīti



call integC_BodesN_fun(Ngrid,r, tools, tools_info,1,  Bess_ik(:,k,lpripri+1,1) *u_all(:,ish)*u    ,integc1)
integc1=integc1*rsfunC(k,2)*Bess_ik(:,k,lpripri+1,2)
call integC_BodesN_fun(Ngrid,r, tools, tools_info,-1, Bess_ik(:,k,lpripri+1,2) *u_all(:,ish)*u    ,integc2)
integc2=integc2*rsfunC(k,2)*Bess_ik(:,k,lpripri+1,1)
rez=rez+rsfunC(k,1)*(integc1+integc2)

enddo

   rsmu=hybx_w(3)
   !rsmu=realpart(rsfunc(Nrsfun+1,1))
   rsrmax=imagpart(rsfunc(Nrsfun+1,1))

   open(11,file='Fock_integral_test.dat',status='replace')
   write(11,*)"Foka apmaiņas integrālis erfc RE(rez) Im(rez) erfc-RE(rez) mu=",rsmu
   do i=1, Ngrid
     write(11,*)r(i), realpart(rez(i)), imagpart(rez(i))
   enddo
   close(11)
!stop
   do i=1, Ngrid
     integ1(i)=2d0*realpart(rez(i))
   enddo

   integ_sr=integ_sr+gc*u_all(:,ish)*(-integ1)



endif
    endif
  enddo
enddo
vx_psi=integ/r
vx_psi_sr=integ_sr/r
vx_psi_lr=vx_psi-vx_psi_sr

!vx_psi=vx_psi_lr

!   open(11,file='Fock_vxpsi.dat',status='replace')
!   
!   write(11,*)"Foka apmaiņas integrālis erfc RE(rez) Im(rez) erfc-RE(rez) mu=",rsmu
!   do i=1, Ngrid
!     write(11,*)r(i), vx_psi(i), vx_psi_sr(i),vx_psi_lr(i)
!   enddo
!   close(11)
!stop




end subroutine
