subroutine get_Fock_ex(Ngrid,r,tools,tools_info,shell,Nshell,shell_l,lmax,psi,psi_all,&
                vx_psi,vx_psi_sr,vx_psi_lr,rs,rsfunC,Nrsfun,hybx_w,Bess_ik)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: tools_info(3),lmax
real(8), intent(in) :: tools(Ngrid,tools_info(1)),hybx_w(5,2)
complex(8), intent(in) ::rsfunC(Nrsfun,2)
integer, intent(in)  :: Nshell,Ngrid,shell,shell_l(Nshell),rs,Nrsfun
real(8), intent(in)  :: psi_all(Ngrid,Nshell),r(Ngrid),psi(Ngrid)
complex(8), intent(in)  :: Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
real(8), intent(out) :: vx_psi(Ngrid),vx_psi_sr(Ngrid),vx_psi_lr(Ngrid)

complex(8) ::rez(Ngrid)
complex(8) :: integc1(Ngrid),integc2(Ngrid)
integer :: k,i,ir,ish,l,lpri,lpripri
real(8) :: gc,u_all(Ngrid,Nshell),u(Ngrid)
real(8) :: integ1(Ngrid),integ2(Ngrid),integ3(Ngrid),integ(Ngrid),integ_sr(Ngrid)
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
      call integ_BodesN_fun(Ngrid,r,tools,tools_info,-1, u_all(:,ish)*u/r**(lpripri+1),integ2)


      integ2=integ2*r**lpripri
      integ3=-integ1-integ2
!      if ((lpripri.eq.2).and.(lpri.eq.1).and.(l.eq.1)) then
!      open(12,file='Fock_C_int_pp2.dat',status='replace')
!         do i=1, Ngrid
!
!      write(12,*)r(i),integ3(i)
!      enddo
!      close(12)
!      !stop 
!     endif
      integ=integ+gc*u_all(:,ish)*integ3

if (rs.eq.1) then !! range seperation short-range component (erfc)/r
!write(*,*)"Foka funkcija ",Nrsfun

rez=cmplx(0d0,0d0,8)*r
do k=1, Nrsfun

!call integC_BodesN_fun(Ngrid,r, tools, tools_info,1,  conjg(Bess_ik(:,k,lpripri+1,1)) *u_all(:,ish)*u    ,integc1)
!integc1=integc1*conjg(rsfunC(k,2))*conjg(Bess_ik(:,k,lpripri+1,2))
!call integC_BodesN_fun(Ngrid,r, tools, tools_info,-1, conjg(Bess_ik(:,k,lpripri+1,2)) *u_all(:,ish)*u    ,integc2)
!integc2=integc2*conjg(rsfunC(k,2))*conjg(Bess_ik(:,k,lpripri+1,1))
!rez=rez+conjg(rsfunC(k,1))*(integc1+integc2)

!varbūt var saīināt ja integrāļi arī ir kompleksi saistīti



call integC_BodesN_fun(Ngrid,r, tools, tools_info,1,  Bess_ik(:,k,lpripri+1,1) *u_all(:,ish)*u    ,integc1)
integc1=integc1*rsfunC(k,2)*Bess_ik(:,k,lpripri+1,2)
call integC_BodesN_fun(Ngrid,r, tools, tools_info,-1, Bess_ik(:,k,lpripri+1,2) *u_all(:,ish)*u    ,integc2)
integc2=integc2*rsfunC(k,2)*Bess_ik(:,k,lpripri+1,1)
rez=rez+rsfunC(k,1)*(integc1+integc2)

enddo


   do i=1, Ngrid
     integ1(i)=-2d0*realpart(rez(i))
   enddo
    integ1=dble(2*lpripri+1)*integ1

   
!   if ((lpripri.eq.2).and.(lpri.eq.1).and.(l.eq.1)) then
!      open(12,file='Fock_RS_int_pp2.dat',status='replace')
!         do i=1, Ngrid
!      write(12,*)r(i),integ1(i)
!      enddo
!      close(12)
!      stop
!      endif




   integ_sr=integ_sr+gc*u_all(:,ish)*integ1



endif !(rs.eq.1)
    endif !(gc.ne.0d0)
  enddo !lpripri
enddo !ish
vx_psi=integ/r
vx_psi_sr=integ_sr/r



!vx_psi=vx_psi_sr



end subroutine
