subroutine get_energy(Ngrid,r,tools,tools_info,vn,Nshell,shell_occ,spin,relativity,v_rel,Nspin,shell_l,hybx_w,&
        vxc,exc,rho,vh,vx_psi,vx_psi_sr,psi,&
        e_kin,e_ext,e_h,e_xc)
use modinteg
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
real(8), intent(in) :: vn(Ngrid)
integer, intent(in) :: Nshell, Nspin
real(8), intent(in) :: shell_occ(Nshell,Nspin)
logical, intent(in) :: spin
logical, intent(in) :: relativity
real(8), intent(in) :: v_rel(Ngrid,Nspin)
integer, intent(in) :: shell_l(Nshell)
real(8), intent(in) :: hybx_w(5,2)
real(8), intent(in) :: vxc(Ngrid,Nspin)
real(8), intent(in) :: exc(Ngrid,Nspin)
real(8), intent(in) :: rho(Ngrid)
real(8), intent(in) :: vh(Ngrid)
real(8), intent(in) :: vx_psi(Ngrid,Nshell,Nspin)
real(8), intent(in) :: vx_psi_sr(Ngrid,Nshell,Nspin)
real(8), intent(in) :: psi(Ngrid,Nshell,Nspin) 
real(8), intent(out) :: e_kin,e_ext,e_h,e_xc


real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)

integer :: ish,isp
real(8) :: e1,e2,e3
real(8) :: ftemp1(Ngrid),ftemp2(Ngrid)

e_kin=0d0
if (.not.relativity) then
do ish=1,Nshell
do isp=1, Nspin
  !call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,ish,isp),ftemp1)
  call deriv_f(Ngrid,r,psi(:,ish,isp),ftemp1)
  call integ_v(Ngrid,r,0.5d0*ftemp1**2*r**2,e1)
  call integ_v(Ngrid,r,0.5d0*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish,isp)**2,e2)
 ! call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*ftemp1**2*r**2,e1)
 ! call integ_BodesN_value(Ngrid,r,tools, tools_info,0.5d0*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish,isp)**2,e2)
  e_kin=e_kin+shell_occ(ish,isp)*(e1+e2)
!  write(*,*)"e1, e2",e1,e2
enddo
enddo
else !relativity
do ish=1,Nshell
do isp=1, Nspin

!    ftemp1=1d0/(1d0-v_rel(:,isp)*alpha2)
!    call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,ish,isp),ftemp2)
!    call rderivative_lagrN(Ngrid,r,tools,tools_info,ftemp1*ftemp2*r**2,ftemp3)
!    ftemp4=ftemp3/r**2
!    call integ_BodesN_value(Ngrid,r,tools,tools_info,-0.5d0*psi(:,ish,isp)*ftemp4*r**2,e1)

    ftemp2=1d0/(1d0-v_rel(:,isp)*alpha2)
   !call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,ish,isp),ftemp1)
   call deriv_f(Ngrid,r,psi(:,ish,isp),ftemp1) 
   call integ_v(Ngrid,r,0.5d0*ftemp2*ftemp1**2*r**2,e1)
   !call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*ftemp2*ftemp1**2*r**2,e1)

    call integ_v(Ngrid,r,0.5d0*ftemp2*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish,isp)**2,e2)
   !call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*ftemp2*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish,isp)**2,e2)
    e_kin=e_kin+shell_occ(ish,isp)*(e1+e2)
!    write(*,*)"e1,e2",e1,e2
enddo
enddo

endif
call integ_v(Ngrid,r, r**2*rho*vn,e_ext)
!call integ_BodesN_value(Ngrid,r,tools,tools_info, r**2*rho*vn,e_ext)

call integ_v(Ngrid,r,r**2*0.5d0*rho*vh,e_h)
!call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*rho*vh,e_h)

e2=0d0
  do isp=1,Nspin
  do ish=1,Nshell
    call integ_v(Ngrid,r,r**2*0.5d0*shell_occ(ish,isp)*psi(:,ish,isp)*&
            (hybx_w(4,1)*vx_psi(:,ish,isp)+hybx_w(5,1)*vx_psi_sr(:,ish,isp)),e1)

    !call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*shell_occ(ish,isp)*psi(:,ish,isp)*&
    !        (hybx_w(4,1)*vx_psi(:,ish,isp)+hybx_w(5,1)*vx_psi_sr(:,ish,isp)),e1)
    e2=e2+e1
  enddo
  enddo
  call integ_v(Ngrid,r,exc(:,1)*rho*r**2,e3)

  !call integ_BodesN_value(Ngrid,r,tools, tools_info,exc(:,1)*rho*r**2,e3)
 e_xc=e2+e3

end subroutine




