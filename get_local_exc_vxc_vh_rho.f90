subroutine get_local_exc_vxc_vh_rho(Ngrid,r,tools,tools_info,Nshell,shell_occ,spin,Nspin,psi,&
         xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func,hybx_w,exc,vxc,vh,rho)
use xc_f03_lib_m
implicit none
integer(8), intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
integer, intent(in) :: Nshell
real(8), intent(in) :: shell_occ(Nshell,Nspin)
logical, intent(in) :: spin
integer, intent(in) :: Nspin
real(8), intent(in) :: psi(Ngrid,Nshell,Nspin) 
integer, intent(in) :: xc1_num,xc2_num,xc3_num
TYPE(xc_f03_func_t),intent(in)  :: xc1_func,xc2_func,xc3_func

real(8), intent(in) :: hybx_w(5,2)
real(8), intent(out) :: exc(Ngrid,Nspin)
real(8), intent(out) :: vxc(Ngrid,Nspin)
real(8), intent(out) :: vh(Ngrid)

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer :: isp, ish, ixc,ir

integer :: xc_num(3)
TYPE(xc_f03_func_t) :: xc_func(3)
TYPE(xc_f03_func_info_t) :: ixc_info

real(8) :: rho_sp(Nspin,Ngrid),rho(Ngrid)

real(8) :: vxc_i(Ngrid), exc_i(Ngrid)
real(8) :: vxc_sp_i(Nspin,Ngrid)

real(8) :: grho(Ngrid),g2rho(Ngrid),grho2(Ngrid),ftemp1(Ngrid),ftemp2(Ngrid),vxcsigma(Ngrid)

real(8) :: grho_sp(2,Ngrid),grho2_sp(3,Ngrid),vxc1_sp(2,Ngrid),vxc2_sp(2,Ngrid),vxc3_sp(2,Ngrid),&
        vxcsigma_sp(3,Ngrid),gvxcsigma_sp(3,Ngrid),g2rho_sp(2,Ngrid)

xc_num(1)=xc1_num
xc_num(2)=xc2_num
xc_num(3)=xc3_num
xc_func(1)=xc1_func
xc_func(2)=xc2_func
xc_func(3)=xc3_func


! Calculate density
  rho=0d0*r
if (spin)then
  rho_sp(1,:)=0d0*r !rho spin up
  rho_sp(2,:)=0d0*r !rho spin down
endif

do ish=1,Nshell
  if (.not.spin)then
    rho=rho+shell_occ(ish,1)*psi(:,ish,1)**2
  else
    rho_sp(1,:)=rho_sp(1,:)+shell_occ(ish,1)*psi(:,ish,1)**2
    rho_sp(2,:)=rho_sp(2,:)+shell_occ(ish,2)*psi(:,ish,2)**2
    rho=rho_sp(1,:)+rho_sp(2,:)
  endif
enddo




rho_sp=rho_sp/(4d0*Pi)
rho=rho/(4d0*Pi)


do isp=1, Nspin
vxc(:,isp)=0d0*r
exc(:,isp)=0d0*r
enddo


do ixc=1,3

vxc_i=0d0*r
exc_i=0d0*r
if (spin) then
vxc_sp_i(1,:)=0d0*r
vxc_sp_i(2,:)=0d0*r
endif
!! EXCHANGE-CORRELATION ixc  !!!

if (xc_num(ixc).gt.0)then
  ixc_info = xc_f03_func_get_info(xc_func(ixc))
  if (xc_f03_func_info_get_family(ixc_info).eq.XC_FAMILY_LDA) then
    if (.not.spin)then 
          call xc_f03_lda_exc_vxc(xc_func(ixc), Ngrid, rho(1), exc_i(1),vxc_i(1))
  else !spin polarised LDA
          call xc_f03_lda_exc_vxc(xc_func(ixc), Ngrid, rho_sp(1,1), exc_i(1),vxc_sp_i(1,1))
  endif
  elseif  ((xc_f03_func_info_get_family(ixc_info).eq.XC_FAMILY_GGA).or.&
                 (xc_f03_func_info_get_family(ixc_info).eq.XC_FAMILY_HYB_GGA)) then
   if (.not.spin)then
     call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
     call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
     g2rho=g2rho*r**(-2)
     grho2=grho**2
     call xc_f03_gga_exc_vxc(xc_func(ixc), Ngrid, rho(1), grho2(1), exc_i(1),vxc_i(1),vxcsigma(1))!

  !!!!! formula 6.0.8 !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vxcsigma,ftemp1)
!  ftemp1=ftemp1*r**(-2)
!  vxc1=vxc1-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!!

!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
      call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma,ftemp1)
     vxc_i=vxc_i-2d0*(grho*ftemp1+vxcsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!!
!
     else! spin polarised GGA

! xc_f03_gga_exc_vxc(func, Ngrid, rho(1), sigma(1), ex(1),vrho(1),vsigma(1)) --libxc documentation
!                     (in)  (in)   (in)     (in)    (out)  (out)   (out)
!                                          grho2       
     call rderivative_lagrN(Ngrid,r,tools,tools_info,rho_sp(1,:),grho_sp(1,:))
     call rderivative_lagrN(Ngrid,r,tools,tools_info,rho_sp(2,:),grho_sp(2,:))
     grho2_sp(1,:)=grho_sp(1,:)**2
     grho2_sp(2,:)=grho_sp(1,:)*grho_sp(2,:)
     grho2_sp(3,:)=grho_sp(2,:)**2
     call xc_f03_gga_exc_vxc(xc_func(ixc),Ngrid,rho_sp(1,1),grho2_sp(1,1)&
            ,exc_i(1),vxc_sp_i(1,1),vxcsigma_sp(1,1))
    
     call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho_sp(1,:),g2rho_sp(1,:))
     call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho_sp(2,:),g2rho_sp(2,:))
     g2rho_sp(1,:)=g2rho_sp(1,:)*r**(-2)
     g2rho_sp(2,:)=g2rho_sp(2,:)*r**(-2)

     call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma_sp(1,:),gvxcsigma_sp(1,:))
     call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma_sp(2,:),gvxcsigma_sp(2,:))
     call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma_sp(3,:),gvxcsigma_sp(3,:))
 
   ! v_xc for spin up is being calculated using this formula:
   ! $V_{xc}^{\uparrow}({\bf r})=\frac{\partial\hat{\epsilon}_{xc}}{\partial
   ! \rho^{\uparrow}({\bf r})}-2\left(\nabla\frac{\partial\hat{\epsilon}_{xc}}
   ! {\partial(\nabla\rho^{\uparrow})^2}\right)\cdot\nabla\rho^{\uparrow}
   ! -2\frac{\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow})^2}\nabla^2
   ! \rho^{\uparrow}-$\\
   ! $-\left(\nabla\frac{\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow}
   ! \cdot\nabla\rho^{\downarrow})}\right)\cdot\nabla\rho^{\downarrow}
   ! -\frac{\partial\hat{\epsilon}_{xc}}{\partial(\nabla\rho^{\uparrow}\cdot
   ! \nabla\rho^{\downarrow})}\nabla^2\rho^{\downarrow},$\\
   
     vxc_sp_i(1,:)=vxc_sp_i(1,:)-2d0*gvxcsigma_sp(1,:)*grho_sp(1,:)-2d0*vxcsigma_sp(1,:)*g2rho_sp(1,:)-&
             gvxcsigma_sp(2,:)*grho_sp(2,:)-vxcsigma_sp(2,:)*g2rho_sp(2,:)

     vxc_sp_i(2,:)=vxc_sp_i(2,:)-2d0*gvxcsigma_sp(3,:)*grho_sp(2,:)-2d0*vxcsigma_sp(3,:)*g2rho_sp(2,:)-&
             gvxcsigma_sp(2,:)*grho_sp(1,:)-vxcsigma_sp(2,:)*g2rho_sp(1,:)
!
!!     write(*,*)"Get spin polarized potential"
!!     do i=1, Ngrid
!!     write(*,*)vxc1_sp(1,i),vxc1_sp(2,i)
!!     enddo
!!    stop
!
     endif !spin 
 
  else
    write(*,*)"1-st XC functional not supported!"
  endif
!elseif (xc1_num.eq.-1) then !Built-in LDA
!   call getvxc(Ngrid,rho,vxc1,exc1)
!elseif (xc1_num.eq.-2) then !Built-in GGA
!  call getxc_pbe(Ngrid,r,tools,tools_info,rho,vxc1,exc1)
endif
!! END EXCHANGE-CORRELATION 1 !!!


exc(:,1)=exc(:,1)+hybx_w(ixc,1)*exc_i

if (.not.spin) then

vxc(:,1)=vxc(:,1)+hybx_w(ixc,1)*vxc_i

else
      vxc(:,1)=vxc(:,1)+hybx_w(ixc,1)*vxc_sp_i(1,:)
      vxc(:,2)=vxc(:,2)+hybx_w(ixc,1)*vxc_sp_i(2,:)
      
endif

enddo ! ixc

rho=rho*(4d0*Pi)
rho_sp=rho_sp*(4d0*Pi)

! Get Coulomb potential

call  integ_BodesN_fun(Ngrid,r,tools,tools_info,1,r**2*rho,ftemp1)
call  integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,r*rho,ftemp2)
vh=ftemp1/r+ftemp2


end subroutine




