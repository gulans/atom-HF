subroutine get_Fock_ex(Ngrid,r,tools,tools_info,shell,Nshell,shell_l,psi,psi_all,vx_psi)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))

integer, intent(in)  :: Nshell,Ngrid,shell,shell_l(Nshell)
real(8), intent(in)  :: psi_all(Ngrid,Nshell),r(Ngrid),psi(Ngrid)
real(8), intent(out) :: vx_psi(Ngrid)


integer :: ir,ish,l,lpri,lpripri
real(8) :: gc,u_all(Ngrid,Nshell),u(Ngrid)
real(8) :: integ1(Ngrid),integ2(Ngrid),integ(Ngrid)

u=psi*r

do ish=1, Nshell
  u_all(:,ish)=psi_all(:,ish)*r
enddo
l=shell_l(shell)
integ=0d0*r
do ish=1, Nshell
  lpri=shell_l(ish)

!  write(*,*)"l'' sum range",abs(l-lpri), l+lpri, "step 2"

  do lpripri=abs(l-lpri),l+lpri,2
    
  
    call wigner3j(l,lpri,lpripri,gc)
    gc=dble(2*lpri+1)*gc**2  !(2*lpri+1) should be replaced with shell_occ 
!    write(*,*)"(l,l',l'') (",l,",",lpri,",",lpripri,")", " Gaunt_coef=",gc
    if (gc.ne.0d0) then
!      call integ_s38_fun(Ngrid,r,   u_all(:,ish)*u*r**lpripri    ,1,integ1)
!      call integ_Bodes6_fun(Ngrid,r,   u_all(:,ish)*u*r**lpripri    ,1,integ1)
     call integ_BodesN_fun(Ngrid,r, tools, tools_info,1,  u_all(:,ish)*u*r**lpripri    ,integ1)

      integ1=integ1/r**(lpripri+1)
!      call integ_s38_fun(Ngrid,r,   u_all(:,ish)*u/r**(lpripri+1),-1,integ2)
!      call integ_Bodes6_fun(Ngrid,r,   u_all(:,ish)*u/r**(lpripri+1),-1,integ2)
      call integ_BodesN_fun(Ngrid,r,tools,tools_info,-1, u_all(:,ish)*u/r**(lpripri+1),integ2)


      integ2=integ2*r**lpripri
      integ=integ+gc*u_all(:,ish)*(-integ1-integ2)
    endif
  enddo
enddo

vx_psi=integ/r

end subroutine
