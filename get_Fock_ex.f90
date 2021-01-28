!call subroutine get_Fock_ex(Ngrid,r,ish,Nshell,shell_l,psi_non_norm(:,ish),psi_all,vx_u)

subroutine get_Fock_ex(Ngrid,r,shell,Nshell,shell_l,psi,psi_all,vx_psi)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0

integer, intent(in)  :: Nshell,Ngrid,shell,shell_l(Nshell)
real(8), intent(in)  :: psi_all(Ngrid,Nshell),r(Ngrid),psi(Ngrid)
real(8), intent(out) :: vx_psi(Ngrid)


integer :: ish,l,lpri,lpripri
real(8) :: gc,u_all(Ngrid,Nshell),u(Ngrid)
real(8) :: integ1(Ngrid),integ2(Ngrid),integ(Ngrid)

u=psi*r

do ish=1, Nshell
  u_all(:,ish)=psi_all(:,ish)*r
enddo
l=shell_l(shell)
integ=0d0*integ
do ish=1, Nshell
  lpri=shell_l(ish)


  do lpripri=abs(l-lpri),l+lpri,1
    
  
    call wigner3j(l,lpri,lpripri,gc)
    gc=dble(2*lpri+1)*gc
    
    write(*,*)"(l,l',l'') (",l,",",lpri,",",lpripri,")", " Gaunt_coef=",gc
    if (gc.ne.0d0) then
      call integ_s38_fun(Ngrid,r,   u_all(:,ish)*u*r**lpripri    ,1,integ1)
      integ1=integ1/r**(lpripri+1)
      call integ_s38_fun(Ngrid,r,   u_all(:,ish)*u/r**(lpripri+1),-1,integ2)
      integ2=integ2*r**lpripri
      integ=integ+gc*u_all(:,ish)*(-integ1-integ2)
    endif
  enddo
enddo
vx_psi=integ/r
end subroutine
