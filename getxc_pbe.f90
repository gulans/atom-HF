subroutine getxc_pbe(Ngrid,r,tools,tools_info,rho,vxc,exc) !,vx,vc,ex,ec)
implicit none

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(inout) :: r(Ngrid),rho(Ngrid) 
real(8), intent(out) :: vxc(Ngrid),exc(Ngrid)
real(8) :: vx(Ngrid),vc(Ngrid),ex(Ngrid),ec(Ngrid)

integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))


integer :: i

real(8) :: kappa,mu,beta,g2rho(Ngrid),test,t11
real(8) :: rhoup(Ngrid),rhodn(Ngrid),grho(Ngrid),gup(Ngrid),gdn(Ngrid),g2up(Ngrid),g2dn(Ngrid),g3rho(Ngrid),g3up(Ngrid),&
     g3dn(Ngrid),vxup(Ngrid),vxdn(Ngrid),vcup(Ngrid),vcdn(Ngrid),temp(Ngrid)

!rho=exp(-r)

kappa=0.80400000000000005d0     
mu=0.21951497276451709d0     
beta=6.6724550603149219d-002

  call rderivative_lagrN(Ngrid,r,tools,tools_info,-rho,grho)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,grho*r**2,g2rho)

g2rho=-g2rho/r**2
 call rderivative_lagrN(Ngrid,r,tools,tools_info,abs(grho),g3rho)


  g3rho=-grho*g3rho



rhoup=rho/2d0
rhodn=rhoup
gup=grho/2d0
gdn=gup
g2up=g2rho/2d0
g2dn=g2up
g3up=g3rho/4d0
g3dn=g3up



call xc_pbe(Ngrid,kappa,mu,beta,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up, &
     g3dn,ex,ec,vxup,vxdn,vcup,vcdn)
exc=ex+ec
vxc=vxup+vcup!+vxdn+vcdn
vx=vxup
vc=vcup
!write(*,*)'rho grho ex ec vx vc'
!do i=1,40
! write(*,*)rho(i),grho(i),ex(i),ec(i),vx(i),vc(i)
!enddo
!stop
end subroutine
