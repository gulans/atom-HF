subroutine getxc_pbe(Ngrid,r,tools,tools_info,rho,vxc,exc)
implicit none

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(inout) :: r(Ngrid),rho(Ngrid) 
real(8), intent(out) :: vxc(Ngrid),exc(Ngrid)

integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))


integer :: i

real(8) :: kappa,mu,beta,g2rho(Ngrid),test,t11
real(8) :: rhoup(Ngrid),rhodn(Ngrid),grho(Ngrid),gup(Ngrid),gdn(Ngrid),g2up(Ngrid),g2dn(Ngrid),g3rho(Ngrid),g3up(Ngrid),&
     g3dn(Ngrid),ex(Ngrid),ec(Ngrid),vxup(Ngrid),vxdn(Ngrid),vcup(Ngrid),vcdn(Ngrid),temp(Ngrid)

!rho=exp(-r)

kappa=0.80400000000000005d0     
mu=0.21951497276451709d0     
beta=6.6724550603149219d-002

!  call rderivative_lagrN(Ngrid,r,tools,tools_info,-rho,grho)
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,grho*r**2,g2rho)
  call rderivative_lagr3(Ngrid,r,-rho,grho)
  call rderivative_lagr3(Ngrid,r,grho*r**2,g2rho)

g2rho=-g2rho/r**2
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,abs(grho),g3rho)

  call rderivative_lagr3(Ngrid,r,abs(grho),g3rho)

  g3rho=-grho*g3rho

!write(*,*)"r rho grho g2rho g3rho"
  do i=1,Ngrid 
!write(*,*)r(i),",",rho(i),",",grho(i),",",g2rho(i),",",g3rho(i)
  enddo
!stop


open(11,file='rho.out',status='replace')
  do i=1,Ngrid
    write(11,*)r(i),",",rho(i),",",grho(i),",",g2rho(i),",",g3rho(i)
  enddo
close(11)

!stop



! call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
! call rderivative_lagrN(Ngrid,r,tools,tools_info,grho,g2rho)
! g2rho=2d0*grho/r+g2rho
! g3rho=grho*g2rho

  
           !Call fderiv (1, nr, r, rho, grho, cf)
! grad^2 rho
           ! Call fderiv (2, nr, r, rho, g2rho, cf)
           ! Do ir = 1, nr
           !    g2rho (ir) = g2rho (ir) + 2.d0 * ri (ir) * grho (ir)
           ! End Do
! approximate (grad rho).(grad |grad rho|)
           ! Do ir = 1, nr
           !    g3rho (ir) = grho (ir) * g2rho (ir)
           ! End Do



!rhoup=rho
!rhodn=rhoup*0d0
!gup=grho
!gdn=gup*0d0
!g2up=g2rho
!g2dn=g2up*0d0
!g3up=g3rho
!g3dn=g3up*0d0

!rho(1)=1d0
!grho(1)=-1d0
!g2rho(1)=-200000d0
!g3rho(1)=1d0


! open (2, file = 'rho_vxc.dat', status = 'old')
!   read(2,*)
!   do i = 1,Ngrid
!      read(2,*) t11, rho(i),grho(i),g2rho(i),g3rho(i)
!   end do

 !  close(2)


rhoup=rho/2d0
rhodn=rhoup
gup=grho/2d0
gdn=gup
g2up=g2rho/2d0
g2dn=g2up
g3up=g3rho/2d0
g3dn=g3up



call xc_pbe(Ngrid,kappa,mu,beta,rhoup,rhodn,grho,gup,gdn,g2up,g2dn,g3rho,g3up, &
     g3dn,ex,ec,vxup,vxdn,vcup,vcdn)
exc=ex+ec
vxc=vxup+vcup!+vxdn+vcdn
!write(*,*)"r(i) vxup(i) vcup(i)"

open(11,file='rho_vxc_atom.dat',status='replace')
write(11,*)"r(i) rho(i) grho(i) g2rho(i) g3rho(i) vxup(i) vcup(i)"
do i=1,Ngrid
write(11,*)r(i),rho(i),grho(i),g2rho(i),g3rho(i),vxup(i),vcup(i),vxup(i)+vcup(i)
enddo
close(11)

do i=1,Ngrid 
! write(*,*)r(i),",",vxup(i),",",vcup(i)
enddo
! write(*,*)"r(i),rho(i),grho(i),g2rho(i),g3rho(i)"

do i=1,Ngrid
! write(*,*)r(i),",",rho(i),",",grho(i),",",g2rho(i),",",g3rho(i)
enddo


!stop

end subroutine
