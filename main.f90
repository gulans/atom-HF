program atomHF
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:),rho(:)
!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl = 20
integer :: ir, iscl

real(8) :: Z
real(8) :: Rmax,hh,e1,e2,e3,energy
integer :: Ngrid



!read input
read(*,*) 
read(*,*) Z, Rmax, Ngrid

allocate(r(Ngrid),vfull(Ngrid),vh(Ngrid),vxc(Ngrid),exc(Ngrid),H(Ngrid,Ngrid),eig(Ngrid),psi(Ngrid),rho(Ngrid))

call gengrid(Ngrid,Rmax,r)
hh=r(2)-r(1)


! Initial guess for potential
!...
do ir=1,Ngrid
  vfull(ir)=-Z/r(ir)
enddo


! ----- version 1 -------
! Self-consistency loop

do iscl=1,maxscl

! Construct Hamiltonian 

  call genham(Ngrid,r,vfull,H)


! Diagonalize
call diasym(Ngrid,H,eig)

#ifdef debug
write(*,*) "eig:"
   do ir=1,4
     write(*,*) eig(ir)
   enddo
#endif

! Normalize WFs
psi=H(:,1)  !tas ir Psi*r
psi=psi/r
call wfnorm(Ngrid,r,psi)
#ifdef debug
  Print *,"PSI1 NORMĒETS? Jābūt 1 ", 4*Pi*sum(psi**2*r**2*hh)
#endif

! Calculate density
rho=2d0*psi**2d0
! Construct the Hartree potential

call getvh(Ngrid,r,-4d0*Pi*rho*r,vh)
! Construct the exchange-correlation potential
call getvxc(Ngrid,rho,vxc,exc)

! Caulculate energy
 e1=2d0*eig(1)
 e2=-0.5d0*4d0*Pi*hh*sum(vh*rho*r**2d0)
 e3=4d0*Pi*hh*sum((exc-vxc)*rho*r**2d0)
 energy=e1+e2+e3
 Print *,"(11) E=",e1,"+",e2,"+",e3,"=",energy


 vfull=-Z/r+vh+vxc



enddo


! ----- version 2 -------
! Self-consistency loop
!do iscl=1,maxscl

! Diagonalize via the shooting method

! Normalize WFs

! Construct the Hartree potential

! Construct the exchange-correlation potential

!enddo


! ----- version 3 -------
! Self-consistency loop
!do iscl=1,maxscl

! integrate the screened Poisson equation
! $\left( \nabla^2+ 2\epsilon \right \psi(\mathbf{r}) = v(\mathbf{r}) \psi(\mathbf{r})$

! Normalize WFs

! Calculate \epsilon as an expectation value

!enddo


! ----- version 4 -------
! Self-consistency loop
!do iscl=1,maxscl

! Davidson method 

!enddo

deallocate(r)


end program

