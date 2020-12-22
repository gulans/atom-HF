program atomHF
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:,:),rho(:)
integer, allocatable :: shell_n(:),shell_l(:)
real(8), allocatable :: shell_occ(:)

!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl = 20
integer :: ir, iscl

real(8) :: Z
real(8) :: Rmax,hh,e1,e2,e3,energy
integer :: Ngrid, Nshell, ish
integer :: i


!read input
read(*,*) 
read(*,*) Z, Rmax, Ngrid
read(*,*)
read(*,*) Nshell
allocate(shell_n(Nshell),shell_l(Nshell),shell_occ(Nshell))
read(*,*)
do ish=1,Nshell
  read(*,*) shell_n(ish), shell_l(ish), shell_occ(ish)
enddo
#ifdef debug
  write(*,*) "shell_n ", shell_n(:)
  write(*,*) "shell_l ", shell_l(:)
  write(*,*) "shell_occ ", shell_occ(:)
#endif


allocate(r(Ngrid),vfull(Ngrid),vh(Ngrid),vxc(Ngrid),exc(Ngrid),H(Ngrid,Ngrid),eig(Ngrid),rho(Ngrid))
allocate(psi(Ngrid,Nshell))

call gengrid(Ngrid,Rmax,r)
hh=r(2)-r(1)


! Initial guess for potential
!...
vfull=-Z/r



! ----- version 1 -------
! Self-consistency loop

do iscl=1,maxscl

  ! Construct Hamiltonian 

  call genham(Ngrid,r,vfull,H)

  ! Diagonalize
  call diasym(Ngrid,H,eig)

#ifdef debug
  write(*,*) "iteration",iscl, ". Eigenvalues:"
     do ir=1,Nshell
       write(*,*) eig(ir)
     enddo
#endif

  ! Normalize WFs
  do ish=1,Nshell
    psi(:,ish)=H(:,ish)  !tas ir Psi*r
    psi(:,ish)=psi(:,ish)/r
    call wfnorm(Ngrid,r,psi(:,ish))
#ifdef debug
    Print *,"PSI",ish," nomalised to:", 4*Pi*sum(psi(:,ish)**2*r**2*hh)
#endif
  enddo
  
  !write 2 wave functions to file 
  open(11,file='wf.dat',status='replace')
  write(11,*)"r psi1 psi2"
   do i = 1,Ngrid 
     write(11,*)r(i),psi(i,1),psi(i,2)
   end do
   close(11)


! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
    rho=rho+shell_occ(ish)*psi(:,ish)**2d0
  enddo
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

