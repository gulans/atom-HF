program atomHF

implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0, alpha=0.5d0
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:,:),rho(:),vfull1(:)
integer, allocatable :: shell_n(:),shell_l(:),count_l(:)
real(8), allocatable :: shell_occ(:),psi_eig(:)

!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl = 20
integer :: ir,il,icl,il_icl,iscl,lmax

real(8) :: Z
real(8) :: Rmax,hh,e1,e2,e3,energy
integer :: Ngrid, Nshell, ish
integer :: i,j,countl0

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

lmax=maxval(shell_l)
allocate(count_l(lmax+1))


!counts how many particular l shells (there must be a shorter way to do it)
do i=0,lmax
 countl0=0
 do j=1,Nshell
 if (i.eq.shell_l(j)) then 
         countl0=countl0+1
 endif
 enddo
 write(*,*)"l=",i," countl=", countl0
 count_l(i+1)=countl0
enddo
 write(*,*) "count_l=",count_l(:)



allocate(r(Ngrid),vfull(Ngrid),vh(Ngrid),vxc(Ngrid),exc(Ngrid),H(Ngrid,Ngrid),eig(Ngrid),rho(Ngrid),vfull1(Ngrid))
allocate(psi(Ngrid,Nshell),psi_eig(Nshell))

call gengrid(Ngrid,Rmax,r)
hh=r(2)-r(1)


! Initial guess for potential
!...
vfull=-Z/r



! ----- version 1 -------
! Self-consistency loop

do iscl=1,maxscl
#ifdef debug
  write(*,*) "iteration: ",iscl
#endif

il_icl=0
do il=1,lmax+1
  ! Construct Hamiltonian 
  
  call genham(Ngrid,r,vfull+0.5d0*(il-1)*il*r**(-2d0),H)

  ! Diagonalize
  call diasym(Ngrid,H,eig)

  ! Normalize WFs
  do icl=1,count_l(il)
    il_icl=il_icl+1
    psi(:,il_icl)=H(:,icl)  !tas ir Psi*r
    psi(:,il_icl)=psi(:,il_icl)/r
    psi_eig(il_icl)=eig(icl)
    call wfnorm(Ngrid,r,psi(:,il_icl))
#ifdef debug
    write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," normalised to:",4*Pi*sum(psi(:,il_icl)**2*r**2*hh),&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)
#endif
  enddo
enddo
  !write 3 wave functions to file 
!  open(11,file='wf.dat',status='replace')
!  write(11,*)"r psi1 psi2"
!   do i = 1,Ngrid 
!     write(11,*)r(i),psi(i,1)!,psi(i,2)!,psi(i,3)
!   end do
!   close(11)


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
  e1=0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
  e2=-0.5d0*4d0*Pi*hh*sum(vh*rho*r**2d0)
  e3=4d0*Pi*hh*sum((exc-vxc)*rho*r**2d0)
  energy=e1+e2+e3
  Print *,"(11) E=",e1,"+",e2,"+",e3,"=",energy
 
  vfull1=-Z/r+vh+vxc
  vfull=alpha*vfull1+(1-alpha)*vfull
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

deallocate(r,vfull,vh,vxc,exc,H,eig,psi,rho,shell_n,shell_l,count_l,shell_occ)


end program

