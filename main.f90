program atomHF

implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0, alpha=0.50d0
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:,:),rho(:),vfull1(:)
integer, allocatable :: shell_n(:),shell_l(:),count_l(:)
real(8), allocatable :: shell_occ(:),psi_eig(:)

!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl = 20
integer :: ir,il,icl,il_icl,iscl,lmax

real(8) :: Z
real(8) :: norm,Rmin,Rmax,hh,e1,e2,e3,energy,energy0
integer :: grid,Ngrid, Nshell, ish
integer :: i,j,countl0, version
   
logical :: file_exists
!read input
read(*,*) 
read(*,*) Z, Rmin, Rmax, Ngrid, version
read(*,*)
read(*,*) grid
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



allocate(r(Ngrid),vfull(Ngrid),vh(Ngrid),vxc(Ngrid),exc(Ngrid),eig(Ngrid),rho(Ngrid),vfull1(Ngrid))
allocate(psi(Ngrid,Nshell),psi_eig(Nshell))

call gengrid(grid,Ngrid,Rmin,Rmax,r)
hh=r(2)-r(1)


! Initial guess for potential
!...
vfull=0*r !vfull is vh+vxc, core potential -Z/r is used seperatly - version 1 has to be modified


if (version.eq.1) then
allocate(H(Ngrid,Ngrid))
! ----- version 1 -------
! Self-consistency loop
write(*,*)"************* version 1 **************"

do iscl=1,1!maxscl
#ifdef debug

  write(*,*) "iteration: ",iscl
#endif

il_icl=0
do il=1,lmax+1
  ! Construct Hamiltonian 
  
  call genham(Ngrid,r,vfull+0.5d0*dble(il-1)*dble(il)*r**(-2d0),H)

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
    write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," normalised to:",4d0*Pi*sum(psi(:,il_icl)**2d0*r**2d0*hh),&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)
#endif
  enddo
enddo
!  write 3 wave functions to file 
!  open(11,file='wf_v1.dat',status='replace')
!  write(11,*)"r psi1 psi2"
!   do i = 1,Ngrid 
!     write(11,*)r(i),psi(i,1),psi(i,2),psi(i,3)
!   end do
!   close(11)


! Calculate density
  rho=0d0*rho

  
  
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2d0
  enddo

  ! Construct the Hartree potential

!  call getvh(Ngrid,r,rho,vh)
  call getvhtrapez(Ngrid,r,rho,vh)
! Construct the exchange-correlation potential
  call getvxc(Ngrid,rho,vxc,exc)

! Caulculate energy
  e1=0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
  e2=0
  e3=0
  do ir=1, Ngrid-1
    hh=r(ir+1)-r(ir)
    e2=e2-0.5d0*4d0*Pi*hh*vh(ir)*rho(ir)*r(ir)**2d0
    e3=e3+4d0*Pi*hh*(exc(ir)-vxc(ir))*rho(ir)*r(ir)**2d0
  enddo
  energy0=energy
  energy=e1+e2+e3

  Print *,"(11) E=",e1,"+",e2,"+",e3,"=",energy
 
  vfull1=vh+vxc
  vfull=alpha*vfull1+(1d0-alpha)*vfull
enddo

deallocate(H)



else if (version.eq.2) then
write(*,*)"********** version 2 *************"
! ----- version 2 -------
! Self-consistency loop
do iscl=1,maxscl
  il_icl=0
  do il=1,lmax+1
    ! Diagonalize via the shooting method
    call diashoot(Ngrid,r,Z,vfull,il-1,count_l(il),il_icl,Nshell,eig,psi)

    ! Normalize WFs
    do icl=1,count_l(il)
      il_icl=il_icl+1
!      psi(:,il_icl)=H(:,icl)
      psi_eig(il_icl)=eig(icl)
      call wfnorm(Ngrid,r,psi(:,il_icl))
#ifdef debug
  norm=0
  do ir=1,Ngrid-1
    hh=r(ir+1)-r(ir)
    norm=norm+4d0*Pi*psi(ir,il_icl)**2d0*r(ir)**2d0*hh
  enddo
  write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," normalised to:",norm,&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)
write(*,*)"PSI: ",psi(1,il_icl),psi(2,il_icl),psi(3,il_icl),psi(4,il_icl)

#endif
    enddo
  enddo

  !write 3 wave functions to file 
!  open(11,file='wf_v2_Newton_new.dat',status='replace')
!  write(11,*)"r psi1 psi2"
!   do i = 1,Ngrid 
!     write(11,*)r(i),psi(i,1),psi(i,2),psi(i,3)
!   end do
!   close(11)

 
 
  ! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2d0
  enddo

  ! Construct the Hartree potential
  call getvhtrapez(Ngrid,r,rho,vh)
  ! Construct the exchange-correlation potential
  call getvxc(Ngrid,rho,vxc,exc)

  ! Caulculate energy
  e1=0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
 
  e2=0
  e3=0
  do ir=1, Ngrid-1
    hh=r(ir+1)-r(ir)
    e2=e2-0.5d0*4d0*Pi*hh*vh(ir)*rho(ir)*r(ir)**2d0
    e3=e3+4d0*Pi*hh*(exc(ir)-vxc(ir))*rho(ir)*r(ir)**2d0
  enddo

 
  energy0=energy
  energy=e1+e2+e3
  Print *,iscl,".iter E=",e1,"+",e2,"+",e3,"=",energy

  vfull1=vh+vxc
  vfull=alpha*vfull1+(1-alpha)*vfull

enddo

endif
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

!  write results 

  inquire(file='results.dat',EXIST=file_exists)
  if (file_exists) then
     open(11,file='results.dat',status='old', access='append')
  else
     open(11,file='results.dat',status='new')
  endif
  write(11,*)"Z version iterations Ngrid Rmax Grid"
  write(11,*)Z, version, iscl-1, Ngrid, Rmax, grid
  write(11,*)"n l eigval"
  il_icl=0
  do il=1,lmax+1
    do icl=1, count_l(il)
      il_icl=il_icl+1
      write(11,*) icl, il, psi_eig(il_icl) 
    enddo
  enddo
  write(11,*)"Tot Energy: ", energy
  write(11,*)"dE        : ", energy-energy0
  write(11,*)"" 
  close(11)

deallocate(r,vfull,vh,vxc,exc,eig,psi,rho,shell_n,shell_l,count_l,shell_occ)


end program

