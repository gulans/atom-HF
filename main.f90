program atomHF
implicit none
integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:)
real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl = 10
integer :: ir, iscl


allocate(r(Ngrid))

!read input
!...

call gengrid(Ngrid,Rmax,r)

#ifdef debug
do ir=1,Ngrid
  write(*,*) r(ir)
enddo
#endif


! Initial guess for potential or WFs(?)
!...

! ----- version 1 -------
! Self-consistency loop
do iscl=1,maxscl

! Construct Hamiltonian 

! Diagonalize

! Normalize WFs

! Construct the Hartree potential

! Construct the exchange-correlation potential

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

