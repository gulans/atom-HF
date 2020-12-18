program atomHF
implicit none
integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:)
real(8), parameter :: Rmax = 10d0

integer :: ir

allocate(r(Ngrid))


call gengrid(Ngrid,Rmax,r)

do ir=1,Ngrid
  write(*,*) r(ir)
enddo

deallocate(r)


end program

