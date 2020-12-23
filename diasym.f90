subroutine diasym(Ngrid,H,eig)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(out) :: H(Ngrid,Ngrid)
real(8), intent(out) :: eig(Ngrid)
integer :: i,j,ir,l,inf
real(8) :: work(Ngrid**2),Aij, AB(2,Ngrid)
real(8) :: eigvec(Ngrid,Ngrid),work1(3*Ngrid-2)


 AB(2,1)=H(1,1)
  
  do j=2,Ngrid
   do i=1,2
     AB(i,j)=H(i+j-2,j)
   enddo
  enddo


 call dsbev('V','U',Ngrid,1,AB,2,eig,eigvec,Ngrid,work1,inf)
 H=eigvec


end subroutine
