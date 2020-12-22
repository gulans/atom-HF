subroutine diasym(Ngrid,H,eig)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: H(Ngrid) 
real(8), intent(out) :: eig(Ngrid)
integer :: ir,l,inf
real(8) :: work(Ngrid**2)

l=Ngrid**2

call dsyev('V','U',Ngrid,H,Ngrid,eig,work,l,inf)
 !print *,'dsyev info: ', inf
 !print *,'dsyev work(1): ', work(1)
 !print *,'dsyev lwork: ', l
end subroutine
