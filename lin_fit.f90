subroutine lin_fit(Ngrid,Nbase,r,psi,base,fitc)
implicit none

integer, intent(in) :: Ngrid,Nbase
real(8), intent(in) :: r(Ngrid),psi(Ngrid),base(Ngrid,Nbase)
real(8), intent(out) :: fitc(Nbase)
integer :: i,ir
real(8) :: A(Ngrid,Nbase), A1(Nbase,Nbase)
real(8) :: sigma(Ngrid),y(Ngrid)

sigma=r**2
!A(:,1)=r**2
!A(:,2)=r**3
!y=4d0*r**2+99d0*r**3

!sigma=1d0
do i=1, Nbase
  A(:,i)=base(:,i)/sigma
enddo
y=psi/sigma

 open(11,file="base_and_fun.dat",status='replace')
 do ir=1, Ngrid
  write(11,*)r(ir),A(ir,1),A(ir,2),y(ir)
 enddo
 close(11)




call inver( matmul(transpose(A),A) ,A1, Nbase )

fitc=matmul( matmul(A1, transpose(A)) ,y ) 
!write(*,*) "fitc",fitc
end subroutine
