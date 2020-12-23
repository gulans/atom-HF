subroutine diasym(Ngrid,H,eig)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(out) :: H(Ngrid,Ngrid)
real(8), intent(out) :: eig(Ngrid)



real(8) :: E(Ngrid), D(Ngrid)
real(8) :: VL,VU
integer, parameter :: IL=1 !the index of the smallest eigenvalue to be returned
integer, parameter :: IU=20 !the index of the largest eigenvalue to be returned
integer, parameter :: M=IU-IL+1
integer :: eigvecfound
logical :: TRYRAC
 real(8) :: eigvec(Ngrid,M),WORK(18*Ngrid)
 integer :: info,i,IWORK(10*Ngrid),ISUPPZ(2*M)



 do i=1,Ngrid
     D(i)=H(i,i) !diagonal
 enddo
  
 do i=1,Ngrid-1
     E(i)=H(i+1,i) !subdiagonal elements
 enddo
 TRYRAC=.TRUE. ! after dstemr its value is FALSE
 call dstemr('V','I',Ngrid,D,E,VL,VU,IL,IU,eigvecfound,eig,eigvec,Ngrid,M,ISUPPZ,TRYRAC,WORK,18*Ngrid,IWORK,10*Ngrid,INFO)
 if (info.ne.0) then
         write(*,*) "dyasym.f90 diagonalisation problem posibility (info not 0)"
 end if

 do i=1,M
     H(:,i)=eigvec(:,i)
 enddo

end subroutine
