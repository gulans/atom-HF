subroutine getvh(Ngrid,r,rho,vh)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),rho(Ngrid) 
real(8), intent(out) :: vh(Ngrid)
integer :: i,j
real(8) :: hh
real(8) :: b(Ngrid), AB(4,Ngrid)
real(8) :: Aij   
real(8) :: Nelec
integer :: info, ipiv(Ngrid)

hh=r(2)-r(1) !RUPJI

!! vh - atrisinājums
!! Matricu vienādojums A0*vh=b

b=-4d0*Pi*rho*r
Nelec=sum(-4d0*Pi*rho*r**2*hh)
write(*,*)'Nelec: ', Nelec 
 
 b(Ngrid)=b(Ngrid)+Nelec/hh**2d0

 do j=1,Ngrid
   do i=1,4
     Aij=0d0
     if (i+j-3==j) then !main diagonal
        Aij=-2d0/(hh**2d0)
     else if ((i+j-3==j+1).OR.(i+j-3==j-1)) then !sub or super diagonal 
        Aij=1d0/(hh**2d0)
     end if
     AB(i,j)=Aij
   enddo
 enddo


vh=b
call dgbsv(Ngrid,1,1,1,AB,4 ,ipiv,vh,Ngrid,info) 
!(dimension,KL=subdiognals,KU=superdiognals,b width,Matrix(LDAB,N),LDAB=2*KL+KU+1, int IPIV(N), b,LDB, info)

 vh=vh/r

end subroutine
