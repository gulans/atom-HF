subroutine getvh(Ngrid,r,rho,vh)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),rho(Ngrid) 
real(8), intent(out) :: vh(Ngrid)
integer :: i,j
real(8) :: hh,norm
real(8) :: A(Ngrid,Ngrid),b(Ngrid)
real(8) :: Aij   
real(8) :: Nelec
integer :: info, ipiv(Ngrid)

hh=r(2)-r(1) !RUPJI

!! vh - atrisinājums
!! Matricu vienādojums A0*vh=b

Nelec=sum(rho*r*hh)
 
 b=rho
 b(Ngrid)=b(Ngrid)+Nelec/hh**2d0
 do i=1,Ngrid      
   do j=1,Ngrid
     Aij=0d0
     if (j==i) then  !GALVENAA DIOGNAALE
       Aij=-2d0/(hh**2d0)
     else if (i+1==j) then !augšējā diogn
       Aij=1d0/(hh**2d0)
     else if (i-1==j) then !apakšējā 
       Aij=1d0/(hh**2d0)
     endif
     A(i,j)=Aij
   enddo
 enddo





!!TAGAD ATRISINAM [A][vh]=[b]"
 vh=b
 call dgesv(Ngrid,1,A,Ngrid,IPIV,vh,Ngrid,info)
 vh=vh/r
end subroutine
