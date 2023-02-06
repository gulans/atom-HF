subroutine gen_vn(Ngrid,Z,r,sigma,vn)
use modinteg
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: Z,r(Ngrid),sigma
real(8), intent(out) :: vn(Ngrid)
integer :: ir

real(8), PARAMETER :: Pi = 3.1415926535897932384d0


real(8) ::test,ksi,rho0,norm, q(Ngrid),vn_sk(Ngrid),ftemp1(Ngrid),ftemp2(Ngrid)


  if (sigma.lt.1d-30)then
    vn=-Z/r 
  else
    vn=(-Z/r)*erf(r/(dsqrt(2d0)*sigma)) 
  endif
write(*,*)"nr",Ngrid
write(*,*)"r",r(1),r(2),r(Ngrid)
write(*,*)vn(1:10)  


ksi=1d0/(2d0*sigma**2)

vn=(-Z/r)*erf(r*dsqrt(ksi))







rho0=4d0*Pi*(ksi/Pi)**(3d0/2d0)
q=-Z*rho0*exp(-ksi*r**2)

q=-Z*4d0*Pi*(ksi/Pi)**(3d0/2d0)*exp(-ksi*r**2)


call integ_v(Ngrid,r,q*r**2,test)
write(*,*)"test total charge",test,"Z",Z






!sigma=1d-3
!q=exp(-0.5d0*r**2/sigma**2)
!call integ_v(Ngrid,r,q*r**2,norm)
!q=q*Z/norm



call  integ_f(Ngrid,r,r**2*q,ftemp1)
call  integ_f_rev(Ngrid,r,r*q,ftemp2)

vn_sk=ftemp1/r+ftemp2









open(11,file='kodols.dat',status='replace')
!write(11,*)"r q vn,vn_sk"
  do ir=1,Ngrid
     write(11,*)r(ir), q(ir),vn(ir),vn_sk(ir)
  enddo
close(11)




!stop
!do ir=1, 60
!  write(*,*)r(ir),vn(ir), -Z/r(ir)
!enddo

end subroutine
