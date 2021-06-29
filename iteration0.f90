
subroutine iteration0 (Ngrid,r,Z,Nshell,shell_l,lmax,count_l,eig,psi)
       Implicit None
      Integer, Intent (In) :: Ngrid
      real(8), Intent (In) :: r(Ngrid)
      real(8), Intent (In) :: Z
      integer, Intent (In) :: Nshell
      integer, intent (In) :: shell_l(Nshell)
      integer, intent (In) :: lmax
      integer, intent (In) :: count_l(lmax+1)
      real(8), Intent (Out) :: eig(Nshell)
      real(8), Intent (Out) :: psi(Ngrid,Nshell)
      
integer :: ish,l,n,ir,nq
integer :: alpha, k !indexes for Laguerre polynomials 
real(8) :: eig1
real(8) :: lag(Ngrid) !Laguerre polynomials
real(8) :: rho(Ngrid) !Laguerre polynomials argument
real(8) :: fact !factoral function
ish=0

do l=0,lmax
  alpha=2*l+1
  do n=1,count_l(l+1)
    k=n-1
    nq=n+l ! quantum number n
    ish=ish+1
    eig(ish)=-0.5d0*Z**2/nq**2
!    write(*,*)ish,". n=",n," l=",l," eig=",eig(ish)," k:",k," alpha: ",alpha

    rho=2d0*Z*r/dble(nq)


    call Laguerre(Ngrid,k,alpha,rho,lag)
 

    if ((k.eq.4).and.(alpha.eq.3)) then
!      open(11,file='Lag43.dat',status='replace')
!      do ir=1,Ngrid
!         write(11,*)r(ir),lag(ir)
!     enddo
!      close(11)
 
    endif

    psi(:,ish)=dsqrt((2d0*Z/dble(nq))**3*fact(nq-l-1)/(2d0*dble(nq)*fact(nq+l)))*&
            exp(-Z*r/dble(nq))*(2d0*Z*r/dble(nq))**l*lag


  enddo
enddo
end subroutine

subroutine Laguerre(Ngrid,k,alpha,rho,lag)
      Integer, Intent (In) :: Ngrid,k,alpha 
      real(8), Intent (In) :: rho(Ngrid)
      real(8), Intent (Out) :: lag(Ngrid)

      real(8) :: lag0(Ngrid),lag1(Ngrid)
  
      lag0=0d0*rho+1
      lag1=dble(1+alpha)-rho

      if (k.eq.0) then
        lag=lag0
      elseif (k.eq.1) then
        lag=lag1
      else
        do i=2,k      
          lag=((dble(2*(i-1)+1+alpha)-rho)*lag1 - dble(i-1+alpha)*lag0)/dble(i)
          lag0=lag1
          lag1=lag
        enddo
      endif




end subroutine

