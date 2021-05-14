function PHInk(n,k,k1,k2)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: n,k
real(8), intent(in) :: k1,k2

real(8) :: PHInk,rez

real(8) :: summa
real(8) :: Dnk,Hn,Fn,facti,dfact

integer :: m,j 


if ((k2.lt.0.4d0).or.(k1.lt.0.5d0).or.(k2.lt.2d0*k1)) then
  summa=0d0         !Taylor expantion
  do j=0,k
  summa=summa+Dnk(n,j,k1)*k2**(n+2*j)/(k1**(n+1))
  enddo
  rez=summa
else
  summa=0d0
  do m=1,n
    summa=summa+Fn(n-m,k1,k2)*Hn(n,k1,k2)*(k1**(2*m)+k2**(2*m))/((k1*k2)**m)
  enddo
  rez=Fn(n,k1,k2)+summa
endif

PHInk=rez
end function 

real(8) function Dnk(n,k,k1)
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: n,k
real(8), intent(in) :: k1
integer :: b1,b2
real(8) :: dfact,fact,summa,b,binomial
if (k.eq.0) then
  summa=0d0
  do m=1, n 
    summa=summa+(2d0*k1**2)**(-m)/dfact(2*n-2*m+1)
  enddo
 !Lehtola 2020:
 ! Dnk=erfc(k1)+dexp(-k1**2)*(2d0*k1**2)**(n+1)*summa/dsqrt(Pi)
 !Angyan 2006:
 Dnk=erfc(k1)+dexp(-k1**2)*2d0**(n+1)*k1**(2*n+1)*summa/dsqrt(Pi)

else
summa=0d0
  do m=1,k
    b1=m-k-1
    b2=m-1
    b=binomial(b1,b2)
    summa=summa+dble(b)*(2d0*k1**2)**(k-m)/dfact(2*n+2*k-2*m+1)
  enddo 

 !Lehtola 2020:
 ! Dnk=dexp(-k1**2)*(2d0*k1**2)**(n+1)*dble(2*n+1)/(dsqrt(Pi)*fact(k)*dble(2*n+2*k+1))*summa
 !Angyan 2006:
  Dnk=dexp(-k1**2)*2d0**(n+1)*k1**(2*n+1)*dble(2*n+1)/(dsqrt(Pi)*fact(k)*dble(2*n+2*k+1))*summa

endif
return
end function


real(8) function Fn(n,k1,k2)
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: n
real(8), intent(in) :: k1,k2
integer :: p
real(8) :: fact,summa
summa=0d0

do p=0, n
 summa=summa + (-1d0/(4d0*k1*k2))**(p+1) * (fact(n+p)/(fact(p)*fact(n-p))) *&
         ((-1)**(n-p)*dexp(-(k2+k2)**2) - dexp(-(k2-k1)**2) ) 
enddo
Fn=(2d0/dsqrt(Pi))*summa
return
end function 

real(8) function Hn(n,k1,k2)
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: n
real(8), intent(in) :: k1,k2
real(8) :: rez,rez1,rez2

rez1=(k1**(2*n+1)+k2**(2*n+1))*erfc(k1+k2)
rez2=(k1**(2*n+1)-k2**(2*n+1))*erfc(k1-k2)
Hn=(rez1-rez2)/(2d0*(k1*k2)**(n+1))
return
end function


real(8) function fact(n)
integer, intent(in) :: n

  if (n < 0) then
    write(*,*) 'ercf_kernel.f90 factorial for negative integer?',n
    stop
  endif
  fact = 1.0d0
  do i = 2, n
    fact = fact * dble(i)
  enddo
return
end function    

real(8) function dfact(n)
integer, intent(in) :: n

  if (n < 0) then
    write(*,*) 'ercf_kernel.f90 dfactorial for negative integer?',n
    stop
  endif
  dfact = 1.0d0
  do i = n, 1,-2
    dfact = dfact * dble(i)
  enddo
return
end function

real(8) function binomial(b1,b2)
integer, intent(in) :: b1,b2
integer :: b
real(8) ::fact
if (b2.ge.0) then
  if (b1.ge.0) then
    binomial=fact(b1)/(fact(b2)*fact(b1-b2))
  else
    b=b2-b1-1
    binomial=(-1d0)**b2*fact(b)/(fact(b2)*fact(b-b2))
  endif
else
   write(*,*) 'ercf_kernel.f90 binomail b2<0, b1=',b1,'b2=',b2
endif
!  write(*,*) 'b1=',b1,'b2=',b2,"rez",binomial

return
end function

