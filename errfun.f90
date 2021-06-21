subroutine errfun(Ngrid,r,n,mu,rsfunC)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid,n
real(8), intent(in) :: r(Ngrid),mu
complex(8),intent(out) :: rsfunC(n,2) 

real(8) ::AA1(n),AA2(n),BB1(n),BB2(n)
complex(8) :: rez(Ngrid)
real(8) :: rmax, temp(Ngrid)
real(8) :: a1,a2
integer :: i,ir,k


write(*,*)n,mu

if (n.eq.8)then

AA1=(/1.2860987030345912d01*0.5d0,-4.9051058740247786d00,-1.9357963384430639d00,1.0117925923974511d00,&
        -9.9013582507628040d-02,-3.3813615287899379d-03,1.0330516726180475d-03,-2.2002738770173973d-05/)

AA2=(/0.0000000000000000d00,8.4204906883706165d00,-3.7444370420962709d00,1.5201429359841603d-01,&
        9.5997444318657268d-02,-1.1743959149248155d-02,1.5238534910174336d-04,5.8451000772024675d-06/)

BB1=(/4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00,&
        4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00,4.5617402703990644d00/)

BB2=(/0.0000000000000000d00,1.0160216352396843d00,2.0413334512076715d00,3.0805139945811777d00,&
        4.1164481828964572d00,5.0646792984847675d00,5.9785623723578114d00,7.2027553450037622d00/)
rmax=5.3d0

endif



rmax=rmax/mu
   do i = 1,n
      rsfunC(i,1)=cmplx(AA1(i),AA2(i),8)
      rsfunC(i,2)=cmplx(BB1(i),BB2(i),8)
   end do

!write(*,*)" A"
!do i=1, n
!write(*,*)rsfunC(i,1)
!enddo
!write(*,*)" B"
!do i=1, n
!write(*,*)rsfunC(i,2)
!enddo


 rsfunC(:,2)=rsfunC(:,2)*mu




   rez=cmplx(0d0,0d0,8)*r
 !  rez=rez+rsfunC(1,1)*exp(-rsfunC(1,2)*r)
   do i=1, n
     rez=rez+rsfunC(i,1)*exp(-rsfunC(i,2)*r)
     rez=rez+conjg(rsfunC(i,1))*exp(-conjg(rsfunC(i,2))*r)
   enddo
   open(11,file='erfc_test.dat',status='replace')
   write(11,*)"erfc RE(rez) Im(rez) erfc-RE(rez) mu=",mu," rmax=",rmax
   do i=1, Ngrid
     write(11,*)r(i), erfc(mu*r(i)), realpart(rez(i)), imagpart(rez(i)),erfc(mu*r(i))-realpart(rez(i))
   enddo
   close(11)

end subroutine

