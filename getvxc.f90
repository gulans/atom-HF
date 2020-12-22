subroutine getvxc(Ngrid,rho,vxc,exc)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384
integer, intent(in) :: Ngrid
real(8), intent(in) :: rho(Ngrid) 
real(8), intent(out) :: vxc(Ngrid),exc(Ngrid)
integer :: i
real(8) :: rhoi, vxc_fun,exc_fun

do i=1, Ngrid
  rhoi=rho(i)
  if (rhoi.lt.(1e-20)) then
    vxc(i)=0
    exc(i)=0
  else
    vxc(i)=vxc_fun(rhoi)
    exc(i)=exc_fun(rhoi) 
  endif   
enddo
end subroutine

real(8) function vxc_fun(n)  !Funkcijas nosaukums jālieto par mainīgo???? Tiešām?
  IMPLICIT NONE
  real(8), intent(in)    :: n ! input
  real(8)  :: A,a1,b1,b2,b3,b4
  real(8), PARAMETER :: Pi = 3.1415927
  A=0.031091
  a1=0.082477
  b1=5.1486
  b2=1.6483
  b3=0.23647
  b4=0.20614
  vxc_fun=        (-3*(3/Pi)**0.3333333333333333)/(4.*(1/n)**0.3333333333333333) - &
       2*A*(1 + (a1*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666)* &
        Log(1 + 1/(2.*A*((b1*(1/n)**0.16666666666666666*(3/Pi)**0.16666666666666666)/2**0.3333333333333333 +  &
               (b2*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666 +  &
               (b3*Sqrt(1/n)*Sqrt(3/Pi))/2. + &
               (3**0.6666666666666666*b4*(1/n)**0.6666666666666666)/&
                (2.*2**0.3333333333333333*Pi**0.6666666666666666)))) + &
       n*(((1 + (a1*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666)*&
             (-(b3*(1/n)**1.5*Sqrt(3/Pi))/4. - &
               (b4*(1/n)**1.6666666666666667)/(6**0.3333333333333333*Pi**0.6666666666666666) - &
               (b2*(1/n)**1.3333333333333333)/(6**0.6666666666666666*Pi**0.3333333333333333) - &
               (b1*(1/n)**1.1666666666666667)/&
                (2.*2**0.3333333333333333*3**0.8333333333333334*Pi**0.16666666666666666)))/&
           ((1 + 1/(2.*A*((b1*(1/n)**0.16666666666666666*(3/Pi)**0.16666666666666666)/2**0.3333333333333333 + &
                    (b2*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666 + &
                    (b3*Sqrt(1/n)*Sqrt(3/Pi))/2. + &
                    (3**0.6666666666666666*b4*(1/n)**0.6666666666666666)/&
                     (2.*2**0.3333333333333333*Pi**0.6666666666666666))))*&
             ((b1*(1/n)**0.16666666666666666*(3/Pi)**0.16666666666666666)/2**0.3333333333333333 + &
                (b2*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666 + &
                (b3*Sqrt(1/n)*Sqrt(3/Pi))/2. + &
                (3**0.6666666666666666*b4*(1/n)**0.6666666666666666)/&
                 (2.*2**0.3333333333333333*Pi**0.6666666666666666))**2) - &
          ((1/n)**0.6666666666666666*(3/Pi)**0.3333333333333333)/4. + &
          (A*a1*(1/n)**1.3333333333333333*(2/Pi)**0.3333333333333333*&
             Log(1 + 1/&
                (2.*A*((b1*(1/n)**0.16666666666666666*(3/Pi)**0.16666666666666666)/2**0.3333333333333333 + &
                    (b2*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666 + &
                    (b3*Sqrt(1/n)*Sqrt(3/Pi))/2. + &
                    (3**0.6666666666666666*b4*(1/n)**0.6666666666666666)/&
                     (2.*2**0.3333333333333333*Pi**0.6666666666666666)))))/3**0.6666666666666666)

  return
end function

real(8) function exc_fun(n)  !Funkcijas nosaukums jālieto par mainīgo???? Tiešām?
  IMPLICIT NONE
  real(8), intent(in)    :: n ! input
  real(8)  :: A,a1,b1,b2,b3,b4
  real(8), PARAMETER :: Pi = 3.1415927
  A=0.031091
  a1=0.082477
  b1=5.1486
  b2=1.6483
  b3=0.23647
  b4=0.20614

  exc_fun=(-3*(3/Pi)**0.3333333333333333)/(4.*(1/n)**0.3333333333333333) - &
       2*A*(1 + (a1*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666)*&
        Log(1 + 1/(2.*A*((b1*(1/n)**0.16666666666666666*(3/Pi)**0.16666666666666666)/2**0.3333333333333333 + &
               (b2*(1/n)**0.3333333333333333*(3/Pi)**0.3333333333333333)/2**0.6666666666666666 + &
               (b3*Sqrt(1/n)*Sqrt(3/Pi))/2. + &
               (3**0.6666666666666666*b4*(1/n)**0.6666666666666666)/&
                (2.*2**0.3333333333333333*Pi**0.6666666666666666))))
  return
end function
