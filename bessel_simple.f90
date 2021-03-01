
subroutine besi_simp (l, x, rez)
       Implicit None
! arguments
      Integer, Intent (In) :: l
      Real (8), Intent (In) :: x
      Real (8), Intent (Out) :: rez
! local variables
real(8), PARAMETER :: Pi = 3.1415926535897932384d0 

if (l.eq.0) then
rez=sinh(x)/x
elseif (l.eq.1) then
rez=(x*cosh(x)-sinh(x))/x**2
elseif (l.eq.2) then
rez=( (x**2+3d0)*dsinh(x)-3d0*x*dcosh(x) )/x**3
endif
  return
end subroutine

subroutine besk_simp (l, x, rez)
       Implicit None
! arguments
      Integer, Intent (In) :: l
      Real (8), Intent (In) :: x
      Real (8), Intent (Out) :: rez
! local variables
real(8), PARAMETER :: Pi = 3.1415926535897932384d0 

if (l.eq.0) then
rez=exp(-x)/x
elseif (l.eq.1) then
rez=exp(-x)*(x+1)/x**2
elseif (l.eq.2) then
rez=exp(-x)*(x**2+3d0*x+3d0)/x**3
endif
  return
end subroutine
