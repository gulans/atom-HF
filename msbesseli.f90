
subroutine msbesseli (lmax, x, il)
       Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: x
      Real (8), Intent (Out) :: il (0:lmax)
! local variables
! staring value for l above lmax (suitable for lmax < 50)
      Integer :: l
      Real (8) :: xi, i0, i1, it, t1, t2, xm
      Real (8) :: cs
      Real (8) :: f
      Real (8) :: f0
      Real (8) :: f1
      Integer :: m
      Integer :: k
      Integer :: msta1
      Integer :: msta2
      If ((lmax .Lt. 0) .Or. (lmax .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(msbesseli): lmax out of range : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      If ((x .Lt. 0.d0) .Or. (x .Gt. 1.d5)) Then
         Write (*,*)
         Write (*, '("Error(msbesseli): x out of range : ", G18.10)') x
         Write (*,*)
         Stop
      End If
      xi = 1.d0 / x
      xm = 1.d-8 
     
      If (x .Lt. xm) Then
         Do l = 0, lmax
            il (l) = 0.0d0
         End do
            il (0) = 1.d0
         Return
       End If
         il (0) = Sinh (x) * xi
         If (lmax .Eq. 0) Return
         il (1) = xi * (Cosh (x) - (Sinh(x)*xi))
         If (lmax .Eq. 1) Return
         i0 = il (0)
      If (lmax .Ge. 2) Then
         m =  msta1(x, 200) 
         If (m .Lt. lmax) Then
             m = lmax
         Else 
             m = msta2(x, lmax, 15)        
         End If
         f0 = 0.d0
         f1 = 1.d0 - 100
         Do l = m, 0, -1
            f = (2.d0 * l + 3.d0) * f1/x + f0
            If (l .Le. lmax) Then
            il(l) = f
            End If 
            f0 = f1 
            f1 = f 
         End Do
      cs = i0/f
      Do l = 0, lmax
         il(l) = cs * il(l)
      End Do
      End If
  Return
End Subroutine

function msta1 ( x, mp )
  implicit none
  real (8) :: a0
  real (8) :: envj
  real (8) :: f
  real (8) :: f0
  real (8) :: f1
  integer (4) :: it
  integer (4) :: mp
  integer (4) :: msta1
  integer (4) :: n0
  integer (4) :: n1
  integer (4) :: nn
  real (8) :: x
  a0 = abs ( x )
  n0 = int ( 1.1D+00 * a0 ) + 1
  f0 = envj ( n0, a0 ) - mp
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - mp
  do it = 1, 20       
    nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )                  
    f = envj ( nn, a0 ) - mp
    if ( abs ( nn - n1 ) .Lt. 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do
  msta1 = nn
  return
end function

function msta2 ( x, n, mp )
  implicit none
  real (8) :: a0
  real (8) :: ejn
  real (8) :: envj
  real (8) :: f
  real (8) :: f0
  real (8) :: f1
  real (8) :: hmp
  integer (4) :: it
  integer (4) :: mp
  integer (4) :: msta2
  integer (4) :: n
  integer (4) :: n0
  integer (4) :: n1
  integer (4) :: nn
  real (8) :: obj
  real (8) :: x
  a0 = abs ( x )
  hmp = 0.5D+00 * mp
  ejn = envj ( n, a0 )
  if ( ejn .Le. hmp ) then
    obj = mp
!
!  Original code:
!
!   n0 = int ( 1.1D+00 * a0 )
!
!  Updated code:
!
    n0 = int ( 1.1D+00 * a0 ) + 1
  else
    obj = hmp + ejn
    n0 = n
  end if
  f0 = envj ( n0, a0 ) - obj
  n1 = n0 + 5
  f1 = envj ( n1, a0 ) - obj
  do it = 1, 20
    nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
    f = envj ( nn, a0 ) - obj
    if ( abs ( nn - n1 ) .Lt. 1 ) then
      exit
    end if
    n0 = n1
    f0 = f1
    n1 = nn
    f1 = f
  end do
  msta2 = nn + 10
  return
end function

function envj ( n, x )
!
  implicit none
  real (8) :: envj
  integer (4) :: n
  real (8) :: x
!
    envj = 0.5D+00 * log10 ( 6.28D+00 * n ) &
      - n * log10 ( 1.36D+00 * x / n )
!
  return
end function
