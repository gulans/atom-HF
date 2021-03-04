!
!BOP
! !ROUTINE: msbesselk
! !INTERFACE:
!
!
Subroutine msbesselk (lmax, x, kl)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   kl   : array of returned values (out,real(0:lmax))
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: lmax
      Real (8), Intent (In) :: x
      Real (8), Intent (Out) :: kl (0:lmax)
! local variables
! staring value for l above lmax (suitable for lmax < 50)
      Integer :: l
      Real (8) :: xi, k0, k1, kt, t1, t2
      If ((lmax .Lt. 0) .Or. (lmax .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(msbesselk): lmax out of range : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      If ((x .Lt. 0.d0) .Or. (x .Gt. 1.d5)) Then
         Write (*,*)
         Write (*, '("Error(msbesselk): x out of range : ", G18.10)') x
         Write (*,*)
         Stop
      End If
! treat x << 1
      xi = 1.d0 / x
      If (x .Lt. 1.d-8) Then
         kl (0) = -xi   
         t1 = 1.d0
         t2 = xi 
         Do l = 1, lmax
            t1 = t1 * dble (2*l-1)
            t2 = t2 * xi
            kl (l) = t2 * t1
         End Do
         Return
      End If
! recursion relation
         kl (0) = Exp (-x) * xi
         If (lmax .Eq. 0) Return
         kl (1) = kl(0) * (1.d0 + xi)
         If (lmax .Eq. 1) Return
         k0 = kl (0)
         k1 = kl (1)
         Do l = 2, lmax
!            kt = k0 + (dble (2*l+1) * xi * k1)
           kt = k0 + (dble (2*l-1) * xi * k1)
            k0 = k1
            k1 = kt
            kl (l) = k1
         End Do
         Return
End Subroutine
!EOC

