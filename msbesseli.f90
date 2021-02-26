
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
