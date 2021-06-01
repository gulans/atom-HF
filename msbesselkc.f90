
subroutine msbesselkc (n, z, rez)
       Implicit None
! arguments
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: z
      Complex (8), Intent (Out) :: rez
      Integer :: i
      Complex (8) :: f(0:n)
      if (n.eq.0) then
      rez=exp(-z)*z**(-1)
      else
      f(0)=z**(-1)
      f(1)=(z+dble(1))*z**(-2)


      do i=2,n
        f(i)=f(i-2)+dble(2*i-1)*f(0)*f(i-1)
      enddo

      rez=exp(-z)*f(n)
      endif
end subroutine 
