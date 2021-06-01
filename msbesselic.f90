
subroutine msbesselic (n, z, rez)
       Implicit None
! arguments
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: z
      Complex (8), Intent (Out) :: rez
      Integer :: i
      Complex (8) :: g(-n-1:n)
      if (n.eq.0) then
      rez=sinh(z)*z**(-1)
      else
      g(0)=z**(-1)
      g(1)=-z**(-2)


      do i=2,n
        g(i)=g(i-2)-dble(2*i-1)*g(0)*g(i-1)
      enddo
      do i=-1,-n-1,-1
        g(i)=g(i+2)+dble(2*i+3)*g(0)*g(i+1)
      enddo

      rez=g(n)*sinh(z)+g(-n-1)*cosh(z)
      endif
end subroutine 
