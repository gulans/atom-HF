subroutine msbesselkc (n, z, rez)
       Implicit None
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: z
      Complex (8), Intent (Out) :: rez
      Integer :: i
      Complex (8) :: f(0:n)

      Complex (8) :: rezi,rezk

Complex(16) :: f16(0:n)
complex(16)  :: rez16
complex(16)  :: y



   if ((n.eq.1) .and. (abs(realpart(z)).lt.1d-3).and.(abs(imagpart(z)).lt.1d-3))  then
     rez=1d0/z**2 - (1d0/2d0)+ (1d0/3d0)*z - (1d0/8d0)*z**2 +  (1d0/30d0)*z**3&
             -(1d0/144d0)*z**4 +(1d0/840d0)*z**5 - (1d0/5760d0)*z**6
        !(1d0/45360d0)*z**7 -(1d0/403200d0)*z**8      
   else

if (.false.)then !double


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

else !!quad

y=cmplx(realpart(z),imagpart(z),16)


      if (n.eq.0) then
      rez16=exp(-y)*y**(-1)
      else
      f16(0)=y**(-1)
      f16(1)=(y+1q0)*y**(-2)


      do i=2,n
        f16(i)=f16(i-2)+real(2*i-1,16)*f16(0)*f16(i-1)
      enddo

      rez16=exp(-y)*f16(n)
      endif
   rez=cmplx(realpart(rez16),imagpart(rez16),8)

endif !!end quad
endif

end subroutine 
