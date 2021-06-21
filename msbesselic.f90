
subroutine msbesselic (n, z, rez)
       Implicit None
! arguments
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: z
      Complex (8), Intent (Out) :: rez
      Integer :: i
      Complex (8) :: g(-n-1:n)


      Complex (16) :: g16(-n-1:n),y,rez16


      Complex (8) :: rezi,rezk
      logical :: cits
      cits=.FALSE.
      !cits=.True.
if (cits) then
      !call bes_ik_c( n, z, rezi,rezk )
      rez=rezi
else
  if ((n.eq.0).and. (abs(realpart(z)).lt.1d-3)) then
          rez=1d0+(1d0/6d0)*z**2 + (1d0/120d0)*z**4 + (1d0/5040d0)*z**6  

  elseif ((n.eq.1) .and. (abs(realpart(z)).lt.1d-1).and.(abs(imagpart(z)).lt.1d-1) ) then
          rez=(1d0/3d0)*z +(1d0/30d0)*z**3 + (1d0/840d0)*z**5 + (1d0/45360d0)*z**7+&
                  (1d0/3991680d0)*z**9

  elseif ((n.eq.2) .and. (abs(realpart(z)).lt.0.5d0).and.(abs(imagpart(z)).lt.0.5d0))  then

          rez=(1d0/15d0)*z**2 +(1d0/210d0)*z**4 + (1d0/7560d0)*z**6 + (1d0/498960d0)*z**8 +&
                  (1d0/51891840d0)*z**10+ (1d0/7783776000d0)*z**12+ (1d0/1587890304000d0)*z**14

  elseif ((n.eq.3).and. (abs(realpart(z)).lt.1d-1)) then
          rez=(1d0/105d0)*z**3 + (1d0/1890d0)*z**5 + (1d0/83160d0)*z**7&
                 + (1d0/6486480d0)*z**9 
  elseif ((n.eq.4).and. (abs(realpart(z)).lt.1d-1)) then
          rez=(1d0/945d0)*z**4 + (1d0/20790d0)*z**6 + (1d0/1081080d0)*z**8&
                +(1d0/97297200d0)*z**10  
  elseif ((n.eq.5).and. (abs(realpart(z)).lt.1d-1)) then
          rez=(1d0/10395d0)*z**5 + (1d0/270270d0)*z**7 + (1d0/16216200d0)*z**9&
                 +(1d0/1654052400d0)*z**11 
  elseif ((n.eq.6).and. (abs(realpart(z)).lt.1d-1)) then
          rez=(1d0/135135d0)*z**6 + (1d0/4054050d0)*z**8 + (1d0/275675400d0)*z**10&
                  +(1d0/31426995600d0)*z**12

  else
  
  if (.false.)then !double


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
 
  else !quad
      y=cmplx(realpart(z),imagpart(z),16)

      if (n.eq.0) then
      rez16=sinh(y)*y**(-1)
      else
      g16(0)=y**(-1)
      g16(1)=-y**(-2)


      do i=2,n
        g16(i)=g16(i-2)-real(2*i-1,16)*g16(0)*g16(i-1)
      enddo
      do i=-1,-n-1,-1
        g16(i)=g16(i+2)+real(2*i+3,16)*g16(0)*g16(i+1)
      enddo

      rez16=g16(n)*sinh(y)+g16(-n-1)*cosh(y)
      endif
   rez=cmplx(realpart(rez16),imagpart(rez16),8)


  endif !end quad

  endif
endif
end subroutine 
