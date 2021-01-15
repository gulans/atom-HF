subroutine gengrid(grid,Ngrid,Rmin,Rmax,r)
implicit none
integer, intent(in) :: grid,Ngrid
real(8), intent(in) :: Rmin,Rmax 
real(8), intent(out) :: r(Ngrid)
integer :: ir

!cubic exciting       spr (ir, is) = sprmin (is)+(dble(ir-1)/dble(nrmt(is)-1))**3*(rmt(is)-sprmin (is))
!exponential exciting   spr (ir, is) = sprmin (is) * Exp (dble(ir-1)*t1*t2)

if (grid.eq.0) then !Equidistant
  do ir=1,Ngrid
    r(ir)=Rmin+(Rmax-Rmin)*dble(ir-1)/dble(Ngrid-1)
  enddo

elseif (grid.eq.1) then !Cubic
  do ir=1,Ngrid
    r(ir)=Rmin+(Rmax-Rmin)*(dble(ir-1)/dble(Ngrid-1))**3
  enddo

elseif (grid.eq.2) then !Exponential
  do ir=1,Ngrid
    r(ir)=Rmin*exp(dble(ir-1)*LOG(Rmax/Rmin)/dble(Ngrid-1))
  enddo
endif
!write(*,*) "Grid:"
!  do ir=1,Ngrid
!    write(*,*) "ri=",ir," r=",r(ir)
!  enddo
 
end subroutine
