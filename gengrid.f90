subroutine gengrid(grid,Z,Ngrid,Rmin,Rmax,r)
implicit none
integer, intent(in) :: grid,Ngrid
real(8), intent(in) :: Rmin,Rmax,Z 
real(8), intent(out) :: r(Ngrid)
integer :: ir

real(8) :: pow, Rmax1, Rmax2
integer :: ir1_start,ir1_stop,ir2_start,ir2_stop,ir3_start,ir3_stop
integer :: ir1, ir2, ir3

real(8) :: m1,m2,m3,p1,p2,p3,norm,Ngrid1,Ngrid2,Ngrid3

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

elseif (grid.eq.3) then !Advanced Cubic grid
  do ir=1,Ngrid
    r(ir)=ir*Rmin+(Rmax-Ngrid*Rmin)*(dble(ir-1)/dble(Ngrid-1))**3
  enddo

elseif (grid.eq.4) then !Advanced power 4 grid
  do ir=1,Ngrid
    r(ir)=ir*Rmin+(Rmax-Ngrid*Rmin)*(dble(ir-1)/dble(Ngrid-1))**4
  enddo
elseif (grid.eq.5) then !Advanced power 5 grid
  do ir=1,Ngrid
    r(ir)=ir*Rmin+(Rmax-Ngrid*Rmin)*(dble(ir-1)/dble(Ngrid-1))**5
  enddo
elseif (grid.eq.6) then !Advanced power 6 grid
  do ir=1,Ngrid
    r(ir)=ir*Rmin+(Rmax-Ngrid*Rmin)*(dble(ir-1)/dble(Ngrid-1))**6
  enddo
elseif (grid.eq.7) then !Advanced power 6 grid
  do ir=1,Ngrid
    r(ir)=ir*Rmin+(Rmax-Ngrid*Rmin)*(dble(ir-1)/dble(Ngrid-1))**7
  enddo


elseif (grid.eq.0) then !optimised polynomial grid
if ((Z.gt.2.5).and.(Z.lt.4.5)) then !Li(3) - Be(4) <1e-7 precision with Ngrid=80
pow=2.7d0
Rmax1=2d-1
Rmax2=2d0
m1=1d0
m2=1.5d0
m3=1d0

elseif ((Z.gt.4.5).and.(Z.lt.10.5)) then !B(5) - Ne(10) <1e-7 precision with Ngrid=130
pow=3.6d0
Rmax1=2d-2
Rmax2=2d0
m1=1d0
m2=1.8d0
m3=1d0

elseif ((Z.gt.10.5).and.(Z.lt.18.5)) then !Na(11) - Ar(18) <1e-6 precision with Ngrid=180 
pow=3.6d0
Rmax1=1d-2
Rmax2=2d0
m1=1d0
m2=2d0
m3=1d0

elseif ((Z.gt.18.5).and.(Z.lt.36.5)) then !K(19) - Kr(36) <1e-6 precision with Ngrid=200 for LDA and GGA, but Ngrid=230 for HF
pow=4.5d0
Rmax1=1d-2
Rmax2=2d0
m1=1d0
m2=2d0
m3=1d0


elseif ((Z.gt.36.5).and.(Z.lt.54.5)) then !Rb(37) - Xe(54) <1e-6 precision with Ngrid=230 
pow=4.4d0
Rmax1=1d-2
Rmax2=2d0
m1=1d0
m2=2d0
m3=1d0

elseif ((Z.gt.54.5).and.(Z.lt.86.5)) then !Cs(55) - Rn(86) <1e-6 precision with Ngrid=300 
pow=5.0d0
Rmax1=5d-3
Rmax2=3d0
m1=1d0
m2=2d0
m3=1d0
else
write(*,*)"There is no optimised grid for element Z= ",Z
pow=5.0d0
Rmax1=Rmax
Rmax2=Rmax
m1=1d0
m2=1d0
m3=1d0

endif
p1=(Rmax1/Rmax)**(1d0/pow)
p2=(Rmax2/Rmax)**(1d0/pow)-p1
p3=1d0-p1-p2
norm=p1*m1+p2*m2+p3*m3
m1=m1/norm
m2=m2/norm
m3=m3/norm

Ngrid1=dble(Ngrid)*m1
Ngrid2=dble(Ngrid)*m2
Ngrid3=dble(Ngrid)*m3

write(*,*)"Ngrids:",Ngrid1,Ngrid2,Ngrid3

ir1_start=1
ir1_stop=int(p1*Ngrid1)

ir2_start=int(p1*Ngrid2)+1
ir2_stop=int((p1+p2)*Ngrid2)

ir3_start=int((p1+p2)*Ngrid3)+1
ir3_stop=ir3_start+(Ngrid-(ir1_stop-ir1_start+1)-(ir2_stop-ir2_start+1))-1

write(*,*)"grid1: ",ir1_start,ir1_stop, "N=",ir1_stop-ir1_start+1
write(*,*)"grid2: ",ir2_start,ir2_stop, "N=",ir2_stop-ir2_start+1
write(*,*)"grid3: ",ir3_start,ir3_stop, "N=",ir3_stop-ir3_start+1
ir=1
do ir1=ir1_start,ir1_stop
r(ir)=Rmin*ir+(Rmax-Ngrid*Rmin)*(dble(ir1-1)/dble(Ngrid1-1))**pow
!write(*,*)"extir=",ir," ir1=",ir1,"r", r(ir)
ir=ir+1
enddo
write(*,*)
do ir2=ir2_start,ir2_stop
r(ir)=Rmin*ir+(Rmax-Ngrid*Rmin)*(dble(ir2-1)/dble(Ngrid2-1))**pow
!write(*,*)"extir=",ir," ir2=",ir2,"r", r(ir)
ir=ir+1
enddo
write(*,*)

do ir3=ir3_start,ir3_stop
r(ir)=Rmin*ir+(Rmax-Ngrid*Rmin)*(dble(ir3-1)/dble(Ngrid3-1))**pow
!write(*,*)"extir=",ir," ir3=",ir3,"r", r(ir)
ir=ir+1
enddo


endif




!write(*,*) "Grid:"
!  do ir=1,Ngrid
!    write(*,*) "ri=",ir," r=",r(ir)
!  enddo
 
end subroutine
