subroutine wfnorm(Ngrid,r,psi)
implicit none
 real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid) 
real(8), intent(out) :: psi(Ngrid)
integer :: i,i_interp
real(8) :: hh,norm
real(8) :: x(4),y(4),psi1,psi2,f0,f1,f2,f3


 
norm=0
i=1

f3=4d0*Pi*r(i)**2d0*psi(i)**2d0 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  



do i=1,Ngrid-1
  hh=(r(i+1)-r(i))/3
  if (i.lt.2) then
    i_interp=1
  else if(i.gt.Ngrid-3) then
    i_interp=Ngrid-3
  else
    i_interp=i-1
  endif
  x=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3)/)
  y=(/psi(i_interp),psi(i_interp+1),psi(i_interp+2),psi(i_interp+3)/)
  call interp(4,x,y,r(i)+hh,psi1)
  call interp(4,x,y,r(i)+2d0*hh,psi2)

  f0=f3 !in numerical recepies (4.1.5) every iteration f0 is previous iterations f3  
  f1=4d0*Pi*(r(i)+hh)**2d0*psi1**2d0
  f2=4d0*Pi*(r(i)+2d0*hh)**2d0*psi2**2d0
  f3=4d0*Pi*r(i+1)**2d0*psi(i+1)**2d0

  norm=norm+hh*(3d0*f0+9d0*f1+9d0*f2+3d0*f3)/8d0
enddo

!!!! This was before !!!!!
!norm=0
!do ir=1, Ngrid-1
!  hh=r(ir+1)-r(ir)
!  norm=norm+4d0*Pi*r(ir)**2d0*psi(ir)**2d0*hh
!enddo

!
!call integ_sph_s38_value(Ngrid,r,psi**2d0,norm) 

psi=psi*norm**(-0.5d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine
