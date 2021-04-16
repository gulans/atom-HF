subroutine getvxc(Ngrid,rho,vxc,exc)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: rho(Ngrid) 
real(8), intent(out) :: vxc(Ngrid),exc(Ngrid)
integer :: i
real(8) :: rhoi, vxc_fun,exc_fun


do i=1, Ngrid
  rhoi=rho(i)
  if (rhoi.lt.(1e-20)) then
    vxc(i)=0
    exc(i)=0
  else
    vxc(i)=vxc_fun(rhoi)
    exc(i)=exc_fun(rhoi) 
  endif   
enddo
end subroutine

real(8) function vxc_fun(n)  !Funkcijas nosaukums jālieto par mainīgo???? Tiešām?
  IMPLICIT NONE
  real(8), intent(in)    :: n ! input
  real(8)  :: A,a1,b1,b2,b3,b4
  real(8) :: rs,vx,vc
  real(8), parameter :: pi = 3.1415926535897932384d0

  A=0.031091d0
  a1=0.21370d0
  b1=7.5957d0
  b2=3.5876d0
  b3=1.6382d0
  b4=0.49294d0

  rs = (3d0/(4d0*pi*n))**(1d0/3d0)
  vx = -(3d0/2d0/pi)**(2d0/3d0)/rs
  vc = -2d0*A*(1d0+2d0*a1*rs/3d0)*Log(1d0+1d0/(2*A*(b1*rs**0.5d0+b2*rs+b3*rs**1.5d0+b4*rs**2d0))) &
       -4d0*pi/9d0*(rs**4)*2d0*A*n*(1d0+a1*rs)/(1d0+1d0/(2*A*(b1*rs**0.5d0+b2*rs+b3*rs**1.5d0+b4*rs**2d0))) &
       *(2*A*(0.5d0*b1/(rs**0.5d0)+b2+1.5d0*b3*rs**0.5d0+2d0*b4*rs))/((2*A*(b1*rs**0.5d0+b2*rs+b3*rs**1.5d0+b4*rs**2d0))**2)

  vxc_fun=vx +vc

  return
end function

real(8) function exc_fun(n)  !Funkcijas nosaukums jālieto par mainīgo???? Tiešām?
  IMPLICIT NONE
  real(8), intent(in)    :: n ! input
  real(8)  :: A,a1,b1,b2,b3,b4
  real(8) :: rs,ex,ec
  real(8), parameter :: pi = 3.1415926535897932384d0

  A=0.031091d0
  a1=0.21370d0
  b1=7.5957d0
  b2=3.5876d0
  b3=1.6382d0
  b4=0.49294d0

  rs = (3d0/(4d0*pi*n))**(1d0/3d0)
  ex = -(9d0/4d0/pi**2)**(1d0/3d0)*3d0/4d0/rs
  ec = -2d0*A*(1d0+a1*rs)*Log(1d0+1d0/(2*A*(b1*rs**0.5d0+b2*rs+b3*rs**1.5d0+b4*rs**2d0)))
  exc_fun = ex + ec

  return
end function
