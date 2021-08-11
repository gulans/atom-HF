subroutine sew_function3(Ngrid,r,tools,tools_info,Z,l,r_sew,eig,psi,psip,psipp,psi_out,psip_out,psipp_out)

implicit none
integer(8), intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid)
integer, intent(in) :: tools_info(3)
real(8), intent(in) :: tools(Ngrid,tools_info(1))
real(8), intent(in) :: Z
integer, intent(in) :: l
real(8), intent(in) :: r_sew
real(8), intent(in) :: eig
real(8), intent(in) :: psi(Ngrid)
real(8), intent(in) :: psip(Ngrid)
real(8), intent(in) :: psipp(Ngrid)
real(8), intent(out) :: psi_out(Ngrid)
real(8), intent(out) :: psip_out(Ngrid)
real(8), intent(out) :: psipp_out(Ngrid)

real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)
real(8), PARAMETER :: lspeed=137.035999084d0
integer :: ir,ir_sew,i
real(8) :: a,psi_sew(Ngrid),psip_sew(Ngrid),psipp_sew(Ngrid)
real(8) :: b,c

do ir=1,Ngrid
if (r(ir).ge.r_sew)then
    ir_sew=ir
    exit
  endif
enddo

a=(-lspeed**2+sqrt(lspeed**4+lspeed**4*dble(l)+lspeed**4*dble(l)**2-lspeed**2*Z**2))/lspeed**2
b=(-2d0*lspeed**4 - eig*Z**2 + 2d0*lspeed**2*(-Z**2 + Sqrt(lspeed**4*(1d0 +dble(l)+dble(l)**2)&
        - lspeed**2*Z**2)))/(Z*(lspeed**2 + 2d0*Sqrt(lspeed**4*(1d0 +dble(l)+dble(l)**2) - lspeed**2*Z**2)))


c=psi(ir_sew)/(r(ir_sew)**a + b*r(ir_sew)**(a+1d0))

psi_sew =  c * (r**a                  + b*r**(a+1d0))
psip_sew=  c * (a*r**(a-1d0)          + b*(a+1d0)*r**a)
psipp_sew= c * (a*(a-1d0)*r**(a-2d0)  + b*(a+1d0)*a*r**(a-1d0))

psi_out=psi
psip_out=psip
psipp_out=psipp

psi_out(1:ir_sew)=psi_sew(1:ir_sew)
psip_out(1:ir_sew)=psip_sew(1:ir_sew)
psipp_out(1:ir_sew)=psipp_sew(1:ir_sew)

end subroutine




