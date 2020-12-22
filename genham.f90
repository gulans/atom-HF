subroutine genham(Ngrid,r,vfull,H)
implicit none
integer, intent(in) :: Ngrid
real(8), intent(in) :: r(Ngrid),vfull(Ngrid) 
real(8), intent(out) :: H(Ngrid,Ngrid)
integer :: i,j
real(8) :: hh,Hij

hh=r(2)-r(1) !PAVIRSHI
do i=1,Ngrid
  do j=1,Ngrid
     if (j==i) then  !GALVENAA DIOGNAALE
       Hij=1d0/(hh**2d0)+vfull(i)
     else if (i+1==j) then !(augšējā diogn)
       Hij=-1d0/(2d0*hh**2d0)
     else if (i-1==j) then !(apakšējā diogn)
       Hij=-1d0/(2d0*hh**2d0)
    else
       Hij=0d0
     endif
     H(i,j)=Hij
  enddo
enddo



end subroutine
