subroutine get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid,Nrsfun,lmax
real(8), intent(in) :: r(Ngrid)
complex(8),intent(in) :: rsfunC(Nrsfun,2) 
complex(8), intent(out) ::Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
complex(8) ::z,zmin,zmax, besi,besk
integer :: ir,l,fun
integer :: ii
real(8) :: remin,immin,remax,immax,zminabs,zmaxabs,zabs
  

  do l=0,lmax*2
  do fun=1,Nrsfun
  do ir = 1,Ngrid
    z=rsfunC(fun,2)*r(ir)
    call msbesselic (l, z, besi)
    call msbesselkc (l, z, besk)
    Bess_ik(ir,fun,l+1,1)=besi
    Bess_ik(ir,fun,l+1,2)=besk

  end do
  end do
  end do
close(11)
close(12)





! Debug
! write all argument values and bessel functions to files and print out range of Bessel function argument
if(.false.) then
       
ii=0
zminabs=1d100
zmaxabs=0d0
open(11,file='besi_all.dat',status='replace')
open(12,file='besk_all.dat',status='replace')
  do l=0,lmax*2
  do fun=1,Nrsfun
  do ir = 1,Ngrid
    ii=ii+1
    z=rsfunC(fun,2)*r(ir)
    zabs=dsqrt(realpart(z)**2+imagpart(z)**2)
    if (zabs.gt.zmaxabs)then
            zmaxabs=zabs
            zmax=z
    endif
    if (zabs.lt.zminabs)then
            zminabs=zabs
            zmin=z
    endif
    write(11,*)ii,l,r(ir),realpart(z),imagpart(z),realpart(Bess_ik(ir,fun,l+1,1)),imagpart(Bess_ik(ir,fun,l+1,1))
    write(12,*)ii,l,r(ir),realpart(z),imagpart(z),realpart(Bess_ik(ir,fun,l+1,2)),imagpart(Bess_ik(ir,fun,l+1,2))


  end do
  end do
  end do
close(11)
close(12)
write(*,*)"Minimum argument for Bessel", zmin
write(*,*)"Maximum argument for Bessel", zmax

endif
! Debug end



end subroutine



