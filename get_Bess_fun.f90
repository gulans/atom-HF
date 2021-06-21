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
open(11,file='bes/besi_all_sm.dat',status='replace')
open(12,file='bes/besk_all_sm.dat',status='replace')
ii=0
zminabs=1d100
zmaxabs=0d0
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

    call msbesselic (l, z, besi)
    write(11,*)ii,l,r(ir),realpart(z),imagpart(z),realpart(besi),imagpart(besi)
    call msbesselkc (l, z, besk)
    write(12,*)ii,l,r(ir),realpart(z),imagpart(z),realpart(besk),imagpart(besk)
    Bess_ik(ir,fun,l+1,1)=besi
    Bess_ik(ir,fun,l+1,2)=besk

  end do
  end do
  end do
close(11)
close(12)

write(*,*)"Minimum argument for Bessel", zmin
write(*,*)"Maximum argument for Bessel", zmax
end subroutine



