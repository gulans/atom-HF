subroutine get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: Ngrid,Nrsfun,lmax
real(8), intent(in) :: r(Ngrid)
complex(8),intent(in) :: rsfunC(Nrsfun+1,2) 
complex(8), intent(out) ::Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)
complex(8) :: besi,besk
integer :: ir,l,fun



  do l=0,lmax*2
  do fun=1,Nrsfun
  do ir = 1,Ngrid

    call msbesselic (l, rsfunC(fun,2)*r(ir), besi)
    call msbesselkc (l, rsfunC(fun,2)*r(ir), besk)
    Bess_ik(ir,fun,l+1,1)=besi
    Bess_ik(ir,fun,l+1,2)=besk
  end do
  end do
  end do

end subroutine

