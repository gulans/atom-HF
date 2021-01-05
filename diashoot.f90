subroutine diashoot(Ngrid,r,vfull,l,num,nummax,eigval,eigfun)
        ! Ngrid
        ! r
        ! vfull - potential (v_n+v_h+v_xc)
        ! l - quantum number l
        ! num - Number Of eigvals and eigfun to solve
        ! nummax - size of eigval array
        ! eigval (OUT) - erigval array
        ! eigfun (OUT) - eigfun array


!--input and output variables--
implicit none

integer, intent(in) :: Ngrid, l, num,nummax
real(8), intent(in) :: r(Ngrid),vfull(Ngrid)

real(8), intent(out) :: eigval(nummax)
real(8), intent(out) :: eigfun(Ngrid,nummax)


!--other variables--

integer :: i,ri,ei, junct,ji,try_dir 
real(8) :: eigtry,psi(Ngrid),eigtry_max,eigtry_min,eigtry_step 
logical ::eigtry_max_OK,eigtry_min_OK

! try_dir: -1 if we have to try smaller eigtry, +1 - larger, 0 - we found the eigenvalue

eigtry=-20d0
eigtry_step=10d0





do ei=1,num 

  eigtry_max_OK=.false.
  eigtry_min_OK=.false.


  if (ei.NE.1) then
   eigtry_min_OK=.TRUE.
   eigtry_min=eigval(ei-1)
   eigtry=eigval(ei-1)+eigtry_step
  endif

  
  !Searching_for bisection minimum and maximum
  do i=1,1000 !to avoid deadlock
    call shoot_using_Euler(Ngrid,r,vfull,l,ei,eigtry,psi,try_dir)
    if (try_dir.EQ.-1) then
      eigtry_max_OK=.TRUE. 
      eigtry_max=eigtry
      eigtry=eigtry-eigtry_step
    else if (try_dir.EQ.1) then
      eigtry_min_OK=.TRUE.
      eigtry_min=eigtry
      eigtry=eigtry+eigtry_step
    end if
    if (eigtry_max_OK.and.eigtry_min_OK) then
      EXIT      
    endif
  enddo
  if (i.gt.1000) write(*,*)"Bisection min/max not found eig_min=",eigtry_min," eig_max=",eigtry_max," ei=",ei," l=",l;
!  write(*,*)"Starting bisection eig_min=",eigtry_min," eig_max=",eigtry_max
  
  !starting_bisection

  do while((try_dir.NE.0))
    eigtry=(eigtry_max+eigtry_min)/2
    call shoot_using_Euler(Ngrid,r,vfull,l,ei,eigtry,psi,try_dir)
    if (try_dir.EQ.-1) then
      eigtry_max=eigtry
    else if (try_dir.EQ.1) then
      eigtry_min=eigtry
    endif
!    write(*,*)"e_max",eigtry_max," eigmin",eigtry_min
    if((eigtry_max-eigtry_min).LT.1e-14) then
      write(*,*)"FAILED to get a solution l=",l," ei=",ei," eigtry=",eigtry,&
              " e_max-e_min=",eigtry_max-eigtry_min," psi(Ngrid)=",psi(Ngrid)
      EXIT
    endif
  enddo
  eigfun(:,ei)=psi
  eigval(ei)=eigtry
!  write(*,*)"ei=",ei," eigval=",eigtry," e_max-e_min=",eigtry_max-eigtry_min," psi(Ngrid)=",psi(Ngrid)
 
enddo
end subroutine

subroutine shoot_using_Euler(Ngrid,r,vfull,l,ei,eigtry,psi,try_dir)
implicit none

integer, intent(in) :: Ngrid,l,ei
real(8), intent(in) :: r(Ngrid),vfull(Ngrid),eigtry

integer, intent(out) :: try_dir 
real(8), intent(out) :: psi(Ngrid)

real(8) :: junct,ji,u(Ngrid),v(Ngrid),vprime(Ngrid),h
integer :: ri
real(8) :: zerro_effective

zerro_effective=1e-9
! u - Psi*r 
! v - du/dr
  

try_dir=2          !something is wrong if it will stay 2 after this subroutine
junct=ei-1       !how many junctions with x-axis WF must have
ji=0
!write(*,*)"junctions=", junct," l=",l," ei=",ei
u(1)=r(1)**dble(l+1d0)   !boundary condition
v(1)=(l+1)*r(1)**l
vprime(1)=(2d0*vfull(1)+l*(l+1)*r(1)**(-2)-2d0*eigtry)*u(1)
psi(1)=u(1)/r(1)

do ri=2,Ngrid
  h=r(ri)-r(ri-1)
  u(ri)=u(ri-1)+v(ri-1)*h
  v(ri)=v(ri-1)+vprime(ri-1)*h
  vprime(ri)=(2d0*vfull(ri)+l*(l+1d0)*r(ri)**(-2)-2d0*eigtry)*u(ri)
  psi(ri)=u(ri)/r(ri)
  !write(*,*)r(ri),u(ri),v(ri),vprime(ri),psi(ri)


!Check if the WF crossed x axis
  if ((psi(ri)*psi(ri-1)).LE.0d0) then
     !if ((psi(ri).EQ.0d0) or (psi(ri-1).EQ.0d0)) !something with low probability but should be cheked
     ji=ji+1
!     write(*,*)"junction r=",r(ri)
     if (ji.GT.junct) then
        try_dir=-1 !we have to try smaller eigtry 
        EXIT
     
     endif
  endif

  !if the wf>1 - then it will not come back to 0
  if (psi(ri).GT.1) then
        try_dir=1 !we have to try larger eigtry 
        EXIT
     endif

enddo

if (try_dir.eq.2)  then
        
  if (ji.EQ.junct) then
    if ((abs(psi(Ngrid)).LT.zerro_effective)) then    !check if psi(rmax) is close enough to zerro
      try_dir=0
    else
      try_dir=1
    endif
  endif
endif
write(*,*)"eigtry=",eigtry," ir=",ri," try_dir=",try_dir," psi(Ngrid)=" ,psi(Ngrid)

end subroutine
