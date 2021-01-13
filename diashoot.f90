subroutine diashoot(Ngrid,r,vfull,l,num,shell0,nummax,eigval,eigfun)
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

integer, intent(in) :: Ngrid, l, num,nummax,shell0
real(8), intent(in) :: r(Ngrid),vfull(Ngrid)

real(8), intent(out) :: eigval(nummax)
real(8), intent(out) :: eigfun(Ngrid,nummax)


!--other variables--

integer :: i,ri,ei, junct,ji,try_dir,shell 
real(8) :: eigtry,psi(Ngrid),eigtry_max,eigtry_min,eigtry_step,psi0NR(Ngrid),psi1NR(Ngrid),psidot(Ngrid),emaxNR,eminNR 
logical ::eigtry_max_OK,eigtry_min_OK,euler
euler=.false.
! try_dir: -1 if we have to try smaller eigtry, +1 - larger, 0 - we found the eigenvalue

eigtry=-20d0
eigtry_step=10d0


!call shoot_using_Euler(Ngrid,r,vfull,0,1,-50.0d0,psi,try_dir,.TRUE.,10.0d0)



do ei=1,num 
  !write(*,*)"shell0=",shell0," ei=",ei
  shell=ei+shell0
  eigtry_max_OK=.false.
  eigtry_min_OK=.false.


  if (ei.NE.1) then
   eigtry_min_OK=.TRUE.
   eigtry_min=eigval(ei-1)
   eigtry=eigval(ei-1)+eigtry_step
  endif

  
  !Searching_for bisection minimum and maximum
  do i=1,1000 !to avoid deadlock
    call shoot_using_Euler(Ngrid,r,vfull,l,ei,eigtry,psi,try_dir,.FALSE.,eigtry_max-eigtry_min,euler)
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
  emaxNR=0
  eminNR=0
  do while((try_dir.NE.0))
    eigtry=(eigtry_max+eigtry_min)/2d0
    call shoot_using_Euler(Ngrid,r,vfull,l,ei,eigtry,psi,try_dir,.FALSE.,eigtry_max-eigtry_min,euler)
    if (try_dir.EQ.-1) then
      eigtry_max=eigtry
    else if (try_dir.EQ.1) then
      eigtry_min=eigtry
    endif
   ! write(*,*)"e_max",eigtry_max," eigmin",eigtry_min

    if(((eigtry_max-eigtry_min).LT.1d-9).AND.(emaxNR.eq.0).AND.(eminNR.eq.0)) then !will be prepared to use Newton–Raphson method
      
      emaxNR=eigtry_max
      eminNR=eigtry_min
 !     emaxNR=-51.914435599464900d0
 !     eminNR=-51.914435598882800d0
     ! write(*,*)"Preparing for Newton–Raphson"
    !for saving time we can put an exit here use Newton–Raphson mathod and remove the next if statement
!    endif
    
!    if((eigtry_max-eigtry_min).LT.1d-13) then
!      write(*,*)"FAILED to get a solution l=",l," ei=",ei," eigtry=",eigtry,&
!              " e_max-e_min=",eigtry_max-eigtry_min," psi(Ngrid)=",psi(Ngrid)

      call shoot_using_Euler(Ngrid,r,vfull,l,ei,eminNR,psi0NR,try_dir,.TRUE.,0d0,euler)

      call shoot_using_Euler(Ngrid,r,vfull,l,ei,emaxNR,psi1NR,try_dir,.TRUE.,0d0,euler)
!This is the new way to get psidot
      call getpsidot_euler(Ngrid,r,vfull,l,psi0NR,eminNR,psidot) !

!This is the old way to get psidot
!      psidot=(psi1NR-psi0NR)/(emaxNR-eminNR)

      call NR_method(Ngrid,eminNR,psi0NR,psi1NR,psidot,eigtry,psi)
      EXIT
    endif
  enddo
  write(*,*)"try_dir=",try_dir," Psi(Ngrid)=",psi(ngrid)," Bisection range etry_max-etry_min=",eigtry_max-eigtry_min
  eigfun(:,shell)=psi
  eigval(ei)=eigtry
!  write(*,*)"ei=",ei," eigval=",eigtry," e_max-e_min=",eigtry_max-eigtry_min," psi(Ngrid)=",psi(Ngrid)
 
enddo
end subroutine

subroutine getpsidot_euler(Ngrid,r,vfull,l,psi0,e0,psidot) 

integer, intent(in) :: Ngrid,l
real(8), intent(in) :: e0,r(Ngrid),vfull(Ngrid),psi0(Ngrid)
real(8), intent(out) :: psidot(Ngrid)

real(8) :: u(Ngrid),v(Ngrid),s(Ngrid),vprime(Ngrid)
integer ri
s(1)=0 !boundary condition
v(1)=0 !it is s'

u=psi0*r

vprime(1)=(2d0*vfull(1)+dble(l)*dble(l+1)*r(1)**(-2d0)-2*e0)*s(1)-2*u(1) 
psidot(1)=s(1)/r(1)

do ri=2,Ngrid
  h=r(ri)-r(ri-1)
  s(ri)=s(ri-1)+v(ri-1)*h
  v(ri)=v(ri-1)+vprime(ri-1)*h
  vprime(ri)=(2d0*vfull(ri)+dble(l)*dble(l+1)*r(ri)**(-2d0)-2*e0)*s(ri)-2*u(ri)
  psidot(ri)=s(ri)/r(ri)
!  write (*,*)ri,s(ri),v(ri),u(ri),vprime(ri),psidot(ri)
enddo


end subroutine



subroutine NR_method(Ngrid,e0,psi0,psi1,psidot,eig,psi)
integer, intent(in) :: Ngrid
real(8), intent(in) :: e0,psi0(Ngrid),psi1(Ngrid),psidot(Ngrid)
real(8), intent(out) :: psi(Ngrid),eig

real(8) :: psi_rmax 
integer :: i_rmax,ri,last_ri
logical :: even_junct

i_rmax=0
psi_rmax=1 !This parameter changes results.
do i=1,Ngrid
  if ((abs(psi0(i)).GT.psi_rmax).or.(abs(psi1(i)).GT.psi_rmax)) then
    i_rmax=i
    EXIT
  endif 
enddo
if (i_rmax.eq.0) then
        i_rmax=Ngrid
endif
eig=e0-psi0(i_rmax)/psidot(i_rmax)

psi=psi0-(psi0(i_rmax)/psidot(i_rmax))*psidot
write(*,*)"Eigval from Newton–Raphson=",eig," New i_rmax=",i_rmax," psi(i_rmax)=",psi(i_rmax) 
last_ri=i_rmax
do ri=last_ri+1, Ngrid
  psi(ri)=0d0
enddo
end subroutine



subroutine shoot_using_Euler(Ngrid,r,vfull,l,ei,eigtry,psi,try_dir,finish,bis_range,Euler)
implicit none

integer, intent(in) :: Ngrid,l,ei
real(8), intent(in) :: r(Ngrid),vfull(Ngrid),eigtry
logical, intent(in) :: finish,Euler
integer, intent(out) :: try_dir 
real(8), intent(out) :: psi(Ngrid)

real(8) :: junct,ji,u(Ngrid),v(Ngrid),vprime(Ngrid),h,r_interp(4),v_interp(4),vfull_interp
integer :: ri,i_interp
real(8) :: bis_range,zerro_effective,minimum_bis_range
real(8) :: k1u,k2u,k3u,k4u,k1v,k2v,k3v,k4v
zerro_effective=1d-4
minimum_bis_range=1d-9
! u - Psi*r 
! v - du/dr
  

try_dir=2          
junct=ei-1       !how many junctions with x-axis WF must have
ji=0
!write(*,*)"junctions=", junct," l=",l," ei=",ei
u(1)=r(1)**dble(l+1d0)   !boundary condition
v(1)=dble(l+1)*r(1)**dble(l)
vprime(1)=(2d0*vfull(1)+dble(l)*dble(l+1)*r(1)**(-2d0)-2d0*eigtry)*u(1)
psi(1)=u(1)/r(1)

!  open(11,file='RK_rez.dat',status='replace')
!   write(11,*)"r v v_interp k1u k2u k3u k4u k1v k2v k3v k4v psi"


do ri=1,Ngrid-1
  h=r(ri+1)-r(ri)

if (Euler) then
!Euler
  u(ri+1)=u(ri)+v(ri)*h
  v(ri+1)=v(ri)+vprime(ri)*h
  vprime(ri+1)=(2d0*vfull(ri+1)+dble(l)*dble(l+1)*r(ri+1)**(-2d0)-2d0*eigtry)*u(ri+1)
else
!RK4
!interpolation to obtain point vfull(ri+1/2) needed
  if (ri.lt.2) then
  i_interp=1
  else if(ri.gt.Ngrid-3) then
  i_interp=Ngrid-3
  else
  i_interp=ri-1
  endif

  r_interp=(/r(i_interp),r(i_interp+1),r(i_interp+2),r(i_interp+3)/)
  v_interp=(/vfull(i_interp),vfull(i_interp+1),vfull(i_interp+2),vfull(i_interp+3)/)
!    call interp(4,r_interp,v_interp,r(ri)+h/2d0,vfull_interp)
  vfull_interp=-4d0/(r(ri)+h/2d0)
  if (ri.lt.6) then
!          write(*,*)"ri ",ri
!          write(*,*)"r_interp ",r_interp
!          write(*,*)"v_interp ",v_interp
!          write(*,*)"v(",r(ri)+h/2,")=",vfull_interp
  endif


  k1u=v(ri)
  k1v=(2d0*vfull(ri)+dble(l)*dble(l+1)*r(ri)**(-2d0)-2d0*eigtry)*u(ri)
  k2u=v(ri)+h*k1v/2d0
  k2v=(2d0*vfull_interp+dble(l)*dble(l+1)*(r(ri)+h/2d0)**(-2d0)-2d0*eigtry)*(u(ri)+h*k1u/2)
  k3u=v(ri)+h*k2v/2d0
  k3v=(2d0*vfull_interp+dble(l)*dble(l+1)*(r(ri)+h/2d0)**(-2d0)-2d0*eigtry)*(u(ri)+h*k2u/2)
  k4u=v(ri)+h*k3v
  k4v=(2d0*vfull(ri+1)+dble(l)*dble(l+1)*(r(ri+1))**(-2d0)-2d0*eigtry)*(u(ri)+h*k3u)
  u(ri+1)=u(ri)+h*(k1u+2d0*k2u+2d0*k3u+k4u)/6d0
  v(ri+1)=v(ri)+h*(k1v+2d0*k2v+2d0*k3v+k4v)/6d0
! End of RK4

endif
psi(ri+1)=u(ri+1)/r(ri+1)
  !write(*,*)r(ri),u(ri),v(ri),vprime(ri),psi(ri)

!  write(11,*)r(ri), vfull(ri), vfull_interp, k1u, k2u, k3u, k4u, k1v, k2v, k3v, k4v, psi(ri)


!Check if the WF crossed x axis
  if ((psi(ri+1)*psi(ri)).LE.0d0) then
     !if ((psi(ri+1).EQ.0d0) or (psi(ri).EQ.0d0)) !something with low probability but should be cheked
     ji=ji+1
!     write(*,*)"junction r=",r(ri+1)
     if (ji.GT.junct) then
        try_dir=-1 !we have to try smaller eigtry 
        
        if (.not.(finish)) then 
                EXIT
        endif
     
     endif
  endif

  !if the wf>1 - then it will not come back to 0
  if (abs(psi(ri+1)).GT.1d0) then 
        try_dir=1 !we have to try larger eigtry 
        if (.not.(finish)) then 
                EXIT
        endif
     endif

enddo

close(11)

if (try_dir.eq.2)  then
        
  if (ji.EQ.junct) then
    if (((abs(psi(Ngrid)).LT.zerro_effective)).AND.(bis_range.LT.minimum_bis_range)) then    !check if psi(rmax) is close enough to zerro
      try_dir=0
    else
      try_dir=1
    endif
  endif
endif
!write(*,*)"eigtry=",eigtry," ri",ri," try_dir=",try_dir," psi(Ngrid)="&
!        ,psi(Ngrid)," bis_range: ",bis_range

end subroutine

subroutine interp(m,x,y,xp,yp)
 implicit none
 integer, intent(in) :: m
 real(8), intent(in) :: x(m),y(m),xp
 real(8), intent(out):: yp

 integer :: i, j
 real(8) :: sk,sauc
 yp=0
 do i=1,m
   sk=1.0d0
   sauc=1.0d0
   do j=1,m
     if (i.ne.j) then
       sk=sk*(xp-x(j))
       sauc=sauc*(x(i)-x(j))

     endif
 enddo
 yp=yp+sk*y(i)/sauc
 enddo


end subroutine
