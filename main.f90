program atomHF

implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0, alpha=0.5d0
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:,:),rho(:),vfull1(:),ftemp1(:),&
        ftemp2(:),vx_phi(:,:),vx_phi1(:,:),vx_psidot(:,:),psidot(:,:),psi_non_norm(:,:),norm_arr(:)
integer, allocatable :: shell_n(:),shell_l(:),count_l(:)
real(8), allocatable :: shell_occ(:),psi_eig(:)

!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl =50 !Maximal itteration number
integer :: ir,il,icl,il_icl,iscl,lmax

real(8) :: Z,rez
real(8) :: norm,Rmin,Rmax,hh,dE_min,e1,e2,e3,energy,energy0,e_kin,e_ext,e_h,e_x,e_pot
integer :: grid,Ngrid, Nshell, ish
integer :: i,j,countl0, version
   
logical :: E_dE_file_exists, file_exists
character(len=1024) :: filename

dE_min=1d-7


!read input
read(*,*) 
read(*,*) Z, Rmin, Rmax, Ngrid, version
read(*,*)
read(*,*) grid
read(*,*) 
read(*,*) Nshell
allocate(shell_n(Nshell),shell_l(Nshell),shell_occ(Nshell),norm_arr(Nshell))
read(*,*)
do ish=1,Nshell
  read(*,*) shell_n(ish), shell_l(ish), shell_occ(ish)
enddo
#ifdef debug
  write(*,*) "shell_n ", shell_n(:)
  write(*,*) "shell_l ", shell_l(:)
  write(*,*) "shell_occ ", shell_occ(:)
#endif

lmax=maxval(shell_l)
allocate(count_l(lmax+1))


!counts how many particular l shells (there must be a shorter way to do it)
do i=0,lmax
 countl0=0
 do j=1,Nshell
 if (i.eq.shell_l(j)) then 
         countl0=countl0+1
 endif
 enddo
 write(*,*)"l=",i," countl=", countl0
 count_l(i+1)=countl0
enddo
 write(*,*) "count_l=",count_l(:)



allocate(r(Ngrid),vfull(Ngrid),vh(Ngrid),vxc(Ngrid),exc(Ngrid),eig(Ngrid),rho(Ngrid),vfull1(Ngrid))
allocate(psi(Ngrid,Nshell),psi_eig(Nshell))
allocate(ftemp1(Ngrid),ftemp2(Ngrid))

call gengrid(grid,Ngrid,Rmin,Rmax,r)
hh=r(2)-r(1)


! Initial guess for potential
!...
vfull=0*r !vfull is vh+vxc, core potential -Z/r is used seperatly - version 1 has to be modified

  inquire(file='norm.out',EXIST=E_dE_file_exists)
  if (E_dE_file_exists) then
     open(11,file='norm.out',status='old', access='append')
  else
     open(11,file='norm.out',status='new')
  endif
  write(11,*)"shell norm"
  close(11)

  inquire(file='E_dE.out',EXIST=E_dE_file_exists)
  if (E_dE_file_exists) then
     open(11,file='E_dE.out',status='old', access='append')
  else
     open(11,file='E_dE.out',status='new')
  endif
  write(11,*)"iter E dE"
  close(11)



if (version.eq.1) then
allocate(H(Ngrid,Ngrid))
! ----- version 1 -------
! Self-consistency loop
write(*,*)"************* version 1 **************"

do iscl=1,1!maxscl
#ifdef debug

  write(*,*) "iteration: ",iscl
#endif

il_icl=0
do il=1,lmax+1
  ! Construct Hamiltonian 
  
  call genham(Ngrid,r,vfull+0.5d0*dble(il-1)*dble(il)*r**(-2d0),H)

  ! Diagonalize
  call diasym(Ngrid,H,eig)

  ! Normalize WFs
  do icl=1,count_l(il)
    il_icl=il_icl+1
    psi(:,il_icl)=H(:,icl)  !tas ir Psi*r
    psi(:,il_icl)=psi(:,il_icl)/r
    psi_eig(il_icl)=eig(icl)
    call wfnorm(Ngrid,r,psi(:,il_icl))
#ifdef debug
    write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," normalised to:",4d0*Pi*sum(psi(:,il_icl)**2d0*r**2d0*hh),&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)
#endif
  enddo
enddo
!  write 3 wave functions to file 
!  open(11,file='wf_v1.dat',status='replace')
!  write(11,*)"r psi1 psi2"
!   do i = 1,Ngrid 
!     write(11,*)r(i),psi(i,1),psi(i,2),psi(i,3)
!   end do
!   close(11)


! Calculate density
  rho=0d0*rho

  
  
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2d0
  enddo

  ! Construct the Hartree potential

!  call getvh(Ngrid,r,rho,vh)
  call getvhtrapez(Ngrid,r,rho,vh)
! Construct the exchange-correlation potential
  call getvxc(Ngrid,rho,vxc,exc)

! Caulculate energy
  e1=0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
  e2=0
  e3=0
  do ir=1, Ngrid-1
    hh=r(ir+1)-r(ir)
    e2=e2-0.5d0*4d0*Pi*hh*vh(ir)*rho(ir)*r(ir)**2d0
    e3=e3+4d0*Pi*hh*(exc(ir)-vxc(ir))*rho(ir)*r(ir)**2d0
  enddo
  energy0=energy
  energy=e1+e2+e3

  Print *,"(11) E=",e1,"+",e2,"+",e3,"=",energy
 
  vfull1=vh+vxc
  vfull=alpha*vfull1+(1d0-alpha)*vfull
enddo

deallocate(H)



else if (version.eq.2) then
        
write(*,*)"********** version 2 *************"
! ----- version 2 -------
! Self-consistency loop
do iscl=1,maxscl
  il_icl=0
  do il=1,lmax+1
    ! Diagonalize via the shooting method
    call diashoot(Ngrid,r,Z,vfull,il-1,count_l(il),il_icl,Nshell,eig,psi)

    ! Normalize WFs
    do icl=1,count_l(il)
      il_icl=il_icl+1
      psi_eig(il_icl)=eig(icl)
      call integ_sph_s38_value(Ngrid,r,psi(:,il_icl)**2d0,norm)
      psi(:,il_icl)=psi(:,il_icl)*norm**(-0.5d0)
      write(*,*)"NORM=",norm
      if (il_icl.eq.2) then
        open(11,file='norm.out',status='old', access='append')
        write(11,*)il_icl, norm, "v",version       
        close(11)
      endif
#ifdef debug
      norm=0
      do ir=1,Ngrid-1
        hh=r(ir+1)-r(ir)
        norm=norm+psi(ir,il_icl)**2d0*r(ir)**2d0*hh
      enddo
      write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," normalised to:",norm,&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)

#endif
    enddo
  enddo

  !write 3 wave functions to file 
!  open(11,file='wf_v2_Newton_new.dat',status='replace')
!  write(11,*)"r psi1 psi2"
!   do i = 1,Ngrid 
!     write(11,*)r(i),psi(i,1),psi(i,2),psi(i,3)
!   end do
!   close(11)

 
 
  ! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2d0
  enddo

  ! Construct the Hartree potential
!  call getvhtrapez(Ngrid,r,rho,vh)
!  call getvhsimp38(Ngrid,r,rho,vh)

  call integ_s38_fun(Ngrid,r,r**2*rho,1,ftemp1)
  call integ_s38_fun(Ngrid,r,r*rho,-1,ftemp2)
  vh=ftemp1/r+ftemp2


  ! Construct the exchange-correlation potential
  call getvxc(Ngrid,rho/(4d0*Pi),vxc,exc)
  ! Caulculate energy
  e1=0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
 
!  e2=0
!  e3=0
!  do ir=1, Ngrid-1
!    hh=r(ir+1)-r(ir)
!    e2=e2-0.5d0*hh*vh(ir)*rho(ir)*r(ir)**2d0
!    e3=e3+hh*(exc(ir)-vxc(ir))*rho(ir)*r(ir)**2d0
!  enddo
  call integ_sph_s38_value(Ngrid,r,-0.5d0*vh*rho,e2)
  call integ_sph_s38_value(Ngrid,r,(exc-vxc)*rho,e3) 
  energy0=energy
  energy=e1+e2+e3


  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0
  close(11)

  Print *,iscl,".iter E=",e1,"+",e2,"+",e3,"=",energy
  if (abs(energy-energy0).LT.dE_min) then
          EXIT
  endif
  vfull1=vh+vxc
  vfull=alpha*vfull1+(1-alpha)*vfull

enddo

elseif(version.eq.3) then
allocate(vx_phi(Ngrid,Nshell),vx_phi1(Ngrid,Nshell),vx_psidot(Ngrid,Nshell),psidot(Ngrid,Nshell)&
        ,psi_non_norm(Ngrid,Nshell))
vx_phi=0d0*vx_phi
vx_psidot=0d0*vx_psidot
write(*,*)"********** version 3 *************"
! ----- version 3 HF_exchange -------
! Self-consistency loop
do iscl=1,maxscl
  il_icl=0
  do il=1,lmax+1
! Diagonalize via the shooting method


  call diashoot2(Ngrid,r,Z,vfull,il-1,count_l(il),il_icl,Nshell,vx_phi,vx_psidot,psidot,psi_eig,psi)



! Normalize WFs
    do icl=1,count_l(il)
      il_icl=il_icl+1
      psi_non_norm(:,il_icl)=psi(:,il_icl)
      call integ_sph_s38_value(Ngrid,r,psi(:,il_icl)**2d0,norm)
      norm_arr(il_icl)=norm
      psi(:,il_icl)=psi(:,il_icl)*norm**(-0.5d0)


#ifdef debug
!  norm=0
!  do ir=1,Ngrid-1
!    hh=r(ir+1)-r(ir)
!    norm=norm+psi(ir,il_icl)**2d0*r(ir)**2d0*hh
!  enddo
  write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," norm:",norm,&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)

#endif
    enddo
  enddo

  !write 3 wave functions to file 
!  write (filename, "(A2,I1)") "Be", iscl
!  print *, trim(filename)
!  open(11,file=filename,status='replace')
!  write(11,*)"r psi1 psi2, psi3"
!   do i = 1,Ngrid
!     write(11,*)r(i),psi_non_norm(i,1),psi_non_norm(i,2)!,psi_non_norm(i,3)
!   end do
!   close(11)




  ! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2d0
  enddo


  ! Construct the Hartree potential

  call integ_s38_fun(Ngrid,r,r**2*rho,1,ftemp1)
  call integ_s38_fun(Ngrid,r,r*rho,-1,ftemp2)
  vh=ftemp1/r+ftemp2


!  Caulculate vx_phi for every orbital  
!   call getvxc(Ngrid,rho/(4d0*Pi),vxc,exc)


  !write vh to file 
!  write (filename, "(A5,I1)") "vbe", iscl
!  print *, trim(filename)
!  open(11,file=filename,status='replace')
!  write(11,*)"r vh"
!   do i = 1,Ngrid
!     write(11,*)r(i),vh(i)
!   end do
!   close(11)



!  vxc=-0.5d0*vh !The line that works for He only
  do ish=1,Nshell
   call get_Fock_ex(Ngrid,r,ish,Nshell,shell_l,psi_non_norm(:,ish),psi,vx_phi1(:,ish))
 !  call get_Fock_ex(Ngrid,r,ish,Nshell,shell_l,psidot(:,ish),psi,vx_psidot(:,ish))



  !write vh to file 
!  write (filename, "(A5,I1,A1,I1)") "wf/vxpsi", iscl,"_",ish
!  print *, trim(filename)
!  open(11,file=filename,status='replace')
!  write(11,*)"r vx*psi vx_psi", iscl, ish 
!   do i = 1,Ngrid
!     write(11,*)r(i),psi_non_norm(i,ish)*vxc(i),vx_phi1(i,ish)
!     write(11,*)r(i),psidot(i,ish)*vxc(i),vx_psidot(i,ish)
!   end do
!   close(11)


!    vx_phi1(:,ish)=psi_non_norm(:,ish)*vxc
!    vx_psidot(:,ish)=psidot(:,ish)*vxc

  enddo
  vx_phi=alpha*vx_phi1+(1d0-alpha)*vx_phi
  !vx_phi=0d0*vx_phi
  !vx_psidot=psidot*0d0

  ! Caulculate energy
!Both give practicly the same result for e_ext
  call get_energy_ext(Ngrid,r,rho,Z,e_ext)
!  call integ_sph_s38_value(Ngrid,r,-rho*Z/r,e_ext)

  write(*,*)"e_ext=",e_ext
  call integ_sph_s38_value(Ngrid,r,0.5*rho*vh,e_h)
  write(*,*)"e_h=",e_h
  e_x=0d0
  do ish=1,Nshell
    !call integ_sph_s38_value(Ngrid,r,vxc*psi(:,ish)**2d0,e2)
    call integ_sph_s38_value(Ngrid,r,0.5d0*shell_occ(ish)*psi(:,ish)*vx_phi(:,ish)*norm_arr(ish)**(-0.5d0),e2)
    e_x=e_x+e2
  enddo
  write(*,*)"e_x=",e_x


  e1=0d0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
  e_kin=e1-e_ext-2d0*e_h-2d0*e_x
  write(*,*)"e_kin=",e_kin
  e_pot=e_ext+e_h+e_x

  write(*,*)"e_pot/e_kin=",e_pot/e_kin, "e_tot=",e_kin+e_pot, " (",iscl,")"

  energy0=energy
  energy=e_kin+e_pot

  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0
  close(11)

  call integ_sph_s38_value(Ngrid,r,-0.5d0*vh*rho,e2)
  call integ_sph_s38_value(Ngrid,r,(exc-vxc)*rho,e3)

  Print *,"Result in case of DFT vxc: ",iscl,".iter E=",e1,"+",e2,"+",e3,"=",e1+e2+e3


  if (abs(energy-energy0).LT.dE_min) then
          EXIT
  endif

  vfull1=vh!+vxc
  vfull=alpha*vfull1+(1d0-alpha)*vfull
!  vfull=vfull1
enddo
deallocate(vx_phi)

deallocate(vx_phi1,vx_psidot,psidot,psi_non_norm)
endif






! ----- version x3 -------
! Self-consistency loop
!do iscl=1,maxscl

! integrate the screened Poisson equation
! $\left( \nabla^2+ 2\epsilon \right \psi(\mathbf{r}) = v(\mathbf{r}) \psi(\mathbf{r})$

! Normalize WFs

! Calculate \epsilon as an expectation value

!enddo

! ----- version x4 -------
! Self-consistency loop
!do iscl=1,maxscl

! Davidson method 

!enddo

!  write results 

  inquire(file='results.dat',EXIST=file_exists)
  if (file_exists) then
     open(11,file='results.dat',status='old', access='append')
  else
     open(11,file='results.dat',status='new')
  endif
  write(11,*)"Z version iterations Ngrid Rmax Grid"
  write(11,*)Z, version, iscl-1, Ngrid, Rmax, grid
  write(11,*)"n l eigval"
  il_icl=0
  do il=1,lmax+1
    do icl=1, count_l(il)
      il_icl=il_icl+1
      write(11,*) icl, il, psi_eig(il_icl) 
    enddo
  enddo
  write(11,*)"Tot Energy: ", energy
  write(11,*)"dE        : ", energy-energy0
  write(11,*)"" 
  close(11)

  inquire(file='energy.out',EXIST=file_exists)
  if (file_exists) then
     open(11,file='energy.out',status='old', access='append')
  else
     open(11,file='energy.out',status='new')
     write(11,*)"Z iterations Ngrid Rmax Grid E dE"
  endif
  write(11,*)Z, iscl-1, Ngrid, Rmax, grid, energy, energy-energy0
  close(11)

  
  
  
  
deallocate(r,vfull,vh,vxc,exc,eig,psi,rho,shell_n,shell_l,count_l,shell_occ,norm_arr)

call wigner3j(3,3,2,rez)
write(*,*)rez
end program

