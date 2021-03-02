program atomHF

implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0, alpha=0.5d0
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:,:),rho(:),vfull1(:),ftemp1(:),&
        ftemp2(:),vx_phi(:,:),vx_phi1(:,:),vx_psidot(:,:),psidot(:,:),psi_non_norm(:,:),norm_arr(:),&
        vhp(:),vxcp(:)
integer, allocatable :: shell_n(:),shell_l(:),count_l(:)
real(8), allocatable :: shell_occ(:),psi_eig(:),psi_eig_temp(:)
real(8), allocatable :: psip(:,:),eigp(:)
!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl =50 !Maximal itteration number
integer :: ir,il,icl,il_icl,iscl,lmax,l_n,inn

real(8) :: Z,rez
real(8) :: norm,Rmin,Rmax,hh,dE_min,e1,e2,e3,energy,energy0,e_kin,e_ext,e_h,e_x,e_pot
real(8) :: e11,e22,e33
integer :: grid,Ngrid, Nshell, ish
integer :: i,j,countl0, version
   
logical :: E_dE_file_exists, file_exists
character(len=1024) :: filename

Real (8)  :: besrez (0:50)


dE_min=1d-8


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
allocate(psi(Ngrid,Nshell),psi_eig(Nshell),psi_eig_temp(Nshell))
allocate(ftemp1(Ngrid),ftemp2(Ngrid))

call gengrid(grid,Ngrid,Rmin,Rmax,r)
hh=r(2)-r(1)


! Initial guess for potential
!...
vfull=0d0*r !vfull is vh+vxc, core potential -Z/r is used seperatly - version 1 has to be modified

  inquire(file='eigval_de.out',EXIST=E_dE_file_exists)
  if (E_dE_file_exists) then
     open(11,file='eigval_de.out',status='old', access='append')
  else
     open(11,file='eigval_de.out',status='new')
  endif
  write(11,*)"eigval de"
  close(11)


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
  write(11,*)"iter E dE E_ext E_h E_x eig_sum Z=",Z

  close(11)

  inquire(file='results_short.dat',EXIST=file_exists)
  if (file_exists) then
     open(11,file='results_short.dat',status='old', access='append')
  else
     open(11,file='results_short.dat',status='new')
     write(11,*)"Z E dE iterations Ngrid Rmax Rmin E_x"
  endif
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
vfull=0d0*r       
write(*,*)"********** version 2 *************"
! ----- version 2 -------
! Self-consistency loop
do iscl=1,maxscl
  psi_eig_temp=psi_eig
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
!to test Laplaciam subroutine alternative energy calculation for Helium
  call laplacian(Ngrid,r,psi(:,1),ftemp1)
  call integ_sph_s38_value(Ngrid,r,2d0*psi(:,1)*(-ftemp1/2d0-psi(:,1)*Z/r),e11)
  call integ_sph_s38_value(Ngrid,r,(vh/2d0+exc)*rho,e22)
  write(*,*)"Alternative Energy calculation:",e11,e22,e11+e22
!end Laplacian test
  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0
  close(11)


  open(11,file='eigval_de.out',status='old', access='append')
  do ish=1,Nshell
    write(11,*)iscl, psi_eig(ish), psi_eig(ish)-psi_eig_temp(ish), shell_occ(ish)
  enddo
  close(11)



  Print *,iscl,".iter E=",e1,"+",e2,"+",e3,"=",energy
  if (abs(energy-energy0).LT.dE_min) then
          EXIT
  endif
  vfull1=vh+vxc
  vfull=alpha*vfull1+(1d0-alpha)*vfull

enddo
  open(11,file='charge_potential.dat',status='replace')
  write(11,*)"r rho vc+vh+vxc"
   do ir = 1,Ngrid
     write(11,*)r(ir), rho(ir), vh(ir)+vxc(ir)+Z/r(ir)
   end do
   close(11)

  open(11,file='H_wf.dat',status='replace')
  write(11,*)"r psi"
   do ir = 1,Ngrid
     write(11,*)r(ir), psi(ir,1)
   end do
   close(11)

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
      call integ_sph_s38_value(Ngrid,r,psi(:,il_icl)**2,norm)
      norm_arr(il_icl)=norm
      psi(:,il_icl)=psi(:,il_icl)/dsqrt(norm)

#ifdef debug
  write(*,*)"PSI l=",il-1," icl=",icl," il_cl=",il_icl," norm:",norm,&
            "eig:",psi_eig(il_icl)," occ:",shell_occ(il_icl)
#endif
    enddo
  enddo
! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2
  enddo


  ! Construct the Hartree potential

  call integ_s38_fun(Ngrid,r,r**2*rho,1,ftemp1)
  call integ_s38_fun(Ngrid,r,r*rho,-1,ftemp2)
  vh=ftemp1/r+ftemp2

!   call getvxc(Ngrid,rho/(4d0*Pi),vxc,exc)
!  Caulculate vx_phi for every orbital  
  do ish=1,Nshell
   call get_Fock_ex(Ngrid,r,ish,Nshell,shell_l,psi_non_norm(:,ish),psi,vx_phi1(:,ish))
  enddo
  vx_phi=alpha*vx_phi1+(1d0-alpha)*vx_phi

  ! Caulculate energy
  call integ_sph_s38_value(Ngrid,r,-rho*Z/r,e_ext)
  write(*,*)"e_ext=",e_ext
  call integ_sph_s38_value(Ngrid,r,0.5d0*rho*vh,e_h)
  write(*,*)"e_h=",e_h
  e_x=0d0
  do ish=1,Nshell
    !call integ_sph_s38_value(Ngrid,r,vxc*psi(:,ish)**2d0,e2)
    call integ_sph_s38_value(Ngrid,r,0.5d0*shell_occ(ish)*psi(:,ish)*vx_phi(:,ish)/dsqrt(norm_arr(ish)),e2)
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
  energy0=energy

  write(*,*)"e_pot/e_kin=",e_pot/e_kin, "e_tot=",e_kin+e_pot, " (",iscl,")"


  !LDA LDA LDA energy
  call integ_sph_s38_value(Ngrid,r,-0.5d0*vh*rho,e2)
  call integ_sph_s38_value(Ngrid,r,(exc-vxc)*rho,e3)
  Print *,"Result in case of DFT vxc: ",iscl,".iter E=",e1,"+",e2,"+",e3,"=",e1+e2+e3
!  energy=e1+e2+e3

  !LDA LDA LDA energy

 energy=e_kin+e_pot

  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0,e_ext,e_h,e_x,e1
  close(11)

  open(11,file='eigval_de.out',status='old', access='append')
  do ish=1,Nshell
    write(11,*)iscl, psi_eig(ish), psi_eig(ish)-psi_eig_temp(ish), shell_occ(ish)
  enddo
  close(11)


  psi_eig_temp=psi_eig

  if (abs(energy-energy0).LT.dE_min) then
          EXIT
  endif

  vfull1=vh!+vxc
  vfull=alpha*vfull1+(1d0-alpha)*vfull
!  vfull=vfull1
enddo


  open(11,file='H_wf.dat',status='replace')
  write(11,*)"r psi"
   do ir = 1,Ngrid
     write(11,*)r(ir), psi(ir,1)
   end do
   close(11)


deallocate(vx_phi)

deallocate(vx_phi1,vx_psidot,psidot,psi_non_norm)

elseif(version.eq.4) then

write(*,*)"Version 4"
write(*,*)"lmax=",lmax



! ----- version x4 -------
! Lippmann–Schwinger iterations and solving screened Poisson equation.
! $\left( \nabla^2+ 2\epsilon \right \psi(\mathbf{r}) = v(\mathbf{r}) \psi(\mathbf{r})$

 
 l_n=0
  do il=1,lmax+1
    call diashoot(Ngrid,r,Z,vfull,il-1,count_l(il),l_n,Nshell,eig,psi)
    do inn=1,count_l(il)
      l_n=l_n+1
      psi_eig(l_n)=eig(inn)
      call integ_sph_s38_value(Ngrid,r,psi(:,l_n)**2d0,norm)
      psi(:,l_n)=psi(:,l_n)/dsqrt(norm)
      norm_arr(l_n)=dsqrt(norm)
      write(*,*)"l_n",l_n,"eig=",psi_eig(l_n),"norm=",norm
#ifdef debug
      write(*,*)"PSI l=",il-1," n=",inn," l_n=",l_n,&
            "eig:",psi_eig(l_n)," occ:",shell_occ(l_n)
#endif
    enddo
  enddo

Allocate(psip(Ngrid,Nshell),eigp(Nshell),vhp(Ngrid),vxcp(Ngrid))
vh=0d0*r
vxc=0d0*r
do iscl=1,20
! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2
  enddo

vhp=vh
vxcp=vxc
! Get exchange potential
  call getvxc(Ngrid,rho/(4d0*Pi),vxc,exc)
! Get Coulomb potential
  call integ_s38_fun(Ngrid,r,r**2*rho,1,ftemp1)
  call integ_s38_fun(Ngrid,r,r*rho,-1,ftemp2)
  vh=ftemp1/r+ftemp2

vh=vh*alpha+(1d0-alpha)*vhp
vxc=vxc*alpha+(1d0-alpha)*vxcp
  write(*,*)"EXTERNAL cycle begining:"
l_n=0
do il=1,lmax+1
  do inn=1,count_l(il)
    l_n=l_n+1
    write(*,*)"l=",il-1," n=",inn," eig=",psi_eig(l_n)
  enddo
enddo
l_n=0

  do il=1,lmax+1
    psip=psi
    eigp=psi_eig
    call LS_iteration(Ngrid,r,Z,il-1,count_l(il),l_n,Nshell,vxc,vh,psip,eigp,norm_arr,psi,psi_eig)
    do inn=1,count_l(il)
       l_n=l_n+1
    enddo

  enddo

  open(11,file='psi_new.dat',status='replace')
  write(11,*)"r psi1 psi2 psi3"
   do ir = 1,Ngrid 
!     write(11,*)r(ir), psi(ir,1),psi(ir,2),psi(ir,3)
   end do  
   close(11)

   enddo

write(*,*)"RESULT::"
l_n=0
do il=1,lmax+1
  do inn=1,count_l(il)
    l_n=l_n+1
    write(*,*)"l=",il-1," n=",inn," eig=",psi_eig(l_n)
  enddo
enddo

  Deallocate(psip,eigp,vxcp,vhp)



! Calculate \epsilon as an expectation value

!enddo

! ----- version x4 -------
! Self-consistency loop
!do iscl=1,maxscl

! Davidson method 

!enddo


endif

!  write results 
  open(11,file='results_short.dat',status='old', access='append')
!  write(11,*)"Z E dE iterations Ngrid Rmax Rmin E_x"
  write(11,*)Z, energy, energy-energy0, iscl-1, Ngrid, Rmax, Rmin, e_x
  close(11)


  inquire(file='results.dat',EXIST=file_exists)
  if (file_exists) then
     open(11,file='results.dat',status='old', access='append')
  else
     open(11,file='results.dat',status='new')
  endif
  write(11,*)"Z E dE iterations Ngrid Rmax Rmin"
  write(11,*)Z, version, iscl-1, Ngrid, Rmax
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


e1=0.8d0
i=0
j=0
do i=1,5
j=2
e1=dble(i)
call msbesselk (j, e1,besrez)
call besk_simp(j,e1,e3)

write(*,*)"x=",e1,"l=",j,"msbeselk=",besrez(j)," simple: ",e3

enddo


do i=1,5
j=2
e1=dble(i)
call msbesseli (j, e1,besrez)
call besi_simp(j,e1,e3)

write(*,*)"x=",e1,"l=",j,"msbeseli=",besrez(j)," simple: ",e3

enddo


end program

