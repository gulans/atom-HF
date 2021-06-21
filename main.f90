program atomHF


!!!!!!!!!! libxc !!!!!!!!!
  use xc_f03_lib_m
  implicit none
  TYPE(xc_f03_func_t) :: x_func,c_func,c1_func
  TYPE(xc_f03_func_info_t) :: x_info,c_info,c1_info
  integer :: vmajor, vminor, vmicro
  character(len=120) :: kind,family,version_text
  logical :: polarised
  real(8) :: hybx_w(4) !hyb exchange weigths  1-HF , 2-SR HF , 3-local 
  !  integer,external :: xc_family_from_id !!tā ir funkcija
  !!!!!!!!!! libxc !!!!!!!!
  !
real(16) :: r16





real(8), PARAMETER :: Pi = 3.1415926535897932384d0, alpha=0.5d0
real(8) :: Dnk,PHInk
!integer,parameter :: Ngrid = 500
real(8), allocatable :: r(:),vfull(:),vh(:),vxc(:),exc(:),H(:,:),eig(:),psi(:,:),rho(:),vfull1(:),ftemp1(:),&
        ftemp2(:),vx_psi(:,:),vx_phi(:,:),vx_phi1(:,:),vx_psidot(:,:),psidot(:,:),psi_non_norm(:,:),norm_arr(:),&
        vhp(:),vxcp(:),grho(:),grho2(:),g2rho(:),g3rho(:),vx_gga(:),vc_gga(:),ex_gga(:),ec_gga(:),vx(:),vc(:),&
        ex(:),ec(:),vxsigma(:),vcsigma(:),rho_pol(:,:),vx_pol(:,:),vc_pol(:,:),grho_pol(:,:),vsigma_pol(:,:),&
        vx_psi_sr(:,:),vx_psi_lr(:,:),vc1(:),ec1(:)
integer, allocatable :: shell_n(:),shell_l(:),count_l(:)
real(8), allocatable :: shell_occ(:),psi_eig(:),psi_eig_temp(:)
real(8), allocatable :: psip(:,:),eigp(:)
complex(8), allocatable :: Bess_ik(:,:,:,:)

!real(8), parameter :: Rmax = 10d0
integer, parameter :: maxscl =50 !Maximal itteration number
integer :: il,icl,il_icl,iscl,lmax,l_n,inn
real(8) :: Z,rez,a1,a2
real(8) :: norm,Rmin,Rmax,hh,dE_min,e1,e2,e3,energy,energy0,e_kin,e_ext,e_h,e_x,e_pot
real(8) :: e11,e22,e33
integer :: grid, Nshell, ish, d_order,i_order
integer :: ir,i,j,countl0, version,tools_info(3)
integer(8) :: Ngrid 
logical :: E_dE_file_exists, file_exists
character(len=1024) :: filename

Real (8)  :: besrez (0:50),time0,time1, t11
complex(8) :: besrezc (0:50)
real(8), allocatable :: tools(:,:)
integer :: x_num,c_num,c1_num

complex(8) :: ac,bc,cc
real(8) :: rsmu
integer :: rs,Nrsfun
complex(8), allocatable :: rsfunC(:,:)



dE_min=1d-8
call timesec(time0)

!read input
read(*,*) 
read(*,*) Z, Rmin, Rmax, Ngrid, version
read(*,*)
read(*,*) x_num, c_num,c1_num
!c1_num=0
read(*,*)
read(*,*) grid
read(*,*) 
read(*,*) Nshell
allocate(shell_n(Nshell),shell_l(Nshell),shell_occ(Nshell),norm_arr(Nshell))
read(*,*)
do ish=1,Nshell
  read(*,*) shell_n(ish), shell_l(ish), shell_occ(ish)
enddo


!!!!!!!!!! libxc !!!!!!!!!
if ((version.eq.5).or.(version.eq.8)) then 
! Print out the version
          call xc_f03_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
!  x_num=1
!  c_num=12
!  c_num=7
!  x_num=101
!  c_num=130

if (x_num.ne.0) then
!!!info Exchange!!!!!!!
write(*,*)"x_num",x_num

  call xc_f03_func_init(x_func, x_num, XC_UNPOLARIZED)

  x_info = xc_f03_func_get_info(x_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(x_info))
  case (XC_EXCHANGE)
    write(kind, '(a)') 'an exchange functional'
  case (XC_CORRELATION)
    write(kind, '(a)') 'a correlation functional'
  case (XC_EXCHANGE_CORRELATION)
    write(kind, '(a)') 'an exchange-correlation functional'
  case (XC_KINETIC)
    write(kind, '(a)') 'a kinetic energy functional'
  case default
    write(kind, '(a)') 'of unknown kind'
  end select
  ! Get the family
write(*,*)"FAMILY:"
  write(*,*)xc_f03_func_info_get_family(x_info)
write(*,*)XC_FAMILY_LDA
  select case (xc_f03_func_info_get_family(x_info))
  case (XC_FAMILY_LDA);
    write(family,'(a)') "LDA"
  case (XC_FAMILY_GGA);
    write(family,'(a)') "GGA"
  case (XC_FAMILY_HYB_GGA);
    write(family,'(a)') "Hybrid GGA"
  case (XC_FAMILY_MGGA);
    write(family,'(a)') "MGGA"
  case (XC_FAMILY_HYB_MGGA);
    write(family,'(a)') "Hybrid MGGA"
  case default;
    write(family,'(a)') "unknown"
  end select
  ! Print out information
  write(*,'("The functional ''", a, "'' is ", a, ", it belongs to the ''", a, "'' family and is defined in the reference(s):")') &
    trim(xc_f03_func_info_get_name(x_info)), trim(kind), trim(family)
  ! Print out references
  i = 0
  !do while(i >= 0)
  !  write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(x_info, i)))
  !end do
  endif
if (c_num.ne.0) then

!!!!info correlation
write(*,*)"c_num",c_num
  call xc_f03_func_init(c_func, c_num, XC_UNPOLARIZED)

  c_info = xc_f03_func_get_info(c_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(c_info))
  case (XC_EXCHANGE)
    write(kind, '(a)') 'an exchange functional'
  case (XC_CORRELATION)
    write(kind, '(a)') 'a correlation functional'
  case (XC_EXCHANGE_CORRELATION)
    write(kind, '(a)') 'an exchange-correlation functional'
  case (XC_KINETIC)
    write(kind, '(a)') 'a kinetic energy functional'
  case default
    write(kind, '(a)') 'of unknown kind'
  end select
  ! Get the family
  select case (xc_f03_func_info_get_family(c_info))
  case (XC_FAMILY_LDA);
    write(family,'(a)') "LDA"
  case (XC_FAMILY_GGA);
    write(family,'(a)') "GGA"
  case (XC_FAMILY_HYB_GGA);
    write(family,'(a)') "Hybrid GGA"
  case (XC_FAMILY_MGGA);
    write(family,'(a)') "MGGA"
  case (XC_FAMILY_HYB_MGGA);
    write(family,'(a)') "Hybrid MGGA"
  case default;
    write(family,'(a)') "unknown"
  end select
  ! Print out information
  write(*,'("The functional ''", a, "'' is ", a, ", it belongs to the ''", a, "'' family and is defined in the reference(s):")') &
    trim(xc_f03_func_info_get_name(c_info)), trim(kind), trim(family)
  ! Print out references
  i = 0
 ! do while(i >= 0)
 !   write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(c_info, i)))
 ! end do
  endif

!!!!info correlation1
if (c1_num.ne.0) then
write(*,*)"c1_num",c1_num

  call xc_f03_func_init(c1_func, c1_num, XC_UNPOLARIZED)

  c1_info = xc_f03_func_get_info(c1_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(c1_info))
  case (XC_EXCHANGE)
    write(kind, '(a)') 'an exchange functional'
  case (XC_CORRELATION)
    write(kind, '(a)') 'a correlation functional'
  case (XC_EXCHANGE_CORRELATION)
    write(kind, '(a)') 'an exchange-correlation functional'
  case (XC_KINETIC)
    write(kind, '(a)') 'a kinetic energy functional'
  case default
    write(kind, '(a)') 'of unknown kind'
  end select
  ! Get the family
  select case (xc_f03_func_info_get_family(c1_info))
  case (XC_FAMILY_LDA);
    write(family,'(a)') "LDA"
  case (XC_FAMILY_GGA);
    write(family,'(a)') "GGA"
  case (XC_FAMILY_HYB_GGA);
    write(family,'(a)') "Hybrid GGA"
  case (XC_FAMILY_MGGA);
    write(family,'(a)') "MGGA"
  case (XC_FAMILY_HYB_MGGA);
    write(family,'(a)') "Hybrid MGGA"
  case default;
    write(family,'(a)') "unknown"
  end select
  ! Print out information
  write(*,'("The functional ''", a, "'' is ", a, ", it belongs to the ''", a, "'' family and is defined in the reference(s):")') &
    trim(xc_f03_func_info_get_name(c1_info)), trim(kind), trim(family)
  ! Print out references
  i = 0
 ! do while(i >= 0)
 !   write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(c1_info, i)))
 ! end do

endif
  
  !!!!!!!!!! libxc !!!!!!!!!
endif



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
allocate(psi(Ngrid,Nshell),psi_eig(Nshell),psi_eig_temp(Nshell),grho2(Ngrid))
allocate(ftemp1(Ngrid),ftemp2(Ngrid),vx(Ngrid),vc(Ngrid),ex(Ngrid),ec(Ngrid),vxsigma(Ngrid),vcsigma(Ngrid),grho(Ngrid))
allocate(rho_pol(2,Ngrid),vx_pol(2,Ngrid),vc_pol(2,Ngrid),grho_pol(3,Ngrid),vsigma_pol(3,Ngrid))
allocate(g2rho(Ngrid),ec1(Ngrid),vc1(Ngrid))
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

  inquire(file='res.dat',EXIST=file_exists)
  if (file_exists) then
     open(11,file='res.dat',status='old', access='append')
  else
     open(11,file='res.dat',status='new')
     write(11,*)"v-text,version,x_num,c_num,d_order,i_order,Ngrid,Rmin,Rmax,Z,energy,time,Nrsfun,rs_mu"
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

d_order=8
i_order=7
        tools_info=(/40,d_order,i_order/)!optimal 8,7
 !info about tools_info array:
 !tools_info(1) - (2-nd dimenstion size of tools array 1-10 - coefficinets for Lagrange interp,
 !                11-20 coefficints for (ri+1/4) interpolation, and 21-30 for (ri+2/4), 31-40 for (ri+3/4) 
 !tools_info(2) - order of intrpolation polynom for Lagrange interpolation when calculating derivative,
 !tools_info(3) - order for inperolation for integration with Bodes folmula)

Allocate(tools(Ngrid,tools_info(1)))
call generate_tools(Ngrid,r,tools,tools_info)




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
!  call getxc_pbe(Ngrid,r,tools,tools_info,rho/(4d0*Pi),vxc,exc)

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
   call get_Fock_ex(Ngrid,r,ish,Nshell,shell_l,psi_non_norm(:,ish),psi,vx_phi1(:,ish),rs,rsfunC,Nrsfun)
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

elseif((version.eq.4).or.(version.eq.5).or.(version.eq.6).or.(version.eq.7).or.(version.eq.8)) then
!4-Fock exchange
!5 - LDA exchange correlation
!6 - GGA-PBE
!7 - PBE0 

Nrsfun=8
allocate(rsfunC(Nrsfun,2))
allocate(Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)) !last atgument 1- 1-st kind (msbesi), 2- 2-nd kind (msbesk) 

!d_order=4
!i_order=3


d_order=8
i_order=7
        tools_info=(/40,d_order,i_order/)!optimal 8,7
 !info about tools_info array:
 !tools_info(1) - (2-nd dimenstion size of tools array 1-10 - coefficinets for Lagrange interp,
 !                11-20 coefficints for (ri+1/4) interpolation, and 21-30 for (ri+2/4), 31-40 for (ri+3/4) 
 !tools_info(2) - order of intrpolation polynom for Lagrange interpolation when calculating derivative,
 !tools_info(3) - order for inperolation for integration with Bodes folmula)

Allocate(tools(Ngrid,tools_info(1)))
call generate_tools(Ngrid,r,tools,tools_info)

if(version.eq.4)then
write(*,*)"Version 4 - Fock"
elseif  (version.eq.5) then
write(*,*)"Version 5 - libxc functionals"
elseif (version.eq.6) then
        if (x_num.eq.1) then
         write(*,*)"Version 6 - Bult-in functional LDA"
        elseif (x_num.eq.2) then
         write(*,*)"Version 6 - Bult-in functional GGA-PBE"
        endif
elseif (version.eq.7) then
write(*,*)"Version 7 - PBE0"
elseif (version.eq.8) then
write(*,*)"Hybrid libxc test"



endif
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
!      write(*,*)"PSI l=",il-1," n=",inn," l_n=",l_n,&
!            "eig:",psi_eig(l_n)," occ:",shell_occ(l_n)
#endif
    enddo
  enddo

Allocate(psip(Ngrid,Nshell),vx_psi(Ngrid,Nshell),vx_psi_sr(Ngrid,Nshell),vx_psi_lr(Ngrid,Nshell),&
        eigp(Nshell),vhp(Ngrid),vxcp(Ngrid))


vh=0d0*r
vxc=0d0*r

if(.true.)then
 rs=0
 rsmu=0.11d0
 if (rs.eq.0)then
         rsmu=0d0
 endif
 call errfun(Ngrid,r,Nrsfun,rsmu,rsfunC)
 call get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
 hybx_w(1)=0d0 !Coulumb interaction weigth in Fock exchange
 hybx_w(2)=0d0 !Short range interaction weigth in Fock exchange
 hybx_w(3)=rsmu !short range parameter
endif

!Some initialisation for libxc Hybrid exchange
if ((version.eq.5).or.(version.eq.8)) then
  if(x_num.ne.0)then
    if  (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_HYB_GGA) then
    write(*,*)"Hybrid Exchange!"
    e1= xc_f03_hyb_exx_coef(x_func) !hyb_exx - exact-exchange - Foka apmaiņas svars
    write(*,*) "Foka apmaiņas svars: ",e1
    call xc_f03_hyb_cam_coef(x_func , rsmu, a1, a2);
    write(*,*) "rs parameter mu:",rsmu
    write(*,*) "rs alpha",a1
    write(*,*) "rs beta",a2
    hybx_w(1)=a1
    hybx_w(2)=a2
    hybx_w(3)=rsmu

    if (rsmu.gt.1d-20) then
         rs=1
         !set array rsfunC - array of complex coeficients for erfc expantion
         call errfun(Ngrid,r,Nrsfun,rsmu,rsfunC)
         !obrain Bessel function Real part for each erfc expantion element 
         call get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
    else
         rs=0
    endif
  endif
  endif
endif



do iscl=1,100
      !open(12,file='r.dat',status='replace')
      !   do i=1, Ngrid
      !write(12,*) i,r(i)
      !enddo
      !close(12)
      !stop


! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2
  enddo

if ((version.eq.5).or.(version.eq.8)) then
  vxcp=vxc
!  call getvxc(Ngrid,rho/(4d0*Pi),vxc,exc)

rho=rho/(4d0*Pi)
vx=0d0*r
vc=0d0*r
ex=0d0*r
ec=0d0*r
vc1=0d0*r
ec1=0d0*r
!! EXCHANGE  !!!

if (x_num.ne.0)then

if (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_LDA) then
        
  call xc_f03_lda_exc_vxc(x_func, Ngrid, rho(1), ex(1),vx(1))

elseif  (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_GGA) then
write(*,*)"GGA"
if(x_num.eq.524)then
 hybx_w(3)=0.11d0
 write(*,*)"RS sparameter for local exchange:", hybx_w(3)
 call xc_f03_func_set_ext_params(x_func,hybx_w(3))

! call xc_f03_func_set_ext_params(x_func,hybx_w(3))
  endif
  call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
  g2rho=g2rho*r**(-2)
  grho2=grho**2
  call xc_f03_gga_exc_vxc(x_func, Ngrid, rho(1), grho2(1), ex(1),vx(1),vxsigma(1))

  !!!!! formula 6.0.8 !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vxsigma,ftemp1)
!  ftemp1=ftemp1*r**(-2)
!  vx=vx-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!


!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
   call rderivative_lagrN(Ngrid,r,tools,tools_info,vxsigma,ftemp1)
   vx=vx-2d0*(grho*ftemp1+vxsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!

elseif  (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_HYB_GGA) then

  call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
  g2rho=g2rho*r**(-2)
  grho2=grho**2

  call xc_f03_gga_exc_vxc(x_func, Ngrid, rho(1), grho2(1), ex(1),vx(1),vxsigma(1))

!  !!!!! formula 6.0.8 !!!!!!
  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vxsigma,ftemp1)
  ftemp1=ftemp1*r**(-2)
  vx=vx-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!


!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,vxsigma,ftemp1)
!  vx=vx-2d0*(grho*ftemp1+vxsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!

else

   write(*,*)"Exchange not supported!"
endif
endif
!! END EXCHANGE  !!!


!! CORRELATION !!!
if (c_num.ne.0) then
if (xc_f03_func_info_get_family(c_info).eq.XC_FAMILY_LDA) then

  call xc_f03_lda_exc_vxc(c_func, Ngrid, rho(1), ec(1),vc(1))

elseif  (xc_f03_func_info_get_family(c_info).eq.XC_FAMILY_GGA) then
  call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
  g2rho=g2rho*r**(-2)
  grho2=grho**2
  call xc_f03_gga_exc_vxc(c_func, Ngrid, rho(1), grho2(1), ec(1),vc(1),vcsigma(1))

 !!!!! formula 6.0.8 !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vcsigma,ftemp1)
!  ftemp1=ftemp1*r**(-2)
!  vc=vc-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!
  
  
!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
   call rderivative_lagrN(Ngrid,r,tools,tools_info,vcsigma,ftemp1)
   vc=vc-2d0*(grho*ftemp1+vcsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!

else
   write(*,*)"Correlation not supported!"
endif
endif
!! END CORRELATION !!

!!!!!! PAPILDUS COREL
if (c1_num.ne.0) then
if (xc_f03_func_info_get_family(c1_info).eq.XC_FAMILY_LDA) then

  call xc_f03_lda_exc_vxc(c1_func, Ngrid, rho(1), ec1(1),vc1(1))

elseif  (xc_f03_func_info_get_family(c1_info).eq.XC_FAMILY_GGA) then
  call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
  g2rho=g2rho*r**(-2)
  grho2=grho**2
  call xc_f03_gga_exc_vxc(c1_func, Ngrid, rho(1), grho2(1), ec1(1),vc1(1),vcsigma(1))

 !!!!! formula 6.0.8 !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vcsigma,ftemp1)
!  ftemp1=ftemp1*r**(-2)
!  vc1=vc1-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!
  
!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
   call rderivative_lagrN(Ngrid,r,tools,tools_info,vcsigma,ftemp1)
   vc1=vc1-2d0*(grho*ftemp1+vcsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!

else
   write(*,*)"Correlation not supported!"
endif
endif

     ! open(12,file='vx_all_exc.dat',status='replace')
     ! write(12,*) 'r -0.25d0*vx(i) , vc(i) , vc1(i) , -0.25d0*vx(i)+vc1(i)+vc(i)'
     ! do i=1, Ngrid
     ! write(12,*) r(i),',',-0.25d0*vx(i),",",vc(i),",",vc1(i),",",-0.25d0*vx(i)+vc1(i)+vc(i)
     ! enddo
     ! close(12)
!stop
!vx=-0.25d0*vx
!ex=-0.25d0*ex

vc=vc!+vc1
ec=ec!+ec1
!!!!!! End Papildus COREL

vxc=vx+vc


exc=ex+ec
rho=rho*(4d0*Pi)

else if (version.eq.6) then
        c_num=0
        if (x_num.eq.1) then
        vxcp=vxc
        vx=0d0*r
        vc=0d0*r
        ex=0d0*r
        ec=0d0*r


        call getvxc(Ngrid,rho/(4d0*Pi),vxc,exc)

        elseif (x_num.eq.2) then
        vxcp=vxc
        call getxc_pbe(Ngrid,r,tools,tools_info,rho/(4d0*Pi),vxc,exc,vx,vc,ex,ec)

        endif
else if (version.eq.7) then
  vxcp=vxc
  call getxc_pbe(Ngrid,r,tools,tools_info,rho/(4d0*Pi),vxc,exc,vx,vc,ex,ec)

endif

! Get exchange potential

if ((version.eq.4).or.(version.eq.7).or.(version.eq.8)) then
 do ish=1,Nshell
   call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,lmax,psi(:,ish),psi,&
           vx_psi(:,ish),vx_psi_sr(:,ish),vx_psi_lr(:,ish),rs,rsfunC,Nrsfun,hybx_w,Bess_ik)
  enddo 

 !open(11,file='Fock_vx.dat',status='replace')
 !write(11,*)
 !  do i=1, Ngrid
 !    write(11,*)r(i), vx_psi(i,1)!,vx_psi(i,2)!,vx_psi(i,3)
 !  enddo
 !  close(11)
 !open(11,file='Fock_vx_rs.dat',status='replace')
 !write(11,*)
 !  do i=1, Ngrid
 !    write(11,*)r(i), vx_psi_sr(i,1)!,vx_psi_sr(i,2)!,vx_psi_sr(i,3)
 !  enddo
 !  close(11)
!stop

endif
  
! Get Coulomb potential

call  integ_BodesN_fun(Ngrid,r,tools,tools_info,1,r**2*rho,ftemp1)
call  integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,r*rho,ftemp2)
vhp=vh
vh=ftemp1/r+ftemp2
vh=vh*alpha+(1d0-alpha)*vhp

!write(*,*)"EXTERNAL cycle begining:"

l_n=0
do il=1,lmax+1
  do inn=1,count_l(il)
    l_n=l_n+1
!    write(*,*)"l=",il-1," n=",inn," eig=",psi_eig(l_n)
  enddo
enddo


!Check convergence

!write(*,*)"External psi_eig-eigp: ",psi_eig-eigp
if((maxval(abs((psi_eig-eigp)/(psi_eig+1d0)))).lt.1d-13)then !1d-13 bija
!        write(*,*)"Arējais cikls konverģē:"
!        write(*,*)"konverģence absolūtā: ",maxval(abs(psi_eig-eigp))," relatīvā: ",maxval(abs((psi_eig-eigp)/(psi_eig+1d0)))

        exit
endif




! Caulculate energy

  e1=0d0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
if (version.eq.4)then
  call integ_BodesN_value(Ngrid,r,tools,tools_info,-r*rho*Z,e_ext)
  write(*,*)"e_ext=",e_ext
  call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*rho*vh,e_h)
  write(*,*)"e_h=",e_h
  e_x=0d0
  do ish=1,Nshell
    call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*shell_occ(ish)*psi(:,ish)*vx_psi(:,ish),e2)
    e_x=e_x+e2
  
  enddo
e3=0d0
  write(*,*)"e_x=",e_x
  e_kin=e1-e_ext-2d0*e_h-2d0*e_x
  write(*,*)"e_kin=",e_kin
  e_pot=e_ext+e_h+e_x
  energy0=energy
  write(*,*)"e_pot/e_kin=",e_pot/e_kin, "e_tot=",e_kin+e_pot, " (",iscl,")"
  energy=e_kin+e_pot


else if((version.eq.5).or.(version.eq.6)) then
  !LDA LDA LDA energy

  call integ_BodesN_value(Ngrid,r,tools, tools_info,-0.5d0*vh*rho*r**2,e2)
  call integ_BodesN_value(Ngrid,r,tools, tools_info,(exc-vxc)*rho*r**2,e3)
  Print *,iscl,".iter E=",e1,"+",e2,"+",e3,"=",e1+e2+e3
  energy0=energy
  energy=e1+e2+e3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

else if(version.eq.7)then
  e_kin=0d0
  do ish=1,Nshell

  call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,ish),ftemp1)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,ftemp1*r**2,ftemp2)
  ftemp2=ftemp2/r**2
  call integ_BodesN_value(Ngrid,r,tools, tools_info,-0.5d0*psi(:,ish)*ftemp2*r**2,e1)
  call integ_BodesN_value(Ngrid,r,tools, tools_info,0.5d0*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish)**2,e2)
  e_kin=e_kin+shell_occ(ish)*(e1+e2)
  enddo
  !write(*,*)"e_kin_new:",e_kin
  call integ_BodesN_value(Ngrid,r,tools,tools_info,-r*rho*Z,e_ext)
  !write(*,*)"e_ext=",e_ext
  call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*rho*vh,e_h)
  !write(*,*)"e_h=",e_h
  e2=0d0
  do ish=1,Nshell
    call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*shell_occ(ish)*psi(:,ish)*0.25d0*vx_psi(:,ish),e1)
    e2=e2+e1
  enddo
  call integ_BodesN_value(Ngrid,r,tools, tools_info,(0.75d0*ex+ec)*rho*r**2,e3)

 e_x=e2+e3
  energy0=energy
  energy=e_kin+e_ext+e_h+e_x

  write(*,*)"Energy_tot:",energy, " (",iscl,")"

else if(version.eq.8)then
  e_kin=0d0
  do ish=1,Nshell

  call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,ish),ftemp1)

!!!šī nedaudz neprecīāka metode
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,ftemp1*r**2,ftemp2)
!  ftemp2=ftemp2/r**2
!  call integ_BodesN_value(Ngrid,r,tools, tools_info,-0.5d0*psi(:,ish)*ftemp2*r**2,e1)
!!ši precīzāka
  call rderivative_lagrN(Ngrid,r,tools,tools_info,ftemp1,ftemp2)
  call integ_BodesN_value(Ngrid,r,tools, tools_info,-0.5d0*psi(:,ish)*(2d0*ftemp1*r+ftemp2*r**2),e1)



call integ_BodesN_value(Ngrid,r,tools, tools_info,0.5d0*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish)**2,e2)
  e_kin=e_kin+shell_occ(ish)*(e1+e2)
  enddo
  write(*,*)"e_kin_new:",e_kin
  call integ_BodesN_value(Ngrid,r,tools,tools_info,-r*rho*Z,e_ext)
  write(*,*)"e_ext=",e_ext
  call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*rho*vh,e_h)
  write(*,*)"e_h=",e_h
  e2=0d0
  do ish=1,Nshell
    call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*shell_occ(ish)*psi(:,ish)*&
            (hybx_w(1)*vx_psi(:,ish)+hybx_w(2)*vx_psi_sr(:,ish)),e1)



    e2=e2+e1
  enddo
  call integ_BodesN_value(Ngrid,r,tools, tools_info,(ex+ec)*rho*r**2,e3)

 e_x=e2+e3
  write(*,*)"e_xc=",e_x
!!!!Alternative kinetic energy:
  e1=0d0
  do ish=1,Nshell
    e1=e1+shell_occ(ish)*psi_eig(ish)
  enddo
e2=e1-e_ext-2d0*e_h-2d0*e_x 



write(*,*)"e kin _alternative=",e2
write(*,*)"e tot _alternative=",e2+e_ext+e_h+e_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 energy0=energy
  energy=e_kin+e_ext+e_h+e_x

  write(*,*)"Energy_tot:",energy, " (",iscl,")"


endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!END LDA LDA LDA energy

  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0,e_ext,e_h,e_x,e1
  close(11)

  open(11,file='eigval_de.out',status='old', access='append')
  do ish=1,Nshell
    write(11,*)iscl, psi_eig(ish), psi_eig(ish)-eigp(ish), shell_occ(ish)
  enddo
  close(11)

!END Caulculate energy





l_n=0
 psip=psi
 eigp=psi_eig

  do il=1,lmax+1
 call LS_iteration(Ngrid,r,tools,tools_info,rs,rsfunC,Nrsfun,hybx_w,version,Z,il-1,shell_l,count_l(il),l_n,&
            Nshell,lmax,vxc,vx,vc,vh,vx_psi,vx_psi_sr,vx_psi_lr,psip,norm_arr,psi,psi_eig,Bess_ik)

!   open(11,file='Fock_vxpsi.dat',status='replace')
!   write(11,*)"r(i), vx_psi(i), vx_psi_sr(i),vx_psi_lr(i)",rsmu
!   do i=1, Ngrid
!     write(11,*)r(i), vx_psi(i,2), vx_psi_sr(i,2),vx_psi_lr(i,2)
!   enddo
!   close(11)
!   stop


!write(*,*)"psi_eig: ",psi_eig
!write(*,*)"eigp: ",eigp


    do inn=1,count_l(il)
       l_n=l_n+1
    enddo

  enddo

!  if (iscl.eq.21)then
!  open(11,file='psi_21.dat',status='replace')
!   do i=1, Ngrid
!     write(11,*)r(i), psi(i,1), psi(i,2),psi(i,3),vx_psi(i,1),vx_psi(i,2),vx_psi(i,3)
!   enddo
!   close(11)
!   stop
!   endif




enddo
!End of self consistent loop

call timesec(time1)

if(.false.)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LIBXC test 101 & 524(w=0)


  call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
  g2rho=g2rho*r**(-2)
  grho2=grho**2
!  call xc_f03_gga_exc_vxc(c1_func, Ngrid, rho(1), grho2(1), ec1(1),vc1(1),vcsigma(1))

 !!!!! formula 6.0.8 !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vcsigma,ftemp1)
!  ftemp1=ftemp1*r**(-2)
!  vc1=vc1-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!

!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
!   call rderivative_lagrN(Ngrid,r,tools,tools_info,vcsigma,ftemp1)
!   vc1=vc1-2d0*(grho*ftemp1+vcsigma*g2rho)

open(11,file='30-punkti.dat',status='replace')
write(11,*)"r=(/"
do i=1,Ngrid,4
write(11,*)r(i),",",r(i+1),",",r(i+2),",",r(i+3)
enddo
write(11,*)"/)"



write(11,*)"rho=(/"
do i=1,Ngrid,4
write(11,*)rho(i),",",rho(i+1),",",rho(i+2),",",rho(i+3)
enddo
write(11,*)"/)"

write(11,*)"grho=(/"
do i=1,Ngrid,4
write(11,*)grho2(i),",",grho2(i+1),",",grho2(i+2),",",grho2(i+3)
enddo
write(11,*)"/)"

close(11)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!LIBXC test 101 & 524(w=0)
endif !!LIBXC test


















write(*,*)"RESULT::"
l_n=0
do il=1,lmax+1
  do inn=1,count_l(il)
    l_n=l_n+1
    write(*,*)"l=",il-1," n=",inn," eig=",psi_eig(l_n)
  enddo
enddo

  Deallocate(vx_psi,psip,eigp,vxcp,vhp,tools)




! Calculate \epsilon as an expectation value

!enddo

! ----- version x4 -------
! Self-consistency loop
!do iscl=1,maxscl

! Davidson method 

!enddo

else if (version.eq.10)then


open(11,file='bestesti_0_small.dat',status='replace')
do i=1,Ngrid
ac=cmplx(r(i),r(i),8)
call msbesselic(0,ac, bc)
write(11,*)realpart(ac),imagpart(ac),realpart(bc),imagpart(bc)
enddo
close(11)

open(11,file='bestestk_0_small.dat',status='replace')
do i=1,Ngrid
ac=cmplx(r(i),r(i),8)
call msbesselkc(0,ac, bc)
write(11,*)realpart(ac),imagpart(ac),realpart(bc),imagpart(bc)
enddo
close(11)



open(11,file='bestesti_1_small.dat',status='replace')
do i=1,Ngrid
ac=cmplx(r(i),r(i),8)
call msbesselic(1,ac, bc)
write(11,*)realpart(ac),imagpart(ac),realpart(bc),imagpart(bc)
enddo
close(11)

open(11,file='bestestk_1_small.dat',status='replace')
do i=1,Ngrid
ac=cmplx(r(i),r(i),8)
call msbesselkc(1,ac, bc)
write(11,*)realpart(ac),imagpart(ac),realpart(bc),imagpart(bc)
enddo
close(11)



open(11,file='bestesti_2_small.dat',status='replace')
do i=1,Ngrid
ac=cmplx(r(i),r(i),8)
call msbesselic(2,ac, bc)
write(11,*)realpart(ac),imagpart(ac),realpart(bc),imagpart(bc)
enddo
close(11)

open(11,file='bestestk_2_small.dat',status='replace')
do i=1,Ngrid
ac=cmplx(r(i),r(i),8)
call msbesselkc(2,ac, bc)
write(11,*)realpart(ac),imagpart(ac),realpart(bc),imagpart(bc)

write(*,*)ac,exp(-ac),(ac+1d0)/ac**2

enddo
close(11)


ac=cmplx(1d0,1d0,8)
call msbesselic(8,ac, bc)
write(*,*)8,ac,bc
call msbesselic(9,ac, bc)
write(*,*)9,ac,bc
call msbesselic(10,ac, bc)
write(*,*)10,ac,bc

ac=cmplx(4.5617402703990642E-010,   1.0160216352396843E-010,8)
call msbesselic(3,ac, bc)
write(*,*)3,ac,bc



do i=-10,1

ac=cmplx(10d0**i,10d0**i,8)
write(*,*)"argum",ac
call msbesselic(3,ac, bc)
write(*,*)"besi",bc
call msbesselkc(3,ac, bc)
write(*,*)"besk",bc
write(*,*)""
enddo



endif



if (version.eq.4) then
version_text='HF'
elseif (version.eq.5) then
  version_text='Libxc'
elseif (version.eq.6) then
        if (x_num.eq.1) then
        version_text='LDA'
        elseif (x_num.eq.2) then
        version_text='PBE'
        endif
elseif (version.eq.7) then
version_text='PBE0'
else
        version_text='unknown'
endif

!  write results 
  open(11,file='res.dat',status='old', access='append')
!  write(11,*)"version Z E dE iterations Ngrid Rmax Rmin E_x, e_pot/e_kin, e_homo"
!  write(11,*)Z, energy, energy-energy0, iscl-1, Ngrid, Rmax, Rmin, e_x, e_pot/e_kin, maxval(psi_eig) 
  !write(11,*)version,",",Z,",",d_order,",",i_order,",",Ngrid,",", Rmin,",", Rmax,",", energy,",", time1-time0
  write(11, '(a8,a1)',advance="no")trim(version_text),","
  write(11, '(i3,a1,i3,a1,i3,a1,i3,a1,i3,a1,i5,a1)',advance="no") version,",",x_num,",",c_num,","&
          ,d_order,",",i_order,",",Ngrid,","
  write(11, '(ES9.2E2,a1)',advance="no")Rmin,","
 write(11, '(f5.2)' ,advance="no") Rmax
 write(11, '(a1)',advance="no")","
  write(11, '(f5.2)',advance="no")Z
  write(11, *)",", energy,",", time1-time0,",",Nrsfun,",",hybx_w(3)
  
  close(11)


  inquire(file='output.dat',EXIST=file_exists)
  if (file_exists) then
     open(11,file='output.dat',status='old', access='append')
  else
     open(11,file='output.dat',status='new')
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


end program

