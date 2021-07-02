program atomHF


!!!!!!!!!! libxc !!!!!!!!!
  use xc_f03_lib_m
  implicit none
  TYPE(xc_f03_func_t) :: x_func,c_func,c1_func
  TYPE(xc_f03_func_info_t) :: x_info,c_info,c1_info
  integer :: vmajor, vminor, vmicro
  character(len=120) :: kind,family,version_text
  logical :: polarised
  real(8) :: hybx_w(5,2) !hyb exchange weigths: (1,1)-vx, (2,1)-vc, (3,1)-vc1, (4,1)-HF, (5,1)-HF_sr
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
integer :: Nrsfun
complex(8), allocatable :: rsfunC(:,:)

logical :: override_libxc_hyb
integer :: param_nr1,param_nr2,param_nr3
real(8), allocatable :: param1(:), param2(:), param3(:)


dE_min=1d-8
call timesec(time0)

!read input
read(*,*) 
read(*,*) Z, Rmin, Rmax, Ngrid, version
read(*,*)
read(*,*) x_num, hybx_w(1,1), param_nr1
allocate(param1(param_nr1))
do i=1, param_nr1
read(*,*) param1(i)
enddo

read(*,*) c_num, hybx_w(2,1), param_nr2
allocate(param2(param_nr2))
do i=1, param_nr2
read(*,*) param2(i)
enddo

read(*,*) c1_num,hybx_w(3,1), param_nr3
allocate(param3(param_nr3))
do i=1, param_nr3
read(*,*) param3(i)
enddo


read(*,*) 
read(*,*) override_libxc_hyb
read(*,*) hybx_w(4,1)
read(*,*) hybx_w(5,1),hybx_w(5,2)
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
if ((version.eq.0).or.(version.eq.1)) then 
! Print out the version
          call xc_f03_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
!  x_num=1
!  c_num=12
!  c_num=7
!  x_num=101
!  c_num=130

if (x_num.gt.0) then
!!!info Exchange!!!!!!!
write(*,*)
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

  i = 0
  if(x_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(x_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(x_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(x_info),  "external parameters."
   
  do i=0, xc_f03_func_info_get_n_ext_params(x_info)-1
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(x_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(x_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(x_info,i))
 enddo
  endif

 write(*,*)

  
if (c_num.gt.0) then

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
  if(c_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(c_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(c_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(c_info),  "external parameters."

  do i=0, xc_f03_func_info_get_n_ext_params(c_info)-1
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(c_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(c_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(c_info,i))
  enddo

  endif
 write(*,*)

!!!!info correlation1
if (c1_num.gt.0) then
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
  if(c1_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(c1_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(c1_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(c1_info),  "external parameters."

  do i=0, xc_f03_func_info_get_n_ext_params(c1_info)-1
  
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(c1_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(c1_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(c1_info,i))
 enddo


endif
 write(*,*)


 if (param_nr1.ne.0)then
   if (xc_f03_func_info_get_n_ext_params(x_info).ne.param_nr1) then
         write(*,*)" Parameter nr for Functional ",x_num, "should be", xc_f03_func_info_get_n_ext_params(x_info),&
                 ", not ",param_nr1
         stop
   else
    call xc_f03_func_set_ext_params(x_func,param1(1))
   endif     
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



  inquire(file='eigval_de.out',EXIST=E_dE_file_exists)
  if (E_dE_file_exists) then
     open(11,file='eigval_de.out',status='old', access='append')
  else
     open(11,file='eigval_de.out',status='new')
  endif
  write(11,*)"eigval de"
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
     write(11,*)"grid,iter,x_num,c_num,c2_num,d_order,i_order,Ngrid,Rmin,Rmax,Z,time,Nrsfuni&
             ,Fock_w,Fock_rs_w,rs_mu,",&
             "dE,energy"
  endif
  close(11)


if (version.eq.1) then

Nrsfun=8
allocate(rsfunC(Nrsfun,2))
allocate(Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)) !last atgument 1- 1-st kind (msbesi), 2- 2-nd kind (msbesk) 


d_order=9
i_order=9
        tools_info=(/40,d_order,i_order/)!optimal 8,7
 !info about tools_info array:
 !tools_info(1) - (2-nd dimenstion size of tools array 1-10 - coefficinets for Lagrange interp,
 !                11-20 coefficints for (ri+1/4) interpolation, and 21-30 for (ri+2/4), 31-40 for (ri+3/4) 
 !tools_info(2) - order of intrpolation polynom for Lagrange interpolation when calculating derivative,
 !tools_info(3) - order for inperolation for integration with Bodes folmula)

Allocate(tools(Ngrid,tools_info(1)))
call generate_tools(Ngrid,r,tools,tools_info)

write(*,*)"lmax=",lmax



! ----- version x4 -------
! Lippmann–Schwinger iterations and solving screened Poisson equation.
! $\left( \nabla^2+ 2\epsilon \right \psi(\mathbf{r}) = v(\mathbf{r}) \psi(\mathbf{r})$



Allocate(psip(Ngrid,Nshell),vx_psi(Ngrid,Nshell),vx_psi_sr(Ngrid,Nshell),vx_psi_lr(Ngrid,Nshell),&
        eigp(Nshell),vhp(Ngrid),vxcp(Ngrid))

call iteration0(Ngrid,r,Z,Nshell,shell_l,lmax,count_l,psi_eig,psi)
do ish=1,Nshell
  call integ_BodesN_value(Ngrid,r,tools, tools_info,psi(:,ish)**2*r**2,norm)
  write(*,*)ish,". "," eig=",psi_eig(ish)," norm=",norm
enddo

eigp=psi_eig*0d0
vh=0d0*r
vxc=0d0*r

if (abs(hybx_w(5,1)).gt.1d-20) then
 call errfun(Ngrid,r,Nrsfun,hybx_w(5,2),rsfunC)
 call get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
endif

!Some initialisation for libxc Hybrid exchange
if(x_num.gt.0)then
if (.not.(override_libxc_hyb)) then
  if  (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_HYB_GGA) then
    write(*,*)"Hybrid Exchange!"
    e1= xc_f03_hyb_exx_coef(x_func) !hyb_exx - exact-exchange - Foka apmaiņas svars
    write(*,*) "Foka apmaiņas svars: ",e1
    call xc_f03_hyb_cam_coef(x_func , rsmu, a1, a2);
    write(*,*) "rs parameter mu:",rsmu
    write(*,*) "rs alpha",a1
    write(*,*) "rs beta",a2
    hybx_w(4,1)=a1
    hybx_w(5,1)=a2
    hybx_w(5,2)=rsmu

    if (abs(hybx_w(5,1)).gt.1d-20) then
         call errfun(Ngrid,r,Nrsfun,hybx_w(5,2),rsfunC)
         call get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
    endif
  endif
endif
endif



do iscl=1,100
!write(*,*)"STARTING LOOP"
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


rho=rho/(4d0*Pi)
vx=0d0*r
vc=0d0*r
ex=0d0*r
ec=0d0*r
vc1=0d0*r
ec1=0d0*r
!! EXCHANGE  !!!

if (x_num.gt.0)then


  if (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_LDA) then
     call xc_f03_lda_exc_vxc(x_func, Ngrid, rho(1), ex(1),vx(1))
  elseif  (xc_f03_func_info_get_family(x_info).eq.XC_FAMILY_GGA) then
     !if(x_num.eq.524)then
     !  hybx_w(3)=0.11d0
     !  write(*,*)"RS sparameter for local exchange:", hybx_w(3)
     !  call xc_f03_func_set_ext_params(x_func,hybx_w(3))
     !endif
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
!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
    call rderivative_lagrN(Ngrid,r,tools,tools_info,vxsigma,ftemp1)
    vx=vx-2d0*(grho*ftemp1+vxsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!
  else
    write(*,*)"Exchange not supported!"
  endif
elseif (x_num.eq.-1) then !Built-in LDA
   call getvxc(Ngrid,rho,vx,ex)
elseif (x_num.eq.-2) then !Built-in GGA
  call getxc_pbe(Ngrid,r,tools,tools_info,rho,vx,ex)
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
!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
    call rderivative_lagrN(Ngrid,r,tools,tools_info,vcsigma,ftemp1)
    vc1=vc1-2d0*(grho*ftemp1+vcsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!
  else
    write(*,*)"Correlation not supported!"
  endif
endif

!!!!!! End Papildus COREL

vxc=hybx_w(1,1)*vx+hybx_w(2,1)*vc+hybx_w(3,1)*vc1
exc=hybx_w(1,1)*ex+hybx_w(2,1)*ec+hybx_w(3,1)*ec1

rho=rho*(4d0*Pi)


! Get non-local exchange

if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
 do ish=1,Nshell
 call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,lmax,psi(:,ish),psi,&
           vx_psi(:,ish),vx_psi_sr(:,ish),rsfunC,Nrsfun,hybx_w,Bess_ik)
  enddo 

endif
  
! Get Coulomb potential

call  integ_BodesN_fun(Ngrid,r,tools,tools_info,1,r**2*rho,ftemp1)
call  integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,r*rho,ftemp2)
vhp=vh
vh=ftemp1/r+ftemp2
vh=vh*alpha+(1d0-alpha)*vhp


!!!!!!!!! Caulculate energy !!!!!!!!!!!!!
e_kin=0d0
do ish=1,Nshell
  call rderivative_lagrN(Ngrid,r,tools,tools_info,psi(:,ish),ftemp1)
  call integ_BodesN_value(Ngrid,r,tools,tools_info,0.5d0*ftemp1**2*r**2,e1)
  call integ_BodesN_value(Ngrid,r,tools, tools_info,0.5d0*dble(shell_l(ish))*dble(shell_l(ish)+1)*psi(:,ish)**2,e2)
  e_kin=e_kin+shell_occ(ish)*(e1+e2)
enddo

call integ_BodesN_value(Ngrid,r,tools,tools_info,-r*rho*Z,e_ext)

call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*rho*vh,e_h)

e2=0d0
  do ish=1,Nshell
    call integ_BodesN_value(Ngrid,r,tools, tools_info,r**2*0.5d0*shell_occ(ish)*psi(:,ish)*&
            (hybx_w(4,1)*vx_psi(:,ish)+hybx_w(5,1)*vx_psi_sr(:,ish)),e1)
    e2=e2+e1
  enddo
  call integ_BodesN_value(Ngrid,r,tools, tools_info,(exc)*rho*r**2,e3)

 e_x=e2+e3
 energy0=energy
 energy=e_kin+e_ext+e_h+e_x

!!!!!!!!! End Caulculate energy !!!!!!!!!!!!!


write(*,*)"Energy_tot:",energy, " (",iscl,")"," dE=",energy-energy0, " e_pot/e_kin+2=",(e_ext+e_h+e_x)/e_kin+2d0


  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0,e_ext,e_h,e_x,e1
  close(11)

  open(11,file='eigval_de.out',status='old', access='append')
  do ish=1,Nshell
    write(11,*)iscl, psi_eig(ish), psi_eig(ish)-eigp(ish), shell_occ(ish)
  enddo
  close(11)



!Check convergence

!write(*,*)iscl,". " ,psi_eig,eigp, "psi_eig eigp"
if(((maxval(abs((psi_eig-eigp)/(psi_eig+1d0)))).lt.1d-13).and.(abs(energy-energy0).lt.1d-6))then !1d-13 bija
        !(abs(energy-energy0).lt.1d-6) ir lai He nekonverģē pēc 1. iterācijas
        write(*,*)"Arējais cikls konverģē:"
        write(*,*)"konverģence absolūtā: ",maxval(abs(psi_eig-eigp))," relatīvā: ",maxval(abs((psi_eig-eigp)/(psi_eig+1d0)))
        exit
endif

l_n=0
psip=psi
eigp=psi_eig
  do il=1,lmax+1
  call LS_iteration(Ngrid,r,tools,tools_info,rsfunC,Nrsfun,hybx_w,Z,il-1,shell_l,count_l(il),l_n,&
            Nshell,lmax,vxc,vh,vx_psi,vx_psi_sr,vx_psi_lr,psip,norm_arr,psi,psi_eig,Bess_ik)
    do inn=1,count_l(il)
       l_n=l_n+1
    enddo
  enddo


enddo
!End of self consistent loop

call timesec(time1)



write(*,*)"RESULT::"
l_n=0
do il=1,lmax+1
  do inn=1,count_l(il)
    l_n=l_n+1
    write(*,*)"l=",il-1," n=",inn," eig=",psi_eig(l_n)
  enddo
enddo

  Deallocate(vx_psi,psip,eigp,vxcp,vhp,tools)






endif



!  write results 
  open(11,file='res.dat',status='old', access='append')
  write(11, '(i1,a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a1,i5,a1)',advance="no") grid,",",iscl,&
          ",",x_num,",",c_num,",",c1_num,",",d_order,",",i_order,",",Ngrid,","
  write(11, '(ES9.2E2,a1)',advance="no")Rmin,","
 write(11, '(f5.2)' ,advance="no") Rmax
 write(11, '(a1)',advance="no")","
  write(11, '(f5.2,a1)',advance="no")Z,","
  write(11, '(f7.2,a1,i1,a1,f5.2,a1,f5.2,a1,f4.2,a1,ES9.2E2)',advance="no") time1-time0,",",Nrsfun,",",hybx_w(4,1),",",&
          hybx_w(5,1),",",hybx_w(5,2),",",energy-energy0
  write(11, *)",",energy
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

