program atomHF


!!!!!!!!!! libxc !!!!!!!!!
  use xc_f03_lib_m
  implicit none
  TYPE(xc_f03_func_t) :: xc1_func,xc2_func,xc3_func
  TYPE(xc_f03_func_info_t) :: xc1_info,xc2_info,xc3_info
  integer :: vmajor, vminor, vmicro
  character(len=120) :: kind,family,version_text
  logical :: polarised
  real(8) :: hybx_w(5,2) !hyb exchange weigths: (1,1)-vx, (2,1)-vc, (3,1)-vc1, (4,1)-HF, (5,1)-HF_sr
!!!!!!!!!! libxc !!!!!!!!





real(8), PARAMETER :: Pi = 3.1415926535897932384d0

real(8), allocatable :: r(:),vh(:),vxc(:),exc(:),psi(:,:),rho(:),ftemp1(:),&
        ftemp2(:),vx_psi(:,:),&
        grho(:),grho2(:),g2rho(:),g3rho(:),vxc1(:),vxc2(:),vxc3(:),&
        exc1(:),exc2(:),exc3(:),vxcsigma(:),&
        vx_psi_sr(:,:)
integer, allocatable :: shell_n(:),shell_l(:),count_l(:)

real(8), allocatable :: shell_occ(:),eig(:)
real(8), allocatable :: psip(:,:),eigp(:)

complex(8), allocatable :: Bess_ik(:,:,:,:)

integer :: il,icl,il_icl,iscl,lmax,l_n,inn, positive_eig_iter
real(8) :: Z,rez,a1,a2
real(8) :: norm,Rmin,Rmax,hh,dE_min,e1,e2,e3,energy,energy0,e_kin,e_ext,e_h,e_x,e_pot
integer :: grid, Nshell, ish, d_order,i_order
integer :: ir,i,j,countl0, version,tools_info(3)
integer(8) :: Ngrid 
logical :: file_exists, sorted
character(len=1024) :: filename

Real (8)  :: time0,time1, t11
real(8), allocatable :: tools(:,:)
integer :: xc1_num,xc2_num,xc3_num

real(8) :: rsmu
integer :: Nrsfun
complex(8), allocatable :: rsfunC(:,:)

logical :: override_libxc_hyb
integer :: param_nr1,param_nr2,param_nr3
real(8), allocatable :: param1(:), param2(:), param3(:)

positive_eig_iter=0

call timesec(time0)

!read input
read(*,*) 
read(*,*) Z, Rmin, Rmax, Ngrid, version
read(*,*)
read(*,*) xc1_num, hybx_w(1,1), param_nr1
allocate(param1(param_nr1))
do i=1, param_nr1
read(*,*) param1(i)
enddo

read(*,*) xc2_num, hybx_w(2,1), param_nr2
allocate(param2(param_nr2))
do i=1, param_nr2
read(*,*) param2(i)
enddo

read(*,*) xc3_num,hybx_w(3,1), param_nr3
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
allocate(shell_n(Nshell),shell_l(Nshell),shell_occ(Nshell))
read(*,*)
do ish=1,Nshell
  read(*,*) shell_n(ish), shell_l(ish), shell_occ(ish)
enddo

!check if list of shells is in correct order by l and sort
sorted=.false.
do while (.not.sorted)
sorted=.true.
  do ish=1,Nshell-1
    if(shell_l(ish).gt.shell_l(ish+1))then
      sorted=.false.
      il=shell_l(ish)
      shell_l(ish)=shell_l(ish+1)
      shell_l(ish+1)=il
      inn=shell_n(ish)
      shell_n(ish)=shell_n(ish+1)
      shell_n(ish+1)=inn
      e1=shell_occ(ish)
      shell_occ(ish)=shell_occ(ish+1)
      shell_occ(ish+1)=e1
    endif
  enddo
enddo


lmax=maxval(shell_l)
allocate(count_l(lmax+1))



!counts how many particular l shells
do il=0,lmax
 countl0=0
 do ish=1,Nshell
 if (il.eq.shell_l(ish)) then
         countl0=countl0+1
 endif
 enddo
 count_l(il+1)=countl0
enddo
write(*,*)"**********Shell configuration info:************"
write(*,*)"Shell count for each l:"
do il=0, lmax
 write(*,*) "l=",il," shell count=",count_l(il+1)
enddo

!check if list of shells is in correct order by n for each l and sort
l_n=1
do il=0, lmax
 if (count_l(il+1).gt.1) then
         sorted=.false.
 else
         sorted=.true.
 endif
 do while (.not.sorted)
   sorted=.true.
   do ish=l_n,l_n+count_l(il+1)-2
     if(shell_n(ish).gt.shell_n(ish+1)) then
      sorted=.false.
      inn=shell_l(ish)
      shell_l(ish)=shell_l(ish+1)
      shell_l(ish+1)=inn
      inn=shell_n(ish)
      shell_n(ish)=shell_n(ish+1)
      shell_n(ish+1)=inn
      e1=shell_occ(ish)
      shell_occ(ish)=shell_occ(ish+1)
      shell_occ(ish+1)=e1
     endif
   enddo  
 enddo
l_n=l_n+count_l(il+1)
enddo

write(*,*)"Shell configuration from input sorted:"
do ish=1,Nshell
  write(*,*) ish,". n: ", shell_n(ish), " l: ",  shell_l(ish), " occ: ",shell_occ(ish)
enddo

write(*,*) "sum(occ)=", sum(shell_occ)," Z=", Z
write(*,*)
write(*,*)"**********XC functional info:************"
!!!!!!!!!! libxc !!!!!!!!!
if ((version.eq.0).or.(version.eq.1)) then 
! Print out libxc version
          call xc_f03_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro

if (xc1_num.gt.0) then
!!!!info 1-st XC functional 
  write(*,*)"1-st XC functional number: ",xc1_num

  call xc_f03_func_init(xc1_func, xc1_num, XC_UNPOLARIZED)

  xc1_info = xc_f03_func_get_info(xc1_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(xc1_info))
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
  select case (xc_f03_func_info_get_family(xc1_info))
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
    trim(xc_f03_func_info_get_name(xc1_info)), trim(kind), trim(family)

  i = 0
  if(xc1_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(xc1_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(xc1_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(xc1_info),  "external parameters."
   
  do i=0, xc_f03_func_info_get_n_ext_params(xc1_info)-1
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(xc1_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(xc1_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(xc1_info,i))
 enddo
 endif

  
if (xc2_num.gt.0) then
!!!!info 2-nd XC functional 
write(*,*)"2-nd XC functional number: ",xc2_num
  call xc_f03_func_init(xc2_func, xc2_num, XC_UNPOLARIZED)

  xc2_info = xc_f03_func_get_info(xc2_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(xc2_info))
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
  select case (xc_f03_func_info_get_family(xc2_info))
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
    trim(xc_f03_func_info_get_name(xc2_info)), trim(kind), trim(family)
  ! Print out references
  i = 0
  if(xc2_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(xc2_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(xc2_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(xc2_info),  "external parameters."

  do i=0, xc_f03_func_info_get_n_ext_params(xc2_info)-1
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(xc2_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(xc2_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(xc2_info,i))
  enddo
  endif
!!!!info 3-rd XC functional  
if (xc3_num.gt.0) then
write(*,*)"3rd XC functional number: ",xc3_num

  call xc_f03_func_init(xc3_func, xc3_num, XC_UNPOLARIZED)

  xc3_info = xc_f03_func_get_info(xc3_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(xc3_info))
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
  select case (xc_f03_func_info_get_family(xc3_info))
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
    trim(xc_f03_func_info_get_name(xc3_info)), trim(kind), trim(family)
  ! Print out references
  i = 0
  if(xc3_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(xc3_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(xc3_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(xc3_info),  "external parameters."

  do i=0, xc_f03_func_info_get_n_ext_params(xc3_info)-1
  
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(xc3_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(xc3_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(xc3_info,i))
 enddo

endif


 if (param_nr1.ne.0)then
   if (xc_f03_func_info_get_n_ext_params(xc1_info).ne.param_nr1) then
         write(*,*)" Parameter nr for Functional ",xc1_num, "should be", xc_f03_func_info_get_n_ext_params(xc1_info),&
                 ", not ",param_nr1
         stop
   else
    call xc_f03_func_set_ext_params(xc1_func,param1(1))
   endif     
 endif

 if (param_nr2.ne.0)then
   if (xc_f03_func_info_get_n_ext_params(xc2_info).ne.param_nr2) then
         write(*,*)" Parameter nr for Functional ",xc2_num, "should be", xc_f03_func_info_get_n_ext_params(xc2_info),&
                 ", not ",param_nr2
         stop
   else
    call xc_f03_func_set_ext_params(xc2_func,param2(1))
   endif
 endif

 if (param_nr3.ne.0)then
   if (xc_f03_func_info_get_n_ext_params(xc3_info).ne.param_nr3) then
         write(*,*)" Parameter nr for Functional ",xc3_num, "should be", xc_f03_func_info_get_n_ext_params(xc3_info),&
                 ", not ",param_nr3
         stop
   else
    call xc_f03_func_set_ext_params(xc3_func,param3(1))
   endif
 endif


  !!!!!!!!!! libxc !!!!!!!!!


write(*,*)
write(*,*)"**********Non-local exchange info************"
!Store non-local parameters in hybx_w(4,1), hybx_w(5,1) and hybx_w(5,2) 
if (override_libxc_hyb) then
   write(*,*)"Non-local exchange parameters taken from input:"
   write(*,*)"Fock exchange with weigth: ",hybx_w(4,1)
   write(*,*)"Fock SR exchange with weigth: ",hybx_w(5,1), " parameter: ", hybx_w(5,2)
else
  if(xc1_num.gt.0)then
    if (xc_f03_func_info_get_family(xc1_info).eq.XC_FAMILY_HYB_GGA) then
      write(*,*)"XC funcitional ",xc1_num , " contains non-local part."
      !e1= xc_f03_hyb_exx_coef(xc1_func) !hyb_exx - exact exchange (Fock weigth)  
      call xc_f03_hyb_cam_coef(xc1_func ,hybx_w(5,2),  hybx_w(4,1), hybx_w(5,1))
      write(*,*)"Fock exchange with weigth: ",hybx_w(4,1)
      write(*,*)"Fock SR exchange with weigth: ",hybx_w(5,1), " parameter: ", hybx_w(5,2)
    else
      write(*,*)"XC funcitional ",xc1_num , " without non-local part."
      hybx_w(5,1)=0d0
      hybx_w(4,1)=0d0
    endif 
  else
    write(*,*)"Calculation will be done without non-local part."
    hybx_w(5,1)=0d0
    hybx_w(4,1)=0d0
  endif 

  if(xc2_num.gt.0)then
    if (xc_f03_func_info_get_family(xc2_info).eq.XC_FAMILY_HYB_GGA) then
      write(*,*)"XC funcitional ",xc2_num ," in 2-nd position of input is hybrid, but non-local part will not be used."
      write(*,*)"1-st position of xc functional is for hybrid."
    endif
  endif

  if(xc3_num.gt.0)then
    if (xc_f03_func_info_get_family(xc3_info).eq.XC_FAMILY_HYB_GGA) then
      write(*,*)"XC funcitional ",xc3_num ," in 3-rd position of input is hybrid, but non-local part will not be used."
      write(*,*)"1-st position of xc functional is for hybrid."
    endif
  endif
endif !override

!End store non-local parameters

endif


if (version.eq.1) then

  inquire(file='eigval_de.out',EXIST=file_exists)
  if (file_exists) then
     open(11,file='eigval_de.out',status='old', access='append')
  else
     open(11,file='eigval_de.out',status='new')
  endif
  write(11,*)"eigval de"
  close(11)

  inquire(file='E_dE.out',EXIST=file_exists)
  if (file_exists) then
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
     write(11,*)"grid,iter,posit_eig_iter,xc1_num,xc2_num,c2_num,d_order,i_order,Ngrid,Rmin,Rmax,&
             Z,time,Nrsfuni&
             ,Fock_w,Fock_rs_w,rs_mu,",&
             "dE,energy"
  endif
  close(11)

allocate(r(Ngrid),vh(Ngrid),rho(Ngrid),vxc(Ngrid),exc(Ngrid),vxc1(Ngrid),exc1(Ngrid),&
        vxc2(Ngrid),exc2(Ngrid),vxc3(Ngrid),exc3(Ngrid),psi(Ngrid,Nshell),eig(Nshell),&
        grho2(Ngrid),ftemp1(Ngrid),ftemp2(Ngrid),vxcsigma(Ngrid),grho(Ngrid),g2rho(Ngrid),&
        psip(Ngrid,Nshell),vx_psi(Ngrid,Nshell),vx_psi_sr(Ngrid,Nshell),eigp(Nshell))

call gengrid(grid,Ngrid,Rmin,Rmax,r)


d_order=9
i_order=9
        tools_info=(/40,d_order,i_order/)
 !info about tools_info array:
 !tools_info(1) - (2-nd dimenstion size of tools array 1-10 - coefficinets for Lagrange interp,
 !                11-20 coefficints for (ri+1/4) interpolation, and 21-30 for (ri+2/4), 31-40 for (ri+3/4) 
 !tools_info(2) - order of intrpolation polynom for Lagrange interpolation when calculating derivative,
 !tools_info(3) - order for inperolation for integration with Bodes folmula)

Allocate(tools(Ngrid,tools_info(1)))
call generate_tools(Ngrid,r,tools,tools_info)

Nrsfun=8 
allocate(rsfunC(Nrsfun,2))
allocate(Bess_ik(Ngrid,Nrsfun,2*lmax+1,2)) !last atgument 1- 1-st kind (msbesi), 2- 2-nd kind (msbesk) 

if (abs(hybx_w(5,1)).gt.1d-20) then
 call errfun(Ngrid,r,Nrsfun,hybx_w(5,2),rsfunC)
 call get_Bess_fun(Ngrid,r,lmax,Nrsfun,rsfunC,Bess_ik)
endif

write(*,*)
write(*,*)"**********Starting self consistent loop************"
! Lippmann–Schwinger iterations and solving screened Poisson equation.

call iteration0(Ngrid,r,Z,Nshell,shell_l,lmax,count_l,eig,psi)
eigp=eig*0d0

! $\left( \nabla^2+ 2\epsilon \right \psi(\mathbf{r}) = v(\mathbf{r}) \psi(\mathbf{r})$

!START self consistent loop
do iscl=1,150

! Calculate density
  rho=0d0*rho
  do ish=1,Nshell
     rho=rho+shell_occ(ish)*psi(:,ish)**2
  enddo


rho=rho/(4d0*Pi)
vxc1=0d0*r
exc1=0d0*r
vxc2=0d0*r
exc2=0d0*r
vxc3=0d0*r
exc3=0d0*r

!! EXCHANGE-CORRELATION 1  !!!

if (xc1_num.gt.0)then


  if (xc_f03_func_info_get_family(xc1_info).eq.XC_FAMILY_LDA) then
     call xc_f03_lda_exc_vxc(xc1_func, Ngrid, rho(1), exc1(1),vxc1(1))
  elseif  ((xc_f03_func_info_get_family(xc1_info).eq.XC_FAMILY_GGA).or.&
                 (xc_f03_func_info_get_family(xc1_info).eq.XC_FAMILY_HYB_GGA)) then
     call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
     call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
     g2rho=g2rho*r**(-2)
     grho2=grho**2
     call xc_f03_gga_exc_vxc(xc1_func, Ngrid, rho(1), grho2(1), exc1(1),vxc1(1),vxcsigma(1))

  !!!!! formula 6.0.8 !!!!!!
!  call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho*vxcsigma,ftemp1)
!  ftemp1=ftemp1*r**(-2)
!  vxc1=vxc1-2d0*ftemp1
  !!!!! end formula 6.0.8 !!!!!!

!   !!!!! formula 6.0.5 (exciting variants) !!!!!!
      call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma,ftemp1)
     vxc1=vxc1-2d0*(grho*ftemp1+vxcsigma*g2rho)
!   !!!!! end formula 6.0.5 (exciting variants) !!!!!!

  else
    write(*,*)"1-st XC functional not supported!"
  endif
elseif (xc1_num.eq.-1) then !Built-in LDA
   call getvxc(Ngrid,rho,vxc1,exc1)
elseif (xc1_num.eq.-2) then !Built-in GGA
  call getxc_pbe(Ngrid,r,tools,tools_info,rho,vxc1,exc1)
endif
!! END EXCHANGE-CORRELATION 1 !!!


!! EXCHANGE-CORRELATION 2 !!!
if (xc2_num.ne.0) then
  if (xc_f03_func_info_get_family(xc2_info).eq.XC_FAMILY_LDA) then
    call xc_f03_lda_exc_vxc(xc2_func, Ngrid, rho(1), exc2(1),vxc2(1))
  elseif  ((xc_f03_func_info_get_family(xc2_info).eq.XC_FAMILY_GGA).or.&
                 (xc_f03_func_info_get_family(xc2_info).eq.XC_FAMILY_HYB_GGA)) then
    call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
    g2rho=g2rho*r**(-2)
    grho2=grho**2
    call xc_f03_gga_exc_vxc(xc2_func, Ngrid, rho(1), grho2(1), exc2(1),vxc2(1),vxcsigma(1))
    call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma,ftemp1)
    vxc2=vxc2-2d0*(grho*ftemp1+vxcsigma*g2rho)
  else
    write(*,*)"2-nd XC functional not supported!"
  endif
endif
!! END EXCHANGE-CORRELATION 2 !!!!!

!! EXCHANGE-CORRELATION 3 !!!
if (xc3_num.ne.0) then
  if (xc_f03_func_info_get_family(xc3_info).eq.XC_FAMILY_LDA) then
    call xc_f03_lda_exc_vxc(xc3_func, Ngrid, rho(1), exc3(1),vxc3(1))
  elseif ((xc_f03_func_info_get_family(xc3_info).eq.XC_FAMILY_GGA).or.&
                 (xc_f03_func_info_get_family(xc3_info).eq.XC_FAMILY_HYB_GGA)) then
    call rderivative_lagrN(Ngrid,r,tools,tools_info,rho,grho)
    call rderivative_lagrN(Ngrid,r,tools,tools_info,r**2*grho,g2rho)
    g2rho=g2rho*r**(-2)
    grho2=grho**2
    call xc_f03_gga_exc_vxc(xc3_func, Ngrid, rho(1), grho2(1), exc3(1),vxc3(1),vxcsigma(1))
    call rderivative_lagrN(Ngrid,r,tools,tools_info,vxcsigma,ftemp1)
    vxc3=vxc3-2d0*(grho*ftemp1+vxcsigma*g2rho)
  else
    write(*,*)"3-rd XC functional not supported!"
  endif
endif

!!!! END EXCHANGE-CORRELATION 3 !!!

vxc=hybx_w(1,1)*vxc1+hybx_w(2,1)*vxc2+hybx_w(3,1)*vxc3
exc=hybx_w(1,1)*exc1+hybx_w(2,1)*exc2+hybx_w(3,1)*exc3

rho=rho*(4d0*Pi)


! Get non-local exchange

if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
 do ish=1,Nshell
 call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ,lmax,psi(:,ish),psi,&
           vx_psi(:,ish),vx_psi_sr(:,ish),rsfunC,Nrsfun,hybx_w,Bess_ik)
  enddo 

endif
  
! Get Coulomb potential

call  integ_BodesN_fun(Ngrid,r,tools,tools_info,1,r**2*rho,ftemp1)
call  integ_BodesN_fun(Ngrid,r,tools,tools_info,-1,r*rho,ftemp2)
vh=ftemp1/r+ftemp2


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

if (maxval(eig).gt.0) then 
  positive_eig_iter=iscl
  write(*,*)"Energy_tot:",energy, " (",iscl,")"," dE=",energy-energy0, " POSITIVE EIGENVALUE! "
else
  write(*,*)"Energy_tot:",energy, " (",iscl,")"," dE=",energy-energy0
endif

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0,e_ext,e_h,e_x,e1
  close(11)

  open(11,file='eigval_de.out',status='old', access='append')
  do ish=1,Nshell
    write(11,*)iscl, eig(ish), eig(ish)-eigp(ish), shell_occ(ish)
  enddo
  close(11)



!Check convergence

!write(*,*)iscl,". " ,eig,eigp, "eig eigp"
if(((maxval(abs((eig-eigp)/(eig-1d0)))).lt.1d-13).and.(abs(energy-energy0).lt.1d-6))then !1d-13 for best result 
        !without (abs(energy-energy0).lt.1d-6) Helium converges after 1-st iteration.
        write(*,*)"Convergence of external cycle reached:"
        write(*,*)"max(eig-eigp) absolute : ",maxval(abs(eig-eigp))," relative: ",maxval(abs(eig-eigp)/(abs(eig-1d0)))
        exit
else
        if (iscl.gt.90)then
         ish=maxloc(abs((eig-eigp)/(eig-1d0)),Nshell)
       !  write(*,*)"eigenvalue with largest relative ", maxval(abs((eig-eigp)/(eig-1d0)))," and absolute ",&
       !          maxval(abs(eig-eigp)), "change:"
       !  write(*,*)"n=",shell_n(ish)," l=",shell_l(ish), " eig=", eig(ish) 
        endif
endif

l_n=0
psip=psi
eigp=eig
  do il=1,lmax+1
  call LS_iteration(Ngrid,r,tools,tools_info,rsfunC,Nrsfun,hybx_w,Z,il-1,shell_l,shell_occ,count_l(il),l_n,&
            Nshell,lmax,vxc,vh,vx_psi,vx_psi_sr,psip,psi,eig,Bess_ik)
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
    write(*,*)"l=",il-1," n=",inn," eig=",eig(l_n)," occ=",shell_occ(l_n)
  enddo
enddo

  Deallocate(tools,rsfunC,Bess_ik)





 
endif



!  write results 
  open(11,file='res.dat',status='old', access='append')
  write(11, '(i1,a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a1,i5,a1)',advance="no") grid,",",iscl,",",positive_eig_iter,&
          ",",xc1_num,",",xc2_num,",",xc3_num,",",d_order,",",i_order,",",Ngrid,","
  write(11, '(ES9.2E2,a1)',advance="no")Rmin,","
 write(11, '(f5.2)' ,advance="no") Rmax
 write(11, '(a1)',advance="no")","
  write(11, '(f6.2,a1)',advance="no")Z,","
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
      write(11,*) icl, il, eig(il_icl) 
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



deallocate(r,vh,rho,vxc,exc,vxc1,exc1,vxc2,exc2,vxc3,exc3,psi,eig,&
        grho2,ftemp1,ftemp2,vxcsigma,grho,g2rho,psip,vx_psi,vx_psi_sr,eigp)  
  
deallocate(shell_n,shell_l,count_l,shell_occ,param1,param2,param3)


end program

