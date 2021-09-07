program atomHF


!!!!!!!!!! libxc !!!!!!!!!
  use xc_f03_lib_m
  implicit none
  TYPE(xc_f03_func_t) :: xc1_func,xc2_func,xc3_func
  TYPE(xc_f03_func_info_t) :: xc1_info,xc2_info,xc3_info
  integer :: vmajor, vminor, vmicro
  character(len=120) :: kind,family,version_text
  logical :: spin
  real(8) :: hybx_w(5,2) !hyb exchange weigths: (1,1)-vx, (2,1)-vc, (3,1)-vc1, (4,1)-HF, (5,1)-HF_sr
!!!!!!!!!! libxc !!!!!!!!


real(8), PARAMETER :: Pi = 3.1415926535897932384d0
real(8), PARAMETER :: alpha2=0.5d0*7.2973525693d-3**2 !1/(2*c^2)

real(8), allocatable :: r(:),vh(:),vhp(:),vxcp(:,:),vxc(:,:),exc(:,:),psi(:,:,:),rho(:),ftemp1(:),&
        ftemp2(:),ftemp3(:),ftemp4(:),ftemp(:),vx_psi(:,:,:),v_rel(:,:),&
        grho(:),grho2(:),g2rho(:),g3rho(:),vxc1(:),vxc2(:),vxc3(:),&
        exc1(:),exc2(:),exc3(:),vxcsigma(:),vx_psip(:,:,:),vx_psi_srp(:,:,:),&
        vx_psi_sr(:,:,:)
!spin polarised variables
real(8), allocatable :: rho_sp(:,:),vxc1_sp(:,:),vxc2_sp(:,:),vxc3_sp(:,:)
real(8), allocatable :: grho_sp(:,:),grho2_sp(:,:),vxcsigma_sp(:,:),gvxcsigma_sp(:,:),g2rho_sp(:,:)

integer, allocatable :: shell_n(:),shell_l(:),count_l(:)

real(8), allocatable :: shell_occ(:,:),eig(:,:)
real(8), allocatable :: psip(:,:,:),eigp(:,:)

complex(8), allocatable :: Bess_ik(:,:,:,:)

integer :: il,icl,il_icl,iscl,lmax,l_n,inn, positive_eig_iter
real(8) :: Z,rez,a1,a2,mixerC
real(8) :: norm,Rmin,Rmax,hh,dE_min,e1,e2,e3,energy,energy0,e_kin,e_ext,e_h,e_xc,e_pot
integer :: grid, Nshell, ish, d_order,i_order
integer :: ir,i,j,countl0, version,tools_info(3)
integer(8) :: Ngrid 
logical :: file_exists, sorted
character(len=1024) :: filename

Real (8)  :: time0,time1, t11
real(8), allocatable :: tools(:,:)
integer :: xc1_num,xc2_num,xc3_num

real(8) :: rsmu
integer,allocatable :: iner_loop(:)
integer :: Nrsfun, Nspin, isp
complex(8), allocatable :: rsfunC(:,:)

logical :: override_libxc_hyb
integer :: param_nr1,param_nr2,param_nr3
real(8), allocatable :: param1(:), param2(:), param3(:)
integer :: sh0,sh1
logical :: relativity
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
read(*,*) relativity
read(*,*)
read(*,*) Nshell
read(*,*)
read(*,*) spin
read(*,*)
if (spin)then
        Nspin=2
else
        Nspin=1
endif

allocate(shell_n(Nshell),shell_l(Nshell),shell_occ(Nshell,Nspin))
if (.not.spin) then
  do ish=1,Nshell
    read(*,*) shell_n(ish), shell_l(ish), shell_occ(ish,1)
  enddo 
else
  do ish=1,Nshell
    read(*,*) shell_n(ish), shell_l(ish), shell_occ(ish,1), shell_occ(ish,2)
  enddo
endif

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
      e1=shell_occ(ish,1)
      shell_occ(ish,1)=shell_occ(ish+1,1)
      shell_occ(ish+1,1)=e1
      if (spin) then
        e1=shell_occ(ish,2)
        shell_occ(ish,2)=shell_occ(ish+1,2)
        shell_occ(ish+1,2)=e1
      endif
    endif
  enddo
enddo


lmax=maxval(shell_l)
allocate(count_l(lmax+1),iner_loop(lmax+1))



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
      e1=shell_occ(ish,1)
      shell_occ(ish,1)=shell_occ(ish+1,1)
      shell_occ(ish+1,1)=e1
      if (spin) then
        e1=shell_occ(ish,2)
        shell_occ(ish,2)=shell_occ(ish+1,2)
        shell_occ(ish+1,2)=e1
      endif
     endif
   enddo  
 enddo
l_n=l_n+count_l(il+1)
enddo

write(*,*)"Shell configuration from input sorted:"
do ish=1,Nshell
if (.not.spin) then
  write(*,*) ish,". n: ", shell_n(ish), " l: ",  shell_l(ish), " occ: ",shell_occ(ish,1)
else
  write(*,*) ish,". n: ", shell_n(ish), " l: ",  shell_l(ish), " occ up: ",shell_occ(ish,1), " occ down: ",shell_occ(ish,2)
endif
enddo

write(*,*) "sum(occ)=", sum(shell_occ)," Z=", Z
write(*,*)
write(*,*)"**********XC functional info:************"
!!!!!!!!!! libxc !!!!!!!!!
if ((version.eq.0).or.(version.eq.1).or.(version.eq.2)) then 
! Print out libxc version
          call xc_f03_version(vmajor, vminor, vmicro)
  write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro

if (xc1_num.gt.0) then
  write(*,*)"1-st XC functional number: ",xc1_num
  if (.not.spin)then
  call xc_f03_func_init(xc1_func, xc1_num, XC_UNPOLARIZED)
  else
  call xc_f03_func_init(xc1_func, xc1_num, XC_POLARIZED)
  endif
  xc1_info = xc_f03_func_get_info(xc1_func)
  call functional_info(xc1_num,xc1_func)
endif

  
if (xc2_num.gt.0) then
!!!!info 2-nd XC functional 
write(*,*)"2-nd XC functional number: ",xc2_num
  if (.not.spin)then
  call xc_f03_func_init(xc2_func, xc2_num, XC_UNPOLARIZED)
  else
  call xc_f03_func_init(xc2_func, xc2_num, XC_POLARIZED)
  endif
  xc2_info = xc_f03_func_get_info(xc2_func)
  call functional_info(xc2_num,xc2_func)
endif

!!!!info 3-rd XC functional  
if (xc3_num.gt.0) then
write(*,*)"3rd XC functional number: ",xc3_num
  if (.not.spin)then
  call xc_f03_func_init(xc3_func, xc3_num, XC_UNPOLARIZED)
  else
  call xc_f03_func_init(xc3_func, xc3_num, XC_POLARIZED)
  endif
  xc3_info = xc_f03_func_get_info(xc3_func)
  call functional_info(xc3_num,xc3_func)
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
   if(abs(hybx_w(5,1)).lt.1d-20)then
           hybx_w(5,2)=0d0 !or else the value in output is a bit confusing
   endif
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


if ((version.eq.1).or.(version.eq.2)) then

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
  write(11,*)"iter E dE E_ext E_h E_xc eig_sum Z=",Z
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

allocate(r(Ngrid),vh(Ngrid),vhp(Ngrid),rho(Ngrid),vxc(Ngrid,Nspin),vxcp(Ngrid,Nspin),exc(Ngrid,Nspin),&
        vxc2(Ngrid),exc2(Ngrid),vxc3(Ngrid),exc3(Ngrid),psi(Ngrid,Nshell,Nspin),eig(Nshell,Nspin),&
        grho2(Ngrid),ftemp1(Ngrid),ftemp2(Ngrid),vxcsigma(Ngrid),grho(Ngrid),g2rho(Ngrid),&
        psip(Ngrid,Nshell,Nspin),vx_psi(Ngrid,Nshell,Nspin),vx_psi_sr(Ngrid,Nshell,Nspin),eigp(Nshell,Nspin),&
        v_rel(Ngrid,Nspin),ftemp3(Ngrid),ftemp4(Ngrid),ftemp(Ngrid),vx_psip(Ngrid,Nshell,Nspin),&
        vx_psi_srp(Ngrid,Nshell,Nspin))
allocate(rho_sp(Nspin,Ngrid),vxc1(Ngrid),exc1(Ngrid)) !differenet order for libxc



allocate(grho_sp(2,Ngrid),grho2_sp(3,Ngrid),vxc1_sp(2,Ngrid),vxc2_sp(2,Ngrid),vxc3_sp(2,Ngrid),&
        vxcsigma_sp(3,Ngrid),gvxcsigma_sp(3,Ngrid),g2rho_sp(2,Ngrid))


call gengrid(grid,Z,Ngrid,Rmin,Rmax,r)


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

if (version.eq.1) then
call iteration0(Ngrid,r,Z,Nshell,shell_l,lmax,count_l,Nspin,eig,psi)
endif

do isp=1, Nspin
vxc(:,isp)=r*0d0
enddo
vxcp=vxc*0d0
vh=0d0*r
vhp=0d0*r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Restart part
if (version.eq.2) then

!input WF from file
  read(*,*) !Header row

do ir=1,Ngrid
    read(*,'(ES25.16E3)',advance="no")e1 !r(ir)
  !  read(*,'(a1)',advance="no") !seperator
  do isp=1,Nspin
  do ish=1,Nshell
    read(*,'(ES25.16E3)',advance="no")psi(ir,ish,isp)
  enddo
  enddo
  read(*,*)!read end of line

  
  enddo  !ir

 ! Get non-local exchange

if ((abs(hybx_w(4,1)).gt.1d-20).or.(abs(hybx_w(5,1)).gt.1d-20)) then
 do ish=1,Nshell
 do isp=1,Nspin

 call get_Fock_ex(Ngrid,r,tools,tools_info,ish,Nshell,shell_l,shell_occ(:,isp),lmax,psi(:,ish,isp),psi(:,:,isp),&
           vx_psi(:,ish,isp),vx_psi_sr(:,ish,isp),rsfunC,Nrsfun,hybx_w,Bess_ik)
   enddo
   enddo
vx_psi=vx_psi*dble(Nspin)
vx_psi_sr=vx_psi_sr*dble(Nspin)

endif !end get Fock exchange


call get_local_exc_vxc_vh_rho(Ngrid,r,tools,tools_info,Nshell,shell_occ,spin,Nspin,psi,&
         xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func,hybx_w,exc,vxc,vh,rho)

vx_psip=vx_psi
vx_psi_srp=vx_psi_sr
do isp=1,Nspin
sh0=0
  do il=1, lmax+1
  sh1=sh0+count_l(il)
  call orthonorm_get_eig(Ngrid,r,tools,tools_info,Z,il-1,count_l(il),relativity,v_rel,hybx_w,&
        vxc(:,isp),vh,vx_psip(:,sh0+1:sh1,isp),vx_psi_srp(:,sh0+1:sh1,isp),&
        psi(:,sh0+1:sh1,isp),eig(sh0+1:sh1,isp),vx_psi(:,sh0+1:sh1,isp),vx_psi_sr(:,sh0+1:sh1,isp))

sh0=sh1
  enddo
enddo
do isp=1,Nspin
do ish=1,Nshell
  write(*,*)ish,eig(ish,isp)
enddo
enddo
endif !version 2
!END Restart part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Total energy after 0-th iteration
rho=0d0*r
do isp=1,Nspin
do ish=1,Nshell
    rho=rho+shell_occ(ish,isp)*psi(:,ish,isp)**2
enddo
enddo
call get_energy(Ngrid,r,tools,tools_info,Z,Nshell,shell_occ,spin,relativity,v_rel,Nspin,shell_l,hybx_w,&
        vxc,exc,rho,vh,vx_psi,vx_psi_sr,psi,&
        e_kin,e_ext,e_h,e_xc)
energy=e_kin+e_ext+e_h+e_xc
if (maxval(eig).gt.0) then
  positive_eig_iter=0
  write(*,*)"Energy_tot:",energy, " (", 0,")"," POSITIVE EIGENVALUE! "
!  write(*,*)"E=Ekin+Eext+Eh+Ex",energy,"=",e_kin,e_ext,e_h,e_xc
else
  write(*,*)"Energy_tot:",energy, " (", 0,")"
!  write(*,*)"E=Ekin+Eext+Eh+Ex",energy,"=",e_kin,e_ext,e_h,e_xc
endif

!End Total energy after 0-th iteration



psip=psi*0d0
eigp=eig*0d0

! $\left( \nabla^2+ 2\epsilon \right \psi(\mathbf{r}) = v(\mathbf{r}) \psi(\mathbf{r})$

!START self consistent loop
do iscl=1,200
vxcp=vxc
vhp=vh
call get_local_exc_vxc_vh_rho(Ngrid,r,tools,tools_info,Nshell,shell_occ,spin,Nspin,psi,&
         xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func,hybx_w,exc,vxc,vh,rho)
 if ((Z.gt.28.5d0).and.(Z.lt.29.5d0))then !Cu case
   mixerC=0.3d0
 else
   mixerC=0.5d0
 endif
vxc=mixerC*vxc+(1d0-mixerC)*vxcp
vh=mixerC*vh+(1d0-mixerC)*vhp
  

psip=psi
eigp=eig
l_n=0
  do il=1,lmax+1
  do isp=1,Nspin
  call LS_iteration(Ngrid,r,tools,tools_info,rsfunC,Nrsfun,hybx_w,Z,il-1,isp,shell_l,shell_occ,count_l(il),l_n,& !in
            Nshell,Nspin,relativity,lmax,vxc,v_rel(:,isp),vh,& !in
            vx_psi(:,l_n+1:l_n+count_l(il),isp),vx_psi_sr(:,l_n+1:l_n+count_l(il),isp),& !inout parameters
            psip,& !in
            psi(:,l_n+1:l_n+count_l(il),isp),eig(l_n+1:l_n+count_l(il),isp),& !inout 
            Bess_ik,iner_loop(il),& !in
            xc1_num,xc2_num,xc3_num,xc1_func,xc2_func,xc3_func) !in
  enddo
    l_n=l_n+count_l(il)
  enddo




call get_energy(Ngrid,r,tools,tools_info,Z,Nshell,shell_occ,spin,relativity,v_rel,Nspin,shell_l,hybx_w,&
        vxc,exc,rho,vh,vx_psi,vx_psi_sr,psi,&
        e_kin,e_ext,e_h,e_xc)

 energy0=energy
 energy=e_kin+e_ext+e_h+e_xc

!!!!!!!!! End Caulculate energy !!!!!!!!!!!!!

if (maxval(eig).gt.0) then 
  positive_eig_iter=iscl
  write(*,*)"Energy_tot:",energy, " (",iscl,")"," dE=",energy-energy0, " POSITIVE EIGENVALUE! ", "(",iner_loop,")"
else
  write(*,*)"Energy_tot:",energy, " (",iscl,")"," dE=",energy-energy0, "(",iner_loop,")"
 ! write(*,*)"E=Ekin+Eext+Eh+Ex",energy,"=",e_kin,e_ext,e_h,e_xc
endif

if (energy.ne.energy)then
write (*,*)"NAN detected in eigenvalues"
exit
endif 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(11,file='E_dE.out',status='old', access='append')
  write(11,*)iscl, energy, energy-energy0,e_ext,e_h,e_xc,e1
  close(11)

  open(11,file='eigval_de.out',status='old', access='append')
  
  do ish=1,Nshell
  do isp=1, Nspin
    write(11,*)iscl, eig(ish,isp), eig(ish,isp)-eigp(ish,isp), shell_occ(ish,isp)
  enddo
  enddo

  close(11)



!Check convergence

!write(*,*)iscl,". " ,eig,eigp, "eig eigp"
if(((maxval(abs((eig-eigp)/(eig-1d0)))).lt.1d-13).and.(abs(energy-energy0).lt.1d-6))then !1d-13 for best result 
        write(*,*)"Convergence of external cycle reached:"
        write(*,*)"max(eig-eigp) absolute : ",maxval(abs(eig-eigp))," relative: ",maxval(abs(eig-eigp)/(abs(eig-1d0)))
        exit
endif


enddo
!End of self consistent loop

call timesec(time1)



write(*,*)"RESULT::"
l_n=0
do il=1,lmax+1
  do inn=1,count_l(il)
    l_n=l_n+1
    if (.not.spin)then
      write(*,*)"l=",il-1," n=",inn," eig=",eig(l_n,1)," occ=",shell_occ(l_n,1)
    else
    write(*,*)"l=",il-1," n=",inn," eig=",eig(l_n,1)," occ=",shell_occ(l_n,1),"up"
    write(*,*)"l=",il-1," n=",inn," eig=",eig(l_n,2)," occ=",shell_occ(l_n,2),"down"
    endif
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
      if (.not.spin)then
      write(11,*) icl, il, eig(il_icl,1) 
      else
      write(11,*) icl, il, eig(il_icl,1),"(occ",shell_occ(il_icl,1) , ")",eig(il_icl,2),"(occ",shell_occ(il_icl,2) , ")"
      endif
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

!write all wave functions to a file
  open(11,file='wave_fun.dat',status='replace')

  write(11,'(a2)',advance="no")"r "
  do isp=1,Nspin
  do ish=1,Nshell
    if (ish.lt.10) then
      write(11,'(a2,i1,a1,i1,a1)',advance="no") "WF",ish,"-",isp," "
    else
      write(11,'(a2,i2,a1,i1,a1)',advance="no") "WF",ish,"-",isp," "
    endif
  enddo
  enddo
  write(11,*)""

  do ir=1,Ngrid
  write(11,'(ES25.16E3)',advance="no")r(ir)
  do isp=1,Nspin
  do ish=1,Nshell
    write(11,'(ES25.16E3)',advance="no") psi(ir,ish,isp)
  enddo
  enddo
  write(11,*)""
  enddo
  close(11)
!end write all wave functions to a file


deallocate(r,vh,rho,vxc,exc,vxc1,exc1,vxc2,exc2,vxc3,exc3,psi,eig,&
        grho2,ftemp1,ftemp2,vxcsigma,grho,g2rho,psip,vx_psi,vx_psi_sr,eigp)  
  
deallocate(shell_n,shell_l,count_l,shell_occ,param1,param2,param3)


end program

