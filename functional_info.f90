subroutine functional_info(xc_num,xc_func)
use xc_f03_lib_m
implicit none

integer, intent(in) :: xc_num
TYPE(xc_f03_func_t),intent(in)  :: xc_func
TYPE(xc_f03_func_info_t) :: xc_info
integer :: i
character(len=120) :: kind,family

  xc_info = xc_f03_func_get_info(xc_func)
  ! Get the type of the functional
  select case(xc_f03_func_info_get_kind(xc_info))
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
  select case (xc_f03_func_info_get_family(xc_info))
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
    trim(xc_f03_func_info_get_name(xc_info)), trim(kind), trim(family)

  i = 0
  if(xc_num.ne.524)then
  do while(i >= 0)
    write(*, '(a,i1,2a)') '[', i+1, '] ', trim(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(xc_info, i)))
  end do
  endif
  write(*,*)"FUCTIONAL: ",trim(xc_f03_func_info_get_name(xc_info))," Supports: ",&
          xc_f03_func_info_get_n_ext_params(xc_info),  "external parameters."
   
  do i=0, xc_f03_func_info_get_n_ext_params(xc_info)-1
  write(*,*)i,". ", trim(xc_f03_func_info_get_ext_params_name(xc_info,i))," default value: ",&
         xc_f03_func_info_get_ext_params_default_value(xc_info,i)," ",&
         trim(xc_f03_func_info_get_ext_params_description(xc_info,i))
 enddo

end subroutine




