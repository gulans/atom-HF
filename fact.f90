real(8) function fact(n)
integer, intent(in) :: n

  if (n < 0) then
    write(*,*) 'ercf_kernel.f90 factorial for negative integer?',n
    stop
  endif
  fact = 1.0d0
  do i = 2, n
    fact = fact * dble(i)
  enddo
return
end function  
