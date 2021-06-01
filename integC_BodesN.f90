subroutine integC_BodesN_fun(Ngrid,r,tools,tools_info,dir,fin,rez)
implicit none
integer, intent(in) :: Ngrid,tools_info(3),dir
real(8), intent(in) :: r(Ngrid),tools(Ngrid,tools_info(1)) 
complex(8), intent(in)  :: fin(Ngrid)
complex(8), intent(out) :: rez(Ngrid)
integer :: ir

real(8) :: finRe(Ngrid),finIm(Ngrid),rezRe(Ngrid),rezIm(Ngrid)

do ir=1, Ngrid
  finRe(ir)=realpart(fin(ir)) 
  finIm(ir)=imagpart(fin(ir))
enddo

call integ_BodesN_fun(Ngrid,r,tools,tools_info,dir,finRe,rezRe)
call integ_BodesN_fun(Ngrid,r,tools,tools_info,dir,finIm,rezIm)

do ir=1, Ngrid
  rez(ir)=cmplx(rezRe(ir),rezIm(ir),8)
enddo



end subroutine
