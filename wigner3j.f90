subroutine wigner3j(j1,j2,j3,rez)
implicit none
real(8), PARAMETER :: Pi = 3.1415926535897932384d0
integer, intent(in) :: j1,j2,j3
real(8), intent(out) :: rez

!integer ::n,i, j1a(23),j2a(23),j3a(23)
!real(8) :: wa(23)
integer ::n,i, j1a(106),j2a(106),j3a(106)
real(8) :: wa(106)

n=106

!j1a=(/0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3/)
!j2a=(/0, 1, 2, 3, 0, 1, 1, 2, 2, 3, 0, 1, 1, 2, 2, 3, 3, 0, 1, 2, 2, 3, 3/)
!j3a=(/0, 1, 2, 3, 1, 0, 2, 1, 3, 2, 2, 1, 3, 0, 2, 1, 3, 3, 2, 1, 3, 0, 2/)
!wa=(/1.0d0, -0.5773502691896258d0, 0.447213595499958d0, -0.3779644730092272d0, -0.5773502691896258d0,&
!        -0.5773502691896258d0, 0.36514837167011077d0, 0.36514837167011077d0, -0.29277002188455997d0,&
!        -0.29277002188455997d0, 0.4472135954999579d0, 0.3651483716701108d0, -0.29277002188456d0,&
!        0.4472135954999579d0, -0.23904572186687875d0, -0.29277002188456d0, 0.1951800145897066d0,&
!        -0.3779644730092272d0, -0.29277002188455997d0, -0.29277002188455997d0, 0.19518001458970663d0,&
!        -0.3779644730092272d0, 0.19518001458970663d0/)

j1a=(/0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,&
        2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,&
        4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6,&
        6, 6, 6, 6, 6, 6, 6, 6, 6/)
j2a=(/0, 1, 2, 3, 4, 5, 6, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5,&
        6, 6, 0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,&
        4, 5, 5, 5, 6, 6, 6, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 0, 1, 2, 2, 3, 3, 4,&
        4, 4, 5, 5, 5, 6, 6, 6, 6/)
j3a=(/0, 1, 2, 3, 4, 5, 6, 1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 2, 1, 3, 0, 2, 4, 1, 3, 5, 2, 4, 6, 3, 5,&
        4, 6, 3, 2, 4, 1, 3, 5, 0, 2, 4, 6, 1, 3, 5, 2, 4, 6, 3, 5, 4, 3, 5, 2, 4, 6, 1, 3, 5, 0, 2, 4,&
        6, 1, 3, 5, 2, 4, 6, 5, 4, 6, 3, 5, 2, 4, 6, 1, 3, 5, 0, 2, 4, 6, 1, 3, 5, 6, 5, 4, 6, 3, 5, 2,&
        4, 6, 1, 3, 5, 0, 2, 4, 6/)
wa=(/1.0d0, -0.5773502691896258d0, 0.447213595499958d0, -0.3779644730092272d0, 0.33333333333333326d0, -0.3015113445777636d0,&
        0.2773500981126146d0, -0.5773502691896258d0, -0.5773502691896258d0, 0.36514837167011077d0, 0.36514837167011077d0,&
        -0.29277002188455997d0, -0.29277002188455997d0, 0.2519763153394848d0, 0.2519763153394848d0, -0.22473328748774732d0,&
        -0.22473328748774732d0, 0.20483662259967567d0, 0.20483662259967567d0, 0.4472135954999579d0, 0.3651483716701108d0,&
        -0.29277002188456d0, 0.4472135954999579d0, -0.23904572186687875d0, 0.2390457218668787d0, -0.29277002188456d0,&
        0.1951800145897066d0, -0.20806259464411975d0, 0.2390457218668787d0, -0.16988239714587514d0, 0.18698939800169145d0,&
        -0.20806259464411975d0, 0.15267620413811483d0, 0.18698939800169145d0, -0.13993005245628823d0, -0.3779644730092272d0,&
        -0.29277002188455997d0, 0.2519763153394849d0, -0.29277002188455997d0, 0.19518001458970663d0, -0.2080625946441198d0,&
        -0.3779644730092272d0, 0.19518001458970663d0, -0.16116459280507608d0, 0.18248296715045298d0, 0.2519763153394849d0,&
        -0.16116459280507608d0, 0.1413506985480439d0, -0.2080625946441198d0, 0.1413506985480439d0, -0.12773807700531709d0,&
        0.18248296715045298d0, -0.12773807700531709d0, 0.3333333333333333d0, 0.2519763153394849d0, -0.22473328748774737d0,&
        0.23904572186687875d0, -0.16988239714587516d0, 0.18698939800169143d0, 0.2519763153394849d0, -0.16116459280507606d0,&
        0.1413506985480439d0, 0.3333333333333333d0, -0.16988239714587516d0, 0.13409704688030227d0, -0.1246595986677943d0,&
        -0.22473328748774737d0, 0.1413506985480439d0, -0.11826247919781653d0, 0.18698939800169143d0, -0.1246595986677943d0,&
        0.10732145112154906d0, -0.30151134457776363d0, -0.22473328748774737d0, 0.2048366225996757d0, -0.20806259464411975d0,&
        0.1526762041381148d0, -0.20806259464411975d0, 0.14135069854804394d0, -0.1277380770053171d0, -0.22473328748774737d0,&
        0.14135069854804394d0, -0.11826247919781652d0, -0.30151134457776363d0, 0.1526762041381148d0, -0.11826247919781652d0,&
        0.10473501197846222d0, 0.2048366225996757d0, -0.1277380770053171d0, 0.10473501197846222d0, 0.2773500981126146d0,&
        0.2048366225996757d0, 0.18698939800169145d0, -0.13993005245628826d0, 0.18248296715045298d0, -0.12773807700531709d0,&
        0.18698939800169145d0, -0.1246595986677943d0, 0.10732145112154912d0, 0.2048366225996757d0, -0.12773807700531709d0,&
        0.10473501197846219d0, 0.2773500981126146d0, -0.13993005245628826d0, 0.10732145112154912d0, -0.09305950021129074d0/)


rez=0d0
do i=1,n
  if ((j1a(i).eq.j1).and.(j2a(i).eq.j2).and.(j3a(i).eq.j3)) then
          rez=wa(i)
          exit
  endif
enddo

end subroutine

