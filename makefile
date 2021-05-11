F90 = gfortran 
F90_OPTS = -O3 -cpp -I /home/janis/atom-HFtest/atom-HF-libxc/include/ 

SRCS := $(wildcard *.f90)
OBJECTS := $(SRCS:%.f90=%.o) -llapack -lblas /home/janis/atom-HFtest/atom-HF-libxc/lib/libxcf03.a /home/janis/atom-HFtest/atom-HF-libxc/lib/libxc.a 

all: ${OBJECTS} link

link:
	$(F90) -o atomHF ${OBJECTS} 

%.o: %.f90
	${F90} $(F90_OPTS) -c  $<

clean:
	rm -rvf *.o ${BINS}
