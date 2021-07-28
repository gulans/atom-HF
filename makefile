F90 = gfortran 
LIBXC=/home/janis/app/libxc/
F90_OPTS = -O3 -cpp -I $(LIBXC)/include/ 

SRCS := $(wildcard *.f90)
OBJECTS := $(SRCS:%.f90=%.o) -llapack -lblas $(LIBXC)/lib/libxcf03.a $(LIBXC)/lib/libxc.a 

all: ${OBJECTS} link

link:
	$(F90) -o atomHF ${OBJECTS} 

%.o: %.f90
	${F90} $(F90_OPTS) -c  $<

clean:
	rm -rvf *.o ${BINS}
