F90 = gfortran 
F90_OPTS = -cpp -Ddebug -fcheck=all -g -fbacktrace 


SRCS := $(wildcard *.f90)
OBJECTS := $(SRCS:%.f90=%.o) -llapack -lblas 


all: ${OBJECTS} link

link:
	$(F90) -o atomHF ${OBJECTS} 

%.o: %.f90
	${F90} $(F90_OPTS) -c  $<

clean:
	rm -rvf *.o ${BINS}
