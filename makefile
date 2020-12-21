F90 = gfortran 
F90_OPTS = -cpp

SRCS := $(wildcard *.f90)
OBJECTS := $(SRCS:%.f90=%.o)

all: ${OBJECTS} link

link:
	$(F90) -o atomHF ${OBJECTS} 

%.o: %.f90
	${F90} $(F90_OPTS) -c $<

clean:
	rm -rvf *.o ${BINS}
