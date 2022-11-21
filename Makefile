
FORT=/usr/local/bin/gfortran
.SUFFIXES: .f95 .mod .f
FFLAGS = -Wall -g -O3 -mtune=native -fbounds-check -ffpe-summary=invalid,zero,overflow -ffree-form  -fPIC
F77FLAGS = -fPIC -O3
LDFLAGS = "-L/usr/local/opt/lapack/lib" -llapack -lblas
CFLAGS = -Wall  -c

F90SRCS = Parameters.f90 integrator.f90 main.f90


OBJS = $(FORTSRCS:.f=.o) $(F90SRCS:.f90=.o) $(CSRCS:.c=.o)
%.o %.mod: %.f90
	$(FORT) -c  $(FFLAGS) $<
%.o: %.f
	$(FORT) -c $(F77FLAGS) $<

all: ckohn
ckohn: ${OBJS} 
	$(FORT) ${LDFLAGS} $(FFLAGS) -o $@ $^
	rm -f *.o *.mod 
clean:
	rm -f *.o *.out ckohn
