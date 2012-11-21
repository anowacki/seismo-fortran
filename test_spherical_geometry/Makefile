FC = gfortran
FCFLAGS = 
MODDIR = /Users/nowacki/Applications/modules/mods
LIBDIR = /Users/nowacki/Applications/modules/lib
LIBS = -lspherical_geometry

default: all

all: test_sphere_sample test_vincenty

test_sphere_sample: test_sphere_sample.o
	$(FC) ${FCFLAGS} -o test_sphere_sample -L$(LIBDIR) ${LIBS} -J$(MODDIR) test_sphere_sample.o

test_vincenty: test_vincenty.o
	$(FC) ${FCFLAGS} -o test_vincenty -L$(LIBDIR) ${LIBS} -J$(MODDIR) test_vincenty.o

%.o: %.f90
	$(FC) ${FCFLAGS} -c -L$(LIBDIR) ${LIBS} -I$(MODDIR) $*.f90

.PHONY: clean

clean:
	/bin/rm -f *.o *.mod