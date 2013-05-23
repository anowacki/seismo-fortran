FC = gfortran-mp-4.6
FCFLAGS = -g -ggdb
MODDIR = /Users/nowacki/Applications/modules/mods
LIBDIR = /Users/nowacki/Applications/modules/lib
LIBS = -lspherical_geometry

default: all

progs = test_sphere_sample \
        test_vincenty \
        test_gcp_points \
        test_step

all: ${progs}

test_sphere_sample: test_sphere_sample.o
	$(FC) ${FCFLAGS} -o test_sphere_sample -L$(LIBDIR) ${LIBS} -J$(MODDIR) test_sphere_sample.o

test_vincenty: test_vincenty.o
	$(FC) ${FCFLAGS} -o test_vincenty -L$(LIBDIR) ${LIBS} -J$(MODDIR) test_vincenty.o

test_gcp_points: test_gcp_points.o
	$(FC) ${FCFLAGS} -o $@ -L$(LIBDIR) ${LIBS} -J$(MODDIR) $^

test_step: test_step.o
	$(FC) ${FCFLAGS} -o $@ -L$(LIBDIR) ${LIBS} -J$(MODDIR) $^

%.o: %.f90
	$(FC) ${FCFLAGS} -c -L$(LIBDIR) ${LIBS} -I$(MODDIR) $*.f90

.PHONY: clean

clean:
	/bin/rm -f *.o *.mod ${progs}