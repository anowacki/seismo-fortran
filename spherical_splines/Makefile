FC = gfortran
FCFLAGS = -O3
LIBS = -framework vecLib -llapack -lspherical_geometry
L = ../lib
M = ../mods
O = ../objs

default: all

all: test $(O)/spherical_splines.o

test: test_spherical_splines test_spherical_splines_AP test_AP_width2h

testSH: test_SHTOOLS

test_spherical_splines_AP: $(O)/spherical_splines.o test_spherical_splines_AP.o
	$(FC) ${FCFLAGS} -o $@ -L$(L) ${LIBS} -J$(M) $^

test_spherical_splines: $(O)/spherical_splines.o test_spherical_splines.o
	$(FC) ${FCFLAGS} -o $@ -L$(L) ${LIBS} -J$(M) $^

test_AP_width2h: $(O)/spherical_splines.o test_AP_width2h.o
	$(FC) ${FCFLAGS} -o $@ -L$(L) ${LIBS} -J$(M) $^

test_SHTOOLS: test_SHTOOLS.f90
#	$(FC) ${FCLAGS} -c -o $(O)/test_SHTOOLS.o -J$(M) $^
	ln -sf SHTOOLS/modules/shtools.mod
	$(FC) ${FCFLAGS} -o $@ -L$(L) -lspherical_geometry -J$(M) -framework vecLib \
        SHTOOLS/lib/libSHTOOLS2.6.a $^

# Compilation rules
$(O)/%.o: %.f90
	$(FC) ${FCFLAGS} ${LIBS} -J$(M) -c -o $(O)/$*.o $*.f90
%.o: %.f90
	$(FC) ${FCFLAGS} -J$(M) -c -o $*.o $*.f90
	
.PHONY: clean

clean:
	/bin/rm -f $(O)/spherical_splines.o *.o *.mod