FC = gfortran
FCOPTS = -Wall -g -ggdb -O3

O = ../../objs
M = ../../mods

LIBS = -L/opt/local/lib -lfftw3 -lfftw3f
INCS = -I/opt/local/include

MODS = $(O)/FFFTW.o

PROGS = test_FFFTW

all: ${PROGS}

default: ${PROGS}

run: ${PROGS}
	for p in ${PROGS}; do \
		echo $$p; \
		time ./$$p; \
	done
	
%: ${MODS} $(O)/%.o
	$(FC) ${FCOPTS} -o $* ${INCS} ${LIBS} -J$(M) ${MODS} $(O)/$*.o

$(O)/%.o: ../%.f03
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} ../$*.f03

$(O)/%.o: ../%.f90
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} ../$*.f90

$(O)/%.o: %.f90
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} $*.f90

$(O)/%.o: %.f03
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} $*.f03

.PHONY: clean

clean:
	/bin/rm -f ${PROGS} ${MODS} ${PROGS:%=$O/%.o}