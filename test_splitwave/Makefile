FC = gfortran
FCOPTS = -Wall -g -ggdb -fcheck=all

O = ../objs
M = ../mods
L = ../lib

LIBS = -L$L -lFFFTW -lf90sac
INCS = 

MODS = $(O)/splitwave.o

PROGS = test_splitwave

all: ${PROGS}

default: ${PROGS}

run: ${PROGS}
	for p in ${PROGS}; do \
		echo $$p; \
		time ./$$p; \
	done
	
%: ${MODS} $(O)/%.o
	$(FC) ${FCOPTS} -o $* ${INCS} ${LIBS} -J$(M) ${MODS} $(O)/$*.o

$(O)/%.o: ../%.f90
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} ../$*.f90

$(O)/%.o: %.f90
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} $*.f90

.PHONY: clean

clean:
	/bin/rm -f ${PROGS} ${MODS} ${PROGS:%=$O/%.o}