FC = gfortran
FCOPTS = -Wall -g -ggdb -fcheck=all -O0

O = ../objs
M = ../mods
L = ../lib

LIBS = -L$L -lFFFTW -lf90sac
INCS = 

MODS = $(O)/splitwave.o

PROGS = test_splitwave test_misfit

all: ${PROGS}

default: ${PROGS}

run: ${PROGS}
	for p in ${PROGS}; do \
		echo $$p; \
		./$$p; \
	done

debug: clean ${PROGS}
	for p in ${PROGS}; do \
		ggdb $$p; \
	done

	
%: ${MODS} $(O)/%.o
	$(FC) ${FCOPTS} -o $* ${INCS} ${LIBS} -J$(M) ${MODS} $(O)/$*.o

$(O)/%.o: ../%.f90
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} ../$*.f90

$(O)/%.o: %.f90
	$(FC) ${FCOPTS} -c -o $(O)/$*.o -J$(M) ${INCS} $*.f90

.SECONDARY: ${MODS} ${PROGS:%=$(O)/%.o}

.PHONY: clean run debug

clean:
	/bin/rm -f ${PROGS} ${MODS} ${PROGS:%=$O/%.o}