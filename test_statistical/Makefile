# Makefile to make tests of module statistical
MOD = ../statistical.o
PROGS = F_cum \
        f_dist \
        fact \
        beta \
        incomp_beta \
        reg_beta \
        f \
        chi2_95 \
        circ_disp \
        von_Mieses \
        small_bootstrap
FC = gfortran

all: ${PROGS}

run: ${PROGS}
	(echo "plotting f-distribution"; ./f_dist 5 3 10 | awk 'BEGIN{print "plot","'\''" "-" "'\'' w l"}{print}' | gnuplot -persist)
	(echo "12 factorial"; ./fact 12)

F_cum: statistical.o F_cum.o 
	$(FC) ${FCOPTS} -o F_cum statistical.o F_cum.o

f_dist: statistical.o f_dist.o
	$(FC) ${FCOPTS} -o f_dist statistical.o f_dist.o

f: statistical.o f.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

fact: statistical.o fact.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

beta: statistical.o beta.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

incomp_beta: statistical.o incomp_beta.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

reg_beta: statistical.o reg_beta.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

chi2_95: statistical.o chi2_95.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

circ_disp: statistical.o circ_disp.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

von_Mieses: statistical.o von_Mieses.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

small_bootstrap: statistical.o small_bootstrap.o
	$(FC) ${FCOPTS} -o $@ statistical.o $@.o

# Make object files
%.o: %.f90
	$(FC) ${FCOPTS} -c $*.f90

# Make module
statistical.o: ../statistical.f90
	$(FC) ${FCOPTS} -c ../statistical.f90

.PHONY: clean

clean:
	/bin/rm *.o *.mod
