# Makefile to make tests of module statistical
MOD = ../objs/EC_grid.o
M = ../mods
PROGS = test_EC_grid_dump_file_bin
FC = gfortran

all: ${PROGS}

run: ${PROGS}
	(echo "plotting f-distribution"; ./f_dist 5 3 10 | awk 'BEGIN{print "plot","'\''" "-" "'\'' w l"}{print}' | gnuplot -persist)
	(echo "12 factorial"; ./fact 12)

test_EC_grid_dump_file_bin: $(MOD) test_EC_grid_dump_file_bin.o
	$(FC) ${FCOPTS} -o $@ -J $M $(MOD) $@.o

# Make object files
%.o: %.f90
	$(FC) ${FCOPTS} -c -J $M $*.f90

.PHONY: clean

clean:
	/bin/rm *.o *.mod
