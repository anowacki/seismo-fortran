# Makefile for wrapper programs for routines provided in modules

FC = gfortran
FFLAGS = -O2 -ggdb -fbounds-check -frange-check
B = bin
L = -L../lib
M = ../mods
O = obj
S = src
LIBS = -lspherical_geometry -lget_args -lglobal_1d_models -lsplitwave -lf90sac -lFFFTW
BINDIR = ~/Applications/bin

PROGS = $(B)/GPa2km \
        $(B)/Vak135 \
        $(B)/Vprem \
		$(B)/ak135 \
        $(B)/azi \
        $(B)/baz \
        $(B)/create_wave \
        $(B)/fdsplitwaveN \
		$(B)/gcp_from_point_azi \
		$(B)/gcp_from_points \
        $(B)/gcp_points \
		$(B)/km2GPa \
		$(B)/midpoint \
		$(B)/prem \
		$(B)/project_to_gcp \
		$(B)/split_misfit \
        $(B)/step

all: progs

install: progs
	cp ${PROGS} $(BINDIR)/

progs: ${PROGS}

$(B)/%: $(O)/%.o
	$(FC) ${FFLAGS} -o $(B)/$* -J$(M) ${L} ${LIBS} $(O)/$*.o

$(O)/%.o: $(S)/%.f90
	$(FC) ${FFLAGS} -c -o $(O)/$*.o -J$(M) $(S)/$*.f90

.SECONDARY: ${PROGS:$(B)/%=$(O)/%.o}

.PHONY: clean

clean:
	/bin/rm -f $(O)/*.o *.mod ${PROGS}
