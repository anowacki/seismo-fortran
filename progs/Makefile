# Makefile for wrapper programs for routines provided in modules

FC = gfortran
FFLAGS = -g -ggdb -fbounds-check -frange-check
B = bin
L = -L../lib
M = ../mods
O = obj
S = src
LIBS = -lspherical_geometry
BINDIR = ~/Applications/bin

PROGS = $(B)/azi \
        $(B)/baz \
		$(B)/gcp_from_point_azi \
		$(B)/gcp_from_points \
        $(B)/gcp_points \
		$(B)/project_to_gcp \
        $(B)/step

all: progs

install: progs
	cp ${PROGS} $(BINDIR)/

progs: ${PROGS}

$(B)/%: $(O)/%.o
	$(FC) ${FFLAGS} -o $(B)/$* -J$(M) ${L} ${LIBS} $(O)/$*.o

$(O)/%.o: $(S)/%.f90
	$(FC) ${FFLAGS} -c -o $(O)/$*.o -J$(M) $(S)/$*.f90


.PHONY: clean

clean:
	/bin/rm -f $(O)/*.o *.mod ${PROGS}
