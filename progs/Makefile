# Makefile for anisotropy_ajn module and test programs

FC = gfortran
FFLAGS = -g -ggdb -fbounds-check -frange-check
B = bin
L = -L../lib
M = ../mods
O = obj
S = src
LIBS = -lget_args -lspherical_geometry
BINDIR = ~/Applications/bin

PROGS = $(B)/gcp_points \
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
