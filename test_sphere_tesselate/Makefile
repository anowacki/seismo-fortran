# Makefile to make tests of module statistical

PROGS = dump_ico \
        generate_list_neighbours \
        save_and_load_list_neighbours \
        load_neighbours_from_cache
FC = gfortran
FCOPTS = -O0 -g -ggdb
O = ../objs
M = ../mods
LIBDIR = ../lib
LIBS = -lspherical_geometry

all: ${PROGS}


dump_ico: $(O)/sphere_tesselate.o dump_ico.o 
	$(FC) ${FCOPTS} -o $@ -L$(LIBDIR) ${LIBS} $(O)/sphere_tesselate.o $@.o 

generate_list_neighbours: $(O)/sphere_tesselate.o generate_list_neighbours.o
	$(FC) ${FCOPTS} -o $@ -L$(LIBDIR) ${LIBS} $(O)/sphere_tesselate.o $@.o 
	
save_and_load_list_neighbours: $(O)/sphere_tesselate.o save_and_load_list_neighbours.o
	$(FC) ${FCOPTS} -o $@ -L$(LIBDIR) ${LIBS} $(O)/sphere_tesselate.o $@.o

load_neighbours_from_cache: $(O)/sphere_tesselate.o load_neighbours_from_cache.o
	$(FC) ${FCOPTS} -o $@ -L$(LIBDIR) ${LIBS} $(O)/sphere_tesselate.o $@.o


# Make object files
%.o: %.f90
	$(FC) ${FCOPTS} -c -J$(M) $*.f90


.PHONY: clean

clean:
	/bin/rm *.o *.mod
