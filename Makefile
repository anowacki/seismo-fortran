# Makefile to create shared objects out of the modules in this folder

# Need the $(CURDIR) so that when link against the libraries, dyld knows where they are.
L = $(CURDIR)/lib
O = objs
M = mods
FC = gfortran
FCOPTS = -O3 -fcheck=all -fbacktrace -g -ggdb
CC = gcc
CCOPTS = -O2
LOPTS = -fpic

# Need to link FFFTW against FFTW3
FFFTWOPTS = -I/opt/local/include -L/opt/local/lib -lfftw3 -lfftw3f
SPLINEOPTS = -framework vecLib -llapack
SPLITWAVEOPTS = -L${L} -lFFFTW -lf90sac -lEmatrixUtils
F90SACOPTS = -DFORCE_BIGENDIAN_SACFILES

MODS = $(O)/constants.o \
       $(O)/density_1d.o \
       $(O)/EC_grid_assumed_int.o \
       $(O)/EC_grid.o \
	   $(O)/f90sac.o \
	   $(O)/FFFTW.o \
       $(O)/functions.o \
       $(O)/global_1d_models.o \
       $(O)/mod_raypaths.o \
       $(O)/spherical_geometry.o \
       $(O)/splitwave.o \
       $(O)/statistical.o 

all: ${MODS} \
     $(O)/anisotropy_ajn.o \
     $(O)/get_args.o \
     $(O)/plate_motion.o \
     $(O)/spherical_splines.o

progs:
	$(MAKE) -C progs

installprogs:
	$(MAKE) -C progs install

$(O)/anisotropy_ajn.o: anisotropy_ajn/anisotropy_ajn.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(M) -o $(O)/anisotropy_ajn.o anisotropy_ajn/anisotropy_ajn.f90
	$(FC) -I$(M) -o $(L)/libanisotropy_ajn.so.1 -shared -W1,-soname,libanisotropy_ajn.so.1 $(O)/anisotropy_ajn.o
	ln -sf $(L)/libanisotropy_ajn.so.1 $(L)/libanisotropy_ajn.so

$(O)/FFFTW.o: FFFTW/FFFTW.f03
	$(FC) ${FCOPTS} ${LOPTS} ${FFFTWOPTS} -c -J$(M) -o $(O)/FFFTW.o FFFTW/FFFTW.f03
	$(FC) -I$(M) ${FFFTWOPTS} -o $(L)/libFFFTW.so.1 -shared -W1,-soname,FFFTW.so.1 $(O)/FFFTW.o
	ln -sf $(L)/libFFFTW.so.1 $(L)/libFFFTW.so

$(O)/get_args.o: get_args/get_args.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(M) -o $(O)/get_args.o get_args/get_args.f90
	$(FC) -I$(M) -o $(L)/libget_args.so.1 -shared -W1,-soname,get_args.so.1 $(O)/get_args.o
	ln -sf $(L)/libget_args.so.1 $(L)/libget_args.so

$(O)/plate_motion.o: plate_motion/plate_motion.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(M) -o $(O)/plate_motion.o plate_motion/plate_motion.f90
	$(FC) -I$(M) -o $(L)/libplate_motion.so.1 -shared -W1,-soname,libplate_motion.so.1 $(O)/plate_motion.o
	ln -sf $(L)/libplate_motion.so.1 $(L)/libplate_motion.so

$(O)/spherical_splines.o: spherical_splines/spherical_splines.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(M) -o $(O)/spherical_splines.o spherical_splines/spherical_splines.f90
	$(FC) -I$(M) ${SPLINEOPTS} -o $(L)/libspherical_splines.so.1 -shared -W1,-soname,libspherical_splines.so.1 $(O)/spherical_splines.o
	ln -sf $(L)/libspherical_splines.so.1 $(L)/libspherical_splines.so

$(O)/splitwave.o: $(O)/f90sac.o $(O)/EmatrixUtils.o $(O)/FFFTW.o splitwave.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(M) -o $(O)/splitwave.o splitwave.f90
	$(FC) -I$(M) ${SPLITWAVEOPTS} -o $(L)/libsplitwave.so.1 -shared -W1,-soname,libsplitwave.so.1 $(O)/splitwave.o
	ln -sf $(L)/libsplitwave.so.1 $(L)/libsplitwave.so

$(O)/f90sac.o: $(O)/f90sac_csubs.o f90sac/f90sac.F90
	$(FC) ${FCOPTS} ${F90SACOPTS} ${LOPTS} -c -J$(M) -o $(O)/f90sac.o f90sac/f90sac.F90
	$(FC) -I$(M) ${F90SACOPTS} -o $(L)/libf90sac.so.1 -shared -W1,-soname,libf90sac.so.1 $(O)/f90sac.o $(O)/f90sac_csubs.o
	ln -sf $(L)/libf90sac.so.1 $(L)/libf90sac.so

$(O)/f90sac_csubs.o: f90sac/f90sac_csubs.c
	$(CC) ${CCOPTS} -c ${LOPTS} -o $(O)/f90sac_csubs.o f90sac/f90sac_csubs.c

$(O)/%.o: %.f90
	$(FC) ${FCOPTS} ${LOPTS} -c $*.f90 -J$(M) -o $(O)/$*.o
	$(FC) -I$(M) -o $(L)/lib$*.so.1 -shared -W1,-soname,lib$*.so.1 $(O)/$*.o
	$(MAKE) --directory=lib OBJ=$*


.PHONY: progs installprogs

clean:
	/bin/rm -f $(M)/*.mod $(O)/*.o
