# Makefile to create shared objects out of the modules in this folder

# Need the $(CURDIR) so that when link against the libraries, dyld knows where they are.
LIBDIR = $(CURDIR)/lib
OBJDIR = objs
MODSDIR = mods
FC = gfortran
FCOPTS = -O3 -fbounds-check -g -ggdb
CC = gcc
CCOPTS = -O2
LOPTS = -fpic

# Need to link FFFTW against FFTW3
FFFTWOPTS = -I/opt/local/include -L/opt/local/lib -lfftw3 -lfftw3f
SPLINEOPTS = -framework vecLib -llapack
SPLITWAVEOPTS = -L${LIBDIR} -lFFFTW -lf90sac
F90SACOPTS = -DFORCE_BIGENDIAN_SACFILES

MODS = $(OBJDIR)/constants.o \
       $(OBJDIR)/density_1d.o \
       $(OBJDIR)/EC_grid_assumed_int.o \
       $(OBJDIR)/EC_grid.o \
	   $(OBJDIR)/FFFTW.o \
       $(OBJDIR)/functions.o \
	   $(OBJDIR)/f90sac.o \
       $(OBJDIR)/EmatrixUtils.o \
       $(OBJDIR)/global_1d_models.o \
       $(OBJDIR)/mod_raypaths.o \
       $(OBJDIR)/spherical_geometry.o \
       $(OBJDIR)/splitwave.o \
       $(OBJDIR)/statistical.o 

all: ${MODS} \
     $(OBJDIR)/anisotropy_ajn.o \
     $(OBJDIR)/get_args.o \
     $(OBJDIR)/plate_motion.o \
     $(OBJDIR)/spherical_splines.o

progs:
	$(MAKE) -C progs

installprogs:
	$(MAKE) -C progs install

$(OBJDIR)/anisotropy_ajn.o: anisotropy_ajn/anisotropy_ajn.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/anisotropy_ajn.o anisotropy_ajn/anisotropy_ajn.f90
	$(FC) -I$(MODSDIR) -o $(LIBDIR)/libanisotropy_ajn.so.1 -shared -W1,-soname,libanisotropy_ajn.so.1 $(OBJDIR)/anisotropy_ajn.o
	ln -sf $(LIBDIR)/libanisotropy_ajn.so.1 $(LIBDIR)/libanisotropy_ajn.so

$(OBJDIR)/FFFTW.o: FFFTW/FFFTW.f03
	$(FC) ${FCOPTS} ${LOPTS} ${FFFTWOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/FFFTW.o FFFTW/FFFTW.f03
	$(FC) -I$(MODSDIR) ${FFFTWOPTS} -o $(LIBDIR)/libFFFTW.so.1 -shared -W1,-soname,FFFTW.so.1 $(OBJDIR)/FFFTW.o
	ln -sf $(LIBDIR)/libFFFTW.so.1 $(LIBDIR)/libFFFTW.so

$(OBJDIR)/get_args.o: get_args/get_args.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/get_args.o get_args/get_args.f90
	$(FC) -I$(MODSDIR) -o $(LIBDIR)/libget_args.so.1 -shared -W1,-soname,get_args.so.1 $(OBJDIR)/get_args.o
	ln -sf $(LIBDIR)/libget_args.so.1 $(LIBDIR)/libget_args.so

$(OBJDIR)/plate_motion.o: plate_motion/plate_motion.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/plate_motion.o plate_motion/plate_motion.f90
	$(FC) -I$(MODSDIR) -o $(LIBDIR)/libplate_motion.so.1 -shared -W1,-soname,libplate_motion.so.1 $(OBJDIR)/plate_motion.o
	ln -sf $(LIBDIR)/libplate_motion.so.1 $(LIBDIR)/libplate_motion.so

$(OBJDIR)/spherical_splines.o: spherical_splines/spherical_splines.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/spherical_splines.o spherical_splines/spherical_splines.f90
	$(FC) -I$(MODSDIR) ${SPLINEOPTS} -o $(LIBDIR)/libspherical_splines.so.1 -shared -W1,-soname,libspherical_splines.so.1 $(OBJDIR)/spherical_splines.o
	ln -sf $(LIBDIR)/libspherical_splines.so.1 $(LIBDIR)/libspherical_splines.so

$(OBJDIR)/splitwave.o: $(OBJDIR)/f90sac.o $(OBJDIR)/FFFTW.o splitwave.f90
	$(FC) ${FCOPTS} ${LOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/splitwave.o splitwave.f90
	$(FC) -I$(MODSDIR) ${SPLITWAVEOPTS} -o $(LIBDIR)/libsplitwave.so.1 -shared -W1,-soname,libsplitwave.so.1 $(OBJDIR)/splitwave.o
	ln -sf $(LIBDIR)/libsplitwave.so.1 $(LIBDIR)/libsplitwave.so

$(OBJDIR)/f90sac.o: $(OBJDIR)/f90sac_csubs.o f90sac/f90sac.F90
	$(FC) ${FCOPTS} ${F90SACOPTS} ${LOPTS} -c -J$(MODSDIR) -o $(OBJDIR)/f90sac.o f90sac/f90sac.F90
	$(FC) -I$(MODSDIR) ${F90SACOPTS} -o $(LIBDIR)/libf90sac.so.1 -shared -W1,-soname,libf90sac.so.1 $(OBJDIR)/f90sac.o $(OBJDIR)/f90sac_csubs.o
	ln -sf $(LIBDIR)/libf90sac.so.1 $(LIBDIR)/libf90sac.so

$(OBJDIR)/f90sac_csubs.o: f90sac/f90sac_csubs.c
	$(CC) ${CCOPTS} -c ${LOPTS} -o $(OBJDIR)/f90sac_csubs.o f90sac/f90sac_csubs.c

$(OBJDIR)/%.o: %.f90
	$(FC) ${FCOPTS} ${LOPTS} -c $*.f90 -J$(MODSDIR) -o $(OBJDIR)/$*.o
	$(FC) -I$(MODSDIR) -o $(LIBDIR)/lib$*.so.1 -shared -W1,-soname,lib$*.so.1 $(OBJDIR)/$*.o
	$(MAKE) --directory=lib OBJ=$*


.PHONY: progs installprogs

clean:
	/bin/rm -f $(MODSDIR)/*.mod $(OBJDIR)/*.o