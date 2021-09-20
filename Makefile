# Makefile to create shared objects out of the modules in this folder

# Need the $(CURDIR) so that when link against the libraries, dyld knows where they are.
L = $(CURDIR)/lib
O = objs
M = mods
FC = gfortran
PIC = -fpic
FCOPTS = -O3 -fcheck=all -fbacktrace -g -ggdb
FCOPTS += $(PIC)
CC = gcc
CCOPTS = -O2
LOPTS = $(PIC)

AR = ar
AROPTS = cru
RANLIB = ranlib

# Get host for machine-specific build instructions
host = $(shell hostname | awk '{if     (/typhon/)  print "typhon"; \
                               else if (/eslogin/) print "archer"; \
			       else if (/uan/)     print "archer2"; \
                               else                print "other"}')

# Need to link FFFTW against FFTW3
# Flags which work using MacPorts on OS X:
ifeq ($(host),other)
    FFFTWOPTS = -I/opt/local/include -L/opt/local/lib -lfftw3 -lfftw3f
    SPLINEOPTS = -framework Accelerate -llapack
endif
# Flags which work on Typhon:
ifeq ($(host),typhon)
    FFFTWOPTS = -I/share/apps/local/include -L/share/apps/local/lib -lfftw3 -lfftw3f
    SPLINEOPTS = -L/share/apps/local/lib -llapack
endif
# Flags which work on ARCHER (must `module load fftw`)
ifeq ($(host),archer)
    FFFTWOPTS = -I${FFTW_INC} -L${FFTW_DIR} -lfftw3 -lfftw3f
endif
ifeq ($(host),archer2)
    FFFTWOPTS = -I${FFTW_INC} -L${FFTW_DIR} -lfftw3 -lfftw3f
endif

SPLITWAVEOPTS = -L${L} -lFFFTW -lf90sac -lanisotropy_ajn
F90SACOPTS = -DFORCE_BIGENDIAN_SACFILES
F90SACLIBS = 
# Options if using the f90sac filtering routines: edit these to point to libxapiir
XAPIIR_DEFS = -DUSE_XAPIIR
XAPIIR_LIBS = -L/usr/local/lib -lxapiir
# Options for sphere_tesselate
TESSOPTS = -L$(L) -lspherical_geometry

MODS = $(O)/constants.o \
       $(O)/density_1d.o \
       $(O)/EC_grid_assumed_int.o \
       $(O)/EC_grid.o \
       $(O)/f90sac.o \
       $(O)/FFFTW.o \
       $(O)/global_1d_models.o \
       $(O)/mod_raypaths.o \
       $(O)/moment_tensor.o \
       $(O)/sphere_tesselate.o \
       $(O)/spherical_geometry.o \
       $(O)/splitwave.o \
       $(O)/statistical.o \
       $(O)/timing.o

all: ${MODS} \
     $(O)/anisotropy_ajn.o \
     $(O)/get_args.o \
     $(O)/plate_motion.o \
     $(O)/spherical_splines.o

progs:
	$(MAKE) -C progs

installprogs:
	$(MAKE) -C progs install

# Automatic macro which defines the name of the module from the object file
nm = $(patsubst $(O)/%,%,$(patsubst %.o,%,$@))

# Add rules to build filtering routines into f90sac if we select the target 'xapiir'
xapiir filt: F90SACOPTS += $(XAPIIR_DEFS)
xapiir filt: F90SACLIBS += $(XAPIIR_LIBS)
xapiir filt: all

$(O)/f90sac.o: f90sac.F90
	$(FC) $(FCOPTS) $(F90SACOPTS) -c -J$(M) -o $(O)/$(nm).o $^
	$(FC) -I$(M) $(F90SACLIBS) -o $(L)/lib$(nm).so.1 -shared $(O)/$(nm).o
	ln -sf $(L)/lib$(nm).so.1 $(L)/lib$(nm).so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o
	$(RANLIB) $(L)/lib$(nm).a

$(O)/anisotropy_ajn.o: anisotropy_ajn/anisotropy_ajn.f90
	$(FC) ${FCOPTS} -c -J$(M) -o $(O)/$(nm).o $^
	$(FC) -I$(M) -o $(L)/lib$(nm).so.1 -shared $(O)/$(nm).o
	ln -sf $(L)/lib$(nm).so.1 $(L)/lib$(nm).so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o
	$(RANLIB) $(L)/lib$(nm).a

$(O)/FFFTW.o: FFFTW/FFFTW.f03
	$(FC) ${FCOPTS} ${FFFTWOPTS} -c -J$(M) -o $(O)/FFFTW.o FFFTW/FFFTW.f03
	$(FC) -I$(M) ${FFFTWOPTS} -o $(L)/libFFFTW.so.1 -shared $(O)/FFFTW.o
	ln -sf $(L)/libFFFTW.so.1 $(L)/libFFFTW.so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o
	$(RANLIB) $(L)/lib$(nm).a

$(O)/get_args.o: get_args/get_args.f90
	$(FC) ${FCOPTS} -c -J$(M) -o $(O)/get_args.o get_args/get_args.f90
	$(FC) -I$(M) -o $(L)/libget_args.so.1 -shared $(O)/get_args.o
	ln -sf $(L)/libget_args.so.1 $(L)/libget_args.so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o
	$(RANLIB) $(L)/lib$(nm).a

$(O)/plate_motion.o: plate_motion/plate_motion.f90
	$(FC) ${FCOPTS} -c -J$(M) -o $(O)/plate_motion.o plate_motion/plate_motion.f90
	$(FC) -I$(M) -o $(L)/libplate_motion.so.1 -shared $(O)/plate_motion.o
	ln -sf $(L)/libplate_motion.so.1 $(L)/libplate_motion.so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o
	$(RANLIB) $(L)/lib$(nm).a

$(O)/sphere_tesselate.o: $(O)/spherical_geometry.o sphere_tesselate.f90
	$(FC) ${FCOPTS} -c -J$(M) -o $(O)/sphere_tesselate.o sphere_tesselate.f90
	$(FC) -I$(M) ${TESSOPTS} -o $(L)/libsphere_tesselate.so.1 -shared $(O)/sphere_tesselate.o
	ln -sf $(L)/libsphere_tesselate.so.1 $(L)/libsphere_tesselate.so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o $(O)/spherical_geometry.o
	$(RANLIB) $(L)/lib$(nm).a
	@[ "${SEISMO_FORTRAN_DATA}" ] || { \
		echo ""; echo "==== NOTICE ===="; \
	    echo "sphere_tesselate: export the following environment variable for caching to work:"; \
		echo "    SEISMO_FORTRAN_DATA=$(CURDIR)/data"; \
		echo "================"; }

$(O)/spherical_splines.o: spherical_splines/spherical_splines.f90
	$(FC) ${FCOPTS} -c -J$(M) -o $(O)/spherical_splines.o spherical_splines/spherical_splines.f90
	$(FC) -I$(M) ${SPLINEOPTS} -o $(L)/libspherical_splines.so.1 -shared $(O)/spherical_splines.o
	ln -sf $(L)/libspherical_splines.so.1 $(L)/libspherical_splines.so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o
	$(RANLIB) $(L)/lib$(nm).a

$(O)/splitwave.o: $(O)/f90sac.o $(O)/anisotropy_ajn.o $(O)/FFFTW.o splitwave.f90
	$(FC) ${FCOPTS} -c -J$(M) -o $(O)/splitwave.o splitwave.f90
	$(FC) -I$(M) ${SPLITWAVEOPTS} -o $(L)/libsplitwave.so.1 -shared $(O)/splitwave.o
	ln -sf $(L)/libsplitwave.so.1 $(L)/libsplitwave.so
	rm -f $(O)/lib$(nm).a
	$(AR) $(AROPTS) $(L)/lib$(nm).a $(O)/$(nm).o $(O)/f90sac.o $(O)/anisotropy_ajn.o $(O)/FFFTW.o
	$(RANLIB) $(L)/lib$(nm).a

# General compilation rules
$(O)/%.o: %.f90
	$(FC) ${FCOPTS} -c $*.f90 -J$(M) -o $(O)/$*.o
	$(FC) -I$(M) -o $(L)/lib$*.so.1 -shared $(O)/$*.o
	ln -sf $(L)/lib$*.so.1 $(L)/lib$*.so
	rm -f $(L)/lib$*.a
	$(AR) ${AROPTS} $(L)/lib$*.a $(O)/$*.o
	$(RANLIB) $(L)/lib$*.a

$(O)/%.o: %.F90
	$(FC) ${FCOPTS} $(DEFINES) -c $*.F90 -J$(M) -o $(O)/$*.o
	$(FC) -I$(M) -o $(L)/lib$*.so.1 -shared $(O)/$*.o
	ln -sf $(L)/lib$*.so.1 $(L)/lib$*.so
	rm -f $(L)/lib$*.a
	$(AR) ${AROPTS} $(L)/lib$*.a $(O)/$*.o
	$(RANLIB) $(L)/lib$*.a

.PHONY: progs installprogs

clean:
	/bin/rm -f $(M)/*.mod $(O)/*.o
