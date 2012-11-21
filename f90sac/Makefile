#-------------------------------------------------------------------------------
#  MAKEFILE
#-------------------------------------------------------------------------------  
#
#  CVS: $Revision: 1.12 $ $Date: 2008/10/17 10:48:03 $
#
#  (C) James Wookey
#  Department of Earth Sciences, University of Bristol
#  Wills Memorial Building, Queen's Road, Bristol, BR8 1RJ, UK
#  j.wookey@bristol.ac.uk
#
#-------------------------------------------------------------------------------
#
#   Makefile for building f90sac
#
#-------------------------------------------------------------------------------
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
#
#-------------------------------------------------------------------------------

F90SAC = f90sac.o f90sac_csubs.o

#-------------------------------------------------------------------------------
# Uncomment to use Intel Compilers (v8+) 
#FC = ifort
#CC = icc
#CFLAGS = -O2
#FFLAGS = -O2 -assume byterecl
#-------------------------------------------------------------------------------
# Uncomment to use Gnu Compilers 
FC = gfortran
CC = gcc
CFLAGS = -O2 # -m32
# AJN addition: In order to link to the tpttsub.o, we need to stick to 32-bit
FFLAGS = -O2 # -m32
#-------------------------------------------------------------------------------
# Uncomment to use g95
#FC = /opt/local/bin/g95
#CC = gcc
#CFLAGS = -O2
#FFLAGS = -O2 -fno-second-underscore
#-------------------------------------------------------------------------------
## Uncomment for sac/sac2000 (binaries sac files are native endian)
#F90SAC_FLAGS = 
#-------------------------------------------------------------------------------
## Uncomment for MacSAC (binary sac files are always big-endian)
F90SAC_FLAGS = -DFORCE_BIGENDIAN_SACFILES
#-------------------------------------------------------------------------------
## Uncomment to use pure Fortran IO (significantly slower)
#F90SAC_FLAGS += -DDISABLE_C_IO
#F90SAC = f90sac.o
#-------------------------------------------------------------------------------
## Options for compilation of libf90sac
LIBDIR = $(CURDIR)/../lib
#-------------------------------------------------------------------------------

all:test_f90sac

test_f90sac: ${F90SAC} test_f90sac.o 
	$(FC) -o test_f90sac ${F90SAC} test_f90sac.o 

libf90sac: f90sac.F90 f90sac_csubs.c
	$(FC) ${FFLAGS} ${F90SAC_FLAGS} -c -fpic -J../mods -o ../objs/f90sac.o f90sac.F90
	$(CC) ${CFLAGS} -c -fpic -o ../objs/f90sac_csubs.o f90sac_csubs.c
	$(FC) -o $(LIBDIR)/libf90sac.so.1 -shared -W1,-soname,libf90sac.so.1 ../objs/f90sac.o ../objs/f90sac_csubs.o
	ln -sf $(LIBDIR)/libf90sac.so.1 $(LIBDIR)/libf90sac.so

f90sac.o: f90sac.F90
	$(FC) $(FFLAGS) $(F90SAC_FLAGS) -c f90sac.F90
	
clean:
	rm -f *.o *.mod

# COMPILE INSTRUCTIONS
%.o: %.F90
	$(FC) $(FFLAGS) -c $*.F90 
%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 
%.o: %.c
	$(CC) $(CFLAGS) -c $*.c
