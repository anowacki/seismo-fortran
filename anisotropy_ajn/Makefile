# Makefile for anisotropy_ajn module and test programs

FC = gfortran
FFLAGS = -O0 -g -ggdb -fbounds-check -frange-check -ffpe-trap=zero,overflow,invalid

MODULE = anisotropy_ajn.o

TESTPROGS = \
      test_CIJ_tandon_and_weng \
      test_CIJ_brow_chev \
      test_read_string \
      test_Reuss_av \
      test_Voigt_av \
      test_axial_average \
      test_flip \
      test_global_VTI \
      test_panning_VTI \
      test_rot_euler \
      test_symm \
      test_disp \
      test_tensor_rot \
      test_transform_M

all: ${MODULE}

test: ${TESTPROGS}

${TESTPROGS}: ${MODULE}

anisotropy_ajn.o: anisotropy_ajn.f90
	$(FC) $(FFLAGS) -c anisotropy_ajn.f90


%: anisotropy_ajn.o %.o
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

.PHONY: clean

clean:
	/bin/rm -f *.o ${TESTPROGS}