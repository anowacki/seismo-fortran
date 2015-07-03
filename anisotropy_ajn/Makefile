# Makefile for anisotropy_ajn module and test programs

FC = gfortran
FFLAGS = -O0 -g -ggdb -fbounds-check -frange-check -ffpe-trap=zero,overflow,invalid

all: anisotropy_ajn.o

anisotropy_ajn.o: anisotropy_ajn.f90
	$(FC) $(FFLAGS) -c anisotropy_ajn.f90

test: test_CIJ_tandon_and_weng \
      test_CIJ_brow_chev \
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

%: anisotropy_ajn.o %.o
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $^

test_CIJ_tandon_and_weng: anisotropy_ajn.o test_CIJ_tandon_and_weng.o
	$(FC) ${FFLAGS} -o test_CIJ_tandon_and_weng test_CIJ_tandon_and_weng.o \
                       anisotropy_ajn.o

test_CIJ_tandon_and_weng.o: test_CIJ_tandon_and_weng.f90
	$(FC) ${FFLAGS} -c test_CIJ_tandon_and_weng.f90


test_CIJ_brow_chev: anisotropy_ajn.o test_CIJ_brow_chev.o
	$(FC) ${FFLAGS} -o test_CIJ_brow_chev anisotropy_ajn.o test_CIJ_brow_chev.o

test_CIJ_brow_chev.o : test_CIJ_brow_chev.f90
	$(FC) ${FFLAGS} -c test_CIJ_brow_chev.f90

test_axial_average: anisotropy_ajn.o test_axial_average.o
	$(FC) ${FFLAGS} -o $@ $^

test_axial_average.o: test_axial_average.f90
	$(FC) ${FFLAGS} -c $^

test_Reuss_av: anisotropy_ajn.o test_Reuss_av.o
	$(FC) ${FFLAGS} -o $@ $^

test_Reuss_av.o: test_Reuss_av.f90
	$(FC) ${FFLAGS} -c $^

test_flip: anisotropy_ajn.o test_flip.o
	$(FC) ${FFLAGS} -o $@ $^

test_flip.o: test_flip.f90
	$(FC) ${FFLAGS} -c $^

test_global_VTI: anisotropy_ajn.o test_global_VTI.o
	$(FC) ${FFLAGS} -o $@ $^

test_global_VTI.o: test_global_VTI.f90
	$(FC) ${FFLAGS} -c $^

test_panning_VTI: anisotropy_ajn.o test_panning_VTI.o
	$(FC) ${FFLAGS} -o $@ $^

test_panning_VTI.o: test_panning_VTI.f90
	$(FC) ${FFLAGS} -c $^

test_Voigt_av: anisotropy_ajn.o test_Voigt_av.o
	$(FC) ${FFLAGS} -o $@ $^

test_Voigt_av.o: test_Voigt_av.f90
	$(FC) ${FFLAGS} -c $^

test_symm: anisotropy_ajn.o test_symm.o
	$(FC) ${FFLAGS} -o $@ $^

test_symm.o: test_symm.f90
	$(FC) ${FFLAGS} -c $^

test_disp: anisotropy_ajn.o test_disp.o
	$(FC) ${FFLAGS} -o $@ $^

test_disp.o: test_disp.f90
	$(FC) ${FFLAGS} -c $^

test_tensor_rot: anisotropy_ajn.o test_tensor_rot.o
	$(FC) ${FFLAGS} -o $@ $^

test_tensor_rot.o: test_tensor_rot.f90
	$(FC) ${FFLAGS} -c $^

test_transform_M: anisotropy_ajn.o test_transform_M.o
	$(FC) ${FFLAGS} -o $@ $^

test_transform_M.o: test_transform_M.f90
	$(FC) ${FFLAGS} -c $^

.PHONY: clean

clean:
	/bin/rm -f *.o