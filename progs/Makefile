# Makefile for wrapper programs for routines provided in modules

FC = gfortran
FFLAGS = -O0 -g -ggdb -fcheck=all -fbacktrace
B = bin
L = -L../lib
M = ../mods
O = obj
S = src
LIBS = \
    -lanisotropy_ajn \
    -lf90sac \
    -lFFFTW \
    -lget_args \
    -lglobal_1d_models \
    -lmoment_tensor \
    -lsphere_tesselate \
    -lspherical_geometry \
    -lsplitwave \
    -lstatistical
BINDIR = ~/Applications/bin

# Extra libraries needed on Typhon
# L += -L/share/apps/local/lib
# LIBS += \
#     -lfftw3 \
#     -lfftw3f

PROGS = $(B)/CIJ_phase_vels \
        $(B)/CIJ_pitl \
        $(B)/CIJ_rotV \
        $(B)/CIJ_rot_azinc \
        $(B)/GPa2km \
        $(B)/MT2rp \
        $(B)/Vak135 \
        $(B)/Vprem \
        $(B)/ak135 \
        $(B)/azi \
        $(B)/baz \
        $(B)/circ_correl \
        $(B)/circ_stats \
        $(B)/circ_test_orient \
        $(B)/create_wave \
        $(B)/delta \
        $(B)/fdsplitwaveN \
        $(B)/fp2ap \
        $(B)/gcp_from_point_azi \
        $(B)/gcp_from_points \
        $(B)/gcp_points \
        $(B)/lonlat2xyz \
        $(B)/km2GPa \
        $(B)/midpoint \
        $(B)/prem \
        $(B)/project_to_gcp \
        $(B)/sdr2MT \
        $(B)/sph \
        $(B)/sph_rand \
        $(B)/sph_tess \
        $(B)/sph_tess_save \
        $(B)/split_misfit \
        $(B)/split_misfit_ecs \
        $(B)/splitwave \
        $(B)/splitwaveN \
        $(B)/step \
        $(B)/triangle_area

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
