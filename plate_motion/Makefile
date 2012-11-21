FC = gfortran
FCOPTS = -O3
BINDIR = ~nowacki/Applications/bin
MODDIR = ~nowacki/Applications/modules/mods
OBJDIR = ~nowacki/Applications/modules/objs
LIBDIR = /Users/nowacki/Applications/modules/lib

default: all

all: progs

test: test_plate_motion test_inout

install: progs
	/bin/cp spread apm plate_name plate_name_list $(BINDIR)/

progs: spread apm plate_name plate_name_list

test_plate_motion: test_plate_motion.o plate_motion.o
	$(FC) ${FCOPTS} -o test_plate_motion plate_motion.o test_plate_motion.o

spread: spread.o $(OBJDIR)/plate_motion.o
	$(FC) ${FCOPTS} -J $(MODDIR) -o spread $(OBJDIR)/plate_motion.o spread.o

apm: apm.o plate_motion.o
	$(FC) ${FCOPTS} -J $(MODDIR) -o apm $(OBJDIR)/plate_motion.o apm.o

plate_name: plate_name.o plate_motion.o
	$(FC) ${FCOPTS} -J $(MODDIR) -o plate_name $(OBJDIR)/plate_motion.o plate_name.o

plate_name_list: plate_name_list.o plate_motion.o
	$(FC) ${FCOPTS} -J $(MODDIR) -o plate_name_list $(OBJDIR)/plate_motion.o plate_name_list.o

test_inout: test_inout.o $(LIBDIR)/libspherical_geometry.so
	$(FC) ${FCOPTS} -J $(MODDIR) -L$(LIBDIR) -lspherical_geometry -o test_inout test_inout.o 

%.o: %.f90
	$(FC) ${FCOPTS} -J $(MODDIR) -c $*.f90

.PHONY: clean

clean:
	/bin/rm -f *.o test_plate_motion