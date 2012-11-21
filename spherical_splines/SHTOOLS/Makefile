###################################################################################
#
#	make all
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". The default is to use "f95".
#
#	make all2
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". LAPACK subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK	libraries with conflicting underscore 
#		conventions. The default is to use "f95".
#
#	make all3
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". FFTW subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK	libraries with conflicting underscore 
#		conventions. The default is to use "f95".
#
#	make install
#		Move the entire directory to usr/local/SHTOOLS2.5. This
#		generally requires root priviledges. In any case, you 
#		should first compile the archive using "make all".
#
#	make remove-doc
#		Remove the man and html-man pages.
#
#	make clean
#		Remove the lib and module files.
#
#	make doc
#		Create the man and html-man pages from input POD files.
#		These are PRE-MADE in the distribution, so it shouldn't
#		be necessary to recreate these unless there is some kind
#		of problem.
#
#
#	Written by Mark Wieczorek (July 2007).
#
#
#	TO DO
#		
#		Write a configure script to determine the underscore conventions
#		used with the FFTW and LAPACK libraries and then determine the necessary
#		F95 compiler flags. Also determine if these libraries are 32 or 64 bit. 
#		If you know how to do this, please let me know!
#
#####################################################################################

SHELL=/bin/tcsh

MAKE = make
DOCDIR = src/doc
SRCDIR = src
F95 = gfortran

.PHONY: all all2 all3 install doc remove-doc clean getflags
	
all: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!

all2: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all2 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!

all3: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all3 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!


getflags:

ifeq ($(F95),f95)
# Default Absoft Pro Fortran flags
# F95FLAGS ?= -O3 -YEXT_NAMES=LCS
#F95FLAGS ?= -O3 -YEXT_NAMES=LCS -YEXT_SFX=_
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_
endif

ifeq ($(F95),g95)
# Default g95 flags
#F95FLAGS ?= -O3 -fno-underscoring
F95FLAGS ?= -O3 -fno-second-underscore
#F95FLAGS ?= -m64 -O3 -fno-second-underscore 
endif

ifeq ($(F95),gfortran)
# Default gfortran flags
F95FLAGS ?= -O3
#F95FLAGS ?= -m64 -O3
endif

ifeq ($(F95),ifort)
# Default intel fortran flags
F95FLAGS ?= -free -O3 -Tf
#F95FLAGS ?= -m64 -free -O3 -Tf
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -O3
# F95FLAGS = -m64 -O3
endif


install:
	@cp -R ../SHTOOLS /usr/local/SHTOOLS2.6/
	@echo
	@echo Entire directory moved to /usr/local/SHTOOLS2.6
	

doc: 
	$(MAKE) -C $(DOCDIR) -f Makefile
	@echo Documentation successfully created


remove-doc:
	@rm -f man/man1/*.1
	@rm -f www/man/*.html
	@echo
	@echo Removed man and html-man files
	

clean:
	$(MAKE) -C $(SRCDIR) -f Makefile clean
	@echo
	@echo Removed lib, module, and object files
	