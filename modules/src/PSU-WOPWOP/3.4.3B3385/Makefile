#.SILENT:
# $Id: Makefile 3372 2017-08-17 12:52:59Z brentner $

# This file is created by the configure script: if you are getting an error, make
# sure you have run ./configure
include .make.config

##################################################################################
# Compilation flags for the compilers we typically use
#ABSOFT_FLAGS = -YEXT_NAMES=LCS -s -B108 -YCFRL=1 -et -Rb -Rc -Rs -m3 -Rp -cpu:host  
#INTEL_FLAGS  = -O2 -g -static -implicitnone -traceback -check all -CB -pad -fpstkchk -ftrapuv -warn
#GNU_FLAGS    = -O2 -g -static -fimplicit-none #-ftrace=full 
#pgf_FLAGS = -g
#PARALLEL_ABSOFT_FLAGS = -s  -et -Rp -Rc -Rs -m3 -O2 
#PARALLEL_INTEL_FLAGS = -openmp-report2 -opt-mem-bandwidth2 -reentrancy threaded -tcheck
#CFLAGS     = -O2 -g $(CFLAGS)
#FFLAGS =  -O2 -g $(FFLAGS)
##################################################################################


##################################################################################
# System rm command location
RM = /bin/rm
##################################################################################


##################################################################################
# Fast Fourier Transform library. Currently FFTW v3.3.3.
#FFT_DIR=${FFT_DIR} # not actually necessary
FFT="$(FFT_DIR)/lib/libfftw3.so"
#FFT_LIB= -L$(FFT_DIR)/lib/ -l"fftw3"
FFT_LIB= -lfftw3
##################################################################################



OBJS =  \
				COB.o \
				$(IO_MODULE) \
				Main.o \
				WOPWOPModelModule.o \
				build_number.o \
				constants.o \
				container.o \
				debug.o \
				four1.o \
				freqAnalysis.o \
				frequencyDomain.o \
				hunt.o \
				interpolate.o \
				loadingData.o \
				mathModule.o \
				$(MPIModule) \
				obsUtility.o \
				observer.o \
				observerContainer.o \
				patch.o \
				broadband.o \
				realfft.o \
				sigma.o \
				strings.o \
				surface.o \
				timeHistory.o \
				toUpper.o \
				wall.o \
				zbrac.o \
				zbrent.o \


######################################################
# The last known good version - used by "make release"
#
# DO NOT UPDATE WITHOUT CONSULTING FULL WOPWOP TEAM !!
#
######################################################
CURRENT_RELEASE_BUILD = 3370
CURRENT_RELEASE_LOCATION = "svn+ssh://128.118.168.9/var/svn/modules/wopwop3/releases/3_4_3/"
######################################################


# Compiler-specific targets
$(F90):
	./update_build_number.sh
	echo "Compiling with F90=$(F90) and CC=$(CC)"
	make compile	TARGET=wopwop3 \
				F77="$(F90)"  F90="$(F90)"  CXX="$(CC)" \
				CC="$(CC)"   LINKER="$(F90)" \
				FFLAGS="$(FFLAGS)"  \
				CFLAGS="$(CFLAGS)"

clean:
	$(RM) -f *.o *.mod core wopwop3

distclean:
	$(RM) -f *.o *.mod *.dat .make.config core wopwop3
	$(RM) -rf $(FFT_DIR)
	cd $(FFT_COMPILE_DIR); make distclean

clean_small:
	$(RM) -f *.o *.mod core

build:
	# Usage: make build BUILD=570 LOCATION=svn+ssh://128.118.168.9/var/svn/wopwop3/trunk/
	./make_build.sh $(BUILD) $(LOCATION)

release:
	make build BUILD=$(CURRENT_RELEASE_BUILD) LOCATION=$(CURRENT_RELEASE_LOCATION)

$(FFT):
#	$(PWD)/compile_fftw.sh "$(FFT_SUFFIX)"
	echo "If no errors follow the \"Linking wopwop3\" line, then any above errors can be safely ignored."


# The real targets: the ones that actually do things:

compile: $(OBJS) $(FFT)
	echo "Linking" $(TARGET)
	$(LINKER)  -o $(TARGET) $(OBJS) $(LIBS) $(FFT_LIB)


.SUFFIXES:.f .f90 .o .mod .c .cpp

.f.o:
	echo "Compiling" $<
	$(F77) $(FFLAGS) -c $<

.f90.o:
	echo "Compiling" $<
	$(F90) $(FFLAGS) -c $*.f90

.f90.mod:
	echo "Compiling" $<
	$(F90) $(FFLAGS) -c $*.f90

.cpp.o:
	echo "Compiling" $<
	$(CXX) $(CFLAGS) -c $<

.c.o:
	echo "Compiling" $<
	$(CC) $(CFLAGS) -c $<


Main.o         : debug.o $(MPIModule) constants.o COB.o wall.o observer.o container.o timeHistory.o observerContainer.o sigma.o strings.o
container.o    : debug.o $(MPIModule) mathModule.o $(IO_MODULE) COB.o patch.o broadband.o
patch.o        : debug.o $(MPIModule) $(IO_MODULE) freqAnalysis.o COB.o surface.o loadingData.o broadband.o
observer.o     : debug.o $(MPIModule) freqAnalysis.o  obsUtility.o $(IO_MODULE) $(ModelModule) strings.o
loadingData.o  : debug.o $(MPIModule) mathModule.o constants.o toUpper.o $(IO_MODULE) $(ModelModule) strings.o
interpolate.o  : debug.o mathModule.o constants.o                    
freqAnalysis.o : debug.o constants.o interpolate.o
COB.o          : debug.o $(MPIModule) constants.o toUpper.o interpolate.o strings.o
wall.o         : debug.o $(MPIModule) mathModule.o
Cylinder.o     : debug.o constants.o freqAnalysis.o
TrailingEdge.o : debug.o constants.o freqAnalysis.o
obsUtility.o   : debug.o mathModule.o COB.o $(IO_MODULE)
constants.o    : debug.o $(MPIModule)
mathModule.o   : debug.o $(MPIModule)
debug.o        : $(MPIModule) strings.o
observer.o : COB.o wall.o timeHistory.o frequencyDomain.o
observerContainer.o : debug.o $(MPIModule) mathModule.o constants.o toUpper.o $(IO_MODULE) interpolate.o freqAnalysis.o obsUtility.o observer.o timeHistory.o frequencyDomain.o
observer.o : debug.o $(MPIModule) mathModule.o constants.o toUpper.o $(IO_MODULE) interpolate.o freqAnalysis.o obsUtility.o timeHistory.o frequencyDomain.o
timeHistory.o : debug.o $(MPIModule) mathModule.o constants.o toUpper.o $(IO_MODULE) interpolate.o
frequencyDomain.o : debug.o $(MPIModule) mathModule.o constants.o toUpper.o $(IO_MODULE) interpolate.o freqAnalysis.o
sigma.o        : observerContainer.o
surface.o      : WOPWOPModelModule.o
broadband.o    : debug.o $(MPIModule) COB.o constants.o mathModule.o strings.o
