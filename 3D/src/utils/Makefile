#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

OBJS = #regridutils.o gtimer_mod.o gutils_mod.o 

CUDA_GHOSTno   =
CUDA_GHOSTyes  = -DDEF_GHOST_CUDA_

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(@:.o=.f90)
regridutils.o:
	$(COMPILE) regridutils.f90
gutils_mod.o:
	$(COMPILE) gutils_mod.f90
gtimer_mod.o:
	$(COMPILE) gtimer_mod.f90
	

edit:
	$(CPP) -DG$(PRECISION)_PRECISION $(CUDA_GHOST$(P_CUDA)) -I../fftp-3  regridutils.fpp -o \
           regridutils.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc* regridutils.f90

