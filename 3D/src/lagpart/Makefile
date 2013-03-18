#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

OBJS = 

CUDA_GHOSTno   =
CUDA_GHOSTyes  = -DDEF_GHOST_CUDA_

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(INCLUDE) $(@:.o=.f90)
gpdbrec_mod.o:
	$(COMPILE) gpdbrec_mod.f90
gpartcomm_mod.o:
	$(COMPILE) gpartcomm_mod.f90
gparts_mod.o:
	$(COMPILE) gparts_mod.f90
gpsplineint_mod.o:
	$(COMPILE) gpsplineint_mod.f90

edit:
	$(CPP) -DG$(PRECISION)_PRECISION $(CUDA_GHOST$(P_CUDA)) -I../fftp-3  regridutils.fpp -o \
           regridutils.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc* regridutils.f90

