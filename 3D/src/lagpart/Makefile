#****************************************************************************
# Makefile for compiling particle modules
# D. Rosenberg - 8/04/2013
#****************************************************************************

OBJS = gpdbrec_mod.o gpartcomm_mod.o gpsplineint_mod.o gparts_mod.o

CUDA_GHOSTno   =
CUDA_GHOSTyes  = -DDEF_GHOST_CUDA_

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(@:.o=.f90)

edit:
	$(CPP) -DG$(PRECISION)_PRECISION $(CUDA_GHOST$(P_CUDA)) -I../fftp-3 \
	regridutils.fpp -o regridutils.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc* regridutils.f90
