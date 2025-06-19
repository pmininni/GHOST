#****************************************************************************
# Makefile for compiling particle modules
# D. Rosenberg - 8/04/2013
#****************************************************************************

OBJS = gpdbrec_mod.o gpartcomm_mod.o gpicsplineint_mod.o gpsplineint_mod.o \
       gparts_mod.o

CUDA_GHOSTno   =
CUDA_GHOSTyes  = -DDEF_GHOST_CUDA_

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(@:.o=.f90)

edit:
	$(CPP) -D$(DO_HYBRID) gpsplineint_mod.fpp -o gpsplineint_mod.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod *.modmic work.pc* ifc* \
        gpsplineint_mod.f90
