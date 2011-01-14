#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

OBJS = gutils_mod.o

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(@:.o=.f90)

edit:
	$(CPP) -DG$(PRECISION)_PRECISION gutils._mod.fpp -o \
           gutils_mod.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*

