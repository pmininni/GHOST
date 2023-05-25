#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

OBJS = binary_mod.o binary_io.o

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(@:.o=.f90)

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod *.modmic work.pc* ifc*
