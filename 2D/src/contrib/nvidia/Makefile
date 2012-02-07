#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

COBJS = setaffinity_for_nvidia.o

all: $(COBJS)

$(COBJS):
	$(CCOMPILE) $(@:.o=.c)

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*

