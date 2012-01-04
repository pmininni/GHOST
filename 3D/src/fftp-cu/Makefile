#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

FOBJS  = fprecision_mod.o fftp_mod.o cuda_bindings_mod.o fftp3D.o 
CUOBJS = cuda_utils.o

all: $(FOBJS) $(CUOBJS)

$(FOBJS):
	$(COMPILE) $(@:.o=.f90)
$(CUOBJS):
	$(CUCOMPILE) $(@:.o=.cu)

edit:
	$(CPP) -DG$(PRECISION)_PRECISION fprecision_mod.fpp -o \
           fprecision_mod.f90
	$(CPP) -DIKIND_=$(IKIND) -DCSIZE_=$(CSIZE) -DNSTRIP_=$(NSTRIP) \
	   -DG$(PRECISION)_PRECISION fftp_mod.fpp -o fftp_mod.f90
	$(CPP) -DG$(PRECISION)_PRECISION fftp3D.fpp -o fftp3D.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	rm -f fprecision_mod.f90 fftp_mod.f90 fftp3D.f90
