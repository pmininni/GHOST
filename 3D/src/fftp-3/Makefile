#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

OBJS = fftp_mod.o fftp.o fftp3D.o

all: $(OBJS)

$(OBJS):
	$(COMPILE) $(@:.o=.f90)

edit:
	$(CPP) -DIKIND_=$(IKIND) -DCSIZE_=$(CSIZE) -DNSTRIP_=$(NSTRIP) \
	   fftp_mod.fpp -o fftp_mod.f90

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	rm -f fftp_mod.f90
