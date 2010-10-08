#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

#****************************************************************************
# Edit library paths, compiler flags, and code options in Makefile.in
#****************************************************************************
include Makefile.in

#****************************************************************************
# Don't edit below this line
#****************************************************************************
CPP         = $(CPP_$(COMPILER))
FC          = $(FC_$(COMPILER))
FFLAGS      = $(FFLAGS_$(COMPILER))
TARGET_ARCH = $(TARGET_$(COMPILER))
MPIINC      = $(MPIINC_$(MPI))
MPILIB      = $(MPILIB_$(MPI))
MPILDFLAGS  = $(MPILD_$(MPI))

GHOME       = $(PWD)
FFTW2SINGLE = s
FFTW2DOUBLE = d
FFTW2FLAGS  = -l$(FFTW2$(PRECISION))rfftw -l$(FFTW2$(PRECISION))fftw -lm
FFTW3SINGLE = f
FFTW3DOUBLE =
FFTW3FLAGS  = -lfftw3$(FFTW3$(PRECISION)) -lm 
FFTWDIR     = $($(FFTP:fftp-%=FFTW%)DIR)
FFTWLDFLAGS = $($(FFTP:fftp-%=FFTW%)FLAGS)
FFTWINC     = $(FFTWDIR)/include
FFTWLIB     = $(FFTWDIR)/lib
LDFLAGS     = $(MPILIB) -L$(FFTWLIB) $(MPILDFLAGS) $(FFTWLDFLAGS)
INCLUDE     = $(MPIINC) -I$(FFTWINC) -I$(GHOME)/$(FFTP) -I$(GHOME)/$(IOLIB) -Iinclude
CPPFLAGS    = -D$(SOLVER)_SOL
PREPROCESS  = $(CPP) $(CPPFLAGS)
COMPILE     = $(FC) $(FFLAGS) $(INCLUDE) $(TARGET_ARCH) -c

HD_OBJS        = pseudospec2D_mod.o module_dns.o pseudospec2D_hd.o
MHD_OBJS       = $(HD_OBJS) pseudospec2D_mhd.o
MHDB_OBJS      = $(MHD_OBJS)
HMHD_OBJS      = $(MHD_OBJS) pseudospec2D_25D.o
OBJS           = $($(SOLVER)_OBJS)
OBJSDIR        = pseudo
SUBDIRS        = posix mpiio fftp-2 fftp-3

all: edit $(FFTP)/fftp2D.o $(IOLIB)/binary_io.o $(OBJS) $(SOLVER) undo

main: all

clean:  undo
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	rm -f $(OBJSDIR)/pseudospec2D_mod.f90
	for dir in $(SUBDIRS); do \
            cd $$dir && $(MAKE) clean && cd .. ; \
        done

dist:
	rm -f *HD* *MHD* *HMHD* *ROTH*
	$(MAKE) clean

edit:
	$(CPP) -DN_=$(N) -DORD_=$(ORD) $(OBJSDIR)/pseudospec2D_mod.fpp -o \
	   $(OBJSDIR)/pseudospec2D_mod.f90
	cd $(FFTP) && $(MAKE) edit CPP='$(CPP)' IKIND='$(IKIND)' \
	   CSIZE='$(CSIZE)' NSTRIP='$(NSTRIP)' PRECISION='$(PRECISION)' \
           DO_HYBRID='DO_HYBRID$(P_HYBRID)'
	if [ ${UNDER} -eq 1 ]; then \
	if [ -e main2D.fpp.bak ]; then \
	 echo 'Files have already been modified!'; else \
	 perl -i.bak -pe 's/fftw_f77\(/fftw_f77${APPEND}\(/g' *.f90 *.fpp ; \
	 for item in `cat external`; do \
	 `echo 'perl -i -pe 's/$${item}[\(]/$${item}${APPEND}\(/g' *.f90'`; \
	 `echo 'perl -i -pe 's/$${item}[\(]/$${item}${APPEND}\(/g' *.fpp'`; \
	 done ; \
	 for dir in include $(OBJSDIR) $(SUBDIRS) $(EXTRASDIR); do \
	  cd $$dir ; \
	  perl -i.bak -pe 's/fftw_f77\(/fftw_f77${APPEND}\(/g' *.f90 ; \
	  for item in `cat ../external`; do \
	  `echo 'perl -i -pe 's/$${item}[\(]/$${item}${APPEND}\(/g' *.f90'`; \
	  done ; cd .. ; \
	 done ; \
	fi ; fi

undo:
	if [ -e main2D.fpp.bak ]; then \
	 for file in *.f??.bak; do \
	  mv $$file $${file%.bak}; \
	 done ; \
	 for dir in include $(OBJSDIR) $(SUBDIRS); do \
	  cd $$dir ; \
	  for file in *.f??.bak; do \
	   mv $$file $${file%.bak}; \
	  done; cd ..; \
	 done ; \
	fi

$(SOLVER):
	$(PREPROCESS) main2D.fpp -o $(SOLVER).f90
	$(COMPILE) $(SOLVER).f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(OBJS) \
           $(SOLVER).o $(LDFLAGS) -o $(SOLVER)

$(FFTP)/fftp2D.o:
	cd $(FFTP) && $(MAKE) COMPILE='$(COMPILE)'

$(IOLIB)/binary_io.o:
	cd $(IOLIB) && $(MAKE) COMPILE='$(COMPILE)' 

$(OBJS) $(TRAN_OBJS) $(STRU_OBJS):
	$(COMPILE) $(OBJSDIR)/$(@:.o=.f90)
