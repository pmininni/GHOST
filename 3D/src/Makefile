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

DO_HYBRIDyes= $($(COMPILER)_OMP)
DO_HYBRIDno =

GHOME       = $(PWD)
FFTW2SINGLE = s
FFTW2DOUBLE = d
FFTW2HYByes = -l$(FFTW2$(PRECISION))rfftw_threads -l$(FFTW2$(PRECISION))fftw_threads
FFTW2HYBno  =
FFTW2FLAGS  = $(FFTW2HYB$(P_HYBRID)) -l$(FFTW2$(PRECISION))rfftw -l$(FFTW2$(PRECISION))fftw -lm
FFTW3SINGLE = f
FFTW3DOUBLE =
FFTW3HYByes = -lfftw3$(FFTW3$(PRECISION))_threads
FFTW3HYBno  =
FFTW3FLAGS  = $(FFTW3HYB$(P_HYBRID)) -lfftw3$(FFTW3$(PRECISION)) -lm 
FFTWDIR     = $($(FFTP:fftp-%=FFTW%)DIR)
FFTWLDFLAGS = $($(FFTP:fftp-%=FFTW%)FLAGS) $(DO_HYBRID$(P_HYBRID))
FFTWINC     = $(FFTWDIR)/include
FFTWLIB     = $(FFTWDIR)/lib
LDFLAGS     = $(MPILIB) -L$(FFTWLIB) $(MPILDFLAGS) $(FFTWLDFLAGS)
INCLUDE     = $(MPIINC) -I$(FFTWINC) -I$(GHOME)/$(FFTP) -I$(GHOME)/$(IOLIB) -I$(GHOME)/utils -Iinclude
CPPFLAGS    = -D$(SOLVER)_SOL
PREPROCESS  = $(CPP) $(CPPFLAGS)
COMPILE     = $(FC) $(FFLAGS) $(DO_HYBRID$(P_HYBRID)) $(INCLUDE) $(TARGET_ARCH) -c

HD_OBJS        = pseudospec3D_mod.o module_dns.o pseudospec3D_hd.o  
PHD_OBJS       = $(HD_OBJS) pseudospec3D_phd.o
MHD_OBJS       = $(HD_OBJS) pseudospec3D_mhd.o
MHDB_OBJS      = $(HD_OBJS) pseudospec3D_mhd.o pseudospec3D_rot.o
HMHD_OBJS      = $(MHD_OBJS)
ROTH_OBJS      = $(HD_OBJS) pseudospec3D_rot.o
PROTH_OBJS     = $(PHD_OBJS) pseudospec3D_rot.o pseudospec3D_prot.o
ALPHA_OBJS     = module_lles.o pseudospec3D_filt.o
CAHD_OBJS      = $(HD_OBJS) $(ALPHA_OBJS) pseudospec3D_lles.o pseudospec3D_cahd.o
LAHD_OBJS      = $(HD_OBJS) $(ALPHA_OBJS) pseudospec3D_lles.o pseudospec3D_lahd.o
LHD_OBJS       = $(HD_OBJS) $(ALPHA_OBJS) pseudospec3D_lhd.o
LAMHD_OBJS     = $(MHD_OBJS) $(ALPHA_OBJS) pseudospec3D_lles.o pseudospec3D_lahd.o
EDQNM_OBJS     = module_edqnm.o pseudospec3D_edqnm.o
EDQNMHD_OBJS   = $(HD_OBJS) $(EDQNM_OBJS)
EDQNMROTH_OBJS = $(ROTH_OBJS) $(EDQNM_OBJS)
TRAN_OBJS      = pseudospec3D_trn.o
STRU_OBJS      = pseudospec3D_shf.o
OBJS           = $($(SOLVER)_OBJS)
UTILS          = utils
OBJSDIR        = pseudo
SUBDIRS        = posix mpiio fftp-2 fftp-3 utils
EXTRASDIR      = analysis

all: edit $(FFTP)/fftp3D.o $(IOLIB)/binary_io.o $(UTILS)/gutils_mod.o $(OBJS) $(SOLVER) undo

main: all

trans triad: edit $(FFTP)/fftp3D.o $(IOLIB)/binary_io.o $(OBJS) $(TRAN_OBJS)
	$(MAKE) $@$(SOLVER)
	$(MAKE) undo

struc:	edit $(FFTP)/fftp3D.o $(IOLIB)/binary_io.o $(OBJS) $(STRU_OBJS)
	$(MAKE) $@$(SOLVER)
	$(MAKE) undo

trunc:	edit uedit $(FFTP)/fftp3D.o $(IOLIB)/binary_io.o $(OBJS) $(UTILS)/gutils_mod.o $(UTILS)/regridutils.o
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I./fftp-3 $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(UTILS)/*.o $(OBJS) \
           $@.o $(LDFLAGS) -o ../bin/$@
	$(MAKE) undo

boots: edit uedit $(FFTP)/fftp3D.o $(IOLIB)/binary_io.o $(OBJS) $(UTILS)/gutils_mod.o $(UTILS)/regridutils.o
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I./fftp-3 $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(UTILS)/*.o $(OBJS) \
           $@.o $(LDFLAGS) -o ../bin/$@
	$(MAKE) undo

clean:  undo
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	rm -f $(OBJSDIR)/pseudospec3D_mod.f90
	for dir in $(SUBDIRS); do \
            cd $$dir && $(MAKE) clean && cd .. ; \
        done

dist:
	rm -f *HD* *MHD* *HMHD* *ROTH*
	$(MAKE) clean

uedit: 
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I./fftp-3 $(UTILS)/regridutils.fpp -o $(UTILS)/regridutils.f90
edit:
	$(CPP) -DN_=$(N) -DORD_=$(ORD) $(OBJSDIR)/pseudospec3D_mod.fpp -o \
	   $(OBJSDIR)/pseudospec3D_mod.f90
	cd $(FFTP) && $(MAKE) edit CPP='$(CPP)' IKIND='$(IKIND)' \
	   CSIZE='$(CSIZE)' NSTRIP='$(NSTRIP)' PRECISION='$(PRECISION)' \
           DO_HYBRID='DO_HYBRID$(P_HYBRID)'	
	if [ ${UNDER} -eq 1 ]; then \
	if [ -e main3D.fpp.bak ]; then \
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
	if [ -e main3D.fpp.bak ]; then \
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
	$(PREPROCESS) main3D.fpp -o $(SOLVER).f90
	$(COMPILE) $(SOLVER).f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(OBJS) \
           $(SOLVER).o $(LDFLAGS) -o ../bin/$(SOLVER)

trans$(SOLVER) triad$(SOLVER):
	$(PREPROCESS) $(EXTRASDIR)/$(@:$(SOLVER)=3D).fpp -o $@.f90
	$(COMPILE) $@.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(OBJS) \
           $(TRAN_OBJS) $@.o $(LDFLAGS) -o ../bin/$@

struc$(SOLVER):
	$(PREPROCESS) $(EXTRASDIR)/$(@:$(SOLVER)=3D).fpp -o $@.f90
	$(COMPILE) $@.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTP)/*.o $(IOLIB)/*.o $(OBJS) \
           $(STRU_OBJS) $@.o $(LDFLAGS) -o ../bin/$@

$(FFTP)/fftp3D.o:
	cd $(FFTP) && $(MAKE) COMPILE='$(COMPILE)'

$(UTILS)/gutils_mod.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)'
$(UTILS)/regridutils.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)'

$(IOLIB)/binary_io.o:
	cd $(IOLIB) && $(MAKE) COMPILE='$(COMPILE)' 

$(OBJS) $(TRAN_OBJS) $(STRU_OBJS):
	$(COMPILE) $(OBJSDIR)/$(@:.o=.f90)
