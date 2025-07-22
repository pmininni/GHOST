#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFT libraries
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
FPSPEC      = $(FPSPEC_$(COMPILER))
TARGET_ARCH = $(TARGET_$(COMPILER))
MPIINC      = $(MPIINC_$(MPI))
MPILIB      = $(MPILIB_$(MPI))
MPILDFLAGS  = $(MPILD_$(MPI))

# Anisotropic box with user-defined size
ARBSIZEDEFyes = -DDEF_ARBSIZE_
ARBSIZEDEFno  =

# OpenMP-MPI hybdridization and GPU offloading
DO_HYBRIDno   =
DO_HYBRIDyes  = $($(COMPILER)_OMP)
DO_HYBRIDoffl = $($(COMPILER)_OFFL)

# CUDA FFT & RT library variables:
CUCC         = $(CUCC_$(CUCOMPILER))
CUFLAGS      = $(CUFLAGS_$(CUCOMPILER)) -DG$(PRECISION)_PRECISION
CULIBDIR     = $(CULIBDIR_$(CUCOMPILER))
CULIB        = $(CULIB_$(CUCOMPILER))
CUINC        = $(CUINC_$(CUCOMPILER)) 
CUGENCODE    = $(CUGENCODE_$(CUCOMPILER))
CUARCH       = $(CUARCH_$(CUCOMPILER))

GHOME        = $(PWD)
CONTRIB      = $(GHOME)/contrib
UTILS        = $(GHOME)/utils
LAGP         = $(GHOME)/lagpart

# Lagrangian particle variables:
LAGPDEFyes   = -DDEF_GHOST_LAGP   
LAGPDEFlag   = $(LAGPDEFyes)
LAGPDEFiner  = -DDEF_GHOST_INERP
LAGPDEFtest  = -DDEF_GHOST_TESTP
LAGPDEFno    = 
LAGPOBJSyes  = $(LAGP)/gpdbrec_mod.o $(LAGP)/gpartcomm_mod.o \
	           $(LAGP)/gpsplineint_mod.o  $(LAGP)/gparts_mod.o
LAGPOBJSlag  = $(LAGPOBJSyes)
LAGPOBJSiner = $(LAGPOBJSyes)
LAGPOBJStest = $(LAGPOBJSyes)
LAGPOBJSno   = 

CUDA_GHOSTno =
CUDA_GHOSTyes= -DDEF_GHOST_CUDA_ -DG_PPN_=$(PPN) $(CUBIND_NVIDIA)
CUDA_OBJSno  =
CUDA_OBJSyes = $(CONTRIB)/nvidia/setaffinity_for_nvidia.o
GDEF         = $(IODEF$(IOLIB)) $(CUDA_GHOST$(P_CUDA)) \
               $(LAGPDEF$(PARTICLES)) $(ARBSIZEDEF$(ARBSIZE))
CONTRIBOBJS  = $(CUDA_OBJS$(P_CUDA))

# IO library variables:
IODEFvdcpio  = -DDEF_VDCPIO_
IODEFmpiioa  = 
IODEFposixio = 

IOLDIRvdcpio = -L$(IOLIBDIR)/lib
IOLDIRmpiio  = 
IOLDIRposixio= 

IOINCvdcpio  = -I$(IOLIBDIR)/include
IOINCmpiio   = 
IOINCposixio = 

IOINCDIR     = $(IOINC$(IOLIB))
IOLIBFLG     = $(IOLDIR$(IOLIB))

# FFTW library variables:
FFTW2SINGLE = s
FFTW2DOUBLE = d
FFTW2HYByes = -l$(FFTW2$(PRECISION))rfftw_threads \
	      -l$(FFTW2$(PRECISION))fftw_threads
FFTW2HYBoffl= $(FFTW2HYByes)
FFTW2HYBno  =
FFTW2FLAGS  = $(FFTW2HYB$(P_HYBRID)) -l$(FFTW2$(PRECISION))rfftw \
	      -l$(FFTW2$(PRECISION))fftw -lm
FFTW3SINGLE = f
FFTW3DOUBLE =
FFTW3HYByes = -lfftw3$(FFTW3$(PRECISION))_threads
FFTW3HYBoffl= $(FFTW3HYByes)
FFTW3HYBno  =
FFTW3FLAGS  = $(FFTW3HYB$(P_HYBRID)) -lfftw3$(FFTW3$(PRECISION)) -lm 
FFTWDIR     = $($(FFTP:fftp-%=FFTW%)DIR)
FFTWLDFLAGS = $($(FFTP:fftp-%=FFTW%)FLAGS) $(DO_HYBRID$(P_HYBRID))
FFTWINC     = -I$(FFTWDIR)/include
FFTWLIBDIR  = -L$(FFTWDIR)/lib

# MKL library variables:
FFTMKLDIR     = $(MKLDIR)
FFTMKLINC     = -I$(MKLDIR)/include
FFTMKLLIBDIR  = -L$(MKLDIR)/lib/intel64
FFTMKLLDno    =
FFTMKLLDyes   =
FFTMKLLDoffl  = -lmkl_sycl
FFTMKLLDFLAGS = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential \
                $(FFTMKLLD$(P_HYBRID))

# FFT Target variables derived from above FFT library quantities:
FFTTYPE3    = W
FFTTYPE2    = W
FFTTYPEmkl  = MKL
FFTTYPEyes  = CU
FFTTYPEno   = $(FFTTYPE$(FFTP:fftp-%=%))
FFTPDIRno   = $(FFTP)
FFTPDIRyes  = fftp-cu
FFTTYPE     = $(FFTTYPE$(P_CUDA))
FFTTLIBDIR  = $(FFT$(FFTTYPE)LIBDIR)
FFTTLIB     = $(FFT$(FFTTYPE)LIB)
FFTTINC     = $(FFT$(FFTTYPE)INC)
FFTTLDFLAGS = $(FFT$(FFTTYPE)LDFLAGS)

FFTPDIR     = $(GHOME)/$(FFTPDIR$(P_CUDA))

FFTCULIBDIR = $(CULIBDIR)
FFTCULIB    = $(CULIB)
FFTCUINC    = $(CUINC) -I$(FFTPDIR)
FFTCULDFLAGS= $(CUFLAGS)

LDFLAGS     = $(MPILIB) $(FFTTLIBDIR) $(FFTTLIB) $(MPILDFLAGS) \
	      $(FFTTLDFLAGS) $(IOLIBFLG)
INCLUDE     = $(MPIINC) $(FFTTINC) -I$(FFTPDIR) -I$(GHOME)/$(IOLIB) \
	      -I$(GHOME)/utils -I$(GHOME)/include $(IOINCDIR) -I$(GHOME) \
	      -I$(LAGP)
CPPFLAGS    = -D$(SOLVER)_SOL -DDO_HYBRID$(P_HYBRID) \
              $(CUDA_GHOST$(P_CUDA)) $(IODEF$(IOLIB))
PREPROCESS  = $(CPP) $(CPPFLAGS) $(GDEF)
COMPILE     = $(FC) $(FFLAGS) $(DO_HYBRID$(P_HYBRID)) $(INCLUDE) \
	      $(TARGET_ARCH) -c
SPECOMP     = $(FC) $(FPSPEC) $(DO_HYBRID$(P_HYBRID)) $(INCLUDE) \
	      $(TARGET_ARCH) -c
CCOMPILE    = $(CC) $(CFLAGS) $(INCLUDE) $(TARGET_ARCH) -c
CUCOMPILE   = $(CUCC) $(CUGENCODE) $(CUARCH) $(CUINC) -I$(FFTPDIR) \
	      $(CUFLAGS) -c
BIN         = ../bin

# Objects
HD_OBJS        = pseudospec3D_mod.o module_dns.o pseudospec3D_hd.o  
PHD_OBJS       = $(HD_OBJS) pseudospec3D_phd.o
MPHD_OBJS      = $(PHD_OBJS)
MHD_OBJS       = $(HD_OBJS) pseudospec3D_mhd.o
MHDB_OBJS      = $(MHD_OBJS) pseudospec3D_rot.o
RMHDB_OBJS     = $(MHDB_OBJS)
HMHD_OBJS      = $(MHD_OBJS)
HMHDB_OBJS     = $(MHDB_OBJS)
COMPRHD_OBJS   = $(PHD_OBJS) pseudospec3D_compr.o
CMHD_OBJS      = $(COMPRHD_OBJS) pseudospec3D_mhd.o
CMHDB_OBJS     = $(CMHD_OBJS) pseudospec3D_rot.o pseudospec3D_prot.o
CHMHD_OBJS     = $(CMHD_OBJS)
CHMHDB_OBJS    = $(CMHDB_OBJS)
ROTH_OBJS      = $(HD_OBJS) pseudospec3D_rot.o
PROTH_OBJS     = $(PHD_OBJS) pseudospec3D_rot.o pseudospec3D_prot.o
MPROTH_OBJS    = $(PROTH_OBJS) 
BOUSS_OBJS     = $(PROTH_OBJS) pseudospec3D_bouss.o lessgs_mod.o
MPBOUSS_OBJS   = $(BOUSS_OBJS)
ROTBOUSS_OBJS  = $(BOUSS_OBJS)
ROTBOUMHDB_OBJS= $(BOUSS_OBJS) pseudospec3D_mhd.o
MPROTBOUSS_OBJS= $(MPROTH_OBJS) pseudospec3D_bouss.o
GPE_OBJS       = $(PHD_OBJS) pseudospec3D_mhd.o pseudospec3D_gpe.o
ARGL_OBJS      = $(GPE_OBJS) pseudospec3D_newt.o
RGPE_OBJS      = $(GPE_OBJS) pseudospec3D_rot.o pseudospec3D_prot.o \
	         pseudospec3D_rgpe.o
RARGL_OBJS     = $(ARGL_OBJS)
ALPHA_OBJS     = module_lles.o pseudospec3D_filt.o
CAHD_OBJS      = $(HD_OBJS) $(ALPHA_OBJS) pseudospec3D_lles.o \
	         pseudospec3D_cahd.o
LAHD_OBJS      = $(HD_OBJS) $(ALPHA_OBJS) pseudospec3D_lles.o \
	         pseudospec3D_lahd.o
LHD_OBJS       = $(HD_OBJS) $(ALPHA_OBJS) pseudospec3D_lhd.o
LAMHD_OBJS     = $(MHD_OBJS) $(ALPHA_OBJS) pseudospec3D_lles.o \
	         pseudospec3D_lahd.o
EDQNM_OBJS     = module_edqnm.o pseudospec3D_edqnm.o
EDQNMHD_OBJS   = $(HD_OBJS) $(EDQNM_OBJS)
EDQNMROTH_OBJS = $(ROTH_OBJS) $(EDQNM_OBJS)
EDQNMBUSS_OBJS = $(BOUSS_OBJS) $(EDQNM_OBJS)
TRAN_OBJS      = pseudospec3D_trn.o
STRU_OBJS      = pseudospec3D_shf.o
OBJS           = $($(SOLVER)_OBJS) 
OBJSDIR        = pseudo
SUBDIRS        = posix mpiio fftp-2 fftp-3 fftp-cu fftp-mkl utils \
                 $(CONTRIB)/nvidia $(LAGP)
EXTRASDIR      = analysis
GTIMER         = $(UTILS)/gtimer_mod.o
GUTILS         = $(UTILS)/gutils_mod.o 
GSGS           = pseudo/lessgs_mod.o 
GREGRID        = $(UTILS)/gtrunc_mod.o $(UTILS)/regridutils.o
GIO            = $(IOLIB)/binary_io.o 

all: edit $(CONTRIBOBJS) $(GTIMER) $(FFTPDIR)/fftp3D.o  \
	  $(OBJS) $(GUTILS) $(GIO) $(LAGPOBJS$(PARTICLES)) $(SOLVER) undo

main: all

trans triad: edit $(GTIMER) $(FFTPDIR)/fftp3D.o $(OBJS) \
        $(GUTILS) $(GIO) $(TRAN_OBJS)
	$(MAKE) $@$(SOLVER)
	$(MAKE) undo

cancel: edit $(GTIMER) $(FFTPDIR)/fftp3D.o $(OBJS) \
           $(GUTILS) $(GIO)
	$(MAKE) $@$(SOLVER)
	$(MAKE) undo

struc:	edit $(GTIMER) $(FFTPDIR)/fftp3D.o $(OBJS) \
	   $(GUTILS) $(GIO) $(STRU_OBJS)
	$(MAKE) $@$(SOLVER)
	$(MAKE) undo

trunc:	edit $(GTIMER) $(FFTPDIR)/fftp3D.o  \
	   $(OBJS) $(GUTILS) $(GIO) $(GREGRID)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o  \
	   $(GTIMER) $(GUTILS) $(GREGRID) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

boots:  edit $(GTIMER) $(FFTPDIR)/fftp3D.o $(OBJS) \
	   $(GUTILS) $(GIO) $(GREGRID)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o \
	   $(GTIMER) $(GUTILS) $(GREGRID) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

lesbouss: edit $(CONTRIBOBJS) $(GTIMER) $(FFTPDIR)/fftp3D.o  \
      $(OBJS) $(GUTILS) $(GIO) $(GREGRID)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o  \
	   $(GTIMER) $(GUTILS) $(GREGRID) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

sgsbouss: edit $(CONTRIBOBJS) $(GTIMER) $(FFTPDIR)/fftp3D.o  \
      $(OBJS) $(GUTILS) $(GIO) $(GREGRID)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o  \
	   $(GTIMER) $(GUTILS) $(GREGRID) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

vt:     edit $(GTIMER) $(FFTPDIR)/fftp3D.o \
	   $(OBJS) $(GIO)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o \
	   $(GTIMER) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

shear:  edit $(GTIMER) $(FFTPDIR)/fftp3D.o $(GMIS) \
	   $(OBJS) $(GIO) $(GUTILS)
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o \
	   $(GTIMER) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

wv:     edit $(GTIMER) $(FFTPDIR)/fftp3D.o $(OBJS) $(GIO) $(GUTILS) 
	$(PREPROCESS) -DG$(PRECISION)_PRECISION -I$(FFTPDIR) \
	   $(EXTRASDIR)/$(@)3D.fpp -o $@.f90
	$(COMPILE) $@.f90 -o $@.o
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o \
	   $(GTIMER) $(OBJS) $@.o $(LDFLAGS) -o $(BIN)/$@
	$(MAKE) undo

clean:  undo
	rm -f *.o *.d *.s *.int *.inf *.mod *.modmic work.pc* ifc*
	rm -f $(OBJSDIR)/pseudospec3D_mod.f90
	for dir in $(SUBDIRS); do \
            cd $$dir && $(MAKE) clean && cd .. ; \
        done

dist:
	rm -f $(BIN)/*HD* $(BIN)/*MHD* $(BIN)/*ROT* $(BIN)/*BOUSS* \
              $(BIN)/*GPE* $(BIN)/*ARGL* $(BIN)/boots $(BIN)/trunc \
              $(BIN)/vt $(BIN)/shear $(BIN)/wv $(BIN)/wt $(BIN)/README_* \
              *HD* *MHD* *ROT* *BOUSS* \
              *GPE* *ARGL* boots trunc vt shear wv
	$(MAKE) clean

edit:
	$(CPP) -DNX_=$(NX) -DNY_=$(NY) -DNZ_=$(NZ) -DORD_=$(ORD) \
	   $(OBJSDIR)/pseudospec3D_mod.fpp -o $(OBJSDIR)/pseudospec3D_mod.f90
	cd $(FFTPDIR) && $(MAKE) edit CPP='$(CPP)' IKIND='$(IKIND)' \
	   CSIZE='$(CSIZE)' NSTRIP='$(NSTRIP)' PRECISION='$(PRECISION)' \
           NSTREAMS='$(NSTREAMS)' DO_HYBRID='DO_HYBRID$(P_HYBRID)'
	cd $(UTILS) && $(MAKE) edit CPP='$(CPP)' IKIND='$(IKIND)' \
	   CSIZE='$(CSIZE)' NSTRIP='$(NSTRIP)' PRECISION='$(PRECISION)' \
	   DO_HYBRID='DO_HYBRID$(P_HYBRID)' P_CUDA='$(P_CUDA)'
	cd $(LAGP) && $(MAKE) edit CPP='$(CPP)' CSIZE='$(CSIZE)' \
	   NSTRIP='$(NSTRIP)' DO_HYBRID='DO_HYBRID$(P_HYBRID)'
	cd ../
	echo "Output for this solver includes the following files:" >> \
	   $(BIN)/README_output.txt
	echo "(note some output files may not be created, depending" >> \
	   $(BIN)/README_output.txt
	echo "on specific solver and options being used)."  >> \
	   $(BIN)/README_output.txt
	echo >> $(BIN)/README_output.txt
	for file in $(OBJS:.o=.f90); do \
	 echo "From subroutines in" $$file "(see file for more details):"  >> \
	   $(BIN)/README_output.txt ; \
	 perl -ne 'print if /Output files contain/ .. /Parameters/' \
	   $(OBJSDIR)/$$file >> $(BIN)/README_output.txt ; \
	 echo >> $(BIN)/README_output.txt ; \
	done
	perl -i -lpe 's/! Parameters//g' $(BIN)/README_output.txt
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
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(DO_HYBRID$(P_HYBRID)) \
           $(FFTPDIR)/*.o $(IOLIB)/*.o $(CONTRIBOBJS) $(GUTILS) \
           $(LAGPOBJS$(PARTICLES)) $(GTIMER) $(OBJS) $(SOLVER).o \
           $(LDFLAGS) -o $(BIN)/$(SOLVER)

trans$(SOLVER) triad$(SOLVER):
	$(PREPROCESS) $(EXTRASDIR)/$(@:$(SOLVER)=3D).fpp -o $@.f90
	$(COMPILE) $@.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o $(OBJS) \
           $(CONTRIBOBJS) $(GUTILS) $(GTIMER) $(TRAN_OBJS) $@.o $(LDFLAGS) \
           -o $(BIN)/$@

cancel$(SOLVER):
	$(PREPROCESS) $(EXTRASDIR)/$(@:$(SOLVER)=3D).fpp -o $@.f90
	$(COMPILE) $@.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o $(OBJS) \
           $(GTIMER) $(GUTILS) $@.o $(LDFLAGS) -o $(BIN)/$@

struc$(SOLVER):
	$(PREPROCESS) $(EXTRASDIR)/$(@:$(SOLVER)=3D).fpp -o $@.f90
	$(COMPILE) $@.f90
	$(FC) $(FFLAGS) $(TARGET_ARCH) $(FFTPDIR)/*.o $(IOLIB)/*.o $(OBJS) \
           $(STRU_OBJS) $(GTIMER) $(GUTILS) $@.o $(LDFLAGS) -o $(BIN)/$@

$(FFTPDIR)/fftp3D.o: 
	cd $(FFTPDIR) && $(MAKE) COMPILE='$(COMPILE)' CUCOMPILE='$(CUCOMPILE)'

$(IOLIB)/binary_io.o:
	cd $(IOLIB) && $(MAKE) COMPILE='$(COMPILE)' 

$(LAGP)/gpdbrec_mod.o:
	cd $(LAGP) && $(MAKE) COMPILE='$(COMPILE)' gpdbrec_mod.o

$(LAGP)/gpartcomm_mod.o:
	cd $(LAGP) && $(MAKE) COMPILE='$(COMPILE)' gpartcomm_mod.o

$(LAGP)/gpsplineint_mod.o:
	cd $(LAGP) && $(MAKE) COMPILE='$(SPECOMP)' gpsplineint_mod.o

$(LAGP)/gparts_mod.o:
	cd $(LAGP) && $(MAKE) COMPILE='$(COMPILE)' gparts_mod.o

$(UTILS)/gutils_mod.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)' gutils_mod.o 

$(UTILS)/gtrunc_mod.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)' gtrunc_mod.o 

$(UTILS)/gtimer_mod.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)' gtimer_mod.o

$(UTILS)/regridutils.o:
	cd $(UTILS) && $(MAKE) COMPILE='$(COMPILE)' regridutils.o

$(CONTRIB)/nvidia/setaffinity_for_nvidia.o:
	cd $(CONTRIB)/nvidia && $(MAKE) CCOMPILE='$(CCOMPILE) -D_GNU_SOURCE' 

$(OBJS) $(TRAN_OBJS) $(STRU_OBJS):
	$(COMPILE) $(OBJSDIR)/$(@:.o=.f90)
