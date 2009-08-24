#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Pablo Daniel Mininni - 10/04/2003
#****************************************************************************

# Compiler options
# For the Intel F95 compiler use -c -w -O when compiling. In RH9 systems, 
# you also need libc-2.2.93 and link with -lc-2.2.93. For the Absoft F90 
# compiler, use -c -YEXT_SFX=__ -YEXT_NAMES=LCS -O to compile and set the 
# variable UNDER to 0. When UNDER is equal to 1 you can use the variable 
# APPEND to add characters to external functions before compiling.
# When EXTEN is equal to 1 you can use the variable F90EXT 
# to change the extension of the F90/95 files.

#
# Desktop computers using the Portland copiler
COMP     = /usr/pgi/linux86/6.1/bin/pgf90 -c -g77libs -fastsse \
    -Msecond_underscore -Minline -Knoieee -Minfo=loop,inline
LINK     = /usr/pgi/linux86/6.1/bin/pgf90 -fastsse -Msecond_underscore \
    -Minline -Knoieee -Minfo=loop,inline
UNDER    = 0 
EXTEN    = 0 
APPEND   = _
F90EXT   = f90
LMPI     = -lmpich -lmpichf90 -lfmpich -lpthread -lrt
IDIR     = -I/home/mininni/library/mpich/include
LDIR     = -L/home/mininni/library/mpich/lib
FFTWLDIR = /home/mininni/library/fftw2/lib
LFFTW    = -lsrfftw -lsfftw -L$(FFTWLDIR)
SHELL    = /bin/bash
MAKE     = /usr/bin/make
PERL     = /usr/bin/perl
ECHO     = echo 

#
# Desktop computers using the GCC4 compiler
#COMP     = /home/mininni/gcc4/bin/gfortran -c -msse3 -O4 -fsecond-underscore
#LINK     = /home/mininni/gcc4/bin/gfortran -msse3 -O4 -fsecond-underscore
#UNDER    = 0 
#EXTEN    = 0 
#APPEND   = _
#F90EXT   = f90
#LMPI     = -llamf77mpi -lmpi -llam -lutil -ldl
#IDIR     = -I/home/mininni/mpilibs/include
#LDIR     = -L/home/mininni/mpilibs/lib
#FFTWLDIR = /home/mininni/fftlibs/lib
#LFFTW    = -lsrfftw -lsfftw -L$(FFTWLDIR)
#SHELL    = /bin/bash
#MAKE     = /usr/bin/make
#PERL     = /usr/bin/perl
#ECHO     = echo 

#
# Desktop computers using the Intel copiler
#COMP     = ifort -c -O3 -ip -pad -unroll -align -axP -w
#LINK     = ifort -O3 -ip -pad -unroll -align -axP -w
#UNDER    = 1 
#EXTEN    = 0 
#APPEND   = _
#F90EXT   = f90
#LMPI     = -llamf77mpi -lmpi -llam
#IDIR     =
#LDIR     =
#FFTWLDIR = /home/mininni/fftlibs/lib
#LFFTW    = -lsrfftw -lsfftw -L$(FFTWLDIR)
#SHELL    = bash
#MAKE     = make
#PERL     = perl
#ECHO     = echo 

default:
	@$(ECHO) "Usage: make <action>"
	@$(ECHO) "Possible actions are: hd2D, mhd2D, mhdu2D, hall25,"
	@$(ECHO) "clean, dist"

fftp.o:
	$(COMP) $(IDIR) fftp.$(F90EXT)

fftp2D.o:
	@$(MAKE) fftp.o
	$(COMP) $(IDIR) fftp_mod.$(F90EXT)
	$(COMP) $(IDIR) $(LDIR) $(LMPI) $(LFFTW) fftp2D.$(F90EXT)

pseudospec2D_inc.o:
	$(COMP) $(IDIR) pseudospec2D_mod.$(F90EXT)
	$(COMP) $(IDIR) pseudospec2D_inc.$(F90EXT)

pseudospec25D_inc.o:
	$(COMP) $(IDIR) pseudospec25D_inc.$(F90EXT)

hd2D:
	@$(MAKE) modify
	@$(MAKE) fftp2D.o
	@$(MAKE) pseudospec2D_inc.o
	$(COMP) $(IDIR) hd2D.para.$(F90EXT)
	$(LINK) fftp.o fftp2D.o fftp_mod.o pseudospec2D_mod.o \
          pseudospec2D_inc.o hd2D.para.o $(LDIR) $(LMPI) $(LFFTW) -o hd2D
	@$(MAKE) undo

mhd2D:
	@$(MAKE) modify
	@$(MAKE) fftp2D.o
	@$(MAKE) pseudospec2D_inc.o
	$(COMP) $(IDIR) mhd2D.para.$(F90EXT)
	$(LINK) fftp.o fftp2D.o fftp_mod.o pseudospec2D_mod.o \
          pseudospec2D_inc.o mhd2D.para.o $(LDIR) $(LMPI) $(LFFTW) -o mhd2D
	@$(MAKE) undo

mhdu2D:
	@$(MAKE) modify
	@$(MAKE) fftp2D.o
	@$(MAKE) pseudospec2D_inc.o
	$(COMP) $(IDIR) pseudospec2D_uni.$(F90EXT)
	$(COMP) $(IDIR) mhd2D.unif.$(F90EXT)
	$(LINK) fftp.o fftp2D.o fftp_mod.o pseudospec2D_mod.o \
          pseudospec2D_inc.o pseudospec2D_uni.o mhd2D.unif.o \
          $(LDIR) $(LMPI) $(LFFTW) -o mhdu2D
	@$(MAKE) undo

hall25:
	@$(MAKE) modify
	@$(MAKE) fftp2D.o
	@$(MAKE) pseudospec2D_inc.o
	@$(MAKE) pseudospec25D_inc.o
	$(COMP) $(IDIR) hallmhd25.para.$(F90EXT)
	$(LINK) fftp.o fftp2D.o fftp_mod.o pseudospec2D_mod.o \
          pseudospec2D_inc.o pseudospec25D_inc.o hallmhd25.para.o \
          $(LDIR) $(LMPI) $(LFFTW) -o hall25
	@$(MAKE) undo

modify:
	if [ ${UNDER} -eq 1 ]; then \
	$(PERL) -i.bak -pe 's/fftw_f77\(/fftw_f77${APPEND}\(/g' *.f90 ; \
	for item in `cat external`; do \
	  `echo '$(PERL) -i -pe 's/$$item/$${item}${APPEND}/g' *.f90'`; \
	done; fi
	if [ ${EXTEN} -eq 1 ]; then \
        for file in *.f90; do \
          if [ ${UNDER} -eq 0 ]; then \
             cp $$file $${file}.bak; \
          fi; \
          mv $$file $${file%.f90}.$(F90EXT); \
        done; fi

undo:
	rm -f *.f
	if [ -a fftp.f90.bak ]; then \
        for file in *.f90.bak; do \
          mv $$file $${file%.bak}; \
        done; fi

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*
	@$(MAKE) undo

dist:
	rm -f hd2D mhd2D mhdu2D hall25
	@$(MAKE) clean
