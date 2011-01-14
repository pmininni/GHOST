!=================================================================
! MODULES for MPI binary I/O
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!=================================================================

  MODULE iovar
      TYPE IOPLAN
         INTEGER :: n,ksta,kend
         INTEGER :: iotype
      END TYPE IOPLAN
      SAVE

  END MODULE iovar
!=================================================================
 
  MODULE iompi
!     INCLUDE 'mpif.h'
      USE commtypes
      INTEGER, SAVE :: ioerr
      INTEGER(KIND=MPI_OFFSET_KIND), SAVE :: disp = 0
      INTEGER, SAVE :: bmangle = 1

  END MODULE iompi
!=================================================================
