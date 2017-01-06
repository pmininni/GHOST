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
         INTEGER      :: ksta,kend
         INTEGER      :: nx,ny,nz
         INTEGER      :: iotype
      END TYPE IOPLAN
      INTEGER :: iswap=0
      INTEGER :: oswap=0
      SAVE

  END MODULE iovar
!=================================================================
 
  MODULE iompi
!     INCLUDE 'mpif.h'
      USE commtypes
      INTEGER, SAVE :: ioerr
      INTEGER(KIND=MPI_OFFSET_KIND), SAVE :: disp = 0
      INTEGER, SAVE :: bmangle = 1
      INTEGER, SAVE :: ihopen,ihread,ihwrite

  END MODULE iompi
!=================================================================
