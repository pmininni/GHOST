!=================================================================
! BINARY_IO subroutines
!
! Subroutines for MPI binary I/O. When these subroutines are 
! called, binary native files are read or written using MPI I/O. 
! This file contains the 2D version of the subroutines.
!
! 2010 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE io_init(myrank,n,jsta,jend,plan)
!-----------------------------------------------------------------
!
! Initializes variables for MPI I/O. Creates plans and MPI 
! derived data types for I/O of the distributed real arrays 
! with the components of the fields.
!
! Parameters
!     myrank: the rank of the processor [IN]
!     n     : the size of the dimensions of the input array [IN]
!     jsta  : start value of the block in the third dimension [IN]
!     jend  : end value of the block in the third dimension [IN]
!     plan  : contains the I/O plan [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE iovar
      USE iompi
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: myrank,n
      INTEGER, INTENT(IN)   :: jsta,jend
      INTEGER, DIMENSION(2) :: sizes,subsizes,starts
      TYPE(IOPLAN), INTENT(OUT) :: plan

      plan%n = n
      plan%jsta = jsta
      plan%jend = jend

      sizes(1) = n
      sizes(2) = n
      subsizes(1) = n
      subsizes(2) = jend-jsta+1
      starts(1) = 0
      starts(2) = jsta-1
      CALL MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,starts, &
           MPI_ORDER_FORTRAN,GC_REAL,plan%iotype,ioerr)
      CALL MPI_TYPE_COMMIT(plan%iotype,ioerr)

      RETURN
      END SUBROUTINE io_init

!*****************************************************************
      SUBROUTINE io_read(unit,dir,fname,nmb,plan,var)
!-----------------------------------------------------------------
!
! Reads field components from individual MPI native 
! binary files, labeled only by the extension 
! indicating the time.
!
! Parameters
!     unit  : file unit [IN]
!     dir   : directory from which the files are read [IN]
!     fname : name of the field component [IN]
!     nmb   : extension with the time label [IN]
!     plan  : I/O plan [IN]
!     var   : the array with the field component [OUT]
!-----------------------------------------------------------------

      USE fprecision
      USE commtypes
      USE iovar
      USE iompi
      IMPLICIT NONE
      
      TYPE(IOPLAN), INTENT(IN)   :: plan
      REAL(KIND=GP), INTENT(OUT) :: var(plan%n,plan%jsta:plan%jend)
      INTEGER, INTENT(IN)        :: unit
      INTEGER                    :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
          '.' // nmb // '.out',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
      CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_READ_ALL(fh,var,plan%n*(plan%jend-plan%jsta+1), &
          GC_REAL,MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)

      RETURN
      END SUBROUTINE io_read

!*****************************************************************
      SUBROUTINE io_write(unit,dir,fname,nmb,plan,var)
!-----------------------------------------------------------------
!
! Writes field components into MPI native binary files, 
! labeled only by the extension indicating the time.
!
! Parameters
!     unit  : file unit [IN]
!     dir   : directory in which the files are written [IN]
!     fname : name of the field component [IN]
!     nmb   : extension with the time label [IN]
!     plan  : I/O plan [IN]
!     var   : the array with the field component [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE commtypes
      USE iovar
      USE iompi
      IMPLICIT NONE

      TYPE(IOPLAN), INTENT(IN)  :: plan
      REAL(KIND=GP), INTENT(IN) :: var(plan%n,plan%jsta:plan%jend)
      INTEGER, INTENT(IN)       :: unit
      INTEGER                   :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
          '.' // nmb // '.out',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,ioerr)
      CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_WRITE_ALL(fh,var,plan%n*(plan%jend-plan%jsta+1), &
          GC_REAL,MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)

      RETURN
      END SUBROUTINE io_write
