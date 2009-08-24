!=================================================================
! BINARY_IO subroutines
!
! Subroutines for MPI binary I/O. When these subroutines are 
! called, binary native files are read or written using MPI I/O.
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE io_init(myrank,n,ksta,kend,plan)
!-----------------------------------------------------------------
!
! Initializes variables for MPI I/O. Creates plans and MPI 
! derived data types for I/O of the distributed real arrays 
! with the components of the fields.
!
! Parameters
!     myrank: the rank of the processor [IN]
!     n     : the size of the dimensions of the input array [IN]
!     ksta  : start value of the block in the third dimension [IN]
!     kend  : end value of the block in the third dimension [IN]
!     plan  : contains the I/O plan [OUT]
!-----------------------------------------------------------------

      USE iovar
      USE iompi
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: myrank,n
      INTEGER, INTENT(IN)   :: ksta,kend
      INTEGER, DIMENSION(3) :: sizes,subsizes,starts
      TYPE(IOPLAN), INTENT(OUT) :: plan


      plan%n = n
      plan%ksta = ksta
      plan%kend = kend

      sizes(1) = n
      sizes(2) = n
      sizes(3) = n
      subsizes(1) = n
      subsizes(2) = n
      subsizes(3) = kend-ksta+1
      starts(1) = 0
      starts(2) = 0
      starts(3) = ksta-1
      CALL MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts, &
           MPI_ORDER_FORTRAN,MPI_REAL,plan%iotype,ioerr)
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

      USE iovar
      USE iompi
      IMPLICIT NONE
      
      TYPE(IOPLAN), INTENT(IN) :: plan
      REAL, INTENT(OUT) :: var(plan%n,plan%n,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)      :: unit
      INTEGER                  :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
          '.' // nmb // '.out',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
      CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_READ_ALL(fh,var, &
          plan%n*plan%n*(plan%kend-plan%ksta+1),MPI_REAL, &
          MPI_STATUS_IGNORE,ioerr)
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

      USE iovar
      USE iompi
      IMPLICIT NONE

      TYPE(IOPLAN), INTENT(IN) :: plan
      REAL, INTENT(IN) :: var(plan%n,plan%n,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)      :: unit
      INTEGER                  :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
          '.' // nmb // '.out',MPI_MODE_CREATE+MPI_MODE_WRONLY, &
          MPI_INFO_NULL,fh,ioerr)
      CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_WRITE_ALL(fh,var, &
          plan%n*plan%n*(plan%kend-plan%ksta+1),MPI_REAL, &
          MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)

      RETURN
      END SUBROUTINE io_write
