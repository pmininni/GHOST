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

      USE commtypes
      USE iovar
      USE iompi
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: myrank,n(3)
      INTEGER, INTENT(IN)   :: ksta,kend
      INTEGER, DIMENSION(3) :: subsizes,starts
      TYPE(IOPLAN), INTENT(OUT) :: plan

      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      plan%ksta = ksta
      plan%kend = kend

      subsizes(1) = n(1)
      subsizes(2) = n(2)
      subsizes(3) = kend-ksta+1
      starts(1) = 0
      starts(2) = 0
      starts(3) = ksta-1
      CALL MPI_TYPE_CREATE_SUBARRAY(3,n,subsizes,starts, &
           MPI_ORDER_FORTRAN,GC_REAL,plan%iotype,ioerr)
      CALL MPI_TYPE_COMMIT(plan%iotype,ioerr)
      CALL GTInitHandle(ihopen ,GT_CPUTIME)
      CALL GTInitHandle(ihread ,GT_CPUTIME)
      CALL GTInitHandle(ihwrite,GT_CPUTIME)

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
!     nmb   : extension with the time label [IN]. May turn off 
!             name-mangling by setting bmangle=0 before call.
!     plan  : I/O plan [IN]
!     var   : the array with the field component [OUT]
!-----------------------------------------------------------------

      USE fprecision
      USE commtypes
      USE iovar
      USE iompi
      USE gtimer
      USE mpivars
      USE gutils
      IMPLICIT NONE
      
      TYPE(IOPLAN),INTENT  (IN)      :: plan
      REAL(KIND=GP),INTENT(OUT) :: var(plan%nx,plan%ny,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)            :: unit
      INTEGER                        :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      CALL GTStart(ihopen)
      IF ( bmangle.EQ.1 ) THEN
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
          '.' // nmb // '.out',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
      IF ( ioerr.NE.MPI_SUCCESS ) THEN
         WRITE(*,*)': io_read: cannot open file for reading: ',      &
              trim(dir) // '/' // fname // '.' // nmb // '.out'
         STOP
      ENDIF
      ELSE

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname  &
                        ,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
      IF ( ioerr.NE.MPI_SUCCESS ) THEN
         WRITE(*,*)': io_read: cannot open file for reading: ',    &
              trim(dir) // '/' // fname 
         STOP
      ENDIF
      ENDIF
      CALL GTStop(ihopen)
      CALL GTStart(ihread)
      CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_READ_ALL(fh,var, &
          plan%nx*plan%ny*(plan%kend-plan%ksta+1),GC_REAL,         &
          MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)
      CALL GTStop(ihread)
      IF ( iswap.gt.0 ) THEN
        CALL rarray_byte_swap(var,plan%nx*plan%ny*(plan%kend-plan%ksta+1))
      ENDIF

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
!     nmb   : extension with the time label [IN]. May turn off
!             name-mangling by setting bmangle=0 before call.
!     plan  : I/O plan [IN]
!     var   : the array with the field component [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE commtypes
      USE iovar
      USE iompi
      USE gtimer
      USE mpivars
      USE gutils
      IMPLICIT NONE

      TYPE(IOPLAN), INTENT(IN)       :: plan
      REAL(KIND=GP), INTENT(INOUT) :: var(plan%nx,plan%ny,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)            :: unit
      INTEGER                        :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      CALL GTStart(ihopen)
      IF ( bmangle.EQ.1 ) THEN
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname // &
          '.' // nmb // '.out',MPI_MODE_CREATE+MPI_MODE_WRONLY,      &
          MPI_INFO_NULL,fh,ioerr)
      ELSE
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,fname            &
                        ,MPI_MODE_CREATE+MPI_MODE_WRONLY &
                        ,MPI_INFO_NULL,fh,ioerr)
      ENDIF
      CALL GTStop(ihopen)
      CALL GTStart(ihwrite)
      CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_WRITE_ALL(fh,var,                      &
          plan%nx*plan%ny*(plan%kend-plan%ksta+1),GC_REAL, &
          MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)
      CALL GTStop(ihwrite)

      RETURN
      END SUBROUTINE io_write

!*****************************************************************
      SUBROUTINE io_initc(myrank,n,ksta,kend,plan)
!-----------------------------------------------------------------
!
! Initializes variables for MPI I/O. Creates plans and MPI 
! derived data types for I/O of the distributed complex arrays 
! with the components of the fields.
!
! Parameters
!     myrank: the rank of the processor [IN]
!     n     : the size of the dimensions of the input array [IN]
!             For Fourier transforms of real fields, this should
!             be a vector (nz,ny,nx/2+1)       
!     ksta  : start value of the block in the third dimension [IN]
!     kend  : end value of the block in the third dimension [IN]
!     plan  : contains the I/O plan [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE iovar
      USE iompi
      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: myrank,n(3)
      INTEGER, INTENT(IN)   :: ksta,kend
      INTEGER, DIMENSION(3) :: subsizes,starts
      TYPE(IOPLAN), INTENT(OUT) :: plan

      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      plan%ksta = ksta
      plan%kend = kend

      subsizes(1) = n(1)
      subsizes(2) = n(2)
      subsizes(3) = kend-ksta+1
      starts(1) = 0
      starts(2) = 0
      starts(3) = ksta-1
      CALL MPI_TYPE_CREATE_SUBARRAY(3,n,subsizes,starts, &
           MPI_ORDER_FORTRAN,GC_COMPLEX,plan%iotype,ioerr)
      CALL MPI_TYPE_COMMIT(plan%iotype,ioerr)

      RETURN
      END SUBROUTINE io_initc

!*****************************************************************
      SUBROUTINE io_readc(unit,dir,fname,nmb,plan,var)
!-----------------------------------------------------------------
!
! Reads field components from individual MPI native 
! complex binary files, labeled only by the extension 
! indicating the time.
!
! Parameters
!     unit  : file unit [IN]
!     dir   : directory from which the files are read [IN]
!     fname : name of the field component [IN]
!     nmb   : extension with the time label [IN]
!     plan  : I/O plan [IN]
!     var   : the array with the complex field component [OUT]
!-----------------------------------------------------------------

      USE fprecision
      USE commtypes
      USE iovar
      USE iompi
      IMPLICIT NONE
      
      TYPE(IOPLAN), INTENT(IN)       :: plan
      COMPLEX(KIND=GP), INTENT(OUT) :: var(plan%nx,plan%ny,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)            :: unit
      INTEGER                        :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      IF ( bmangle.EQ.1 ) THEN
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname //  &
          '.' // nmb // '.out',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
      ELSE

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname     &
                        ,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
      ENDIF
      CALL MPI_FILE_SET_VIEW(fh,disp,GC_COMPLEX,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_READ_ALL(fh,var,                           &
          plan%nx*plan%ny*(plan%kend-plan%ksta+1),GC_COMPLEX,  &
          MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)

      RETURN
      END SUBROUTINE io_readc

!*****************************************************************
      SUBROUTINE io_writec(unit,dir,fname,nmb,plan,var)
!-----------------------------------------------------------------
!
! Writes field components into MPI native complex binary files, 
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

      TYPE(IOPLAN), INTENT(IN)       :: plan
      COMPLEX(KIND=GP), INTENT(IN) :: var(plan%nx,plan%ny,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)            :: unit
      INTEGER                        :: fh
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      IF ( bmangle.EQ.1 ) THEN
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname //  &
          '.' // nmb // '.out',MPI_MODE_CREATE+MPI_MODE_WRONLY,       &
          MPI_INFO_NULL,fh,ioerr)
      ELSE
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(dir) // '/' // fname     &
                        ,MPI_MODE_CREATE+MPI_MODE_WRONLY &
                        ,MPI_INFO_NULL,fh,ioerr)
      ENDIF
      CALL MPI_FILE_SET_VIEW(fh,disp,GC_COMPLEX,plan%iotype,'native', &
          MPI_INFO_NULL,ioerr)
      CALL MPI_FILE_WRITE_ALL(fh,var,                         &
          plan%nx*plan%ny*(plan%kend-plan%ksta+1),GC_COMPLEX, &
          MPI_STATUS_IGNORE,ioerr)
      CALL MPI_FILE_CLOSE(fh,ioerr)

      RETURN
      END SUBROUTINE io_writec
    
