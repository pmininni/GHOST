!=================================================================
! BINARY_IO subroutines
!
! Subroutines for POSIX binary I/O. When these subroutines are 
! called, binary native files are read or written by each 
! individual processor, and numbered according to the processor 
! rank and a time label.
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE io_init(myrank,n,ksta,kend,plan)
!-----------------------------------------------------------------
!
! Initializes variables for POSIX binary unformatted I/O.
!
! Parameters
!     myrank: the rank of the processor [IN]
!     n     : the size of the dimensions of the input array [IN]
!     ksta  : start value of the block in the third dimension [IN]
!     kend  : end value of the block in the third dimension [IN]
!     plan  : contains the I/O plan [OUT]
!-----------------------------------------------------------------

      USE iovar
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: myrank,n
      INTEGER, INTENT(IN) :: ksta,kend
      TYPE(IOPLAN), INTENT(OUT) :: plan

      plan%n = n
      plan%ksta = ksta
      plan%kend = kend
      WRITE(plan%node, fmtnod) myrank

      RETURN
      END SUBROUTINE io_init

!*****************************************************************
      SUBROUTINE io_read(unit,dir,fname,nmb,plan,var)
!-----------------------------------------------------------------
!
! Reads field components from individual POSIX binary 
! unformatted files, numbered by the node rank and by 
! an extension indicating the time.
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
      IMPLICIT NONE

      TYPE(IOPLAN), INTENT(IN) :: plan
      REAL, INTENT(OUT) :: var(plan%n,plan%n,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)      :: unit
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      OPEN(unit,file=trim(dir) // '/' // fname // '.' // & 
           plan%node // '.' // nmb // '.out',form='unformatted')
      READ(unit) var
      CLOSE(unit)

      RETURN
      END SUBROUTINE io_read

!*****************************************************************
      SUBROUTINE io_write(unit,dir,fname,nmb,plan,var)
!-----------------------------------------------------------------
!
! Writes field components into individual POSIX binary 
! unformatted files, numbered by the node rank and by 
! an extension indicating the time.
!
! Parameters
!     unit    : file handler [INOUT]
!     dir   : directory in which the files are written [IN]
!     fname : name of the field component [IN]
!     nmb   : extension with the time label [IN]
!     plan  : I/O plan [IN]
!     var   : the array with the field component [IN]
!-----------------------------------------------------------------

      USE iovar
      IMPLICIT NONE

      TYPE(IOPLAN), INTENT(IN) :: plan
      REAL, INTENT(IN) :: var(plan%n,plan%n,plan%ksta:plan%kend)
      INTEGER, INTENT(IN)      :: unit
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*), INTENT(IN)   :: nmb
      CHARACTER(len=*), INTENT(IN)   :: fname

      OPEN(unit,file=trim(dir) // '/' // fname // '.' // & 
           plan%node // '.' // nmb // '.out',form='unformatted')
      WRITE(unit) var
      CLOSE(unit)

      RETURN
      END SUBROUTINE io_write
