!=================================================================
! GHOST suite Fortran timing utilities 
!
! 2013 D. Rosenberg
!       ORNL
!
! 30 Jan 2013: Initial version
!=================================================================
MODULE gtimer
      IMPLICIT NONE
!  Timer data, if any:
!
! ...
! end, member data
      INTEGER, PARAMETER                :: MAXLEVELS=128, GTNULLHANDLE=-1
      INTEGER, PARAMETER                :: GT_WTIME=0, GT_CPUTIME=1, GT_OMPTIME=2

      INTEGER                           :: ihandle_(0:MAXLEVELS-1)
      INTEGER                           :: itype_(0:MAXLEVELS-1)
      DOUBLE PRECISION                  :: t0_(0:MAXLEVELS-1) 
      DOUBLE PRECISION                  :: t1_(0:MAXLEVELS-1) 

!
! Methods:
      PUBLIC                            :: GTStart, GTStop, GTAcc
      PUBLIC                            :: GTGetTime, GTFreeHandle
      PRIVATE                           :: GTbasic, GTGetHandle
!
!
! Methods:
      CONTAINS

      
      SUBROUTINE GTStart(handle,itype)
!-----------------------------------------------------------------
! Initializes time level, returns handle to it if successful
!
! Parameters
!     handle : (returned), Integer handle to timer level (-1 if unsuccessful)
!     itype  : timer type to use. See GT PARAMETERS for valid types
!-----------------------------------------------------------------
      USE mpivars
      USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(OUT)                       :: handle
      INTEGER, INTENT (IN)                       :: itype

      INTEGER                                    :: ih

      IF ( itype .NE. GT_WTIME .and. itype .NE. GT_CPUTIME ) THEN
        WRITE(*,*) 'GTStart: invalid time type'
        STOP
      ENDIF

      ih = GTGetHandle()
      IF ( ih .LT. 0 ) THEN
        WRITE(*,*) 'GTStart: invalid handle. Increase no. levels?'
        STOP
      ENDIF
      ihandle_(ih) = ih
      itype_  (ih) = itype
      t0_     (ih) = GTbasic(itype_(ih))      

      RETURN

      END SUBROUTINE GTStart
!
!
      SUBROUTINE GTStop(ih)
!-----------------------------------------------------------------
! Finalizes time level, and computes elapsed time for specified handle..
! The handle is not deleted with this call; caller is responsible for
! deleting it. 
!
! Parameters
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!-----------------------------------------------------------------
      USE mpivars
      USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(INOUT)                     :: ih

      IF ( ih .GE. 0 .and. ih .LT. MAXLEVELS ) THEN
        t1_(ih) = GTbasic(itype_(ih)) - t0_(ih);
      ELSE
        WRITE(*,*) 'GTStop: invalid handle: ', ih
        STOP
      ENDIF

      RETURN

      END SUBROUTINE GTStop
!
!

      SUBROUTINE GTAcc(ih)
!-----------------------------------------------------------------
! Accumulates time. GTStart must still be called to initialze handle
! With each call, the total elapsed time is computed, and made
! available to caller. Caller must call GTFreeHandle explicitly.
! Do not call GTStop after this call, or the time history will end, 
! and accumulation will include only the time between the last
! GTAcc call, and the GTStop call.
!
! Parameters
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!-----------------------------------------------------------------
      USE mpivars
      USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(IN)                        :: ih

      IF ( ih .GE. 0 .and. ih .LT. MAXLEVELS ) THEN
        t1_(ih) += GTbasic(itype_(ih)) - t0_(ih)
        t0_(ih)  = GTbasic(itype_(ih))
      ELSE
        WRITE(*,*) 'GTStop: invalid handle: ', handle
        STOP
      ENDIF
      
      RETURN

      END SUBROUTINE GTAcc
!
!
      DOUBLE PRECISION FUNCTION GTGetTime(ih)
!-----------------------------------------------------------------
! Returns elapsed time for handle ih, after call to _either_ GTStop
! or GTAcc.
!
! Parameters
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!-----------------------------------------------------------------

      INTEGER  :: ih

      IF ( ih .LT. 0 .OR. ih .GE. MAXLEVELS ) THEN
        WRITE(*,*) 'GTGetTime: invalid handle: ', ih
        STOP
      ENDIF
      GTGetTime = t1_(ih);

      END FUNCTION GTGetTime
!
!
      INTEGER FUNCTION GTGetHandle()
!-----------------------------------------------------------------
! Returns valid handle, or GTNULLHANDLE if unsuccessful
!
!-----------------------------------------------------------------

      INTEGER  :: ih, j

      j = 0
      GTGetHandle = GTNULLHANDLE

      DO WHILE ( j.LT. MAXLEVELS .AND. ihandle_(j) .LT. 0 )
        j = j + 1
      ENDDO 
      IF ( j.LT.MAXLEVELS ) THEN
        GTGetHandle = ih
      ENDIF

      END FUNCTION GTGetHandle
!
!
      INTEGER FUNCTION GTFreeHandle(ih)
!-----------------------------------------------------------------
! Frees up all data associated with handle ih, and nullifies it
!
! Parameters
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!
!-----------------------------------------------------------------

      INTEGER, INTENT(INOUT)  :: ih

      IF ( ih .LT. 0 .OR. ih .GE. MAXLEVELS ) THEN
        WRITE(*,*) 'GTFreeHandle: invalid handle: ', ih
        STOP
      ENDIF
      ihandle_(ih) = GTNULLHANDLE
      t0_     (ih) = 0.0D0
      t1_     (ih) = 0.0D0
      itype_  (ih) = 0

      END FUNCTION GTFreeHandle
!
!
      DOUBLE PRECISION FUNCTION GTbasic(itype)
!-----------------------------------------------------------------
! Returns time in seconds from arbitrary time in past
!
! Parameters
!     itype : time type
!-----------------------------------------------------------------
      INCLUDE 'mpif.h'
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: itype
      DOUBLE PRECISION     :: tt

      SELECTCASE(itype)
        CASE(0)
          GTbasic = MPI_WTIME()
        CASE(1)
          CALL CPU_Time(GTbasic)
        CASE(2)
         tt = 0.0D0
!$       tt = omp_get_wtime()
        CASE DEFAULT
          WRITE(*,*) 'GTbasic: invalid time type'
          STOP
      END SELECT

      END FUNCTION GTbasic
!
!
END MODULE gtimer
