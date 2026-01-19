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
!  Timer data
!
      INTEGER, PARAMETER                :: MAXLEVELS=128
      INTEGER, PARAMETER, PUBLIC        :: GTNULLHANDLE=-1
      INTEGER, PARAMETER, PUBLIC        :: GT_WTIME=0, GT_CPUTIME=1, GT_OMPTIME=2
      INTEGER, PARAMETER, PUBLIC        :: GTERR_GOOD_HANDLE=0, GTERR_INVAL_HANDLE=1, GTERR_STALE_HANDLE=2

      INTEGER                           :: ihandle_(0:MAXLEVELS-1)=GTNULLHANDLE
      INTEGER                           :: itype_(0:MAXLEVELS-1)=GT_WTIME
      DOUBLE PRECISION                  :: t0_(0:MAXLEVELS-1)=0.0D0 
      DOUBLE PRECISION                  :: t1_(0:MAXLEVELS-1)=0.0D0 
! end, member data
      
!
! Methods:
      PUBLIC                            :: GTStop, GTStart, GTInitHandle, GTAcc, GTFree
      PUBLIC                            :: GTGetTime, GTValidHandle, GTReset
      PUBLIC                            :: GTGetHandle

!
!
!
! Methods:
      CONTAINS

      
      SUBROUTINE GTInitHandle(ih,itype)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Initializes time level, returns handle to it if 
! successful. Handle is assumed not to be in use, and a new handle will
! be found for this interface, even if the old is not stale. 
!
! ARGUMENTS  :
!     ih     : (returned), Integer handle to timer level
!              (GTNULLHANDLE if unsuccessful)
!     itype  : timer type to use. See GT PARAMETERS for valid types
!-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)                     :: ih
      INTEGER, INTENT   (IN)                     :: itype
      INTEGER                                    :: ierr


      IF ( itype .LT. GT_WTIME .OR. itype .GT. GT_OMPTIME ) THEN
        WRITE(*,*) 'GTInitHandle: invalid time type: ', itype 
        STOP
      ENDIF

!     if handle is ok (currently in use), then just
!     reset time (time type is checked)
!!    ierr = GTValidHandle(ih)
!!    IF ( ierr.EQ.GTERR_GOOD_HANDLE ) THEN
!!      IF ( itype.NE.itype_(ih) ) THEN
!!        WRITE(*,*) 'GTInitHandle: invalid time type for valid handle: ', itype 
!!        STOP
!!      ENDIF
!!      t0_     (ih) = GTbasic(itype_(ih))      
!!      RETURN
!!    ENDIF

      ih   = GTGetHandle()
      ierr = GTValidHandle(ih)
      CALL GTHandleCatch(ih,ierr,'GTInitHandle')
      ihandle_(ih) = ih
      itype_  (ih) = itype
      t0_     (ih) = GTbasic(itype_(ih))      

      RETURN

      END SUBROUTINE GTInitHandle
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE GTStart(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Starts timer for specified handle. Handle must be valid, 
! created by call to GTInitHandle method.
!
! ARGUMENTS  :
!     ih     : (returned), Integer handle to timer level
!              (GTNULLHANDLE if unsuccessful)
!-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)                     :: ih
      INTEGER                                    :: ierr


!     If handle is ok (currently in use), then just
!     reset time (time type is checked)
      ierr = GTValidHandle(ih)
      IF ( ierr.NE.GTERR_GOOD_HANDLE ) THEN
        WRITE(*,*) 'GTStart: Invalide handle. Make sure GT'
      ENDIF
      t0_     (ih) = GTbasic(itype_(ih))

      RETURN

      END SUBROUTINE GTStart
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!
!
      SUBROUTINE GTStop(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Finalizes time level, and computes cumulative elapsed time 
! for specified ! handle. The handle is not deleted with this call; caller 
! is responsible for deleting it. 
!
! ARGUMENTS  :
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)                     :: ih
      INTEGER                                    :: ierr

   
      ierr = GTValidHandle(ih)
      IF ( ierr.NE.GTERR_GOOD_HANDLE ) THEN
        CALL GTHandleCatch(ih,ierr,'GTStop')
      ENDIF

      t1_(ih) = t1_(ih) + GTbasic(itype_(ih)) - t0_(ih);

      RETURN

      END SUBROUTINE GTStop
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!

      SUBROUTINE GTAcc(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Accumulates time. GTStart must still be called to initialze 
! handle. With each call, the total elapsed time is computed, and made
! available to caller. Caller must call GTFree explicitly.
! Do not call GTStop after this call, or the time history will end, 
! and accumulation will include only the time between the last
! GTAcc call, and the GTStop call.
!
! ARGUMENTS  :
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!-----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: ih
      INTEGER              :: ierr

      ierr = GTValidHandle(ih)
      IF ( ierr.NE.GTERR_GOOD_HANDLE ) THEN
        CALL GTHandleCatch(ih,ierr,'GTAcc')
      ENDIF

      t1_(ih) = t1_(ih) + (GTbasic(itype_(ih)) - t0_(ih))
      t0_(ih) = GTbasic(itype_(ih))
      
      RETURN

      END SUBROUTINE GTAcc
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      DOUBLE PRECISION FUNCTION GTGetTime(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Returns elapsed time for handle ih, after call to _either_ GTStop
! or GTAcc.
!
! ARGUMENTS  :
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!-----------------------------------------------------------------
      INTEGER, INTENT(IN)  :: ih
      INTEGER              :: ierr

      ierr = GTValidHandle(ih)
      IF ( ierr.NE.GTERR_GOOD_HANDLE ) THEN
        CALL GTHandleCatch(ih,ierr,'GTGetTime')
      ENDIF

      GTGetTime = t1_(ih);

      END FUNCTION GTGetTime
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      INTEGER FUNCTION GTGetHandle()
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Returns valid handle, or GTNULLHANDLE if unsuccessful
!
!-----------------------------------------------------------------

      INTEGER  :: j

      j = 0
      GTGetHandle = GTNULLHANDLE

      DO WHILE ( j.LT.MAXLEVELS .AND. ihandle_(j) .GE. 0 )
        j = j + 1
      ENDDO 
      IF ( j.LT.MAXLEVELS ) THEN
        GTGetHandle = j
      ENDIF

      END FUNCTION GTGetHandle
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE GTFree(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Frees up all data associated with handle ih, and nullifies it
! Nothing is done if handle is bad.
!
! ARGUMENTS  :
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!
!-----------------------------------------------------------------

      INTEGER, INTENT(INOUT)  :: ih

      IF ( GTValidHandle(ih).NE.GTERR_GOOD_HANDLE ) THEN
        RETURN
      ENDIF

      ihandle_(ih) = GTNULLHANDLE
      t0_     (ih) = 0.0D0
      t1_     (ih) = 0.0D0
      itype_  (ih) = 0

      END SUBROUTINE GTFree
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      INTEGER FUNCTION GTValidHandle(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Checks for handle validity, both that it's within bounds,
! and that it's not stale. If ok, then return GTERR_GOOD_HANDLE. If out
! of bounds, return GTERR_INVAL_HANDLE, if stale, return 
! GTERR_STALE_HANDLE
!
! ARGUMENTS  :
!     ih : Input integer handle to timer level (checked for validity,
!          but this won't affect elapsed time).
!
!-----------------------------------------------------------------

      INTEGER, INTENT(IN)  :: ih

      IF ( ih .LT. 0 .OR. ih .GE. MAXLEVELS ) THEN
        GTValidHandle = GTERR_INVAL_HANDLE
        RETURN
      ELSE IF ( ih .EQ. GTNULLHANDLE ) THEN
        GTValidHandle = GTERR_STALE_HANDLE
        RETURN
      ELSE
        GTValidHandle = GTERR_GOOD_HANDLE
      ENDIF

      END FUNCTION GTValidHandle
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE GTHandleCatch(ih,ierr,scaller)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! 
! DESCRIPTION: Catches handle errors, given error id, ierr, and reports
! errors & halts program
!
! ARGUMENTS  :
!     ih     : Input integer handle to timer level (checked for validity,
!              but this won't affect elapsed time).
!     ierr   : Input integer error code
!     scaller: error message: should be caller method name
!
!-----------------------------------------------------------------

      INTEGER, INTENT(IN)       :: ih
      INTEGER, INTENT(IN)       :: ierr
      CHARACTER*(*), INTENT(IN) :: scaller

      SELECTCASE(ierr)
        CASE(GTERR_INVAL_HANDLE)
          WRITE(*,*)trim(scaller),': invalid handle: ', ih
          STOP
        CASE(GTERR_STALE_HANDLE)
          WRITE(*,*)trim(scaller),': stale handle: ', ih
          STOP
      END SELECT

      END SUBROUTINE GTHandleCatch
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE GTReset(ih)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Resets time ranges. Useful for clearing times obtained
! from GTAcc call.
!
! ARGUMENTS  :
!     ih : input handle; must be valid 
!-----------------------------------------------------------------
      INTEGER, INTENT(INOUT)  :: ih

      IF ( GTValidHandle(ih).NE.GTERR_GOOD_HANDLE ) THEN
        RETURN
      ENDIF

      t0_     (ih) = 0.0D0
      t1_     (ih) = 0.0D0

      END SUBROUTINE GTReset
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      DOUBLE PRECISION FUNCTION GTbasic(itype)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
! DESCRIPTION: Returns time in seconds from arbitrary time in past
!
! ARGUMENTS  :
!     itype : time type
!-----------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'mpif.h'
!$    DOUBLE PRECISION, EXTERNAL :: omp_get_wtime
      INTEGER, INTENT(IN)  :: itype
      DOUBLE PRECISION     :: tt

      SELECTCASE(itype)
        CASE(0)
          tt = MPI_WTIME()
        CASE(1)
          CALL CPU_Time(tt)
        CASE(2)
         tt = 0.0D0
!$       tt = omp_get_wtime()
        CASE DEFAULT
          WRITE(*,*) 'GTbasic: invalid time type'
          STOP
      END SELECT
      GTBasic = tt

      END FUNCTION GTbasic
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
END MODULE gtimer
