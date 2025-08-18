#if !defined(DEF_GHOST_CUDA_)
!#include "fftw_wrappers.h"
#endif
!*****************************************************************
      SUBROUTINE fftp3d_create_trblock(n,nt,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 3D matrix between processors for truncation
! operator. The data types are used to transpose the matrix during 
! the FFT.
!
! Parameters
!     n      : the size of the dimensions of the (largest) input array [IN]
!     nt     : the truncation size [IN]
!     nprocs : the number of processors for untruncated data [IN]
!     myrank : the rank of the processor [IN]
!     itype1 : contains a derived data type for sending [OUT]
!     itype2 : contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------
      USE commtypes

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3),nt(3),nprocs
      INTEGER, INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      INTEGER :: itemp1,itemp2

      CALL trrange(1,n(3),nt(3),nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         CALL trrange(1,n(1)/2+1,nt(1)/2+1,nprocs,irank,ista,iend)
         CALL block3d(1,nt(1)/2+1,1,nt(2),ksta,ista,iend,1,nt(2), &
                      ksta,kend,GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL trrange(1,n(1)/2+1,nt(1)/2+1,nprocs,myrank,ista,iend)
      iend = min(iend,ista+nt(1))
      DO krank = 0,nprocs-1
         CALL trrange(1,n(3),nt(3),nprocs,krank,ksta,kend)
         CALL block3d(ista,iend,1,nt(2),1,ista,iend,1,nt(2),      &
                      ksta,kend,GC_COMPLEX,itemp2)
         itype2(krank) = itemp2
      END DO

      RETURN
      END SUBROUTINE fftp3d_create_trblock

!*****************************************************************
      SUBROUTINE trrange(n1,n2,nt,nprocs,irank,itsta,itend)
!-----------------------------------------------------------------
!
! Soubroutine for computing the local coordinate range 
! when splitting the original array among cores, for the
! truncation operation.
!
! Parameters
!     n1      : the minimum value in the splitted dimension [IN]
!     n2      : the maximum value in the splitted dimension [IN]
!     nt      : the maximum truncated value in the splitted dimension [IN]
!     nprocs  : the number of processors for untruncated data [IN]
!     irank   : the rank of the processor [IN]
!     itsta   : start value for the local truncation coordinate [OUT]
!     itend   : end value for the local truncation coordinate [OUT]
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n1,n2,nt
      INTEGER, INTENT(IN)  :: nprocs,irank
      INTEGER, INTENT(OUT) :: itsta,itend

      INTEGER :: iwork

      CALL range(n1,n2,nprocs,irank,itsta,iwork)
      IF ( itsta.gt.nt ) THEN
        itend = itsta-1
      ELSE
        itend = min(iwork,nt)
      ENDIF
      
      RETURN
      END SUBROUTINE trrange

!*****************************************************************
      SUBROUTINE fftp3d_create_trplan_comm(plan,n,nt,fftdir,flags,comm)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
! NOTE: itsta/itend & ktsta/ktend for the truncated grid must be set
!       outside of this call, since they are set in fftplans module.
!       Method trrange can be used to compute these.
!
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     nt     : the size of the dimensions of the truncated array [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_MEASURE (optimal but slower) or 
!              FFTW_ESTIMATE (sub-optimal but faster)
!     comm   : MPI communicator
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
!$    USE threads
      IMPLICIT NONE

      INCLUDE 'mpif.h'

      INTEGER, INTENT(IN) :: n(3), nt(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      INTEGER, INTENT(IN) :: comm
      TYPE(FFTPLAN), INTENT(OUT) :: plan


      IF ( comm .NE. MPI_COMM_NULL ) THEN
        CALL MPI_COMM_DUP(comm, plan%comm, ierr)
        IF ( ierr .NE. MPI_SUCCESS .OR. plan%comm .EQ. MPI_COMM_NULL ) THEN
          write(*,*) 'fftp3d_create_trplan_comm: MPI_COMM_DUP failed'
          STOP
        ENDIF
      ENDIF

      plan%nx   = nt(1)
      plan%ny   = nt(2)
      plan%nz   = nt(3)
      plan%ksta = ktsta
      plan%kend = ktend
      plan%ista = itsta
      plan%iend = itend
      plan%nprocs = 0
      plan%myrank = 0

      IF ( comm .EQ. MPI_COMM_NULL ) THEN
              RETURN
      ENDIF

      ALLOCATE ( plan%ccarr(nt(3),nt(2),plan%ista:plan%iend)    )
      ALLOCATE ( plan%carr(nt(1)/2+1,nt(2),plan%ksta:plan%kend) )
      ALLOCATE ( plan%rarr(nt(1),nt(2),plan%ksta:plan%kend)     )

      CALL MPI_COMM_SIZE(plan%comm,plan%nprocs,ierr)
      CALL MPI_COMM_RANK(plan%comm,plan%myrank,ierr)


#if !defined(DEF_GHOST_CUDA_)
      IF ( fftdir.EQ.FFTW_REAL_TO_COMPLEX ) THEN
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planr,2,(/nt(1),nt(2)/),plan%kend-plan%ksta+1,  &
                       plan%rarr,(/nt(1),nt(2)*(plan%kend-plan%ksta+1)/),1,nt(1)*nt(2),     &
                       plan%carr,(/nt(1)/2+1,nt(2)*(plan%kend-plan%ksta+1)/),1,             &
		       (nt(1)/2+1)*nt(2),flags)
      ELSE
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planr,2,(/nt(1),nt(2)/),plan%kend-plan%ksta+1,  &
                       plan%carr,(/nt(1)/2+1,nt(2)*(plan%kend-plan%ksta+1)/),1,             &
		       (nt(1)/2+1)*nt(2),plan%rarr,(/nt(1),nt(2)*(plan%kend-plan%ksta+1)/), &
		       1,nt(1)*nt(2),flags)
      ENDIF
      CALL GPMANGLE(plan_many_dft)(plan%planc,1,nt(3),nt(2)*(plan%iend-plan%ista+1),        &
                       plan%ccarr,(plan%iend-plan%ista+1)*nt(2)*nt(3),1,nt(3),              &
                       plan%ccarr,(plan%iend-plan%ista+1)*nt(2)*nt(3),1,nt(3),fftdir,flags)
#endif

      ! NOTE: In principle, mpivars::nprocs should be 
      !       the _full_ MPI_COMM_WORLD.
      
      ALLOCATE( plan%itype1(0:plan%nprocs-1) )
      ALLOCATE( plan%itype2(0:plan%nprocs-1) )
      CALL fftp3d_create_trblock(n,nt,plan%nprocs,plan%myrank,plan%itype1,plan%itype2)

      RETURN
      END SUBROUTINE fftp3d_create_trplan_comm


!*****************************************************************
      SUBROUTINE create_trcomm(n, nt, oldcomm, newcomm, newgrp)
!-----------------------------------------------------------------
!
! Create communicator for truncated grid operations.
!
! Parameters
!     n      : the size of the dimensions of the (largest) input array [IN]
!     nt     : the truncation size [IN]
!     oldcomm: original communicator
!     newcomm: new communicator
!     newgrp : new comm group
!-----------------------------------------------------------------
      USE commtypes

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: newcomm, newgrp
      INTEGER, INTENT(IN)  :: n(3),nt(3), oldcomm
     
      INTEGER              :: ierr, grpworld
      INTEGER              :: nmt, np, nprocs, ntprocs
      INTEGER              :: iExclude(3,1), iInclude(3,1)
      INTEGER              :: rlen
      INTEGER, PARAMETER   :: MAX_ERR_STRING_LEN = MPI_MAX_ERROR_STRING
      CHARACTER(LEN=MAX_ERR_STRING_LEN) :: serr


      CALL MPI_COMM_SIZE(oldcomm, nprocs, ierr)
      IF (ierr .NE. MPI_SUCCESS) THEN
        CALL MPI_Error_string(ierr, serr, rlen, ierr)
        WRITE(*,*) 'create_trcomm:', TRIM(serr(:rlen))
        STOP
      END IF
      np      = n(1) / nprocs
      ntprocs = nt(1) / np
      nmt     = mod(nt(1),np) 
      IF ( nmt .GT. 0 ) ntprocs = ntprocs + 1
      ntprocs = min(ntprocs, nprocs)

      CALL MPI_COMM_GROUP(oldcomm, grpworld, ierr)
      IF (ierr .NE. MPI_SUCCESS) THEN
        CALL MPI_Error_string(ierr, serr, rlen, ierr)
        WRITE(*,*) 'create_trcomm:', TRIM(serr(:rlen))
        STOP
      END IF
      newcomm  = MPI_COMM_NULL
      newgrp = MPI_GROUP_NULL
      IF ( ntprocs .LT. nprocs ) THEN
        iExclude(1,1) = ntprocs
        iExclude(2,1) = nprocs-1
        iExclude(3,1) = 1
        write(*,*) 'create_trcomm: 0: ntprocs=', ntprocs, ' nprocs=', nprocs
        CALL MPI_GROUP_RANGE_EXCL(grpworld, 1, iExclude, newgrp, ierr)   
        IF (ierr .NE. MPI_SUCCESS) THEN
          CALL MPI_Error_string(ierr, serr, rlen, ierr)
          WRITE(*,*) 'create_trcomm:', TRIM(serr(:rlen))
          STOP
        END IF
        IF (ierr .EQ. MPI_ERR_GROUP) THEN
          WRITE(*,*) 'create_trcomm: MPI_GROUP_RANGE_EXCL failed'
          STOP
        END IF
        CALL MPI_COMM_CREATE(oldcomm, newgrp, newcomm, ierr)
        IF (ierr .NE. MPI_SUCCESS) THEN
          CALL MPI_Error_string(ierr, serr, rlen, ierr)
          WRITE(*,*) 'create_trcomm:', TRIM(serr(:rlen))
          STOP
        END IF
        IF (ierr .EQ. MPI_ERR_COMM) THEN
          WRITE(*,*) 'create_trcomm: MPI_ERR_COMM'
          STOP
        END IF
        IF (ierr .EQ. MPI_ERR_GROUP) THEN
          WRITE(*,*) 'create_trcomm: MPI_ERR_GROUP'
          STOP
        END IF
      ELSE IF ( ntprocs .EQ. nprocs ) THEN
        write(*,*) 'create_trcomm: 1'
        CALL MPI_COMM_DUP(oldcomm,newcomm,ierr)
        IF (ierr .NE. MPI_SUCCESS) THEN
          CALL MPI_Error_string(ierr, serr, rlen, ierr)
          WRITE(*,*) 'create_trcomm:', TRIM(serr(:rlen))
          STOP
        END IF
        IF (ierr .EQ. MPI_ERR_COMM) THEN
          WRITE(*,*) 'create_trcomm: MPI_ERR_COMM'
          STOP
        END IF
        CALL MPI_COMM_GROUP(newcomm,newgrp,ierr)
        IF (ierr .NE. MPI_SUCCESS) THEN
          CALL MPI_Error_string(ierr, serr, rlen, ierr)
        END IF
      ELSE
        WRITE(*,*) 'create_trcomm: ntprocs invalid: ntprocs= ', ntprocs, ' nprocs=', nprocs
        STOP
      ENDIF

      RETURN
      END SUBROUTINE create_trcomm
