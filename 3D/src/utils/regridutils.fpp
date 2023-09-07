#if !defined(DEF_GHOST_CUDA_)
#include "fftw_wrappers.h"
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

      INTEGER, INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n(3),nt(3),nprocs
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
      SUBROUTINE fftp3d_create_trplan(plan,n,nt,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
! NOTE: itsta/itend & ktsta/ktend for the truncated grid must be set
!       outside of this call, since they are set in fftplans module.
!       Method trrange can be used to compute these.
!
!       In principle, nprocs should be the _full_ MPI_COMM_WORLD.
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
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3), nt(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      ALLOCATE ( plan%ccarr(nt(3),nt(2),itsta:itend)    )
      ALLOCATE ( plan%carr(nt(1)/2+1,nt(2),ktsta:ktend) )
      ALLOCATE ( plan%rarr(nt(1),nt(2),ktsta:ktend)     )

#if !defined(DEF_GHOST_CUDA_)
      IF ( fftdir.EQ.FFTW_REAL_TO_COMPLEX ) THEN
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planr,2,(/nt(1),nt(2)/),ktend-ktsta+1,  &
                       plan%rarr,(/nt(1),nt(2)*(ktend-ktsta+1)/),1,nt(1)*nt(2),     &
                       plan%carr,(/nt(1)/2+1,nt(2)*(ktend-ktsta+1)/),1,             &
		       (nt(1)/2+1)*nt(2),flags)
      ELSE
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planr,2,(/nt(1),nt(2)/),ktend-ktsta+1,  &
                       plan%carr,(/nt(1)/2+1,nt(2)*(ktend-ktsta+1)/),1,             &
		       (nt(1)/2+1)*nt(2),plan%rarr,(/nt(1),nt(2)*(ktend-ktsta+1)/), &
		       1,nt(1)*nt(2),flags)
      ENDIF
      CALL GPMANGLE(plan_many_dft)(plan%planc,1,nt(3),nt(2)*(itend-itsta+1),        &
                       plan%ccarr,(itend-itsta+1)*nt(2)*nt(3),1,nt(3),              &
                       plan%ccarr,(itend-itsta+1)*nt(2)*nt(3),1,nt(3),fftdir,flags)
#endif
      plan%nx = nt(1)
      plan%ny = nt(2)
      plan%nz = nt(3)
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_trblock(n,nt,nprocs,myrank,plan%itype1,plan%itype2)

      RETURN
      END SUBROUTINE fftp3d_create_trplan

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

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: newcomm, newgrp
      INTEGER, INTENT(IN)  :: n(3),nt(3), oldcomm
     
      INTEGER              :: ierr, grpworld
      INTEGER              :: nmt, np, nprocs, ntprocs
      INTEGER              :: iExclude(3,1), iInclude(3,1)

      CALL MPI_COMM_SIZE(oldcomm, nprocs, ierr)
      np      = nx / nprocs
      ntprocs = nxt / np
      nmt     = mod(nxt,np) 
      IF ( nmt .GT. 0 ) ntprocs = ntprocs + 1
      ntprocs = min(ntprocs, nprocs)

      CALL MPI_COMM_GROUP(oldcomm, grpworld, ierr)
      newcomm  = MPI_COMM_NULL
      newgrp = MPI_GROUP_NULL
      IF ( ntprocs .LT. nprocs ) THEN
        iExclude(1,1) = ntprocs
        iExclude(2,1) = nprocs-1
        iExclude(3,1) = 1
        CALL MPI_GROUP_RANGE_EXCL(grpworld, 1, iExclude, newgrp, ierr)   
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, newgrp, newcomm, ierr)
      ELSE IF ( ntprocs .EQ. nprocs ) THEN
        CALL MPI_COMM_DUP(MPI_COMM_WORLD,newcomm,ierr)
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD,newgrp,ierr)
      ENDIF

      RETURN
      END SUBROUTINE create_trcomm
