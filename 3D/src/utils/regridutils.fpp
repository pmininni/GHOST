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
! NOTE: ista/iend & ksta/kend for the truncated grid must be set
!       outside of this call, since they are set in fftplans module.
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

      ALLOCATE ( plan%ccarr(nt(3),nt(2),ista:iend)    )
      ALLOCATE ( plan%carr(nt(1)/2+1,nt(2),ksta:kend) )
      ALLOCATE ( plan%rarr(nt(1),nt(2),ksta:kend)     )

#if !defined(DEF_GHOST_CUDA_)
      IF ( fftdir.EQ.FFTW_REAL_TO_COMPLEX ) THEN
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planr,2,(/nt(1),nt(2)/),kend-ksta+1,  &
                       plan%rarr,(/nt(1),nt(2)*(kend-ksta+1)/),1,nt(1)*nt(2),     &
                       plan%carr,(/nt(1)/2+1,nt(2)*(kend-ksta+1)/),1,             &
		       (nt(1)/2+1)*nt(2),flags)
      ELSE
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planr,2,(/nt(1),nt(2)/),kend-ksta+1,  &
                       plan%carr,(/nt(1)/2+1,nt(2)*(kend-ksta+1)/),1,             &
		       (nt(1)/2+1)*nt(2),plan%rarr,(/nt(1),nt(2)*(kend-ksta+1)/), &
		       1,nt(1)*nt(2),flags)
      ENDIF
      CALL GPMANGLE(plan_many_dft)(plan%planc,1,nt(3),nt(2)*(iend-ista+1),        &
                       plan%ccarr,(iend-ista+1)*nt(2)*nt(3),1,nt(3),              &
                       plan%ccarr,(iend-ista+1)*nt(2)*nt(3),1,nt(3),fftdir,flags)
#endif
      plan%nx = nt(1)
      plan%ny = nt(2)
      plan%nz = nt(3)
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_trblock(n,nt,nprocs,myrank,plan%itype1,plan%itype2)

      RETURN
      END SUBROUTINE fftp3d_create_trplan

