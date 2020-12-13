!=================================================================
! FFTP2D v3
! Parallel Fast Fourier Transform in 2D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the FFTW 3.x library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that calls any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp2d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp2d_create_block'.
!
! 2010 Pablo D. Mininni
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!
! References:
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================
#include "fftw_wrappers.h"

!*****************************************************************
      SUBROUTINE fftp2d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
!
! Parameters
!     plan   : contains the parallel 2D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_MEASURE (optimal but slower) or 
!              FFTW_ESTIMATE (sub-optimal but faster)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      ALLOCATE ( plan%ccarr(n,ista:iend)    )
      ALLOCATE ( plan%carr(n/2+1,jsta:jend) )
      ALLOCATE ( plan%rarr(n,jsta:jend)     )
      IF (fftdir.eq.FFTW_REAL_TO_COMPLEX) THEN
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planr,1,n,jend-jsta+1, &
                              plan%rarr,n*(jend-jsta+1),1,n,       &
                              plan%carr,(n/2+1)*(jend-jsta+1),1,   &
                              n/2+1,flags)
      ELSE
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planr,1,n,jend-jsta+1, &
                         plan%carr,(n/2+1)*(jend-jsta+1),1,n/2+1,  &
                         plan%rarr,n*(jend-jsta+1),1,n,flags)
      ENDIF
      CALL GPMANGLE(plan_many_dft)(plan%planc,1,n,iend-ista+1,     &
                         plan%ccarr,n*(iend-ista+1),1,n,           &
                         plan%ccarr,n*(iend-ista+1),1,n,fftdir,flags)
      plan%n = n
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp2d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)

      RETURN
      END SUBROUTINE fftp2d_create_plan

!*****************************************************************
      SUBROUTINE fftp2d_destroy_plan(plan)
!-----------------------------------------------------------------
!
! Destroys FFTW plans in each node.
!
! Parameters
!     plan : the parallel 2D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      CALL GPMANGLE(destroy_plan)(plan%planr)
      CALL GPMANGLE(destroy_plan)(plan%planc)
      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr   )
      DEALLOCATE( plan%rarr   )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

      RETURN
      END SUBROUTINE fftp2d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp2d_create_block(n,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 2D matrix between processors. The data 
! types are used to transpose the matrix during the FFT.
!
! Parameters
!     n      : the size of the dimensions of the input array [IN]
!     nprocs : the number of processors [IN]
!     myrank : the rank of the processor [IN]
!     itype1 : contains a derived data type for sending [OUT]
!     itype2 : contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------

      USE commtypes
      IMPLICIT NONE

      INTEGER, INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n,nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: jsta,jend
      INTEGER :: irank,jrank
      INTEGER :: itemp1,itemp2

      CALL range(1,n,nprocs,myrank,jsta,jend)
      DO irank = 0,nprocs-1
         CALL range(1,n/2+1,nprocs,irank,ista,iend)
         CALL block2d(1,n/2+1,jsta,ista,iend,jsta,jend, &
                     GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      DO jrank = 0,nprocs-1
         CALL range(1,n,nprocs,jrank,jsta,jend)
         CALL block2d(ista,iend,1,ista,iend,jsta,jend,  &
                     GC_COMPLEX,itemp2)
         itype2(jrank) = itemp2
      END DO

      RETURN
      END SUBROUTINE fftp2d_create_block

!*****************************************************************
      SUBROUTINE fftp2d_real_to_complex(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 2D real-to-complex FFT in parallel. The 
! complex output has the same structure than the output 
! of the 2D FFTW, but the output is transposed.
!
! Parameters
!     plan : the 2D plan created with fftp2d_create_plan [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%n,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n)              :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%n,jsta:jend)     :: in

      DOUBLE PRECISION                    :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j
      INTEGER :: ii,jj
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

      CALL CPU_TIME(t0)

!
! 1D real-to-complex FFT in each node using the FFTW library
!
      CALL GPMANGLE(execute_dft_r2c)(plan%planr,in,plan%carr)
      CALL CPU_TIME(t1); ffttime = ffttime + t1-t0

      CALL CPU_TIME(t0)
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs

            CALL MPI_IRECV(c1,1,plan%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(plan%carr,1,plan%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
!
! Cache friendly transposition
!
      DO ii = ista,iend,csize
         DO jj = 1,plan%n,csize
            DO i = ii,min(iend,ii+csize-1)
            DO j = jj,min(plan%n,jj+csize-1)
               out(j,i) = c1(i,j)
            END DO
            END DO
         END DO
      END DO
      CALL CPU_TIME(t1); tratime = tratime + t1-t0
!
! 1D FFT in each node using the FFTW library
!
      CALL CPU_TIME(t0)
      CALL GPMANGLE(execute_dft)(plan%planc,out,out)
      CALL CPU_TIME(t1); ffttime = ffttime + t1-t0

      RETURN
      END SUBROUTINE fftp2d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp2d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 2D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 2D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 2D plan created with fftp2d_create_plan [IN]
!     in   : complex input array [IN]
!     out  : real output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%n,ista:iend) :: in
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n)             :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%n,jsta:jend)   :: out

      DOUBLE PRECISION                    :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j
      INTEGER :: ii,jj
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

      CALL CPU_TIME(t0)
!
! 1D FFT in each node using the FFTW library
!
      CALL GPMANGLE(execute_dft)(plan%planc,in,in)
      CALL CPU_TIME(t1); ffttime = ffttime + t1-t0

      CALL CPU_TIME(t0)
!
! Cache friendly transposition
!
      DO ii = ista,iend,csize
         DO jj = 1,plan%n,csize
            DO i = ii,min(iend,ii+csize-1)
            DO j = jj,min(plan%n,jj+csize-1)
               c1(i,j) = in(j,i)
            END DO
            END DO
         END DO
      END DO
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs

            CALL MPI_IRECV(plan%carr,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%itype2(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL CPU_TIME(t1); tratime = tratime + t1-t0

      CALL CPU_TIME(t0)
!
! 1D FFT in each node using the FFTW library
!
      CALL GPMANGLE(execute_dft_c2r)(plan%planr,plan%carr,out)
      CALL CPU_TIME(t1); ffttime = ffttime + t1-t0

      RETURN
      END SUBROUTINE fftp2d_complex_to_real

!*****************************************************************
      SUBROUTINE block2d(imin,imax,jmin,ista,iend,jsta,jend, &
                        ioldtype,inewtype)
!-----------------------------------------------------------------
!
! Soubroutine for defining derived data types in 2D.
!
! Parameters
!     imin : the minimum value in the first dimension [IN]
!     imax : the maximum value in the first dimension [IN]
!     jmin : the minimum value in the second dimension [IN]
!     ista : start value of the block in the first dimension [IN]
!     iend : end value of the block in the first dimension [IN]
!     jsta : start value of the block in the second dimension [IN]
!     jend : end value of the block in the second dimension [IN]
!     ioldtype: data type of the elements in the block [IN]
!     inewtype: the derived data type for the block [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE fftplans
      IMPLICIT NONE

      INTEGER, DIMENSION (2) :: iblock,idisp,itype

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin
      INTEGER, INTENT(IN)  :: ioldtype
      INTEGER, INTENT(OUT) :: inewtype

      INTEGER :: ilen,jlen
      INTEGER :: ilb,isize,idist
      INTEGER :: itemp
      INTEGER :: ierr

      CALL MPI_TYPE_GET_EXTENT(ioldtype,ilb,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      iblock(1) = 1
      iblock(2) = 1
      idisp(1) = 0
      idisp(2) = ((imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      itype(1) = MPI_LB
      itype(2) = itemp
      CALL MPI_TYPE_STRUCT(2,iblock,idisp,itype,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block2d