!=================================================================
! FFTP3D v3
! Parallel Fast Fourier Transform in 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the FFTW 3.x library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that calls any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!
! 16 Feb 2004: Performs complex FFTs in place.
!  8 Jul 2004: itype pointers only used to store datatypes of 
!              blocks in the row and column each processor is.
!  9 Jul 2004: Transposition uses data cache blocking.
! 13 Feb 2007: Transposition uses strip mining (rreddy@psc.edu)
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
!  3 Jan 2017: Anisotropic boxes (P. Mininni)
!
! References:
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================
#include "fftw_wrappers.h"

!*****************************************************************
      SUBROUTINE fftp3d_init_threads(err)
!-----------------------------------------------------------------
!
! Initializes FFTW threads.
!
! Parameters
!     err : if zero, the initialization failed
!-----------------------------------------------------------------

!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: err

!$    CALL GPMANGLE(init_threads)(err)
      IF (err.eq.0) PRINT *,'FFTP threads initialization failed!'

      RETURN
      END SUBROUTINE fftp3d_init_threads

!*****************************************************************
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_ESTIMATE (sub-optimal but faster)
!              FFTW_MEASURE (optimal but slower to create plans)
!              FFTW_PATIENT AND FFTW_EXHAUSTIVE are also available
!              for extra performance, but may take a long time to
!              create plans (specially when using OpenMP)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
!$    USE threads
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      ALLOCATE ( plan%ccarr(n(3),n(2),ista:iend)    )
      ALLOCATE ( plan%carr(n(1)/2+1,n(2),ksta:kend) )
      ALLOCATE ( plan%rarr(n(1),n(2),ksta:kend)     )
!$    CALL GPMANGLE(plan_with_nthreads)(nth)

      IF (fftdir.eq.FFTW_REAL_TO_COMPLEX) THEN
      CALL GPMANGLE(plan_many_dft_r2c)(plan%planr,2,(/n(1),n(2)/),    &
                         kend-ksta+1,plan%rarr,                       &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),     &
                         plan%carr,(/n(1)/2+1,n(2)*(kend-ksta+1)/),1, &
                         (n(1)/2+1)*n(2),flags)
      ELSE
      CALL GPMANGLE(plan_many_dft_c2r)(plan%planr,2,(/n(1),n(2)/),    &
                         kend-ksta+1,plan%carr,                       &
                         (/n(1)/2+1,n(2)*(kend-ksta+1)/),1,           &
                         (n(1)/2+1)*n(2),plan%rarr,                   &
                         (/n(1),n(2)*(kend-ksta+1)/),1,n(1)*n(2),flags)
      ENDIF
      CALL GPMANGLE(plan_many_dft)(plan%planc,1,n(3),n(2)*(iend-ista+1), &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         plan%ccarr,(iend-ista+1)*n(2)*n(3),1,n(3),      &
                         fftdir,flags)
      plan%nx = n(1)
      plan%ny = n(2)
      plan%nz = n(3)
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)

      CALL GTInitHandle(hcom,GT_WTIME)
      CALL GTInitHandle(hfft,GT_WTIME)
      CALL GTInitHandle(htra,GT_WTIME)
      CALL GTInitHandle(htot,GT_WTIME)


      RETURN
      END SUBROUTINE fftp3d_create_plan

!*****************************************************************
      SUBROUTINE fftp3d_destroy_plan(plan)
!-----------------------------------------------------------------
!
! Destroys FFTW plans in each node.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      CALL GPMANGLE(destroy_plan)(plan%planr)
      CALL GPMANGLE(destroy_plan)(plan%planc)
      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr   )
      DEALLOCATE( plan%rarr   )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

      CALL GTFree(hcom)
      CALL GTFree(hfft)
      CALL GTFree(htra)
      CALL GTFree(htot)

      RETURN
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
!
! Defines derived data types for sending and receiving 
! blocks of the 3D matrix between processors. The data 
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
      INTEGER, INTENT(IN) :: n(3),nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      INTEGER :: itemp1,itemp2

      CALL range(1,n(3),nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         CALL range(1,n(1)/2+1,nprocs,irank,ista,iend)
         CALL block3d(1,n(1)/2+1,1,n(2),ksta,ista,iend,1,n(2), &
                     ksta,kend,GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL range(1,n(1)/2+1,nprocs,myrank,ista,iend)
      DO krank = 0,nprocs-1
         CALL range(1,n(3),nprocs,krank,ksta,kend)
         CALL block3d(ista,iend,1,n(2),1,ista,iend,1,n(2),     &
                     ksta,kend,GC_COMPLEX,itemp2)
         itype2(krank) = itemp2
      END DO

      RETURN
      END SUBROUTINE fftp3d_create_block

!*****************************************************************
      SUBROUTINE fftp3d_real_to_complex(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D real-to-complex FFT in parallel. The 
! complex output has the same structure than the output 
! of the 3D FFTW, but the output is transposed. 
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : real input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%nz,plan%ny,ista:iend) :: out 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)          :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%nx,plan%ny,ksta:kend) :: in

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(htot)

      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_r2c)(plan%planr,in,plan%carr)
      CALL GTStop(hfft); 

!
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
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
      CALL GTStop(hcom); 
!
! Cache friendly transposition
!
      CALL GTStart(htra)
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  out(k,j,i) = c1(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra)
!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planc,out,out)
      CALL GTStop(hfft); 

      CALL GTStop(htot); 
     
      ! Update local accumulated timers:
      ffttime = GTGetTime(hfft)
      tratime = GTGetTime(htra)
      comtime = GTGetTime(hcom)
      tottime = GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp3d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 3D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 3D FFTW, but should be transposed. The real 
! output has the same order than the output of the FFTW.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 3D plan created with fftp3d_create_plan [IN]
!     in   : complex input array [IN]
!     out  : real output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%nz,plan%ny,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%ny,plan%nz)           :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%nx,plan%ny,ksta:kend) :: out

      DOUBLE PRECISION                    :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

!
! 1D FFT in each node using the FFTW library

      CALL GTStart(htot)

      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft)(plan%planc,in,in)
      CALL GTStop(hfft); 

!
! Cache friendly transposition
!
      CALL GTStart(htra)
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%ny,csize
            DO kk = 1,plan%nz,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%ny,jj+csize-1)
               DO k = kk,min(plan%nz,kk+csize-1)
                  c1(i,j,k) = in(k,j,i)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra)
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
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
      CALL GTStop(hcom); 
!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GPMANGLE(execute_dft_c2r)(plan%planr,plan%carr,out)
      CALL GTStop(hfft); 
      
      CALL GTStop(htot); 

      ! Update local accumulated timers:
      ffttime = GTGetTime(hfft)
      tratime = GTGetTime(htra)
      comtime = GTGetTime(hcom)
      tottime = GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp3d_complex_to_real

!*****************************************************************
      SUBROUTINE block3d(imin,imax,jmin,jmax,kmin,ista,iend, &
                        jsta,jend,ksta,kend,ioldtype,inewtype)
!-----------------------------------------------------------------
!
! Soubroutine for defining derived data types in 3D.
!
! Parameters
!     imin : the minimum value in the first dimension [IN]
!     imax : the maximum value in the first dimension [IN]
!     jmin : the minimum value in the second dimension [IN]
!     jmax : the maximum value in the second dimension [IN]
!     kmin : the minimum value in the third dimension [IN]
!     ista : start value of the block in the first dimension [IN]
!     iend : end value of the block in the first dimension [IN]
!     jsta : start value of the block in the second dimension [IN]
!     jend : end value of the block in the second dimension [IN]
!     ksta : start value of the block in the third dimension [IN]
!     kend : end value of the block in the third dimension [IN]
!     ioldtype: data type of the elements in the block [IN]
!     inewtype: the derived data type for the block [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE fftplans
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin,jmax,kmin
      INTEGER, INTENT(IN)  :: ioldtype
      INTEGER, INTENT(OUT) :: inewtype

      INTEGER(KIND=MPI_ADDRESS_KIND) :: ilb,isize,idist
      INTEGER :: ilen,jlen,klen
      INTEGER :: itemp,itemp2
      INTEGER :: ierr

      CALL MPI_TYPE_GET_EXTENT(ioldtype,ilb,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      klen = kend-ksta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      idist = (imax-imin+1)*(jmax-jmin+1)*isize
      CALL MPI_TYPE_CREATE_HVECTOR(klen,1,idist,itemp,itemp2,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      idist = ((imax-imin+1)*(jmax-jmin+1)*(ksta-kmin) &
              +(imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      CALL MPI_TYPE_CREATE_STRUCT(1,1,idist,itemp2,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block3d
