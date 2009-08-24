!=================================================================
! FFTP2D
! Parallel Fast Fourier Transform in 2D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the FFTW library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that call any of the subroutines in this 
! file. Also, you must create plans for the FFT using the 
! 'fftp2d_create_plan' subroutine.
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!
! 16 Feb 2004: Performs complex FFTs in place.
!  8 Jul 2004: itype pointers only store datatypes of blocks 
!              in the row and column each processor is.
!  9 Jul 2004: Transposition uses data cache blocking.
!
! References:
! Mininni PD, Montgomery DC, Pouquet A, Phys.Fluids 17, 035112 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================

!*****************************************************************
      SUBROUTINE fftp2d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTW in each node.
!
! Parameters
!     plan   : contains the parallel 2D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!              FFTW_FORWARD or FFTW_REAL_TO COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_MEASURE (sub-optimal but faster) or 
!              FFTW_ESTIMATE (optimal but slower)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
      IMPLICIT NONE

      TYPE (FFTPLAN) :: plan
      INTEGER :: n
      INTEGER :: fftdir
      INTEGER :: flags

      CALL rfftwnd_f77_create_plan(plan%planr,1,n,fftdir,flags)
      CALL fftw_f77_create_plan(plan%planc,n,fftdir, &
                               flags+FFTW_IN_PLACE)
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
!     plan : the parallel 2D plan [IN]
!-----------------------------------------------------------------

      USE fftplans
      IMPLICIT NONE

      TYPE (FFTPLAN) :: plan

      CALL rfftwnd_f77_destroy_plan(plan%planr)
      CALL fftw_f77_destroy_plan(plan%planc)
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
!     itype1 : contains a derived data type for sending [OUT]
!     itype2 : contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'mpif.h'
      
      INTEGER :: n,nprocs
      INTEGER, DIMENSION(0:nprocs-1) :: itype1,itype2

      INTEGER :: myrank
      INTEGER :: ista,iend
      INTEGER :: jsta,jend
      INTEGER :: irank,jrank
      INTEGER :: itemp1,itemp2

      CALL range(1,n,nprocs,myrank,jsta,jend)
      DO irank = 0,nprocs-1
         CALL range(1,n/2+1,nprocs,irank,ista,iend)
         CALL block2d(1,n/2+1,jsta,ista,iend,jsta,jend, &
                     MPI_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      DO jrank = 0,nprocs-1
         CALL range(1,n,nprocs,jrank,jsta,jend)
         CALL block2d(ista,iend,1,ista,iend,jsta,jend,  &
                     MPI_COMPLEX,itemp2)
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

      USE mpivars
      USE fftplans
      IMPLICIT NONE

      TYPE (FFTPLAN) :: plan

      COMPLEX, DIMENSION(plan%n,ista:iend)      :: out 
      COMPLEX, DIMENSION(plan%n/2+1,jsta:jend)  :: c1
      COMPLEX, DIMENSION(ista:iend,plan%n)      :: c2
      COMPLEX, DIMENSION(plan%n,ista:iend)      :: c3
      REAL, DIMENSION(plan%n,jsta:jend)         :: in
      INTEGER, DIMENSION(0:nprocs-1)            :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE)       :: istatus

      INTEGER :: comm
      INTEGER :: i,j
      INTEGER :: ii,jj
      INTEGER :: irank

!
! 1D real-to-complex FFT in each node using the FFTW library
!
      CALL rfftwnd_f77_real_to_complex(plan%planr,jend-jsta+1,in,1, &
                                      plan%n,c1,1,plan%n/2+1)
!
! Transposes the result between nodes
!
      DO irank = 0,nprocs-1
         IF (irank.ne.myrank) THEN
            CALL MPI_ISEND(c1,1,plan%itype1(irank),irank,1, &
                          comm,ireq1(irank),ierr)
            CALL MPI_IRECV(c2,1,plan%itype2(irank),irank,1, & 
                          comm,ireq2(irank),ierr)
         ELSE
            c2(ista:iend,jsta:jend) = c1(ista:iend,jsta:jend)
         ENDIF
      END DO
      DO irank = 0,nprocs-1
         IF (irank.ne.myrank) THEN
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         ENDIF
      END DO
!
! Cache friendly transposition
!
      DO ii = ista,iend,csize
         DO jj = 1,plan%n,csize
            DO i = ii,min(iend,ii+csize-1)
            DO j = jj,min(plan%n,jj+csize-1)
               out(j,i) = c2(i,j)
            END DO
            END DO
         END DO
      END DO
!
! 1D FFT in each node using the FFTW library
!
      CALL fftw_f77(plan%planc,iend-ista+1,out,1,plan%n,c2,1,plan%n)

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

      USE mpivars
      USE fftplans
      IMPLICIT NONE

      TYPE (FFTPLAN) :: plan

      COMPLEX, DIMENSION(plan%n,ista:iend)      :: in 
      COMPLEX, DIMENSION(plan%n/2+1,jsta:jend)  :: c1
      COMPLEX, DIMENSION(ista:iend,plan%n)      :: c2
      REAL, DIMENSION(plan%n,jsta:jend)         :: out
      INTEGER, DIMENSION(0:nprocs-1)            :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE)       :: istatus

      INTEGER :: comm
      INTEGER :: i,j
      INTEGER :: ii,jj
      INTEGER :: irank

!
! 1D FFT in each node using the FFTW library
!
      CALL fftw_f77(plan%planc,iend-ista+1,in,1,plan%n,c2,1,plan%n)
!
! Cache friendly transposition
!
      DO ii = ista,iend,csize
         DO jj = 1,plan%n,csize
            DO i = ii,min(iend,ii+csize-1)
            DO j = jj,min(plan%n,jj+csize-1)
               c2(i,j) = in(j,i)
            END DO
            END DO
         END DO
      END DO
!
! Transposes the result between nodes
!
      DO irank = 0,nprocs-1
         IF (irank.ne.myrank) THEN
            CALL MPI_ISEND(c2,1,plan%itype2(irank),irank,1, &
                          comm,ireq1(irank),ierr)
            CALL MPI_IRECV(c1,1,plan%itype1(irank),irank,1, & 
                          comm,ireq2(irank),ierr)
         ELSE 
            c1(ista:iend,jsta:jend) = c2(ista:iend,jsta:jend)
         ENDIF
      END DO
      DO irank = 0,nprocs-1
         IF (irank.ne.myrank) THEN
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         ENDIF
      END DO
!
! 1D FFT in each node using the FFTW library
!
      CALL rfftwnd_f77_complex_to_real(plan%planr,jend-jsta+1,c1,1, &
                                      plan%n/2+1,out,1,plan%n)

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

      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, DIMENSION (2) :: iblock,idisp,itype

      INTEGER :: ista,iend
      INTEGER :: jsta,jend
      INTEGER :: imin,imax
      INTEGER :: jmin
      INTEGER :: ioldtype,inewtype

      INTEGER :: ilen,jlen
      INTEGER :: isize,idist
      INTEGER :: itemp,itemp2
      INTEGER :: ierr

      CALL MPI_TYPE_EXTENT(ioldtype,isize,ierr)
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
      CALL MPI_TYPE_FREE(itemp)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block2d
