!=================================================================
! FFTP3D v2
! Parallel Fast Fourier Transform in 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the FFTW 2.x library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that calls any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!
! 16 Feb 2004: Performs complex FFTs in place.
!  8 Jul 2004: itype pointers only used to store datatypes of 
!              blocks in the row and column each processor is.
!  9 Jul 2004: Transposition uses data cache blocking.
! 13 Feb 2007: Transposition uses strip mining (rreddy@psc.edu)
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
!
! References:
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================

!*****************************************************************
      SUBROUTINE fftp3d_init_threads(err)
!-----------------------------------------------------------------
!
! Initializes FFTW threads.
!
! Parameters
!     err : if not zero, the initialization failed
!-----------------------------------------------------------------

!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: err

!$    CALL fftw_f77_threads_init(err)
      IF (err.ne.0) PRINT *,'FFTP threads initialization failed!'

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
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX(KIND=GP) (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_MEASURE (optimal but slower) or 
!              FFTW_ESTIMATE (sub-optimal but faster)
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      CALL rfftw2d_f77_create_plan(plan%planr,n,n,fftdir,flags)
      CALL fftw_f77_create_plan(plan%planc,n,fftdir, &
                               flags+FFTW_IN_PLACE)
      plan%n = n
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,nprocs,myrank,plan%itype1, &
                              plan%itype2)

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
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      CALL rfftwnd_f77_destroy_plan(plan%planr)
      CALL fftw_f77_destroy_plan(plan%planc)
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

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
      INTEGER, INTENT(IN) :: n,nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: ista,iend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      INTEGER :: itemp1,itemp2

      CALL range(1,n,nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         CALL range(1,n/2+1,nprocs,irank,ista,iend)
         CALL block3d(1,n/2+1,1,n,ksta,ista,iend,1,n,ksta, &
                     kend,GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      DO krank = 0,nprocs-1
         CALL range(1,n,nprocs,krank,ksta,kend)
         CALL block3d(ista,iend,1,n,1,ista,iend,1,n,ksta, &
                     kend,GC_COMPLEX,itemp2)
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

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%n,plan%n,ista:iend) :: out 
      COMPLEX(KIND=GP), DIMENSION(plan%n/2+1,plan%n,ksta:kend)          :: c1
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n,plan%n)              :: c2
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%n,plan%n,ksta:kend)     :: in

      DOUBLE PRECISION                    :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc
      INTEGER :: jfft,htra

!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft,GT_WTIME)
#ifdef DO_HYBRIDyes
      CALL rfftwnd_f77_threads_real_to_complex(nth,plan%planr,kend-ksta+1, &
                          in,1,plan%n*plan%n,c1,1,plan%n*(plan%n/2+1))
#else
      CALL rfftwnd_f77_real_to_complex(plan%planr,kend-ksta+1,in,          &
                          1,plan%n*plan%n,c1,1,plan%n*(plan%n/2+1))
#endif
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(htra,GT_WTIME)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs

            CALL MPI_IRECV(c2,1,plan%itype2(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c1,1,plan%itype1(isendTo),isendTo, &
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
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%n,csize
            DO kk = 1,plan%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%n,jj+csize-1)
               DO k = kk,min(plan%n,kk+csize-1)
                  out(k,j,i) = c2(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)

!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
#ifdef DO_HYBRIDyes
      CALL fftw_f77_threads(nth,plan%planc,plan%n*(iend-ista+1),out,1, &
                   plan%n,c2,1,plan%n)
#else
      CALL fftw_f77(plan%planc,plan%n*(iend-ista+1),out,1,plan%n,      &
                   c2,1,plan%n)
#endif
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL  GTFree(hfft); CALL GTFree(htra)


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

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%n,plan%n,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(plan%n/2+1,plan%n,ksta:kend)         :: c1
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n,plan%n)             :: c2
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%n,plan%n,ksta:kend)   :: out

      DOUBLE PRECISION                    :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc
      INTEGER :: hfft,htra

!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft,GT_WTIME)
#ifdef DO_HYBRIDyes
      CALL fftw_f77_threads(nth,plan%planc,plan%n*(iend-ista+1),in,1, &
                   plan%n,c2,1,plan%n)
#else
      CALL fftw_f77(plan%planc,plan%n*(iend-ista+1),in,1,plan%n, &
                   c2,1,plan%n)
#endif
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
!
! Cache friendly transposition
!

      CALL GTStart(htra,GT_WTIME)
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%n,csize
            DO kk = 1,plan%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%n,jj+csize-1)
               DO k = kk,min(plan%n,kk+csize-1)
                  c2(i,j,k) = in(k,j,i)
               END DO
               END DO
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

            CALL MPI_IRECV(c1,1,plan%itype1(igetFrom),igetFrom, & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(c2,1,plan%itype2(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)

!
! 2D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
#ifdef DO_HYBRIDyes
      CALL rfftwnd_f77_threads_complex_to_real(nth,plan%planr,kend-ksta+1, &
                         c1,1,plan%n*(plan%n/2+1),out,1,plan%n*plan%n)
#else
      CALL rfftwnd_f77_complex_to_real(plan%planr,kend-ksta+1,c1,          &
                         1,plan%n*(plan%n/2+1),out,1,plan%n*plan%n)
#endif
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      CALL  GTFree(hfft); CALL GTFree(htra)

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

      INTEGER, DIMENSION (2) :: iblock,idisp,itype

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin,jmax,kmin
      INTEGER, INTENT(IN)  :: ioldtype
      INTEGER, INTENT(OUT) :: inewtype

      INTEGER :: ilen,jlen,klen
      INTEGER :: isize,idist
      INTEGER :: itemp,itemp2
      INTEGER :: ierr

      CALL MPI_TYPE_EXTENT(ioldtype,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      klen = kend-ksta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      idist = (imax-imin+1)*(jmax-jmin+1)*isize
      CALL MPI_TYPE_HVECTOR(klen,1,idist,itemp,itemp2,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      iblock(1) = 1
      iblock(2) = 1
      idisp(1) = 0
      idisp(2) = ((imax-imin+1)*(jmax-jmin+1)*(ksta-kmin) &
                 +(imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      itype(1) = MPI_LB
      itype(2) = itemp2
      CALL MPI_TYPE_STRUCT(2,iblock,idisp,itype,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block3d
