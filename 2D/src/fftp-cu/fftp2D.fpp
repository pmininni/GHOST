!=================================================================
! FFTP2D v3
! Parallel Fast Fourier Transform in 2D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the CUDA library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that calls any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp2d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp2d_create_block'.
!
! 2011 Duane Rosenberg & Pablo D. Mininni
!      National Center for Atmospheric Research
!      e-mail: mininni@df.uba.ar 
!
! References:
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================

!*****************************************************************
      SUBROUTINE fftp2d_create_plan(plan,n,fftdir,flags)
!-----------------------------------------------------------------
!
! Creates plans for the FFTCU in each node.
!
! Parameters
!     plan   : contains the parallel 2D plan [OUT]
!     n      : the size of the dimensions of the input array [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTCU_REAL_TO_COMPLEX (-1)
!              or FFTCU_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTCU [IN]
!              Curently unused.
!-----------------------------------------------------------------

      USE mpivars
      USE fftplans
      USE fprecision
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      INTEGER             :: iret, i1(2), i2(2), i3(2)


      plan%szccd_= 2* n     *(iend-ista+1)*GP
      plan%szcd_ = 2*(n/2+1)*(jend-jsta+1)*GP
      plan%szrd_ =    n     *(jend-jsta+1)*GP

      iret = cudaMallocHost ( plan%pccarr_,  plan%szccd_ )
      iret = cudaMallocHost ( plan%pcarr_, plan%szcd_ )
      iret = cudaMallocHost ( plan%prarr_, plan%szrd_ )

      i1(1) = n    ; i1(2) = iend-ista+1;
      i2(1) = n/2+1; i2(2) = jend-jsta+1;
      i3(1) = n    ; i3(2) = jend-jsta+1;
      call c_f_pointer ( plan%pccarr_, plan%ccarr, i1 )
      call c_f_pointer ( plan%pcarr_ , plan%carr , i2 )
      call c_f_pointer ( plan%prarr_ , plan%rarr , i3 )

!     ALLOCATE ( plan%ccarr(n,ista:iend)    )
!     ALLOCATE ( plan%carr(n/2+1,jsta:jend) )
!     ALLOCATE ( plan%rarr(n,jsta:jend)     )
!
      IF (fftdir.eq.FFTCU_REAL_TO_COMPLEX) THEN
      iret = cudaMalloc(plan%cu_ccd_, plan%szccd_)
      iret = cudaMalloc(plan%cu_cd_ , plan%szcd_ )
      iret = cudaMalloc(plan%cu_rd_ , plan%szrd_ )
      ! reverse order of rank for C-calling order:
!     CALL c_f_pointer(plan%icuplanr_, plan%icuplanr_)
      iret = cufftPlan1d(plan%icuplanr_, n, CUFFT_R2C, jend-jsta+1)
      iret = cufftPlan1d(plan%icuplanc_, n, CUFFT_C2C, iend-ista+1)
      ELSE
      iret = cudaMalloc(plan%cu_ccd_, plan%szccd_)
      iret = cudaMalloc(plan%cu_cd_ , plan%szcd_ )
      iret = cudaMalloc(plan%cu_rd_ , plan%szrd_ )
      iret = cufftPlan1d(plan%icuplanr_, n    , CUFFT_C2R, jend-jsta+1)
      iret = cufftPlan1d(plan%icuplanc_, n, CUFFT_C2C, iend-ista+1)
      ENDIF
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
! Destroys FFTCU plans in each node.
!
! Parameters
!     plan : the parallel 2D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE cuda_bindings
      USE cutypes

      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      INTEGER                      :: iret

      iret = cufftDestroy(plan%icuplanr_)
      iret = cufftDestroy(plan%icuplanc_)
      iret = cudaFreeHost (plan%pccarr_)
      iret = cudaFreeHost (plan%pcarr_)
      iret = cudaFreeHost (plan%prarr_)
      iret = cudaFree(plan%cu_ccd_)
      iret = cudaFree(plan%cu_cd_)
      iret = cudaFree(plan%cu_rd_)
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
      USE cuda_bindings
      USE cutypes
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
! of the 2D FFTCU, but the output is transposed.
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
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%n,ista:iend) :: out
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n)              :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%n,jsta:jend)     :: in

      DOUBLE PRECISION                    :: t0, t1
      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,iret,j
      INTEGER :: ii,jj
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

!
! 1D real-to-complex FFT in each node using the FFTCU library
      
      DO j = jsta,jend
         DO i = 1,plan%n
            plan%rarr(i,j) = in(i,j)
         END DO
      END DO


! data sent to cufftXXXXXX must reside on device:
      call CPU_TIME(t0)
      iret = cudaMemCpyHost2Dev(plan%cu_rd_, plan%prarr_, plan%szrd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0
 
      call CPU_TIME(t0)
      IF ( GP.EQ. 4 ) THEN
        iret = cufftExecR2C(plan%icuplanr_, plan%cu_rd_, plan%cu_cd_)
      ELSE
        iret = cufftExecD2Z(plan%icuplanr_, plan%cu_rd_, plan%cu_cd_)
      ENDIF
      call CPU_TIME(t1); ffttime = ffttime + t1-t0
      call CPU_TIME(t0)
      iret = cudaMemCpyDev2Host(plan%pcarr_, plan%cu_cd_, plan%szcd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0

! NOTE: If nrocs = 1, then we can carry out the transpose directly
!       on the CUDA device.
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      call CPU_TIME(t0)
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
               plan%ccarr(j,i) = c1(i,j)
            END DO
            END DO
         END DO
      END DO
      call CPU_TIME(t1); tratime = tratime + t1-t0
!
! 1D FFT in each node using the FFTCU library
!
      call CPU_TIME(t0); 
      iret = cudaMemCpyHost2Dev(plan%cu_ccd_, plan%pccarr_, plan%szccd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0
      call CPU_TIME(t0); 
      IF ( GP.EQ. 4 ) THEN
      iret = cufftExecC2C(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_REAL_TO_COMPLEX)
      ELSE
      iret = cufftExecZ2Z(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_REAL_TO_COMPLEX)
      ENDIF
      call CPU_TIME(t1); ffttime = ffttime + t1-t0
      call CPU_TIME(t0); 
      iret = cudaMemCpyDev2Host(plan%pccarr_, plan%cu_ccd_, plan%szccd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0

      DO j = ista,iend
         DO i = 1,plan%n
             out(i,j) = plan%ccarr(i,j)
         END DO
      END DO

      RETURN
      END SUBROUTINE fftp2d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp2d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 2D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 2D FFTCU, but should be transposed. The real 
! output has the same order than the output of the FFTCU.
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
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%n,ista:iend) :: in
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n)             :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%n,jsta:jend)   :: out

      DOUBLE PRECISION                                          :: t0, t1

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,iret,j
      INTEGER :: ii,jj
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc

!
! 1D FFT in each node using the FFTCU library
!
!     plan%ccarr = in
      DO j = ista,iend
         DO i = 1,plan%n
             plan%ccarr(i,j) = in(i,j)
         END DO
      END DO
      call CPU_TIME(t0);
      iret = cudaMemCpyHost2Dev(plan%cu_ccd_, plan%pccarr_, plan%szccd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0
      call CPU_TIME(t0);
      IF ( GP.EQ. 4 ) THEN
      iret = cufftExecC2C(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_COMPLEX_TO_REAL)
      ELSE
      iret = cufftExecZ2Z(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_COMPLEX_TO_REAL)
      ENDIF
      call CPU_TIME(t1); ffttime = ffttime + t1-t0
      call CPU_TIME(t0);
      iret = cudaMemCpyDev2Host(plan%pccarr_, plan%cu_ccd_, plan%szccd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0

      call CPU_TIME(t0);
!
! Cache friendly transposition
!
      DO ii = ista,iend,csize
         DO jj = 1,plan%n,csize
            DO i = ii,min(iend,ii+csize-1)
            DO j = jj,min(plan%n,jj+csize-1)
               c1(i,j) = plan%ccarr(j,i)
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

      call CPU_TIME(t1); tratime = tratime + t1-t0
!
! 1D FFT in each node using the FFTCU library
!

      call CPU_TIME(t0);
      iret = cudaMemCpyHost2Dev(plan%cu_cd_, plan%pcarr_, plan%szcd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0
      call CPU_TIME(t0);
      IF ( GP.EQ. 4 ) THEN
      iret = cufftExecC2R(plan%icuplanr_, plan%cu_cd_, plan%cu_rd_)
      ELSE
      iret = cufftExecZ2D(plan%icuplanr_, plan%cu_cd_, plan%cu_rd_)
      ENDIF
      call CPU_TIME(t1); ffttime = ffttime + t1-t0
      call CPU_TIME(t0);
      iret = cudaMemCpyDev2Host(plan%prarr_, plan%cu_rd_, plan%szrd_ )
      call CPU_TIME(t1); memtime = memtime + t1-t0
      out = plan%rarr

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
      CALL MPI_TYPE_FREE(itemp,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block2d
