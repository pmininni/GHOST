!=================================================================
! FFTP3D v3
! Parallel Fast Fourier Transform in 2D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the CUDA library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'fftp_mod.f90') 
! in each program that CALLs any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
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
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir,flags)
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

      INTEGER             :: istr, idist, ostr, odist, nrank
      INTEGER             :: iret, i1(3), i2(3), i3(3)
      INTEGER             :: na(2), pinembed(2), ponembed(2)
      INTEGER             :: iplan


      plan%szccd_= 2* n     *n*(iend-ista+1)*GP
      plan%szcd_ = 2*(n/2+1)*n*(kend-ksta+1)*GP
      plan%szrd_ =    n     *n*(kend-ksta+1)*GP

      iret = cudaHostAlloc ( plan%pccarr_, plan%szccd_, cudaHostAllocPortable) 
      iret = cudaHostAlloc ( plan%pcarr_ , plan%szcd_ , cudaHostAllocPortable) 
      iret = cudaHostAlloc ( plan%prarr_ , plan%szrd_ , cudaHostAllocPortable)

      i1(1) = n    ; i1(2) = n; i1(3) = iend-ista+1;
      i2(1) = n/2+1; i2(2) = n; i2(3) = kend-ksta+1; 
      i3(1) = n    ; i3(2) = n; i3(3) = kend-ksta+1; 
      CALL c_f_pointer ( plan%pccarr_, plan%ccarr, i1 )
      CALL c_f_pointer ( plan%pcarr_ , plan%carr , i2 )
      CALL c_f_pointer ( plan%prarr_ , plan%rarr , i3 )


      iret = cudaMalloc(plan%cu_ccd_, plan%szccd_)
      iret = cudaMalloc(plan%cu_cd_ , plan%szcd_ )
      iret = cudaMalloc(plan%cu_rd_ , plan%szrd_ )
      IF (fftdir.eq.FFTCU_REAL_TO_COMPLEX) THEN
        nrank= 2
        na      (1) = n            ; na      (2) = n              ;
        pinembed(1) = n            ; pinembed(2) = n*(kend-ksta+1);            
        istr        = 1            ; idist       = n*n            ;
        ponembed(1) = n/2+1        ; ponembed(2) = n*(kend-ksta+1);            
        ostr        = 1            ; odist       = n*(n/2+1)      ;
        iret = cufftPlanMany(plan%icuplanr_,2,na,pinembed,istr,idist,&
                             ponembed,ostr,odist,CUFFT_R2C,kend-ksta+1);
        nrank       = 1
        na      (1) = n                ; 
        pinembed(1) = n*n*(iend-ista+1);
        istr        = 1                ; idist       = n           ;
        ponembed(1) = n*n*(iend-ista+1); 
        ostr        = 1                ; odist       = n           ; 
        iret = cufftPlanMany(plan%icuplanc_,nrank,na,pinembed,istr,idist,&
                             ponembed,ostr,odist,CUFFT_C2C,n*(iend-ista+1));
      ELSE
        nrank= 2;
        na      (1) = n            ; na      (2) = n              ;
        pinembed(1) = n/2+1        ; pinembed(2) = n*(kend-ksta+1);
        istr        = 1            ; idist       = n*(n/2+1)      ;
        ponembed(1) = n            ; ponembed(2) = n*(kend-ksta+1);
        ostr        = 1            ; odist       = n*n            ; 
        iret = cufftPlanMany(plan%icuplanr_,nrank,na,pinembed,istr,idist,&
                             ponembed,ostr,odist,CUFFT_C2R,kend-ksta+1);
        nrank       = 1
        na      (1) = n                ; 
        pinembed(1) = n*n*(iend-ista+1);
        istr        = 1                ; idist       = n           ;
        ponembed(1) = n*n*(iend-ista+1); 
        ostr        = 1                ; odist       = n           ; 
        iret = cufftPlanMany(plan%icuplanc_,nrank,na,pinembed,istr,idist,&
                             ponembed,ostr,odist,CUFFT_C2C,n*(iend-ista+1));
      ENDIF
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
! Destroys FFTCU plans in each node.
!
! Parameters
!     plan : the parallel 2D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE cuda_bindings
      USE cutypes
      USE threads

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
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nprocs,myrank,itype1,itype2)
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
! of the 3D FFTCU, but the output is transposed.
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
      USE, INTRINSIC :: iso_c_binding
      USE cuda_bindings
      USE cutypes
      USE threads
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%n,plan%n,ista:iend) :: out 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n,plan%n)              :: c1
      REAL(KIND=GP), INTENT(IN), DIMENSION(plan%n,plan%n,ksta:kend)     :: in

      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,iret,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc
      INTEGER :: hfft,htra,hmem

!
! 2D real-to-complex FFT in each node using the FFTCU library
     
      plan%rarr = in
!
! Data sent to cufftXXXXXX must reside on device:
      CALL GTStart(hmem,GT_WTIME)
      iret = cudaMemCpyHost2Dev(plan%cu_rd_, plan%prarr_, plan%szrd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)
 
      CALL GTStart(hfft,GT_WTIME)
      IF ( GP.EQ. 4 ) THEN
        iret = cufftExecR2C(plan%icuplanr_, plan%cu_rd_, plan%cu_cd_ )
      ELSE
        iret = cufftExecD2Z(plan%icuplanr_, plan%cu_rd_, plan%cu_cd_ )
      ENDIF

      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStart(hmem)
      iret = cudaMemCpyDev2Host(plan%pcarr_, plan%cu_cd_, plan%szcd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)

! NOTE: If nrocs = 1, then we can carry out the transpose directly
!       on the CUDA device.
!
! Transposes the result between nodes using 
! strip mining when nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(htra,GT_WTIME)
      DO iproc = 0, nprocs-1, nstrip
         DO istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            IF ( isendTo .GE. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            IF ( igetFrom .LT. 0 ) igetFrom = igetFrom + nprocs

            CALL MPI_IRECV(c1,1,plan%itype2(igetFrom),igetFrom,      & 
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(plan%carr,1,plan%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         ENDDO

         DO istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         ENDDO
      ENDDO
!
! Cache friendly transposition
!
!
!$omp parallel do if ((iend-ista)/csize.ge.nth) private (jj,kk,i,j,k)
      DO ii = ista,iend,csize
!$omp parallel do if ((iend-ista)/csize.lt.nth) private (kk,i,j,k)
         DO jj = 1,plan%n,csize
            DO kk = 1,plan%n,csize
               DO i = ii,min(iend,ii+csize-1)
               DO j = jj,min(plan%n,jj+csize-1)
               DO k = kk,min(plan%n,kk+csize-1)
                 !ReCALL that ccarr is dimensioned (:,:), starting at (1,1):
                  plan%ccarr(k,j,i-ista+1) = c1(i,j,k)
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)
!
! 1D FFT in each node using the FFTCU library
!
      CALL GTStart(hmem)
      iret = cudaMemCpyHost2Dev(plan%cu_ccd_, plan%pccarr_, plan%szccd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)
      CALL GTStart(hfft)
      IF ( GP.EQ. 4 ) THEN
      iret = cufftExecC2C(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_REAL_TO_COMPLEX)
      ELSE
      iret = cufftExecZ2Z(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_REAL_TO_COMPLEX)
      ENDIF
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStart(hmem)
      iret = cudaMemCpyDev2Host(plan%pccarr_, plan%cu_ccd_, plan%szccd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)

      CALL  GTFree(hfft); CALL GTFree(htra); CALL GTFree(hmem)

      out = plan%ccarr

      RETURN
      END SUBROUTINE fftp3d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
!
! Computes the 2D complex-to-real FFT in parallel. The 
! complex input has the same structure than the input 
! of the 2D FFTCU, but should be transposed. The real 
! output has the same order than the output of the FFTCU.
! The input data is destroyed during the computation.
!
! Parameters
!     plan : the 2D plan created with fftp3d_create_plan [IN]
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
      USE threads
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%n,plan%n,ista:iend) :: in 
      COMPLEX(KIND=GP), DIMENSION(ista:iend,plan%n,plan%n)             :: c1
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%n,plan%n,ksta:kend)   :: out


      INTEGER, DIMENSION(0:nprocs-1)      :: ireq1,ireq2
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
      INTEGER, INTENT(IN)                 :: comm
      INTEGER :: i,iret,j,k
      INTEGER :: ii,jj,kk
      INTEGER :: irank
      INTEGER :: isendTo, igetFrom
      INTEGER :: istrip,iproc
      INTEGER :: hfft,htra,hmem

!
! 1D FFT in each node using the FFTCU library
!
      plan%ccarr = in
     
      CALL GTStart(hmem,GT_WTIME);
      iret = cudaMemCpyHost2Dev(plan%cu_ccd_, plan%pccarr_, plan%szccd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)
      CALL GTStart(hfft);
      IF ( GP.EQ. 4 ) THEN
        iret = cufftExecC2C(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_COMPLEX_TO_REAL)
!       iret = cufftExecC2R(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_)
      ELSE
        iret = cufftExecZ2Z(plan%icuplanc_, plan%cu_ccd_, plan%cu_ccd_, FFTCU_COMPLEX_TO_REAL)
      ENDIF

      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStart(hmem);
      iret = cudaMemCpyDev2Host(plan%pccarr_, plan%cu_ccd_, plan%szccd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)

      CALL GTStart(htra,GT_WTIME);
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
                 !ReCALL that ccarr is dimensioned (:,:,:), starting at (1,1,1):
                  c1(i,j,k) = plan%ccarr(k,j,i-ista+1)
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

      CALL GTStop(htra); tratime = tratime + GTGetTime(htra)
!
! 2D FFT in each node using the FFTCU library
!
      CALL GTStart(hmem);
      iret = cudaMemCpyHost2Dev(plan%cu_cd_, plan%pcarr_, plan%szcd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)
      CALL GTStart(hfft);
      IF ( GP.EQ. 4 ) THEN
      iret = cufftExecC2R(plan%icuplanr_, plan%cu_cd_, plan%cu_rd_)
      ELSE
      iret = cufftExecZ2D(plan%icuplanr_, plan%cu_cd_, plan%cu_rd_)
      ENDIF
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStart(hmem)
      iret = cudaMemCpyDev2Host(plan%prarr_, plan%cu_rd_, plan%szrd_ )
      CALL GTStop(hmem); memtime = memtime + GTGetTime(hmem)
      out = plan%rarr

      CALL  GTFree(hfft); CALL GTFree(htra); CALL GTFree(hmem)

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
