!=================================================================
      PROGRAM BOOTS3D
!=================================================================
! BOOTS3D code (part of the GHOST suite)
!
! Reads dataset, and computes dataset that is prolongated in 
! wavenumber space, and outputs to real space. The files
! to be read are specified on input, as is the max linear size of the
! old grid; the output files are the originial files appended with with the 
! prolongation size 'PXXX'. This is the 'boostrap regridding'
! procedure (BOOTS).
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2010 D. Rosenberg
!      NCAR
!
! 17 Nov 2010: Initial version
!=================================================================

!
! Definitions for conditional compilation

! Modules

      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE iompi
      USE grid
      USE fft
      USE ali
      USE var
      USE kes
      USE fftplans
      USE threads
      USE gutils
      IMPLICIT NONE

!
! Arrays for the fields and structure functions

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1t
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: B1


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vvt
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: br
      REAL(KIND=GP)                                 :: fact
!
! Auxiliary variables

      INTEGER :: nfiles,nmt,np,npt,nt,ntprocs,npkeep
      INTEGER :: i,ib,ie,ind,iswap,itsta,itend,j,k,ktsta,ktend
      INTEGER :: istak,iendk,kstak,kendk
      INTEGER :: mykrank,mytrank

      INTEGER :: commtrunc, fh, groupworld, flags, grouptrunc, iExclude(3,1), iInclude(3,1)

      TYPE(IOPLAN)  :: planio, planiot
      TYPE(FFTPLAN) :: planrct

      CHARACTER(len=8)   :: suff
      CHARACTER(len=100) :: odir,idir
      CHARACTER(len=256) :: fname, fout, msg
      CHARACTER(len=1024):: fnlist
!
      NAMELIST / regrid / idir, odir, fnlist, iswap, nt

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)

!
! Reads from the external file 'boots.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     fnlist : file list to prolongate, separated by ';'
!     iswap  : do endian swap?
!     nt     : original linear size of the old grid

      IF (myrank.eq.0) THEN
         OPEN(1,file='boots.txt',status='unknown',form="formatted")
         READ(1,NML=regrid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(fnlist,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nt    ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      IF ( nt .GT. n .OR. nt .LT. 1 ) THEN
        STOP 'MAIN: prolongation specification incorrect; input nt must be less than N' 
      ENDIF

      WRITE(suff,'(a2,i5.5)') '_P', n
!
!  n, here refers to the prolongated grid, and nt
!  is the grid size of the 'old grid'
     
! Create new communicator for the truncated data:
      np      = n / nprocs
      ntprocs = nt / np
      nmt     = mod(nt,np) 
      IF ( nmt .GT. 0 ) ntprocs = ntprocs + 1
      ntprocs = min(ntprocs, nprocs)
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD, groupworld, ierr)
      commtrunc  = MPI_COMM_NULL
      grouptrunc = MPI_GROUP_NULL
      IF ( ntprocs .LT. nprocs ) THEN
        iExclude(1,1) = ntprocs
        iExclude(2,1) = nprocs-1
        iExclude(3,1) = 1
        CALL MPI_GROUP_RANGE_EXCL(groupworld, 1, iExclude, grouptrunc, ierr)   
        CALL MPI_COMM_CREATE(MPI_COMM_WORLD, grouptrunc, commtrunc, ierr)
      ELSE IF ( ntprocs .EQ. nprocs ) THEN
        CALL MPI_COMM_DUP(MPI_COMM_WORLD,commtrunc,ierr)
        CALL MPI_COMM_GROUP(MPI_COMM_WORLD,grouptrunc,ierr)
      ENDIF
      CALL trrange(1,n    ,nt    ,nprocs,myrank,ktsta,ktend)
      CALL trrange(1,n/2+1,nt/2+1,nprocs,myrank,itsta,itend)


      ALLOCATE( vvt(nt,nt,ktsta:ktend) )
      ALLOCATE( C1t(nt,nt,itsta:itend) )
      ALLOCATE( B1(n,n,ista:iend) )
      ALLOCATE( ka(nt),ka2(nt,nt,itsta:itend) )
      ALLOCATE( br(n,n,ksta:kend) )

      IF ( myrank .LT. ntprocs ) THEN
        npkeep = nprocs; nprocs = ntprocs
        CALL io_init(myrank,nt,ktsta,ktend,planiot)
        nprocs = npkeep
      ENDIF
!
      flags  = FFTW_MEASURE

      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, flags)
      CALL trrange(1,n    ,nt    ,nprocs,myrank,ksta,kend)
      CALL trrange(1,n/2+1,nt/2+1,nprocs,myrank,ista,iend)
write(*,*)'main: creating trplan...'
      CALL fftp3d_create_trplan(planrct,n,nt,FFTW_REAL_TO_COMPLEX,flags)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)

#if 0
      IF ( myrank .LT. ntprocs ) THEN
!       CALL MPI_COMM_RANK(commtrunc,myrank,ierr)
        npkeep = nprocs ; istak  = ista ; iendk  = iend ; kstak  = ksta ; kendk  = kend ;
        nprocs = ntprocs; ista   = itsta; iend   = itend; ksta   = ktsta; kend   = ktend;
        IF ( ierr .NE. MPI_SUCCESS ) THEN
          WRITE(*,*) 'main: MPI_COMM_RANK failed on commtrunc; myrank= ', myrank
          STOP
        ENDIF
        CALL fftp3d_create_plan(planrct,nt, FFTW_REAL_TO_COMPLEX, FFTW_MEASURE)
        nprocs = npkeep ; ista   = istak; iend   = iendk; ksta   = kstak; kend   = kendk;
!       CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      ENDIF
#endif
!
      kmax = (real(nt,kind=GP)/3.)**2
      DO i = 1,nt/2
         ka(i) = real(i-1,kind=GP)
         ka(i+nt/2) = real(i-nt/2-1,kind=GP)
      END DO

      DO i = itsta,itend
         DO j = 1,nt
            DO k = 1,nt
               ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
            END DO
         END DO
      END DO

     
      ib = 1
      ie = len(fnlist)
      DO WHILE ( len(trim(fnlist(ib:ie))) .GT. 0 )  

         ind = index(fnlist(ib:ie),";")
         IF ( ind .eq. 0 ) THEN
           fname = trim(adjustl(fnlist(ib:ie)))
           ib = ie + 1
         ELSE
           fname = trim(adjustl(fnlist(ib:(ib+ind-2))))
           ib = ib + ind
         ENDIF
!
! Read the external binary files in real space:
         IF ( myrank .LT. ntprocs ) THEN
           CALL MPI_FILE_OPEN(commtrunc,trim(idir) // '/' // trim(fname), &
                 MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error opening file: ', trim(idir) // '/' // trim(fname)
             STOP
           ENDIF
           CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,planiot%iotype,'native', &
               MPI_INFO_NULL,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error with MPI_FILE_SET_VIEW'
             STOP
           ENDIF
           CALL MPI_FILE_READ_ALL(fh,vvt, &
               planiot%n*planiot%n*(planiot%kend-planiot%ksta+1),GC_REAL, &
               MPI_STATUS_IGNORE,ioerr)
           IF ( ioerr .NE. MPI_SUCCESS ) THEN
             WRITE(*,*) 'main: Error with MPI_FILE_READ_ALL'
             STOP
           ENDIF
           CALL MPI_FILE_CLOSE(fh,ioerr)

           IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(vvt, nt)
           ENDIF

         ENDIF

!
! Compute FT of variable:
!        IF ( myrank .LT. ntprocs ) THEN
#if 0
           npkeep = nprocs ; istak  = ista ; iendk  = iend ; kstak  = ksta ; kendk  = kend ;
           nprocs = ntprocs; ista   = itsta; iend   = itend; ksta   = ktsta; kend   = ktend;
           CALL MPI_COMM_RANK(commtrunc,myrank,ierr)
           CALL fftp3d_real_to_complex(planrct,vvt,C1t,MPI_COMM_WORLD)
           nprocs = npkeep ; ista   = istak; iend   = iendk; ksta   = kstak; kend   = kendk;
           CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
#endif
if ( myrank.lt.ntprocs ) write(*,*)'main: C1t=', C1t(1:10,1,1)
if ( myrank.lt.ntprocs ) write(*,*)'main: vvt=', vvt(1,1,1:10)
           CALL trrange(1,n    ,nt    ,nprocs,myrank,ksta,kend)
           CALL trrange(1,n/2+1,nt/2+1,nprocs,myrank,ista,iend)
           ista   = itsta; iend   = itend; ksta   = ktsta; kend   = ktend;
           CALL fftp3d_real_to_complex(planrct,vvt,C1t,MPI_COMM_WORLD)
           CALL range(1,n/2+1,nprocs,myrank,ista,iend)
           CALL range(1,n,nprocs,myrank,ksta,kend)

!        ENDIF
!
! Prolongate in Fourier space:
!
!
         B1 = 0.0
! Compute inverse FT of truncated variable, and store in
! larger storage for output:
         fact = 1.0_GP / real(nt,kind=GP)**3
         IF ( myrank .LT. ntprocs ) THEN
           DO i = itsta,itend
              DO j = 1,nt/2
                 DO k = 1,nt/2
                    B1(k,j,i) = C1t(k,j,i) * fact
                 END DO
                 DO k = nt/2+1,nt
                    B1(k,j,i) = C1t(k-nt+n,j,i) * fact
                 END DO
              END DO
              DO j = nt/2+1,nt
                 DO k = 1,nt/2
                    B1(k,j,i) = C1t(k,j-nt+n,i) * fact
                 END DO
                 DO k = nt/2+1,nt
                    B1(k,j,i) = C1t(k-nt+n,j-nt+n,i) * fact
                 END DO
              END DO
           END DO

         ENDIF
!
! Compute inverse FT of prolongated variable:
         CALL fftp3d_complex_to_real(plancr,B1,br,MPI_COMM_WORLD)
!
! Put to disk:
         fout = trim(odir) // '/' // trim(fname) // trim(suff)
         bmangle = 0
         CALL io_write(1,idir,fout,1,planio,br)
         bmangle = 1
#if 0
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD, fout, &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, &
           MPI_INFO_NULL,fh,ioerr)
         CALL MPI_FILE_SET_VIEW(fh,disp,GC_REAL,planio%iotype,'native', &
           MPI_INFO_NULL,ioerr)
         CALL MPI_FILE_WRITE_ALL(fh,br, &
           planio%n*planio%n*(planio%kend-planio%ksta+1),GC_REAL, &
           MPI_STATUS_IGNORE,ioerr)
         CALL MPI_FILE_CLOSE(fh,ioerr)
#endif
         IF ( myrank .EQ. 0 ) THEN
         WRITE(*,*) 'main: ',trim(fout),' written.'
         ENDIF
      ENDDO
!
      CALL fftp3d_destroy_plan(plancr)
      IF ( myrank .LT. ntprocs ) THEN
        CALL fftp3d_destroy_plan(planrct)
        CALL MPI_COMM_FREE(commtrunc, ierr)
        CALL MPI_GROUP_FREE(grouptrunc, ierr)
      ENDIF


      DEALLOCATE ( vvt )
      DEALLOCATE ( C1t )
      DEALLOCATE ( B1 )
      DEALLOCATE ( ka)
      DEALLOCATE ( ka2)
      DEALLOCATE ( br )

      CALL MPI_FINALIZE(ierr)

      END PROGRAM BOOTS3D


      SUBROUTINE rarray_byte_swap( Rin, nin )
!-----------------------------------------------------------------
!
! Performs endian conversion on array of floats
!
! Parameters
!     Rin : 3D array whose values will be endian-swapped, returned
!     nin : linear dimension of rin
!-----------------------------------------------------------------
      USE mpivars
      USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(IN)                                        :: nin
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nin,nin,ksta:kend) :: Rin

      INTEGER :: i, ie0, ie1, j, k, m, nb

      nb = 8

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,nin
            DO i = 1,nin
               ie0 = TRANSFER(Rin(i,j,k), 0)
               DO m = 1, GP
                  CALL MVBITS( ie0, (GP-m)*nb, nb, ie1, (m-1)*nb  )
               END DO
               Rin(i,j,k) = TRANSFER(ie1, 0.0)
            END DO
         END DO
       END DO
      
      RETURN

      END SUBROUTINE

