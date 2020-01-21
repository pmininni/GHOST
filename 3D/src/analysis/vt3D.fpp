!=================================================================
      PROGRAM VT3D
!=================================================================
! VT3D code (part of the GHOST suite)
!
! Reads velocity binaries and computes all or specified
! components of velocity gradient tensor in Fourier space, and
! then outputs these to disk in real space. This tool ONLY works
! with cubic data in (2.pi)^3 domains.
!
! 2011 D. Rosenberg
!      NCAR
!
! 17 Sep 2011: Initial version
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
! Arrays for the fields

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vc
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: dvc


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: rv
      REAL(KIND=GP)                                 :: tmp,rmp,rmq,rms
!
! Auxiliary variables

      INTEGER :: n
      INTEGER :: i,ic,iir,ind,ir,it,j,jc,jjc,k
      INTEGER :: irow(3),jcol(3),istat(1024), nrow,ncol,nstat

      TYPE(IOPLAN) :: planio
      CHARACTER(len=8)   :: pref
      CHARACTER(len=256) :: odir,idir
      CHARACTER(len=256) :: fout
      CHARACTER(len=100 ):: srow, scol
      CHARACTER(len=4096):: stat
!
      NAMELIST / vt / idir, odir, srow, scol, stat, iswap, oswap

!
! Verifies proper compilation of the tool

      IF ( (nx.ne.ny).or.(ny.ne.nz) ) THEN
        IF (myrank.eq.0) &
           PRINT *,'This tool only works with cubic data in (2.pi)^3 domains'
        STOP
      ENDIF
      n = nx
!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/n,n,n/),ksta,kend,planio)
      idir   = '.'
      odir   = '.'
      srow   = '1'
      scol   = '1'
      stat  = '0'
      iswap  = 0
      oswap  = 0
      irow   = 0
      jcol   = 0
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     stat  : time index for which to compute VT, or a ';--separated list
!     srow   : list of rows to do (';' separated)
!     scol   : list of columns to do (';' separated)
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?

      IF (myrank.eq.0) THEN
         OPEN(1,file='vt.txt',status='unknown',form="formatted")
         READ(1,NML=vt)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(srow  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(scol  ,100 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

      pref = 'VT'
!
     
      ALLOCATE( vc(n,n,ista:iend) )
      ALLOCATE( dvc(n,n,ista:iend) )
      ALLOCATE( kx(n), ky(n), kz(n) )
      ALLOCATE( kn2(n,n,ista:iend), kk2(n,n,ista:iend) )
      ALLOCATE( rv(n,n,ksta:kend) )
!

      CALL fftp3d_create_plan(planrc,(/n,n,n/),FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,(/n,n,n/),FFTW_COMPLEX_TO_REAL, &
          FFTW_MEASURE)

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,n/2
         kx(i) = real(i-1,kind=GP)
         kx(i+n/2) = real(i-n/2-1,kind=GP)
      END DO
      ky = kx
      kz = kx
      rmp = 1.0_GP/real(nx,kind=GP)**2
      rmq = 1.0_GP/real(ny,kind=GP)**2
      rms = 1.0_GP/real(nz,kind=GP)**2
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               kn2(k,j,i) = rmp*kx(i)**2+rmq*ky(j)**2+rms*kz(k)**2
            END DO
         END DO
      END DO
      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
            END DO
         END DO
      END DO

      CALL parseind(srow, ';', irow , 3   , nrow) 
      CALL parseind(scol, ';', jcol , 3   , ncol) 
      CALL parseind(stat,';', istat , 1024, nstat) 

      tmp = 1./REAL(n,KIND=GP)**3
      DO it = 1,nstat
        WRITE(ext, fmtext) istat(it)
        DO ir = 1, nrow
          iir = irow(ir)
! read in appropriate file:
          IF ( iir.EQ.1 ) THEN
            CALL io_read(1,idir,'vx',ext,planio,rv)
          ELSE IF ( iir.EQ.2 ) THEN
            CALL io_read(1,idir,'vy',ext,planio,rv)
          ELSE IF ( iir.EQ.3 ) THEN
            CALL io_read(1,idir,'vz',ext,planio,rv)
          ELSE
            WRITE(*,*) 'main: invalid row specified: ',i
            STOP 
          ENDIF
!
! Byte-swap on input:
          IF ( iswap .NE. 0 ) THEN
            CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
!
! take FFT of component:
          CALL fftp3d_real_to_complex(planrc,rv,vc,MPI_COMM_WORLD)
          DO jc = 1, ncol
            jjc = jcol(jc)
            IF ( jjc.LT.1 .OR. jjc.GT.3 ) THEN
              WRITE(*,*) 'main: invalid column specified: ',j
              STOP 
            ENDIF
            CALL derivk3(vc, dvc, jjc)
            DO i = ista,iend
              DO j = 1,n
                DO k = 1,n
                  dvc(k,j,i) = dvc(k,j,i)*tmp
                END DO
              END DO
            END DO
            CALL fftp3d_complex_to_real(plancr,dvc,rv,MPI_COMM_WORLD)
!
! Byte-swap on output:
            IF ( oswap .NE. 0 ) THEN
              CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
            ENDIF
            WRITE(fout,'(a2,i01,i01)') trim(pref),iir,jjc
            fout = trim(odir) // '/' // trim(fout) // '.' // ext // '.out' 
            bmangle = 0
            CALL io_write(1,odir,fout,0,planio,rv)
            bmangle = 1
            IF ( myrank.EQ.0 ) THEN
              WRITE(*,*)'main: ', trim(fout), ' written.'
            ENDIF
          ENDDO ! tensor column loop
        ENDDO   ! tensor row loop
      ENDDO     ! time index loop
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)


      DEALLOCATE ( vc,dvc)
      DEALLOCATE ( rv)
      DEALLOCATE ( kx,ky,kz)
      DEALLOCATE ( kn2,kk2)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM VT3D

