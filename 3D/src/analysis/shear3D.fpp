!=================================================================
      PROGRAM SHEAR3D
!=================================================================
! SHEAR3D code (part of the GHOST suite)
!
! Reads velocity binaries and computes all required
! components of strain rate tensor in Fourier space, computes
! eigenvalues, and isotropic spectra of the eigenvallue fields.
! Note: this utility is _strictly_ for incompressible flows!
!
!
! 2013 D. Rosenberg
!      ORNL/NCCS
!
! 28 May 2013: Initial version
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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ctmp,sij


      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: lamb,R1,R2,R3,R4,R5
      REAL   (KIND=GP)                                 :: btrunc,ktrunc,tmp
!
! Auxiliary variables

      INTEGER :: i,ic,iir,ind,ir,iswap,it,j,jc,jjc,k
      INTEGER :: inorm,istat(1024), nstat

      TYPE(IOPLAN) :: planio
      CHARACTER(len=10)   :: suff
      CHARACTER(len=128)  :: pref
      CHARACTER(len=256)  :: odir,idir
      CHARACTER(len=1024) :: fnout
      CHARACTER(len=2048) :: fntmp
      CHARACTER(len=4096) :: stat
!
      NAMELIST / shear / btrunc, ktrunc, iswap, idir, odir, pref, stat

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)
     
      kmax   = real(n/2+1,kind=GP)

      idir   = '.'
      odir   = '.'
      stat  = '0'
      iswap  = 0
      btrunc = 0
      ktrunc = kmax
      pref = 'ksplambda'
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     stat  : time index for which to compute SHEAR, or a ';--separated list
!     ktrunc: if set to  < n/2+1, will truncate strain in computing spectrum.
!     iswap  : do endian swap on input?

      IF (myrank.eq.0) THEN
         OPEN(1,file='shear.txt',status='unknown',form="formatted")
         READ(1,NML=shear)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pref  ,128 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat  ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(btrunc,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktrunc,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)

!

      ALLOCATE( ctmp(n,n,ista:iend) )
      ALLOCATE( sij (n,n,ista:iend) )
      ALLOCATE( vx(n,n,ista:iend) )
      ALLOCATE( vy(n,n,ista:iend) )
      ALLOCATE( vz(n,n,ista:iend) )
      ALLOCATE( ka(n) )
      ALLOCATE( ka2(n,n,ista:iend) )
      ALLOCATE( R1(n,n,ksta:kend) )
      ALLOCATE( R2(n,n,ksta:kend) )
      ALLOCATE( R3(n,n,ksta:kend) )
      ALLOCATE( R4(n,n,ksta:kend) )
      ALLOCATE( R5(n,n,ksta:kend) )
      ALLOCATE( lamb(n,n,ksta:kend) )
!

      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
          FFTW_MEASURE)
!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      tiny = 1e-5
      if ( btrunc.eq.0 ) ktrunc = kmax

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,n/2
         ka(i) = REAL(i-1,KIND=GP)
         ka(i+n/2) = REAL(i-n/2-1,KIND=GP)
      END DO
      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
            END DO
         END DO
      END DO

      CALL parseind(stat,';', istat , 1024, nstat) 

      tmp = 1.0_GP/REAL(n,KIND=GP)**3
      DO it = 1,nstat
        WRITE(ext, fmtext) istat(it)
! read in appropriate file:
        CALL io_read(1,idir,'vx',ext,planio,R1)
        CALL io_read(1,idir,'vy',ext,planio,R2)
        CALL io_read(1,idir,'vz',ext,planio,R3)
!
! Byte-swap on input:
        IF ( iswap .NE. 0 ) THEN
          CALL rarray_byte_swap(R1, n*n*(kend-ksta+1))
          CALL rarray_byte_swap(R2, n*n*(kend-ksta+1))
          CALL rarray_byte_swap(R3, n*n*(kend-ksta+1))
        ENDIF
!
! take FFT of component:
        CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
        CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
        CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)

! Compute required strain rate components:
        inorm = 1
        CALL Strain(vx,vy,vz,1,1,ktrunc,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R1,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,1,2,ktrunc,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R2,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,1,3,ktrunc,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R3,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,2,2,ktrunc,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R4,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,2,3,ktrunc,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R5,MPI_COMM_WORLD)

! Compute required eigenvalue field of strain:
        CALL Eigen(R1,R2,R3,R4,R5,lamb)
        CALL fftp3d_real_to_complex(planrc,lamb,ctmp,MPI_COMM_WORLD)
!
! Compute power spectrum and output it:
        IF ( btrunc.EQ.0 ) THEN
          fnout = trim(pref) // '.' // ext  // '.txt';
        ELSE
          WRITE(suff,'(a2,i5.5)') '_T', int(ktrunc)
          fnout = trim(pref) // '.' // ext //  trim(suff) // '.txt'
        ENDIF
        fntmp = trim(odir) // '/' // trim(fnout)

        IF ( myrank.EQ. 0 ) THEN
          write(*,*)'main: fntmp=',trim(fntmp),' ktrunc=',ktrunc
          write(*,*)'main: time index ', trim(ext), ' done.'
        ENDIF
        CALL pspectrum(ctmp,fntmp,int(ktrunc))

      ENDDO   ! time (stat) loop
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)
      CALL MPI_FINALIZE(ierr)

      IF ( ALLOCATED(ctmp) ) DEALLOCATE(ctmp)
      IF ( ALLOCATED (sij) ) DEALLOCATE (sij)
      IF ( ALLOCATED  (vx) ) DEALLOCATE  (vx)
      IF ( ALLOCATED  (vy) ) DEALLOCATE  (vy)
      IF ( ALLOCATED  (vz) ) DEALLOCATE  (vz)
      IF ( ALLOCATED(lamb) ) DEALLOCATE(lamb)
      IF ( ALLOCATED  (R1) ) DEALLOCATE  (R1)
      IF ( ALLOCATED  (R2) ) DEALLOCATE  (R2)
      IF ( ALLOCATED  (R3) ) DEALLOCATE  (R3)
      IF ( ALLOCATED  (R4) ) DEALLOCATE  (R4)
      IF ( ALLOCATED  (R5) ) DEALLOCATE  (R5)
      IF ( ALLOCATED  (ka) ) DEALLOCATE  (ka)
!!    IF ( ALLOCATED (ka2) ) DEALLOCATE (ka2)
!!write(*,*)'main: after ka2'

      END PROGRAM SHEAR3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE Strain(vx,vy,vz,ir,jc,ktrunc,inorm,sij,ctmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the complex strain rate component 
!
! Parameters
!     vi   : input velocities
!     sij  : returned complex component of strain rate tensor 
!     ir,jc: the row and col of sij
!     krtunc: truncaton wavenumber for spherical truncation
!     inorm : normalize (1), or not (0)
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: ctmp,sij
      REAL   (KIND=GP), INTENT   (IN)                           :: ktrunc
      REAL   (KIND=GP)                                          :: ktrunc2,tmp
!
      INTEGER         , INTENT   (IN)                           :: inorm,ir,jc
      INTEGER                                                   :: i,j,k

      IF ( ir.NE.1 .AND. ir.NE.2 .AND. ir.NE.3 &
      .AND.jc.NE.1 .AND. jc.NE.2 .AND. jc.NE.3 ) THEN
        WRITE(*,*)'Strain: invalid row/column specification: ', ir, jc
        STOP
      ENDIF

      ktrunc2 = ktrunc**2

      IF ( ir.EQ.1 ) THEN
        CALL derivk3(vx, sij, jc)
        SELECT CASE (jc)
          CASE(1)
            DO i = ista,iend
              DO j = 1,n
                DO k = 1,n
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
          CASE(2)
            CALL derivk3(vy, ctmp, 1)
          CASE(3)
            CALL derivk3(vz, ctmp, 1)
        END SELECT
      ELSE IF ( ir.EQ.2 ) THEN
        CALL derivk3(vy, sij, jc)
        SELECT CASE (jc)
          CASE(1)
            CALL derivk3(vx, ctmp, 2)
          CASE(2)
            DO i = ista,iend
              DO j = 1,n
                DO k = 1,n
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
          CASE(3)
            CALL derivk3(vz, ctmp, 2)
        END SELECT
      ELSE IF ( ir.EQ.3 ) THEN
        CALL derivk3(vz, sij, jc)
        SELECT CASE (jc)
          CASE(1)
            CALL derivk3(vx, ctmp, 3)
          CASE(2)
            CALL derivk3(vy, ctmp, 3)
          CASE(3)
            DO i = ista,iend
              DO j = 1,n
                DO k = 1,n
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
        END SELECT
      ENDIF

      DO i = ista,iend
        DO j = 1,n
          DO k = 1,n
            sij(k,j,i) = 0.50_GP*(sij(k,j,i)+ctmp(k,j,i)) 
          END DO
        END DO
      END DO


      ! truncate spherically:
      DO i = ista,iend
        DO j = 1,n
          DO k = 1,n
            IF ((ka2(k,j,i).gt.ktrunc2 ).and.(ka2(k,j,i).ge.tiny)) THEN
              sij(k,j,i) = 0.0_GP
            ENDIF
          END DO
        END DO
      END DO


      IF ( inorm.GT.0 ) THEN
        
        tmp = 1.0_GP/REAL(n,KIND=GP)**3
        DO i = ista,iend
          DO j = 1,n
            DO k = 1,n
              sij(k,j,i) = sij(k,j,i)*tmp
            END DO
          END DO
        END DO

      ENDIF

      END SUBROUTINE Strain
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE Eigen(S11,S12,S13,S22,S23,lamb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the eigenvalue field for strain rate tensor, whos
! required components are specified
!
! Parameters
!     SIJ  : strain rate components (real)
!     lamb : returned eigenvalue field
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: S11,S12,S13,S22,S23
      REAL   (KIND=GP), INTENT  (OUT), DIMENSION(n,n,ksta:kend) :: lamb
      REAL   (KIND=GP)                                          :: sa,sb,sc,sd,se,sf
      REAL   (KIND=GP)                                          :: a,b,c,d
      REAL   (KIND=GP)                                          :: del0,del1,ll(3),lmax
      COMPLEX(KIND=GP)                                          :: CC,u(3)
      COMPLEX(KIND=GP)                                          :: D0,D1

!
      INTEGER                                                   :: i,j,k,l

      u(1) = cmplx(1.0_GP , 0.0_GP)
      u(2) = cmplx(-0.5_GP, 0.5_GP*sqrt(3.0_GP))
      u(3) = cmplx(-0.5_GP,-0.5_GP*sqrt(3.0_GP))
      DO k = ksta,kend
        DO j = 1,n
          DO i = 1,n
            sa = S11(i,j,k); sb = S12(i,j,k); sc = S13(i,j,k); 
            sd = S22(i,j,k); se = S23(i,j,k); sf = -(sa + sd);
            a    = 1.0_GP
            b    = 0.0_GP ! sa + sd + sf
            c    = sa*sd+sa*sf+sd*sf-sb**2 - sc**2 - se**2
            d    = sa*se**2 + sf*sb**2 + sd*sc**2 - sa*sd*sf - 2.0_GP*sb*sc*se
!           del0 = b**2 - 3.0_GP*a*c
!           del1 = 2.0_GP*b**3 - 9.0_GP*a*b*c + 27.0_GP*d*a**2
            del0 = -3.0_GP*a*c
            del1 =  27.0_GP*d*a**2
            D0   = cmplx(del0,0)
            D1   = cmplx(del1,0)
            CC   = (0.5*D1 + 0.5*sqrt(D1**2 -4.0_GP*D0**3) )**(1.0_GP/3.0_GP)
            lmax = 0.0_GP
            DO l = 1,3
!             ll(l) = real(-( b + u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
              ll(l) = real(-( u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
              lmax = max(lmax,abs(ll(l)))
            ENDDO
!if ( i.eq.10 .and. j.eq.10 .and. k.gt.10 .and. k.lt.15) then
!!write(*,*)'Eigen: sa=',sa,' sb=',sb,' sc=',sc,' sd=',sd, 'se=',se
!write(*,*)'Eigen: u=',u,' CC=', CC,' del0=',del0,' del1=',del1,' lmax=',lmax, ' ll=',ll
!endif

            lamb(i,j,k) = lmax
          END DO
        END DO
      END DO


      END SUBROUTINE Eigen

