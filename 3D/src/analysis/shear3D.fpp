      PROGRAM SHEAR3D
!=================================================================
! SHEAR3D code (part of the GHOST suite)
!
! Reads velocity binaries and computes all required
! components of strain rate tensor in Fourier space, computes
! eigenvalues, and isotropic spectra of the eigenvallue fields.
! Note: this utility is _strictly_ for incompressible flows!
!
! Additionally, if isolve > 0 , the utility will also compute the
! eigenvector at each spatial location corresponding to the max
! |eigenvalue|. This gives what might be termed the 'principle-principle'
! axis representing the direction of the largest amount of 'shear'
! This computation relied heavily on the following article:
!
!  'Eigensystems for 3x3 Symmetric Matrices (Revisited)'
!   David Eberly, Geometric Tools, LLC
!   http://www.geometrictools.com
!   (c) 1998-2012, May 2011 
!
! Finally, if jpdf=1, jpoint and 1d pdfs of energy dissipation and 
! entrophy density, energy diss and lambda, and energy diss and (relative?)
! helicity will be computed and written.
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
      USE random
      USE threads
      USE gutils
      IMPLICIT NONE

!
! Arrays for the fields and structure functions

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: ctmp,sij
      COMPLEX(KIND=GP)                                 :: cdump,jdump



      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: lamb,R1,R2,R3,R4,R5
      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: evx,evy,evz
      REAL   (KIND=GP)                                 :: btrunc,sg,sl,tmp
      REAL   (KIND=GP)                                 :: ensmin,ensmax,dismin,dismax
      REAL   (KIND=GP)                                 :: krmin,krmax,phase
      REAL   (KIND=GP)                                 :: krmin2,krmax2 
      REAL   (KIND=GP)                                 :: ktmin ,ktmax 
      REAL   (KIND=GP)                                 :: ktmin2,ktmax2

!
! Auxiliary variables

      INTEGER :: i,ic,iir,ind,ir,isolve,iswap,it,j,jc,jjc,k,oswap
      INTEGER :: demean,dolog,ilamb,inorm,istat(1024),jpdf,nstat
      INTEGER :: irand,nbinx,nbiny,nbins(2),seed
!$    INTEGER, EXTERNAL :: omp_get_max_threads

      TYPE(IOPLAN) :: planio
      CHARACTER(len=16)   :: suff
      CHARACTER(len=128)  :: pref
      CHARACTER(len=256)  :: odir,idir
      CHARACTER(len=1024) :: fnout
      CHARACTER(len=2048) :: fntmp
      CHARACTER(len=2048) :: fntmp1,fntmp2,fntmp3,fntmp4,fntmp5
      CHARACTER(len=4096) :: stat
!
      NAMELIST / shear / demean,ilamb,isolve,iswap
      NAMELIST / shear / dolog,oswap,idir,odir,pref,stat
      NAMELIST / shear / dismin,dismax,ensmin,ensmax,jpdf,nbinx,nbiny
      NAMELIST / shear / irand,krmin,krmax,seed
      NAMELIST / shear / btrunc,ktmin,ktmax

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)
      nth = 1
!$    nth = omp_get_max_threads()
!$    CALL fftp3d_init_threads(ierr)
!     tiny: minimum truncation for dealiasing
      tiny = 1e-5
     
      kmax   = real(n/2+1,kind=GP)

      idir   = '.'
      odir   = '.'
      stat  = '0'
      iswap  = 0
      oswap  = 0
      btrunc = 0
      demean = 0
      dolog  = 1
      ilamb  = 0
      irand  = 0
      isolve = 0
      jpdf   = 1
      krmin  = tiny
      krmax  = kmax
      ktmin  = tiny
      ktmax  = kmax
      seed   = 1000
      pref = 'ksplambda'
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     stat   : time index for which to compute SHEAR, or a ';--separated list
!     btrunc : if == 1, truncate spectral range to [ktmin,ktmax]
!     ktmin  : min wavenumber for truncation if btrunc=1
!     ktmax  : max wavenumber for truncation if btrunc=1
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     isolve : 0==>just max |eigenvale| field; 1==>max |eigenvalue| field plus
!              corresponding eigenvector inevx,evy,evz
!     ilamb  : 1==>write eigenvalue field to disk; 0==>don't
!     irand  : randomize phases between [krmin,krmax] if 1; else, don't
!     krmin  : min wavenumber for randomization if irand=1
!     krmax  : max wavenumber for randomization if irand=1
!     jpdf   : 1==>do joint pdf of energy diss and other things; 0==>don't
!     demean : demean the eigenvalue field?
!     dolog  : compute PDFs in log=space?

      IF (myrank.eq.0) THEN
write(*,*)'main: opening shear.txt...'
         OPEN(1,file='shear.txt',status='unknown',form="formatted")
         READ(1,NML=shear)
         CLOSE(1)
write(*,*)'main: shear.txt read.'
      ENDIF
      CALL MPI_BCAST(idir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pref  ,128 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat  ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(btrunc,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(demean,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(dolog ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(isolve,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilamb ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(irand ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(jpdf  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktmin ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ktmax ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbinx ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nbiny ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(krmin ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(krmax ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
if (myrank.eq.0) write(*,*)'main: broadcast done.'
!
      nbins(1) = nbinx
      nbins(2) = nbiny
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
      
      ALLOCATE( evx(n,n,ksta:kend) )
      IF ( isolve.GT.0 ) THEN
      ALLOCATE( evy(n,n,ksta:kend) )
      ALLOCATE( evz(n,n,ksta:kend) )
      ENDIF
!

if (myrank.eq.0) write(*,*)'main: creating plans...'
      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
          FFTW_MEASURE)
if (myrank.eq.0) write(*,*)'main: plans done.'
!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
      IF ( irand.eq.0 ) THEN
         krmin  = tiny
         krmax  = kmax
      ENDIF
      IF ( btrunc.eq.0 ) THEN
         ktmin  = tiny
         ktmax  = kmax
      ENDIF
      krmin2 = krmin**2
      krmax2 = krmax**2
      ktmin2 = ktmin**2
      ktmax2 = ktmax**2
if (myrank.eq.0) write(*,*)'main: d-extrema done.'

!
! Builds the wave number and the square wave 
! number matrixes

      DO i = 1,n/2
         ka(i) = REAL(i-1,KIND=GP)
         ka(i+n/2) = REAL(i-n/2-1,KIND=GP)
      END DO
if (myrank.eq.0) write(*,*)'main: ka done.'

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
            END DO
         END DO
      END DO
if (myrank.eq.0) write(*,*)'main: ka2 done.'

      CALL parseind(stat,';', istat , 1024, nstat) 
if (myrank.eq.0) write(*,*)'main: index parsing done: nstat=',nstat

      tmp = 1.0_GP/REAL(n,KIND=GP)**3
      DO it = 1,nstat
if (myrank.eq.0) write(*,*)'main: writing ext: fmtext:', fmtext,' istat=',istat(it)
        WRITE(ext, fmtext) istat(it)
! read in appropriate file:
if (myrank.eq.0) write(*,*)'main: reading vx...'
        CALL io_read(1,idir,'vx',ext,planio,R1)
if (myrank.eq.0) write(*,*)'main: reading vy...'
        CALL io_read(1,idir,'vy',ext,planio,R2)
if (myrank.eq.0) write(*,*)'main: reading vz...'
        CALL io_read(1,idir,'vz',ext,planio,R3)
if (myrank.eq.0) write(*,*)'main: reading done.'


!
! Byte-swap on input:
        IF ( iswap .NE. 0 ) THEN
          CALL rarray_byte_swap(R1, n*n*(kend-ksta+1))
          CALL rarray_byte_swap(R2, n*n*(kend-ksta+1))
          CALL rarray_byte_swap(R3, n*n*(kend-ksta+1))
        ENDIF
!
! take FFT of component:
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for vx...'
        CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for vy...'
        CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: real 2 cmplex for vz...'
        CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)
if (myrank.eq.0) write(*,*)'main: real 2 cmplex done.'
        IF ( irand.GT.0 ) THEN
          IF (myrank.eq.0) phase = 2*pi*randu(seed)
          CALL MPI_BCAST(phase,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
          cdump = COS(phase)+im*SIN(phase)
          jdump = conjg(cdump)
          CALL Randomize(vx,vy,vz,krmin,krmax,phase)
        ENDIF

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,n
            DO k = 1,n
              IF ((ka2(k,j,i).lt.ktmin2).or.(ka2(k,j,i).gt.ktmax2)) THEN
                vx(k,j,i) = 0.0_GP
                vy(k,j,i) = 0.0_GP
                vz(k,j,i) = 0.0_GP
              ENDIF
            END DO
          END DO
        END DO

write(*,*)'main: ktmin2=',ktmin2, ' ktmax2=',ktmax2

! Compute required strain rate components:
        inorm = 1
        CALL Strain(vx,vy,vz,1,1,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R1,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,1,2,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R2,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,1,3,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R3,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,2,2,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R4,MPI_COMM_WORLD)
        CALL Strain(vx,vy,vz,2,3,ktmin,ktmax,inorm,sij,ctmp)
        CALL fftp3d_complex_to_real(plancr,sij,R5,MPI_COMM_WORLD)

! Compute required eigenvalue field of strain:
        IF ( isolve.EQ.0 ) THEN
          CALL EigenValMax(R1,R2,R3,R4,R5,lamb)
        ELSE
          CALL EigenSolveMax(R1,R2,R3,R4,R5,lamb,evx,evy,evz)
        ENDIF
!
! Remove mean, if desired:

        IF ( demean.GT.0 ) THEN
        sl = 0.0_GP
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,n
            DO i = 1,n
!$omp atomic
              sl = sl + lamb(i,j,k)
             ENDDO
           ENDDO
         ENDDO
         call MPI_ALLREDUCE(sl,sg,1,GC_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
         sg = sg / real(n,kind=GP)**3

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
          DO j = 1,n
            DO i = 1,n
!$omp atomic
              lamb(i,j,k) = lamb(i,j,k) - sg
             ENDDO
           ENDDO
         ENDDO
         ENDIF

        CALL fftp3d_real_to_complex(planrc,lamb,ctmp,MPI_COMM_WORLD)
!
! Compute power spectrum of e-value  and output it:
        IF ( btrunc.EQ.0 ) THEN
          fnout = trim(pref) // '.' // ext  // '.txt';
        ELSE
          WRITE(suff,'(a2,i5.5,a1,i5.5)') '_T', int(ktmin),'_',int(ktmax)
          fnout = trim(pref) // '.' // ext //  trim(suff) // '.txt'
        ENDIF
        fntmp = trim(odir) // '/' // trim(fnout)

        CALL pspectrum(ctmp,fntmp,int(kmax))

! 
! Prepare eignesystem for output if necessary
! (don't forget to undo swaps for later):
        IF ( ilamb .GT. 0 ) THEN
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(lamb, n*n*(kend-ksta+1))
          ENDIF
          CALL io_write(1,odir,'lmb',ext,planio,lamb)
        ENDIF
        IF ( isolve .GT. 0 ) THEN
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(evx, n*n*(kend-ksta+1))
            CALL rarray_byte_swap(evy, n*n*(kend-ksta+1))
            CALL rarray_byte_swap(evz, n*n*(kend-ksta+1))
          ENDIF
          CALL io_write(1,odir,'evx',ext,planio,evx)
          CALL io_write(1,odir,'evy',ext,planio,evy)
          CALL io_write(1,odir,'evz',ext,planio,evz)
        ENDIF
!  
!  
! Compute and print joint pdfs for epsilon/(2\nu), and enstropohy, helicity, lambda:
        IF ( jpdf.GT.0 ) THEN
          ! Note, the strain rate components RI & lambda will be modified on exit:
          fnout  = trim(pref) // '.' // ext  // '.txt';
          fntmp1 = trim(odir) // '/'// 'jpdf_diss_enst.'// ext // '.txt'
          fntmp2 = trim(odir) // '/'// 'jpdf_diss_lamb.'// ext // '.txt'
          fntmp3 = trim(odir) // '/'// 'jpdf_enst_lamb.'// ext // '.txt'
          fntmp4 = trim(odir) // '/'// 'jpdf_diss_hel.' // ext // '.txt'
          fntmp5 = trim(odir) // '/'// 'jpdf_enst_hel.' // ext // '.txt'
          IF ( oswap .NE. 0 ) THEN
            CALL rarray_byte_swap(lamb, n*n*(kend-ksta+1))
          ENDIF
          CALL DoDissJPDF(R1,R2,R3,R4,R5,vx,vy,vz,lamb,ext    , &
               fntmp1,fntmp2,fntmp3,fntmp4,fntmp5,nbins,dolog , &
               trim(odir) // '/disspdf.' // ext // '.txt' , &
               trim(odir) // '/enstpdf.' // ext // '.txt' , & 
               trim(odir) // '/lambpdf.' // ext // '.txt' , &
               trim(odir) // '/helpdf.'  // ext // '.txt' , &
               trim(odir) // '/flat'            // '.txt' , &
               trim(odir) // '/skew'            // '.txt' , &
               ctmp,sij,evx)
        ENDIF

        IF ( myrank.EQ. 0 ) THEN
          write(*,*)'main: fntmp=',trim(fntmp),' ktmin=',ktmin,' ktmax=',ktmax
          write(*,*)'main: time index ', trim(ext), ' done.'
        ENDIF

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
      IF ( ALLOCATED (evx) ) DEALLOCATE (evx)
      IF ( ALLOCATED (evy) ) DEALLOCATE (evy)
      IF ( ALLOCATED (evz) ) DEALLOCATE (evz)
      IF ( ALLOCATED  (ka) ) DEALLOCATE  (ka)
!!    IF ( ALLOCATED (ka2) ) DEALLOCATE (ka2)

      END PROGRAM SHEAR3D
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE Strain(vx,vy,vz,ir,jc,ktmin,ktmax,inorm,sij,ctmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the complex strain rate component 
!
! Parameters
!     vi   : input velocities
!     sij  : returned complex component of strain rate tensor 
!     ir,jc: the row and col of sij
!     ktmin: truncaton min wavenumber for spherical truncation
!     ktmax: truncaton max wavenumber for spherical truncation
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
      REAL   (KIND=GP), INTENT   (IN)                           :: ktmin,ktmax
      REAL   (KIND=GP)                                          :: ktmin2,ktmax2,tmp
!
      INTEGER         , INTENT   (IN)                           :: inorm,ir,jc
      INTEGER                                                   :: i,j,k

      IF ( ir.NE.1 .AND. ir.NE.2 .AND. ir.NE.3 &
      .AND.jc.NE.1 .AND. jc.NE.2 .AND. jc.NE.3 ) THEN
        WRITE(*,*)'Strain: invalid row/column specification: ', ir, jc
        STOP
      ENDIF

      ktmin2 = ktmin**2
      ktmax2 = ktmax**2

      IF ( ir.EQ.1 ) THEN
        CALL derivk3(vx, sij, jc)
        SELECT CASE (jc)
          CASE(1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
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
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
              DO j = 1,n
                DO k = 1,n
                  ctmp(k,j,i) = sij(k,j,i)
                END DO
              END DO
            END DO
        END SELECT
      ENDIF

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            sij(k,j,i) = 0.50_GP*(sij(k,j,i)+ctmp(k,j,i)) 
          END DO
        END DO
      END DO


      ! truncate spherically:
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            IF ((ka2(k,j,i).lt.ktmin2 ).or.(ka2(k,j,i).gt.ktmax2)) THEN
              sij(k,j,i) = 0.0_GP
            ENDIF
          END DO
        END DO
      END DO


      IF ( inorm.GT.0 ) THEN
        
        tmp = 1.0_GP/REAL(n,KIND=GP)**3
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
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


      SUBROUTINE EigenValMax(S11,S12,S13,S22,S23,lamb)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the max |eigenvalue| field for strain rate tensor, whos
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
!$omp parallel do if (kend-ksta.ge.nth) private (j,i,l,sa,sb,sc,sd,se,sf,a,b,c,d,del0,del1,D0,D1,CC,lmax)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,l,sa,sb,sc,sd,se,sf,a,b,c,d,del0,del1,D0,D1,CC,lmax)
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


      END SUBROUTINE EigenValMax
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE EigenSolveMax(S11,S12,S13,S22,S23,lamb,evx,evy,evz)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the max |eigenvalue| field for strain rate tensor, whos
! required components are specified. Also computes the compnents of
! the eigenvector corresponding to this max eigenvalue at each point.
! If at a point the strain rate tensor is not of rank 2, meaning 3
! distinct e-values, and 3 e-vectors, then return (0,0,0) for the 
! e-vector of the max e-value.
!
! Parameters
!     SIJ  : strain rate components (real)
!     lamb : returned eigenvalue field
!     evx,
!     evy,
!     evz  : field of eigenvector components corresp. to lamb, returned
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT   (IN), DIMENSION(n,n,ksta:kend) :: S11,S12,S13,S22,S23
      REAL   (KIND=GP), INTENT  (OUT), DIMENSION(n,n,ksta:kend) :: lamb,evx,evy,evz
      REAL   (KIND=GP)                                          :: s   (3,3),m1  (3,3),m2  (3,3),m3  (3,3)
      REAL   (KIND=GP)                                          :: eval(3)  ,evec(3)  ,row1(3)  ,row2(3)
      REAL   (KIND=GP)                                          :: ev1 (3)  ,ev2 (3)
      INTEGER                                                   :: i,j,k,l,m
      INTEGER                                                   :: GetRank
!$    INTEGER,EXTERNAL                                          :: omp_get_thread_num

!$omp parallel do if (kend-ksta.ge.nth) private (j,i,s,m1,m2,m3,evec,eval,ev1,ev2,row1,row2)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i,s,m1,m2,m3,evec,eval,ev1,ev2,row1,row2)
        DO j = 1,n
          DO i = 1,n
            s(1,1)=S11(i,j,k); s(1,2)=S12(i,j,k); s(1,3)=S13(i,j,k)
            s(2,1)=s(1,2); s(2,2)=S22(i,j,k); s(2,3)=S23(i,j,k)
            s(3,1)=s(1,3); s(3,2)=s(2,3); s(3,3)=-(s(1,1)+s(2,2))
            CALL GetRoots(s,eval)
            IF ( GetRank(s).LT.2 ) THEN
              evec = 0.0_GP
              CYCLE
            ENDIF
            m1 = s;
            m1(1,1) = s(1,1)-eval(1); m1(2,2) = s(2,2)-eval(1); m1(3,3) = s(3,3)-eval(1)
            IF ( GetRank(m1).LT.2 ) THEN
              evec = 0.0_GP
              CYCLE
            ENDIF
            row1(1:3) = m1(1,1:3)
            row2(1:3) = m1(2,1:3)
            CALL GetComplement1(row1,row2,ev1)

            m2 = s;
            m2(1,1) = s(1,1)-eval(2); m2(2,2) = s(2,2)-eval(2); m2(3,3) = s(3,3)-eval(2)
            IF ( GetRank(m2).LT.2 ) THEN
              evec(1:3) = 0.0_GP
              CYCLE
            ENDIF
            row1(1:3) = m2(1,1:3)
            row2(1:3) = m2(2,1:3)
            CALL GetComplement1(row1,row2,ev2)

            m3 = s;
            m3(1,1) = s(1,1)-eval(3); m3(2,2) = s(2,2)-eval(3); m3(3,3) = s(3,3)-eval(3)
            IF ( GetRank(m3).LT.2 ) THEN
              evec(1:3) = 0.0_GP
              CYCLE
            ENDIF
            CALL GetComplement1(ev1,ev2,evec)
            evx(i,j,k)=evec(1); evy(i,j,k)=evec(2); evz(i,j,k)=evec(3); lamb(i,j,k)=abs(eval(3))
          ENDDO
        ENDDO
      ENDDO


      END SUBROUTINE EigenSolveMax
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      RECURSIVE SUBROUTINE GetComplement1(u,v,w)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the normalized cross product of 3-vectors u,v, and
! stores result in w
!
! Parameters
!     u,v  : input 3-vectors
!     w    : returned result
!
      USE fprecision
      USE ali
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3) :: u,v
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3) :: w
      REAL   (KIND=GP)                              :: cmag
      INTEGER                                       :: j  

      ! compute cross product:
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)

      cmag = 0.0_GP
      DO j = 1,3
        cmag = cmag + w(j)*w(j)
      ENDDO
      cmag = sqrt(cmag) 
      IF ( cmag.GT.tiny ) cmag = 1.0_GP/cmag
      DO j = 1,3
        w(j) = w(j)*cmag
      ENDDO

      END SUBROUTINE GetComplement1
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      RECURSIVE INTEGER FUNCTION GetRank(m)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the rank of the 3x3 matrix, m, and returns it
!
! Parameters
!     m  : input 3x3 matrix
!
      USE fprecision
      USE ali
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3,3) :: m
      REAL   (KIND=GP)                                :: fabs,fmax,fmaxi,fsave
      INTEGER                                         :: i,j,maxrow,maxcol

      maxrow = 0
      maxcol = 0

      fmax = 0.0_GP
      DO i = 1, 3
        DO j = i, 3
          fabs = abs(m(i,j))
          IF ( fabs .GT. fmax ) THEN
            fmax = fabs
            maxrow = i
            maxcol = j
          ENDIF
        ENDDO
      ENDDO

      IF ( fmax .LT. tiny ) THEN
        GetRank = 0 ! e-value has multiplicity three
        RETURN
      ENDIF
     
      ! rank is >= 1. Swap row containing the max-magnitude entry
      ! with wow 1:

      IF ( maxrow .NE. 1 ) THEN
       DO j = 1, 3
          DO i = 1, 3
            fsave = m(1,j)
            m(1,j) = m(maxrow,j)
            m(maxrow,j) = fsave
          ENDDO
        ENDDO
      ENDIF

      ! row-reduce matrix:

      ! scale row 1 to generate a 1-valued pivot:
      fmaxi = 0.0_GP
      IF ( abs(m(1,maxcol)).GT.tiny ) fmaxi = 1.0_GP / m(1,maxcol)
      m(1,1) = m(1,1)*fmaxi
      m(1,2) = m(1,2)*fmaxi
      m(1,3) = m(1,3)*fmaxi
     
      ! eliminate maxcol column entries in rows 2 & 2:
      IF ( maxcol .EQ. 1 ) THEN 
        m(2,2) = m(2,2) - m(2,1)*m(1,2)
        m(2,3) = m(2,3) - m(2,1)*m(1,3)
        m(3,2) = m(3,2) - m(3,1)*m(1,2)
        m(3,3) = m(3,3) - m(3,1)*m(1,3)
        m(2,1) = 0.0_GP
        m(3,1) = 0.0_GP
      ELSE IF ( maxcol .EQ. 2 ) THEN
        m(2,1) = m(2,1) - m(2,2)*m(1,1)
        m(2,3) = m(2,3) - m(2,2)*m(1,3)
        m(3,1) = m(3,1) - m(3,2)*m(1,1)
        m(3,3) = m(3,3) - m(3,2)*m(1,3)
        m(2,2) = 0.0_GP
        m(3,2) = 0.0_GP
      ELSE
        m(2,1) = m(2,1) - m(2,3)*m(1,1)
        m(2,2) = m(2,2) - m(2,3)*m(1,2)
        m(3,1) = m(3,1) - m(3,3)*m(1,1)
        m(3,2) = m(3,2) - m(3,3)*m(1,2)
        m(2,3) = 0.0_GP
        m(3,3) = 0.0_GP
      ENDIF

      ! compute max-magnitude entry of the last 2 rows
      ! of row-reduced matrix:
      fmax = -1.0_GP
      maxrow = 0
      maxcol = 0

      DO j = 1, 3
        DO i = 2, 3
          fabs = abs(m(i,j))
          IF ( fabs .GT. fmax ) THEN
            fmax = fabs
            maxrow = i
            maxcol = j
          ENDIF
        ENDDO
      ENDDO

      IF ( fmax .lt. tiny ) THEN
        GetRank = 1 ! e-value has multiplicity 2
        RETURN
      ENDIF

      ! if row 2 has the max-magnitude entry, swap it with 
      ! row 1:
      IF ( maxrow .EQ. 3 ) THEN
        DO j = 1, 3
          fsave = m(2,j)
          m(2,j) = m(3,j)
          m(3,j) = fsave
        ENDDO
      ENDIF

      ! scale row 1 to generate a 1-vaued point:
      fmaxi = 0.0_GP
      IF ( abs(m(2,maxcol)).GT.tiny ) fmaxi = 1.0_GP / m(2,maxcol)
      m(2,1) = m(2,1)*fmaxi
      m(2,2) = m(2,2)*fmaxi
      m(2,3) = m(2,3)*fmaxi
 
      GetRank = 2 ! e-value has multiplicity 1

      END FUNCTION GetRank
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      RECURSIVE SUBROUTINE GetRoots(m,eval)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the ordered eigen values of input 3x3 symmetric matrix, m,
! returns in eval vector from smallest to largest. It is assumed that
! the input matrix, m, is trace free, which applies to the straing
! rate matrix for incompressible flows only!
!
! Parameters
!     m    : input 3x3 matrix
!     eval : eigenvalues from smallest to largest returned 
!
      USE fprecision
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(3,3) :: m
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION  (3) :: eval
      REAL   (KIND=GP)                                :: sa,sb,sc,sd,se,sf
      REAL   (KIND=GP)                                :: a,b,c,d
      REAL   (KIND=GP)                                :: del0,del1
      COMPLEX(KIND=GP)                                :: CC,u(3)
      COMPLEX(KIND=GP)                                :: D0,D1
      INTEGER                                         :: j,l

      u(1) = cmplx(1.0_GP , 0.0_GP)
      u(2) = cmplx(-0.5_GP, 0.5_GP*sqrt(3.0_GP))
      u(3) = cmplx(-0.5_GP,-0.5_GP*sqrt(3.0_GP))

      sa   = m(1,1); sb = m(1,2); sc = m(1,3); 
      sd   = m(2,2); se = m(2,3); sf = -(sa+sd)
      a    = 1.0_GP
      b    = 0.0_GP ! sa + sd + sf; trace-free
      c    = sa*sd+sa*sf+sd*sf-sb**2 - sc**2 - se**2
      d    = sa*se**2 + sf*sb**2 + sd*sc**2 - sa*sd*sf - 2.0_GP*sb*sc*se
!     del0 = b**2 - 3.0_GP*a*c
!     del1 = 2.0_GP*b**3 - 9.0_GP*a*b*c + 27.0_GP*d*a**2
      del0 = -3.0_GP*a*c
      del1 =  27.0_GP*d*a**2
      D0   = cmplx(del0,0)
      D1   = cmplx(del1,0)
      CC   = (0.5*D1 + 0.5*sqrt(D1**2 -4.0_GP*D0**3) )**(1.0_GP/3.0_GP)
      DO l = 1,3
!       eval(l) = real(-( b + u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
        eval(l) = real(-( u(l)*CC + D0/(u(l)*CC))/(3.0_GP*a),kind=GP)
      ENDDO

      ! reorder from smallest to largest:in magnitude
      DO j = 1,2
        DO l = j+1,3
          IF ( abs(eval(l)).LT.abs(eval(j)) ) THEN
            a       = eval(j)
            eval(j) = eval(l)
            eval(l) = a
          ENDIF
        ENDDO
      ENDDO

      END SUBROUTINE GetRoots
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE DoDissJPDF(S11,S12,S13,S22,S23,vx,vy,vz,lambda,ext, &
                            fnjde,fnjdl,fnjel,fnjdh,fnjeh, &
                            nbins,dolog,fndiss,fnenst,fnlamb,fnhel,fnflat,fnskew,ctmp,vtmp,rtmp)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes the joint PDF of energy dissipation/(2*nu) and ( enstrophy
! lambda, helicity), ! and output to files. 1d PDFs of dissipation, 
! enstrophy density, lambda, and helicity are also written.
!
! Note: after this call, the data should be expected to be overwritten.
!
! Parameters
!     SIJ   : strain rate components (real)
!     vx,
!     vy,
!     vz    : complex velocities
!     lambda: e-value field
!     ext   : string time index
!     fnjde : output file name for diss-enstr JPDF
!     fnjdl : output file name for diss-lambda JPDF
!     fnjel : output file name for diss-enst JPDF
!     fnjdh : output file name for diss-helicity JPDF
!     fnjeh : output file name for enstrophy-helicity JPDF
!     nbins : 2-elem array with no. bins for (enstrophy,energy diss)
!     dolog : flag to do (1) or not (0) logs of ranges when computing bins
!     fndiss: file name for 1d PDF of dissipation  
!     fnenst: file name for 1d PDF of enstrophy density
!     fnlamb: file name for 1d PDF of lambda
!     fnflat: file name for flatness quantities
!     fnskew: file name for skewness quantities
!     fnhel : file name for 1d PDF of helicity
!     ctmp  : complex temp array of size vx,vy,vz
!     vtmp  : complex temp array of size vx,vy,vz
!     rtmp  : real temp array of size lambda
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      USE fft
      USE var
      USE fftplans
      USE ali
      USE gutils
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: lambda,rtmp,S11,S12,S13,S22,S23
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: ctmp,vtmp
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      REAL   (KIND=GP)                                          :: fact,fmin(2),fmax(2),xnorm,xnormi
      REAL   (KIND=GP)                                          :: ss11,ss12,ss13,ss22,ss23,ss33
      REAL   (KIND=GP)                                          :: sf11,sf12,sf13,sf22,sf23,sf33
      REAL   (KIND=GP)                                          :: s2,s3,s4
      REAL   (KIND=GP)                                          :: ws1,ws2,ws3
      REAL   (KIND=GP)                                          :: wf1,wf2,wf3
      REAL   (KIND=GP)                                          :: fdiss,sdiss
      REAL   (KIND=GP)                                          :: fenst,senst
      INTEGER         , INTENT   (IN)                           :: dolog,nbins(2)
      INTEGER                                                   :: i,j,k,nin,sr
      CHARACTER(len=*), INTENT   (IN)                           :: fnjde,fnjdl,fnjel,fnjdh,fnjeh,fndiss
      CHARACTER(len=*), INTENT   (IN)                           :: fnenst,fnlamb,fnhel,fnflat,fnskew,ext


      nin = n*n*(kend-ksta+1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            vx(k,j,i) = vx(k,j,i)/real(n,kind=GP)**3
            vy(k,j,i) = vy(k,j,i)/real(n,kind=GP)**3
            vz(k,j,i) = vz(k,j,i)/real(n,kind=GP)**3
          END DO
        END DO
      END DO

      ! compute S33 for flatness/skewness:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            rtmp(i,j,k) =  -(S11(i,j,k)+S22(i,j,k)) 
          ENDDO
        ENDDO
      ENDDO

      CALL skewflat(S11 ,nin,n,ss11,sf11,s2,s3,s4)
if ( myrank.eq.0 ) write(*,*)' time=',ext,' s11_s2=',s2,' s11_s3=',s3,' s11_s4=',s4
      CALL skewflat(S12 ,nin,n,ss12,sf12,s2,s3,s4)
if ( myrank.eq.0 ) write(*,*)' time=',ext,' s12_s2=',s2,' s12_s3=',s3,' s12_s4=',s4
      CALL skewflat(S13 ,nin,n,ss13,sf13,s2,s3,s4)
if ( myrank.eq.0 ) write(*,*)' time=',ext,' s13_s2=',s2,' s13_s3=',s3,' s12_s4=',s4
      CALL skewflat(S22 ,nin,n,ss22,sf22,s2,s3,s4)
if ( myrank.eq.0 ) write(*,*)' time=',ext,' s22_s2=',s2,' s22_s3=',s3,' s22_s4=',s4
      CALL skewflat(S23 ,nin,n,ss23,sf23,s2,s3,s4)
if ( myrank.eq.0 ) write(*,*)' time=',ext,' s23_s2=',s2,' s23_s3=',s3,' s23_s4=',s4
      CALL skewflat(rtmp,nin,n,ss33,sf33,s2,s3,s4)
if ( myrank.eq.0 ) write(*,*)' time=',ext,' s33_s2=',s2,' s33_s3=',s3,' s33_s4=',s4
      ! Compute, write, 1d Sij pdfs:
      CALL dopdfr(S11 ,nin,n,'s11pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 
      CALL dopdfr(S12 ,nin,n,'s12pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 
      CALL dopdfr(S13 ,nin,n,'s13pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 
      CALL dopdfr(S22 ,nin,n,'s22pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 
      CALL dopdfr(S23 ,nin,n,'s23pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 
      CALL dopdfr(rtmp,nin,n,'s33pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 

      ! Compute normalized energy dissipation field, store in ! S11:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            S11(i,j,k) = 2.0_GP*( S11(i,j,k)*S11(i,j,k) &
                                 +S22(i,j,k)*S22(i,j,k) &
                                 +S22(i,j,k)*S11(i,j,k) &
                                 +S12(i,j,k)*S12(i,j,k) &
                                 +S13(i,j,k)*S13(i,j,k) &
                                 +S23(i,j,k)*S23(i,j,k) &
                                )
          ENDDO
        ENDDO
      ENDDO

      ! Compute joint PDF for energy diss and lambda (order 
      ! switched, so that lambda is on x-axis, and energy diss on y axis:
      CALL dojpdfr(lambda,'lambda',S11,'diss',nin,n,fnjdl,nbins,0,fmin,fmax,[dolog,dolog])

      ! Do 1d and jooint pdfs for diss and lambda:
      CALL dopdfr(S11   ,nin,n,fndiss,nbins(2),0,fmin(2),fmax(2),dolog) 
      ! Compute, write, 1d lambda pdf:
      CALL dopdfr(lambda,nin,n,fnlamb,nbins(1),0,fmin(1),fmax(1),dolog) 

      ! Compute enstrophy density,helicity field, store in S22, S13:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            S22 (i,j,k) = 0.0_GP ! Enstrophy density: omega^2
            S12 (i,j,k) = 0.0_GP ! v^2
            S13 (i,j,k) = 0.0_GP ! v.omega
          ENDDO
        ENDDO
      ENDDO

!     fact  = 1.0_GP/real(n,kind=GP)**6 ! no  longer needed; is in v_i now
      CALL rotor3(vy,vz,ctmp,1)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23   ,MPI_COMM_WORLD)
      vtmp = vx
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            S22 (i,j,k) = S22 (i,j,k) + S23   (i,j,k)**2         ! Enstrophy dens
            S13 (i,j,k) = S13 (i,j,k) + rtmp  (i,j,k)*S23(i,j,k) ! v.omega
            S12 (i,j,k) = S12 (i,j,k) + rtmp  (i,j,k)**2         ! v^2
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(rtmp,nin,n,ws1,wf1,s2,s3,s4)
      ! Compute, write, 1d \omega_i pdfs:
      CALL dopdfr(rtmp,nin,n,'w1pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 

      CALL rotor3(vz,vx,ctmp,2)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23   ,MPI_COMM_WORLD)
      vtmp = vy
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp  ,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            S22 (i,j,k) = S22 (i,j,k) + S23   (i,j,k)**2         ! Enstrophy dens
            S13 (i,j,k) = S13 (i,j,k) + rtmp  (i,j,k)*S23(i,j,k) ! v.omega
            S12 (i,j,k) = S12 (i,j,k) + rtmp  (i,j,k)**2         ! v^2
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(rtmp,nin,n,ws2,wf2,s2,s3,s4)
      ! Compute, write, 1d \omega_i pdfs:
      CALL dopdfr(rtmp,nin,n,'w2pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 

      CALL rotor3(vx,vy,ctmp,3)
      CALL fftp3d_complex_to_real(plancr,ctmp,S23,MPI_COMM_WORLD)
      vtmp = vz
      CALL fftp3d_complex_to_real(plancr,vtmp,rtmp,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            S22 (i,j,k) = S22 (i,j,k) + S23   (i,j,k)**2         ! Enstrophy dens
            S13 (i,j,k) = S13 (i,j,k) + rtmp  (i,j,k)*S23(i,j,k) ! v.omega
            S12 (i,j,k) = S12 (i,j,k) + rtmp  (i,j,k)**2         ! v^2
          ENDDO
        ENDDO
      ENDDO
      CALL skewflat(rtmp,nin,n,ws3,wf3,s2,s3,s4)
      CALL skewflat(S11 ,nin,n,sdiss,fdiss,s2,s3,s4)
      CALL skewflat(S22 ,nin,n,senst,fenst,s2,s3,s4)
      ! Compute, write, 1d \omega_i pdfs:
      CALL dopdfr(rtmp,nin,n,'w3pdf.'//ext//'.txt',nbins(1),0,fmin(1),fmax(1),0) 

      ! Print out skewness and flatness data:
      IF ( myrank.EQ.0 ) THEN
        OPEN(1,file=trim(fnskew),position='append')
        WRITE(1,*)ext,ss11,ss12,ss13,ss22,ss23,ss33,ws1,ws2,ws3,sdiss,senst
        CLOSE(1)
        OPEN(1,file=trim(fnflat),position='append')
        WRITE(1,*)ext,sf11,sf12,sf13,sf22,sf23,sf33,wf1,wf2,wf3,fdiss,fenst
        CLOSE(1)
      ENDIF

!     Compute relative helicity:
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
        DO j = 1,n
          DO i = 1,n
            sr = sign(1.0,S13(i,j,k));
            xnorm = sqrt( S22(i,j,k)*S12(i,j,k) )
            xnormi = 0.0_GP
            IF ( xnorm .GT. tiny ) THEN
              xnormi = 1.0_GP/xnorm
            ENDIF
            S13 (i,j,k) = S13(i,j,k) * xnormi
            S13 (i,j,k) = sr*min(abs(S13(i,j,k)),1.0)
          ENDDO
        ENDDO
      ENDDO

      ! Compute, write, 1d enstrophy density pdf:
      CALL dopdfr(S22,nin,n,fnenst,nbins(1),0,fmin(1),fmax(1),dolog) 
      ! Compute, write, 1d relative helicity pdf:
      CALL dopdfr(S13,nin,n,fnhel ,nbins(1),0,fmin(1),fmax(1),0) 

      ! Compute joint PDF for enstrophy and diss (order switched, 
      ! so that enstrophy is on x-axis, and energy diss on y axis:
      CALL dojpdfr(S22,'enst',S11,'diss',nin,n,fnjde,nbins,0,fmin,fmax,[dolog,dolog])
      ! Compute joint PDF for rel. helicity and diss :
      CALL dojpdfr(S13,'rhel',S11,'diss',nin,n,fnjdh,nbins,0,fmin,fmax,[0,dolog])
      ! Compute joint PDF for rel. helicity and enstroph. density:
      CALL dojpdfr(S13,'rhel',S22,'enst',nin,n,fnjeh,nbins,0,fmin,fmax,[0,dolog])


      CALL dojpdfr(lambda,'lambda',S22,'enst',nin,n,fnjdl,nbins,0,fmin,fmax,[dolog,dolog])

      END SUBROUTINE DoDissJPDF
!
!
      SUBROUTINE skewflat(fx,nin,n,skew,flat,s2,s3,s4)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes skewnewss and flatness of input random variable
! 
!
! Parameters
!     fx    : input real random variable
!     nin   : total size of input array
!     n     : linear size of input array (s.t. n^3 is global size)
!     skew  : skewness, valid only for MPI taks 0
!     flat  : flatness/kurtosis
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      IMPLICIT NONE

      REAL(KIND=GP), INTENT (IN), DIMENSION(nin) :: fx
      REAL(KIND=GP), INTENT(OUT)                 :: skew,flat
      REAL(KIND=GP), INTENT(OUT)                 :: s2,s3,s4
      REAL(KIND=GP)                              :: avg,gs(3),s(3),xnorm
      INTEGER      , INTENT (IN)                 :: nin,n
      INTEGER                                    :: ierr,j

      xnorm = 1.0_GP/real(n,kind=GP)**3

      s2 = 0.0_GP
!$omp parallel do private(j) reduction(+:s2)
      DO j = 1, nin
        s2 = s2 + fx(j)
      ENDDO

      CALL MPI_ALLREDUCE(s2, avg, 1, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      avg = avg*xnorm

      s2 = 0.0_GP
!$omp parallel do private(j) reduction(+:s2)
      DO j = 1, nin
        s2 = s2 + (fx(j)-avg)**2
      ENDDO

      s3 = 0.0_GP
!$omp parallel do private(j) reduction(+:s3)
      DO j = 1, nin
        s3 = s3 + (fx(j)-avg)**3
      ENDDO

      s4 = 0.0_GP
!$omp parallel do private(j) reduction(+:s4)
      DO j = 1, nin
        s4 = s4 + (fx(j)-avg)**4
      ENDDO

      s(1)=s2; s(2)=s3; s(3)=s4
      CALL MPI_ALLREDUCE(s, gs, 3, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      if ( ierr.ne.MPI_SUCCESS ) then
        write(*,*)'skewflat: final allreduce failed'
        stop
      endif
      s2=gs(1)*xnorm; s3=gs(2)*xnorm; s4=gs(3)*xnorm

      skew = s3 / s2**(1.5)
      flat = s4 / s2**(2.0)

      END SUBROUTINE SkewFlat
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!
!
      SUBROUTINE Randomize(fx,fy,fz,krmin,krmax,phase)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Randomizes phases of vector field, (vx,vy,vz) between
! wavenumbers krmin,krmax
!
! Parameters
!     vx,
!     vy,
!     vz    : complex velocities
!     krmin : minimum wavenumber
!     krmax : maximum wavenumber
!     phase : phase 
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE threads
      USE fft
      USE var
      USE fftplans
      USE gutils
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: fx,fy,fz
      REAL   (KIND=GP), INTENT   (IN)                           :: krmin,krmax,phase
      REAL   (KIND=GP)                                          :: krmin2,krmax2
      COMPLEX(KIND=GP)                                          :: cdump,jdump
      INTEGER                                                   :: i,j,k

      cdump = COS(phase)+im*SIN(phase)
      jdump = conjg(cdump)

      krmin2 = krmin**2
      krmax2 = krmax**2
      IF (ista.eq.1) THEN

      DO j = 2,n/2+1
        IF ( ka2(1,j,1).gt.krmin2 .and. ka2(1,j,1).lt.krmax2 ) THEN
        fx    (1,j,1) = fx    (1,j,1)*cdump
        fx(1,n-j+2,1) = fx(1,n-j+2,1)*jdump
        fy    (1,j,1) = fy    (1,j,1)*cdump
        fy(1,n-j+2,1) = fy(1,n-j+2,1)*jdump
        fz    (1,j,1) = fz    (1,j,1)*cdump
        fz(1,n-j+2,1) = fz(1,n-j+2,1)*jdump
        ENDIF
      ENDDO

      !$omp parallel do
      DO k = 2,n/2+1
        IF ( ka2(k,1,1).gt.krmin2 .and. ka2(k,1,1).lt.krmax2 ) THEN
        fx    (k,1,1) = fx    (k,1,1)*cdump
        fx(n-k+2,1,1) = fx(n-k+2,1,1)*jdump
        fy    (k,1,1) = fy    (k,1,1)*cdump
        fy(n-k+2,1,1) = fy(n-k+2,1,1)*jdump
        fz    (k,1,1) = fz    (k,1,1)*cdump
        fz(n-k+2,1,1) = fz(n-k+2,1,1)*jdump
        ENDIF
      END DO
      
!$omp parallel do private (k)
      DO j = 2,n
        DO k = 2,n/2+1
          IF ( ka2(k,j,1).gt.krmin2 .and. ka2(k,j,1).lt.krmax2 ) THEN
          fx        (k,j,1) = fx        (k,j,1)*cdump
          fx(n-k+2,n-j+2,1) = fx(n-k+2,n-j+2,1)*jdump
          fy        (k,j,1) = fy        (k,j,1)*cdump
          fy(n-k+2,n-j+2,1) = fy(n-k+2,n-j+2,1)*jdump
          fz        (k,j,1) = fz        (k,j,1)*cdump
          fz(n-k+2,n-j+2,1) = fz(n-k+2,n-j+2,1)*jdump
          ENDIF
        ENDDO
      ENDDO


!$omp parallel do if (iend-2.ge.nth) private (j,k)
      DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k)
        DO j = 1,n
          DO k = 1,n
            IF ( ka2(k,j,i).gt.krmin2 .and. ka2(k,j,i).lt.krmax2 ) THEN
            fx(k,j,i) = fx(k,j,i)*cdump
            fy(k,j,i) = fy(k,j,i)*cdump
            fz(k,j,i) = fz(k,j,i)*cdump
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ELSE
    
      DO  i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
           DO k = 1,n
             IF ( ka2(k,j,i).gt.krmin2 .and. ka2(k,j,i).lt.krmax2 ) THEN
             fx(k,j,i) = fx(k,j,i)*cdump
             fy(k,j,i) = fy(k,j,i)*cdump
             fz(k,j,i) = fz(k,j,i)*cdump
            ENDIF
           ENDDO
         ENDDO
       ENDDO
       ENDIF
  
      END SUBROUTINE Randomize
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
