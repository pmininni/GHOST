!=================================================================
      PROGRAM WV3D
!=================================================================
! WV3D code (part of the GHOST suite)
!
! Reads velocity and scalar binaries and computes wave-vortical
! decomposition of the flow. Utility can average all modes in
! Fourier space first before decomposing, or it can operate  
! separately on file corresponding to a list of time-stamps.
!
! Decomposition utilizes that specified 
!  C. Herbert et al. JFM 758:374 (2014) which is closely 
! connected to that in Bartello: J. Atm. Sci.52(24):4410 (1995), 
! except that the propagation matrix is Hermitian in the former.
! This allows us to store only a single set of wave modes.
!
! 2014 D. Rosenberg
!      ORNL
!
! 14 Nov 2014: Initial version
! 11 Feb 2015: Corrections for potential enstrophy computations (CH)
!              Modified PVn, wgradt and wvzspectrum
!              Added functions specprodg and specprodaxig for cubic term
! 12 Mar 2015: Added projections wvproj0, wvprojw (CH)
!              Added transfer functions wvtrans (CH)
!              and tmp aux functions: entransc, entrans2Dc, sctransc sctrans2Dc
!              Ultimately put in pseudo?
! 02 May 2015: Performance improvements for transfer functions (CH)
! 06 May 2015: Added wvtransfull to compute detailed transfer functions (CH)
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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz,th,a0,am,ap
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: c1,c2,c3,c4
      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: r1,r2,r3
      REAL   (KIND=GP)                                 :: omega,bvfreq
      REAL   (KIND=GP)                                 :: tmp
!
! Auxiliary variables

      INTEGER :: i,i2d,it,iavg,ic,ind,ivec,putnm,j,k,trans
      INTEGER :: istat(1024),nfiles,nstat

      TYPE(IOPLAN)        :: planio,planioc
      CHARACTER(len=8)    :: pref
      CHARACTER(len=2048) :: odir,idir
      CHARACTER(len=2048) :: fout
      CHARACTER(len=4096) :: stat
      CHARACTER(len=4)    :: ext1
!
      NAMELIST / wv / idir, odir, stat, iswap, oswap, iavg, putnm, ivec, omega, bvfreq, i2d, trans

      tiny  = 1e-5_GP
      tinyf = 1e-15_GP

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init (myrank,n,ksta,kend,planio)
      CALL io_initc(myrank,n,ista,iend,planioc)
      idir   = '.'
      odir   = '.'
      stat   = '0'
      iswap  = 0
      oswap  = 0
      iavg   = 0 
      i2d    = 0 
      putnm  = 1 
      ivec   = 0
      trans  = 0
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     stat  : time index for which to compute WV, or a ';--separated list
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     iavg   : time average spectrally the time index list?
!     i2d    : write 2d spectra? (1,0; 0 is default)
!     ivec   : 0==don't print vectors, ==1; print velocity comp only;
!              2==print magnitude only; ==3, print velocity comps and mag

!              default is 0 (magnitudes only)
!     trans  : write transfer functions? (1,0; 0 is default)
!     putnm  : write normal mode fields (0==don't; 1==in wave space; 2==in real space, 3==both)
!     omega  : rotation rate
!     bvfreq : Brunt-Vaisalla frequency
      IF (myrank.eq.0) THEN
         OPEN(1,file='wv.txt',status='unknown',form="formatted")
         READ(1,NML=wv)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,2048,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,2048,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat  ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iavg  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(i2d   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ivec  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(putnm ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(trans ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omega ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bvfreq,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)

      pref = 'WV'
!
     
      kmax   = (real(n,kind=GP)/3.)**2

      ALLOCATE( a0(n,n,ista:iend) )
      ALLOCATE( am(n,n,ista:iend) )
      ALLOCATE( ap(n,n,ista:iend) )
      ALLOCATE( vx(n,n,ista:iend) )
      ALLOCATE( vy(n,n,ista:iend) )
      ALLOCATE( vz(n,n,ista:iend) )
      ALLOCATE( th(n,n,ista:iend) )
      ALLOCATE( c1(n,n,ista:iend) )
      ALLOCATE( c2(n,n,ista:iend) )
      ALLOCATE( c3(n,n,ista:iend) )
      ALLOCATE( c4(n,n,ista:iend) )
      ALLOCATE( r1(n,n,ksta:kend) )
      ALLOCATE( r2(n,n,ksta:kend) )
      ALLOCATE( r3(n,n,ksta:kend) )
      ALLOCATE( ka(n),ka2(n,n,ista:iend) )
!

      CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
          FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
          FFTW_MEASURE)
!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing

      kmax = (REAL(n,KIND=GP)/3.)**2
      tiny = 1e-5

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
      ext1 = ''
      tmp = 1.0/REAL(n,KIND=GP)**3
      IF ( iavg.EQ.1 ) THEN
         CALL DoSpAvg(vx,vy,vz,th,istat,nstat,idir,planrc,planio,r1,c1)
         CALL WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq)
         WRITE(ext , fmtext) istat(1)     
         WRITE(ext1, fmtext) istat(nstat) 
         IF ( ibits(putnm,0,1).EQ.1 ) THEN ! output cmplx norm. modes
            c1 = a0
            WRITE(fout,'(a,i4.4,a,i4.4,a)'),'a0av_',istat(1),'_',istat(nstat),'.out'
            CALL PutNormModes(1,odir,trim(fout),planioc,c1)

            c1 = am
            WRITE(fout,'(a,i4.4,a,i4.4,a)'),'amav_',istat(1),'_',istat(nstat),'.out'
            CALL PutNormModes(1,odir,trim(fout),planioc,c1)
            
            c1 = ap
            WRITE(fout,'(a,i4.4,a,i4.4,a)'),'apav_',istat(1),'_',istat(nstat),'.out'
            CALL PutNormModes(1,odir,trim(fout),planioc,c1)
         ENDIF
          IF ( ibits(putnm,1,1).EQ.1 ) THEN ! output real fields
            CALL PutRealFields(1,odir,ext,planio,a0,am,ap,vx,vy,vz,th,omega,bvfreq,c1,c2,c3,c4,r1,r2,ivec,'av')
          ENDIF
         IF ( trans.LE.10 ) THEN
!            CALL WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq) ! this seems unnecessary
            CALL wvspectrum(a0,am,ap,omega,bvfreq,odir,ext,ext1,i2d)
            CALL wvzspectrum(a0,vx,vy,vz,th,omega,bvfreq,odir,ext,ext1,i2d,c1,c2,c3,r1,r2,r3)
         ENDIF
         IF ( MOD(trans,10).EQ.1 ) THEN
            CALL wvtrans(a0,am,ap,vx,vy,vz,th,omega,bvfreq,odir,ext,ext1,i2d,c1,c2,c3)
         ELSEIF ( MOD(trans,10).EQ.2 ) THEN
            CALL wvtransfull(a0,am,ap,vx,vy,vz,th,omega,bvfreq,odir,ext,ext1,i2d,c1,c2,c3)
         ENDIF
      ELSE
         DO it = 1,nstat
          WRITE(ext, fmtext) istat(it)
! read in appropriate file:
!if ( myrank.eq.0 ) write(*,*)'main: reading vx...',ext
          CALL io_read(1,idir,'vx',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
!if ( myrank.eq.0 ) write(*,*)'main: fft_rc on vx...',ext
          CALL fftp3d_real_to_complex(planrc,r1,vx,MPI_COMM_WORLD)
!if ( myrank.eq.0 ) write(*,*)'main: reading vy...',ext
          CALL io_read(1,idir,'vy',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
!if ( myrank.eq.0 ) write(*,*)'main: fft_rc on vy...', ext
          CALL fftp3d_real_to_complex(planrc,r1,vy,MPI_COMM_WORLD)
!if ( myrank.eq.0 ) write(*,*)'main: reading vz...',ext
          CALL io_read(1,idir,'vz',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
!if ( myrank.eq.0 ) write(*,*)'main: fft_rc on vz...',ext
          CALL fftp3d_real_to_complex(planrc,r1,vz,MPI_COMM_WORLD)
!if ( myrank.eq.0 ) write(*,*)'main: reading th...',ext
          CALL io_read(1,idir,'th',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
!if ( myrank.eq.0 ) write(*,*)'main: fft_rc on th...',ext
          CALL fftp3d_real_to_complex(planrc,r1,th,MPI_COMM_WORLD)
!if ( myrank.eq.0 ) write(*,*)'main: calling WVNormal ...'
          CALL WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq)
          IF ( trans.LE.10 ) THEN
!if ( myrank.eq.0 ) write(*,*)'main: calling wvspectrum...'
             CALL wvspectrum(a0,am,ap,omega,bvfreq,odir,ext,'',i2d)
!if ( myrank.eq.0 ) write(*,*)'main: calling wvzspectrum...'
             CALL wvzspectrum(a0,vx,vy,vz,th,omega,bvfreq,odir,ext,'',i2d,c1,c2,c3,r1,r2,r3)
!if ( myrank.eq.0 ) write(*,*)'main: wvzspectrum done'
          ENDIF
          IF ( MOD(trans,10).EQ.1 ) THEN
!if ( myrank.eq.0 ) write(*,*)'main: calling wvtrans...'
             CALL wvtrans(a0,am,ap,vx,vy,vz,th,omega,bvfreq,odir,ext,'',i2d,c1,c2,c3)
          ELSEIF ( MOD(trans,10).EQ.2 ) THEN
!if ( myrank.eq.0 ) write(*,*)'main: calling wvtransful...'
             CALL wvtransfull(a0,am,ap,vx,vy,vz,th,omega,bvfreq,odir,ext,'',i2d,c1,c2,c3)             
          ENDIF

!if ( myrank.eq.0 ) write(*,*)'main: doing binary output...'

          ! Do output:
          IF ( ibits(putnm,0,1).EQ.1 ) THEN ! output norm. mode coeffs
!if ( myrank.eq.0 ) write(*,*)'main: doing PutNormModes...'
            c1 = a0
            fout = 'a0.' // ext // '.out'
            CALL PutNormModes(1,odir,trim(fout),planioc,c1)

            c1 = am
            fout = 'am.' // ext // '.out'
            CALL PutNormModes(1,odir,trim(fout),planioc,c1)

            c1 = ap
            fout = 'ap.' // ext // '.out'
            CALL PutNormModes(1,odir,trim(fout),planioc,c1)
          ENDIF
          IF ( ibits(putnm,1,1).EQ.1 ) THEN ! output real fields
!if ( myrank.eq.0 ) write(*,*)'main: doing PutRealFields...'
            CALL PutRealFields(1,odir,ext,planio,a0,am,ap,vx,vy,vz,th,omega,bvfreq,c1,c2,c3,c4,r1,r2,ivec,'')
          ENDIF
        ENDDO ! end of it loop
      ENDIF
!
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)


      DEALLOCATE (a0,am,ap)
      DEALLOCATE (vx,vy,vz,th)
      DEALLOCATE (c1,c2,c3,c4)
      DEALLOCATE (r1,r2,r3)
      DEALLOCATE (ka)
      DEALLOCATE (ka2)

      CALL MPI_FINALIZE(ierr)

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END PROGRAM WV3D


!*****************************************************************
      SUBROUTINE WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Carries out wave-vortical decomposition tasks.
! Computes from Herbert, et al. JFM 758:374 (2014)
! eq A14,A16.
!
! Note: after this call, the data should be expected to be overwritten.
!
! Parameters
!     a0,p,m: complex temp array of size vx,vy,vz, containing normal modes, returned
!     vx,
!     vy,
!     vz    : complex velocities, overwritten with normal mode fields
!     th    : complex potential temperature
!     omega : rotation rate
!     bvfreq: Brunt-Vaisalla frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: a0,am,ap
      COMPLEX(KIND=GP)                                          :: ic
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig,tmp
      INTEGER                                                   :: i,j,k

      a0 = 0.0_GP
      am = 0.0_GP
      ap = 0.0_GP

      ic = cmplx(0.0_GP,1.0_GP);
      f = 2.0_GP * omega
!if ( myrank.eq.0 ) write(*,*)'bvfreq=',bvfreq,' omega=',omega,' f=',f,' tiny=',tiny
      tmp = 1.0_GP/sqrt(2.0_GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
     
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A16:
                    a0(k,j,i) = -th(k,j,i)
                    am(k,j,i) =  tmp*( vy(k,j,i) - ic*vx(k,j,i) )
                    ap(k,j,i) =  tmp*( vy(k,j,i) + ic*vx(k,j,i) )
                  ELSE
                    ! From decomp from Herbert, eq. A14:
                    ! A^0
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    a0(k,j,i) = ( -bvfreq*( ka(j)*vx(k,j,i) - ka(i)*vy(k,j,i) ) &
                              - f*ka(k)*th(k,j,i) ) / (ks*sig)
                    ! A^-
                    am(k,j,i) = ( f*ka(k)*( ka(j)*vx(k,j,i) - ka(i)*vy(k,j,i) ) &
                              - ic*ka2(k,j,i)*sig*vz(k,j,i) - bvfreq*kp**2*th(k,j,i) ) / (sqrt(2.0)*ks*kp*sig)
                    ! A^+
                    ap(k,j,i) = ( f*ka(k)*( ka(j)*vx(k,j,i) - ka(i)*vy(k,j,i) ) &
                              + ic*ka2(k,j,i)*sig*vz(k,j,i) - bvfreq*kp**2*th(k,j,i) ) / (sqrt(2.0)*ks*kp*sig)
                  
                 ENDIF

               END DO
            END DO
         END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE WVNormal


!*****************************************************************
      SUBROUTINE wvproj0(a0,vx,vy,vz,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the slow mode coefficients from the normal mode decomposition,
! computes the slow component of the fields
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
! Parameters
!     a0    : complex array of size vx,vy,vz, containing normal modes
!     vx,
!     vy    : complex slow velocities, returned
!             (no need to compute vz, it always vanishes for the slow modes)
!     th    : complex slow potential temperature, returned
!     omega : rotation rate
!     bvfreq: Brunt-Vaisala frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig
      INTEGER                                                   :: i,j,k

      vx = 0.0_GP
      vy = 0.0_GP
      th = 0.0_GP
      vz = 0.0_GP

      f = 2.0_GP * omega
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
 
                  vz(k,j,i) = 0.0_GP
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A15:
                    vx(k,j,i) = 0.0_GP
                    vy(k,j,i) = 0.0_GP                    
                    th(k,j,i) = -a0(k,j,i)
                  ELSE
                    ! From decomp from Herbert, eq. A11:
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    vx(k,j,i) = -bvfreq*ka(j)*a0(k,j,i)/(ks*sig)
                    vy(k,j,i) =  bvfreq*ka(i)*a0(k,j,i)/(ks*sig)
                    th(k,j,i) = -f*ka(k)*a0(k,j,i)/(ks*sig)
                 ENDIF

               END DO
            END DO
         END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvproj0

!*****************************************************************
      SUBROUTINE wvprojv0(a0,vx,vy,vz,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the slow mode coefficients from the normal mode decomposition,
! computes the slow component of the velocity field.
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
! Parameters
!     a0    : complex array of size vx,vy,vz, containing normal modes
!     vx,
!     vy    : complex slow velocities, returned
!             (no need to compute vz, it always vanishes for the slow modes)
!     omega : rotation rate
!     bvfreq: Brunt-Vaisala frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig
      INTEGER                                                   :: i,j,k

      vx = 0.0_GP
      vy = 0.0_GP
      vz = 0.0_GP

      f = 2.0_GP * omega
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
 
                  vz(k,j,i) = 0.0_GP
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A15:
                    vx(k,j,i) = 0.0_GP
                    vy(k,j,i) = 0.0_GP                    
                  ELSE
                    ! From decomp from Herbert, eq. A11:
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    vx(k,j,i) = -bvfreq*ka(j)*a0(k,j,i)/(ks*sig)
                    vy(k,j,i) =  bvfreq*ka(i)*a0(k,j,i)/(ks*sig)
                 ENDIF

               END DO
            END DO
         END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvprojv0

!*****************************************************************
      SUBROUTINE wvprojt0(a0,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the slow mode coefficients from the normal mode decomposition,
! computes the slow component of the temperature field.
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
! Parameters
!     a0    : complex array of size th, containing normal modes
!     th    : complex slow temperature, returned
!     omega : rotation rate
!     bvfreq: Brunt-Vaisala frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: th
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig
      INTEGER                                                   :: i,j,k

      th = 0.0_GP
      f = 2.0_GP * omega
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A15:
                    th(k,j,i) = -a0(k,j,i)
                  ELSE
                    ! From decomp from Herbert, eq. A11:
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    th(k,j,i) = -f*ka(k)*a0(k,j,i)/(ks*sig)
                 ENDIF

               END DO
            END DO
         END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvprojt0
      

!*****************************************************************
      SUBROUTINE wvprojw(am,ap,vx,vy,vz,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the wave mode coefficients from the normal mode decomposition,
! computes the projection of the fields on the wave manifold.
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
! Parameters
!     am,ap : complex array of size vx,vy,vz, containing normal modes
!     vx,
!     vy,
!     vz    : complex velocities, projected on the wave manifold, returned
!     th    : complex potential temperature projected on the wave manifold, returned
!     omega : rotation rate
!     bvfreq: Brunt-Vaisala frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: am,ap
      COMPLEX(KIND=GP)                                          :: ic
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig,tmp
      INTEGER                                                   :: i,j,k

      vx = 0.0_GP
      vy = 0.0_GP
      vz = 0.0_GP
      th = 0.0_GP

      ic = cmplx(0.0_GP,1.0_GP);
      f = 2.0_GP * omega
      tmp = 1.0_GP/sqrt(2.0_GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
     
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A15:
                    vx(k,j,i) = ic*tmp*(am(k,j,i)-ap(k,j,i))
                    vy(k,j,i) = tmp*(am(k,j,i)+ap(k,j,i))  
                    vz(k,j,i) = 0.0_GP
                    th(k,j,i) = 0.0_GP
                  ELSE
                    ! From decomp from Herbert, eq. A11:
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    vx(k,j,i) = tmp*ka(k)*(f*ka(j)*(ap(k,j,i)+am(k,j,i))+ic*ka(i)*sig*(ap(k,j,i)-am(k,j,i)))/(ks*kp*sig)
                    vy(k,j,i) = tmp*ka(k)*(-f*ka(i)*(ap(k,j,i)+am(k,j,i))+ic*ka(j)*sig*(ap(k,j,i)-am(k,j,i)))/(ks*kp*sig)
                    vz(k,j,i) = -tmp*ic*kp*(ap(k,j,i)-am(k,j,i))/ks
                    th(k,j,i) = -tmp*kp*bvfreq*(ap(k,j,i)+am(k,j,i))/(ks*sig)
                 ENDIF

               END DO
            END DO
         END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvprojw

!*****************************************************************
      SUBROUTINE wvprojvw(am,ap,vx,vy,vz,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the wave mode coefficients from the normal mode decomposition,
! computes the projection of the velocity field on the wave manifold.
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
! Parameters
!     am,ap : complex array of size vx,vy,vz, containing normal modes
!     vx,
!     vy,
!     vz    : complex velocities, projected on the wave manifold, returned
!     omega : rotation rate
!     bvfreq: Brunt-Vaisala frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: am,ap
      COMPLEX(KIND=GP)                                          :: ic
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig,tmp
      INTEGER                                                   :: i,j,k

      vx = 0.0_GP
      vy = 0.0_GP
      vz = 0.0_GP

      ic = cmplx(0.0_GP,1.0_GP);
      f = 2.0_GP * omega
      tmp = 1.0_GP/sqrt(2.0_GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
     
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A15:
                    vx(k,j,i) = ic*tmp*(am(k,j,i)-ap(k,j,i))
                    vy(k,j,i) = tmp*(am(k,j,i)+ap(k,j,i))  
                    vz(k,j,i) = 0.0_GP
                  ELSE
                    ! From decomp from Herbert, eq. A11:
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    vx(k,j,i) = tmp*ka(k)*(f*ka(j)*(ap(k,j,i)+am(k,j,i))+ic*ka(i)*sig*(ap(k,j,i)-am(k,j,i)))/(ks*kp*sig)
                    vy(k,j,i) = tmp*ka(k)*(-f*ka(i)*(ap(k,j,i)+am(k,j,i))+ic*ka(j)*sig*(ap(k,j,i)-am(k,j,i)))/(ks*kp*sig)
                    vz(k,j,i) = -tmp*ic*kp*(ap(k,j,i)-am(k,j,i))/ks
                 ENDIF

               END DO
            END DO
         END DO

      END SUBROUTINE wvprojvw

!*****************************************************************
      SUBROUTINE wvprojtw(am,ap,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the wave mode coefficients from the normal mode decomposition,
! computes the projection of the temperature field on the wave manifold.
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
! Parameters
!     am,ap : complex array of size vx,vy,vz, containing normal modes
!     th    : complex potential temperature projected on the wave manifold, returned
!     omega : rotation rate
!     bvfreq: Brunt-Vaisala frequency
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: th
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: am,ap
      COMPLEX(KIND=GP)                                          :: ic
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig,tmp
      INTEGER                                                   :: i,j,k

      th = 0.0_GP
      ic = cmplx(0.0_GP,1.0_GP);
      f = 2.0_GP * omega
      tmp = 1.0_GP/sqrt(2.0_GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
     
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A15:
                    th(k,j,i) = 0.0_GP
                  ELSE
                    ! From decomp from Herbert, eq. A11:
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    th(k,j,i) = -tmp*kp*bvfreq*(ap(k,j,i)+am(k,j,i))/(ks*sig)
                 ENDIF

               END DO
            END DO
         END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvprojtw
      


!*****************************************************************
      SUBROUTINE wvspectrum(a0,am,ap,omega,bvfreq,dir,nmb,nmb1,i2d)
!-----------------------------------------------------------------
!
! Computes the energy and helicity power 
! spectrum. The output is written to a 
! file by the first node.
!
! Parameters
!     a0 : input matrix: vortical modes
!     am : input matrix: wave modes (-)
!     ap : input matrix: wave modes (+)
!     bvfreq: Brunt-Vaisalla freq
!     omega : rotation rate
!     nmb: the extension used when writting the file
!     nmb1: if lenght>0, used to specify range
!     i2d : do 2D spectra (>0), or not (<=0)
!
!-----------------------------------------------------------------
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: E0k,EWk,EV0k,EP0k,EVWk,EPWk
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a0,am,ap
      REAL   (KIND=GP),             DIMENSION(n/2+1,n/2+1)   :: F0axi,FWaxi
      REAL   (KIND=GP), INTENT(IN)                           :: bvfreq,omega
      INTEGER         , INTENT(IN)                           :: i2d
      INTEGER                      :: j
      CHARACTER(len=*), INTENT(IN) :: dir,nmb,nmb1

!
!-----------------------------------------------------------------
!                        isotropic spectra:
!-----------------------------------------------------------------
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,1,1,E0k ,EWk  ) ! total energy
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,3,1,EV0k,EVWk ) ! KE (only EV0)
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,4,1,EP0k,EPWk ) ! PE 
      DO j = 1,n/2+1
        EVWk(j) = E0k(j)+EWk(j)-EV0k(j)-EP0k(j)-EPWk(j)
      ENDDO
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvekspectrum.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvekspectrum.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(6(E23.15,1X))') E0k(j),EWk(j),EV0k(j),EVWk(j),EP0k(j),EPWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,2,1,E0k ,EWk  ) ! Helicity 
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvhkspectrum.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvhkspectrum.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(E23.15,1X,E23.15)') E0k(j), EWk(j)
         ENDDO
         CLOSE(1)
      ENDIF

      IF ( i2d.GT.0 ) THEN
      CALL wvspecaxic(a0,am,ap,omega,bvfreq,1,F0axi,FWaxi) ! 2D axisymm spec E0, EW
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/kspec2DE0.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DE0.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEW.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEW.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) FWaxi
        CLOSE(1)
      ENDIF

      CALL wvspecaxic(a0,am,ap,omega,bvfreq,3,F0axi,FWaxi) ! 2D axisymm spec EV0 (only F0 filled)
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEV0.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEV0.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      CALL wvspecaxic(a0,am,ap,omega,bvfreq,4,F0axi,FWaxi) ! 2D axisymm spec EP0, EPW
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEP0.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEP0.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEPW.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DEPW.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) FWaxi
        CLOSE(1)
      ENDIF

      CALL wvspecaxic(a0,am,ap,omega,bvfreq,2,F0axi,FWaxi) ! 2D axisymm spec H0, HW
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DH0W.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DH0W.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DHW.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file=trim(dir) // '/' // 'kspec2DHW.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) FWaxi
        CLOSE(1)
      ENDIF
      ENDIF
!
!-----------------------------------------------------------------
!                        perpendicular spectra:
!-----------------------------------------------------------------
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,1,2,E0k ,EWk  ) ! total energy
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,3,2,EV0k,EVWk ) ! KE (only E0)
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,4,2,EP0k,EPWk ) ! PE 
      DO j = 1,n/2+1
        EVWk(j) = E0k(j)+EWk(j)-EV0k(j)-EP0k(j)-EPWk(j)
      ENDDO
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvekspecperp.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvekspecperp.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(6(E23.15,1X))') E0k(j),EWk(j),EV0k(j),EVWk(j),EP0k(j),EPWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,2,2,E0k ,EWk  ) ! Helicity 
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvhkspecperp.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvhkspecperp.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(E23.15,1X,E23.15)') E0k(j), EWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
!-----------------------------------------------------------------
!                        parallel spectra:
!-----------------------------------------------------------------
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,1,3,E0k ,EWk  ) ! total energy
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,3,3,EV0k,EVWk ) ! KE (only E0)
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,4,3,EP0k,EPWk ) ! PE 
      DO j = 1,n/2+1
        EVWk(j) = E0k(j)+EWk(j)-EV0k(j)-EP0k(j)-EPWk(j)
      ENDDO
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvekspecpara.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvekspecpara.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(6(E23.15,1X))') E0k(j),EWk(j),EV0k(j),EVWk(j),EP0k(j),EPWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      CALL wvspectrumc(a0,am,ap,omega,bvfreq,2,3,E0k ,EWk  ) ! Helicity 
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvhkspecpara.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvhkspecpara.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(E23.15,1X,E23.15)') E0k(j), EWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvspectrum


!*****************************************************************
      SUBROUTINE wvspecaxic(a0,am,ap,omega,bvfreq,kin,F0k,FWk)
!-----------------------------------------------------------------
!
! Computes axisymmetric (2d) wave and vortical spectra for various quantities, 
! returns them
!
! Parameters
!     a0   : input matri: vortical modes
!     am   : input matrix: wave modes (-) 
!     ap   : input matrix: wave modes (+) 
!     bvfreq: Brunt-Vaisalla frequency
!     f    : rotation parameter = 2\Omega
!     kin  : = 1 computes the total energy spectra wave/vortical
!            = 2 computes the helicity spectra alone wave/vortical 
!            = 3 computes the kinetic energy; vortical only; wave KE must be
!                computed from the outside as EVWk = E0k+EWk-EV0k-EP0k-EPWk,
!                using multiple calls to this routine; only F0k is filled then
!            = 4 computes the potential energy; wave/vortical
!     F0k  : vortical spectrum, returned
!     FWk  : wave spectrum, returned
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(OUT), DIMENSION(n/2+1,n/2+1)   :: F0k,FWk
      REAL   (KIND=GP),              DIMENSION(n/2+1,n/2+1)   :: E0k,EWk
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a0,am,ap
      COMPLEX(KIND=8 )                                        :: adel
      REAL   (KIND=GP), INTENT (IN)                           :: bvfreq,omega
      INTEGER, INTENT(IN)          :: kin
!     REAL   (KIND=GP)             :: ks,sig,tm0,tmi,tmr,tmw
      DOUBLE PRECISION             :: ks,kp2,sig,tm0,tmi,tmr,tmw
      REAL(KIND=GP)                :: f,tmp
      INTEGER                      :: i,ibeg,j,k,km
      INTEGER                      :: kperp,kpara

      f = 2.0*omega
      km      = n/2+1
      E0k     = 0.0
      EWk     = 0.0
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = ista
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.1 ) THEN ! Total energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kperp,kpara,tm0,tmw)
            DO j = 1,n
               kperp = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
               IF ( (kperp.lt.1).or.(kperp.gt.km) ) CYCLE
               DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  tm0     =     (abs(a0(k,j,1))**2 ) * tmp
                  tmw     =     (abs(am(k,j,1))**2 \
                          +      abs(ap(k,j,1))**2 ) * tmp
!$omp critical
                  E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                  EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-ibeg.ge.nth) private (j,k,kperp,kpara,tm0,tmw)
          DO i = ibeg,iend
!$omp parallel do if (iend-ibeg.lt.nth) private (k,kperp,kpara,tm0,tmw)
             DO j = 1,n
                kperp = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                IF ( (kperp.lt.1).or.(kperp.gt.km) ) CYCLE
                DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  tm0     = 2.0 * (abs(a0(k,j,i))**2 ) * tmp
                  tmw     = 2.0 * (abs(am(k,j,i))**2 \
                          +        abs(ap(k,j,i))**2 ) * tmp
!$omp critical
                  E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                  EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical

               END DO 
             END DO
          END DO
!
          CALL MPI_REDUCE(E0k,F0k,(n/2+1)*(n/2+1),GC_REAL,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(EWk,FWk,(n/2+1)*(n/2+1),GC_REAL,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.2 ) THEN ! Helicity spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,ks,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
            DO j = 1,n
               kp2   = ka(1)**2+ka(j)**2
               kperp = int(sqrt(kp2)+1.501)
               IF ( (kperp.gt.km) ) CYCLE
               DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0/(sqrt(2.0*(1.0+tmr)))
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     ks      = sqrt( ka2(k,j,1) )
                     adel    = conjg(ap(k,j,1))-conjg(am(k,j,1))
                     tm0     = dble( ks*a0(k,j,1) * adel *tmp*tmi )
                     tmw     = f*ka(k)*ks*( abs(am(k,j,1))**2 -   \
                               abs(ap(k,j,1))**2 )/sig * tmp
                             
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                     EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical
                  ELSE
                     tmw     = ka(k)*( abs(am(k,j,1))**2-abs(ap(k,j,1))**2 )*tmp
!$omp critical
                     EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,ks,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,ks,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
             DO j = 1,n
                kp2   = ka(i)**2+ka(j)**2
                kperp = int(sqrt(kp2)+1.501)
                IF ( (kperp.gt.km) ) CYCLE
                DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0/(sqrt(2.0*(1.0+tmr)))
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     ks      = sqrt ( ka2(k,j,i) ) 
                     adel    = conjg(ap(k,j,i))-conjg(am(k,j,i))
                     tm0     = 2.0D0*real( ks*a0(k,j,i) * adel *tmp*tmi )
                     tmw     = 2.0D0*f*ka(k)*ks*( abs(am(k,j,i))**2 -   \
                               abs(ap(k,j,i))**2 )/sig * tmp
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                     EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical
                  ELSE
                     tmw     = 2.0*ka(k)*( abs(am(k,j,i))**2-abs(ap(k,j,i))**2 )*tmp
!$omp critical
                     EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO

!         Compute reduction between nodes
!
          CALL MPI_REDUCE(E0k,F0k,(n/2+1)*(n/2+1),GC_REAL,      &    
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(EWk,FWk,(n/2+1)*(n/2+1),GC_REAL,      &    
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.3 ) THEN ! Kinetic energy spectra
         FWk = 0.0
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kperp,kpara,kp2,tm0,tmr)
            DO j = 1,n
               kp2   = ka(1)**2+ka(j)**2
               kperp = int(sqrt(kp2)+1.501)
               IF ( (kperp.gt.km) ) CYCLE
               DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tm0     = abs(a0(k,j,1))**2*tmp/(1.0+tmr)
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kperp,kpara,kp2,tm0,tmr)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kperp,kpara,kp2,tm0,tmr)
             DO j = 1,n
                kp2   = ka(i)**2+ka(j)**2
                kperp = int(sqrt(kp2)+1.501)
                IF ( (kperp.gt.km) ) CYCLE
                DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tm0     = 2.0*abs(a0(k,j,i))**2*tmp/(1.0+tmr)
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_REDUCE(E0k,F0k,(n/2+1)*(n/2+1),GC_REAL,      &    
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.4 ) THEN  ! Potential energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kperp,kpara,kp2,tm0,tmw,tmi,tmr)
            DO j = 1,n
               kp2   = ka(1)**2+ka(j)**2
               kperp = int(sqrt(kp2)+1.501)
               IF ( (kperp.gt.km) ) CYCLE
               DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0/(1.0+tmr)
                     tm0     = tmr*abs(a0(k,j,1))**2*tmp*tmi
                     tmw     = 0.50*abs(am(k,j,1)+ap(k,j,1))**2*tmp*tmi
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                     EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
                  ELSE
                     tm0     = abs(a0(k,j,1))**2*tmp
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
!$omp end critical
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kperp,kpara,kp2,tm0,tmw,tmi,tmr)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kperp,kpara,kp2,tm0,tmw,tmi,tmr)
             DO j = 1,n
                kp2   = ka(i)**2+ka(j)**2
                kperp = int(sqrt(kp2)+1.501)
                IF ( (kperp.gt.km) ) CYCLE
                DO k = 1,n
                  kpara = int(abs(ka(k))+1)
                  IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0/(1.0+tmr)
!$omp critical
                     tm0     = 2.0*tmr*abs(a0(k,j,i))**2*tmp*tmi
                     tmw     =         abs(am(k,j,i)+ap(k,j,i))**2*tmp*tmi
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
                     EWk(kperp,kpara) = EWk(kperp,kpara)+tmw
!$omp end critical
                  ELSE
                     tm0     = 2.0*abs(a0(k,j,i))**2*tmp
!$omp critical
                     E0k(kperp,kpara) = E0k(kperp,kpara)+tm0
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_REDUCE(E0k,F0k,(n/2+1)*(n/2+1),GC_REAL,      &    
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(EWk,FWk,(n/2+1)*(n/2+1),GC_REAL,      &    
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 4
!
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvspecaxic




!*****************************************************************
      SUBROUTINE wvspectrumc(a0,am,ap,omega,bvfreq,kin,kgeo,F0k,FWk)
!-----------------------------------------------------------------
!
! Computes 1d wave and vortical spectra for various quantities, returns them
!
! Parameters
!     a0   : input matrix: vortical modes
!     am   : input matrix: wave modes (-) 
!     ap   : input matrix: wave modes (+) 
!     bvfreq: Brunt-Vaisalla frequency
!     f    : rotation parameter = 2\Omega
!     kin  : = 1 computes the total energy spectra wave/vortical
!            = 2 computes the helicity spectra alone wave/vortical 
!            = 3 computes the kinetic energy; vortical only; wave KE must be
!                computed from the outside as EVWk = E0k+EWk-EV0k-EP0k-EPWk,
!                using multiple calls to this routine; only F0k is filled then
!            = 4 computes the potential energy; wave/vortical
!     kgeo : = 1, then do isotropic spectral dependence; 2, then do perpendicular 
!              spectral dependence; =3, then do parallel spectral dependence.
!     F0k  : vortical spectrum, returned
!     FWk  : wave spectrum, returned
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)         :: F0k,FWk
      DOUBLE PRECISION,  DIMENSION(n/2+1)                     :: E0k,EWk
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a0,am,ap
      REAL   (KIND=GP), INTENT (IN)                           :: bvfreq,omega
      INTEGER, INTENT(IN)          :: kin,kgeo
      DOUBLE PRECISION             :: ks,sig,tm0,tmi,tmr,tmw
      REAL(KIND=GP)                :: f,kp2,tmp
      INTEGER                      :: i,ibeg,j,k,kmsk(3)
      INTEGER                      :: kmn,kiso,kperp,kpara

      kmsk(1:3) = 0
      IF ( kgeo.lt. 1 .OR. kgeo.gt.3 ) THEN
        WRITE(*,*)'wvspectrumc: geometry parameter invalid.'
        STOP
      ENDIF
      kmsk(kgeo) = 1
      f = 2.0*omega
!
      E0k     = 0.0D0
      EWk     = 0.0D0
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = ista
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.1 ) THEN ! Total energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,tm0,tmr)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                  tm0     =     (abs(a0(k,j,1))**2 ) * tmp
                  tmw     =     (abs(am(k,j,1))**2 \
                          +      abs(ap(k,j,1))**2 ) * tmp
!$omp critical
                  E0k(kmn) = E0k(kmn)+tm0
                  EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-ibeg.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tm0,tmw)
          DO i = ibeg,iend
!$omp parallel do if (iend-ibeg.lt.nth) private (k,kmn,kiso,kperp,kpara,tm0,tmw)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE 
                  tm0     = 2.0D0 * (abs(a0(k,j,i))**2 ) * tmp
                  tmw     = 2.0D0 * (abs(am(k,j,i))**2 \
                          +        abs(ap(k,j,i))**2 ) * tmp
!$omp critical
                  E0k(kmn) = E0k(kmn)+tm0
                  EWk(kmn) = EWk(kmn)+tmw
!$omp end critical

               END DO 
             END DO
          END DO
!
          CALL MPI_REDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.2 ) THEN ! Helicity spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,ks,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kp2   = ka(1)**2+ka(j)**2
                  kperp = int(sqrt(kp2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0D0/(sqrt(2.0*(1.0D0+tmr)))
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     ks      = sqrt( ka(k)**2+kp2 )
                     tm0     = real( ks*a0(k,j,1) *  \
                             ( conjg(ap(k,j,1))-conjg(am(k,j,1)) ) )*tmp*tmi
                     tmw     = f*ka(k)*ks*( abs(am(k,j,1))**2 -   \
                               abs(ap(k,j,1))**2 )/sig * tmp
                             
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ELSE
                     tmw     = ka(k)*( abs(am(k,j,1))**2-abs(ap(k,j,1))**2 )*tmp
!$omp critical
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,ks,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,ks,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
             DO j = 1,n
                DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,i))+.501)
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kp2   = ka(i)**2+ka(j)**2
                  kperp = int(sqrt(kp2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE

                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0D0/(sqrt(2.0*(1.0D0+tmr)))
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     ks      = sqrt ( ka(k)**2+kp2 ) 
                     tm0     = 2.0D0*real( ks*a0(k,j,i) * \
                             ( conjg(ap(k,j,i))-conjg(am(k,j,i)) ) )*tmp*tmi
                     tmw     = 2.0D0*f*ka(k)*ks*( abs(am(k,j,i))**2 -   \
                               abs(ap(k,j,i))**2 )/sig * tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ELSE
                     tmw     = 2.0D0*ka(k)*( abs(am(k,j,i))**2-abs(ap(k,j,i))**2 )*tmp
!$omp critical
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO

!         Compute reduction between nodes
!
          CALL MPI_REDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.3 ) THEN ! Kinetic energy spectra
         FWk = 0.0
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmr)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kp2   = ka(1)**2+ka(j)**2
                  kperp = int(sqrt(kp2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tm0     = abs(a0(k,j,1))**2*tmp/(1.0D0+tmr)
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,tmr)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmr)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kp2   = ka(i)**2+ka(j)**2
                  kperp = int(sqrt(kp2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE

                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tm0     = 2.0D0*abs(a0(k,j,i))**2*tmp/(1.0D0+tmr)
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_REDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.4 ) THEN  ! Potential energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kp2   = ka(1)**2+ka(j)**2
                  kperp = int(sqrt(kp2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE

                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0D0/(1.0D0+tmr)
                     tm0     = tmr*abs(a0(k,j,1))**2*tmp*tmi
                     tmw     = 0.50D0*abs(am(k,j,1)+ap(k,j,1))**2*tmp*tmi
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
                  ELSE
                     tm0     = abs(a0(k,j,1))**2*tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
!$omp end critical
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kp2   = ka(i)**2+ka(j)**2
                  kperp = int(sqrt(kp2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE

                  IF ((kp2.gt.0.0)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0D0/(1.0D0+tmr)
!$omp critical
                     tm0     = 2.0D0*tmr*abs(a0(k,j,i))**2*tmp*tmi
                     tmw     =     abs(am(k,j,i)+ap(k,j,i))**2*tmp*tmi
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ELSE
                     tm0     = 2.0D0*abs(a0(k,j,i))**2*tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_REDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,0,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 4
!
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvspectrumc

!*****************************************************************
      SUBROUTINE sctransc(a,b,kgeo,Ektot)
!-----------------------------------------------------------------
!
! Computes the scalar transfer in Fourier space 
! in 3D and returns it.
!
! Parameters
!     a  : scalar
!     b  : nonlinear term
!     kgeo : = 1, then do isotropic spectral dependence; 2, then do perpendicular 
!              spectral dependence; =3, then do parallel spectral dependence.
!     Ektot: transfer function (output)
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX (KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      DOUBLE PRECISION,              DIMENSION(n/2+1) :: Ek
      DOUBLE PRECISION,  INTENT(OUT),DIMENSION(n/2+1) :: Ektot
      DOUBLE PRECISION :: tmq
      REAL(KIND=GP)    :: tmp
      INTEGER, INTENT(IN) :: kgeo
      INTEGER          :: i,j,k,kmsk(3)
      INTEGER          :: kmn,kiso,kperp,kpara

      kmsk(1:3) = 0
      IF ( kgeo.lt.1 .OR. kgeo.gt.3 ) THEN
        WRITE(*,*)'sctransc: geometry parameter invalid.'
        STOP
      ENDIF
      kmsk(kgeo) = 1
!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the scalar transfer
!
      tmp = 1./real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,tmq)
         DO j = 1,n
            DO k = 1,n
               kiso  = int(sqrt(ka2(k,j,1))+.501)
               kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
               kpara = int(abs(ka(k))+1)
               kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
!$omp atomic
                  Ek(kmn) = Ek(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,tmq)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,kiso,kperp,kpara,tmq)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE sctransc

!*****************************************************************
      SUBROUTINE sctrans2Dc(a,b,Ektot)
!-----------------------------------------------------------------
!
! Computes the axysimmetric scalar transfer in Fourier space
! and returns it.
! The spectrum is angle-averaged in the azimuthal direction, 
! and depends on two wavenumbers, kperp=0,...,n/2 and 
! kpara=0,....,n/2. The output is written to a binary file by 
! the first node.
!
! Parameters
!     a  : scalar
!     b  : nonlinear term
!     Ektot: transfer function (output)
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b
      REAL(KIND=GP),                DIMENSION(n/2+1,n/2+1)   :: Ek
      REAL(KIND=GP),    INTENT(OUT),DIMENSION(n/2+1,n/2+1)   :: Ektot
      REAL(KIND=GP)       :: tmq,tmp
      INTEGER             :: i,j,k
      INTEGER             :: kmn,kz

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         DO j = 1,n/2+1
            Ek(i,j) = 0.0_GP
         END DO
      END DO
! 
! Computes the scalar transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kz,kmn,tmq)
         DO j = 1,n
            kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               DO k = 1,n
                  kz = int(abs(ka(k))+1)
                  IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN                                          
                  tmq = tmp*real(a(k,j,1)*conjg(b(k,j,1)))
!$omp atomic
                  Ek(kmn,kz) = Ek(kmn,kz)+tmq
                  ENDIF
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kz,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     kz = int(abs(ka(k))+1)
                     IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
!$omp atomic
                     Ek(kmn,kz) = Ek(kmn,kz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kz,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     kz = int(abs(ka(k))+1)
                     IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                     tmq = 2*tmp*real(a(k,j,i)*conjg(b(k,j,i)))
!$omp atomic
                     Ek(kmn,kz) = Ek(kmn,kz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,(n/2+1)*(n/2+1),GC_REAL,            &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE sctrans2Dc


!*****************************************************************
      SUBROUTINE entransc(a,b,c,d,e,f,kin,kgeo,Ektot)
!-----------------------------------------------------------------
!
! Computes the energy (or cross-helicity) transfer 
! in Fourier space in 3D and returns it.
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force work (energy transfer)
!          =3 computes the magnetic cross-helicity transfer
!          =4 computes the kinetic cross-helicity transfer
!     kgeo : = 1, then do isotropic spectral dependence; 2, then do perpendicular 
!              spectral dependence; =3, then do parallel spectral dependence.
!     Ektot: transfer function, output
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION,             DIMENSION(n/2+1) :: Ek
      DOUBLE PRECISION, INTENT(OUT),DIMENSION(n/2+1) :: Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin,kgeo
      INTEGER             :: i,j,k,kmsk(3)
      INTEGER             :: kmn,kiso,kperp,kpara

      kmsk(1:3) = 0
      IF ( kgeo.lt.1 .OR. kgeo.gt.3 ) THEN
        WRITE(*,*)'entransc: geometry parameter invalid.'
        STOP
      ENDIF
      kmsk(kgeo) = 1
!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (kin.ge.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,tmq)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+            &
                            real(b(k,j,1)*conjg(e(k,j,1)))+            &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kiso  = int(sqrt(ka2(k,j,i))+.501)
                     kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                     kpara = int(abs(ka(k))+1)
                     kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,kiso,kperp,kpara,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kiso  = int(sqrt(ka2(k,j,i))+.501)
                     kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                     kpara = int(abs(ka(k))+1)
                     kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy transfer
!
      ELSE
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,tmq)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     tmq = ka2(k,j,1)*                                 &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+            &
                            real(b(k,j,1)*conjg(e(k,j,1)))+            &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kiso  = int(sqrt(ka2(k,j,i))+.501)
                     kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                     kpara = int(abs(ka(k))+1)
                     kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*ka2(k,j,i)*                            &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+         &
                               real(b(k,j,i)*conjg(e(k,j,i)))+         &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic 
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,kiso,kperp,kpara,tmq)
               DO j = 1,n
                  DO k = 1,n
                     kiso  = int(sqrt(ka2(k,j,i))+.501)
                     kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                     kpara = int(abs(ka(k))+1)
                     kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                     IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                        tmq = 2*ka2(k,j,i)*                            &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+         &
                               real(b(k,j,i)*conjg(e(k,j,i)))+         &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE entransc

!*****************************************************************
      SUBROUTINE entrans2Dc(a,b,c,d,e,f,kin,Ektot)
!-----------------------------------------------------------------
!
! Computes the axysimmetric energy transfer in Fourier space
! and returns it. 
! The spectrum is angle-averaged in the azimuthal direction, 
! and depends on two wavenumbers, kperp=0,...,n/2 and 
! kpara=0,....,n/2. The output is written to a binary file by 
! the first node.
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force work (energy transfer)
!          =3 computes the magnetic cross-helicity transfer
!          =4 computes the kinetic cross-helicity transfer
!     Ektot: transfer function, output
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: d,e,f
      REAL(KIND=GP),                DIMENSION(n/2+1,n/2+1)   :: Ek
      REAL(KIND=GP),    INTENT(OUT),DIMENSION(n/2+1,n/2+1)   :: Ektot
      REAL(KIND=GP)       :: tmq,tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn,kz

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         DO j = 1,n/2+1
            Ek(i,j) = 0.0_GP
         END DO
      END DO
! 
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/real(n,kind=GP)**6
      IF (kin.ge.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     kz = int(abs(ka(k))+1)
                     IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN                                          
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+            &
                            real(b(k,j,1)*conjg(e(k,j,1)))+            &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn,kz) = Ek(kmn,kz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+       &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy transfer
!
      ELSE
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kz,kmn,tmq)
            DO j = 1,n
               kmn = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  DO k = 1,n
                     kz = int(abs(ka(k))+1)
                     IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                     tmq = ka2(k,j,1)*                                 &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+            &
                            real(b(k,j,1)*conjg(e(k,j,1)))+            &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn,kz) = Ek(kmn,kz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*ka2(k,j,i)*                            &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+         &
                               real(b(k,j,i)*conjg(e(k,j,i)))+         &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kz,kmn,tmq)
               DO j = 1,n
                  kmn = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
                  IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                     DO k = 1,n
                        kz = int(abs(ka(k))+1)
                        IF ((kz.gt.0).and.(kz.le.n/2+1)) THEN
                        tmq = 2*ka2(k,j,i)*                            &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+         &
                               real(b(k,j,i)*conjg(e(k,j,i)))+         &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kz) = Ek(kmn,kz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(Ek,Ektot,(n/2+1)*(n/2+1),GC_REAL,            &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE entrans2Dc


!*****************************************************************
      SUBROUTINE wvtrans(a0,am,ap,vx,vy,vz,th,omega,bvfreq,dir,nmb,nmb1,i2d,c1,c2,c3)
!-----------------------------------------------------------------
!
! Computes the energy (and helicity) transfer 
! functions. The output is written to a 
! file by the first node.
!
! Parameters
!     a0 : input matrix: vortical modes
!     am : input matrix: wave modes (-)
!     ap : input matrix: wave modes (+)
!     bvfreq: Brunt-Vaisalla freq
!     omega : rotation rate
!     nmb: the extension used when writting the file
!     nmb1: if lenght>0, used to specify range
!     i2d : do 2D spectra (>0), or not (<=0)
!     c1,c2,c3: tmp complex arrays
!
!-----------------------------------------------------------------
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TV0k,TP0k,TVWk,TPWk,GWk,TV0kperp,TP0kperp,TVWkperp,TPWkperp,GWkperp,TV0kpar,TP0kpar,TVWkpar,TPWkpar,GWkpar
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a0,am,ap,vx,vy,vz,th
      COMPLEX(KIND=GP),             DIMENSION(n,n,ista:iend) :: vxt,vyt,vzt
      COMPLEX(KIND=GP), INTENT(INOUT),DIMENSION(n,n,ista:iend) :: c1,c2,c3
      REAL   (KIND=GP),             DIMENSION(n/2+1,n/2+1)   :: Taxi
      REAL   (KIND=GP), INTENT(IN)                           :: bvfreq,omega
      INTEGER         , INTENT(IN)                           :: i2d
      INTEGER                      :: j,i
      CHARACTER(len=*), INTENT(IN) :: dir,nmb,nmb1

      CALL prodre3(vx,vy,vz,vxt,vyt,vzt)
      CALL nonlhd3(vxt,vyt,vzt,c1,1)
      CALL nonlhd3(vxt,vyt,vzt,c2,2)
      CALL nonlhd3(vxt,vyt,vzt,c3,3)

!      CALL wvproj0(a0,vxt,vyt,vzt,tht,omega,bvfreq)
      CALL wvprojv0(a0,vxt,vyt,vzt,omega,bvfreq)
      CALL entransc(vxt,vyt,vzt,c1,c2,c3,1,1,TV0k)      ! Isotropic transfer function     TV0
      CALL entransc(vxt,vyt,vzt,c1,c2,c3,1,2,TV0kperp)  ! Perpendicular transfer function TV0
      CALL entransc(vxt,vyt,vzt,c1,c2,c3,1,3,TV0kpar)   ! Parallel transfer function      TV0
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(vxt,vyt,vzt,c1,c2,c3,1,Taxi)   ! 2D transfer function            TV0
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvev0ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvev0ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL advect3(vx,vy,vz,th,vxt)
      CALL wvprojt0(a0,vyt,omega,bvfreq)
      CALL sctransc(vyt,vxt,1,TP0k)                     ! Isotropic transfer function     TP0
      CALL sctransc(vyt,vxt,2,TP0kperp)                 ! Perpendicular transfer function TP0
      CALL sctransc(vyt,vxt,3,TP0kpar)                  ! Parallel transfer function      TP0
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(vyt,vxt,Taxi)                  ! 2D transfer function            TP0
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvep0ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvep0ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

!      CALL wvprojw(am,ap,vxt,vyt,vzt,tht,omega,bvfreq)
      CALL wvprojvw(am,ap,vxt,vyt,vzt,omega,bvfreq)
      CALL entransc(vxt,vyt,vzt,c1,c2,c3,1,1,TVWk)      ! Isotropic transfer function     TVW
      CALL entransc(vxt,vyt,vzt,c1,c2,c3,1,2,TVWkperp)  ! Perpendicular transfer function TVW
      CALL entransc(vxt,vyt,vzt,c1,c2,c3,1,3,TVWkpar)   ! Parallel transfer function      TVW
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(vxt,vyt,vzt,c1,c2,c3,1,Taxi)   ! 2D transfer function            TVW
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvevwktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvevwktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL advect3(vx,vy,vz,th,c1)
      CALL wvprojtw(am,ap,c2,omega,bvfreq)
      CALL sctransc(c2,c1,1,TPWk)                     ! Isotropic transfer function     TPW
      CALL sctransc(c2,c1,2,TPWkperp)                 ! Perpendicular transfer function TPW
      CALL sctransc(c2,c1,3,TPWkpar)                  ! Parallel transfer function      TPW
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c2,c1,Taxi)                  ! 2D transfer function            TPW
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvepwktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvepwktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

! Compute the wave exchange terms between kinetic and potential energy:
! (sign convention is that of the EP budget)
      CALL specprodg(vzt,th,1,GWk)     ! Isotropic
      CALL specprodg(vzt,th,2,GWkperp) ! Perpendicular
      CALL specprodg(vzt,th,3,GWkpar)  ! Parallel

      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
            OPEN(1,file='wvektransfer.' // nmb // '_' // trim(nmb1) //'.txt')
         else
            OPEN(1,file='wvektransfer.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
            WRITE(1,FMT='(7(E23.15,1X))') TV0k(j)+TP0k(j),TVWk(j)+TPWk(j),TV0k(j),TVWk(j),TP0k(j),TPWk(j),bvfreq*GWk(j)
         ENDDO
         CLOSE(1)

         if ( len_trim(nmb1).gt.0 ) then
            OPEN(1,file='wvektranperp.' // nmb // '_' // trim(nmb1) //'.txt')
         else
            OPEN(1,file='wvektranperp.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
            WRITE(1,FMT='(7(E23.15,1X))') TV0kperp(j)+TP0kperp(j),TVWkperp(j)+TPWkperp(j),TV0kperp(j),TVWkperp(j),TP0kperp(j),TPWkperp(j),bvfreq*GWkperp(j)
         ENDDO
         CLOSE(1)

         if ( len_trim(nmb1).gt.0 ) then
            OPEN(1,file='wvektranpara.' // nmb // '_' // trim(nmb1) //'.txt')
         else
            OPEN(1,file='wvektranpara.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
            WRITE(1,FMT='(7(E23.15,1X))') TV0kpar(j)+TP0kpar(j),TVWkpar(j)+TPWkpar(j),TV0kpar(j),TVWkpar(j),TP0kpar(j),TPWkpar(j),bvfreq*GWkpar(j)
         ENDDO
         CLOSE(1)
      ENDIF

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvtrans

      
!*****************************************************************
      SUBROUTINE vector3(a,b,c,d,e,f,x,y,z)
!-----------------------------------------------------------------
!
! Computes the product AxB in real space. The components 
! of the vector fields A and B are given by the matrixes 
! a,b,c,d,e and f, following the right hand convention.
!
! Parameters
!     a: input matrix with A_x
!     b: input matrix with A_y
!     c: input matrix with A_z
!     d: input matrix with B_x
!     e: input matrix with B_y
!     f: input matrix with B_z
!     x: at the output contains (AxB)_x in Fourier space
!     y: at the output contains (AxB)_y in Fourier space
!     z: at the output contains (AxB)_z in Fourier space
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend)  :: d,e,f
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,n,ista:iend) :: x,y,z
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r3,r4
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r5,r6
      REAL(KIND=GP), DIMENSION(n,n,ksta:kend) :: r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               x(k,j,i) = a(k,j,i)
               y(k,j,i) = b(k,j,i)
               z(k,j,i) = c(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,x,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,y,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,z,r3,MPI_COMM_WORLD)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,n
            DO k = 1,n
               x(k,j,i) = d(k,j,i)
               y(k,j,i) = e(k,j,i)
               z(k,j,i) = f(k,j,i)
            END DO
         END DO
      END DO
      CALL fftp3d_complex_to_real(plancr,x,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,y,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(plancr,z,r6,MPI_COMM_WORLD)

      tmp = 1.0_GP/real(n,kind=GP)**6
!$omp parallel do if (iend-ista.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (iend-ista.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,r7,x,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r3,y,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,r1,z,MPI_COMM_WORLD)

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE vector3
      
      
!*****************************************************************
      SUBROUTINE wvtransfull(a0,am,ap,vx,vy,vz,th,omega,bvfreq,dir,nmb,nmb1,i2d,c1,c2,c3)
!-----------------------------------------------------------------
!
! Computes the energy (and helicity) transfer 
! functions, distinguishing the components by types of triads.
! The output is written to a file by the first node.
!
! Parameters
!     a0 : input matrix: vortical modes
!     am : input matrix: wave modes (-)
!     ap : input matrix: wave modes (+)
!     bvfreq: Brunt-Vaisalla freq
!     omega : rotation rate
!     nmb: the extension used when writting the file
!     nmb1: if lenght>0, used to specify range
!     i2d : do 2D spectra (>0), or not (<=0)
!     c1,c2,c3: tmp complex arrays
!     CAUTION: vx,vy,vz,th are overwritten!        
!
!-----------------------------------------------------------------
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: GW0k,GW0kperp,GW0kpar,GWWk,GWWkperp,GWWkpar
      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TV000k,    TV00Wk,    TV0W0k,    TV0WWk,    TVW00k,    TVW0Wk,    TVWW0k,    TVWWWk
      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TV000kperp,TV00Wkperp,TV0W0kperp,TV0WWkperp,TVW00kperp,TVW0Wkperp,TVWW0kperp,TVWWWkperp
      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TV000kpar, TV00Wkpar, TV0W0kpar, TV0WWkpar, TVW00kpar, TVW0Wkpar, TVWW0kpar, TVWWWkpar
      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TP000k,    TP00Wk,    TP0W0k,    TP0WWk,    TPW00k,    TPW0Wk,    TPWW0k,    TPWWWk
      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TP000kperp,TP00Wkperp,TP0W0kperp,TP0WWkperp,TPW00kperp,TPW0Wkperp,TPWW0kperp,TPWWWkperp
      DOUBLE PRECISION,             DIMENSION(n/2+1)         :: TP000kpar, TP00Wkpar, TP0W0kpar, TP0WWkpar, TPW00kpar, TPW0Wkpar, TPWW0kpar, TPWWWkpar
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a0,am,ap
      COMPLEX(KIND=GP), INTENT(INOUT),DIMENSION(n,n,ista:iend) :: vx,vy,vz,th      
      COMPLEX(KIND=GP),             DIMENSION(n,n,ista:iend) :: c11,c22,c33,c44,c55,c66,c77
      COMPLEX(KIND=GP), INTENT(INOUT),DIMENSION(n,n,ista:iend) :: c1,c2,c3
      REAL   (KIND=GP),             DIMENSION(n/2+1,n/2+1)   :: Taxi
      REAL   (KIND=GP), INTENT(IN)                           :: bvfreq,omega
      INTEGER         , INTENT(IN)                           :: i2d
      INTEGER                      :: j,i
      CHARACTER(len=*), INTENT(IN) :: dir,nmb,nmb1

      CALL wvprojv0(a0,vx,vy,vz,omega,bvfreq) ! u0

      CALL rotor3(vy,vz,c1,1) ! omega0
      CALL rotor3(vx,vz,c2,2)
      CALL rotor3(vx,vy,c3,3)

      CALL vector3(c1,c2,c3,vx,vy,vz,c33,c44,th) ! nonlinear term 00
      CALL nonlhd3(c33,c44,th,c11,1)
      CALL nonlhd3(c33,c44,th,c22,2)
      CALL nonlhd3(c33,c44,th,c33,3)
      
      CALL entransc(vx,vy,vz,c11,c22,c33,1,1,TV000k)      ! Isotropic transfer function     TV000
      CALL entransc(vx,vy,vz,c11,c22,c33,1,2,TV000kperp)  ! Perpendicular transfer function TV000
      CALL entransc(vx,vy,vz,c11,c22,c33,1,3,TV000kpar)   ! Parallel transfer function      TV000
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(vx,vy,vz,c11,c22,c33,1,Taxi)   ! 2D transfer function            TV000
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvev000ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvev000ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL wvprojvw(am,ap,c44,c55,c66,omega,bvfreq) ! uw
      CALL entransc(c44,c55,c66,c11,c22,c33,1,1,TVW00k)      ! Isotropic transfer function     TVW00
      CALL entransc(c44,c55,c66,c11,c22,c33,1,2,TVW00kperp)  ! Perpendicular transfer function TVW00
      CALL entransc(c44,c55,c66,c11,c22,c33,1,3,TVW00kpar)   ! Parallel transfer function      TVW00
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(c44,c55,c66,c11,c22,c33,1,Taxi)   ! 2D transfer function            TVW00
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvevw00ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvevw00ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF
      
      CALL vector3(c1,c2,c3,c44,c55,c66,c77,th,c33) ! nonlinear term 0W
      CALL nonlhd3(c77,th,c33,c11,1)
      CALL nonlhd3(c77,th,c33,c22,2)
      CALL nonlhd3(c77,th,c33,c33,3)

      CALL entransc(vx,vy,vz,c11,c22,c33,1,1,TV00Wk)      ! Isotropic transfer function     TV00W
      CALL entransc(vx,vy,vz,c11,c22,c33,1,2,TV00Wkperp)  ! Perpendicular transfer function TV00W
      CALL entransc(vx,vy,vz,c11,c22,c33,1,3,TV00Wkpar)   ! Parallel transfer function      TV00W
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(vx,vy,vz,c11,c22,c33,1,Taxi)   ! 2D transfer function            TV00W
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvev00wktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvev00wktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF


      CALL entransc(c44,c55,c66,c11,c22,c33,1,1,TVW0Wk)      ! Isotropic transfer function     TVW0W
      CALL entransc(c44,c55,c66,c11,c22,c33,1,2,TVW0Wkperp)  ! Perpendicular transfer function TVW0W
      CALL entransc(c44,c55,c66,c11,c22,c33,1,3,TVW0Wkpar)   ! Parallel transfer function      TVW0W
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(c44,c55,c66,c11,c22,c33,1,Taxi)   ! 2D transfer function            TVW0W
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvevw0wktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvevw0wktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      
! now I don't need omega0 anymore so arrays c1,c2,c3 are free to be used elsewhere
      CALL rotor3(c55,c66,c1,1) ! omegaW
      CALL rotor3(c44,c66,c2,2)
      CALL rotor3(c44,c55,c3,3)

      CALL vector3(c1,c2,c3,vx,vy,vz,c77,th,c33) ! nonlinear term W0
      CALL nonlhd3(c77,th,c33,c11,1)
      CALL nonlhd3(c77,th,c33,c22,2)
      CALL nonlhd3(c77,th,c33,c33,3)
      
      CALL entransc(vx,vy,vz,c11,c22,c33,1,1,TV0W0k)      ! Isotropic transfer function     TV0W0
      CALL entransc(vx,vy,vz,c11,c22,c33,1,2,TV0W0kperp)  ! Perpendicular transfer function TV0W0
      CALL entransc(vx,vy,vz,c11,c22,c33,1,3,TV0W0kpar)   ! Parallel transfer function      TV0W0
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(vx,vy,vz,c11,c22,c33,1,Taxi)   ! 2D transfer function            TV0W0
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvev0w0ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvev0w0ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF
      
      CALL entransc(c44,c55,c66,c11,c22,c33,1,1,TVWW0k)      ! Isotropic transfer function     TVWW0
      CALL entransc(c44,c55,c66,c11,c22,c33,1,2,TVWW0kperp)  ! Perpendicular transfer function TVWW0
      CALL entransc(c44,c55,c66,c11,c22,c33,1,3,TVWW0kpar)   ! Parallel transfer function      TVWW0
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(c44,c55,c66,c11,c22,c33,1,Taxi)   ! 2D transfer function            TVWW0
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvevww0ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvevww0ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF
      
      CALL vector3(c1,c2,c3,c44,c55,c66,c77,th,c33) ! nonlinear term WW
      CALL nonlhd3(c77,th,c33,c11,1)
      CALL nonlhd3(c77,th,c33,c22,2)
      CALL nonlhd3(c77,th,c33,c33,3)

      CALL entransc(vx,vy,vz,c11,c22,c33,1,1,TV0WWk)      ! Isotropic transfer function     TV0WW
      CALL entransc(vx,vy,vz,c11,c22,c33,1,2,TV0WWkperp)  ! Perpendicular transfer function TV0WW
      CALL entransc(vx,vy,vz,c11,c22,c33,1,3,TV0WWkpar)   ! Parallel transfer function      TV0WW
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(vx,vy,vz,c11,c22,c33,1,Taxi)   ! 2D transfer function            TV0WW
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvev0wwktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvev0wwktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL entransc(c44,c55,c66,c11,c22,c33,1,1,TVWWWk)      ! Isotropic transfer function     TVWWW
      CALL entransc(c44,c55,c66,c11,c22,c33,1,2,TVWWWkperp)  ! Perpendicular transfer function TVWWW
      CALL entransc(c44,c55,c66,c11,c22,c33,1,3,TVWWWkpar)   ! Parallel transfer function      TVWWW
      IF ( i2d.GT.0 ) THEN
         CALL entrans2Dc(c44,c55,c66,c11,c22,c33,1,Taxi)   ! 2D transfer function            TVWWW
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvevwwwktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvevwwwktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF


! Now we do the potential energy transfer functions
      CALL wvprojt0(a0,th,omega,bvfreq) ! th0
      
      CALL advect3(vx,vy,vz,th,c1)
      CALL sctransc(c1,th,1,TP000k)                     ! Isotropic transfer function     TP000
      CALL sctransc(c1,th,2,TP000kperp)                 ! Perpendicular transfer function TP000
      CALL sctransc(c1,th,3,TP000kpar)                  ! Parallel transfer function      TP000
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,th,Taxi)                  ! 2D transfer function            TP000
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvep000ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvep000ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF
      
      CALL wvprojtw(am,ap,c77,omega,bvfreq) ! thW
      CALL sctransc(c1,c77,1,TPW00k)                     ! Isotropic transfer function     TPW00
      CALL sctransc(c1,c77,2,TPW00kperp)                 ! Perpendicular transfer function TPW00
      CALL sctransc(c1,c77,3,TPW00kpar)                  ! Parallel transfer function      TPW00
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,c77,Taxi)                  ! 2D transfer function            TPW00
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvepw00ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvepw00ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF


      CALL advect3(c44,c55,c66,th,c1)
      CALL sctransc(c1,th,1,TP0W0k)                     ! Isotropic transfer function     TP0W0
      CALL sctransc(c1,th,2,TP0W0kperp)                 ! Perpendicular transfer function TP0W0
      CALL sctransc(c1,th,3,TP0W0kpar)                  ! Parallel transfer function      TP0W0
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,th,Taxi)                  ! 2D transfer function            TP0W0
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvep0w0ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvep0w0ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL sctransc(c1,c77,1,TPWW0k)                     ! Isotropic transfer function     TPWW0
      CALL sctransc(c1,c77,2,TPWW0kperp)                 ! Perpendicular transfer function TPWW0
      CALL sctransc(c1,c77,3,TPWW0kpar)                  ! Parallel transfer function      TPWW0
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,c77,Taxi)                  ! 2D transfer function            TPWW0
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvepww0ktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvepww0ktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF


      CALL advect3(vx,vy,vz,c77,c1)
      CALL sctransc(c1,th,1,TP00Wk)                     ! Isotropic transfer function     TP00W
      CALL sctransc(c1,th,2,TP00Wkperp)                 ! Perpendicular transfer function TP00W
      CALL sctransc(c1,th,3,TP00Wkpar)                  ! Parallel transfer function      TP00W
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,th,Taxi)                  ! 2D transfer function            TP00W
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvep00wktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvep00wktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL sctransc(c1,c77,1,TPW0Wk)                     ! Isotropic transfer function     TPW0W
      CALL sctransc(c1,c77,2,TPW0Wkperp)                 ! Perpendicular transfer function TPW0W
      CALL sctransc(c1,c77,3,TPW0Wkpar)                  ! Parallel transfer function      TPW0W
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,c77,Taxi)                  ! 2D transfer function            TPW0W
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvepw0wktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvepw0wktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF


      CALL advect3(c44,c55,c66,c77,c1)
      CALL sctransc(c1,th,1,TP0WWk)                     ! Isotropic transfer function     TP0WW
      CALL sctransc(c1,th,2,TP0WWkperp)                 ! Perpendicular transfer function TP0WW
      CALL sctransc(c1,th,3,TP0WWkpar)                  ! Parallel transfer function      TP0WW
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,th,Taxi)                  ! 2D transfer function            TP0WW
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvep0wwktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvep0wwktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

      CALL sctransc(c1,c77,1,TPWWWk)                     ! Isotropic transfer function     TPWWW
      CALL sctransc(c1,c77,2,TPWWWkperp)                 ! Perpendicular transfer function TPWWW
      CALL sctransc(c1,c77,3,TPWWWkpar)                  ! Parallel transfer function      TPWWW
      IF ( i2d.GT.0 ) THEN
         CALL sctrans2Dc(c1,c77,Taxi)                  ! 2D transfer function            TPWWW
         IF (myrank.eq.0) THEN 
            if ( len_trim(nmb1).gt.0 ) then 
               OPEN(1,file=trim(dir) // '/' // 'wvepwwwktrans2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
            else 
               OPEN(1,file=trim(dir) // '/' // 'wvepwwwktrans2D.' // nmb // '.out',form='unformatted',access='stream')
            endif
            WRITE(1) Taxi
            CLOSE(1)
         ENDIF
      ENDIF

! Compute the wave exchange terms between kinetic and potential energy:
! (sign convention is that of the EP budget)
      CALL specprodg(c66,th,1,GW0k)     ! Isotropic
      CALL specprodg(c66,th,2,GW0kperp) ! Perpendicular
      CALL specprodg(c66,th,3,GW0kpar)  ! Parallel
      CALL specprodg(c66,c77,1,GWWk)     ! Isotropic
      CALL specprodg(c66,c77,2,GWWkperp) ! Perpendicular
      CALL specprodg(c66,c77,3,GWWkpar)  ! Parallel


      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
            OPEN(1,file='wvektransfer.' // nmb // '_' // trim(nmb1) //'.txt')
         else
            OPEN(1,file='wvektransfer.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
            WRITE(1,FMT='(18(E23.15,1X))') TV000k(j),TV00Wk(j),TV0W0k(j),TV0WWk(j),TVW00k(j),TVW0Wk(j),TVWW0k(j),TVWWWk(j), &
                 TP000k(j),TP00Wk(j),TP0W0k(j),TP0WWk(j),TPW00k(j),TPW0Wk(j),TPWW0k(j),TPWWWk(j), &
                 bvfreq*GW0k(j),bvfreq*GWWk(j)
         ENDDO
         CLOSE(1)

         if ( len_trim(nmb1).gt.0 ) then
            OPEN(1,file='wvektranperp.' // nmb // '_' // trim(nmb1) //'.txt')
         else
            OPEN(1,file='wvektranperp.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
            WRITE(1,FMT='(18(E23.15,1X))') TV000kperp(j),TV00Wkperp(j),TV0W0kperp(j),TV0WWkperp(j),TVW00kperp(j),TVW0Wkperp(j),TVWW0kperp(j),TVWWWkperp(j), &
                 TP000kperp(j),TP00Wkperp(j),TP0W0kperp(j),TP0WWkperp(j),TPW00kperp(j),TPW0Wkperp(j),TPWW0kperp(j),TPWWWkperp(j), &
                 bvfreq*GW0kperp(j),bvfreq*GWWkperp(j)
         ENDDO
         CLOSE(1)

         if ( len_trim(nmb1).gt.0 ) then
            OPEN(1,file='wvektranpara.' // nmb // '_' // trim(nmb1) //'.txt')
         else
            OPEN(1,file='wvektranpara.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
            WRITE(1,FMT='(18(E23.15,1X))') TV000kpar(j),TV00Wkpar(j),TV0W0kpar(j),TV0WWkpar(j),TVW00kpar(j),TVW0Wkpar(j),TVWW0kpar(j),TVWWWkpar(j), &
                 TP000kpar(j),TP00Wkpar(j),TP0W0kpar(j),TP0WWkpar(j),TPW00kpar(j),TPW0Wkpar(j),TPWW0kpar(j),TPWWWkpar(j), &
                 bvfreq*GW0kpar(j),bvfreq*GWWkpar(j)
         ENDDO
         CLOSE(1)
      ENDIF

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvtransfull
      

!*****************************************************************
      SUBROUTINE DoSpAvg(vx,vy,vz,th,istat,nstat,idir,planrc,planio,rv,c1)
!-----------------------------------------------------------------
!
! Computes time average of fluctuation fields.
!
! Parameters
!     vx,vy,vz,th  : complex arrays averaged over istat time stamps
!     istat        : array of time stamp integers
!     nstat        : no. istat timestamps
!     iswap        : do byte swap on input?
!     planrc       : real to complex plan for FFT
!     planio       : IO plan
!     rv           : real tmp array
!     c1           : complex tmp array, size of vx, etc.
!
!-----------------------------------------------------------------
      USE fftplans
      USE iovar
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      USE filefmt
      USE gutils
!$    USE threads
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(IN)                                 :: planrc
      TYPE (IOPLAN), INTENT(IN)                                 :: planio
      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: rv
      REAL   (KIND=GP)                                          :: tmp
      INTEGER         , INTENT   (IN), DIMENSION            (*) :: istat
      INTEGER         , INTENT   (IN)                           :: nstat
      INTEGER                                                   :: i,j,k
      CHARACTER(len=*), INTENT   (IN)                           :: idir

      vx = 0.0_GP
      vy = 0.0_GP
      vz = 0.0_GP
      th = 0.0_GP
      DO i = 1,nstat
        WRITE(ext, fmtext) istat(i)
! read in appropriate file:
        CALL io_read(1,idir,'vx',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,c1,MPI_COMM_WORLD)
        CALL padd(vx,c1,1.0)

        CALL io_read(1,idir,'vy',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,c1,MPI_COMM_WORLD)
        CALL padd(vy,c1,1.0)

        CALL io_read(1,idir,'vz',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,c1,MPI_COMM_WORLD)
        CALL padd(vz,c1,1.0)

        CALL io_read(1,idir,'th',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,c1,MPI_COMM_WORLD)
        CALL padd(th,c1,1.0)
      ENDDO
     
       IF ( nstat.EQ. 1) return 
       tmp = 1.0 / real(nstat,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
       DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,n
             DO k = 1,n
                vx(k,j,i) = vx(k,j,i)*tmp
                vy(k,j,i) = vy(k,j,i)*tmp
                vz(k,j,i) = vz(k,j,i)*tmp
                th(k,j,i) = th(k,j,i)*tmp
             END DO
          END DO
       END DO


!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE DoSpAvg


!*****************************************************************
      SUBROUTINE padd(v,p,c)
!-----------------------------------------------------------------
!
! Computes v = v + c*p
!
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: v
      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: p
      REAL   (KIND=GP), INTENT   (IN)                           :: c
      INTEGER                                                   :: i,j,k

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
       DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
          DO j = 1,n
             DO k = 1,n
!$omp atomic
                v(k,j,i) = v(k,j,i) + c*p(k,j,i)
             END DO
          END DO
       END DO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE padd

!
      SUBROUTINE PVn(gn,a0,vx,vy,vz,th,omega,bvfreq,nt,c1,r1,r2,r3)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes nt-order contribution to potential vorticity
!
! Parameters
!     gn: nt-order potential vorticity component (returned)
!     a0: complex temp array of size vx,vy,vz, containing vortical modes
!     vx,
!     vy,
!     vz    : complex velocities, overwritten with normal mode fields
!     th    : complex potential temperature
!     omega : rotation rate
!     bvfreq: Brunt-Vaisalla frequency
!     nt    : order of component (1,2). 
!     c1  : complex tmp array
!     r1-3  : real tmp array
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: gn
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1,r2,r3
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ksigk,tmp
      COMPLEX(KIND=GP)                                          :: ic
      INTEGER         , INTENT (IN)                             :: nt
      INTEGER                                                   :: i,j,k

      IF ( nt.lt.1 .or. nt.gt.2 ) THEN
        STOP 'PVn: invalid order'
      ENDIF


      f = 2.0*omega
      gn = (0.0_GP,0.0_GP)
      ic = cmplx(0.0_GP,1.0_GP);

      IF ( nt.eq.1 ) THEN
        tmp = 1.0_GP/real(n,kind=GP)**6
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,n
              kp  = sqrt(ka(i)**2+ka(j)**2)
              DO k = 1,n
                ksigk = sqrt(f**2*ka(k)**2+bvfreq**2*kp**2)
                gn(k,j,i) = -ic*ksigk*a0(k,j,i)
              END DO
           END DO
        END DO
      ENDIF

      IF ( nt.eq.2 ) THEN
        CALL wgradt(gn,vx,vy,vz,th,c1,r1,r2,r3)
      ENDIF

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE PVn


      SUBROUTINE wgradt(gn,vx,vy,vz,th,c1,r1,r2,r3)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes omega.Grad theta, where
!
!  omega = curl(v), using little memory
!
! Parameters
!     gn: product (returned)
!     th    : complex temp.
!     vx,
!     vy,
!     vz    : complex velocities
!     c1    : complex tmp array
!     r1-3  : real tmp array
!
!-----------------------------------------------------------------
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


      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: gn
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1,r2,r3
      REAL   (KIND=GP)                                          :: tmp
      INTEGER                                                   :: i,j,k

      tmp = 1.0_GP/real(n,kind=GP)**6
      CALL derivk3(th,gn,1)
      CALL fftp3d_complex_to_real(plancr,gn,r1,MPI_COMM_WORLD)
      CALL rotor3(vy,vz,c1,1)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = r1(i,j,k)*r2(i,j,k)*tmp
            END DO
         END DO
      END DO

      CALL derivk3(th,gn,2)
      CALL fftp3d_complex_to_real(plancr,gn,r1,MPI_COMM_WORLD)
      CALL rotor3(vx,vz,c1,2)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = r3(i,j,k) + r1(i,j,k)*r2(i,j,k)*tmp
            END DO
         END DO
      END DO

      CALL derivk3(th,gn,3)
      CALL fftp3d_complex_to_real(plancr,gn,r1,MPI_COMM_WORLD)
      CALL rotor3(vx,vy,c1,3)
      CALL fftp3d_complex_to_real(plancr,c1,r2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
      DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
            DO i = 1,n
               r3(i,j,k) = r3(i,j,k) + r1(i,j,k)*r2(i,j,k)*tmp
            END DO
         END DO
      END DO

       CALL fftp3d_real_to_complex(planrc,r3,gn,MPI_COMM_WORLD)

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wgradt


!*****************************************************************
      SUBROUTINE wvzspectrum(a0,vx,vy,vz,th,omega,bvfreq,dir,nmb,nmb1,i2d,c1,c2,c3,r1,r2,r3)
!-----------------------------------------------------------------
!
! Computes the spectra for the enstropy constributions with
! isotropic, parallel, and perpendicular dependencies. 
! The output is written to a file by the first node.
!
! Parameters
!     a0    : input matrix: vortical modes
!     vx-vz : velocity components, complex
!     th    : potl temp., complex
!     bvfreq: Brunt-Vaisalla freq
!     omega : rotation rate
!     dir   : output directory (for 2d spectra)
!     nmb   : the time stamp used when writing the file
!     nmb1  : if lenght>0, used to specify range
!     i2d   : do 2D spectra (>0) or not (<=0)
!     c1-3  : complex tmp array
!     r1-3  : real tmp array
!-----------------------------------------------------------------
      USE filefmt
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: E2k  ,E3k  ,E4k
      DOUBLE PRECISION, DIMENSION(n/2+1) :: E2kpa,E3kpa,E4kpa
      DOUBLE PRECISION, DIMENSION(n/2+1) :: E2kpr,E3kpr,E4kpr
      REAL   (KIND=GP),                DIMENSION(n/2+1,n/2+1)   :: eaxi
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0,vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1,c2,c3
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: r1,r2,r3
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      INTEGER,          INTENT   (IN)                           :: i2d
      INTEGER                      :: j
      CHARACTER(len=*), INTENT(IN) :: dir,nmb,nmb1

!
! isotropic spectra:
      CALL PVn(c2,a0,vx,vy,vz,th,omega,bvfreq,1,c1,r1,r2,r3) 
      CALL spectrumg(c2,1,E2k)   ! isotropic spectrum of Z2
      CALL spectrumg(c2,2,E2kpr) ! reduced perp spectrum of Z2
      CALL spectrumg(c2,3,E2kpa) ! reduced para spectrum of Z2
      IF ( i2d.GT.0 ) THEN
      CALL specaxig (c2,   eaxi) ! axisymm. 2d spectra for Z2
      IF (myrank.eq.0) THEN
        if ( len_trim(nmb1).gt.0 ) then
        OPEN(1,file=trim(dir) // '/' // 'kspec2DZ2.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else
        OPEN(1,file=trim(dir) // '/' // 'kspec2DZ2.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) eaxi
        CLOSE(1)
      ENDIF
      ENDIF

      CALL wgradt(c3,vx,vy,vz,th,c1,r1,r2,r3)
      CALL spectrumg(c3,1,E4k)   ! isotropic spectrum of Z4
      CALL spectrumg(c3,2,E4kpr) ! reduced perp spectrum of Z4
      CALL spectrumg(c3,3,E4kpa) ! reduced para spectrum of Z4
      IF ( i2d.GT.0 ) THEN
      CALL specaxig (c3,   eaxi) ! axisymm. 2d spectra for Z4
      IF (myrank.eq.0) THEN
        if ( len_trim(nmb1).gt.0 ) then
        OPEN(1,file=trim(dir) // '/' // 'kspec2DZ4.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else
        OPEN(1,file=trim(dir) // '/' // 'kspec2DZ4.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) eaxi
        CLOSE(1)
      ENDIF
      ENDIF

!      CALL PVn(c3,a0,vx,vy,vz,th,omega,bvfreq,3,c1,c2,r1,r2,r3)
if(myrank.eq.0) write(*,*)'wvzspectrum: calling specprodg iso....'
      CALL specprodg(c2,c3,1,E3k)   ! isotropic spectrum of Z3
if(myrank.eq.0) write(*,*)'wvzspectrum: calling specprodg perp....'
      CALL specprodg(c2,c3,2,E3kpr) ! reduced perp spectrum of Z3
if(myrank.eq.0) write(*,*)'wvzspectrum: calling specprodg para....'
      CALL specprodg(c2,c3,3,E3kpa) ! reduced para spectrum of Z3
if(myrank.eq.0) write(*,*)'wvzspectrum: specprodg para done.'
      IF ( i2d.GT.0 ) THEN
if(myrank.eq.0) write(*,*)'wvzspectrum: calling specprodaxig iso....'
      CALL specprodaxig (c2,c3,eaxi) ! axisymm. 2d spectra for Z3
if(myrank.eq.0) write(*,*)'wvzspectrum: specprodaxig iso done.'
      IF (myrank.eq.0) THEN
        if ( len_trim(nmb1).gt.0 ) then
        OPEN(1,file=trim(dir) // '/' // 'kspec2DZ3.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else
        OPEN(1,file=trim(dir) // '/' // 'kspec2DZ3.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) eaxi
        CLOSE(1)
      ENDIF
      ENDIF


if(myrank.eq.0) write(*,*)'wvzspectrum: writing z iso spectra....'
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvzkspectrum.' // nmb // '_' // trim(nmb1) // '.txt' )
         else
         OPEN(1,file='wvzkspectrum.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(3(E23.15,1X))') E2k(j),E3k(j),E4k(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
! perp spectra:
if(myrank.eq.0) write(*,*)'wvzspectrum: writing z perp spectra....'
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvzkspecperp.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvzkspecperp.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(3(E23.15,1X))') E2kpr(j),E3kpr(j),E4kpr(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
! para spectra:
if(myrank.eq.0) write(*,*)'wvzspectrum: writing z para spectra....'
      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvzkspecpara.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvzkspecpara.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(3(E23.15,1X))') E2kpa(j),E3kpa(j),E4kpa(j)
         ENDDO
         CLOSE(1)
      ENDIF
if(myrank.eq.0) write(*,*)'wvzspectrum: done.'
!
      CALL mpi_barrier(MPI_COMM_WORLD,ierr)

      RETURN
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE wvzspectrum


!*****************************************************************
      SUBROUTINE spectrumg(a,kgeo,F0k)
!-----------------------------------------------------------------
!
! Computes 1d iso or aniso spectrum of input quantity, and
! outputs to F0k on all tasks.
!
! Parameters
!     a    : input quantity, complex
!     kgeo : = 1, then do isotropic spectral dependence; 2, then do perpendicular 
!              spectral dependence; =3, then do parallel spectral dependence.
!     F0k  : spectrum, returned
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)         :: F0k
      DOUBLE PRECISION,  DIMENSION(n/2+1)                     :: E0k
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a
      INTEGER, INTENT(IN)          :: kgeo
      DOUBLE PRECISION             :: ks,tmr,tmp
      INTEGER                      :: i,ibeg,j,k
      INTEGER                      :: kmn,kiso,kperp,kpara,kmsk(3)

      kmsk(1:3) = 0
      IF ( kgeo.lt. 1 .OR. kgeo.gt.3 ) THEN
        WRITE(*,*)'wvspectrumc: geometry parameter invalid.'
        STOP
      ENDIF
      kmsk(kgeo) = 1
!
      E0k     = 0.0D0
      F0k     = 0.0D0
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = ista
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,tmr)
          DO j = 1,n
             DO k = 1,n
                kiso  = int(sqrt(ka2(k,j,1))+.501)
                kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                kpara = int(abs(ka(k))+1)
                kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                tmr     =      abs(a(k,j,1))**2  * tmp
!$omp critical
                E0k(kmn) = E0k(kmn)+tmr
!$omp end critical
             END DO
          END DO
        ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmr)
        DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,tmr)
           DO j = 1,n
              DO k = 1,n
                kiso  = int(sqrt(ka2(k,j,i))+.501)
                kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                kpara = int(abs(ka(k))+1)
                kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                tmr     = 2.0D0* abs(a(k,j,i))**2 * tmp
!$omp critical
                E0k(kmn) = E0k(kmn)+tmr
!$omp end critical

             END DO
           END DO
        END DO
!
        CALL MPI_REDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)

        RETURN

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE spectrumg

!*****************************************************************
      SUBROUTINE specprodg(a,b,kgeo,F0k)
!-----------------------------------------------------------------
!
! Computes 1d iso or aniso spectrum of input quantity, and
! outputs to F0k on all tasks.
!
! Parameters
!     a    : input quantity, complex
!     b    : input quantity, complex
!     kgeo : = 1, then do isotropic spectral dependence; 2, then do perpendicular 
!              spectral dependence; =3, then do parallel spectral dependence.
!     F0k  : spectrum, returned
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT), DIMENSION(n/2+1)         :: F0k
      DOUBLE PRECISION,  DIMENSION(n/2+1)                     :: E0k
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a,b
      INTEGER, INTENT(IN)          :: kgeo
      DOUBLE PRECISION             :: ks,tmr,tmp
      INTEGER                      :: i,ibeg,j,k
      INTEGER                      :: kmn,kiso,kperp,kpara,kmsk(3)

      kmsk(1:3) = 0
      IF ( kgeo.lt. 1 .OR. kgeo.gt.3 ) THEN
        WRITE(*,*)'wvspectrumc: geometry parameter invalid.'
        STOP
      ENDIF
      kmsk(kgeo) = 1
!
      E0k     = 0.0D0
      F0k     = 0.0D0
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = ista
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,tmr)
          DO j = 1,n
             DO k = 1,n
                kiso  = int(sqrt(ka2(k,j,1))+.501)
                kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                kpara = int(abs(ka(k))+1)
                kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                tmr     =      real(a(k,j,1)*conjg(b(k,j,1)))*tmp
!$omp critical
                E0k(kmn) = E0k(kmn)+tmr
!$omp end critical
             END DO
          END DO
        ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,tmr)
        DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,tmr)
           DO j = 1,n
              DO k = 1,n
                kiso  = int(sqrt(ka2(k,j,i))+.501)
                kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                kpara = int(abs(ka(k))+1)
                kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                IF ( (kmn.lt.1).or.(kmn.gt.n/2+1) ) CYCLE
                tmr     = 2.0D0* real(a(k,j,i)*conjg(b(k,j,i)))*tmp
!$omp critical
                E0k(kmn) = E0k(kmn)+tmr
!$omp end critical

             END DO
           END DO
        END DO
!
        CALL MPI_REDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
             MPI_SUM,0,MPI_COMM_WORLD,ierr)

        RETURN

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE specprodg


!*****************************************************************
      SUBROUTINE specaxig(a,F0k)
!-----------------------------------------------------------------
!
! Computes 2d axisymmetric spectrum of input quantity, and
! outputs to F0k on all tasks.
!
! Parameters
!     a    : input quantity, complex
!     F0k  : spectrum, returned
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(OUT), DIMENSION(n/2+1,n/2+1)   :: F0k
      REAL   (KIND=GP),              DIMENSION(n/2+1,n/2+1)   :: E0k
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a
      REAL   (KIND=GP)             :: ks,tmr,tmp
      INTEGER                      :: i,ibeg,j,k,km
      INTEGER                      :: kperp,kpara

      km      = n/2+1
      E0k     = 0.0
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = ista
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ista.eq.1) THEN
!$omp parallel do private (k,kperp,kpara,tmr)
          DO j = 1,n
             kperp = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
             IF ( (kperp.lt.1).or.(kperp.gt.km) ) CYCLE
             DO k = 1,n
                kpara = int(abs(ka(k))+1)
                IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                tmr     =      abs(a(k,j,1))**2  * tmp
!$omp critical
                E0k(kperp,kpara) = E0k(kperp,kpara)+tmr
!$omp end critical
             END DO
          END DO
        ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kperp,kpara,tmr)
        DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kperp,kpara,tmr)
           DO j = 1,n
              kperp = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
              IF ( (kperp.lt.1).or.(kperp.gt.km) ) CYCLE
              DO k = 1,n
                kpara = int(abs(ka(k))+1)
                IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                tmr   = 2.0*abs(a(k,j,i))**2 * tmp
!$omp critical
                E0k(kperp,kpara) = E0k(kperp,kpara)+tmr
!$omp end critical

             END DO
           END DO
        END DO
!
        CALL MPI_REDUCE(E0k,F0k,(n/2+1)*(n/2+1),GC_REAL,&    
             MPI_SUM,0,MPI_COMM_WORLD,ierr)


        RETURN

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE specaxig


!*****************************************************************
      SUBROUTINE specprodaxig(a,b,F0k)
!-----------------------------------------------------------------
!
! Computes 2d axisymmetric spectrum of input quantity, and
! outputs to F0k on all tasks.
!
! Parameters
!     a    : input quantity, complex
!     b    : input quantity, complex
!     F0k  : spectrum, returned
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE ali
!$    USE threads
      IMPLICIT NONE

      REAL   (KIND=GP), INTENT(OUT), DIMENSION(n/2+1,n/2+1)   :: F0k
      REAL   (KIND=GP),              DIMENSION(n/2+1,n/2+1)   :: E0k
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a,b
      REAL   (KIND=GP)             :: ks,tmr,tmp
      INTEGER                      :: i,ibeg,j,k,km
      INTEGER                      :: kperp,kpara

      km      = n/2+1
      E0k     = 0.0
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = ista
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       IF (ista.eq.1) THEN
!$omp parallel do private (k,kperp,kpara,tmr)
          DO j = 1,n
             kperp = int(sqrt(ka(1)**2+ka(j)**2)+1.501)
             IF ( (kperp.lt.1).or.(kperp.gt.km) ) CYCLE
             DO k = 1,n
                kpara = int(abs(ka(k))+1)
                IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                tmr     =      real(a(k,j,1)*conjg(b(k,j,1)))*tmp
!$omp critical
                E0k(kperp,kpara) = E0k(kperp,kpara)+tmr
!$omp end critical
             END DO
          END DO
        ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kperp,kpara,tmr)
        DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kperp,kpara,tmr)
           DO j = 1,n
              kperp = int(sqrt(ka(i)**2+ka(j)**2)+1.501)
              IF ( (kperp.lt.1).or.(kperp.gt.km) ) CYCLE
              DO k = 1,n
                kpara = int(abs(ka(k))+1)
                IF ( (kpara.lt.1).or.(kpara.gt.km ) ) CYCLE
                tmr   = 2.0*real(a(k,j,i)*conjg(b(k,j,i)))*tmp
!$omp critical
                E0k(kperp,kpara) = E0k(kperp,kpara)+tmr
!$omp end critical

             END DO
           END DO
        END DO
!
        CALL MPI_REDUCE(E0k,F0k,(n/2+1)*(n/2+1),GC_REAL,&    
             MPI_SUM,0,MPI_COMM_WORLD,ierr)


        RETURN

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE specprodaxig


!*****************************************************************
      SUBROUTINE PutNormModes(iunit,odir,fout,planio,c1)
!-----------------------------------------------------------------
!
! Writes a complex field to disk.
!
! Parameters
!   iunit   : logical I/O unit
!   odir    : output directory
!   fout    : file name (minus directory name, that is)
!   planio  : IO plan -- for complex I/O
!   c1      : complex array to output
!
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE iovar
      USE iompi
      USE ali
      USE gutils
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1
      INTEGER         , INTENT   (IN)                           :: iunit
      INTEGER                                                   :: bmanghold,i,j,k
      CHARACTER*(*)   , INTENT   (IN)                           :: fout,odir
      TYPE(IOPLAN)    , INTENT(INOUT)                           :: planio


     
      bmanghold = bmangle
      IF ( oswap.NE.0 ) THEN
        CALL carray_byte_swap(c1, n*n*(iend-ista+1))
      ENDIF

      bmangle = 0
      CALL io_writec(iunit,odir,trim(fout),'',planio,c1)
      bmangle = bmanghold

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE PutNormModes


!*****************************************************************
      SUBROUTINE PutRealFields(iunit,odir,ext,planio,a0,am,ap,vx,vy,vz,th, &
                               omega,bvfreq,c1,c2,c3,c4,r1,r2,ivec,suff)
!-----------------------------------------------------------------
!
! Writes a complex field to disk.
!
! Parameters
!   iunit   : logical I/O unit
!   odir    : output directory
!   ext     : time index (char form)
!   planio  : IO plan -- for complex I/O
!   a0,m,p  : complex norm. mode coeffs, input
!   vx,y,z,th 
!   omega   : rotation rate
!   bvfreq  : Brunt-Vaisalla freq.
!   c1-4    : complex tmp arrays
!   r1-2    : real tmp arrays.,
!   ivec: > 0 ==> output vector mag only (th still printed, though)
!   ssuff   : file prefix's suffix (may be null: '')
!
!-----------------------------------------------------------------
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
      USE iompi
      USE iovar
      USE gutils
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0,am,ap
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1,c2,c3,c4
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL   (KIND=GP), INTENT   (IN)                           :: omega,bvfreq
      REAL   (KIND=GP)                                          :: tmp
      INTEGER         , INTENT   (IN)                           :: iunit,ivec
      INTEGER                                                   :: i,j,k,m
      CHARACTER*(*)   , INTENT   (IN)                           :: odir,ext,suff
      TYPE(IOPLAN)    , INTENT(INOUT)                           :: planio

      DO m = 0, 2 ! 0: vortical; 1: +modes; 2: -modes
        CALL DoExpProj(vx,vy,vz,th,a0,ap,am,omega,bvfreq,m)
        CALL printprim(vx,vy,vz,th,odir,ext,planio,c1,r1,r2,m,ivec,suff)
        CALL printdprim(vx,vy,vz,th,odir,ext,planio,c1,c2,c3,c4,r1,r2,m,ivec,suff)
      ENDDO

        ! Compute and output total fields:
      vx = 0.0;
      vy = 0.0;
      vz = 0.0;
      th = 0.0;
      DO m = 0, 2
!if ( myrank.eq.0 ) write(*,*)'main: Total DoExpProj, j=',m, ' time=',ext
        CALL DoExpProj(c1,c2,c3,c4,a0,ap,am,omega,bvfreq,m)
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
          DO j = 1,n
            DO k = 1,n
              vx(k,j,i) = vx(k,j,i) + c1(k,j,i)
              vy(k,j,i) = vy(k,j,i) + c2(k,j,i)
              vz(k,j,i) = vz(k,j,i) + c3(k,j,i)
              th(k,j,i) = th(k,j,i) + c4(k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDDO ! m-loop
      CALL printprim (vx,vy,vz,th,odir,ext,planio,c1,r1,r2,3,ivec,suff)
      CALL printdprim(vx,vy,vz,th,odir,ext,planio,c1,c2,c3,c4,r1,r2,3,ivec,suff)

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE PutRealFields


!*****************************************************************
      SUBROUTINE DoExpProj(vx,vy,vz,th,a0,ap,am,omega,bvfreq,kout)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Based on the coefficients from the normal mode decomposition,
! computes the slow mode, + and - contribution to the velocity and theta
! projections ('expanded' projections).
! Ref: Herbert, et al. JFM 758:374 (2014) Eqs A11, A15
!
!
! Parameters
!     vx,
!     vy,
!     vz,
!     th    : complex output fields (either vel. or vorticity), returned;
!     a0,p,m: complex array of size a,b,c, containing normal modes, input
!     omega : rotation rate
!     bvfreq: Brunt-Vaisalla frequency
!     kout  : 0: vortical modes; 1: '+'-modes; 2: '-' modes
!
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(n,n,ista:iend) :: a0,am,ap
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP)                                          :: ic
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      INTEGER         , INTENT   (IN)                           :: kout
      REAL   (KIND=GP)                                          :: f,kp,ks,sig,sq2,sq2i,tmp
      INTEGER                                                   :: i,j,k


      ic  = cmplx(0.0_GP,1.0_GP);
      f   = 2.0_GP * omega
      sq2i= 1.0_GP/sqrt(2.0_GP)     
      sq2 = sqrt(2.0_GP)     
!if ( myrank.eq.0 ) write(*,*)'bvfreq=',bvfreq,' omega=',omega,' f=',f,' tiny=',tiny
      IF ( kout .eq. 0 ) THEN  ! vortical modes
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig,tmp)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))
     
                  IF ( kp.lt.tiny) THEN
                    ! From decomp from Herbert, eq. A16:
                    vx(k,j,i) = 0.0
                    vy(k,j,i) = 0.0
                    vz(k,j,i) = 0.0
                    th(k,j,i) = -a0(k,j,i)
                  ELSE
                    ! From decomp from Herbert, eq. A14:
                    ! A^0
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    tmp = 1.0/(ks*kp*sig)
                    vx(k,j,i) = -tmp*bvfreq*ka(j)*kp * a0(k,j,i)
                    vy(k,j,i) =  tmp*bvfreq*ka(i)*kp * a0(k,j,i)
                    vz(k,j,i) =  0.0
                    th(k,j,i) = -tmp*f*ka(k)*kp*a0(k,j,i)
                 ENDIF

               END DO
            END DO
         END DO
      ENDIF

      IF ( kout .EQ. 1 ) THEN   ! '+' waves

!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig,tmp)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))

                  IF ( kp.lt.tiny) THEN
                    vx(k,j,i) = -sq2i*ic*ap(k,j,i)
                    vy(k,j,i) =  sq2i*   ap(k,j,i)
                    vz(k,j,i) = 0.0
                    th(k,j,i) = 0.0
                  ELSE

                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    tmp = 1.0/(sq2*ks*kp*sig)
                    vx(k,j,i) = tmp*( f*ka(j)*ka(k)+ic*ka(i)*ka(k)*sig) * ap(k,j,i)
                    vy(k,j,i) = tmp*(-f*ka(i)*ka(k)+ic*ka(j)*ka(k)*sig) * ap(k,j,i)
                    vz(k,j,i) = -tmp*ic*(kp**2)*sig*ap(k,j,i)
                    th(k,j,i) = -tmp*bvfreq*(kp**2)*ap(k,j,i)
                 ENDIF

               END DO
            END DO
         END DO
      ENDIF

                  
      IF ( kout .EQ. 2 ) THEN  ! '-' waves

!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig,tmp)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig,tmp)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(j)**2)
                  ks  = sqrt(ka2(k,j,i))

                  IF ( kp.lt.tiny) THEN
                    vx(k,j,i) =  sq2i*ic*am(k,j,i)
                    vy(k,j,i) =  sq2i*   am(k,j,i)
                    vz(k,j,i) = 0.0
                    th(k,j,i) = 0.0
                  ELSE
                  
                    sig = sqrt((f**2*ka(k)**2+bvfreq**2*kp**2)/ka2(k,j,i))
                    tmp = 1.0/(sq2*ks*kp*sig)
                    vx(k,j,i) = tmp*( f*ka(j)*ka(k)-ic*ka(i)*ka(k)*sig) * am(k,j,i)
                    vy(k,j,i) = tmp*(-f*ka(i)*ka(k)-ic*ka(j)*ka(k)*sig) * am(k,j,i)
                    vz(k,j,i) =  tmp*ic*(kp**2)*sig*am(k,j,i)
                    th(k,j,i) = -tmp*bvfreq*(kp**2)*am(k,j,i)
                 ENDIF

               END DO
            END DO
         END DO
      ENDIF

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE DoExpProj


!*****************************************************************
      SUBROUTINE printprim(vx,vy,vz,th,odir,ext,planio,c1,r1,r2,kout,ivec,prend)
!-----------------------------------------------------------------
!   vx, vy, vz, th: primitive fields
!   odir    : output directory
!   ext     : string time stamp
!   planio  : IO plan
!   c1      : complex tmp array
!   r1-2    : real tmp arrays
!   kout    :  0 == vort. modes, 1 = + modes, 2 = - modes 3 = entire
!   field
!   ivec    : ==0, don't print vectors, ==1; print velocity comp only; 2:
!             print magnitude only; ==3, print velocity comps and mag
!   prend   : file prefix's suffix (may be null: '')
!-----------------------------------------------------------------
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
      USE iompi
      USE iovar
      USE gutils
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL   (KIND=GP)                                          :: tmp,rmin,rmax,rmean,rrms,rstd
      INTEGER         , INTENT   (IN)                           :: kout,ivec
      INTEGER                                                   :: i,j,k
      CHARACTER*(*)   , INTENT   (IN)                           :: odir,ext,prend
      CHARACTER*1                                               :: suff(4)
      CHARACTER*128                                             :: sout
      TYPE(IOPLAN)                                              :: planio
      
      suff(1) = '0' ! 0-modes
      suff(2) = 'p' ! +-modes
      suff(3) = 'm' ! --modes
      suff(4) = 't' ! total field (sum of all modes)
      IF ( kout.lt.0 .or. kout.gt.3 ) THEN
        WRITE(*,*)'printprim: kout=',kout,' must be 0 <= kout < 4'
        STOP
      ENDIF     
 
      tmp = 1.0_GP/real(n,kind=GP)**3

      c1 = th*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'th'   // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)
   

      IF ( ibits(ivec,1,1).eq.1 ) THEN
        CALL VecMagC2R(r1,vx,vy,vz,c1,r2)
        CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
        IF ( oswap.NE.0 ) THEN
          CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
        ENDIF
        sout = 'vmag'  // trim(prend) // suff(kout+1) 
        CALL io_write(1,odir,trim(sout),ext,planio,r1)
        CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)
      ENDIF

      IF ( ibits(ivec,0,1) .eq. 0 ) RETURN
 
      c1 = vx*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'vx'   // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)

      c1 = vy*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'vy'  // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)

      c1 = vz*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'vz' // trim(prend) //  suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)


!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE printprim


!*****************************************************************
      SUBROUTINE printdprim(vx,vy,vz,th,odir,ext,planio,c1,c2,c3,c4,r1,r2,kout,ivec,prend)
!-----------------------------------------------------------------
!   vx, vy, vz, th: primitive fields
!   odir    : output directory
!   ext     : string time stamp
!   planio  : IO plan
!   c1-4    : complex tmp arrays
!   r1-2    : real tmp arrays
!   kout    :  0 == vort. modes, 1 = + modes, 2 = - modes, 3 = entire
!   field
!   ivec: > 0 ==> output vector mag only (th still printed, though)
!   prend   : file prefix's suffix (may be null: '')
!-----------------------------------------------------------------
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
      USE iompi
      USE iovar
      USE gutils
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz,th
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1,c2,c3,c4
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1,r2
      REAL   (KIND=GP)                                          :: tmp,rmin,rmax,rmean,rrms,rstd
      INTEGER         , INTENT   (IN)                           :: kout,ivec
      INTEGER                                                   :: i,j,k
      CHARACTER*(*)   , INTENT   (IN)                           :: odir,ext,prend
      CHARACTER*1                                               :: suff(4)
      CHARACTER*128                                             :: sout
      TYPE(IOPLAN)    , INTENT(INOUT)                           :: planio

      
      suff(1) = '0' ! 0-modes
      suff(2) = 'p' ! +-modes
      suff(3) = 'm' ! --modes
      suff(4) = 't' ! total field (sum of all modes)

      tmp = 1.0_GP/real(n,kind=GP)**3

      c1 = th*tmp
      CALL derivk3(th,c1,1)
      CALL derivk3(th,c2,2)
      CALL derivk3(th,c3,3)
      CALL VecMagC2R(r1,c1,c2,c3,c4,r2)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'dth' // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
   
      IF ( ibits(ivec,0,1).eq.1 ) THEN
        CALL rotor3(vy,vz,c1,1)
        CALL rotor3(vx,vz,c2,2)
        CALL rotor3(vx,vy,c3,3)
      ENDIF

      IF ( ibits(ivec,1,1).eq.1 ) THEN
        CALL VecMagC2R(r1,c1,c2,c3,c4,r2)
        CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)

        IF ( oswap.NE.0 ) THEN
          CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
        ENDIF
        sout = 'wmag' // trim(prend) // suff(kout+1)
        CALL io_write(1,odir,trim(sout),ext,planio,r1)
        CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)
      ENDIF

      IF ( ibits(ivec,0,1).eq.0 ) RETURN

      c1 = c1*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'wx'  // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)

      c1 = c2*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'wy'  // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)

      c1 = c3*tmp
      CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL rarray_props(r1,n*n*(kend-ksta+1),n,rmin,rmax,rmean,rrms,rstd)
      IF ( oswap.NE.0 ) THEN
        CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
      ENDIF
      sout = 'wz'  // trim(prend) // suff(kout+1)
      CALL io_write(1,odir,trim(sout),ext,planio,r1)
      CALL printprops('primprops',trim(sout),ext,rmin,rmax,rmean,rrms,rstd)

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE printdprim



!*****************************************************************
      SUBROUTINE VecMagC2R(rmag,vx,vy,vz,c1,r1)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
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

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: vx,vy,vz
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1
      REAL   (KIND=GP), INTENT  (OUT), DIMENSION(n,n,ksta:kend) :: rmag
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1
      REAL   (KIND=GP)                                          :: tmp
      INTEGER                                                   :: i,j,k

      tmp = 1.0_GP/real(n,kind=GP)**3

       c1 = vx
       CALL fftp3d_complex_to_real(plancr,c1,rmag,MPI_COMM_WORLD)
       c1 = vy
       CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
       DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
           DO i = 1,n
             rmag(i,j,k) = rmag(i,j,k)**2 + r1(i,j,k)**2
           ENDDO
         ENDDO
       ENDDO
       c1 = vz
       CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
       DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
         DO j = 1,n
           DO i = 1,n
             rmag(i,j,k) = sqrt(rmag(i,j,k) + (r1(i,j,k)**2)) * tmp
           ENDDO
         ENDDO
       ENDDO

!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE VecMagC2R


!*****************************************************************
      SUBROUTINE printprops(fname,vname,ext,rmin,rmax,rmean,rrms,rstd)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
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

      REAL   (KIND=GP), INTENT   (IN)                           :: rmin,rmax,rmean,rrms,rstd
      INTEGER                                                   :: i,j,k
      CHARACTER*(*)   , INTENT   (IN)                           :: fname,vname,ext

      IF ( myrank.NE.0 ) RETURN

      OPEN(1,file=fname,position='append')
      WRITE(1,*) vname,ext,rmin,rmax,rmean,rrms,rstd
      CLOSE(1)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
      END SUBROUTINE printprops
