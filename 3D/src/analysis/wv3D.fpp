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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz,th,a0,am,cv


      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: rv
      REAL(KIND=GP)                                 :: omega,bvfreq
      REAL(KIND=GP)                                 :: tmp
!
! Auxiliary variables

      INTEGER :: i,it,iavg,ic,ind,inm,j,k
      INTEGER :: istat(1024),nfiles,nstat

      TYPE(IOPLAN)       :: planio
      CHARACTER(len=8)   :: pref
      CHARACTER(len=256) :: odir,idir
      CHARACTER(len=256) :: fout
      CHARACTER(len=4096):: stat
!
      NAMELIST / wv / idir, odir, stat, iswap, oswap, iavg, inm, omega, bvfreq

      tiny  = 1e-5_GP
      tinyf = 1e-15_GP

!
! Initializes the MPI and I/O libraries
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)
      idir   = '.'
      odir   = '.'
      stat   = '0'
      iswap  = 0
      oswap  = 0
      iavg   = 0 
      inm    = 0 
!
! Reads from the external file 'vt`.txt' the 
! parameters that will be used to compute the transfer
!     idir   : directory for unformatted input (field components)
!     odir   : directory for unformatted output (prolongated data)
!     stat  : time index for which to compute WV, or a ';--separated list
!     iswap  : do endian swap on input?
!     oswap  : do endian swap on output?
!     iavg   : time average spectrally the time index list?
!     inm    : write normal mode fields (0==don't; 1==in real space; 2==in wavenumber space)
!     omega  : rotation rate
!     bvfreq : Brunt-Vaisalla frequency
      IF (myrank.eq.0) THEN
         OPEN(1,file='wv.txt',status='unknown',form="formatted")
         READ(1,NML=wv)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir  ,256 ,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat  ,4096,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oswap ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(iavg  ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(inm   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(omega ,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bvfreq,1   ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)

      pref = 'WV'
!
     
      kmax   = (real(n,kind=GP)/3.)**2

      ALLOCATE( a0(n,n,ista:iend) )
      ALLOCATE( am(n,n,ista:iend) )
      ALLOCATE( vx(n,n,ista:iend) )
      ALLOCATE( vy(n,n,ista:iend) )
      ALLOCATE( vz(n,n,ista:iend) )
      ALLOCATE( th(n,n,ista:iend) )
      ALLOCATE( cv(n,n,ista:iend) )
      ALLOCATE( rv(n,n,ksta:kend) )
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

      tmp = 1.0/REAL(n,KIND=GP)**3
      IF ( iavg.EQ.1 ) THEN
write(*,*)'main: calling DoSpAvg...'
        CALL DoSpAvg(vx,vy,vz,th,istat,nstat,idir,planrc,planio,rv,cv)
        IF ( inm.EQ.1 ) THEN
          bmangle = 0
          CALL fftp3d_complex_to_real(plancr,a0,rv,MPI_COMM_WORLD)
          IF ( oswap.NE.0 ) THEN
            CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
          WRITE(fout,'(a,i4.4,a,i4.4,a)'),'a0av_',istat(1),'_',istat(nstat),'.out'
          CALL io_write(1,odir,trim(fout),ext,planio,rv)
          CALL fftp3d_complex_to_real(plancr,am,rv,MPI_COMM_WORLD)
          IF ( oswap.NE.0 ) THEN
            CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
          WRITE(fout,'(a,i4.4,a,i4.4,a)'),'amav_',istat(1),'_',istat(nstat),'.out'
          CALL io_write(1,odir,fout,ext,planio,rv)
          bmangle = 1
        ENDIF
        WRITE(ext, fmtext) istat(nstat)
        CALL WVNormal(a0,am,vx,vy,vz,th,omega,bvfreq)
        CALL wvspectrum(a0,am,bvfreq,omega,ext,1)
      ELSE
        DO it = 1,nstat
          WRITE(ext, fmtext) istat(it)
! read in appropriate file:
write(*,*)'main: reading vx...'
          CALL io_read(1,idir,'vx',ext,planio,rv)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on vx...'
          CALL fftp3d_real_to_complex(planrc,rv,vx,MPI_COMM_WORLD)
write(*,*)'main: reading vy...'
          CALL io_read(1,idir,'vy',ext,planio,rv)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on vy...'
          CALL fftp3d_real_to_complex(planrc,rv,vy,MPI_COMM_WORLD)
write(*,*)'main: reading vz...'
          CALL io_read(1,idir,'vz',ext,planio,rv)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on vz...'
          CALL fftp3d_real_to_complex(planrc,rv,vz,MPI_COMM_WORLD)
write(*,*)'main: reading th...'
          CALL io_read(1,idir,'th',ext,planio,rv)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on th...'
          CALL fftp3d_real_to_complex(planrc,rv,th,MPI_COMM_WORLD)
write(*,*)'main: calling WVNormal ...'
          CALL WVNormal(a0,am,vx,vy,vz,th,omega,bvfreq)
          CALL wvspectrum(a0,am,bvfreq,omega,ext,1)

          ! Do output:
          IF ( inm.EQ.1 ) THEN
write(*,*)'main: fft_cr on a0...'
            CALL fftp3d_complex_to_real(plancr,a0,rv,MPI_COMM_WORLD)
            IF ( oswap.NE.0 ) THEN
              CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
            ENDIF
write(*,*)'main: writing a0...'
            CALL io_write(1,odir,'a0',ext,planio,rv)
write(*,*)'main: fft_cr on a0...'
            CALL fftp3d_complex_to_real(plancr,am,rv,MPI_COMM_WORLD)
            IF ( oswap.NE.0 ) THEN
              CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
            ENDIF
write(*,*)'main: writing am...'
            CALL io_write(1,odir,'am',ext,planio,rv)
          ENDIF
        ENDDO
      ENDIF
!
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)


      DEALLOCATE (a0,am)
      DEALLOCATE (vx,vy,vz,th)
      DEALLOCATE (rv,cv)
      DEALLOCATE (ka)
      DEALLOCATE (ka2)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM WV3D


      SUBROUTINE WVNormal(a0,am,vx,vy,vz,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Carries out wave-vortical decomposition tasks.
! Computes from Bartello: J. Atm. Sci.  ! 52(24):4410 (1995) 
! eq 6a-b, except we divide a0, a+- there by kperp.
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
      COMPLEX(KIND=GP), INTENT  (OUT), DIMENSION(n,n,ista:iend) :: a0,am
      COMPLEX(KIND=GP)                                          :: ic
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,ks,sig,tmp
      INTEGER                                                   :: i,j,k

      a0 = 0.0_GP
      am = 0.0_GP

write(*,*)'WVNormal 0...'
      ic = cmplx(0.0_GP,1.0_GP);
      f = 2.0_GP * omega
      tmp = 1.0_GP/sqrt(2.0_GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kp,ks,sig)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kp,ks,sig)
            DO j = 1,n
               DO k = 1,n
                  kp  = sqrt(ka(i)**2+ka(k)**2)
                  ks  = sqrt(ka2(k,j,i))
                  sig = sqrt(f**2*ka(k)**2+bvfreq**2*kp**2)
     
                    IF ( kp.lt.tiny) THEN
                      a0(k,j,i) = -th(k,j,i)
                      am(k,j,i) = tmp*( vy(k,j,i) - ic*vx(k,j,i) )
                    ELSE
                      ! if A is the decomp from Bartello, eq. 6a-b, then
                      ! A^0/kperp:
                      a0(k,j,i) = ( -bvfreq*( ka(j)*vx(k,j,i) - ka(i)*vy(k,j,i) ) &
                                  - f*ka(k)*th(k,j,i) ) / sig
                      ! A^-/kperp:
                      
                      am(k,j,i) = ( f*ka(k)*( ka(j)*vx(k,j,i) - ka(i)*vy(k,j,i) ) &
                                - ka2(k,j,i)*ic*sig*vz(k,j,i) - bvfreq*kp**2*th(k,j,i) ) / (kp*sig)
                  
                      ! In this basis, A^+ = conjg(A^-)
                  ENDIF

               END DO
            END DO
         END DO

      END SUBROUTINE WVNormal


!*****************************************************************
      SUBROUTINE wvspectrum(a0,am,bvfreq,omega,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy and helicity power 
! spectrum. The output is written to a 
! file by the first node.
!
! Parameters
!     a0 : input matrix: vortical modes
!     am : input matrix: wave modes (-)
!     bvfreq: Brunt-Vaisalla freq
!     omega : rotation rate
!     nmb: the extension used when writting the file
!     kin: = 1 (setting bit 0) computes the energy spectra alone
!          = 2 (setting bit 1) computes the helicity spectra alone
!          = 4 (setting bit 2) computes the enstrophy spectra alone
!          = 8 (setting bit 3) computes V and P energy spectra 
!          Setting kin=3 computes E & H; kin=5 computes E &Enst; kin=7
!          computes E & H & Enst; kin=15, computes everything
!
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(n/2+1) :: E0k,EWk,EV0k,EP0k,EVWk,EPWk
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,n,ista:iend) :: a0,am
      REAL   (KIND=GP), INTENT(IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                       :: f
      INTEGER, INTENT(IN)          :: kin
      INTEGER                      :: j
      CHARACTER(len=*), INTENT(IN) :: nmb

!
      f = 2.0*omega
!
! Puts the spectra to files:
      CALL wvspectrumc(a0,am,bvfreq,f,1,1,E0k ,EWk  ) ! total energy
      CALL wvspectrumc(a0,am,bvfreq,f,4,1,EV0k,EVWk ) ! KE (only E0)
      CALL wvspectrumc(a0,am,bvfreq,f,5,1,EP0k,EPWk ) ! PE 
      DO j = 1,n/2+1
        EVWk(j) = E0k(j)+EWk(j)-EV0k(j)-EP0k(j)-EPWk(j)
      ENDDO
      IF (myrank.eq.0) THEN
         OPEN(1,file='wvekspectrum.' // nmb // '.txt')
         DO j=1,n/2+1
           WRITE(1,FMT='(6(E23.15,1X))') E0k(j),EWk(j),EV0k(j),EVWk(j),EP0k(j),EPWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      CALL wvspectrumc(a0,am,bvfreq,f,3,1,E0k ,EWk  ) ! Enstrophy (only E0)
      IF (myrank.eq.0) THEN
         OPEN(1,file='wvenskspectrum.' // nmb // '.txt')
         DO j=1,n/2+1
           WRITE(1,FMT='(E23.15)') E0k(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      CALL wvspectrumc(a0,am,bvfreq,f,2,1,E0k ,EWk  ) ! Helicity 
      IF (myrank.eq.0) THEN
         OPEN(1,file='wvhkspectrum.' // nmb // '.txt')
         DO j=1,n/2+1
           WRITE(1,FMT='(E23.15,1X,E23.15)') E0k(j), EWk(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
      RETURN
      END SUBROUTINE wvspectrum


!*****************************************************************
      SUBROUTINE wvspectrumc(a0,am,bvfreq,f,kin,kgeo,F0k,FWk)
!-----------------------------------------------------------------
!
! Computes 1d wave and vortical spectra for various quantities, returns them
!
! Parameters
!     a0   : input matri: vortical modes
!     am   : input matrix: wave modes (-) 
!     bvfreq: Brunt-Vaisalla frequency
!     f    : rotation parameter = 2\Omega
!     kin  : = 1 computes the total energy spectra wave/vortical
!            = 2 computes the helicity spectra alone wave/vortical 
!            = 3 computes the enstrophy spectra alone wavevortical; only F0k is filled
!            = 4 computes the kinetic energy; vortical only; wave KE must be
!                computed from the outside as EVWk = E0k+EWk-EV0k-EP0k-EPWk,
!                using multiple calls to this routine; only F0k is filled then
!            = 5 computes the potential energy; wave/vortical
!     kgeo : = 1, then do isotropic spectral dependence; 2, then do perpendicular 
!              spectral dependence; =3, then do parallel spectral dependence.
!     F0k  : vortical spectrum, returned
!     FWk  : wave spectrum, returned
!
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
      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(n,n,ista:iend) :: a0,am
      COMPLEX(KIND=GP)                                        :: apc
      REAL   (KIND=GP), INTENT (IN)                           :: bvfreq,f
      INTEGER, INTENT(IN)          :: kin,kgeo
      DOUBLE PRECISION             :: ks,sig,tm0,tmi,tmr,tmw
      REAL(KIND=GP)                :: tmp
      INTEGER                      :: i,ibeg,j,k
      INTEGER                      :: kmn,kiso,kperp,kpara,kp2,kmsk(3)

      kmsk(1:3) = 0
      IF ( kgeo.lt. 1 .OR. kgeo.gt.3 ) THEN
        WRITE(*,*)'wvspectrumc: geometry parameter invalid.'
        STOP
      ENDIF
      kmsk(kgeo) = 1

!
      E0k     = 0.0_GP
      EWk     = 0.0_GP
      tmp     = 1.0_GP/real(n,kind=GP)**6
      ibeg    = 1
      IF (ista.eq.1) ibeg = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.1 ) THEN ! Total energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmr,apc)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ( (kmn.le.tiny).or.(kmn.gt.n/2+1) ) CYCLE
                  apc     = conjg(am(k,j,1))
                  tm0     =     (abs(a0(k,j,1))**2 \
                          +      abs(      apc)**2 ) * tmp
                  tmw     =     (abs(am(k,j,1))**2 \
                          +      abs(      apc)**2 ) * tmp
!$omp critical
                  E0k(kmn) = E0k(kmn)+tm0
                  EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,sig,apc)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,sig,apc)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  IF ((kmn.le.tiny).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,i))
                  tm0     = 2.0* (abs(a0(k,j,i))**2 \
                          +       abs(      apc)**2 ) * tmp
                  tmw     = 2.0 * (abs(am(k,j,1))**2 \
                          +        abs(      apc)**2 ) * tmp
!$omp critical
                  E0k(kmn) = E0k(kmn)+tm0
                  EWk(kmn) = EWk(kmn)+tmw
!$omp end critical

               END DO
             END DO
          END DO
!
          CALL MPI_ALLREDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.2 ) THEN ! Helicity spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.tiny).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,1))
                  IF ((kperp.gt.tiny)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0/(sqrt(2.0)*(1.0+tmr))
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     ks      = sqrt( ka(k)**2+kp2 )
                     tm0     = real( ks*a0(k,j,1) *  \
                             ( conjg(apc)-conjg(am(k,j,1)) ) )*tmp*tmi
                     tmw     = f*ka(k)*ks*( abs(am(k,j,1))**2 -   \
                               abs(apc)**2 )/sig * tmp
                             
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ELSE
                     tmw     = ka(k)*( abs(am(k,j,1))**2-abs(apc)**2 )*tmp
!$omp critical
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmw,tmi,tmr,sig)
             DO j = 1,n
                DO k = 1,n
                  kmn = int(sqrt(ka2(k,j,i))+.501)
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.tiny).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,i))

                  IF ((kp2.gt.tiny)) THEN
                     tmr     = (f*ka(k)/bvfreq)**2/kp2
                     tmi     = 1.0/(sqrt(2.0)*(1.0+tmr))
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     ks      = sqrt ( ka(k)**2+kp2 ) 
                     tm0     = 2.0*real( ks*a0(k,j,i) * \
                             ( conjg(apc)-conjg(am(k,j,i)) ) )*tmp*tmi
                     tmw     = 2.0*f*ka(k)*ks*( abs(am(k,j,i))**2 -   \
                               abs(apc)**2 )/sig * tmp
write(*,*)'wvsectrumc: f=',f,' ks=',ks,' tmp=',tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ELSE
                     tmw     = 2.0*ka(k)*( abs(am(k,j,i))**2-abs(apc)**2 )*tmp
!$omp critical
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO

!         Compute reduction between nodes
!
          CALL MPI_ALLREDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.3 ) THEN ! Enstropy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,sig,apc)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.0).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,1))
                  sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                  tm0     =     sig**2*abs(a0(k,j,1))**2*tmp
!$omp critical
                  E0k(kmn) = E0k(kmn)+tm0
!$omp end critical
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,sig,apc)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,sig,apc)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.0).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,i))
                  sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                  tm0     =     2.0*sig**2*abs(a0(k,j,i))**2*tmp
!$omp critical
                  E0k(kmn) = E0k(kmn)+tm0
!$omp end critical

               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_ALLREDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.4 ) THEN ! Kinetic energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmr,apc,sig)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.0).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,1))
                  IF ((kp2.gt.0)) THEN
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     tmr     = (bvfreq**2)*kp2/ sig**2
                     tm0     = abs(a0(k,j,1))**2*tmr*tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                  ENDIF
               END DO
            END DO
          ENDIF

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,tmr,apc,sig)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmr,apc,sig)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.0).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,i))

                  IF ((kp2.gt.0)) THEN
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     tmr     = (bvfreq**2)*kp2/ sig**2
                     tmr     = (bvfreq**2)*kp2/ sig**2
                     tm0     = 2.0*abs(a0(k,j,i))**2*tmr*tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_ALLREDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF ( kin.eq.5 ) THEN  ! Potential energy spectra
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmi,tmr,apc,sig)
            DO j = 1,n
               DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,1))+.501)
                  kperp = int(sqrt(ka(1)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.0).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,1))

                  IF ((kp2.gt.0)) THEN
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     tmr     = (f*ka(k)/(bvfreq*kperp))**2
                     tmi     = 1.0/(1.0+tmr)
                     tm0     = tmr*abs(a0(k,j,1))**2*tmp*tmi
                     tmw     = 0.5*abs(am(k,j,1)+apc)**2*tmp*tmi
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

!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,kiso,kperp,kpara,kp2,tm0,tmi,tmr,apc,sig)
          DO i = ibeg,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,kiso,kperp,kpara,kp2,tm0,tmi,tmr,apc,sig)
             DO j = 1,n
                DO k = 1,n
                  kiso  = int(sqrt(ka2(k,j,i))+.501)
                  kperp = int(sqrt(ka(i)**2+ka(j)**2)+.501)
                  kpara = int(abs(ka(k))+1)
                  kmn   = kmsk(1)*kiso + kmsk(2)*kperp + kmsk(3)*kpara
                  kp2   = kperp**2
                  IF ((kmn.le.0).or.(kmn.gt.n/2+1)) CYCLE
                  apc     = conjg(am(k,j,i))

                  IF ((kp2.gt.0)) THEN
                     sig     = sqrt(f**2 * ka(k)**2+ bvfreq**2 * kp2)
                     tmr     = (f*ka(k)/(bvfreq*kperp))**2
                     tmi     = 1.0/(1.0+tmr)
                     tm0     = tmr*abs(a0(k,j,i))**2*tmp*tmi
                     tmw     =     abs(am(k,j,i)+apc)**2*tmp*tmi
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
                     EWk(kmn) = EWk(kmn)+tmw
!$omp end critical
                  ELSE
                     tm0     = 2.0*abs(a0(k,j,i))**2*tmp
!$omp critical
                     E0k(kmn) = E0k(kmn)+tm0
!$omp end critical
                 ENDIF
               END DO
             END DO
          END DO
!         Compute reduction between nodes
!
          CALL MPI_ALLREDUCE(E0k,F0k,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)
          CALL MPI_ALLREDUCE(EWk,FWk,n/2+1,MPI_DOUBLE_PRECISION,      &
               MPI_SUM,MPI_COMM_WORLD,ierr)

          RETURN
      ENDIF ! kin = 5
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END SUBROUTINE wvspectrumc


!*****************************************************************
      SUBROUTINE DoSpAvg(vx,vy,vz,th,istat,nstat,idir,planrc,planio,rv,cv)
!-----------------------------------------------------------------
!
! Computes 1d wave and vortical spectra for various quantities, returns them
!
! Parameters
!     vx,vy,vz,th  : complex arrays averaged over istat time stamps
!     istat        : array of time stamp integers
!     nstat        : no. istat timestamps
!     iswap        : do byte swap on input?
!     planrc       : real to complex plan for FFT
!     planio       : IO plan
!     rv           : real tmp array
!     cv           : complex tmp array, size of vx, etc.
!
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
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: cv
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
        CALL fftp3d_real_to_complex(planrc,rv,cv,MPI_COMM_WORLD)
        CALL padd(vx,cv,1.0)

        CALL io_read(1,idir,'vy',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,cv,MPI_COMM_WORLD)
        CALL padd(vy,cv,1.0)

        CALL io_read(1,idir,'vz',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,cv,MPI_COMM_WORLD)
        CALL padd(vz,cv,1.0)

        CALL io_read(1,idir,'th',ext,planio,rv)
        IF ( iswap .NE. 0 ) THEN
           CALL rarray_byte_swap(rv, n*n*(kend-ksta+1))
        ENDIF
        CALL fftp3d_real_to_complex(planrc,rv,cv,MPI_COMM_WORLD)
        CALL padd(th,cv,1.0)
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


      END SUBROUTINE DoSpAvg


!*****************************************************************
      SUBROUTINE padd(v,p,c)
!-----------------------------------------------------------------
!
! Computes v = v + c*p
!
!
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

      END SUBROUTINE padd
