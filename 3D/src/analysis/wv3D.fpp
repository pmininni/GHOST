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

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz,th,a0,am,ap
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: c1,c2,c3
      REAL   (KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: r1,r2,r3
      REAL   (KIND=GP)                                 :: omega,bvfreq
      REAL   (KIND=GP)                                 :: tmp
!
! Auxiliary variables

      INTEGER :: i,i2d,it,iavg,ic,ind,putnm,j,k
      INTEGER :: istat(1024),nfiles,nstat

      TYPE(IOPLAN)       :: planio
      CHARACTER(len=8)   :: pref
      CHARACTER(len=256) :: odir,idir
      CHARACTER(len=256) :: fout
      CHARACTER(len=4096):: stat
      CHARACTER(len=4)   :: ext1
!
      NAMELIST / wv / idir, odir, stat, iswap, oswap, iavg, putnm, omega, bvfreq, i2d

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
      i2d    = 0 
      putnm  = 0 
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
!     putnm  : write normal mode fields (0==don't; 1==in real space; 2==in wavenumber space)
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
      CALL MPI_BCAST(i2d   ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(putnm ,1   ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
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
        IF ( putnm.EQ.1 ) THEN
          bmangle = 0
          c1 = a0*tmp
          CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
          IF ( oswap.NE.0 ) THEN
            CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
          WRITE(fout,'(a,i4.4,a,i4.4,a)'),'a0av_',istat(1),'_',istat(nstat),'.out'
          CALL io_write(1,odir,trim(fout),ext,planio,r1)
          c1 = am*tmp
          CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
          IF ( oswap.NE.0 ) THEN
            CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
          WRITE(fout,'(a,i4.4,a,i4.4,a)'),'amav_',istat(1),'_',istat(nstat),'.out'
          CALL io_write(1,odir,fout,ext,planio,r1)
          c1 = ap*tmp
          CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
          IF ( oswap.NE.0 ) THEN
            CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
          WRITE(fout,'(a,i4.4,a,i4.4,a)'),'apav_',istat(1),'_',istat(nstat),'.out'
          CALL io_write(1,odir,fout,ext,planio,r1)
          bmangle = 1
        ENDIF
        CALL WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq)
        CALL wvspectrum(a0,am,ap,omega,bvfreq,ext,ext1,i2d)
        CALL wvzspectrum(a0,vx,vy,vz,th,omega,bvfreq,ext,ext1,i2d,c1,c2,c3,r1,r2,r3)
      ELSE
        DO it = 1,nstat
          WRITE(ext, fmtext) istat(it)
! read in appropriate file:
write(*,*)'main: reading vx...'
          CALL io_read(1,idir,'vx',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on vx...'
          CALL fftp3d_real_to_complex(planrc,r1,vx,MPI_COMM_WORLD)
write(*,*)'main: reading vy...'
          CALL io_read(1,idir,'vy',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on vy...'
          CALL fftp3d_real_to_complex(planrc,r1,vy,MPI_COMM_WORLD)
write(*,*)'main: reading vz...'
          CALL io_read(1,idir,'vz',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on vz...'
          CALL fftp3d_real_to_complex(planrc,r1,vz,MPI_COMM_WORLD)
write(*,*)'main: reading th...'
          CALL io_read(1,idir,'th',ext,planio,r1)
          IF ( iswap .NE. 0 ) THEN
             CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
          ENDIF
write(*,*)'main: fft_rc on th...'
          CALL fftp3d_real_to_complex(planrc,r1,th,MPI_COMM_WORLD)
write(*,*)'main: calling WVNormal ...'
          CALL WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq)
          CALL wvspectrum(a0,am,ap,omega,bvfreq,ext,'',i2d)
          CALL wvzspectrum(a0,vx,vy,vz,th,omega,bvfreq,ext,'',i2d,c1,c2,c3,r1,r2,r3)

          ! Do output:
          IF ( putnm.EQ.1 ) THEN
            c1 = a0*tmp
            CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
            IF ( oswap.NE.0 ) THEN
              CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
            ENDIF
            CALL io_write(1,odir,'a0',ext,planio,r1)
            c1 = am*tmp
            CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
            IF ( oswap.NE.0 ) THEN
              CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
            ENDIF
            CALL io_write(1,odir,'am',ext,planio,r1)
            c1 = ap*tmp
            CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
            IF ( oswap.NE.0 ) THEN
              CALL rarray_byte_swap(r1, n*n*(kend-ksta+1))
            ENDIF
            CALL io_write(1,odir,'ap',ext,planio,r1)
          ENDIF
        ENDDO
      ENDIF
!
!
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)


      DEALLOCATE (a0,am,ap)
      DEALLOCATE (vx,vy,vz,th)
      DEALLOCATE (c1,c2,c3)
      DEALLOCATE (r1,r2,r3)
      DEALLOCATE (ka)
      DEALLOCATE (ka2)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM WV3D


      SUBROUTINE WVNormal(a0,am,ap,vx,vy,vz,th,omega,bvfreq)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Carries out wave-vortical decomposition tasks.
! Computes from Herbert, et al. JFM 758:374 (2014)
! eq A14-15.
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
write(*,*)'bvfreq=',bvfreq,' omega=',omega,' f=',f,' tiny=',tiny
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

      END SUBROUTINE WVNormal


!*****************************************************************
      SUBROUTINE wvspectrum(a0,am,ap,omega,bvfreq,nmb,nmb1,i2d)
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
      CHARACTER(len=*), INTENT(IN) :: nmb,nmb1

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
        OPEN(1,file='kE02D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kE02D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file='kEW2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kEW2D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) FWaxi
        CLOSE(1)
      ENDIF

      CALL wvspecaxic(a0,am,ap,omega,bvfreq,3,F0axi,FWaxi) ! 2D axisymm spec EV0 (only F0 filled)
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file='kEV02D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kEV02D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      CALL wvspecaxic(a0,am,ap,omega,bvfreq,4,F0axi,FWaxi) ! 2D axisymm spec P0, PW
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file='kP02D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kP02D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file='kPW2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kPW2D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) FWaxi
        CLOSE(1)
      ENDIF

      CALL wvspecaxic(a0,am,ap,omega,bvfreq,2,F0axi,FWaxi) ! 2D axisymm spec H0, HW
      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file='kH02D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kH02D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) F0axi
        CLOSE(1)
      ENDIF

      IF (myrank.eq.0) THEN 
        if ( len_trim(nmb1).gt.0 ) then 
        OPEN(1,file='kHW2D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else 
        OPEN(1,file='kHW2D.' // nmb // '.out',form='unformatted',access='stream')
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END SUBROUTINE wvspecaxic




!*****************************************************************
      SUBROUTINE wvspectrumc(a0,am,ap,omega,bvfreq,kin,kgeo,F0k,FWk)
!-----------------------------------------------------------------
!
! Computes 1d wave and vortical spectra for various quantities, returns them
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END SUBROUTINE wvspectrumc


!*****************************************************************
      SUBROUTINE DoSpAvg(vx,vy,vz,th,istat,nstat,idir,planrc,planio,rv,c1)
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

      END SUBROUTINE padd

!
      SUBROUTINE PVn(gn,a0,vx,vy,vz,th,omega,bvfreq,nt,c1,c2,r1,r2,r3)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes nt-order contribution to potential vorticity
!
! Parameters
!     gn: nt-order enstrophy component (returned)
!     a0: complex temp array of size vx,vy,vz, containing vortical modes
!     vx,
!     vy,
!     vz    : complex velocities, overwritten with normal mode fields
!     th    : complex potential temperature
!     omega : rotation rate
!     bvfreq: Brunt-Vaisalla frequency
!     nt    : order of component (2,3,4). If 1, then quadratic term is
!             computed.
!     c1-2  : complex tmp array
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
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(n,n,ista:iend) :: c1,c2
      REAL   (KIND=GP), INTENT(INOUT), DIMENSION(n,n,ksta:kend) :: r1,r2,r3
      REAL   (KIND=GP), INTENT   (IN)                           :: bvfreq,omega
      REAL   (KIND=GP)                                          :: f,kp,sig,tmp
      INTEGER         , INTENT (IN)                             :: nt
      INTEGER                                                   :: i,j,k

      IF ( nt.lt.1 .or. nt.gt.4 ) THEN
        STOP 'PVn: invalid order'
      ENDIF


      f = 2.0*omega
      gn = (0.0_GP,0.0_GP)

      IF ( nt.le.2 ) THEN
        tmp = 1.0_GP/real(n,kind=GP)**6
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,n
              kp  = sqrt(ka(i)**2+ka(j)**2)
              DO k = 1,n
              ! ks  = sqrt(ka2(k,j,i))
                sig = sqrt(f**2*ka(k)**2+bvfreq**2*kp**2)
                gn(k,j,i) = sig*a0(k,j,i)
              END DO
           END DO
        END DO
!        c1 = a0 
!        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
!        DO k = ksta,kend
!!$omp parallel do if (kend-ksta.lt.nth) private (i)
!           DO j = 1,n
!              DO i = 1,n
!                 r1(i,j,k) = r1(i,j,k)*r1(i,j,k)*tmp
!              END DO
!           END DO
!        END DO
!        CALL fftp3d_real_to_complex(planrc,r1,gn,MPI_COMM_WORLD)
      ENDIF

      IF ( nt.eq.3 ) THEN
        tmp = 1.0_GP/real(n,kind=GP)**6
        CALL wgradt(gn,vx,vy,vz,th,c1,r1,r2,r3)
        CALL derivk3(th,c1,1)
        CALL rotor3(vy,vz,c2,3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
        DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
           DO j = 1,n
              DO k = 1,n
                 c1(k,j,i) = f*c1(k,j,i)-bvfreq*c2(k,j,i)
              END DO
           END DO
        END DO
        CALL fftp3d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
        CALL fftp3d_complex_to_real(plancr,gn,r2,MPI_COMM_WORLD)

!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
                 r1(i,j,k) = r1(i,j,k)*r2(i,j,k)*tmp
              END DO
           END DO
        END DO
        CALL fftp3d_real_to_complex(planrc,r1,gn,MPI_COMM_WORLD)
      ENDIF

      IF ( nt.eq.4 ) THEN
        tmp = 1.0_GP/real(n,kind=GP)**6
        CALL wgradt(gn,vx,vy,vz,th,c1,r1,r2,r3)
        CALL fftp3d_complex_to_real(plancr,gn,r2,MPI_COMM_WORLD)
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
        DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
           DO j = 1,n
              DO i = 1,n
                 r1(i,j,k) = r1(i,j,k)*r1(i,j,k)*tmp
              END DO
           END DO
        END DO
        CALL fftp3d_real_to_complex(planrc,r1,gn,MPI_COMM_WORLD)
      ENDIF

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
               r3(i,j,k) = r3(i,j,k) + r1(i,j,k)*r2(i,j,k)*tmp
            END DO
         END DO
      END DO

      CALL derivk3(th,gn,2)
      CALL fftp3d_complex_to_real(plancr,gn,r1,MPI_COMM_WORLD)
      CALL rotor3(vz,vx,c1,2)
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

      CALL derivk3(th,gn,2)
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

      END SUBROUTINE wgradt


!*****************************************************************
      SUBROUTINE wvzspectrum(a0,vx,vy,vz,th,omega,bvfreq,nmb,nmb1,i2d,c1,c2,c3,r1,r2,r3)
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
!     nmb   : the extension used when writting the file
!     nmb1  : if lenght>0, used to specify range
!     i2d   : do 2D spectra (>0) or not (<=0)
!     c1-3  : complex tmp array
!     r1-3  : real tmp array
!-----------------------------------------------------------------
      USE kes
      USE grid
      USE mpivars
      USE filefmt
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
      CHARACTER(len=*), INTENT(IN) :: nmb,nmb1

!
! isotropic spectra:
      CALL PVn(c3,a0,vx,vy,vz,th,omega,bvfreq,2,c1,c2,r1,r2,r3) 
      CALL spectrumg(c3,1,E2k)   ! isotropic spectrum of Z2
      CALL spectrumg(c3,2,E2kpr) ! reduced perp spectrum of Z2
      CALL spectrumg(c3,3,E2kpa) ! reduced para spectrum of Z2
      IF ( i2d.GT.0 ) THEN
      CALL specaxig (c3,   eaxi) ! axisymm. 2d spectra for Z2
      IF (myrank.eq.0) THEN
        if ( len_trim(nmb1).gt.0 ) then
        OPEN(1,file='kz22D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else
        OPEN(1,file='kz22D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) eaxi
        CLOSE(1)
      ENDIF
      ENDIF

      CALL PVn(c3,a0,vx,vy,vz,th,omega,bvfreq,3,c1,c2,r1,r2,r3)
      CALL spectrumg(c3,1,E3k)   ! isotropic spectrum of Z3
      CALL spectrumg(c3,2,E3kpr) ! reduced perp spectrum of Z3
      CALL spectrumg(c3,3,E3kpa) ! reduced para spectrum of Z3
      IF ( i2d.GT.0 ) THEN
      CALL specaxig (c3,   eaxi) ! axisymm. 2d spectra for Z3
      IF (myrank.eq.0) THEN
        if ( len_trim(nmb1).gt.0 ) then
        OPEN(1,file='kz32D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else
        OPEN(1,file='kz32D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) eaxi
        CLOSE(1)
      ENDIF
      ENDIF
      CALL PVn(c3,a0,vx,vy,vz,th,omega,bvfreq,4,c1,c2,r1,r2,r3)
      CALL spectrumg(c3,1,E4k)   ! isotropic spectrum of Z4
      CALL spectrumg(c3,2,E4kpr) ! reduced perp spectrum of Z4
      CALL spectrumg(c3,3,E4kpa) ! reduced para spectrum of Z4
      IF ( i2d.GT.0 ) THEN
      CALL specaxig (c3,   eaxi) ! axisymm. 2d spectra for Z4
      IF (myrank.eq.0) THEN
        if ( len_trim(nmb1).gt.0 ) then
        OPEN(1,file='kz42D.' // nmb // '_' // trim(nmb1) // '.out',form='unformatted',access='stream')
        else
        OPEN(1,file='kz42D.' // nmb // '.out',form='unformatted',access='stream')
        endif
        WRITE(1) eaxi
        CLOSE(1)
      ENDIF
      ENDIF

      IF (myrank.eq.0) THEN
         if ( len_trim(nmb1).gt.0 ) then
         OPEN(1,file='wvzkspectrum.' // nmb // '_' // trim(nmb1) //'.txt')
         else
         OPEN(1,file='wvzkspectrum.' // nmb // '.txt')
         endif
         DO j=1,n/2+1
           WRITE(1,FMT='(3(E23.15,1X))') E2k(j),E3k(j),E4k(j)
         ENDDO
         CLOSE(1)
      ENDIF
!
! parallel spectra:
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
! perp spectra:
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
!
      RETURN
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

      END SUBROUTINE spectrumg


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
             kperp = int(sqrt(ka(1)**2+ka(j)**2)+0.501)
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
              kperp = int(sqrt(ka(i)**2+ka(j)**2)+0.501)
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

      END SUBROUTINE specaxig
