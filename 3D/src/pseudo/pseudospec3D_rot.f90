!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Extra subroutines to compute the energy spectrum and energy 
! transfer functions in the HD, MHD, and Hall-MHD equations 
! in 3D in the anisotropic case (e.g. in the rotating frame 
! or with an imposed external magnetic field). Quantities 
! parallel and perpendicular to the z direction in Fourier 
! space can be computed. You should use the FFTPLANS and 
! MPIVARS modules (see the file 'fftp_mod.f90') in each 
! program that calls any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar
!
! 7 Aug 2010: New specperp and 2D spectrum (T. Teitelbaum)
!=================================================================

!*****************************************************************
      SUBROUTINE specpara(a,b,c,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the reduced energy and helicity power spectrum 
! in the direction parallel to the preferred direction 
! (rotation, gravity, or uniform magnetic field, assumed
! to be in the z-direction). As a result, the k-shells
! are planes with normal (0,0,kz), kz = Dkz*(0,...,nz/2).
! Normalization of the reduced spectrum is such that
! E = sum[E(kz).Dkz], where Dkz is the width of the Fourier
! shells in kz. The output is written to a file by the
! first node.        
!
! Output files contain:
! 'kspecpara.XXX.txt': kz, Ev(kz), Ev_perp(kz), Ev_z(kz)
!   [Ev: kinetic energy; Ev_perp: energy in vx,vy; Ev_z: same in vz]        
! 'mspecpara.XXX.txt': kz, Eb(kz), Eb_perp(kz), Eb_z(kz)
! 'khelipara.XXX.txt': kz, Hv(kz), Hv_perp(kz), Hv_z(kz)
!   [Hv: kinetic helicity; Hv_perp: v_perp.w_perp, Hv_z: vz.wz]
! 'mhelipara.XXX.txt': kz, Hb(kz), Hb_perp(kz), Hb_z(kz)
! 'ghelipara.XXX.txt': kz, G(kz)  ,G_perp(kz) , G_z(kz) 
!   [Generalized helicity in Hall-MHD]
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     kin: =2 skips energy spectrum computation
!          =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Ek,Ektot
      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Ekh,Ekhtot
      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Ekv,Ekvtot
      DOUBLE PRECISION    :: tmq,tmr
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO k = 1,nz/2+1
         Ek (k) = 0.0D0
         Ekh(k) = 0.0D0
         Ekv(k) = 0.0D0
      END DO
!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2)*tmp
                     tmr = (abs(c(k,j,1))**2)*tmp
!$omp critical
                     Ek (kmn) = Ek (kmn)+tmq+tmr                       
                     Ekh(kmn) = Ekh(kmn)+tmq
                     Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2)*tmp
                        tmr = 2*(abs(c(k,j,i))**2)*tmp
!$omp critical
                        Ek (kmn) = Ek (kmn)+tmq+tmr
                        Ekh(kmn) = Ekh(kmn)+tmq
                        Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2)*tmp
                        tmr = 2*(abs(c(k,j,i))**2)*tmp
!$omp critical
                        Ek (kmn) = Ek (kmn)+tmq+tmr
                        Ekh(kmn) = Ekh(kmn)+tmq
                        Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2)*tmp
                     tmr = (abs(c3(k,j,1))**2)*tmp
!$omp critical
                     Ek (kmn) = Ek (kmn)+tmq+tmr
                     Ekh(kmn) = Ekh(kmn)+tmq
                     Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical

                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2)*tmp
                        tmr = 2*(abs(c3(k,j,i))**2)*tmp
!$omp critical
                        Ek (kmn) = Ek (kmn)+tmq+tmr
                        Ekh(kmn) = Ekh(kmn)+tmq
                        Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                     ENDIF
                  END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2)*tmp
                        tmr = 2*(abs(c3(k,j,i))**2)*tmp
!$omp critical
                        Ek (kmn) = Ek (kmn)+tmq+tmr
                        Ekh(kmn) = Ekh(kmn)+tmq
                        Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      IF (kin.le.1) THEN
         CALL MPI_REDUCE(Ek ,Ektot ,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekh,Ekhtot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekv,Ekvtot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kspecpara.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mspecpara.' // nmb // '.txt')
            ENDIF
            DO k = 1,nz/2+1
               WRITE(1,FMT='(E13.6,E23.15,E23.15,E23.15)') &
                              Dkz*(k-1),.5_GP*Ektot(k)*Lz, &
                    .5_GP*Ekhtot(k)*Lz,.5_GP*Ekvtot(k)*Lz
            END DO
            CLOSE(1)
         ENDIF
      END IF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO k = 1,nz/2+1
            Ek (k) = 0.0D0
            Ekh(k) = 0.0D0
            Ekv(k) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+            &
                            real(b(k,j,1)*conjg(c2(k,j,1))))*tmp
                     tmr = (real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp critical
                     Ek (kmn) = Ek (kmn)+tmq+tmr
                     Ekh(kmn) = Ekh(kmn)+tmq
                     Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(c2(k,j,i))))*tmp
                        tmr = 2*(real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp critical
                        Ek(kmn) = Ek(kmn)+tmq+tmr
                        Ekh(kmn) = Ekh(kmn)+tmq
                        Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+       &
                                 real(b(k,j,i)*conjg(c2(k,j,i))))*tmp
                        tmr = 2*(real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp critical
                        Ek (kmn) = Ek (kmn)+tmq+tmr
                        Ekh(kmn) = Ekh(kmn)+tmq
                        Ekv(kmn) = Ekv(kmn)+tmr
!$omp end critical
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek ,Ektot ,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekh,Ekhtot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekv,Ekvtot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                         0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='khelipara.' // nmb // '.txt')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file='mhelipara.' // nmb // '.txt')
            ELSE
               OPEN(1,file='ghelipara.' // nmb // '.txt')
            ENDIF
            DO k = 1,nz/2+1
               WRITE(1,FMT='(E13.6,E23.15,E23.15,E23.15)') &
                     Dkz*(k-1),Ektot(k)*Lz,Ekhtot(k)*Lz,Ekvtot(k)*Lz
            ENDDO
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE specpara

!*****************************************************************
      SUBROUTINE specperp(a,b,c,nmb,kin,hel)
!-----------------------------------------------------------------
!
! Computes the reduced energy and helicity power spectrum 
! in the direction perpendicular to the preferred direction 
! (rotation, gravity, or uniform magnetic field, assumed to
! be in the z-direction. The k-shells are cylindrical
! surfaces with kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2).
! It also computes the spectrum of 2D modes with kz = 0 of
! (x,y)-field components and the z-field component separately.
! Normalization of the reduced spectrum is such that
! E = sum[E(kperp).Dkk], where Dkk is the width of the Fourier
! shells. The output is written to a file with three columns
! by the first node.
!
! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'kspecperp.XXX.txt': kp, E_v(kp), ev_x,y(kp,kz=0), ev_z(kp,kz=0)
!   [E_v: kin. ener.; ev_x,y: 2D spec. for vx,vy; ev_z: same for vz]
! 'mspecperp.XXX.txt': kp, E_b(kp), eb_x,y(kp,kz=0), eb_z(kp,kz=0)
! 'kheliperp.XXX.txt': kp, H_v(kp), hv_x,y(kp,kz=0), hv_z(kp,kz=0)
! 'mheliperp.XXX.txt': kp, H_b(kp), hb_x,y(kp,kz=0), hb_z(kp,kz=0)
! 'gheliperp.XXX.txt': kp, G(kz),   g(kp,kz=0),      g(kp,kz=0)
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     kin: =2 skips energy spectrum computation
!          =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ek,Ektot
      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ekp,Eptot
      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ekz,Eztot
      DOUBLE PRECISION    :: tmq,tmr
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,nmaxperp/2+1
         Ek(i) = 0.0D0
         Ekp(i) = 0.0D0
         Ekz(i) = 0.0D0
      END DO
!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq,tmr)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  tmq = (abs(a(1,j,1))**2+abs(b(1,j,1))**2)*tmp
                  tmr = (abs(c(1,j,1))**2)*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ekz(kmn) = Ekz(kmn)+tmr
                  Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                  DO k = 2,nz
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2+          &
                            abs(c(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     tmq = 2*(abs(a(1,j,i))**2+abs(b(1,j,i))**2)*tmp
                     tmr = 2*abs(c(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,nz
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+     &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
          ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     tmq = 2*(abs(a(1,j,i))**2+abs(b(1,j,i))**2)*tmp
                     tmr = 2*abs(c(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,nz
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+     &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
          ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq,tmr)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  tmq = (abs(c1(1,j,1))**2+abs(c2(1,j,1))**2)*tmp
                  tmr = abs(c3(1,j,1))**2*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ekz(kmn) = Ekz(kmn)+tmr
                  Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                  DO k = 2,nz
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+        &
                            abs(c3(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     tmq = 2*(abs(c1(1,j,i))**2+abs(c2(1,j,i))**2)*tmp
                     tmr = 2*abs(c3(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,nz
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+   &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     tmq = 2*(abs(c1(1,j,i))**2+abs(c2(1,j,i))**2)*tmp
                     tmr = 2*abs(c3(1,j,i))**2*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,nz
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+   &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      IF (kin.le.1) THEN
         CALL MPI_REDUCE(Ek,Ektot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,    &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekp,Eptot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekz,Eztot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kspecperp.' // nmb // '.txt')
            ELSE
               OPEN(1,file='mspecperp.' // nmb // '.txt')
            ENDIF
            DO j = 1,nmaxperp/2+1
               WRITE(1,FMT='(E13.6,E23.15,E23.15,E23.15)') &
                             Dkk*(j-1),.5_GP*Ektot(j)/Dkk, &
                    .5_GP*Eptot(j)/Dkk,.5_GP*Eztot(j)/Dkk
            END DO
            CLOSE(1)
         ENDIF
      END IF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,nmaxperp/2+1
            Ek(i) = 0.0D0
            Ekp(i) = 0.0D0
            Ekz(i) = 0.0D0
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq,tmr)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  tmq = (real(a(1,j,1)*conjg(c1(1,j,1)))+              &
                         real(b(1,j,1)*conjg(c2(1,j,1))))*tmp
                  tmr = real(c(1,j,1)*conjg(c3(1,j,1)))*tmp
!$omp critical
                  Ekp(kmn) = Ekp(kmn)+tmq
                  Ekz(kmn) = Ekz(kmn)+tmr
                  Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                  DO k = 2,nz
                     tmq = (real(a(k,j,1)*conjg(c1(k,j,1)))+           &
                            real(b(k,j,1)*conjg(c2(k,j,1)))+           &
                            real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     tmq = 2*(real(a(1,j,i)*conjg(c1(1,j,i)))+         &
                           real(b(1,j,i)*conjg(c2(1,j,i))))*tmp
                     tmr = 2*real(c(1,j,i)*conjg(c3(1,j,i)))*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,nz
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+      &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+      &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq,tmr)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq,tmr)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     tmq = 2*(real(a(1,j,i)*conjg(c1(1,j,i)))+         &
                           real(b(1,j,i)*conjg(c2(1,j,i))))*tmp
                     tmr = 2*real(c(1,j,i)*conjg(c3(1,j,i)))*tmp
!$omp critical
                     Ekp(kmn) = Ekp(kmn)+tmq
                     Ekz(kmn) = Ekz(kmn)+tmr
                     Ek(kmn) = Ek(kmn)+tmq+tmr
!$omp end critical
                     DO k = 2,nz
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+      &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+      &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,   &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekp,Eptot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,  &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         CALL MPI_REDUCE(Ekz,Eztot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,  &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file='kheliperp.' // nmb // '.txt')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file='mheliperp.' // nmb // '.txt')
            ELSE
               OPEN(1,file='gheliperp.' // nmb // '.txt')
            ENDIF
            DO j = 1,nmaxperp/2+1
               WRITE(1,FMT='(E13.6,E23.15,E23.15,E23.15)') Dkk*(j-1),   &
                              Ektot(j)/Dkk,Eptot(j)/Dkk,Eztot(j)/Dkk
            END DO
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE specperp

!****************************************************************
      SUBROUTINE entpara(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy transfer in the direction parallel 
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are planes with 
! normal (0,0,kz), kz = Dkz*(0,...,nz/2). Normalization of
! the transfer function is such that the flux is
! Pi = -sum[T(kz).Dkz], where Dkz is the width of the
! Fourier shells. The output is written to a file by the
! first node.
!
! Output files contain:
! 'ktranpara.XXX.txt': kz, Tv(kz) (kinetic energy transfer function)
! 'mtranpara.XXX.txt': kz, Tb(kz) (magnetic energy transfer)
! 'jtranpara.XXX.txt': kz, Hj(kz) (Lorentz force work)
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO k = 1,nz/2+1
         Ek(k) = 0.0D0
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
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
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = kk2(k,j,1) *                             &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  ENDIF
               END DO
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*kk2(k,j,i) *                        &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     ENDIF
                 END DO
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  DO k = 1,nz
                     kmn = int(abs(kz(k))*Lz+1)
                     IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                        tmq = 2*kk2(k,j,i)*                         &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
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
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='mtranpara.' // nmb // '.txt')
         ELSEIF (kin.eq.1) THEN
            OPEN(1,file='ktranpara.' // nmb // '.txt')
         ELSE
            OPEN(1,file='jtranpara.' // nmb // '.txt')
         ENDIF
         DO k = 1,nz/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkz*(k-1),Ektot(k)*Lz
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE entpara

!*****************************************************************
      SUBROUTINE entperp(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy transfer in the direction perpendicular
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are cylindrical 
! surfaces with kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2).
! Normalization of the transfer function is such that the
! flux is Pi = -sum[T(kperp).Dkk], where Dkk is the width of the
! Fourier shells. The output is written to a file by the first
! node.
!
! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'ktranperp.XXX.txt': kp, Tv(kp) (kinetic energy transfer function)
! 'mtranperp.XXX.txt': kp, Tb(kp) (magnetic energy transfer)
! 'jtranperp.XXX.txt': kp, Hj(kp) (Lorentz force work)
!
! Parameters
!     a  : field component in the x-direction
!     b  : field component in the y-direction
!     c  : field component in the z-direction
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic energy transfer
!          =1 computes the kinetic energy transfer
!          =2 computes the Lorentz force transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads      
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Ek,Ektot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,nmaxperp/2+1
         Ek(i) = 0.0D0
      END DO
!
! Computes the kinetic energy transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF ((kin.eq.1).or.(kin.eq.2)) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     tmq = (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        tmq = 2*(real(a(k,j,i)*conjg(d(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(e(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
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
!$omp parallel do private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     tmq = kk2(k,j,1) *                             &
                           (real(a(k,j,1)*conjg(d(k,j,1)))+         &
                            real(b(k,j,1)*conjg(e(k,j,1)))+         &
                            real(c(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                     Ek(kmn) = Ek(kmn)+tmq
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        tmq = 2*kk2(k,j,i) *                        &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        tmq = 2*kk2(k,j,i)*                         &
                              (real(a(k,j,i)*conjg(d(k,j,i)))+      &
                               real(b(k,j,i)*conjg(e(k,j,i)))+      &
                               real(c(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn) = Ek(kmn)+tmq
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,   &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='mtranperp.' // nmb // '.txt')
         ELSEIF (kin.eq.1) THEN
            OPEN(1,file='ktranperp.' // nmb // '.txt')
         ELSE
            OPEN(1,file='jtranperp.' // nmb // '.txt')
         ENDIF
         DO j = 1,nmaxperp/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(j-1),Ektot(j)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE entperp

!*****************************************************************
      SUBROUTINE heltpara(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the helicity transfer in the direction parallel 
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are planes with 
! normal (0,0,kz), kz = Dkz*(0,...,nz/2). Normalization of
! the transfer function is such that the flux is
! Pi = -sum[T(kz).Dkz], where Dkz is the width of the
! Fourier shells. The output is written to a file by the
! first node.
!
! Output files contain:
! 'hktranpara.XXX.txt': kz, TH_v(kz) (kinetic helicity transfer)
! 'hmtranpara.XXX.txt': kz, TH_b(kz) (magnetic helicity transfer)
!
! Parameters
!     a  : field component in the x-direction (v or a)
!     b  : field component in the y-direction (v or a)
!     c  : field component in the z-direction (v or a)
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic helicity transfer
!          =1 computes the kinetic helicity transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nz/2+1) :: Hk,Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO k = 1,nz/2+1
         Hk(k) = 0.0D0
      END DO
!
! Computes the helicity transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            DO k = 1,nz
               kmn = int(abs(kz(k))*Lz+1)
               IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                  tmq = (real(c1(k,j,1)*conjg(d(k,j,1)))+          &
                         real(c2(k,j,1)*conjg(e(k,j,1)))+          &
                         real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                  Hk(kmn) = Hk(kmn)+tmq
               ENDIF
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  ENDIF
              END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               DO k = 1,nz
                  kmn = int(abs(kz(k))*Lz+1)
                  IF ((kmn.gt.0).and.(kmn.le.nz/2+1)) THEN
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  ENDIF
               END DO
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Hk,Hktot,nz/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='hmtranpara.' // nmb // '.txt')
         ELSE
            OPEN(1,file='hktranpara.' // nmb // '.txt')
         ENDIF
         DO k = 1,nz/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkz*(k-1),Hktot(k)*Lz
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE heltpara

!*****************************************************************
      SUBROUTINE heltperp(a,b,c,d,e,f,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the helicity transfer in the direction perpendicular
! to the preferred direction (rotation or uniform magnetic 
! field) in 3D Fourier space. The k-shells are cylindrical 
! surfaces with kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2).
! Normalization of the transfer function is such that the
! flux is Pi = -sum[T(kperp).Dkk], where Dkk is the width of the
! Fourier shells. The output is written to a file by the first
! node.
!
! Output files contain [kp = Dkk*sqrt(kx**2+ky**2)]:
! 'hktranperp.XXX.txt': kp, TH_v(kp) (kinetic helicity transfer)
! 'hmtranperp.XXX.txt': kp, TH_b(kp) (magnetic helicity transfer)
!
! Parameters
!     a  : field component in the x-direction (v or a)
!     b  : field component in the y-direction (v or a)
!     c  : field component in the z-direction (v or a)
!     d  : nonlinear term in the x-direction
!     e  : nonlinear term in the y-direction
!     f  : nonlinear term in the z-direction
!     nmb: the extension used when writting the file
!     kin: =0 computes the magnetic helicity transfer
!          =1 computes the kinetic helicity transfer
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmaxperp/2+1) :: Hk,Hktot
      DOUBLE PRECISION    :: tmq
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP)       :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      INTEGER             :: kmn
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Hk to zero
!
      DO i = 1,nmaxperp/2+1
         Hk(i) = 0.0D0
      END DO
!
! Computes the helicity transfer
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      CALL rotor3(b,c,c1,1)
      CALL rotor3(a,c,c2,2)
      CALL rotor3(a,b,c3,3)
      IF (ista.eq.1) THEN
!$omp parallel do private (k,kmn,tmq)
         DO j = 1,ny
            kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
            IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
               DO k = 1,nz
                  tmq = (real(c1(k,j,1)*conjg(d(k,j,1)))+          &
                         real(c2(k,j,1)*conjg(e(k,j,1)))+          &
                         real(c3(k,j,1)*conjg(f(k,j,1))))*tmp
!$omp atomic
                  Hk(kmn) = Hk(kmn)+tmq
               END DO
            ENDIF
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmn,tmq)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq
                  END DO
               ENDIF
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmn,tmq)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     tmq = 2*(real(c1(k,j,i)*conjg(d(k,j,i)))+     &
                              real(c2(k,j,i)*conjg(e(k,j,i)))+     &
                              real(c3(k,j,i)*conjg(f(k,j,i))))*tmp
!$omp atomic
                     Hk(kmn) = Hk(kmn)+tmq                       
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Hk,Hktot,nmaxperp/2+1,MPI_DOUBLE_PRECISION,  &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='hmtranperp.' // nmb // '.txt')
         ELSE
            OPEN(1,file='hktranperp.' // nmb // '.txt')
         ENDIF
         DO j = 1,nmaxperp/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*(j-1),Hktot(j)/Dkk
         END DO
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE heltperp

!*****************************************************************
      SUBROUTINE spec2D(a,b,c,nmb,dir,kin,hel)
!-----------------------------------------------------------------
!
! Computes the axysimmetric energy and helicity power spectrum. 
! The spectrum is angle-averaged in the azimuthal direction, 
! and depends on two wavenumbers, (kperp,kpara), with
! kperp = Dkk*(0,...,max{nx*Dkx/Dkk,nyDky/Dkk}/2) and
! kpara = Dkz*(0,...,nz/2). The actual values of these wavenumbers
! can be found in the first column of 'kspecpara' and 'kspecperp'
! files. This spectrum is not normalized (i.e., the energy is not
! divided by 2, it is not divided by the area of the bins in
! Fourier space, nor by sin(theta) to obtain curcles in the
! isotropic case). The output is written to a binary file by the
! first node.
!
! Output files contain:
! 'odir/kspec2D.XXX.out': kinetic energy 2D spectrum ev(kperp,kpara)
! 'odir/mspec2D.XXX.out': magnetic energy spectrum   eb(kperp,kpara)
! 'odir/kheli2D.XXX.out': kinetic helicity spectrum  hv(kperp,kpara)
! 'odir/mheli2D.XXX.out': magnetic helicity spectrum hv(kperp,kpara)
! 'odir/gheli2D.XXX.out': generalized helicity spec. g(kperp,kpara)
!        
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     nmb: the extension used when writting the file
!     dir: directory where the files are written
!     kin: =2 skips energy spectrum computation
!          =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!     hel: =1 computes the helicity spectrum
!          =0 skips helicity spectrum computation
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), DIMENSION(nz,ny,ista:iend)          :: c1,c2,c3
      REAL(KIND=GP),    DIMENSION(nmaxperp/2+1,nz/2+1)      :: Ek,Ektot
      REAL(KIND=GP)       :: tmq,tmp
      INTEGER, INTENT(IN) :: kin,hel
      INTEGER             :: i,j,k
      INTEGER             :: kmn,kmz
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*),   INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,nmaxperp/2+1
         DO k = 1,nz/2+1
            Ek(i,k) = 0.0_GP
         END DO
      END DO
!
! Computes the curl of the field if needed
!
      IF ((kin.eq.0).or.(hel.eq.1)) THEN
         CALL rotor3(b,c,c1,1)
         CALL rotor3(a,c,c2,2)
         CALL rotor3(a,b,c3,3)
      ENDIF
!
! Computes the kinetic energy spectrum
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (kin.eq.1) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmz,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     kmz = int(abs(kz(k))*Lz+1)
                     IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                     tmq = (abs(a(k,j,1))**2+abs(b(k,j,1))**2+        &
                            abs(c(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmz,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        kmz = int(abs(kz(k))*Lz+1)
                        IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmz,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        kmz = int(abs(kz(k))*Lz+1)
                        IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                        tmq = 2*(abs(a(k,j,i))**2+abs(b(k,j,i))**2+   &
                                 abs(c(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the magnetic energy spectrum
!
      ELSE IF (kin.eq.0) THEN
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmz,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                     kmz = int(abs(kz(k))*Lz+1)
                     IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                     tmq = (abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+      &
                            abs(c3(k,j,1))**2)*tmp
!$omp atomic
                     Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                     ENDIF
                  END DO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmz,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        kmz = int(abs(kz(k))*Lz+1)
                        IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmz,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        kmz = int(abs(kz(k))*Lz+1)
                        IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                        tmq = 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                                 abs(c3(k,j,i))**2)*tmp
!$omp atomic
                        Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      IF (kin.le.1) THEN
         CALL MPI_REDUCE(Ek,Ektot,(nmaxperp/2+1)*(nz/2+1),GC_REAL,    &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file=trim(dir) // '/' // 'kspec2D.' // nmb //   &
                    '.out',form='unformatted')
            ELSE
               OPEN(1,file=trim(dir) // '/' // 'mspec2D.' // nmb //   &
                    '.out',form='unformatted')
            ENDIF
            WRITE(1) Ektot
            CLOSE(1)
         ENDIF
      ENDIF
!
! Computes the helicity spectrum
!
      IF (hel.eq.1) THEN
         DO i = 1,nmaxperp/2+1
            DO k = 1,nz/2+1
               Ek(i,k) = 0.0_GP
            END DO
         END DO
         IF (ista.eq.1) THEN
!$omp parallel do private (k,kmz,kmn,tmq)
            DO j = 1,ny
               kmn = int(sqrt(kx(1)**2+ky(j)**2)/Dkk+1)
               IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                  DO k = 1,nz
                  kmz = int(abs(kz(k))*Lz+1)
                  IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                  tmq = (real(a(k,j,1)*conjg(c1(1,j,1)))+       &
                         real(b(k,j,1)*conjg(c2(1,j,1)))+       &
                         real(c(k,j,1)*conjg(c3(k,j,1))))*tmp
!$omp atomic
                  Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                  ENDIF
                  ENDDO
               ENDIF
            END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k,kmz,kmn,tmq)
            DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k,kmz,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        kmz = int(abs(kz(k))*Lz+1)
                        IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k,kmz,kmn,tmq)
            DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k,kmz,kmn,tmq)
               DO j = 1,ny
                  kmn = int(sqrt(kx(i)**2+ky(j)**2)/Dkk+1)
                  IF ((kmn.gt.0).and.(kmn.le.nmaxperp/2+1)) THEN
                     DO k = 1,nz
                        kmz = int(abs(kz(k))*Lz+1)
                        IF ((kmz.gt.0).and.(kmz.le.nz/2+1)) THEN
                        tmq = 2*(real(a(k,j,i)*conjg(c1(k,j,i)))+    &
                                 real(b(k,j,i)*conjg(c2(k,j,i)))+    &
                                 real(c(k,j,i)*conjg(c3(k,j,i))))*tmp
!$omp atomic
                        Ek(kmn,kmz) = Ek(kmn,kmz)+tmq
                        ENDIF
                     END DO
                  ENDIF
               END DO
            END DO
         ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
         CALL MPI_REDUCE(Ek,Ektot,(nmaxperp/2+1)*(nz/2+1),GC_REAL,   &
                         MPI_SUM,0,MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (kin.eq.1) THEN
               OPEN(1,file=trim(dir) // '/' // 'kheli2D.' // nmb //  &
                    '.out',form='unformatted')
            ELSE IF (kin.eq.0) THEN
               OPEN(1,file=trim(dir) // '/' // 'mheli2D.' // nmb //  &
                    '.out',form='unformatted')
            ELSE
               OPEN(1,file=trim(dir) // '/' // 'gheli2D.' // nmb //  &
                    '.out',form='unformatted')
            ENDIF
            WRITE(1) Ektot
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE spec2d

!*****************************************************************
      SUBROUTINE write_fourier(a,fname,nmb,dir)
!-----------------------------------------------------------------
!
! Writes slices with the Fourier modes in the planes kx=0,
! ky=0, and kz=0. The output is written to a binary file
! by the first node.
!
! Output files contain [field = vx,vy,...]:
! 'odir/field_kx': field(kz,ky,kx=0) with size (nz,ny)
! 'odir/field_ky': field(kz,ky=0,kx) with size (nz,nx/2+1)
! 'odir/field_kz': field(kz=0,ky,kx) with size (ny,nx/2+1)
!   [Wavenumbers are k_i=(0,1,...)/L_i (i=x,y,z)]
!
! Parameters
!     a    : input matrix
!     fname: name of the field component
!     nmb  : the extension used when writting the file
!     dir  : directory where the files are written
!
      USE fprecision
      USE commtypes
      USE kes
      USE grid
      USE mpivars
      USE filefmt
      USE boxsize
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), DIMENSION(ny,nx/2+1)             :: uk1,ut1
      COMPLEX(KIND=GP), DIMENSION(nz,nx/2+1)             :: uk2,ut2
      COMPLEX(KIND=GP), DIMENSION(nz,ny)                 ::     ut3
      REAL                :: tmp
      INTEGER             :: i,j,k
      CHARACTER(len=100), INTENT(IN) :: dir
      CHARACTER(len=*),   INTENT(IN) :: nmb,fname

      tmp = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
!
! Computes the kz=0 slice
!
!$omp parallel do private (j)
      DO i = ista,iend
         DO j = 1,ny
            uk1(j,i) = a(1,j,i)*tmp
         END DO
      END DO
!
! Computes the ky=0 slice
!
!$omp parallel do private (k)
      DO i = ista,iend
         DO k = 1,nz
            uk2(k,i) = a(k,1,i)*tmp
         END DO
      END DO
!
! Computes the kx=0 slice
!
      IF (myrank.eq.0) THEN
!$omp parallel do private (k)
         DO j = 1,ny
            DO k = 1,nz
               ut3(k,j) = a(k,j,1)*tmp
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(uk1,ut1,ny*(nx/2+1),GC_COMPLEX,          &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(uk2,ut2,nz*(nx/2+1),GC_COMPLEX,            &
                      MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(dir) // '/' // trim(fname) // '_kz.' // &
              nmb // '.out',form='unformatted')
         WRITE(1) ut1
         CLOSE(1)
         OPEN(1,file=trim(dir) // '/' // trim(fname) // '_ky.' // &
              nmb // '.out',form='unformatted')
         WRITE(1) ut2
         CLOSE(1)
         OPEN(1,file=trim(dir) // '/' // trim(fname) // '_kx.' // &
              nmb // '.out',form='unformatted')
         WRITE(1) ut3
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE write_fourier
