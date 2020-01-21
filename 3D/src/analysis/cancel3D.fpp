!=================================================================
      PROGRAM CANCEL3D
!=================================================================
! CANCEL3D code (part of the GHOST suite)
!
! Numerically computes the cancellation exponent in 3D 
! simulations with the GHOST code. This tool ONLY works
! with cubic data in (2.pi)^3 domains.
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2011 Ernesto Horne and Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!
! 30 Nov 2011: Main program for all decompositions
!=================================================================

!
! Modules

      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE iovar
      USE grid
      USE fft
      USE ali
      USE var
      USE kes
      IMPLICIT NONE

!
! Arrays for the fields

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: v1,v2,v3
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: r
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: chi

!
! Auxiliary variables

      REAL(KIND=GP)    :: tmp,mean,norm
      REAL(KIND=GP)    :: total_box,total_loc

      INTEGER :: n
      INTEGER :: sini
      INTEGER :: cold
      INTEGER :: curl
      INTEGER :: comp
      INTEGER :: stat
      INTEGER :: nboxes,boxsize
      INTEGER :: domain_end
      INTEGER :: xbox_sta,xbox_end
      INTEGER :: ybox_sta,ybox_end
      INTEGER :: zbox_sta,zbox_end
      INTEGER :: mysta,myend
      INTEGER :: i,j,k,ni,nj,nk

      TYPE(IOPLAN) :: planio

      CHARACTER(len=100) :: odir,idir

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
      CALL io_init(myrank,n,ksta,kend,planio)

!
! Reads from the external file 'cancel.txt' the 
! parameters that will be used to compute the transfer
!     idir : directory for unformatted input (field components)
!     odir : directory for unformatted output (exponent)
!     stat : number of the file to analyze
!     cold : =0 restart a previous computation
!            =1 start a new computation
!     curl : =0 computes cancellation of the field
!            =1 computes cancellation of the vorticity
!            =2 computes cancellation of the helicity
!     comp : 1,2,3 field component (ignored if curl=2)

      IF (myrank.eq.0) THEN
         OPEN(1,file='cancel.txt',status='unknown')
         READ(1,'(a100)') idir
         READ(1,'(a100)') odir
         READ(1,*) stat
         READ(1,*) cold
         READ(1,*) curl
         READ(1,*) comp
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cold,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(curl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(comp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Allocates memory for distributed blocks and the output

      ALLOCATE( v1(n,n,ksta:kend) )
      ALLOCATE( r(n) )
      ALLOCATE( chi(n) )

!
! If continuing a previous run (cold=0), reads 'sini.txt'
! and the arrays with the partition functions

      sini = 1
      WRITE(ext, fmtext) stat

 RI : IF (cold.eq.0) THEN

      IF (myrank.eq.0) THEN
         OPEN(1,file='sini.txt',status='unknown')
         READ(1,*) sini
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/increment.' &
              // ext // '.out',form='unformatted')
         READ(1) r
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/chi.'       &
              // ext // '.out',form='unformatted')
         READ(1) chi
         CLOSE(1)
      ENDIF

      ENDIF RI

!
! Reads the external binary files with the 
! fields components, and computes auxiliary
! fields if needed.
! If required, initializes the FFT library and 
! arrays with wavenumbers. Use FFTW_ESTIMATE in 
! short runs and FFTW_MEASURE in long runs.

      IF (curl.eq.0) THEN
         IF (comp.eq.1) THEN
            CALL io_read(1,idir,'vx',ext,planio,v1)
         ElSE IF (comp.eq.2) THEN
            CALL io_read(1,idir,'vy',ext,planio,v1)
         ELSE
            CALL io_read(1,idir,'vz',ext,planio,v1)
         ENDIF
      ELSE IF ((curl.eq.1).and.(sini.eq.1)) THEN
         ALLOCATE( C1(n,n,ista:iend) )
         ALLOCATE( C2(n,n,ista:iend) )
         CALL fftp3d_create_plan(planrc,(/n,n,n/),FFTW_REAL_TO_COMPLEX, &
             FFTW_ESTIMATE)
         CALL fftp3d_create_plan(plancr,(/n,n,n/),FFTW_COMPLEX_TO_REAL, &
             FFTW_ESTIMATE)
         ALLOCATE( kx(n), ky(n), kz(n) )
         ALLOCATE( kn2(n,n,ista:iend), kk2(n,n,ista:iend) )
         kmax = (real(n,kind=GP)/3.)**2
         tiny = 1e-5
         DO i = 1,n/2
            kx(i) = real(i-1,kind=GP)
            kx(i+n/2) = real(i-n/2-1,kind=GP)
         END DO
	 ky = kx
	 kz = kx
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  kk2(k,j,i) = kx(i)**2+ky(j)**2+kz(k)**2
               END DO
            END DO
         END DO
	 kn2 = kk2
         IF (comp.eq.1) THEN
            CALL io_read(1,idir,'vy',ext,planio,v1)
            CALL fftp3d_real_to_complex(planrc,v1,C1,MPI_COMM_WORLD)
            CALL io_read(1,idir,'vz',ext,planio,v1)
            CALL fftp3d_real_to_complex(planrc,v1,C2,MPI_COMM_WORLD)
         ELSE IF (comp.eq.2) THEN
            CALL io_read(1,idir,'vx',ext,planio,v1)
            CALL fftp3d_real_to_complex(planrc,v1,C1,MPI_COMM_WORLD)
            CALL io_read(1,idir,'vz',ext,planio,v1)
            CALL fftp3d_real_to_complex(planrc,v1,C2,MPI_COMM_WORLD)
         ELSE
            CALL io_read(1,idir,'vx',ext,planio,v1)
            CALL fftp3d_real_to_complex(planrc,v1,C1,MPI_COMM_WORLD)
            CALL io_read(1,idir,'vy',ext,planio,v1)
            CALL fftp3d_real_to_complex(planrc,v1,C2,MPI_COMM_WORLD)
         ENDIF
         CALL rotor3(C1,C2,C1,comp)
         CALL fftp3d_complex_to_real(plancr,C1,v1,MPI_COMM_WORLD)
         tmp = 1./real(n,kind=GP)**3
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  v1(i,j,k) = v1(i,j,k)*tmp
                  v1(i,j,k) = v1(i,j,k)*tmp
                  v1(i,j,k) = v1(i,j,k)*tmp
               END DO
            END DO
         END DO
         IF (comp.eq.1) THEN
            CALL io_write(1,idir,'wx',ext,planio,v1)
         ELSE IF (comp.eq.2) THEN
            CALL io_write(1,idir,'wy',ext,planio,v1)
         ELSE
            CALL io_write(1,idir,'wz',ext,planio,v1)
         ENDIF
         DEALLOCATE ( C1,C2 )
      ELSE IF ((curl.eq.1).and.(sini.gt.1)) THEN
         IF (comp.eq.1) THEN
            CALL io_read(1,idir,'wx',ext,planio,v1)
         ELSE IF (comp.eq.2) THEN
            CALL io_read(1,idir,'wy',ext,planio,v1)
         ELSE
            CALL io_read(1,idir,'wz',ext,planio,v1)
         ENDIF
      ELSE IF ((curl.eq.2).and.(sini.eq.1)) THEN
         ALLOCATE( v2(n,n,ksta:kend) )
         ALLOCATE( v3(n,n,ksta:kend) )
         CALL io_read(1,idir,'vx',ext,planio,v2)
         CALL io_read(1,idir,'wx',ext,planio,v3)
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  v1(i,j,k) = v2(i,j,k)*v3(i,j,k)
               END DO
            END DO
         END DO
         CALL io_read(1,idir,'vy',ext,planio,v2)
         CALL io_read(1,idir,'wy',ext,planio,v3)
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  v1(i,j,k) = v1(i,j,k)+v2(i,j,k)*v3(i,j,k)
               END DO
            END DO
         END DO
         CALL io_read(1,idir,'vz',ext,planio,v2)
         CALL io_read(1,idir,'wz',ext,planio,v3)
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  v1(i,j,k) = v1(i,j,k)+v2(i,j,k)*v3(i,j,k)
               END DO
            END DO
         END DO
         DEALLOCATE ( v2,v3 )
         CALL io_write(1,idir,'h',ext,planio,v1)
      ELSE IF ((curl.eq.2).and.(sini.gt.1)) THEN
         CALL io_read(1,idir,'h',ext,planio,v1)
      ENDIF

!
! Computes the mean value of the field in the entire box

      total_loc = 0.0_GP
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
               total_loc = total_loc+v1(i,j,k)
            END DO
         END DO
      END DO
      CALL MPI_ALLREDUCE(total_loc,mean,1,GC_REAL,MPI_SUM,  &
                        MPI_COMM_WORLD,ierr)
      mean = mean/float(n)**3

!
! Computes the partition functions

 BX:  DO boxsize= sini,n

      chi(boxsize) = 0.0_GP
      r(boxsize) = 2*pi*float(boxsize)/float(n)
      nboxes = n/boxsize 
      domain_end = nboxes*boxsize
      IF ( domain_end.gt.kend ) THEN
         myend = kend
      ELSE
         myend = domain_end
      ENDIF

!
! Computes the normalization for the partition function

      total_loc = 0.0_GP
      DO k = ksta,myend
         DO j = 1,domain_end
            DO i = 1,domain_end
               total_loc = total_loc+abs(v1(i,j,k)-mean)
            END DO
         END DO
      END DO
      CALL MPI_ALLREDUCE(total_loc,norm,1,GC_REAL,MPI_SUM, &
                        MPI_COMM_WORLD,ierr)

 CA:  DO ni = 1,nboxes
      DO nj = 1,nboxes
      DO nk = 1,nboxes

         xbox_sta=(ni-1)*boxsize+1
         xbox_end=ni*boxsize
         ybox_sta=(nj-1)*boxsize+1
         ybox_end=nj*boxsize
         zbox_sta=(nk-1)*boxsize+1
         zbox_end=nk*boxsize

         IF ( zbox_end.gt.kend ) THEN
            myend = kend
         ELSE
            myend = zbox_end
         ENDIF
         IF ( zbox_sta.le.ksta ) THEN
            mysta = ksta
         ELSE
            mysta = zbox_sta
         ENDIF

         total_loc = 0.0_GP
         DO k = mysta,myend
            DO j = ybox_sta,ybox_end
               DO i = xbox_sta,xbox_end
                  total_loc = total_loc+v1(i,j,k)-mean
               END DO
            END DO
         END DO
         CALL MPI_ALLREDUCE(total_loc,total_box,1,GC_REAL,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
         chi(boxsize) = chi(boxsize)+abs(total_box/norm)

      END DO
      END DO
      END DO CA

 RO : IF (myrank.eq.0) THEN
         OPEN(1,file='sini.txt',status='unknown')
         WRITE(1,*) boxsize+1
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/increment.'  &
              // ext // '.out' ,form='unformatted')
         WRITE(1) r
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/chi.'        &
              // ext // '.out' ,form='unformatted')
         WRITE(1) chi
         CLOSE(1)
      ENDIF RO

      END DO BX

!
! End of CANCEL3D

      CALL MPI_FINALIZE(ierr)
      DEALLOCATE( v1 )
      DEALLOCATE( r,chi )

      END PROGRAM CANCEL3D
