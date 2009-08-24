!=================================================================
      PROGRAM STRUCTURE3D
!=================================================================
! STRUCTURE3D code (part of the GHOST suite)
!
! Numerically computes the structure functions in 3D HD, 
! MHD, and Hall-MHD simulations with the GHOST code. An 
! SO(2) or SO(3) decomposition is used to compute field 
! increments depending on the solver.
!
! Conditional compilation options:
!           HD_SOL   Hydrodynamic SO(3) structure functions
!           ROTH_SOL SO(2) HD structure functions (rotation)
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2008 Luis Martin and Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!
! 2 Jul 2008: Main program for all decompositions and solvers
! 4 Feb 2009: Option for increments of helicity and vorticity
!=================================================================

!
! Definitions for conditional compilation

#ifdef HD_SOL
#define SO3_
#endif

#ifdef ROTH_SOL
#define SO2_
#endif

!
! Modules

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
! Arrays for the fields and structure functions

      COMPLEX, ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2,C3
      COMPLEX, ALLOCATABLE, DIMENSION (:,:,:) :: C4,C5

      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: vxl,vyl,vzl
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: wx,wy,wz
      REAL, ALLOCATABLE, DIMENSION (:,:,:) :: wxl,wyl,wzl

      REAL, ALLOCATABLE, DIMENSION (:,:)   :: sp
      REAL, ALLOCATABLE, DIMENSION (:)     :: r

!
! Auxiliary variables

      REAL    :: tmp,norm,spaux

#ifdef SO3_
      INTEGER, PARAMETER        :: ngen = 146
#endif
#ifdef SO2_
      INTEGER, PARAMETER        :: ngen = 26
#endif
      INTEGER, DIMENSION (ngen) :: j1,j2,j3
      INTEGER :: gini
      INTEGER :: curl
      INTEGER :: stat
      INTEGER :: pmax
      INTEGER :: maxd,dis
      INTEGER :: i,j,k,l,d,p,q

      TYPE(IOPLAN) :: planio

      CHARACTER(len=100) :: odir,idir

!
! Initializes the MPI and I/O libraries

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,n/2+1,nprocs,myrank,ista,iend)
      CALL range(1,n,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,n,ksta,kend,planio)

!
! Reads from the external file 'structure.txt' the 
! parameters that will be used to compute the transfer
!     idir : directory for unformatted input (field components)
!     odir : directory for unformatted output (transfers)
!     stat : number of the file to analyze
!     gini : =1 start a new computation
!            >1 restart a previous computation from this generator
!     curl : =0 computes increments of the field
!            =1 computes increments of the helicity of the field
!            =2 computes increments of the vorticity of the field
!     pmax : maximum order p computed

      IF (myrank.eq.0) THEN
         OPEN(1,file='structure.txt',status='unknown')
         READ(1,'(a100)') idir
         READ(1,'(a100)') odir
         READ(1,*) stat
         READ(1,*) gini
         READ(1,*) curl
         READ(1,*) pmax
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gini,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(curl,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(pmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Allocates memory for distributed blocks

      IF (curl.gt.0) THEN
         ALLOCATE( C1(n,n,ista:iend), C2(n,n,ista:iend) )
         ALLOCATE( C3(n,n,ista:iend), C4(n,n,ista:iend) )
         ALLOCATE( C5(n,n,ista:iend) )
         ALLOCATE( wx(n,n,ksta:kend) )
         ALLOCATE( wy(n,n,ksta:kend) )
         ALLOCATE( wz(n,n,ksta:kend) )
      ENDIF
      ALLOCATE( vx(n,n,ksta:kend) )
      ALLOCATE( vy(n,n,ksta:kend) )
      ALLOCATE( vz(n,n,ksta:kend) )
      ALLOCATE( sp(n/2-1,pmax) )
      ALLOCATE( r(n/2-1) )

!
! Reads the external binary files with the 
! fields components, and computes auxiliary
! fields if needed.
! If required, initializes the FFT library and 
! arrays with wavenumbers. Use FFTW_ESTIMATE in 
! short runs and FFTW_MEASURE in long runs

      WRITE(ext, fmtext) stat
      CALL io_read(1,idir,'vx',ext,planio,vx)
      CALL io_read(1,idir,'vy',ext,planio,vy)
      CALL io_read(1,idir,'vz',ext,planio,vz)
      IF (curl.gt.0) THEN
         CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
             FFTW_ESTIMATE)
         CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
             FFTW_ESTIMATE)
         ALLOCATE( ka(n), ka2(n,n,ista:iend) )
         kmax = (float(n)/3.)**2
         tiny = 1e-5
         DO i = 1,n/2
            ka(i) = float(i-1)
            ka(i+n/2) = float(i-n/2-1)
         END DO
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n
                  ka2(k,j,i) = ka(i)**2+ka(j)**2+ka(k)**2
               END DO
            END DO
         END DO
         CALL fftp3d_real_to_complex(planrc,vx,C1,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,vy,C2,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,vz,C3,MPI_COMM_WORLD)
         CALL rotor3(C2,C3,C4,1)
         CALL rotor3(C1,C3,C5,2)
         CALL rotor3(C1,C2,C3,3)
         CALL fftp3d_complex_to_real(plancr,C4,wx,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C5,wy,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(plancr,C3,wz,MPI_COMM_WORLD)
         DEALLOCATE ( C1,C2,C3,C4,C5,ka,ka2 )
         CALL fftp3d_destroy_plan(plancr)
         CALL fftp3d_destroy_plan(planrc)
      ENDIF
      ALLOCATE ( vxl(n,n,ksta:kend) )
      ALLOCATE ( vyl(n,n,ksta:kend) )
      ALLOCATE ( vzl(n,n,ksta:kend) )
      IF (curl.eq.1) THEN
         tmp = 1./float(n)**3
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  wx(i,j,k) = wx(i,j,k)*tmp
                  wy(i,j,k) = wy(i,j,k)*tmp
                  wz(i,j,k) = wz(i,j,k)*tmp
               END DO
            END DO
         END DO
         ALLOCATE ( wxl(n,n,ksta:kend) )
         ALLOCATE ( wyl(n,n,ksta:kend) )
         ALLOCATE ( wzl(n,n,ksta:kend) )
      ENDIF
      IF (curl.eq.2) THEN
         tmp = 1./float(n)**3
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  vx(i,j,k) = wx(i,j,k)*tmp
                  vy(i,j,k) = wy(i,j,k)*tmp
                  vz(i,j,k) = wz(i,j,k)*tmp
               END DO
            END DO
         END DO
         DEALLOCATE ( wx,wy,wz )
      ENDIF

!
! Defines the generators

#ifdef SO3_
j1 = (/1,1,1,2,2,2,3,3,0,0,0,1,-1,0,-1,-1,1,1,1,0,0,2,1,-2,-1,0,0,-2,   &
      -1,1,1,-2,2,2,-1,1,1,-1,1,1,1,2,-2,2,2,-1,1,1,-2,2,2,1,0,0,3,1,   &
      -3,-1,0,0,-3,-1,1,1,-3,3,3,-1,1,1,-1,1,1,-1,-1,-1,-2,-2,-2,-3,    &
      -3,-0,-0,-0,-1,1,0,1,1,-1,-1,-1,0,0,-2,-1,2,1,0,0,2,1,-1,-1,2,    &
      -2,-2,1,-1,-1,1,-1,-1,-1,-2,2,-2,-2,1,-1,-1,2,-2,-2,-1,0,0,-3,    &
      -1,3,1,0,0,3,1,-1,-1,3,-3,-3,1,-1,-1,1,-1,-1/)
j2 = (/0,1,1,1,1,2,1,1,1,0,1,0,1,-1,0,1,-1,1,2,2,1,0,0,1,2,-2,-1,0,0,   &
      2,1,1,-1,1,2,-2,2,1,-1,1,2,1,2,-2,2,2,-2,2,1,-1,1,3,3,1,0,0,1,3,  &
      -3,-1,0,0,3,1,1,-1,1,3,-3,3,1,-1,1,0,-1,-1,-1,-1,-2,-1,-1,-1,0,   &
      -1,0,-1,1,0,-1,1,-1,-2,-2,-1,0,0,-1,-2,2,1,0,0,-2,-1,-1,1,-1,-2,  &
      2,-2,-1,1,-1,-2,-1,-2,2,-2,-2,2,-2,-1,1,-1,-3,-3,-1,0,0,-1,-3,3,  &
      1,0,0,-3,-1,-1,1,-1,-3,3,-3,-1,1,-1/)
j3 = (/0,0,1,0,1,1,0,1,0,1,1,1,0,1,1,1,1,-1,0,1,2,1,2,0,0,1,2,1,2,1,2,  &
      1,1,-1,1,1,-1,2,2,-2,2,2,1,1,-1,2,2,-2,2,2,-2,0,1,3,1,3,0,0,1,3,  &
      1,3,1,3,1,1,-1,1,1,-1,3,3,-3,0,0,-1,0,-1,-1,0,-1,0,-1,-1,-1,0,-1, &
      -1,-1,-1,1,0,-1,-2,-1,-2,0,0,-1,-2,-1,-2,-1,-2,-1,-1,1,-1,-1,1,   &
      -2,-2,2,-2,-2,-1,-1,1,-2,-2,2,-2,-2,2,0,-1,-3,-1,-3,0,0,-1,-3,-1, &
      -3,-1,-3,-1,-1,1,-1,-1,1,-3,-3,3/)
#endif
#ifdef SO2_
j1 = (/1,1,2,3,0,-1,1,-2,-1,1,-3,-1,0,-1,-1,-2,-3,0,1,-1,2,1,-1,3,1,0/)
j2 = (/0,1,1,1,1,1,2,1,2,3,1,3,0,0,-1,-1,-1,-1,-1,-2,-1,-2,-3,-1,-3,0/)
j3 = (/0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,-1/)
#endif

!
! Computes the structure functions

 CU : IF (curl.eq.1) THEN               ! helicity increments

      DO d = gini,ngen/2
         WRITE(ext, fmtext) d
         DO l = 1,n/2-1
            r(l) = 0.
         END DO
         DO p = 1,pmax
            DO l = 1,n/2-1
               sp(l,p) = 0.
            END DO
         END DO
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  vxl(i,j,k) = vx(i,j,k)
                  vyl(i,j,k) = vy(i,j,k)
                  vzl(i,j,k) = vz(i,j,k)
                  wxl(i,j,k) = wx(i,j,k)
                  wyl(i,j,k) = wy(i,j,k)
                  wzl(i,j,k) = wz(i,j,k)
               END DO
            END DO
         END DO
         norm = 1./SQRT(float(j1(d)**2+j2(d)**2+j3(d)**2))
         maxd = max(abs(j1(d)),abs(j2(d)))
         maxd = 2*max(maxd,abs(j3(d)))
         DO l = 1,INT(float(n-1)/maxd)
            IF ( abs(j3(d)).gt.0 ) THEN
               dis = j3(d)/abs(j3(d))
               DO i = 1,abs(j3(d))
                  CALL SHIFTZ(vxl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(vyl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(vzl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(wxl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(wyl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(wzl,dis,MPI_COMM_WORLD)
               END DO
            ENDIF
            IF ( abs(j1(d)).gt.0 ) THEN
               CALL SHIFTX(vxl,j1(d))
               CALL SHIFTX(vyl,j1(d))
               CALL SHIFTX(vzl,j1(d))
               CALL SHIFTX(wxl,j1(d))
               CALL SHIFTX(wyl,j1(d))
               CALL SHIFTX(wzl,j1(d))
            ENDIF
            IF ( abs(j2(d)).gt.0 ) THEN
               CALL SHIFTY(vxl,j2(d))
               CALL SHIFTY(vyl,j2(d))
               CALL SHIFTY(vzl,j2(d))
               CALL SHIFTY(wxl,j2(d))
               CALL SHIFTY(wyl,j2(d))
               CALL SHIFTY(wzl,j2(d))
            ENDIF
            r(l) = 2*pi*float(l)/(norm*n)
            DO p = 1,pmax
               spaux = 0
               DO k = ksta,kend
                  DO j = 1,n
                     DO i = 1,n
                        spaux = spaux+                                       &
                        abs((vxl(i,j,k)-vx(i,j,k))*(wxl(i,j,k)-wx(i,j,k))+   &
                        (vyl(i,j,k)-vy(i,j,k))*(wyl(i,j,k)-wy(i,j,k))+       &
                        (vzl(i,j,k)-vz(i,j,k))*(wzl(i,j,k)-wz(i,j,k)))**p
                     END DO
                  END DO
               END DO
               CALL MPI_REDUCE(spaux,tmp,1,MPI_REAL,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
               sp(l,p) = tmp/float(n)**3
            END DO
         END DO
         IF (myrank.eq.0) THEN
            OPEN(1,file=trim(odir) // '/hincrement.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) r
            CLOSE(1)
            OPEN(1,file=trim(odir) // '/hstructure.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) sp
            CLOSE(1)
         ENDIF
      END DO

      ELSE                              ! velocity or vorticity increments

      DO d = gini,ngen/2
         WRITE(ext, fmtext) d
         DO l = 1,n/2-1
            r(l) = 0.
         END DO
         DO p = 1,pmax
            DO l = 1,n/2-1
               sp(l,p) = 0.
            END DO
         END DO
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  vxl(i,j,k) = vx(i,j,k)
                  vyl(i,j,k) = vy(i,j,k)
                  vzl(i,j,k) = vz(i,j,k)
               END DO
            END DO
         END DO
         norm = 1./SQRT(float(j1(d)**2+j2(d)**2+j3(d)**2))
         maxd = max(abs(j1(d)),abs(j2(d)))
         maxd = 2*max(maxd,abs(j3(d)))
         DO l = 1,INT(float(n-1)/maxd)
            IF ( abs(j3(d)).gt.0 ) THEN
               dis = j3(d)/abs(j3(d))
               DO i = 1,abs(j3(d))
                  CALL SHIFTZ(vxl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(vyl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(vzl,dis,MPI_COMM_WORLD)
               END DO
            ENDIF
            IF ( abs(j1(d)).gt.0 ) THEN
               CALL SHIFTX(vxl,j1(d))
               CALL SHIFTX(vyl,j1(d))
               CALL SHIFTX(vzl,j1(d))
            ENDIF
            IF ( abs(j2(d)).gt.0 ) THEN
               CALL SHIFTY(vxl,j2(d))
               CALL SHIFTY(vyl,j2(d))
               CALL SHIFTY(vzl,j2(d))
            ENDIF
            r(l) = 2*pi*float(l)/(norm*n)
            DO p = 1,pmax
               spaux = 0
               DO k = ksta,kend
                  DO j = 1,n
                     DO i = 1,n
                        spaux = spaux+abs((((vxl(i,j,k)-vx(i,j,k))*j1(d)+ &
                         (vyl(i,j,k)-vy(i,j,k))*j2(d)+ &
                         (vzl(i,j,k)-vz(i,j,k))*j3(d))*norm)**p)
                     END DO
                  END DO
               END DO
               CALL MPI_REDUCE(spaux,tmp,1,MPI_REAL,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
               sp(l,p) = tmp/float(n)**3
            END DO
         END DO
         IF (myrank.eq.0) THEN
            OPEN(1,file=trim(odir) // '/increment.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) r
            CLOSE(1)
            OPEN(1,file=trim(odir) // '/structure.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) sp
            CLOSE(1)
         ENDIF
      END DO

      ENDIF CU

!
! End of STRUCTURE3D

      CALL MPI_FINALIZE(ierr)
      DEALLOCATE( vx,vxl )
      IF ((curl.eq.0).or.(curl.eq.2)) THEN
         DEALLOCATE ( vy,vz,vyl,vzl )
      ENDIF
      DEALLOCATE( sp,r )

      END PROGRAM STRUCTURE3D
