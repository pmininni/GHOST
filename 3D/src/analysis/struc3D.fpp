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
!           HD_SOL    Hydrodynamic SO(3) structure functions
!           ROTH_SOL  SO(2) HD structure functions (rotation)
!           PHD_SOL   Hydrodynamic and passive scalar with SO(3)
!           PROTH_SOL Hydrodynamic and passive scalar with SO(2)
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
!  2 Jul 2008: Main program for all decompositions and solvers
!  4 Feb 2009: Option for increments of helicity and vorticity
! 30 Nov 2009: Option for passive scalar (with P. Rodriguez Imazio)
!=================================================================

!
! Definitions for conditional compilation

#ifdef HD_SOL
#define SO3_
#endif

#ifdef ROTH_SOL
#define SO2_
#endif

#ifdef PHD_SOL
#define SO3_
#define SCALAR_
#endif

#ifdef PROTH_SOL
#define SO2_
#define SCALAR_
#endif

#ifdef ROTBOUSS_SOL
#define SO2_
#define SCALAR_
#endif

#ifdef BOUSS_SOL
#define SO2_
#define SCALAR_
#endif


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
      USE threads
      IMPLICIT NONE

!
! Arrays for the fields and structure functions

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C1,C2,C3
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: C4,C5

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vx,vy,vz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: vxl,vyl,vzl
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: wx,wy,wz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: wxl,wyl,wzl
#ifdef SCALAR_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:,:) :: th,thl
#endif

      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:)   :: sp
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: r
#ifdef SCALAR_
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:,:)   :: zp
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:)     :: flux
#endif

!
! Auxiliary variables

      REAL(KIND=GP)    :: tmp,norm,spaux

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
!     curl : =0 computes increments of the field (and the scalar)
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

      ALLOCATE( vx(n,n,ksta:kend) )
      ALLOCATE( vy(n,n,ksta:kend) )
      ALLOCATE( vz(n,n,ksta:kend) )
#ifdef SCALAR_
      ALLOCATE( th(n,n,ksta:kend) )
      ALLOCATE( zp(n/2-1,pmax) )
      ALLOCATE( flux(n/2-1) )
#endif
      ALLOCATE( sp(n/2-1,pmax) )
      ALLOCATE( r(n/2-1) )

!
! Reads the external binary files with the 
! fields components, and computes auxiliary
! fields if needed.
! If required, initializes the FFT library and 
! arrays with wavenumbers. Use FFTW_ESTIMATE in 
! short runs and FFTW_MEASURE in long runs.

      WRITE(ext, fmtext) stat
if ( myrank.eq. 0 ) write(*,*)'main: opening vx stat=',stat
      CALL io_read(1,idir,'vx',ext,planio,vx)
if ( myrank.eq. 0 ) write(*,*)'main: opening vy stat=',stat
      CALL io_read(1,idir,'vy',ext,planio,vy)
if ( myrank.eq. 0 ) write(*,*)'main: opening vz stat=',stat
      CALL io_read(1,idir,'vz',ext,planio,vz)
#ifdef SCALAR_
if ( myrank.eq. 0 ) write(*,*)'main: opening th stat=',stat
      CALL io_read(1,idir,'th',ext,planio,th)
#endif
if ( myrank.eq. 0 ) write(*,*)'main: stat=',stat,' done.'
      IF (curl.gt.0) THEN
         ALLOCATE( C1(n,n,ista:iend), C2(n,n,ista:iend) )
         ALLOCATE( C3(n,n,ista:iend), C4(n,n,ista:iend) )
         ALLOCATE( C5(n,n,ista:iend) )
         ALLOCATE( wx(n,n,ksta:kend) )
         ALLOCATE( wy(n,n,ksta:kend) )
         ALLOCATE( wz(n,n,ksta:kend) )
         CALL fftp3d_create_plan(planrc,n,FFTW_REAL_TO_COMPLEX, &
             FFTW_ESTIMATE)
         CALL fftp3d_create_plan(plancr,n,FFTW_COMPLEX_TO_REAL, &
             FFTW_ESTIMATE)
         ALLOCATE( ka(n), ka2(n,n,ista:iend) )
         kmax = (real(n,kind=GP)/3.)**2
         tiny = 1e-5
         DO i = 1,n/2
            ka(i) = real(i-1,kind=GP)
            ka(i+n/2) = real(i-n/2-1,kind=GP)
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
      IF (curl.eq.1) THEN
         tmp = 1./real(n,kind=GP)**3
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
         tmp = 1./real(n,kind=GP)**3
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
      ALLOCATE ( vxl(n,n,ksta:kend) )
      ALLOCATE ( vyl(n,n,ksta:kend) )
      ALLOCATE ( vzl(n,n,ksta:kend) )
#ifdef SCALAR_
      ALLOCATE( thl(n,n,ksta:kend) )
#endif

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
         norm = 1./SQRT(real(j1(d)**2+j2(d)**2+j3(d)**2,kind=GP))
         maxd = max(abs(j1(d)),abs(j2(d)))
         maxd = 2*max(maxd,abs(j3(d)))
         DO l = 1,INT(real(n-1,kind=GP)/maxd)
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
            r(l) = 2*pi*real(l,kind=GP)/(norm*n)
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
               CALL MPI_REDUCE(spaux,tmp,1,GC_REAL,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
               sp(l,p) = tmp/real(n,kind=GP)**3
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
if ( myrank.eq. 0 ) write(*,*)'main: doing velocity increments... d=',d
         WRITE(ext, fmtext) d
         DO l = 1,n/2-1
            r(l) = 0.
#ifdef SCALAR_
            flux(l) = 0.
#endif
         END DO
         DO p = 1,pmax
            DO l = 1,n/2-1
               sp(l,p) = 0.
#ifdef SCALAR_
               zp(l,p) = 0.
#endif
            END DO
         END DO
if ( myrank.eq. 0 ) write(*,*)'main: setting vxl...'
         DO k = ksta,kend
            DO j = 1,n
               DO i = 1,n
                  vxl(i,j,k) = vx(i,j,k)
                  vyl(i,j,k) = vy(i,j,k)
                  vzl(i,j,k) = vz(i,j,k)
#ifdef SCALAR_
                  thl(i,j,k) = th(i,j,k)
#endif
               END DO
            END DO
         END DO
if ( myrank.eq. 0 ) write(*,*)'main: setting of vxl,, thl, ... done.'
         norm = 1./SQRT(real(j1(d)**2+j2(d)**2+j3(d)**2,kind=GP))
         maxd = max(abs(j1(d)),abs(j2(d)))
         maxd = 2*max(maxd,abs(j3(d)))
         DO l = 1,INT(real(n-1,kind=GP)/maxd)
if ( myrank.eq. 0 ) write(*,*)'main: shifting loop, l=',l,' maxd=',maxd
            IF ( abs(j3(d)).gt.0 ) THEN
               dis = j3(d)/abs(j3(d))
               DO i = 1,abs(j3(d))
                  CALL SHIFTZ(vxl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(vyl,dis,MPI_COMM_WORLD)
                  CALL SHIFTZ(vzl,dis,MPI_COMM_WORLD)
#ifdef SCALAR_
                  CALL SHIFTZ(thl,dis,MPI_COMM_WORLD)
#endif
               END DO
            ENDIF
            IF ( abs(j1(d)).gt.0 ) THEN
               CALL SHIFTX(vxl,j1(d))
               CALL SHIFTX(vyl,j1(d))
               CALL SHIFTX(vzl,j1(d))
#ifdef SCALAR_
               CALL SHIFTX(thl,j1(d))
#endif
            ENDIF
            IF ( abs(j2(d)).gt.0 ) THEN
               CALL SHIFTY(vxl,j2(d))
               CALL SHIFTY(vyl,j2(d))
               CALL SHIFTY(vzl,j2(d))
#ifdef SCALAR_
               CALL SHIFTY(thl,j2(d))
#endif
            ENDIF
if ( myrank.eq. 0 ) write(*,*)'main: shifting loop, l=',l,' done.'
            r(l) = 2*pi*real(l,kind=GP)/(norm*n)
 LP :       DO p = 1,pmax
               spaux = 0.               ! vector field increments
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
               DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                  DO j = 1,n
                     DO i = 1,n
                        spaux = spaux+abs((((vxl(i,j,k)-vx(i,j,k))*j1(d)+ &
                         (vyl(i,j,k)-vy(i,j,k))*j2(d)+ &
                         (vzl(i,j,k)-vz(i,j,k))*j3(d))*norm)**p)
                     END DO
                  END DO
               END DO
               CALL MPI_REDUCE(spaux,tmp,1,GC_REAL,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
               sp(l,p) = tmp/real(n,kind=GP)**3
if ( myrank.eq. 0 ) write(*,*)'main: velocity increments done.'
#ifdef SCALAR_
               IF (curl.eq.0) THEN      ! passive scalar increments
               spaux = 0.
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
               DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                  DO j = 1,n
                     DO i = 1,n
                        spaux = spaux+abs(thl(i,j,k)-th(i,j,k))**p
                     END DO
                  END DO
               END DO
               CALL MPI_REDUCE(spaux,tmp,1,GC_REAL,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
               zp(l,p) = tmp/real(n,kind=GP)**3
               ENDIF
if ( myrank.eq. 0 ) write(*,*)'main: scalar increments done.'
#endif
if ( myrank.eq. 0 ) write(*,*)'main: increments p=', p, ' done.'
            END DO LP
#ifdef SCALAR_
            IF (curl.eq.0) THEN         ! passive scalar theorem
               spaux = 0.
!$omp parallel do if (kend-ksta.ge.nth) private (j,i)
               DO k = ksta,kend
!$omp parallel do if (kend-ksta.lt.nth) private (i)
                  DO j = 1,n
                     DO i = 1,n
                        spaux = spaux+abs(((vxl(i,j,k)-vx(i,j,k))*j1(d)+ &
                         (vyl(i,j,k)-vy(i,j,k))*j2(d)+ &
                         (vzl(i,j,k)-vz(i,j,k))*j3(d))*norm)* &
                         (thl(i,j,k)-th(i,j,k))**2
                     END DO
                  END DO
               END DO
               CALL MPI_REDUCE(spaux,tmp,1,GC_REAL,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
               flux(l) = tmp/real(n,kind=GP)**3
            ENDIF
#endif
         END DO
         IF (myrank.eq.0) THEN
if ( myrank.eq. 0 ) write(*,*)'main: writing increment: ', ext

            OPEN(1,file=trim(odir) // '/increment.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) r
            CLOSE(1)
            OPEN(1,file=trim(odir) // '/structure.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) sp
            CLOSE(1)
         ENDIF
#ifdef SCALAR_
         IF ((myrank.eq.0).and.(curl.eq.0)) THEN
            OPEN(1,file=trim(odir) // '/scalarstr.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) zp
            CLOSE(1)
            OPEN(1,file=trim(odir) // '/sctheorem.' &
              // ext // '.out' ,form='unformatted')
            WRITE(1) flux
            CLOSE(1)
         ENDIF
#endif         
if ( myrank.eq. 0 ) write(*,*)'main: increment ', ext, ' written.'
      END DO

      ENDIF CU

!
! End of STRUCTURE3D

      CALL MPI_FINALIZE(ierr)
      IF ((curl.eq.0).or.(curl.eq.2)) THEN
         DEALLOCATE ( vy,vz,vyl,vzl )
      ENDIF
      DEALLOCATE( vx,vxl )
      DEALLOCATE( sp,r )
#ifdef SCALAR_
      DEALLOCATE( th,thl )
      DEALLOCATE( zp,flux )
#endif         

      END PROGRAM STRUCTURE3D
