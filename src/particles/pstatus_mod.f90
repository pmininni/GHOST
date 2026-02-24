!=================================================================
! General MODULES for particle methods
!
! 2026 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

  MODULE pstatus
!
! General status flags for all particle integrators
! Reads a "&pstatus" namelist with:    
!  maxparts    : Maximum number of particles
!  lgmult      : Multiplier for part output (must divide tstep evenly)
!  ilgintrptype: Interpolation type: only GPINTRP_CSPLINE currently
!  ilgexchtype : Boundary exchange type: GPEXCHTYPE_VDB (voxel db)
!                or GPEXCHTYPE_NN (nearest-neighbor)
!  ilgouttype  : Particle output type: 0=binary; 1=ASCII
!  ilgwrtunit  : Units for part position write: 0=box units; 1=grid units
!  lgseedfile  : Name of seed file if using user initialization
!  ilgcoll     : 1=binary collective I/O; 0=task 0 binary (posix) I/O
!  ilgwrtunit  : write particle positions in box units (==1) (i.e.,
!            x,y,z in [0,2.pi]), or in grid units (==0) (x,y,z in [0,N]).
      INTEGER             :: pstep, tind, timep
      INTEGER             :: maxparts
      INTEGER             :: lgmult
      INTEGER             :: ilgintrptype, ilgexchtype
      INTEGER             :: ilgouttype, ilgwrtunit
      INTEGER             :: ilgcoll, ilgfpfiletype
      CHARACTER(len=1024) :: slgfpfile, lgseedfile
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE pstatus_init(infile_)
!
! Initializes general global status parameters for the particle
! solvers
      USE particlebase_mod
      USE commtypes
      USE mpivars
      USE status
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: infile_
      NAMELIST / pstatus / maxparts,lgmult,ilgintrptype
      NAMELIST / pstatus / ilgexchtype,ilgouttype,ilgwrtunit
      NAMELIST / pstatus / ilgcoll,ilgfpfiletype,slgfpfile,lgseedfile

      maxparts     = 1000
      ilginittype  = GPINIT_RANDLOC
      ilgintrptype = GPINTRP_CSPLINE
      ilgexchtype  = GPEXCHTYPE_VDB
      ilgouttype   = 0
      ilgwrtunit   = 0
      lgmult       = 1
      lgseedfile   = 'gplag.dat'
      ilgcoll      = 1
      rbal         = 0.0
      nwpart       = 0
      ilgfpfiletype= 0
      slgfpfile    = 'xlgInitRndSeed.000.txt'
      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=pstatus)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(maxparts     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lgmult       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgintrptype ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgexchtype  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgouttype   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgwrtunit   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgcoll      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ilgfpfiletype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(slgfpfile ,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(lgseedfile,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      IF ( mod(tstep,lgmult).NE.0 ) THEN
        WRITE(*,*)'main: lgmult must divide tstep evenly'
        STOP
      ENDIF
      pstep = tstep/lgmult
      END SUBROUTINE pstatus_init

  END MODULE pstatus
