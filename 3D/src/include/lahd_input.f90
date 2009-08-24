! Extra input files for the LAHD solver

! Reads from the external file 'hall.txt' 
! parameters for Hall-MHD runs
!     ep  : amplitude of the Hall effect
!     gspe: = 0 skips generalized helicity spectrum computation
!           = 1 computes the spectrum of generalized helicity

!
! Reads from the external file 'alpha.txt'
! the value of kinetic alpha
!     alpk: kinetic alpha

      IF (myrank.eq.0) THEN
         OPEN(1,file='alpha.txt',status='unknown')
         READ(1,*) alpk
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(alpk,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
