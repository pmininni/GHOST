! External mechanical forcing.
! This file contains the expression used for the external
! mechanical forcing. You can use temporary real arrays
! R1-R3 of size (1:nx,1:ny,ksta:kend) and temporary complex
! arrays C1-C8 of size (nz,ny,ista:iend) to do intermediate
! computations. The variable f0 should control the global
! amplitude of the forcing, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays fx, fy, and fz.

! Reads next nudging field
!     kdn : minimum wave number
!     kup : maximum wave number
!     fstep : interpolation time (invoked in main when using rand=2)

      ! Get next field
      IF (ndginitialcall.eq.1) THEN
         ndginitialcall = 0
      ELSE
         WRITE(ndgext, ndgfmtext) (t-1)+fstep
         CALL io_read(1,idir,'ndgx',ndgext,planio,R1)
         CALL io_read(1,idir,'ndgy',ndgext,planio,R2)
         CALL io_read(1,idir,'ndgz',ndgext,planio,R3)
         CALL fftp3d_real_to_complex(planrc,R1,fxold,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R2,fyold,MPI_COMM_WORLD)
         CALL fftp3d_real_to_complex(planrc,R3,fzold,MPI_COMM_WORLD)
      ENDIF

      ! Change flag
      ndgfilesloaded = 1
