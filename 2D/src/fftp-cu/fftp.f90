!=================================================================
! FFTP
! Parallel Fast Fourier Transform in 2D and 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the CUDA library in each node. This file contains 
! subroutines common to the 2D and 3D versions of FFTP.
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE range(n1,n2,nprocs,irank,ista,iend)
!-----------------------------------------------------------------
!
! Soubroutine for computing the local coordinate range 
! when splitting the original array into the nodes
!
! Parameters
!     n1     : the minimum value in the splitted dimension [IN]
!     n2     : the maximum value in the splitted dimension [IN]
!     nprocs : the number of processors [IN]
!     irank  : the rank of the processor [IN]
!     ista   : start value for the local coordinate [OUT]
!     iend   : end value for the local coordinate [OUT]
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n1,n2
      INTEGER, INTENT(IN)  :: nprocs,irank
      INTEGER, INTENT(OUT) :: ista,iend

      INTEGER :: myrank
      INTEGER :: iwork1,iwork2

      iwork1 = (n2-n1+1)/nprocs
      iwork2 = MOD(n2-n1+1,nprocs)
      ista = irank*iwork1+n1+MIN(irank,iwork2)
      iend = ista+iwork1-1
      IF (iwork2.gt.irank) iend = iend+1

      RETURN
      END SUBROUTINE range
