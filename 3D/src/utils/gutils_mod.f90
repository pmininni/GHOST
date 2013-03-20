!=================================================================
! GHOST suite Fortran utilities 
!
! 2010 D. Rosenberg
!      NCAR
!
! 17 Nov 2010: Initial version
!=================================================================
MODULE gutils
! Utils data, if any:
!
! ...
! end, member data
      REAL, DIMENSION(:), ALLOCATABLE   :: fpdf_
      REAL, DIMENSION(:), ALLOCATABLE   :: gpdf_
      INTEGER                           :: nbins_=50
!
!
! Methods:
      CONTAINS

      INTEGER FUNCTION imaxarrp(Iin, nin) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Performax max function over integer array for positive-def array
!
! Parameters
!     Iin : integer array to find max of
!     nin : dimension of Iin
!-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN),DIMENSION(:) :: Iin
      INTEGER,INTENT(IN)              :: nin
      INTEGER                         :: i

      imaxarrp = 0
      DO i = 1, nin
        imaxarrp = max(imaxarrp,Iin(i))
      ENDDO

      END FUNCTION imaxarrp
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      INTEGER FUNCTION iavgarrp(Iin, nin) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Performax averaging function over integer array for positive-def array
!
! Parameters
!     Iin : integer array to find max of
!     nin : dimension of Iin
!-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN),DIMENSION(:) :: Iin
      INTEGER,INTENT(IN)              :: nin
      INTEGER                         :: i

      iavgarrp = 0
      DO i = 1, nin
        iavgarrp = iavgarrp + Iin(i)
      ENDDO
      iavgarrp = int(iavgarrp/nin)

      END FUNCTION iavgarrp
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      INTEGER FUNCTION isumarrp(Iin, nin)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Performa  sum function over integer array for positive-def array
!
! Parameters
!     Iin : integer array to find max of
!     nin : dimension of Iin
!-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN),DIMENSION(:) :: Iin
      INTEGER,INTENT(IN)              :: nin
      INTEGER                         :: i

      isumarrp = 0
      DO i = 1, nin
        isumarrp = isumarrp + Iin(i)
      ENDDO

      END FUNCTION isumarrp
!-----------------------------------------------------------------
!-----------------------------------------------------------------


      SUBROUTINE rarray_byte_swap(Rin, nin)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Performs endian conversion on array of floats
!
! Parameters
!     Rin : artibrary-rank array whose values will be endian-swapped, returned
!     nin : dimension of Rin
!-----------------------------------------------------------------
      USE mpivars
      USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(IN)                        :: nin
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(*) :: Rin

      INTEGER(KIND=GP) :: ie0,ie1
      INTEGER          :: i, j, k, m, nb

      nb = 8  ! no. bits per byte

      ie1 = 0
      DO k = 1, nin
          ie0 = TRANSFER(Rin(k), 0_GP)
          DO m = 1, GP
             CALL MVBITS( ie0, (GP-m)*nb, nb, ie1, (m-1)*nb  )
          END DO
          Rin(k) = TRANSFER(ie1, 0.0_GP)
       END DO

      RETURN

      END SUBROUTINE rarray_byte_swap
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!
!
      SUBROUTINE parseind(sind, sep, ind, nmax, nind) 
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Parses string of integers that are ';' separated, and stores, and
! stores them in the integer array, specifying how many integers
! were found in the string. A test is made against the specified
! max number of integer indices, nmax, before attempting to add
! any more indices to this array.
!
! Parameters
!     sind : ';'-separated string of integers (IN)
!     sep  : separator string (IN)
!     ind  : integer array contining integers found in sind (IN)
!     nmax : max size of 'ind' array (OUT)
!     nind : num. integers found in 'sind' (OUT)
!-----------------------------------------------------------------
      INTEGER, INTENT(OUT)         :: ind(*), nind
      INTEGER, INTENT(IN)          :: nmax
      CHARACTER(len=*), INTENT(IN) :: sind, sep

      INTEGER                      :: i
      CHARACTER(len=1024)          :: sint

      ib = 1;
      ie = len(sind)
      nind = 0
      DO WHILE ( len(trim(sind(ib:ie))) .GT. 0 )
        i = index(sind(ib:ie),sep)
        IF ( i .eq. 0 ) THEN
          sint = trim(adjustl(sind(ib:ie)))
          ib = ie + 1
        ELSE
          sint = trim(adjustl(sind(ib:(ib+i-2))))
          ib = ib + i
        ENDIF
        nind = nind + 1
        IF ( nind.GT.nmax ) RETURN
        READ(sint,'(I10)') ind(nind)
      ENDDO
      RETURN

      END SUBROUTINE parseind
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!
!
      SUBROUTINE dopdfr(Rin, nin, spref, nmb, nbins, ifixdr, fmin, fmax, dolog)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes pdf of real quantity Rin, and outputs it to a file whose
! prefix is specified.
!
! Parameters
!     Rin   : real array (kind GP) of dimension nin
!     nin   : integer size of array
!     spref : pdf output file name prefix
!     nmb   : output interval id extension
!     nbins : number of bins to use for PDF (>0)
!     ifixdr: if > 0, will use dr0(1) as bounds for dynamic
!             range; else it will compute them dynamically
!     fmin  : highest dynamic range value, returned if ifixdr>0
!     fmax  : lowest dynamic range value, returned if ifixdr>0
!     dolog : take log of quantity magnitude when creating bins (0, 1)
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(IN), DIMENSION(*)    :: Rin
      REAL(KIND=GP), INTENT(INOUT)               :: fmin,fmax
      INTEGER, INTENT(IN)                        :: nin,nbins,ifixdr,dolog
      CHARACTER(len=*), INTENT(IN)               :: spref, nmb

      REAL(KIND=GP)                              :: del,gmin,gmax
      REAL(KIND=GP)                              :: smin,smax,sr
      INTEGER                                    :: i,ibin
      CHARACTER(len=1024)                        :: shead

      IF ( nbins.GT.nbins_       .OR. &
          .NOT. ALLOCATED(fpdf_) .OR. &
          .NOT. ALLOCATED(gpdf_)      &
          ) THEN
        ! Re-allocate if necessary:
        IF ( ALLOCATED(fpdf_) ) DEALLOCATE(fpdf_)
        IF ( ALLOCATED(gpdf_) ) DEALLOCATE(gpdf_)
        ALLOCATE(fpdf_(nbins))
        ALLOCATE(gpdf_(nbins))
        nbins_ = nbins
      ENDIF
      fpdf_ = 0.0
      gpdf_ = 0.0

      IF ( ifixdr .le. 0 ) THEN
        ! Compute dynamic range of PDF
        fmin = MINVAL(Rin(1:nin))
        fmax = MAXVAL(Rin(1:nin))
        CALL MPI_ALLREDUCE(fmin,gmin,1, GC_REAL,      &
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(fmax,gmax,1, GC_REAL,      &
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        fmin = gmin
        fmax = gmax
      ENDIF
  
      ! Compute local PDF:

      IF ( dolog .GT. 0 ) THEN
        smin = sign(1.0,fmin)
        smax = sign(1.0,fmax)
        del = ABS(log10(fmax)-log10(fmin+TINY(1.0_GP)))/dble(nbins)
        DO i = 1, nin
          sr   = sign(1.0,Rin(i))
          ibin = NINT( ( sr*log10(abs(Rin(i))) - smin*log10(abs(fmin+tiny(1.0_GP))) )/del )
          ibin = MIN(MAX(ibin,1),nbins_)
          fpdf_(ibin) = fpdf_(ibin) + 1.0
        ENDDO
      ELSE
        del = (fmax-fmin)/dble(nbins)
        DO i = 1, nin
          ibin = NINT( (Rin(i) - fmin)/del )
          ibin = MIN(MAX(ibin,1),nbins_)
          fpdf_(ibin) = fpdf_(ibin) + 1.0
        ENDDO
      ENDIF

      ! Compute global reduction between MPI tasks:
      CALL MPI_REDUCE(fpdf_, gpdf_, nbins_, GC_REAL, &
                      MPI_SUM, 0, MPI_COMM_WORLD,ierr)

      ! Write PDF to disk:
      IF ( myrank.eq.0 ) THEN
         WRITE(shead,'(A9,E23.15,A1,E23.15,A8,I4,A7,I1)') '# range=[', fmin, ',' , fmax, &
           ']; nbin=', nbins, '; blog=', dolog   
         OPEN(1,file=trim(spref)// '.' // nmb // '.txt')
         WRITE(1,'(A)') trim(shead)
         WRITE(1,40) gpdf_
   40    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF
 
      RETURN

      END SUBROUTINE dopdfr
!-----------------------------------------------------------------
!-----------------------------------------------------------------


END MODULE gutils
