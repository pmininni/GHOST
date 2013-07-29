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
      REAL, DIMENSION  (:), ALLOCATABLE   :: fpdf_
      REAL, DIMENSION  (:), ALLOCATABLE   :: gpdf_
      REAL, DIMENSION(:,:), ALLOCATABLE   :: fpdf2_
      REAL, DIMENSION(:,:), ALLOCATABLE   :: gpdf2_
      INTEGER                             :: nbins_=50,nbins2_(2)=(/50,50/)
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

      SUBROUTINE parsestr(strin, sep, astr, astrlen, nmax, nind)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Parses string of names, separated by 'sep' delimiter, and
! returns them in an array.
!
! Parameters
!     strin: 'sep'-separated string of strings(IN)
!     sep  : separator string (IN)
!     astr : array containing strings found in 'strin' (IN)
!     nmax : max size of 'astr' array (OUT)
!     nstr : num. strings found in 'strin' (OUT)
!-----------------------------------------------------------------
      INTEGER               , INTENT(IN)  :: astrlen,nmax
      CHARACTER(len=*)      , INTENT (IN) :: strin,sep
      CHARACTER(len=astrlen), INTENT(OUT) :: astr(nmax)
      CHARACTER(len=astrlen)              :: tstr
      INTEGER                             :: i, ib, ie

      ib = 1;
      ie = len(strin)
      nind = 0
      DO WHILE ( len(trim(strin(ib:ie))) .GT. 0 )
        i = index(strin(ib:ie),sep)
        IF ( i .eq. 0 ) THEN
          astr(nind+1) = trim(adjustl(strin(ib:ie)))
          ib = ie + 1
        ELSE
          astr(nind+1) = trim(adjustl(strin(ib:(ib+i-2))))
          ib = ib + i
        ENDIF
        nind = nind + 1
        IF ( nind.GE.nmax ) RETURN
      ENDDO
      RETURN

      END SUBROUTINE parsestr
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE dopdfr(Rin, nin, n, fname, nbins, ifixdr, fmin, fmax, dolog)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes pdf of real quantity Rin, and outputs it to a file whose
! prefix is specified.
!
! Parameters
!     Rin   : real array (kind GP) of dimension nin
!     nin   : integer size of array
!     n     : full linear domain size
!     fname : full file name
!     nbins : number of bins to use for PDF (>0)
!     ifixdr: if > 0, will use fmin/fmax as bounds for dynamic
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
      INTEGER, INTENT(IN)                        :: nin,n,nbins,ifixdr,dolog
      CHARACTER(len=*), INTENT(IN)               :: fname

      REAL(KIND=GP)                              :: del,gmin,gmax
      REAL(KIND=GP)                              :: gavg,sumr
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

      DO i = 1, nbins
       fpdf_(i) = 0.0_GP
       gpdf_(i) = 0.0_GP
      ENDDO

      IF ( ifixdr .le. 0 ) THEN
        ! Compute dynamic range of PDF
        fmin = MINVAL(Rin(1:nin))
        fmax = MAXVAL(Rin(1:nin))
        IF ( dolog .GT. 0 ) THEN
          fmin = log10(fmin+tiny(1.0_GP))
          fmax = log10(fmax+tiny(1.0_GP))
        ENDIF
        CALL MPI_ALLREDUCE(fmin,gmin,1, GC_REAL,      &
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(fmax,gmax,1, GC_REAL,      &
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        fmin = gmin
        fmax = gmax
      ENDIF
  
      sumr = 0.0_GP
!$omp parallel do 
      DO i = 1, nin
!$omp atomic
        sumr  = sumr +  Rin(i)
      ENDDO

      CALL MPI_ALLREDUCE(sumr, gavg, 1, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      gavg = gavg/real(n,kind=GP)**3

      del = ABS(fmax-fmin)/dble(nbins)
      ! Compute local PDF:
      IF ( dolog .GT. 0 ) THEN
!$omp parallel do private (ibin)
        DO i = 1, nin
          ibin = NINT( ( log10(abs(Rin(i)+tiny(1.0_GP))) - fmin )/del )
          ibin = MIN(MAX(ibin,1),nbins_)
!$omp atomic
          fpdf_(ibin) = fpdf_(ibin) + 1.0_GP
        ENDDO
      ELSE
!$omp parallel do private (ibin)
        DO i = 1, nin
          ibin = NINT( (Rin(i) - fmin)/del )
          ibin = MIN(MAX(ibin,1),nbins_)
!$omp atomic
          fpdf_(ibin) = fpdf_(ibin) + 1.0_GP
        ENDDO
      ENDIF

      ! Compute global reduction between MPI tasks:
      CALL MPI_REDUCE(fpdf_, gpdf_, nbins_, MPI_FLOAT, &
                      MPI_SUM, 0, MPI_COMM_WORLD,ierr)
!     CALL MPI_ALLREDUCE(fpdf_, gpdf_, nbins_, MPI_FLOAT, &
!                     MPI_SUM, MPI_COMM_WORLD,ierr)

      
      ibin = 0
      DO i = 1, nbins
        ibin = ibin + fpdf_(i) 
      ENDDO
      if (ibin.ne.nin) then
        write(*,*)'dopdfr: inconsistent data: expected: ',nin, ' found: ',ibin
        stop
      endif

      ! Write PDF to disk:
      IF ( myrank.eq.0 ) THEN
          IF ( dolog .GT. 0 ) THEN
            fmin = 10.0_GP**fmin
            fmax = 10.0_GP**fmax
          ENDIF
         WRITE(shead,'(A9,E16.8,A1,E16.8,A7,E16.8,A7,I4,A7,I1)') '# range=[', fmin, ',' , fmax, &
           ']; avg=',gavg,'; nbin=', nbins, '; blog=', dolog   
         OPEN(1,file=trim(fname))
         WRITE(1,'(A)') trim(shead)
         WRITE(1,40) gpdf_
   40    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF
 
      RETURN

      END SUBROUTINE dopdfr
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
!
      SUBROUTINE dojpdfr(R1, sR1, R2, sR2, nin, n, fname, nbins, ifixdr, fmin, fmax, dolog)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes joint pdf of quantities R1, R2, and outputs it to a file whose
! prefix is specified.
!
! Parameters
!     R1/R2  : real arrays (kind GP) of dimension nin
!     sR1/sR2: descriptions of quantities R1 & R2; only 4 chars are used
!     nin    : integer size of arrays R2/R2
!     n     : full linear domain size
!     fname  : output interval id extension
!     nbins  : 2D array giving number of bins to use for PDF (>0) 
!              in each 'direction'
!     ifixdr : if > 0, will use dr0(1) as bounds for dynamic
!              range; else it will compute them dynamically
!     fmin   : highest dynamic range value, returned if ifixdr>0
!     fmax   : lowest dynamic range value, returned if ifixdr>0
!     dolog  : 2-element vectory indicating whether take log of quantity magnitude 
!              in 1-, 2- directions when creating bins (0, 1)
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      REAL(KIND=GP)   , INTENT   (IN), DIMENSION(nin):: R1,R2
      REAL(KIND=GP)   , INTENT(INOUT)                :: fmin(2),fmax(2)
      INTEGER         , INTENT   (IN)                :: nin,n,nbins(2),ifixdr,dolog(2)
      CHARACTER(len=*), INTENT   (IN)                :: sR1,sR2,fname

      REAL(KIND=GP)                                  :: del (2),fck,gmin(2),gmax(2)
      REAL(KIND=GP)                                  :: aa,gavg(2),sumr(2)
      INTEGER                                        :: i,j,jx,jy
      CHARACTER(len=2048)                            :: shead

      IF ( nbins(1).GT.nbins2_(1) .OR. &
           nbins(2).GT.nbins2_(2) .OR. &
          .NOT. ALLOCATED(fpdf2_) .OR. &
          .NOT. ALLOCATED(gpdf2_)      ) THEN

        ! Re-allocate if necessary:
        IF ( ALLOCATED(fpdf2_) ) DEALLOCATE(fpdf2_)
        IF ( ALLOCATED(gpdf2_) ) DEALLOCATE(gpdf2_)
        ALLOCATE(fpdf2_(nbins(1),nbins(2)))
        ALLOCATE(gpdf2_(nbins(1),nbins(2)))
        nbins2_(1:2) = nbins(1:2)
      ENDIF
      fpdf2_(1:nbins(1),1:nbins(2)) = 0.0_GP
      gpdf2_(1:nbins(1),1:nbins(2)) = 0.0_GP

      IF ( ifixdr .le. 0 ) THEN
        ! Compute dynamic range of PDF
        fmin(1) = MINVAL(R1(1:nin))
        fmax(1) = MAXVAL(R1(1:nin))
        fmin(2) = MINVAL(R2(1:nin))
        fmax(2) = MAXVAL(R2(1:nin))
        DO j = 1, 2
          IF ( dolog(j) .GT. 0 ) fmin(j) = log10(fmin(j)+tiny(1.0_GP))
          IF ( dolog(j) .GT. 0 ) fmax(j) = log10(fmax(j)+tiny(1.0_GP))
        ENDDO
        CALL MPI_ALLREDUCE(fmin,gmin,2, GC_REAL,      &
                           MPI_MIN,MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(fmax,gmax,2, GC_REAL,      &
                           MPI_MAX,MPI_COMM_WORLD,ierr)
        fmin(1:2) = gmin(1:2)
        fmax(1:2) = gmax(1:2)
      ENDIF
  

      aa = 0.0_GP
!$omp parallel do reduction(+:aa)
      DO i = 1, nin
        aa  = aa +  R1(i)
      ENDDO
      sumr(1) = aa

      aa = 0.0_GP
!$omp parallel do reduction(+:aa)
      DO i = 1, nin
        aa  = aa +  R2(i)
      ENDDO
      sumr(2) = aa

      CALL MPI_ALLREDUCE(sumr, gavg, 2, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      gavg = gavg/real(n,kind=GP)**3

      ! Compute local PDF:
      DO j = 1, 2
        del (j) = ABS(fmax(j)-fmin(j))/dble(nbins(j))
      ENDDO

      ! Compute local PDF:
      IF        ( dolog(1).GT.0 .AND. dolog(2).GT.0 ) THEN
!$omp parallel do private (jx,jy) 
        DO i = 1, nin
          jx  = NINT( ( log10(abs(R1(i))+tiny(1.0_GP)) - fmin(1) )/del(1) )
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( ( log10(abs(R2(i))+tiny(1.0_GP)) - fmin(2) )/del(2) )
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
        ENDDO
      ELSE IF ( dolog(1).GT.0 .AND. dolog(2).LE.0 ) THEN
!$omp parallel do private (jx,jy) 
        DO i = 1, nin
          jx  = NINT( ( log10(abs(R1(i))+tiny(1.0_GP)) - fmin(1) )/del(1) )
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( (R2(i) - fmin(2))/del(2))
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
        ENDDO
      ELSE IF ( dolog(1).LE.0 .AND. dolog(2).GT.0 ) THEN
!$omp parallel do private (jx,jy) 
        DO i = 1, nin
          jx  = NINT( (R1(i) - fmin(1))/del(1)) 
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( ( log10(abs(R2(i))+tiny(1.0_GP)) - fmin(2) )/del(2) )
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
        ENDDO
      ELSE IF ( dolog(1).LE.0 .AND. dolog(2).LE.0 ) THEN
!$omp parallel do private (jx,jy) 
        DO i = 1, nin
          jx  = NINT( (R1(i) - fmin(1))/del(1)) 
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( (R2(i) - fmin(2))/del(2))
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
        ENDDO
      ENDIF

      ! Compute global reduction between MPI tasks:
      CALL MPI_REDUCE(fpdf2_, gpdf2_, nbins(1)*nbins(2), MPI_FLOAT, &
                      MPI_SUM, 0, MPI_COMM_WORLD,ierr)

      ! Write PDF to disk:
      IF ( myrank.eq.0 ) THEN
! Check:
        fck = 0.0_GP
        DO i = 1, nbins(1)
          DO j = 1, nbins(2)
            fck = fck + fpdf2_(i,j)
          ENDDO
        ENDDO
        IF ( fck .ne. real(nin,kind=GP) ) THEN
          WRITE(*,*)'dojpdfr: inconsistent data: nin=',nin, ' no. samples: ', fck
          STOP
        ENDIF
        DO j = 1, 2
          IF ( dolog(j) .GT. 0 ) fmin(j) = 10.0_GP**fmin(j)
          IF ( dolog(j) .GT. 0 ) fmax(j) = 10.0_GP**fmax(j)
        ENDDO

        WRITE(shead,'(A1,2(A4,A6,E16.8,A1,E16.8,A3,A4,A5,E16.8,A2),A6,I4,A1,I4,A9,I1,A1,I1,A1)') '#',&
        sR1(1:4),'_rng=[', fmin(1), ',' , fmax(1),']; ', sR1(1:4),'_avg=',gavg(1),'; ',&
        sR2(1:4),'_rng=[', fmin(2), ',' , fmax(2), ']; ',sR2(1:4),'_avg=',gavg(2),'; ',&
        'nbin=[', nbins(1),',',nbins(2), ']; blog=[', dolog(1),',',dolog(2),']'   
         OPEN(1,file=fname)
         WRITE(1,'(A)') trim(shead)
         WRITE(1,40) gpdf2_
   40    FORMAT( E23.15 )
         CLOSE(1)
      ENDIF
 
      RETURN

      END SUBROUTINE dojpdfr
!-----------------------------------------------------------------
!-----------------------------------------------------------------

END MODULE gutils
