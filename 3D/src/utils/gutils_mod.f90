!=================================================================
! GHOST suite Fortran utilities 
!
! 2010 D. Rosenberg
!      NCAR
!
! 17 Nov 2010: Initial version
!=================================================================
MODULE gutils
!
! Utils data, if any:
!
! ...
! end, member data
      REAL, DIMENSION  (:), ALLOCATABLE   :: fpdf_
      REAL, DIMENSION  (:), ALLOCATABLE   :: gpdf_
      REAL, DIMENSION(:,:), ALLOCATABLE   :: fpdf2_
      REAL, DIMENSION(:,:), ALLOCATABLE   :: gpdf2_
      INTEGER                             :: nbins_=1000,nbins2_(2)=(/1000,1000/)
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

      REAL(KIND=GP)                              :: del,fbin,gmin,gmax,tmin,tmax,test
      REAL(KIND=GP)                              :: gavg,sig,sumr,xnorm
      REAL(KIND=GP)                              :: fkeep, gkeep
      INTEGER                                    :: i,ibin,nkeep
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

      ! Compute dynamic range of PDF
      fmin = MINVAL(Rin(1:nin),nin)
      fmax = MAXVAL(Rin(1:nin),nin)
      CALL MPI_ALLREDUCE(fmin,gmin,1, GC_REAL,      &
                         MPI_MIN,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(fmax,gmax,1, GC_REAL,      &
                         MPI_MAX,MPI_COMM_WORLD,ierr)
      IF ( ifixdr .le. 0 ) THEN
        fmin = gmin
        fmax = gmax
      ENDIF
      IF ( dolog .GT. 0 ) THEN
        fmin = log10(fmin+tiny(1.0_GP))
        fmax = log10(fmax+tiny(1.0_GP))
      ENDIF
  
      tmin = fmin
      tmax = fmax
      IF ( dolog .GT. 0 ) THEN
        tmin = 10.0_GP**(fmin)
        tmax = 10.0_GP**(fmax)
      ENDIF

      sumr = 0.0_GP
!$omp parallel do private(i)
      DO i = 1, nin
        IF ( Rin(i).GE.tmin .AND. Rin(i).LE.tmax ) THEN
!$omp atomic
          sumr  = sumr +  Rin(i)
        ENDIF
      ENDDO

      fkeep = 0.0
!$omp parallel do private(i)
      DO i = 1, nin
        IF ( Rin(i).GE.tmin .AND. Rin(i).LE.tmax ) THEN
!$omp atomic
          fkeep = fkeep + 1.0
        ENDIF
      ENDDO
      CALL MPI_ALLREDUCE(fkeep, gkeep, 1, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
!write(*,*)'dopdf: tmin=',tmin,' tmax=',tmax,' nkeep=',fkeep,' gkeep=',gkeep

      IF ( gkeep .LE. 0.0 .AND.myrank.EQ.0 ) THEN
        WRITE(*,*) myrank,': dopdfr: no samples: fmin,fmax=',fmin,fmax,'gmin,gmax=',gmin,gmax
        RETURN
      ENDIF        
      xnorm = 1.0_GP / gkeep
      CALL MPI_ALLREDUCE(sumr, gavg, 1, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      gavg = gavg*xnorm
 
      ! Compute standard deviation:
      sumr = 0.0_GP
!$omp parallel do 
      DO i = 1, nin
        IF ( Rin(i).GE.tmin.AND. Rin(i).LE.tmax ) THEN
!$omp atomic
          sumr  = sumr +  (Rin(i)-gavg)**2
        ENDIF
      ENDDO
      CALL MPI_ALLREDUCE(sumr, sig, 1, GC_REAL, &
                        MPI_SUM, MPI_COMM_WORLD,ierr)
      sig = sqrt(sig*xnorm)

      del = ABS(fmax-fmin)/dble(nbins)
      ! Compute local PDF:
      IF ( dolog .GT. 0 ) THEN
!$omp parallel do private (ibin,test)
        DO i = 1, nin
          test = log10(abs(Rin(i))+tiny(1.0_GP))
          IF ( test .GE.tmin .AND. test.LE.tmax ) THEN
            ibin = NINT( ( test - fmin )/del+1 )
            ibin = MIN(MAX(ibin,1),nbins)
!$omp atomic
            fpdf_(ibin) = fpdf_(ibin) + 1.0_GP
          ENDIF
        ENDDO
      ELSE
!$omp parallel do private (ibin,test)
        DO i = 1, nin
          test = Rin(i)
          IF ( test.GE.tmin .AND. test.LE.tmax ) THEN
            ibin = NINT( ( test - fmin )/del+1 )
            ibin = MIN(MAX(ibin,1),nbins)
!$omp atomic
            fpdf_(ibin) = fpdf_(ibin) + 1.0_GP
          ENDIF
        ENDDO
      ENDIF

      ! Compute global reduction between MPI tasks:
      CALL MPI_REDUCE(fpdf_, gpdf_, nbins, MPI_REAL, &
                      MPI_SUM, 0, MPI_COMM_WORLD,ierr)
!     CALL MPI_ALLREDUCE(fpdf_, gpdf_, nbins, MPI_REAL, &
!                     MPI_SUM, MPI_COMM_WORLD,ierr)

      ! First, do a sanity check:
      fbin = 0.0_GP
!$omp parallel do reduction(+:fbin)
      DO i = 1, nbins
        fbin = fbin + fpdf_(i)
      ENDDO
      IF (fbin.NE.fkeep ) THEN
        WRITE (*,*)myrank, ': dopdfr: inconsistent data: expected: ',fkeep, ' found: ',fbin
        WRITE (*,*)myrank,': dopdfr: fmin=',fmin,' fmax=',fmax,' nbins=',nbins,' dolog=',dolog,' ifixdr=',ifixdr,' del=',del
        WRITE (*,*)myrank,': dopdfr: file ', fname, ' not written.'
      ENDIF
     
      ! Write PDF to disk:
      IF ( myrank.eq.0 ) THEN
        IF ( dolog .GT. 0 ) THEN
          fmin = 10.0_GP**fmin
          fmax = 10.0_GP**fmax
        ENDIF
        WRITE(shead,'(A9,E16.8,A1,E16.8,A7,E16.8,A6,E16.8,A7,I7,A7,I1,A8,F12.0)') '# range=[', & 
        fmin, ',' , fmax, ']; avg=',gavg,'; sig=',sig,'; nbin=', nbins, '; blog=', dolog,   &
        '; nkeep=',gkeep   
        OPEN(1,file=trim(fname))
        WRITE(1,'(A)') trim(shead)
        WRITE(1,40) gpdf_(1:nbins)
   40   FORMAT( E23.15 )
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

      REAL(KIND=GP)                                  :: del (2),fck,gmin(2),gmax(2),tmin(2),tmax(2),test(2)
      REAL(KIND=GP)                                  :: aa,gavg(2),sumr(2),sig(2),xnorm(2)
      REAL(KIND=GP)                                  :: fkeep(2),gkeep(2),fkeep2,gkeep2
      INTEGER                                        :: i,j,jx,jy,nkeep(2)
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

      ! Compute dynamic range of PDF
      fmin(1) = MINVAL(R1(1:nin),nin)
      fmax(1) = MAXVAL(R1(1:nin),nin)
      fmin(2) = MINVAL(R2(1:nin),nin)
      fmax(2) = MAXVAL(R2(1:nin),nin)
      DO j = 1, 2
        IF ( dolog(j) .GT. 0 ) fmin(j) = log10(fmin(j)+tiny(1.0_GP))
        IF ( dolog(j) .GT. 0 ) fmax(j) = log10(fmax(j)+tiny(1.0_GP))
      ENDDO
      CALL MPI_ALLREDUCE(fmin,gmin,2, GC_REAL,      &
                         MPI_MIN,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(fmax,gmax,2, GC_REAL,      &
                         MPI_MAX,MPI_COMM_WORLD,ierr)
      IF ( ifixdr .le. 0 ) THEN
        fmin(1:2) = gmin(1:2)
        fmax(1:2) = gmax(1:2)
      ENDIF
  
      tmin(1:2) = fmin(1:2)
      tmax(1:2) = fmax(1:2)
      IF ( dolog(1).GT.0 ) THEN
        tmin(1) = 10.0_GP**fmin(1)
        tmin(2) = 10.0_GP**fmin(2)
        tmax(1) = 10.0_GP**fmax(1)
        tmax(2) = 10.0_GP**fmax(2)
      ENDIF

      fkeep(1)  = 0.0
      aa = 0.0_GP
!$omp parallel do private(i)
      DO i = 1, nin
        IF ( R1(i).GE.tmin(1) .AND. R1(i).LE.tmax(1) ) THEN
!$omp critical
          aa  = aa +  R1(i)
          fkeep(1) = fkeep(1) + 1.0
!$omp end critical
        ENDIF
      ENDDO
      sumr(1) = aa
  
      fkeep(2)  = 0.0
      aa = 0.0_GP
!$omp parallel do 
      DO i = 1, nin
        IF ( R2(i).GE.tmin(2) .AND. R2(i).LE.tmax(2) ) THEN
!$omp critical
          aa  = aa +  R2(i)
          fkeep(2) = fkeep(2) + 1.0
!$omp end critical
        ENDIF
       ENDDO
       sumr(2) = aa

      CALL MPI_ALLREDUCE(fkeep, gkeep, 2, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      IF ( gkeep(1) .LE. 0.0 .OR. gkeep(1).LE.0.0  ) THEN
        WRITE(*,*) 'dopdfr: no samples, fmin=',fmin(1:2), ' fmax=',fmax(1:2), ' gmin==',gmin(1:2),' gmax=',gmax(1:2)
        STOP
      ENDIF        

      xnorm(1:2) = 1.0_GP/gkeep(1:2)
      CALL MPI_ALLREDUCE(sumr, gavg, 2, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      gavg(1:2) = gavg(1:2)*xnorm(1:2)

      aa = 0.0_GP
!$omp parallel do 
      DO i = 1, nin
        IF ( R1(i).GE.tmin(1) .AND. R1(i).LE.tmax(1) ) THEN
!$omp atomic
          aa  = aa +  (R1(i)-gavg(1))**2
        ENDIF
      ENDDO
      sumr(1) = aa

      aa = 0.0_GP
!$omp parallel do 
      DO i = 1, nin
        IF ( R2(i).GE.tmin(2) .AND. R2(i).LE.tmax(2) ) THEN
!$omp atomic
         aa  = aa +  (R2(i)-gavg(2))**2
        ENDIF
      ENDDO
      sumr(2) = aa
      CALL MPI_ALLREDUCE(sumr, sig, 2, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      sig(1:2) = sqrt(sig(1:2)*xnorm(1:2))


      ! Compute local PDF:
      DO j = 1, 2
        del (j) = ABS(fmax(j)-fmin(j))/dble(nbins(j))
      ENDDO

      fkeep2 = 0.0
      ! Compute local PDF:
      IF ( dolog(1).GT.0 .AND. dolog(2).GT.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nin
          test(1) = log10(abs(R1(i))+tiny(1.0_GP))
          test(2) = log10(abs(R2(i))+tiny(1.0_GP))
          IF ( test(1).GE.tmin(1) .AND. test(1).LE.tmax(1) .AND. &
               test(2).GE.tmin(2) .AND. test(2).LE.tmax(2) ) THEN
            jx  = NINT( ( test(1) - fmin(1) )/del(1) )
            jx  = MIN(MAX(jx,1),nbins(1))
            jy  = NINT( ( test(2) - fmin(2) )/del(2) )
            jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
            fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
!$omp atomic
            fkeep2 = fkeep2 + 1.0
          ENDIF
        ENDDO
      ELSE IF ( dolog(1).GT.0 .AND. dolog(2).LE.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nin
          test(1) = log10(abs(R1(i))+tiny(1.0_GP))
          test(2) = R2(i)
          IF ( test(1).GE.tmin(1) .AND. test(1).LE.tmax(1) .AND. &
               test(2).GE.tmin(2) .AND. test(2).LE.tmax(2) ) THEN
            jx  = NINT( ( test(1) - fmin(1) )/del(1) )
            jx  = MIN(MAX(jx,1),nbins(1))
            jy  = NINT( ( test(2) - fmin(2) )/del(2) )
            jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
            fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
!$omp atomic
            fkeep2 = fkeep2 + 1.0
          ENDIF
        ENDDO
      ELSE IF ( dolog(1).LE.0 .AND. dolog(2).GT.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nin
          test(1) = R1(i)
          test(2) = log10(abs(R2(i))+tiny(1.0_GP))
          IF ( test(1).GE.tmin(1) .AND. test(1).LE.tmax(1) .AND. &
               test(2).GE.tmin(2) .AND. test(2).LE.tmax(2) ) THEN
            jx  = NINT( ( test(1) - fmin(1) )/del(1) )
            jx  = MIN(MAX(jx,1),nbins(1))
            jy  = NINT( ( test(2) - fmin(2) )/del(2) )
            jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
            fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
!$omp atomic
            fkeep2 = fkeep2 + 1.0
          ENDIF
        ENDDO
      ELSE IF ( dolog(1).LE.0 .AND. dolog(2).LE.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nin
          test(1) = R1(i)
          test(2) = R2(i)
          IF ( test(1).GE.tmin(1) .AND. test(1).LE.tmax(1) .AND. &
               test(2).GE.tmin(2) .AND. test(2).LE.tmax(2) ) THEN
            jx  = NINT( ( test(1) - fmin(1) )/del(1) )
            jx  = MIN(MAX(jx,1),nbins(1))
            jy  = NINT( ( test(2) - fmin(2) )/del(2) )
            jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
            fpdf2_(jx,jy) = fpdf2_(jx,jy) + 1.0_GP
!$omp atomic
            fkeep2 = fkeep2 + 1.0
          ENDIF
        ENDDO
      ENDIF
     
      CALL MPI_ALLREDUCE(fkeep2, gkeep2, 1, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)

      ! Compute global reduction between MPI tasks:
      CALL MPI_REDUCE(fpdf2_, gpdf2_, nbins(1)*nbins(2), MPI_REAL, &
                      MPI_SUM, 0, MPI_COMM_WORLD,ierr)

      ! Write PDF to disk:
      IF ( myrank.eq.0 ) THEN
! Check:
        fck = 0.0_GP
!$omp parallel do private(j) reduction(+:fck)
        DO i = 1, nbins(1)
          DO j = 1, nbins(2)
            fck = fck + fpdf2_(i,j)
          ENDDO
        ENDDO
        IF ( fck .ne. fkeep2 ) THEN
          WRITE(*,*)'dojpdfr: inconsistent data: expected: fkeep=',fkeep2, ' found: ', fck
          STOP
        ENDIF
        DO j = 1, 2
          IF ( dolog(j) .GT. 0 ) fmin(j) = 10.0_GP**fmin(j)
          IF ( dolog(j) .GT. 0 ) fmax(j) = 10.0_GP**fmax(j)
        ENDDO

        WRITE(shead,'(A1,2(A4,A6,E16.8,A1,E16.8,A3,A4,A5,E16.8,A2,A4,A5,E16.8,A2),A6,I7,A1,I7,A9,I1,A1,I1,A10,F12.0,A1,F12.0,A1)') '#',&
        sR1(1:4),'_rng=[',fmin(1),',',fmax(1),']; ',sR1(1:4),'_avg=',gavg(1),'; ',sR1(1:4),'_sig=',sig(1),'; ',&
        sR2(1:4),'_rng=[',fmin(2),',',fmax(2),']; ',sR2(1:4),'_avg=',gavg(2),'; ',sR2(1:4),'_sig=',sig(2),'; ',&
        'nbin=[', nbins(1),',',nbins(2), ']; blog=[', dolog(1),',',dolog(2),']; nkeep=',gkeep2
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
