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
      REAL   , DIMENSION  (:), ALLOCATABLE   :: fpdf_
      REAL   , DIMENSION  (:), ALLOCATABLE   :: gpdf_
      REAL   , DIMENSION  (:), ALLOCATABLE   :: fpdf2_
      REAL   , DIMENSION  (:), ALLOCATABLE   :: gpdf2_
      INTEGER, DIMENSION  (:), ALLOCATABLE   :: ikeep_
      INTEGER  :: nbins_=2500,nbins2_(2)=(/2500,2500/),nikeep_=0
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
!$    USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(IN)                        :: nin
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(*) :: Rin

      INTEGER(KIND=GP) :: ie0,ie1
      INTEGER          :: i, j, k, m, nb

      nb = 8  ! no. bits per byte

      ie1 = 0
! Note using GP like this is not good practice....
!$omp parallel do if (nin.ge.nth) private(ie0,ie1,m)
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
      SUBROUTINE carray_byte_swap(Cin, nin)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Performs endian conversion on array of complex numbers
!
! Parameters
!     Cin : aribtrary-rank array whose values will be endian-swapped, returned
!     nin : dimension of Cin
!-----------------------------------------------------------------
      USE mpivars
!$    USE threads
      USE fprecision

      IMPLICIT NONE

      INTEGER, INTENT(IN)                           :: nin
      REAL   (KIND=GP)                              :: c,r
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(*) :: Cin

      INTEGER(KIND=GP) :: ir0,ic0,ie1,ie2
      INTEGER          :: i, j, k, m, nb

      nb = 8  ! no. bits per byte

      ie1 = 0
! Note using GP like this is not good practice....
!$omp parallel do if (nin.ge.nth) private(ir0,ic0,ie1,ie2,m,r,c)
      DO k = 1, nin
          ir0 = TRANSFER(real(Cin(k)), 0_GP)
          DO m = 1, GP
             CALL MVBITS( ir0, (GP-m)*nb, nb, ie1, (m-1)*nb  )
          END DO
          ic0 = TRANSFER(aimag(Cin(k)), 0_GP)
          DO m = 1, GP
             CALL MVBITS( ic0, (GP-m)*nb, nb, ie2, (m-1)*nb  )
          END DO
          r      = TRANSFER(ie1, 0.0_GP)
          c      = TRANSFER(ie2, 0.0_GP)
          Cin(k) = cmplx(r,c,kind=GP)
       END DO

      RETURN

      END SUBROUTINE carray_byte_swap
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
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(IN), DIMENSION(*)    :: Rin
      REAL(KIND=GP), INTENT(INOUT)               :: fmin,fmax
      INTEGER, INTENT(IN)                        :: nin,n,nbins,ifixdr,dolog
      CHARACTER(len=*), INTENT(IN)               :: fname

      REAL(KIND=GP)                              :: del,gmin,gmax,tmin,tmax,test
      REAL(KIND=GP)                              :: fmin1,fmax1
      REAL(KIND=GP)                              :: gavg,sig,sumr,xnorm
      REAL                                       :: fbin,fkeep,gkeep
      INTEGER                                    :: hwrite,i,ibin,nkeep
      CHARACTER(len=1024)                        :: shead

      IF ( .NOT. ALLOCATED(ikeep_) .OR. nin.GT.nikeep_ ) THEN
        IF ( ALLOCATED(ikeep_ ) ) DEALLOCATE(ikeep_)
        ALLOCATE(ikeep_(nin))
        nikeep_ = nin
      ENDIF
      IF ( nbins.GT.nbins_        .OR. &
          .NOT. ALLOCATED(fpdf_)  .OR. &
          .NOT. ALLOCATED(gpdf_)       &
          ) THEN
        ! Re-allocate if necessary:
        IF ( ALLOCATED (fpdf_) ) DEALLOCATE(fpdf_)
        IF ( ALLOCATED (gpdf_) ) DEALLOCATE(gpdf_)
        ALLOCATE(fpdf_(nbins))
        ALLOCATE(gpdf_(nbins))
        nbins_ = nbins
      ENDIF

      DO i = 1, nbins_
       fpdf_(i) = 0.0_GP
       gpdf_(i) = 0.0_GP
      ENDDO

      ! Compute dynamic range of PDF
      fmin1 = MINVAL(Rin(1:nin),nin)
      fmax1 = MAXVAL(Rin(1:nin),nin)
      CALL MPI_ALLREDUCE(fmin1,gmin,1, GC_REAL,      &
                         MPI_MIN,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(fmax1,gmax,1, GC_REAL,      &
                         MPI_MAX,MPI_COMM_WORLD,ierr)
      IF ( ifixdr .le. 0 ) THEN
        fmin = gmin
        fmax = gmax
      ENDIF
      IF ( dolog .GT. 0 ) THEN
        fmin = log10(abs(fmin)+tiny(1.0_GP))
        fmax = log10(abs(fmax)+tiny(1.0_GP))
      ENDIF
  
      tmin = fmin
      tmax = fmax
      IF ( dolog .GT. 0 ) THEN
        tmin = 10.0_GP**(fmin)
        tmax = 10.0_GP**(fmax)
      ENDIF
!
! Find indices that meet dyn. range criterion:
      fkeep = 0.0
!$omp parallel do 
      DO i = 1, nin
        IF ( Rin(i).GE.tmin.AND. Rin(i).LE.tmax ) THEN
!$omp critical
          fkeep = fkeep + 1.0
          ikeep_(int(fkeep)) = i
!$omp end critical
        ENDIF
      ENDDO
!$omp end parallel do 

      nkeep = int(fkeep)

! Check global samples:
      CALL MPI_ALLREDUCE(fkeep, gkeep, 1, MPI_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)

      IF ( gkeep .LE. 0.0 ) THEN
        IF ( myrank.eq.0 ) THEN
          WRITE(*,*) 'dopdfr: no samples: fmin,fmax=',fmin,fmax,'gmin,gmax=',gmin,gmax
          WRITE(*,*) 'dopdfr: file not written: ',trim(fname)
        ENDIF
        RETURN
      ENDIF        
      xnorm = 1.0_GP / gkeep

      sumr = 0.0_GP
!$omp parallel do default(shared) private(i) reduction(+:sumr)
      DO i = 1, nkeep
        sumr  = sumr +  Rin(ikeep_(i))
      ENDDO
!$omp end parallel do 

      CALL MPI_ALLREDUCE(sumr, gavg, 1, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      gavg = gavg*xnorm
 
      ! Compute standard deviation:
      sumr = 0.0_GP
!$omp parallel do default(shared) private(i) reduction(+:sumr)
      DO i = 1, nkeep
        sumr  = sumr +  (Rin(ikeep_(i))-gavg)**2
      ENDDO
!$omp end parallel do 
      CALL MPI_ALLREDUCE(sumr, sig, 1, GC_REAL, &
                        MPI_SUM, MPI_COMM_WORLD,ierr)
      sig = sqrt(sig*xnorm)

      del = ABS(fmax-fmin)/dble(nbins)
      ! Compute local PDF:
      IF ( dolog .GT. 0 ) THEN
!$omp parallel do private (ibin,test)
        DO i = 1, nkeep
          test = log10(abs(Rin(ikeep_(i)))+tiny(1.0_GP))
          ibin = NINT( ( test - fmin )/del+1 )
          ibin = MIN(MAX(ibin,1),nbins)
!$omp atomic
          fpdf_(ibin) = fpdf_(ibin) + 1.0_GP
        ENDDO
!$omp end parallel do 
      ELSE
!$omp parallel do private (ibin,test)
        DO i = 1, nkeep
          test = Rin(ikeep_(i))
          ibin = NINT( ( test - fmin )/del+1 )
          ibin = MIN(MAX(ibin,1),nbins)
!$omp atomic
          fpdf_(ibin) = fpdf_(ibin) + 1.0_GP
        ENDDO
!$omp end parallel do 
      ENDIF

      ! First, do a sanity check:
      fbin = 0.0
!$omp parallel do default(shared) reduction(+:fbin)
      DO i = 1, nbins
        fbin = fbin + fpdf_(i)
      ENDDO
!$omp end parallel do 
      IF (fbin.NE.fkeep ) THEN
        WRITE (*,*)myrank,': dopdfr: inconsistent data: expected: ',fkeep, ' found: ',fbin
        WRITE (*,*)myrank,': dopdfr: fmin=',fmin,' fmax=',fmax,' nbins=',nbins,' dolog=',dolog,' ifixdr=',ifixdr,' del=',del
        WRITE (*,*)myrank,': dopdfr: file ', fname, ' not written.'
      ENDIF

      ! Compute global reduction between MPI tasks:
      CALL MPI_ALLREDUCE(fpdf_, gpdf_, nbins_, MPI_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      IF ( ierr.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)'dopdf: error in fpdf reduction'
        STOP
      ENDIF

      ! Write PDF to disk:
      CALL GTStart(hwrite,GT_WTIME)
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
        CALL GTStop(hwrite)
        write(*,*)'dopdf: file: ',trim(fname),': write time: ',GTGetTime(hwrite)
      ENDIF
      CALL GTFree(hwrite)
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
!     n      : full linear domain size
!     fname  : output interval id extension
!     nbins  : 2D array giving number of bins to use for PDF (>0) 
!              in each 'direction'. Must be the same for all MPI tasks on entry.
!     ifixdr : if > 0, will use fmin(*),fmax(*) as bounds for dynamic
!              range; else it will compute them dynamically
!     fmin   : highest dynamic range value, returned if ifixdr>0
!     fmax   : lowest dynamic range value, returned if ifixdr>0
!     dolog  : 2-element vectory indicating whether take log of quantity magnitude 
!              in 1-, 2- directions when creating bins (0, 1)
!-----------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE mpivars
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      REAL(KIND=GP)   , INTENT   (IN), DIMENSION(nin):: R1,R2
      REAL(KIND=GP)   , INTENT(INOUT)                :: fmin(2),fmax(2)
      INTEGER         , INTENT   (IN)                :: nin,n,nbins(2),ifixdr(2),dolog(2)
      CHARACTER(len=*), INTENT   (IN)                :: sR1,sR2,fname

      REAL(KIND=GP)                                  :: del (2),gmin(2),gmax(2),tmin(2),tmax(2),test(2)
      REAL(KIND=GP)                                  :: fmin1(2),fmax1(2)
      REAL(KIND=GP)                                  :: aa,gavg(2),sumr(2),sig(2),xnorm(2)
      REAL                                           :: fbin,fkeep,gkeep
      INTEGER                                        :: hwrite,i,j,jx,jy,nkeep
      CHARACTER(len=2048)                            :: shead

!if ( myrank.eq. 0 ) write(*,*)'dojpdf: sR1=',sR1,' sR2=',sR2, ' nbins=',nbins, ' ifixdr=',ifixdr, ' dolog=',dolog

      IF ( .NOT. ALLOCATED(ikeep_) .OR. nin.GT.nikeep_ ) THEN
        IF ( ALLOCATED(ikeep_ ) )DEALLOCATE(ikeep_)
        ALLOCATE(ikeep_(nin))
        nikeep_ = nin
      ENDIF
!if ( myrank.eq. 0 ) write(*,*)'dojpdf: checking need for allocation:'
      IF ( nbins(1).NE.nbins2_(1) .OR. &
           nbins(2).NE.nbins2_(2) .OR. &
          .NOT. ALLOCATED(fpdf2_) .OR. &
          .NOT. ALLOCATED(gpdf2_)      &
         ) THEN
!if ( myrank.eq. 0 ) write(*,*)'dojpdf: allocating fpdf2,gpdf2:'

        ! Re-allocate if necessary:
        IF ( ALLOCATED(fpdf2_) ) DEALLOCATE(fpdf2_)
        IF ( ALLOCATED(gpdf2_) ) DEALLOCATE(gpdf2_)
        ALLOCATE(fpdf2_(nbins(1)*nbins(2)))
        ALLOCATE(gpdf2_(nbins(1)*nbins(2)))
        nbins2_(1:2) = nbins(1:2)
      ENDIF
      fpdf2_(1:nbins2_(1)*nbins2_(2)) = 0.0_GP
      gpdf2_(1:nbins2_(1)*nbins2_(2)) = 0.0_GP

      ! Compute dynamic range of PDF
      fmin1(1) = MINVAL(R1(1:nin),nin)
      fmax1(1) = MAXVAL(R1(1:nin),nin)
      fmin1(2) = MINVAL(R2(1:nin),nin)
      fmax1(2) = MAXVAL(R2(1:nin),nin)
!if ( myrank.eq. 0 ) write(*,*)'dojpdf: reductions for min/max:'
      CALL MPI_ALLREDUCE(fmin1,gmin,2, GC_REAL,      &
                         MPI_MIN,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(fmax1,gmax,2, GC_REAL,      &
                         MPI_MAX,MPI_COMM_WORLD,ierr)
      IF ( ifixdr(1) .le. 0 ) THEN
        fmin(1) = gmin(1)
        fmax(1) = gmax(1)
      ENDIF
      IF ( ifixdr(2) .le. 0 ) THEN
        fmin(2) = gmin(2)
        fmax(2) = gmax(2)
      ENDIF
      DO j = 1, 2
        IF ( dolog(j) .GT. 0 ) fmin(j) = log10(fmin(j)+tiny(1.0_GP))
        IF ( dolog(j) .GT. 0 ) fmax(j) = log10(fmax(j)+tiny(1.0_GP))
      ENDDO
  
      tmin(1:2) = fmin(1:2)
      tmax(1:2) = fmax(1:2)
      IF ( dolog(1).GT.0 ) THEN
        tmin(1) = 10.0_GP**fmin(1)
        tmin(2) = 10.0_GP**fmin(2)
        tmax(1) = 10.0_GP**fmax(1)
        tmax(2) = 10.0_GP**fmax(2)
      ENDIF

!
! Find indices that meet dyn. range criterion:
      fkeep = 0.0
!$omp parallel do 
      DO i = 1, nin
        IF ( R1(i).GE.tmin(1).AND. R1(i).LE.tmax(1) .AND. &
             R2(i).GE.tmin(2).AND. R2(i).LE.tmax(2) ) THEN
!$omp critical
          fkeep = fkeep + 1.0
          ikeep_(int(fkeep)) = i
!$omp end critical
        ENDIF
      ENDDO
!$omp end parallel do 

!if ( myrank.eq. 0 ) write(*,*)'dojpdf: reductions for gkeep:'
      CALL MPI_ALLREDUCE(fkeep, gkeep, 1, MPI_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      IF ( gkeep.LE.0.0 ) THEN
        IF ( myrank.eq.0 ) THEN
          WRITE(*,*) 'dojpdfr: no samples, fmin=',fmin(1:2), ' fmax=',fmax(1:2), ' gmin==',gmin(1:2),' gmax=',gmax(1:2)
          WRITE(*,*) 'dojpdfr: file not written: ',trim(fname)
        ENDIF
        RETURN
      ENDIF        
      nkeep = int(fkeep)


      aa = 0.0_GP
!$omp parallel do default(shared) private(i) reduction(+:aa)
      DO i = 1, nkeep
        aa  = aa +  R1(ikeep_(i))
      ENDDO
!$omp end parallel do 
      sumr(1) = aa
  
      aa = 0.0_GP
!$omp parallel do default(shared) private(i) reduction(+:aa)
      DO i = 1, nkeep
        aa  = aa +  R2(ikeep_(i))
       ENDDO
!$omp end parallel do 
       sumr(2) = aa

      xnorm(1:2) = 1.0_GP/gkeep
!if ( myrank.eq. 0 ) write(*,*)'dojpdf: reductions for avg:'
      CALL MPI_ALLREDUCE(sumr, gavg, 2, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      gavg(1:2) = gavg(1:2)*xnorm(1:2)

      aa = 0.0_GP
!$omp parallel do default(shared) private(i) reduction(+:aa)
      DO i = 1, nkeep
        aa  = aa +  (R1(ikeep_(i))-gavg(1))**2
      ENDDO
!$omp end parallel do 
      sumr(1) = aa

      aa = 0.0_GP
!$omp parallel do default(shared) private(i) reduction(+:aa)
      DO i = 1, nkeep
        aa  = aa +  (R2(ikeep_(i))-gavg(2))**2
      ENDDO
!$omp end parallel do 
      sumr(2) = aa

!if ( myrank.eq. 0 ) write(*,*)'dojpdf: reductions for var:'
      CALL MPI_ALLREDUCE(sumr, sig, 2, GC_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      DO i = 1,2
        sig(i) = sqrt(sig(i)*xnorm(i))
      ENDDO


      ! Compute local PDF:
      DO j = 1, 2
        del (j) = ABS(fmax(j)-fmin(j))/dble(nbins(j))
      ENDDO

!if ( myrank.eq. 0 ) write(*,*)'dojpdf: compute local pdf:'
      ! Compute local PDF:
      IF ( dolog(1).GT.0 .AND. dolog(2).GT.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nkeep
          test(1) = log10(abs(R1(ikeep_(i)))+tiny(1.0_GP))
          test(2) = log10(abs(R2(ikeep_(i)))+tiny(1.0_GP))
          jx  = NINT( ( test(1) - fmin(1) )/del(1) )
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( ( test(2) - fmin(2) )/del(2) )
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx+(jy-1)*nbins(1)) = fpdf2_(jx+(jy-1)*nbins(1)) + 1.0_GP
        ENDDO
!$omp end parallel do 
      ELSE IF ( dolog(1).GT.0 .AND. dolog(2).LE.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nkeep
          test(1) = log10(abs(R1(ikeep_(i)))+tiny(1.0_GP))
          test(2) = R2(ikeep_(i))
          jx  = NINT( ( test(1) - fmin(1) )/del(1) )
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( ( test(2) - fmin(2) )/del(2) )
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx+(jy-1)*nbins(1)) = fpdf2_(jx+(jy-1)*nbins(1)) + 1.0_GP
        ENDDO
!$omp end parallel do 
      ELSE IF ( dolog(1).LE.0 .AND. dolog(2).GT.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nkeep
          test(1) = R1(ikeep_(i))
          test(2) = log10(abs(R2(ikeep_(i)))+tiny(1.0_GP))
          jx  = NINT( ( test(1) - fmin(1) )/del(1) )
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( ( test(2) - fmin(2) )/del(2) )
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx+(jy-1)*nbins(1)) = fpdf2_(jx+(jy-1)*nbins(1)) + 1.0_GP
        ENDDO
!$omp end parallel do 
      ELSE IF ( dolog(1).LE.0 .AND. dolog(2).LE.0 ) THEN
!$omp parallel do private (jx,jy,test) 
        DO i = 1, nkeep
          test(1) = R1(ikeep_(i))
          test(2) = R2(ikeep_(i))
          jx  = NINT( ( test(1) - fmin(1) )/del(1) )
          jx  = MIN(MAX(jx,1),nbins(1))
          jy  = NINT( ( test(2) - fmin(2) )/del(2) )
          jy  = MIN(MAX(jy,1),nbins(2))
!$omp atomic
          fpdf2_(jx+(jy-1)*nbins(1)) = fpdf2_(jx+(jy-1)*nbins(1)) + 1.0_GP
        ENDDO
!$omp end parallel do 
      ENDIF

! Check local data:
        fbin = 0.0
!$omp parallel do default(shared) private(i) reduction(+:fbin)
        DO j = 1, nbins(2)
          DO i = 1, nbins(1)
            fbin = fbin + fpdf2_(i+(j-1)*nbins(1))
          ENDDO
        ENDDO
!$omp end parallel do 

        IF ( fbin .ne. fkeep ) THEN
          WRITE(*,*) myrank,': dojpdfr: inconsistent data: expected: fkeep=',fkeep, ' found: ', &
                     fbin,' nbins=',nbins,' fmin=',fmin, ' fmax=',fmax,' dolog=',dolog,' ifixdr=',ifixdr
          STOP
        ENDIF
     
!if ( myrank.eq. 0 ) write(*,*)'dojpdf: compute global pdf: nbins2_=',nbins2_
      ! Compute global reduction between MPI tasks:
!     CALL MPI_REDUCE(fpdf2_, gpdf2_, nbins2_(1)*nbins2_(2), MPI_REAL, &
!                     MPI_SUM, 0, MPI_COMM_WORLD,ierr)
      CALL MPI_ALLREDUCE(fpdf2_, gpdf2_, nbins2_(1)*nbins2_(2), MPI_REAL, &
                      MPI_SUM, MPI_COMM_WORLD,ierr)
      IF ( ierr.NE.MPI_SUCCESS ) THEN
        WRITE(*,*)'dojpdf: error in fpdf reduction'
        STOP
      ENDIF

!if ( myrank.eq. 0 ) write(*,*)'dojpdf: writing to disk...'
      ! Write PDF to disk:
      CALL GTStart(hwrite,GT_WTIME)
      IF ( myrank.eq.0 ) THEN
        DO j = 1, 2
          IF ( dolog(j) .GT. 0 ) fmin(j) = 10.0_GP**fmin(j)
          IF ( dolog(j) .GT. 0 ) fmax(j) = 10.0_GP**fmax(j)
        ENDDO

        WRITE(shead, &
        '(A1,2(A4,A6,E16.8,A1,E16.8,A3,A4,A5,E16.8,A2,A4,A5,E16.8,A2),A6,I7,A1,I7,A9,I1,A1,I1,A10,F12.0,A1,F12.0,A1)')&
        '#',&
        trim(sR1),'_rng=[',fmin(1),',',fmax(1),']; ',trim(sR1),'_avg=',gavg(1),'; ',trim(sR1),'_sig=',sig(1),'; ',&
        trim(sR2),'_rng=[',fmin(2),',',fmax(2),']; ',trim(sR2),'_avg=',gavg(2),'; ',trim(sR2),'_sig=',sig(2),'; ',&
        'nbin=[', nbins(1),',',nbins(2), ']; blog=[', dolog(1),',',dolog(2),']; nkeep=',gkeep
         OPEN(1,file=fname,iostat=ierr)
         IF ( ierr.ne.0 ) THEN
           WRITE(*,*) 'dojpdf2: error opening file: ', fname
           STOP
         ENDIF
         WRITE(1,'(A)') trim(shead)
         WRITE(1,40) gpdf2_
   40    FORMAT( E23.15 )
         CLOSE(1)
         CALL GTStop(hwrite)
         write(*,*)'dojpdf: file: ',trim(fname),': write time: ',GTGetTime(hwrite)
      ENDIF
!if ( myrank.eq. 0 ) write(*,*)'dojpdf: writing to disk done.'
      CALL GTFree(hwrite)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 

      RETURN

      END SUBROUTINE dojpdfr
!-----------------------------------------------------------------
!-----------------------------------------------------------------

      SUBROUTINE rarray_props(Rin,nin,n,rmin,rmax,mean,rrms,rstd)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Find max of input array
!
! Parameters
!     Rin : arbitrary-rank array 
!     nin : local dimension of Rin
!     n   : linear dimension
!     rmin: min value found
!     rmax: max value found
!     mean: mean
!     rrms: rms
!     std : std deviation from the mean
!-----------------------------------------------------------------
      USE mpivars
!$    USE threads
      USE fprecision
      USE commtypes

      IMPLICIT NONE

      INTEGER, INTENT(IN)                        :: n,nin 
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(*) :: Rin 
      REAL(KIND=GP), INTENT  (OUT)               :: rmin,rmax,rrms,mean,rstd
      REAL(KIND=GP)                              :: loc(2),glo(2),ng

      INTEGER          :: k

      ng = 1.0_GP/real(n,kind=GP)**3
      rmax   = -huge(rmax)
      rmin   =  huge(rmax)
      loc(1) = 0.0_GP
      loc(2) = 0.0_GP
!$omp parallel do if (nin.ge.nth) 
      DO k = 1, nin 
!$omp critical
          rmax   = max(rmax,Rin(k))
          rmin   = min(rmin,Rin(k))
          loc(1) = loc(1) + Rin(k)**2 ! rms
          loc(2) = loc(2) + Rin(k)    ! sum/mean
!$omp end critical
       END DO
       CALL MPI_ALLREDUCE(loc,glo,2,GC_REAL,MPI_SUM, &
                          MPI_COMM_WORLD,ierr)
       rrms = sqrt(glo(1)*ng)
       mean = glo(2)*ng

       loc(1) = rmin
       CALL MPI_ALLREDUCE(loc,glo,1,GC_REAL,MPI_MIN, &
                          MPI_COMM_WORLD,ierr)
       rmin = glo(1)

       loc(1) = rmax
       CALL MPI_ALLREDUCE(loc,glo,1,GC_REAL,MPI_MAX, &
                          MPI_COMM_WORLD,ierr)
       rmax = glo(1)

      loc = 0.0_GP
!$omp parallel do if (nin.ge.nth) 
      DO k = 1, nin 
!$omp atomic
          loc(1)= loc(1) + (Rin(k)-mean)**2
       END DO
  
       CALL MPI_ALLREDUCE(loc,glo,1,GC_REAL,MPI_SUM, &
                          MPI_COMM_WORLD,ierr)
       rstd = sqrt(glo(1)*ng)

      RETURN

      END SUBROUTINE rarray_props
!-----------------------------------------------------------------
!-----------------------------------------------------------------

END MODULE gutils
