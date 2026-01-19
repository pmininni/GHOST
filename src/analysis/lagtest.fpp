!=================================================================
      program lagtest
!=================================================================
! Lagrangian particle test stub. This tool ONLY works with cubic
! data in (2.pi)^3 domains.
!
! 2013 D. Rosenberg
!      NCCS:ORNL
!=================================================================

!
! Definitions for conditional compilation

! Modules

      use fprecision
      use commtypes
      use mpivars
      use filefmt
      use iovar
      use iompi
      use grid
      use threads
      use gutils
      use random
      use class_GPartComm
      use class_GPSplineInt
      use class_GPart

      implicit none

      real(KIND=GP),ALLOCATABLE,DIMENSION (:,:,:) :: bot,tmp1,tmp2,top,v,vb,ve,vt
      real(KIND=GP),ALLOCATABLE,DIMENSION   (:,:) :: vdb,ptmp
      real(KIND=GP),ALLOCATABLE,DIMENSION     (:) :: dv,px,py,pz,vlag,valag
      real(KIND=GP),ALLOCATABLE,DIMENSION     (:) :: gpx,gpy,gpz
      real(KIND=GP)                               :: del,r,suml,sumel,sumg,sumeg,time,x,y,z
      real(KIND=GP)                               :: kx,ky,kz,pi,tiny,xbnds(3,2),xi,xj,xk,xm,xn
      integer                                     :: n
      integer      ,ALLOCATABLE,DIMENSION     (:) :: id,gid
      integer                                     :: irank,i,ir,j,jm,k,ki,km,kbsta,kbend,ktsta,ktend,no
      integer                                     :: maxparts,nd(3),nparts,intorder
      integer                                     :: ibnds(3,2),idims(3),iseed
      integer                                     :: obnds(3,2)
      integer                                     :: csize,nstrip,tind
      integer                                     :: fh,kmax
      type(GPSplineInt)                           :: gpspline
      type  (GPartComm)                           :: gpcomm
      type      (GPart)                           :: lagpart

      character(len=8)   :: suff
      character(len=100) :: odir,idir
      character(len=256) :: fname,fout,msg
      character(len=1024):: fnlist

!
! Verifies proper compilation of the tool

      IF ( (nx.ne.ny).or.(ny.ne.nz) ) THEN
        IF (myrank.eq.0) &
           PRINT *,'This tool only works with cubic data in (2.pi)^3 domains'
        STOP
      ENDIF
      n = nx
!    
!
! Initializes the MPI and I/O libraries
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!     call range(1,n/2+1,nprocs,myrank,ista,iend)
      call range(1,n,nprocs,myrank,ista,iend)
      call range(1,n,nprocs,myrank,ksta,kend)
 
      write(*,*) myrank,' main: ksta=',ksta, ' kend=',kend

      tiny = 1.0e-8
      pi   = 4.0_GP * atan(1.0_GP)

      idir   = '.'
      odir   = '.'
      iseed  = 1000

      tind   = 1
      write(ext, fmtext) tind

      intorder = 3
      csize    = 8
      nstrip   = 1

      maxparts = 32
!
      nd  (1)  = n
      nd  (2)  = n
      nd  (3)  = n
!
      ibnds(1,1) = 1
      ibnds(1,2) = n
      ibnds(2,1) = 1
      ibnds(2,2) = n
      ibnds(3,1) = ksta
      ibnds(3,2) = kend

      obnds(1,1) = 1
      obnds(1,2) = n
      obnds(2,1) = 1
      obnds(2,2) = n
      obnds(3,1) = ista
      obnds(3,2) = iend
!
      xbnds(1,1) = 0.0
      xbnds(1,2) = real(n-1,kind=GP)
      xbnds(2,1) = 0.0
      xbnds(2,2) = real(n-1,kind=GP)
      xbnds(3,1) = real(ksta-1,kind=GP)-0.5_GP+epsilon(1.0_GP)
      xbnds(3,2) = real(kend-1,kind=GP)+0.5_GP-epsilon(1.0_GP)
!
      allocate( ve   (n,n,1:(kend-ksta+1+2*(intorder-1))) )
      allocate( v    (n,n,ksta:kend) )
      allocate( vb   (n,n,ksta:kend) )
      allocate( tmp1 (ista:iend,n,n) )
      allocate( tmp2 (ista:iend,n,n) )
      allocate( vt   (n,n,ista:iend) )
      allocate( top (n,n,intorder-1) )
      allocate( bot (n,n,intorder-1) )
      allocate( vdb     (3,maxparts) )
      allocate( ptmp    (3,maxparts) )
      allocate( id        (maxparts) )
      allocate( px        (maxparts) )
      allocate( py        (maxparts) )
      allocate( pz        (maxparts) )
      allocate( gid       (maxparts) )
      allocate( gpx       (maxparts) )
      allocate( gpy       (maxparts) )
      allocate( gpz       (maxparts) )
      allocate( vlag      (maxparts) )
      allocate( valag     (maxparts) )
      allocate( dv        (maxparts) )

     call gpcomm%GPartComm_ctor(0,maxparts, &
          nd,intorder-1,MPI_COMM_WORLD)
!
! set v with linear index made up from _global_ indices:
      xn = real(n,kind=GP)
      do k = ksta,kend
         xk = real(k,kind=GP)
         do j = 1,n
           xj = real(j,kind=GP)
           do i = 1,n
             xi = real(i,kind=GP)
             v  (i,j,k) = xi + (xj-1.0)*xn + (xk-1.0)*xn*xn 
            enddo
         enddo
      enddo
!
! fill top check buffer for SlabExchenage check:
      xn = real(n,kind=GP)
      no = intorder-1
      if ( nprocs .gt. 1 ) then
        jm = 0; ir = 1
        do while ( jm.lt.no)
          irank = mod(myrank+ir,nprocs)
          call range(1,n,nprocs,irank,ktsta,ktend)
          do k = ktsta,min(ktsta+no-1,ktend)
            if ( jm .ge. no ) cycle
            jm = jm + 1
            xk = real(k,kind=GP)
            do j = 1,n
              xj = real(j,kind=GP)
              do i = 1,n
                xi = real(i,kind=GP)
                top(i,j,jm) = xi + (xj-1.0)*xn + (xk-1.0)*xn*xn 
               enddo
            enddo
          enddo
          ir = ir + 1
        enddo
      else
        do k = ksta,ksta+no-1
           xk = real(k,kind=GP)
           do j = 1,n
             xj = real(j,kind=GP)
             do i = 1,n
               xi = real(i,kind=GP)
               top(i,j,k-ksta+1) = xi + (xj-1.0)*xn + (xk-1.0)*xn*xn 
              enddo
           enddo
        enddo
      endif
!
! fill bottom check buffer for SlabExchenage check:
      if ( nprocs .gt. 1 ) then
        jm = 0; ir = 1;
        do while ( jm.lt.no)
          irank = mod(myrank-ir+nprocs,nprocs)
          call range(1,n,nprocs,irank,kbsta,kbend)
          do k = kbend,max(kbend-no+1,kbsta),-1
            if ( jm .ge. no ) cycle
            jm = jm + 1
            xk = real(k,kind=GP)
            do j = 1,n
              xj = real(j,kind=GP)
              do i = 1,n
                xi = real(i,kind=GP)
                bot(i,j,no-jm+1) = xi + (xj-1.0)*xn + (xk-1.0)*xn*xn 
               enddo
            enddo
          enddo
          ir = ir + 1
        enddo
      else
        do k = kend-no+1,kend
          xk = real(k,kind=GP)
          do j = 1,n
            xj = real(j,kind=GP)
              do i = 1,n
                xi = real(i,kind=GP)
                bot(i,j,k-kend+no) = xi + (xj-1.0)*xn + (xk-1.0)*xn*xn 
              enddo
          enddo
        enddo
      endif
!
! Set particle positions:
      del  = 1.0_GP/(maxparts-1)
      do j = 1,maxparts
        gid(j) = j-1
        r      = 0.5*(randu(iseed)+1.0)
        gpx(j) = r*real(n,kind=GP)
        gpy(j) = r*real(n,kind=GP)
        gpz(j) = r*real(n,kind=GP)
!       gpx(j) = n*real(j,kind=GP)*del
!       gpy(j) = n*real(j,kind=GP)*del
!       gpz(j) = n*real(j,kind=GP)*del
        if ( myrank.eq.0 ) then
          write(*,*)'main: gpdb(',j,')=',gpx(j),gpy(j),gpz(j)
        endif
      enddo
      
      nparts = 0
      do j = 1,maxparts
        if ( gpz(j).GE.xbnds(3,1).AND.gpz(j).LE.xbnds(3,2) ) then
          nparts = nparts + 1
          id(nparts) = gid(j)
          px(nparts) = gpx(j)
          py(nparts) = gpy(j)
          pz(nparts) = gpz(j)
        endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Particle comm. interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( myrank .eq. 0 ) then 
        write(*,*) 'Checking VDBSynch...'
      endif
      call gpcomm%VDBSynch(vdb,maxparts,id,px,py,pz,nparts,ptmp)
#if 0
      write(*,*) 'myrank=',myrank,' nparts=',nparts,' local particles:'
      do j = 1,nparts
        write(*,100) myrank, '> ', id(j), ': pos=', px(j),', ',py(j),', ',pz(j)
 100    format(i2,a1,i4,a,f7.1,a,f7.1,a,f7.1)
      enddo
      
      if ( myrank .eq. 0 ) then
        write(*,*) 'globalparticles:'
        do j = 1,maxparts
          write(*,*) j-1, ': pos=', vdb(1,j),', ',vdb(2,j),', ',vdb(3,j)
        enddo
      endif
#endif

      suml  = 0.0
      sumel = 0.0
      do j = 1,nparts
        if ( id(j).lt. 0 .or. id(j).gt.maxparts-1 ) then
          write(*,*) myrank,': j=',j,' id=',id(j)
          sumel = 1.0
        endif
        suml = suml + vdb(1,id(j)+1) - px(j)
        suml = suml + vdb(2,id(j)+1) - py(j)
        suml = suml + vdb(3,id(j)+1) - pz(j)
      enddo

      call MPI_ALLREDUCE(suml ,sumg ,1,GC_REAL,   &
                       MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(sumel,sumeg,1,GC_REAL,   &
                       MPI_SUM,MPI_COMM_WORLD,ierr)

      if ( myrank .eq. 0 ) then 
        if ( sumeg .gt. tiny ) then 
          write(*,*) 'VDBSynch memory access error: id ; error=', sumeg
          stop
        endif
        if ( sumg .gt. tiny ) then
          write(*,*) '............................................VDBSynch FAILURE'
          write(*,*) 'VDBSynch error=', sumg
          stop
        else
          write(*,*) '............................................VDBSynch passed'
        endif
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Field comm. interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( myrank .eq. 0 ) then 
        write(*,*) 'Checking Transpose...'
      endif

      call gpcomm%GTranspose (vt,obnds,v ,ibnds,3,tmp1) ! make z-y complete
      call gpcomm%GITranspose(vb,ibnds,vt,obnds,3,tmp1) ! make x-y complete again
      suml = 0.0
      do k = ksta,kend
         do j = 1,n
           do i = 1,n
               suml  = suml + abs(v(i,j,k)-vb(i,j,k))
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(suml,sumg,1,GC_REAL,   &
                         MPI_SUM,MPI_COMM_WORLD,ierr)

      if ( myrank .eq. 0 ) then 
        if ( sumg .gt. tiny ) then
          write(*,*) '............................................Transpose FAILURE'
          write(*,*) 'Back Transpose error=', sumg
          stop
        else
          write(*,*) '............................................Transpose passed'
        endif
      endif
!
!
      call gpcomm%SlabDataExchangeSF(ve,v)

      if ( myrank .eq. 0 ) then 
        write(*,*) 'Checking Bottom SlabExch...'
      endif
      suml = 0.0
      do k = 1,no
         do j = 1,n
           do i = 1,n
               suml  = suml + abs(ve(i,j,k)-bot(i,j,k))
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(suml,sumg,1,GC_REAL,   &
                         MPI_SUM,MPI_COMM_WORLD,ierr)
      if ( myrank .eq. 0 ) then 
        if ( sumg .gt. tiny ) then
          write(*,*) '............................................SlabExch BOTTOM FAILURE'
          write(*,*) 'SlabExch error=', sumg
          stop
        else
          write(*,*) '............................................SlabExch BOTTOM passed'
        endif
      endif

      if ( myrank .eq. 0 ) then 
        write(*,*) 'Checking TOP SlabExch...'
      endif
      no = intorder-1
      suml = 0.0
      do k = 1,no
         do j = 1,n
           do i = 1,n
               suml  = suml + abs(ve(i,j,no+kend-ksta+1+k)-top(i,j,k))
            enddo
         enddo
      enddo

      call MPI_ALLREDUCE(suml,sumg,1,GC_REAL,   &
                         MPI_SUM,MPI_COMM_WORLD,ierr)

      if ( myrank .eq. 0 ) then 
        if ( sumg .gt. tiny ) then
          write(*,*) '............................................SlabExch TOP FAILURE'
          write(*,*) 'SlabExch error=', sumg
          stop
        else
          write(*,*) '............................................SlabExch TOP passed'
        endif
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Interpolation interfaces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call gpspline%GPSplineInt_ctor(3,nd,ibnds,maxparts,gpcomm)
      ! Set scalar according te 'xbnds'  grid:
      kx = 2.0_GP*pi/(nd(1)) !-1.0_GP)
      ky = 2.0_GP*pi/(nd(2)) !-1.0_GP)
      kz = 2.0_GP*pi/(nd(3)) !-1.0_GP)
      do k = ksta,kend
         z = real(k-1,kind=GP)
         do j = 1,n
           y = real(j-1,kind=GP)
           do i = 1,n
             x = real(i-1,kind=GP)
             v  (i,j,k) = sin(kx*x)*sin(ky*y)*sin(kz*z)
            enddo
         enddo
      enddo
      CALL gpspline%PartUpdate3D(px,py,pz,nparts)
      CALL gpspline%CompSpline3D(v,tmp1,tmp2)
      CALL gpspline%DoInterp3D(vlag,nparts)

      if ( myrank .eq. 0 ) then 
        write(*,*) 'Checking Splin_int...'
      endif
      suml = 0.0
      kmax = 1
      do k = 1, nparts
        valag(k) = sin(kx*px(k))*sin(ky*py(k))*sin(kz*pz(k))
        dv   (k) = abs(vlag(k)-valag(k))
        if ( suml .lt. dv(k) ) then
          suml = dv(k)
          kmax = k
        endif 
      enddo
      call MPI_ALLREDUCE(suml,sumg,1,GC_REAL,   &
                         MPI_MAX,MPI_COMM_WORLD,ierr)
write(*,*)'kmax =',kmax
write(*,*)myrank,': dv_kmax=',dv(kmax),' px=',px(kmax),' py=',py(kmax),' pz=',pz(kmax)
      if ( myrank .eq. 0 ) then 
        open(1,file='splerr.dat',position='append')
        write(1,*) n, sumg, nprocs, maxparts 
        close(1)
        if ( sumg .gt. 1.0e-4 ) then
          write(*,*) '............................................GSplineInt FAILURE'
          write(*,*) 'GSplineInt error=', sumg
          stop
        else
          write(*,*) '............................................GSplineInt passed'
        endif
      endif

      deallocate  ( dv,v,vb,ve,vlag,vt )
      deallocate                 ( bot )
      deallocate                 ( top )
      deallocate                 ( vdb )
      deallocate                ( ptmp )
      deallocate         ( id,px,py,pz )
      deallocate     ( gid,gpx,gpy,gpz )
      deallocate                 ( tmp1)

      call MPI_FINALIZE(ierr)

      end program lagtest
