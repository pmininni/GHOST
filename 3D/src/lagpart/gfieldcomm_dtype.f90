!=================================================================
! GHOST GPartComm field communicator subclass
!
! 2023 F. Pugliese
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!=================================================================
   INTEGER  :: UNPACK_REP = 0, UNPACK_SUM = 1

   TYPE, PUBLIC, EXTENDS ( GPartComm ) :: GFieldComm
        PRIVATE
        ! Member data
        INTEGER                                :: nbret_ ,ntret_
        INTEGER, ALLOCATABLE, DIMENSION(:)     :: ibretp_,itretp_,ibretnz_ ,itretnz_
        INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: ibret_ ,itret_ ,ibretdst_,itretdst_
      CONTAINS
        ! Public methods:
        PROCEDURE,PUBLIC :: GFieldComm_ctor
        FINAL            :: GFieldComm_dtor
        PROCEDURE,PUBLIC :: SlabDataReturnMF => GFieldComm_SlabDataExchangeMF
        PROCEDURE,PUBLIC :: SlabDataReturnSF => GFieldComm_SlabDataExchangeSF
      END TYPE GFieldComm
