!
MODULE INVERT_PARAM
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) :: INV_MAS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:)     :: INV_SLA
  !
  INTEGER, PUBLIC                                 :: MSGSIZE
  !
  INTEGER, PUBLIC                                 :: NFREQ, NSTKINV
  INTEGER, PUBLIC                                 :: NFREEV, NFREEP
  !
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE      :: NSLAB_PER_FREEV
  !
  INTEGER, PUBLIC                                 :: TFREEP, INFREEPMAX
  INTEGER, PUBLIC                                 :: CURIC
  INTEGER, PUBLIC                                 :: MAXITER
  INTEGER, PUBLIC                                 :: STOP_CRIT
  !
  REAL(SP), PUBLIC, ALLOCATABLE, DIMENSION(:)     :: YTOFIT, YGUESS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:)     :: ISIGMAP
  !
  REAL(SP), PUBLIC, ALLOCATABLE, DIMENSION(:)     :: MPERTURBATION
  !REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) :: DMODEL
  !
  LOGICAL, PUBLIC, DIMENSION(8)                   :: INV_ATMPAR
  !
  LOGICAL, PUBLIC, DIMENSION(4)                   :: INV_STK
  !
  INTEGER, PUBLIC, DIMENSION(8)                   :: NSLB_MAX
  !
  REAL(SP), PUBLIC, DIMENSION(4)                  :: WSTK
  !
  LOGICAL, PUBLIC                                 :: AUTOWEIGHT
  !
  REAL(DP), PUBLIC                                :: ISIGMA
  INTEGER, PUBLIC                                 :: CYCPOW
  INTEGER, PUBLIC                                 :: MAXSTEPS
  !
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:)    :: AM_I_DONE
  REAL(SP), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   :: IMASK
  !
  LOGICAL, PUBLIC                                 :: ASSIST_T, ASSIST_P, ASSIST_B, ASSIST_V
  LOGICAL, PUBLIC, DIMENSION(4)                   :: ASSIST_ATMPAR
  !
  REAL(DP), PUBLIC                                :: RSVDTOL, RSIGMAP
  !
  REAL(DP), PUBLIC, DIMENSION(8), TARGET    :: ATM_FACTOR
  !
   REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: JACOB
   REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:),TARGET   :: DELTA
   REAL(DP), PUBLIC                              :: LAMBDA
  !
   INTEGER, PUBLIC                               :: NJEVALS
   INTEGER, PUBLIC                               :: IFREEP
   INTEGER, PUBLIC                               :: SVDTOL
  !
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:)    :: INTERVALS
  !
   REAL(DP), PUBLIC                              :: MPERT
  !
  INTEGER, PUBLIC                                :: STEPS_WITHOUT_IMPROVEMENT
  REAL(DP), PUBLIC                               :: CHICUR, CHIPRE
  INTEGER, PUBLIC                                :: PRECYC
  !
  REAL(DP), PUBLIC                              :: INU
  INTEGER, PUBLIC                               :: TOFFSET
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  !------------------------------------------------
  !
  !
  !================================================
  !
END MODULE INVERT_PARAM
!
