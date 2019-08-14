!
MODULE PHYS_PARAM
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  ! PHYS1D_PARAM
  !
  REAL(SP), ALLOCATABLE, PUBLIC               :: OBS1D(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC               :: DSYN1D(:,:,:,:), SYN1D(:,:)
  REAL(SP), ALLOCATABLE, PUBLIC               :: FS_SYN1D(:,:,:)
  REAL(DP), ALLOCATABLE, PUBLIC               :: RSYN1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC               :: BEST_RSYN1D(:)
  REAL(DP), ALLOCATABLE, PUBLIC               :: ROBS1D(:)
  REAL(SP), ALLOCATABLE, PUBLIC               :: EQVDSYN(:,:)
  REAL(SP), ALLOCATABLE, PUBLIC               :: BEST_EQVDSYN(:,:)
  REAL(SP), ALLOCATABLE, PUBLIC, TARGET       :: MODEL1D_RCV(:,:)
  REAL(SP), ALLOCATABLE, PUBLIC, TARGET       :: MODEL2D_RCV(:,:,:)
  !
  REAL(DP), ALLOCATABLE, PUBLIC               :: EVOLG(:,:,:,:)
  !
  !INTEGER, ALLOCATABLE, PUBLIC                :: CALC_RFSP(:)
  !
  ! PHYS3D_PARAM
  !
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: TEM3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BX3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BY3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BZ3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: RHO3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: VX3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: VY3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: VZ3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: PG3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: MW3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: PEL3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: TAU3D5(:,:,:)
  !
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: ETEM3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EBX3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EBY3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EBZ3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: ERHO3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EVX3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EVY3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EVZ3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EPG3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EMW3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: EPEL3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: ETAU3D5(:,:,:)
  !
  REAL(SP), PUBLIC, ALLOCATABLE            :: OBS(:,:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE            :: MODEL1D_SND(:,:)
  REAL(SP), PUBLIC, ALLOCATABLE            :: MODEL2D_SND(:,:,:)
  REAL(DP), PUBLIC, ALLOCATABLE            :: OBS3D(:,:,:)
  !
  ! If Non-LTE
  ! REAL(SP),    ALLOCATABLE    :: SYN(:,:,:,:,:)
  ! If LTE
  REAL(SP), PUBLIC, ALLOCATABLE            :: SYN(:,:,:,:)
  REAL(DP), PUBLIC, ALLOCATABLE            :: SYN3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE            :: SYN5D(:,:,:,:,:)
  REAL(DP), PUBLIC, ALLOCATABLE            :: ISIGMAP3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE            :: DSYN(:,:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE            :: WDSYN(:,:,:,:,:,:)
  !
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_TEM3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_BX3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_BY3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_BZ3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_RHO3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_VX3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_VY3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_VZ3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_PG3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_MW3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_PEL3D(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE, TARGET    :: BEST_TAU3D5(:,:,:)
  !
  REAL(DP), PUBLIC, ALLOCATABLE            :: BEST_SYN(:,:,:)
  REAL(SP), PUBLIC, ALLOCATABLE            :: BEST_DSYN(:,:,:,:)
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
END MODULE PHYS_PARAM
!
