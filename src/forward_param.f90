!
MODULE FORWARD_PARAM
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: DP, SP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, PUBLIC                               :: NUML_DATABASE
  INTEGER, PUBLIC, ALLOCATABLE                  :: LINE_POS(:), LINE_NUM(:), LINE_ZN(:)
  INTEGER, PUBLIC, ALLOCATABLE                  :: LINE_ION(:), OETRANSITION(:,:)
  REAL(DP), PUBLIC, ALLOCATABLE                 :: SL(:), LL(:), JL(:), SU(:), LU(:), JU(:), LINE_L0(:)
  REAL(DP), PUBLIC, ALLOCATABLE                 :: LOGGF(:), EPLOW(:), ALPHA(:), SIGMA(:)
  INTEGER, PUBLIC                               :: NUML, NUMW
  INTEGER, PUBLIC, ALLOCATABLE                  :: IND_LINE(:), NUMWAVE(:)
  INTEGER, PUBLIC, ALLOCATABLE                  :: PIXEL_INI(:), PIXEL_END(:)
  REAL(DP), PUBLIC, ALLOCATABLE                 :: WAVE(:), DELTA_LAMBDA(:), WAVE_INI(:)
  REAL(DP), PUBLIC, ALLOCATABLE                 :: INDEX(:)
  REAL(DP), PUBLIC, ALLOCATABLE                 :: KC5(:), KC(:), KLIN(:,:)
  REAL(DP), PUBLIC, ALLOCATABLE                 :: KLINTAU(:,:)
  REAL(SP), PUBLIC, ALLOCATABLE                 :: TAULIN(:,:)
  !
  REAL(SP), PUBLIC, ALLOCATABLE                 :: GAIN1D(:,:)
  ! TEM, PG, RHO, PEL, MW, BX, BY. BZ, VX, VY, VZ, TAU5
  INTEGER, PUBLIC, PARAMETER                    :: ATM_ARGS = 13 
  ! D()/DT_RHO, D()/DPG_TEM, D()/DRHO_TEM, D()/DBX, D()/DBY, D()/DBZ, D()/DVZ
  INTEGER, PUBLIC, PARAMETER                    :: DER_ARGS = 8
  !
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:)  :: BLENDSID
  INTEGER, PUBLIC                               :: BLENDSMAX
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: BLENDSDIFF
  !
  INTEGER, PUBLIC                               :: N_FWD_MODELLING
  !
  LOGICAL, PUBLIC                               :: HYDRO_TOP
  !
  REAL(DP), PUBLIC,   ALLOCATABLE               :: LSF_SIGMA(:), LSF_W0(:), LSF_KERNEL(:)
  CHARACTER*800, PUBLIC, ALLOCATABLE            :: LSF_FILE(:)
  LOGICAL, PUBLIC, ALLOCATABLE, DIMENSION(:)    :: LSF_VALID
  !
  LOGICAL, PUBLIC                               :: FULL_STOKES
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
END MODULE FORWARD_PARAM
!
