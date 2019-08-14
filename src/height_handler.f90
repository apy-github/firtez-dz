!
MODULE HEIGHT_HANDLER
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: TM_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: PG_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: P0_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: RH_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: BX_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: BY_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: BZ_COEFS
  REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:),TARGET    :: VZ_COEFS
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
END MODULE HEIGHT_HANDLER
!
