!
MODULE COUPLED_PARAM
  !
  USE CONS_PARAM, ONLY: DP, SP
  !
  IMPLICIT NONE
  !
  INTEGER                                   :: COU_BLCKSZ
  INTEGER                                   :: COU_NMTHRD
  INTEGER                                   :: COU_PSFRAD
  !
  INTEGER                                   :: COU_NPX
  INTEGER                                   :: COU_NPY
  !
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: COU_PSF
  CHARACTER*800                             :: COU_PSFFNAME
  !
  INTEGER                                   :: COU_NSPX
  INTEGER                                   :: COU_NSPY
  !
  INTEGER                                   :: COU_NPGRP
  !
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: COU_SPSF
  !
  ! BETA2DC vars:
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: COU_ALPHA
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: COU_BETA
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: COU_CH
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: COU_CHISQ
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: COU_OCHISQ
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: COU_OBS
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: COU_YMOD
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: COU_DYDA
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: COU_SIG
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: COU_ID
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: COU_IWL1
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: COU_IWLN
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: COU_NFACT
  !
  ! DA2DC vars:
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: COU_DA
  REAL(DP), DIMENSION(:,:), ALLOCATABLE     :: COU_SVDLAMBDA
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: COU_STEPS
  !REAL(DP)                                  :: COU_SVDLAMBDA
  !REAL(DP)                                  :: COU_STEPS
  REAL(DP)                                  :: COU_PRECYC
  ! Dummy var, it is an input but we have commented the lines of DA2DC...
  ! ... where it is used. Anyhow, it is yet an input.
  INTEGER, DIMENSION(:), ALLOCATABLE        :: COU_LISTA
  INTEGER, DIMENSION(:), ALLOCATABLE        :: COU_IPGRP
  !
  ! For MPI slaves:
  REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET       :: RCV_DA1D
  REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET       :: RCV_NFACT1D
  !
END MODULE COUPLED_PARAM
!
