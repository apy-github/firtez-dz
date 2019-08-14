!
MODULE GRID_PARAM
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: DP, SP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, PUBLIC                  :: NX, NY, NZ
  REAL(DP), PUBLIC                 :: DX, DY, DZ
  !REAL(SP), ALLOCATABLE    :: XX(:,:,:), YY(:,:,:), ZZ(:,:,:)
  REAL(DP), PUBLIC, ALLOCATABLE    :: XX(:), YY(:), ZZ(:)
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
END MODULE GRID_PARAM
!
