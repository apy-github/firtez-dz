!
MODULE ATM_PARAM
  !
  !================================================
  !
  USE FORWARD_PARAM, ONLY: ATM_ARGS
  USE CONS_PARAM, ONLY: SP
  USE GRID_PARAM, ONLY: ZZ, NZ
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(SP), POINTER, PUBLIC   :: TEM(:)
  REAL(SP), POINTER, PUBLIC   :: PG(:)
  REAL(SP), POINTER, PUBLIC   :: RHO(:)
  REAL(SP), POINTER, PUBLIC   :: PEL(:)
  REAL(SP), POINTER, PUBLIC   :: MW(:)
  REAL(SP), POINTER, PUBLIC   :: BX(:)
  REAL(SP), POINTER, PUBLIC   :: BY(:)
  REAL(SP), POINTER, PUBLIC   :: BZ(:)
  REAL(SP), POINTER, PUBLIC   :: VX(:)
  REAL(SP), POINTER, PUBLIC   :: VY(:)
  REAL(SP), POINTER, PUBLIC   :: VZ(:)
  REAL(SP), POINTER, PUBLIC   :: TAU5(:)
  !
  PUBLIC :: SPLIT_MODEL
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! split_model
  !
  !------------------------------------------------
  !
  SUBROUTINE SPLIT_MODEL(MODEL)
    !
    REAL(SP), INTENT(IN), TARGET  :: MODEL(NZ,ATM_ARGS)
    INTEGER                          :: K
    !
    DO K=1,ATM_ARGS-1
       SELECT CASE (K)
       CASE(1)
          TEM => MODEL(:,K)
       CASE(2)
          PG => MODEL(:,K)
       CASE(3)
          RHO => MODEL(:,K)
       CASE(4)
          PEL => MODEL(:,K)
       CASE(5)
          MW => MODEL(:,K)
       CASE(6)
          BX => MODEL(:,K)
       CASE(7)
          BY => MODEL(:,K)
       CASE(8)
          BZ => MODEL(:,K)
       CASE(9)
          VX => MODEL(:,K)
       CASE(10)
          VY => MODEL(:,K)
       CASE(11)
          VZ => MODEL(:,K)
       CASE(12)
          TAU5 => MODEL(:,K)
       ENDSELECT
!print*, K, ':', SUM(MODEL(:,K))
    ENDDO
    ZZ(:)=MODEL(:,13)
!print*, 13, ':', SUM(MODEL(:,13))
    !
  END SUBROUTINE SPLIT_MODEL
  !
  !================================================
  !
END MODULE ATM_PARAM
!
