!
SUBROUTINE JOIN_MODEL(I,J)
  !
  USE CONS_PARAM, ONLY: SP
  USE GRID_PARAM, ONLY: ZZ
  USE PHYS_PARAM
  USE FORWARD_PARAM, ONLY: ATM_ARGS
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, INTENT(IN)   :: I, J
  INTEGER               :: K
  REAL(SP), POINTER     :: ARR(:,:,:)
  !
  DO K=1,ATM_ARGS-1
     SELECT CASE (K)
     CASE(1)
        ARR => TEM3D
     CASE(2)
        ARR => PG3D
     CASE(3)
        ARR => RHO3D
     CASE(4)
        ARR => PEL3D
     CASE(5)
        ARR => MW3D
     CASE(6)
        ARR => BX3D
     CASE(7)
        ARR => BY3D
     CASE(8)
        ARR => BZ3D
     CASE(9)
        ARR => VX3D
     CASE(10)
        ARR => VY3D
     CASE(11)
        ARR => VZ3D
     CASE(12)
        ARR => TAU3D5
     ENDSELECT
     MODEL1D_SND(:,K) = ARR (:,J,I)
  !
     NULLIFY(ARR)
  ENDDO
  MODEL1D_SND(:,13) = ZZ(:)
  !
END SUBROUTINE JOIN_MODEL
!
