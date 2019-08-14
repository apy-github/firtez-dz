!
PROGRAM MAIN
  !
  !================================================
  !
  USE user_mpi, ONLY: START_MPI, END_MPI
  USE pre_post_duties, ONLY: START_STEP, END_STEP
  USE inversion, ONLY: INVERT3D
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  !
  !************************************************
  !
  CALL START_MPI()
  !
  CALL START_STEP()
  !
  CALL INVERT3D()
  !
  CALL END_STEP()
  !
  CALL END_MPI()
  !
  !================================================
  !
END PROGRAM MAIN
!
