!
MODULE TIME_PARAM
  !
  !================================================
  !
  ! J M Borrero
  ! Jan 8, 2007
  ! HAO-NCAR for HMI-Stanford
  !
  USE CONS_PARAM, ONLY: SP, DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, PUBLIC               :: TSTART, TEND, TRATE, TMAX
  REAL(SP), PUBLIC              :: TIMEUSED
  REAL(DP), PUBLIC              :: DSTART
  REAL(DP), PUBLIC              :: TIME1, TIME2
  REAL(DP), PUBLIC              :: TIME1B, TIME2B
  REAL(DP), PUBLIC              :: TIME1C, TIME2C
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
END MODULE TIME_PARAM
!
