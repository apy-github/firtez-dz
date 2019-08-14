!
MODULE DERIVVAR
  !
  !================================================
  !
  ! J M Borrero
  ! Sep 16, 2017
  ! KIS, Freiburg
  !
  USE CONS_PARAM, ONLY: DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP), PUBLIC  :: DNHYDDTEMP_PG, DNHYDDTEMP_RHO, DNHYDDPG_TEMP
  REAL(DP), PUBLIC  :: DNHYDDRHO_TEMP
  REAL(DP), PUBLIC  :: DDAMDTEMP_PG, DDAMDTEMP_RHO, DDAMDPG_TEMP
  REAL(DP), PUBLIC  :: DDAMDRHO_TEMP
  REAL(DP), PUBLIC  :: DNELECDPG_TEMP, DNELECDTEMP_RHO, DNELECDTEMP_PG
  REAL(DP), PUBLIC  :: DNELECDTEMP_PGAS, DNELECDRHO_TEMP
  REAL(DP), PUBLIC  :: DKLDTEMP_RHO, DKLDTEMP_PG, DKLDPG_TEMP, DKLDRHO_TEMP
  REAL(DP), PUBLIC  :: DKCDTEMP_RHO, DKCDTEMP_PG, DKCDPG_TEMP, DKCDRHO_TEMP
  REAL(DP), PUBLIC  :: DNIRDPG_TEMP, DNIRDRHO_TEMP
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
END MODULE DERIVVAR
!
