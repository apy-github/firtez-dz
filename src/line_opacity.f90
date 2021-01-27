!
MODULE LINE_OPACITY
  !
  !================================================
  !
  ! A. Pastor Yabar & J M Borrero
  ! June 2017
  ! KIS Freiburg
  !
  USE CONS_PARAM, ONLY: DP, SP, DPI, EVOLT, HPLA, KBOL &
      , MELE, LIGHT, MAMU, QELE
  USE ATOM_DATABASE, ONLY: ABUND, MATOM, PARTITION_FUNCTION
  USE DERIVVAR
  USE CHEMICAL_EQUILIBRIUM, ONLY: ION_CALC_PGTEMP, ION_CALC_TEMPX
  USE CODE_MODES, ONLY: MRESPFUNCT
  USE FORWARD_PARAM, ONLY: LINE_ZN, LINE_ION, LINE_L0, EPLOW, LOGGF, POPU, POPL
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: OPAC_LINE
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! opac_line
  !
  !------------------------------------------------
  !
  ! Obsolete with nlte departure coefficients:SUBROUTINE OPAC_LINE(INDEX,TEMP,NELEC,NHYD,KLIN,DKLDT_RHO,DKLDT_PG,DKLDPG_TEM)
  SUBROUTINE OPAC_LINE(BL,BU,INDEX,TEMP,NELEC,NHYD,KLIN,DKLDT_RHO,DKLDT_PG,DKLDPG_TEM)
    !
    IMPLICIT NONE
    ! Input
    REAL(DP), INTENT(IN)               :: BL, BU
    !
    INTEGER, INTENT(IN)                :: INDEX
    REAL(SP), INTENT(IN)               :: TEMP
    REAL(DP), INTENT(IN)               :: NELEC, NHYD
    ! Output
    REAL(DP), INTENT(INOUT)              :: KLIN
    ! Optional output
    REAL(DP), INTENT(INOUT), OPTIONAL    :: DKLDT_RHO, DKLDT_PG, DKLDPG_TEM
    ! Internal
    REAL(DP)                           :: DNIRDT_PG, DNIRDT_RHO, N
    REAL(DP)                           :: TERM1, TERM2, TERM3, TERM4, TERM5
    REAL(DP)                           :: DTERM2, DTERM3, DTERM4, DTERM5
    REAL(DP)                           :: UI, UII, UIII, DUI, DUII, DUIII
    REAL(DP)                           :: U_PI, DU_PI
    REAL(DP)                           :: N_PI
    !
    REAL(DP)                            :: FI, FII, FIII
    REAL(DP)                            :: DFI, DFII, DFIII
    REAL(DP)                            :: FIR, DFIR
    !
    REAL(DP)                           :: T8
    !
    T8=DBLE(TEMP)
    ! Number of atoms per cm^3 of the species whose spectral line...
    ! ...we are interested in
    N=NHYD*10.0D0**(ABUND(LINE_ZN(INDEX))-12.0D0)
    CALL ION_CALC_TEMPX(INDEX,T8,NELEC,N,N_PI,DNIRDT_PG,DNIRDT_RHO)
    !
    ! N_PI IS THE DENSITY OF THE IONISATION OF THE SPECIE DIRECTLY!!!!!
    CALL PARTITION_FUNCTION(LINE_ZN(INDEX),T8,UI,UII,UIII,DUI,DUII,DUIII)
    !
    ! Note: ION_CALC_TEMPX already returns the ionisation
    !       state density as well as its derivatives with respect to 
    !       the temperature at constant RHO and PG. That is the 
    !       reason why this SELECT CASE estatement just deals
    !       with the partition functions and their derivatives.
    !
    SELECT CASE (LINE_ION(INDEX))
       CASE(1)
          U_PI = DBLE(UI)
          DU_PI = DBLE(DUI)
       CASE(2)
          U_PI = DBLE(UII)
          DU_PI = DBLE(DUII)
       CASE(3)
          U_PI = DBLE(UIII)
          DU_PI = DBLE(DUIII)
       CASE DEFAULT
          PRINT*, 'THE IONISATION STATE OF THE ELEMENTS '
          PRINT*, 'MUST BE BETWEEN I and III!!!'
          PRINT*, 'CURRENT VALUE: ', LINE_ION(INDEX)
          STOP
    END SELECT
    !
    ! KLIN itself
    !
    TERM1=SQRT(DPI)*QELE**2D0/(MELE*LIGHT)*10D0**(LOGGF(INDEX))&
        *LINE_L0(INDEX)*1D-8
    !
    TERM2=SQRT(MATOM(LINE_ZN(INDEX))*MAMU/(2D0*KBOL*T8))
    !
    ! Obsolete nlte dep. coef.:TERM3=N_PI/U_PI!(N_PI/U_PI)
    TERM3=BL*N_PI/U_PI
    !
    TERM4=DEXP(-EPLOW(INDEX)*EVOLT/(KBOL*T8))
    !
    ! Obsolete nlte dep. coef.:TERM5=1.0D0-DEXP(-HPLA*LIGHT/(1.0D-8*LINE_L0(INDEX)*KBOL*T8))
    TERM5=1.0D0-(BU/BL)*DEXP(-HPLA*LIGHT/(1.0D-8*LINE_L0(INDEX)*KBOL*T8))
    !
    KLIN=TERM1*TERM2*TERM3*TERM4*TERM5
    !
    ! DERIVATIVES CALCULATIONS:
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      !
      ! DNIRDPG_TEMP
      ! Calculate neutrals, single ionized and double ionized atoms...
      ! ...of species J
      CALL ION_CALC_PGTEMP(LINE_ZN(INDEX),T8,NELEC,DNELECDPG_TEMP&
        ,FI,FII,FIII,DFI,DFII,DFIII)
      !
      SELECT CASE (LINE_ION(INDEX))
         CASE(1)
            DFIR=DFI
            FIR=FI
         CASE(2)
            DFIR=DFII
            FIR=FII
         CASE(3)
            DFIR=DFIII
            FIR=FIII
         CASE DEFAULT
            PRINT*, 'THE IONISATION STATE OF THE ELEMENTS '
            PRINT*, 'MUST BE BETWEEN I and III!!!'
            PRINT*, 'CURRENT VALUE: ', LINE_ION(INDEX)
            STOP
      END SELECT
      DNIRDPG_TEMP=N/NHYD*FIR*DNHYDDPG_TEMP&
          +N*DFIR
      ! END TEST DNIRDPG_TEMP
!  I AM NO WRONG, THE ONLY DIFFERENCE BETWEEN d/dpg_t and d/drho_t should be the...
! . use of dne/drho_t instead of dne/dpg_t:
      ! DNIRDRHO_TEMP
      ! Calculate neutrals, single ionized and double ionized atoms...
      ! ...of species J
      CALL ION_CALC_PGTEMP(LINE_ZN(INDEX),T8,NELEC,DNELECDRHO_TEMP&
        ,FI,FII,FIII,DFI,DFII,DFIII)
      !
      SELECT CASE (LINE_ION(INDEX))
         CASE(1)
            DFIR=DFI
            FIR=FI
         CASE(2)
            DFIR=DFII
            FIR=FII
         CASE(3)
            DFIR=DFIII
            FIR=FIII
         CASE DEFAULT
            PRINT*, 'THE IONISATION STATE OF THE ELEMENTS '
            PRINT*, 'MUST BE BETWEEN I and III!!!'
            PRINT*, 'CURRENT VALUE: ', LINE_ION(INDEX)
            STOP
      END SELECT
      DNIRDRHO_TEMP=N/NHYD*FIR*DNHYDDRHO_TEMP&
          +N*DFIR
      ! 
      ! Derivatives of KL with T at constant PG and RHO
      !
      ! TERM1 is constant with T -> DTERM1 = 0.
      ! Derivative of TERM2 
      DTERM2=-0.5D0/TEMP*TERM2
      !
      ! Derivative of TERM3 at constant PG
      ! Obsolete nlte dep. coef.:DTERM3=(DNIRDT_PG-N_PI*DU_PI)/(U_PI)
      DTERM3=BL*(DNIRDT_PG-N_PI*DU_PI)/(U_PI)
      !
      ! Derivative of TERM4
      DTERM4=EPLOW(INDEX)*EVOLT/(KBOL*TEMP*TEMP)*TERM4
      !
      ! Derivative of TERM5
      DTERM5=-HPLA*LIGHT/(1.0D-8*LINE_L0(INDEX)*KBOL*T8*T8)*(1.0D0-TERM5)
      !
      DKLDTEMP_PG=TERM1*(&
           DTERM2*TERM3*TERM4*TERM5+&
           TERM2*DTERM3*TERM4*TERM5+&
           TERM2*TERM3*DTERM4*TERM5+&
           TERM2*TERM3*TERM4*DTERM5)
      !
      ! At constant RHO
      !
      ! Derivative of TERM3 at constant RHO
      ! Obsolete nlte dep. coef.:DTERM3=(DNIRDT_RHO-N_PI*DU_PI)/(U_PI)
      DTERM3=BL*(DNIRDT_RHO-N_PI*DU_PI)/(U_PI)
      !
      DKLDTEMP_RHO=TERM1*(&
           DTERM2*TERM3*TERM4*TERM5+&
           TERM2*DTERM3*TERM4*TERM5+&
           TERM2*TERM3*DTERM4*TERM5+&
           TERM2*TERM3*TERM4*DTERM5)
      !
      ! Derivative of KL with respect to gas pressure at constant temperature
      !
      DKLDPG_TEMP=DBLE(KLIN/N_PI*DNIRDPG_TEMP)
      !
      ! Derivative of KL with respect to gas pressure at constant temperature
      !
      DKLDRHO_TEMP=DBLE(KLIN/N_PI*DNIRDRHO_TEMP)
      ! Optional output
      IF(PRESENT(DKLDT_PG).AND.PRESENT(DKLDT_RHO).AND.PRESENT(DKLDPG_TEM)) THEN
         DKLDT_PG   = DKLDTEMP_PG
         DKLDT_RHO  = DKLDTEMP_RHO
         DKLDPG_TEM = DKLDPG_TEMP
      ENDIF
      !
      ! End derivatives of KLIN   
      !
    ENDIF
    !
  END SUBROUTINE OPAC_LINE
  !
  !================================================
  !
END MODULE LINE_OPACITY
!
