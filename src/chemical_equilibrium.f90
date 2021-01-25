!
MODULE CHEMICAL_EQUILIBRIUM
  !
  !================================================
  !
  ! J M Borrero
  ! April 2, 2013
  ! KIS, Freiburg
  !
  USE CONS_PARAM, ONLY: SP, DP, KBOL, CSAHA1, CSAHA2 &
      , MELE, MAMU
  USE ATOM_DATABASE, ONLY: ABUND, XI, XII, NELEM &
      , DENO_SA, DENO_SAM, PARTITION_FUNCTION
  USE DERIVVAR
  USE FORWARD_PARAM, ONLY: LINE_ION, LINE_ZN
  !
  IMPLICIT NONE
  !
  REAL(DP), DIMENSION(NELEM)  :: VALPHAI,VALPHAII,VNT,VNI,VNII,VNIII
  REAL(DP), DIMENSION(NELEM)  :: VDALPHAI,VDALPHAII,VDNT,VDNI,VDNII,VDNIII
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: GET_RHO
  PUBLIC :: GET_DNHYDDTEMP
  PUBLIC :: GET_NHNEMU_G
  PUBLIC :: ION_CALC_TEMPRHO
  PUBLIC :: ION_CALC_PGTEMP
  PUBLIC :: ION_CALC_TEMPX
  !
  PRIVATE :: GET_PG
  PRIVATE :: GET_NE_RHO
  PRIVATE :: GET_NE_RHO_NIT
  PRIVATE :: GET_NE_RHO_NIT_SNE
  PRIVATE :: GET_NE_TEMP
  PRIVATE :: GET_NE_TEMP_SNE
  PRIVATE :: GET_NE_TEMP_SNE_ALT
  PRIVATE :: GET_MUTEM
  PRIVATE :: GET_NE_PG
  PRIVATE :: GET_NE_DNEDRHO_TEMP_NIT
  PRIVATE :: GET_NE_DNEDRHO_TEMP
  PRIVATE :: ION_CALC_NETEMP
  PRIVATE :: ION_CALC_NETEMP_NEW
  PRIVATE :: ION_CALC_TEMPPGAS
  PRIVATE :: ION_CALC_TEMPRHO_NIT
  PRIVATE :: GET_NHNEMU
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! get_pg
  ! get_rho
  ! get_ne_rho
  ! get_ne_rho_nit
  ! get_ne_rho_nit_sne
  ! get_ne_temp
  ! get_ne_temp_sne
  ! get_ne_temp_sne_alt
  ! get_mutem
  ! get_ne_pg
  ! get_ne_dnedrho_temp_nit
  ! get_ne_dnedrho_temp
  ! ion_calc_netemp
  ! ion_calc_netemp_new
  ! ion_calc_pgtemp
  ! ion_calc_tempx
  ! ion_calc_temprho
  ! ion_calc_temppgas
  ! ion_calc_temprho_nit
  ! get_nhnemu
  ! get_dnhyddtemp
  ! get_nhnemu_g
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_PG(TEMP,DENS,NHYD,NELEC,PATOM,PELEC,MOLECW)
    ! This routine obtains the total pressure (electrons+atoms)
    ! from the density and temperature
    !
    INTEGER                     :: I, J, ITER
    REAL(DP),   INTENT(IN)      :: TEMP, DENS
    REAL(DP),   INTENT(INOUT)     :: NHYD, NELEC, PELEC, PATOM, MOLECW
    REAL(DP)                    :: NT, NTI, NTII, NTIII
    REAL(DP)                    :: NELEC_NEW, NELEC_OLD, NELEC_EST, ERROR
    REAL(DP)                    :: NUME

    MOLECW=0D0
    PATOM=0D0
    PELEC=0D0
    !
    ! Estimation of electron density
    !
    IF ((NELEC.LT.1D0).OR.(NELEC.NE.NELEC)) NELEC = DENS/(150D0*MAMU)
    !
    ERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1E-4.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       ! Estimation of the hydrogen density
       NUME=DENS-NELEC_NEW*MELE
       NHYD=NUME/DENO_SAM
       NELEC_EST=0D0
       ! Initialize arrays
       NT=0D0
       NTI=0D0
       NTII=0D0
       NTIII=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NHYD*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms 
          ! of specie J
          CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,0D0,NT,NTI&
              ,NTII,NTIII)
          ! In the case of hydrogen (J=1) NTI corresponds to H- (electron 
          ! capture)
          IF (J.GT.1) NELEC_EST=NELEC_EST+NTII+NTIII
          IF (J.EQ.1) NELEC_EST=NELEC_EST-NTI+NTIII
       ENDDO
       !
       IF (NELEC_EST.LT.0.D0) NELEC_EST=0.D0
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       ITER = ITER + 1
    ENDDO
    !
    ! Calculate pressure due to electrons and due to atoms
    ! Calculate also mean molecuar weight
    !
    NELEC=NELEC_NEW
    PELEC=NELEC*KBOL*TEMP

    NUME=DENS-NELEC*MELE
    NHYD=NUME/DENO_SAM

    !
    PATOM=0.0D0
    DO I=1,NELEM
       PATOM=PATOM+10D0**(ABUND(I)-12D0)
    ENDDO
    PATOM=PATOM*NHYD*KBOL*TEMP
    MOLECW=((DENS*KBOL*TEMP)/MAMU)*(PELEC+PATOM)**(-1D0)
    !
  END SUBROUTINE GET_PG
  !
  !------------------------------------------------
  !








  SUBROUTINE GET_RHO(TEMP,PG,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW,GETDER)
    ! This routine obtains the density and electron/atom pressure
    ! from the total pressure and temperature
    !
    REAL(DP),   INTENT(IN)      :: TEMP, PG
    REAL(DP),   INTENT(INOUT)     :: NHYD, PELEC, PATOM, MOLECW, DENS
    REAL(DP),   INTENT(INOUT)     :: NELEC
    LOGICAL, OPTIONAL :: GETDER
    !
    INTEGER                     :: J, ITER
    REAL(DP)                    :: NELEC_NEW, NELEC_OLD, NELEC_EST, ERROR
    REAL(DP)                    :: NUME
    LOGICAL :: USEDER
    !
    USEDER=.FALSE.
    IF (PRESENT(GETDER)) USEDER=.TRUE.
    !
    DENS=0.0D0
    MOLECW=0.0D0
    PATOM=0.0D0
    NELEC_NEW=0.0D0
    NELEC_OLD=0.0D0
    NELEC_EST=0.0D0
    ERROR=0.0D0
    NUME=0.0D0
    NHYD=0.0D0
    !
    ! Estimation of electron density
    !
    IF ((NELEC.LT.1D0).OR.(NELEC.NE.NELEC)) NELEC=PG/(101D0*KBOL*TEMP)
    PELEC=NELEC*KBOL*TEMP
    !
    ! If our estimation for NELEC is way above the maximum possible...
    ! ...value, set it to give half of the pressure we have
    IF (PELEC.GT.PG) NELEC=(PG*0.5)/(KBOL*TEMP) 
    !
    ! Pre calculate some stuff needed to ...
    ! ...calculate neutrals, single ionized and double ionized atoms of
    !  all species
    CALL SHORT_PRE_ION_CALC_TEMPRHO(TEMP,NELEC,USEDER)
    !
    ERROR=1.0D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1.0D-4.AND.ITER.LE.100.AND.NHYD.GE.0.0D0)
       NELEC_OLD = NELEC_NEW
       ! Estimation of the hydrogen density
       NUME=PG-NELEC_NEW*KBOL*TEMP
       NHYD=NUME/(DENO_SA*KBOL*TEMP)
       NELEC_EST=0.0D0

       VNT(:)=NHYD*10.0D0**(ABUND(:)-12.0D0)
       CALL SHORT_ION_CALC_TEMPRHO(TEMP,NELEC_NEW,0.0D0,USEDER)

       ! In the case of hydrogen (J=1) NTI corresponds to H- 
       ! (electron capture)
       NELEC_EST=VNIII(1)-VNI(1)
       NELEC_EST=NELEC_EST+SUM(VNII(2:))+SUM(VNIII(2:))
       !
!!!       PELEC=NELEC_EST*KBOL*TEMP
!!!       IF (PELEC.GT.PG) THEN
!!!         PRINT*, NELEC_EST, PELEC, PG
!!!         NELEC_EST=(PG*0.9)/(KBOL*TEMP)
!!!       PELEC=NELEC_EST*KBOL*TEMP
!!!         PRINT*, NELEC_EST, PELEC, PG, (PG-NELEC_EST*KBOL*TEMP)/(DENO_SA*KBOL*TEMP)
!!!       ENDIF
       !
       NELEC_NEW=(MAX(NELEC_EST,0.0D0)+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       ITER = ITER + 1
    ENDDO
    !
    IF ( (ITER.GE.100) .OR. (NHYD.LT.0.D0) ) THEN
      PRINT*, ' Either ', NHYD, ' < 0.e0 or ', ITER, ' > 100'
      PRINT*, TEMP, PG, PATOM+PELEC, DENS, NELEC, NHYD
      PRINT*, PG-NELEC_NEW*KBOL*TEMP, NELEC_NEW*KBOL*TEMP, NELEC_EST, NELEC_NEW
      STOP
    ENDIF
    !
    ! Calculate pressure due to electrons and due to atoms
    !
    NELEC=NELEC_NEW
    !
    ! Calculate density
    NUME=PG-NELEC*KBOL*TEMP
    NHYD=NUME/(DENO_SA*KBOL*TEMP)
    DENS=NHYD*DENO_SAM+NELEC*MELE
    MOLECW=((DENS*KBOL*TEMP)/MAMU)*PG**(-1D0)
    !
    IF ((NHYD.NE.NHYD).OR.(NELEC_NEW.NE.NELEC_NEW)) THEN
      PRINT*, 'I am ce.f90, Either Nh or Ne are NaN!'
      PRINT*, TEMP, PG, NHYD, DENS, NELEC
    ENDIF
    IF ( (NELEC.LT.0) .OR. (NHYD.LE.0) ) THEN
      PRINT*, 'I am ce.f90, Either Nh or Ne are below 0!'
      PRINT*, TEMP, PG, PATOM+PELEC, DENS, NELEC, NHYD
      STOP
    ENDIF
!!!!!    IF (DENS.LT.0) THEN
!!!!!      PRINT*, ' *get_rho* '
!!!!!      PRINT*, ' Density is negative! '
!!!!!      PRINT*, 'TEMP', 'PG', 'NHYD', 'DENS', 'NELEC'
!!!!!      PRINT*, TEMP, PG, NHYD, DENS, NELEC
!!!!!      STOP
!!!!!    ENDIF

    !
  END SUBROUTINE GET_RHO
  !
  !------------------------------------------------
  !




  SUBROUTINE OLD_GET_RHO(TEMP,PG,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW,PRINTTRUE)
    ! This routine obtains the density and electron/atom pressure
    ! from the total pressure and temperature
    INTEGER                     :: J, ITER
    REAL(DP),   INTENT(IN)      :: TEMP, PG
    LOGICAL, OPTIONAL            :: PRINTTRUE
    REAL(DP),   INTENT(INOUT)     :: NHYD, PELEC, PATOM, MOLECW, DENS
    REAL(DP),   INTENT(INOUT)     :: NELEC
    REAL(DP)                    :: NT, NTI, NTII, NTIII
    REAL(DP)                    :: NELEC_NEW, NELEC_OLD, NELEC_EST, ERROR
    REAL(DP)                    :: NUME
    !
    DENS=0.0D0
    MOLECW=0.0D0
    PATOM=0.0D0
    PELEC=0.0D0
    NT=0.0D0
    NTI=0.0D0
    NTII=0.0D0
    NTIII=0.0D0
    NELEC_NEW=0.0D0
    NELEC_OLD=0.0D0
    NELEC_EST=0.0D0
    ERROR=0.0D0
    NUME=0.0D0
    !
    ! Estimation of electron density
    !
    IF ((NELEC.LT.1D0).OR.(NELEC.NE.NELEC)) NELEC = PG/(101D0*KBOL*TEMP)
    IF (PRESENT(PRINTTRUE)) WRITE(*,*) 'I am ce.f90', NELEC, PG, TEMP
    !
    ERROR=1.0D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1.0D-4.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       ! Estimation of the hydrogen density
       NUME=PG-NELEC_NEW*KBOL*TEMP
       NHYD=NUME/(DENO_SA*KBOL*TEMP)
       NELEC_EST=0.0D0
       ! Initialize arrays
       NT=0.0D0
       NTI=0.0D0
       NTII=0.0D0
       NTIII=0.0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NHYD*10.0D0**(ABUND(J)-12.0D0)
          ! Calculate neutrals, single ionized and double ionized atoms of
          !  specie J
          CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,0.0D0,NT,NTI&
              ,NTII,NTIII)
          ! In the case of hydrogen (J=1) NTI corresponds to H- 
          ! (electron capture)
          IF (J.GT.1) NELEC_EST=NELEC_EST+NTII+NTIII
          IF (J.EQ.1) NELEC_EST=NELEC_EST-NTI+NTIII
       ENDDO
       !
       IF (NELEC_EST.LT.0.D0) NELEC_EST=0.D0
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       ITER = ITER + 1
    ENDDO
    !
!PRINT*, 'NITER: ', ITER
    ! Calculate pressure due to electrons and due to atoms
    !
    NELEC=NELEC_NEW
    PELEC=NELEC*KBOL*TEMP
    !
    PATOM=DENO_SA*NHYD*KBOL*TEMP
    !
    ! Calculate density
    !
    NUME=PG-NELEC_NEW*KBOL*TEMP
    NHYD=NUME/(DENO_SA*KBOL*TEMP)
    DENS=NHYD*DENO_SAM+NELEC*MELE
    MOLECW=((DENS*KBOL*TEMP)/MAMU)*(PELEC+PATOM)**(-1D0)
    !
    IF ((NHYD.NE.NHYD).OR.(NELEC_NEW.NE.NELEC_NEW)) THEN
      PRINT*, 'I am ce.f90, Either Nh or Ne are NaN!'
      PRINT*, TEMP, PG, NHYD, DENS, NELEC
    ENDIF
    !
    IF (PRESENT(PRINTTRUE)) WRITE(*,*) 'I am ce.f90', NELEC, MOLECW, ITER, ERROR,NELEC_NEW,NELEC_OLD
    !
  END SUBROUTINE OLD_GET_RHO
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_RHO(TEMP,RHO,NE, DNE)
    ! Determines the electron pressure and its derivative with respect to T
    ! at constant density: DNELECDTEMP_RHO
    !
    ! Input
    REAL(DP), INTENT(IN)                  :: TEMP, RHO
    ! Output
    REAL(DP), INTENT(INOUT)               :: NE
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                               :: J, ITER
    REAL(DP)                              :: NH, ERROR, DERROR, NELEC&
        , NELEC_NEW, NELEC_OLD
    REAL(DP)                              :: NELEC_EST, DNELEC_NEW&
        , DNELEC_OLD, DNELEC_EST
    REAL(DP)                              :: NT, NTI, NTII, NTIII, DNTI&
        , DNTII, DNTIII
    !
    IF((NE.LT.1).OR.(NE.NE.NE)) THEN
      ! Determine NH
      NH=RHO/DENO_SAM
      ! Estimate number of electrons per cm^3
      NELEC=0.1*NH*DENO_SA
    ELSE
      NELEC=NE*1.D0
      NH=(RHO-NE*MELE)/DENO_SAM
    ENDIF
    !
    ERROR=1D0
    DERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    DNELEC_NEW=0D0
    DNELEC_OLD=0D0
    !
    ! Start loop until convergence
    !
    DO WHILE (DERROR.GT.1E-4.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       DNELEC_OLD = DNELEC_NEW
       !
       NELEC_EST=0D0
       DNELEC_EST=0D0
       ! Initialize arrays
       NT=0D0
       NTI=0D0
       NTII=0D0
       NTIII=0D0
       DNTI=0D0
       DNTII=0D0
       DNTIII=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms
          ! of specie J
          CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,DNELEC_NEW,NT&
              ,NTI,NTII,NTIII,DNTI,DNTII,DNTIII)
          ! In the case of hydrogen (J=1) NTI corresponds to H-
          ! (electron capture)
          IF (J.GT.1) THEN
             NELEC_EST=NELEC_EST+NTII+NTIII
             DNELEC_EST=DNELEC_EST+DNTII+DNTIII
          ENDIF
          IF (J.EQ.1) THEN
             NELEC_EST=NELEC_EST-NTI+NTIII
             DNELEC_EST=DNELEC_EST-DNTI+DNTIII
          ENDIF
       ENDDO
       !
       IF (NELEC_EST.LT.0.D0) NELEC_EST=0.D0
       !
       !NELEC_NEW=NELEC_OLD
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       !
       DNELEC_NEW=(DNELEC_EST+DNELEC_OLD)/2D0
       DERROR=ABS((DNELEC_NEW-DNELEC_OLD)/DNELEC_NEW)
       !
       ! UPDATE NH
       NH=(RHO-NELEC_NEW*MELE)/DENO_SAM
       !
       ITER = ITER + 1
       !
    ENDDO
!PRINT*, 'ITER #: ', ITER, ' ; ', NELEC, NELEC_NEW

    !
    NE=NELEC_NEW
    !
    ! Above calculation for Ne is NOT the one with the new dne:
    ! we repeat an additional iteration to make them consistent
    DNELEC_EST=0D0
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    DNTI=0D0
    DNTII=0D0
    DNTIII=0D0
    ! Loop through all elements
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms 
       ! of specie J
       CALL ION_CALC_TEMPRHO(J,TEMP,NELEC_NEW,DNELEC_NEW&
           ,NT,NTI,NTII,NTIII,DNTI,DNTII,DNTIII)
       ! In the case of hydrogen (J=1) NTI corresponds to H-
       ! (electron capture)
       IF (J.GT.1) THEN
          DNELEC_EST=DNELEC_EST+DNTII+DNTIII
       ENDIF
       IF (J.EQ.1) THEN
          DNELEC_EST=DNELEC_EST-DNTI+DNTIII
       ENDIF
    ENDDO
    !
    DNELECDTEMP_RHO=DNELEC_EST
    IF (PRESENT(DNE)) DNE=DNELECDTEMP_RHO
    !
  END SUBROUTINE GET_NE_RHO
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_RHO_NIT(TEMP,RHO,NE,DNE)
    !
    ! Determines the electron pressure and its derivative with respect to T
    ! at constant density: DNELECDTEMP_RHO
    !
    ! Input
    REAL(DP), INTENT(IN)                  :: TEMP, RHO
    REAL(DP), INTENT(IN)                  :: NE
    ! Output
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                               :: J
    REAL(DP)                              :: NH
    REAL(DP)                              :: NT, NTI, NTII, NTIII
    !
    REAL(DP)                              :: DNTI_N, DNTII_N, DNTIII_N
    REAL(DP)                              :: DNTI_D, DNTII_D, DNTIII_D
    !
    REAL(DP)                              :: DNT_N, DNT_D, DNELEC
    !
    NH=(RHO-NE*MELE)/DENO_SAM
    !
    ! Above calculation for Ne is NOT the one with the new dne:
    ! we repeat an additional iteration to make them consistent
    DNT_N=0.D0
    DNT_D=0.D0
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    DNTI_N=0D0
    DNTII_N=0D0
    DNTIII_N=0D0
    DNTI_D=0D0
    DNTII_D=0D0
    DNTIII_D=0D0
    ! Loop through all elements
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms 
       ! of specie J
       CALL ION_CALC_TEMPRHO_NIT(J,TEMP,NE,NT,NTI,NTII,NTIII&
           ,DNTI_N,DNTII_N,DNTIII_N,DNTI_D,DNTII_D,DNTIII_D)
       ! In the case of hydrogen (J=1) NTI corresponds to H-
       ! (electron capture)
       IF (J.GT.1) THEN
          DNT_N=DNT_N+DNTII_N+DNTIII_N
          DNT_D=DNT_D+DNTII_D+DNTIII_D
       ENDIF
       IF (J.EQ.1) THEN
          DNT_N=DNT_N-DNTI_N+DNTIII_N
          DNT_D=DNT_D-DNTI_D+DNTIII_D
       ENDIF
    ENDDO
    !
    DNELEC=DNT_N/(1.D0-DNT_D)
    !
    DNELECDTEMP_RHO=DNELEC
    IF (PRESENT(DNE)) DNE=DNELECDTEMP_RHO
    !
  END SUBROUTINE GET_NE_RHO_NIT
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_RHO_NIT_SNE(TEMP,RHO,NE,DNE)
    !
    ! Determines the electron pressure and its derivative with respect to T
    ! at constant density: DNELECDTEMP_RHO
    !
    ! Input
    REAL(DP), INTENT(IN)                  :: TEMP, RHO
    ! Output
    REAL(DP), INTENT(INOUT)               :: NE
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                               :: J
    REAL(DP)                              :: NH
    REAL(DP)                              :: NT, NTI, NTII, NTIII
    !
    REAL(DP)                              :: DNTI_N, DNTII_N, DNTIII_N
    REAL(DP)                              :: DNTI_D, DNTII_D, DNTIII_D
    !
    REAL(DP)                              :: DNT_N, DNT_D, DNELEC
    !
    NH=(RHO-NE*MELE)/DENO_SAM
    !
    ! Above calculation for Ne is NOT the one with the new dne:
    ! we repeat an additional iteration to make them consistent
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    DNTI_N=0D0
    DNTII_N=0D0
    DNTIII_N=0D0
    DNTI_D=0D0
    DNTII_D=0D0
    DNTIII_D=0D0
    ! Loop through all elements
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms 
       ! of specie J
       CALL ION_CALC_TEMPRHO_NIT(J,TEMP,NE,NT,NTI,NTII,NTIII&
           ,DNTI_N,DNTII_N,DNTIII_N,DNTI_D,DNTII_D,DNTIII_D)
       ! In the case of hydrogen (J=1) NTI corresponds to H-
       ! (electron capture)
       IF (J.GT.1) THEN
          DNT_N=DNT_N+DNTII_N+DNTIII_N
          DNT_D=DNT_D+DNTII_D+DNTIII_D
       ENDIF
       IF (J.EQ.1) THEN
          DNT_N=DNT_N-DNTI_N+DNTIII_N
          DNT_D=DNT_D-DNTI_D+DNTIII_D
       ENDIF
    ENDDO
    !
    DNELEC=DNT_N/(1.D0-DNT_D)
    !
    DNELECDTEMP_RHO=DNELEC
    IF (PRESENT(DNE)) DNE=DNELECDTEMP_RHO
    !
  END SUBROUTINE GET_NE_RHO_NIT_SNE
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_TEMP(TEMP,PGAS,L,NE,DNE)
    ! Determines the number of electrons per cm^3 and its derivative
    ! with respect to the pressure at constant temperature: DNELECDPG_TEMP
    ! Input
    REAL(DP), INTENT(IN)                :: TEMP, PGAS
    INTEGER, INTENT(IN)                 :: L
    ! Output
    REAL(DP), INTENT(INOUT)               :: NE
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                             :: J, ITER
    REAL(DP)                            :: NH, ERROR, NELEC, NELEC_NEW&
        , NELEC_OLD
    REAL(DP)                            :: NELEC_EST, NUM, DEN, DN1DNE
    REAL(DP)                            :: NT, NTI, NTII, NTIII
    REAL(DP)                            :: FI, FII, FIII, DFI, DFII, DFIII
    REAL(DP)                            :: FIR, DFIR, NTR
    REAL(DP)                            :: DNIDNH, DFIRDNH, DNIRDNE, DNIRDNH
    !
    ! Estimate number of electrons per cm^3
    NELEC=NE*1.D0
    IF ((NELEC.LT.1).OR.(NELEC.NE.NELEC)) NELEC=0.1*PGAS/(KBOL*TEMP)
    ! Estimate number of hydrogen atoms per cm^3
    NH=(PGAS-NELEC*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
    !
    ERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       NELEC_EST=0D0
       ! Initialize arrays
       NT=0D0
       NTI=0D0
       NTII=0D0
       NTIII=0D0
       FI=0D0
       FII=0D0
       FIII=0D0
       DFI=0D0
       DFII=0D0
       DFIII=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms
          ! of specie J
          CALL ION_CALC_NETEMP(J,TEMP,NELEC_NEW,FI,FII,FIII&
              ,DFI,DFII,DFIII)
          NTI=NT*FI
          NTII=NT*FII
          NTIII=NT*FIII
          ! In the case of hydrogen (J=1) NTI corresponds to H-
          ! (electron capture)
          IF (J.GT.1) THEN
             NELEC_EST=NELEC_EST+NTII+NTIII
          ENDIF
          IF (J.EQ.1) THEN
             NELEC_EST=NELEC_EST-NTI+NTIII
          ENDIF
       ENDDO
       !
       IF (NELEC_EST.LT.0.D0) NELEC_EST=0.D0
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       ! Update number of hydrogen atoms per cm^3
       NH=(PGAS-NELEC_NEW*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
       !
       ITER = ITER + 1
    ENDDO
!PRINT*, 'ITER #: ', ITER, ' ; ', NELEC, NELEC_NEW
    !
    NE=NELEC_NEW
    !
    NTR=0D0
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    FI=0D0
    FII=0D0
    FIII=0D0
    DFI=0D0
    DFII=0D0
    DFIII=0D0
    ! Loop through all elements
    DN1DNE=0D0
    NUM=0D0
    DEN=0D0
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms of
       ! specie J
       CALL ION_CALC_NETEMP(J,TEMP,NE,FI,FII,FIII&
           ,DFI,DFII,DFIII)
       NTI=NT*FI
       NTII=NT*FII
       NTIII=NT*FIII
       IF (J.EQ.1) THEN
          NUM=NUM+NT*(-DFI+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(-FI+FIII)
       ENDIF
       IF (J.GT.1) THEN
          NUM=NUM+NT*(DFII+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(FII+FIII)
       ENDIF
! FOR THE NEXT STEP IT IS CONVENIENT TO STORE SOME PARAMETERS:
       IF (J.EQ.LINE_ZN(L)) THEN
          NTR=NT
          SELECT CASE (LINE_ION(L))
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
                PRINT*, 'CURRENT VALUE: ', LINE_ION(L)
                STOP
          END SELECT
       ENDIF
    ENDDO
    !
    DN1DNE=(1D0-NUM)/DEN
    !
    DNELECDPG_TEMP=1D0/((DN1DNE*DENO_SA+1D0) * KBOL*TEMP)
    !
    IF (PRESENT(DNE)) DNE=DNELECDPG_TEMP
    !
    DNHYDDPG_TEMP=1D0/((DENO_SA + 1.D0/DN1DNE)*KBOL*TEMP)
    ! THERE EXISTS A VARIABLE THAT ALREADY IS USED FOR THIS! IT IS IN
    ! DERIVVAR MODULE! DNEDPG_TEMP=DNELECDPG_TEMP
    !
    ! We need the derivatives of DNIRDPG_TEMP
    DNIDNH=10.D0**(ABUND(LINE_ZN(L))-12.D0)
    DFIRDNH=DFIR/DN1DNE
    DNIRDNH=DNIDNH*FIR+NTR*DFIRDNH
    DNIRDNE=DN1DNE*10D0**(ABUND(LINE_ZN(L))-12D0)*FIR+NH*DNIDNH*DFIR
    !
    DNIRDPG_TEMP=1.D0/((DENO_SA/DNIRDNH+1.D0/DNIRDNE)*KBOL*TEMP)
    !
  END SUBROUTINE GET_NE_TEMP
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_TEMP_SNE(TEMP,PGAS,L,NE,DNE)
    ! Determines the number of electrons per cm^3 and its derivative
    ! with respect to the pressure at constant temperature: DNELECDPG_TEMP
    ! Input
    REAL(DP), INTENT(IN)                :: TEMP, PGAS
    INTEGER, INTENT(IN)                 :: L
    REAL(DP), INTENT(IN)                :: NE
    ! Output
    REAL(DP), INTENT(INOUT), OPTIONAL   :: DNE
    ! Internal
    INTEGER                             :: J
    REAL(DP)                            :: NH
    REAL(DP)                            :: NUM, DEN, DN1DNE
    REAL(DP)                            :: NT, NTI, NTII, NTIII
    REAL(DP)                            :: FI, FII, FIII, DFI, DFII, DFIII
    REAL(DP)                            :: FIR, DFIR, NTR
    REAL(DP)                            :: DNIDNH, DFIRDNH, DNIRDNE, DNIRDNH
    !
    NH=(PGAS-NE*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
    !
    !
    NTR=0D0
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    FI=0D0
    FII=0D0
    FIII=0D0
    DFI=0D0
    DFII=0D0
    DFIII=0D0
    ! Loop through all elements
    DN1DNE=0D0
    NUM=0D0
    DEN=0D0
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms of
       ! specie J
       CALL ION_CALC_NETEMP(J,TEMP,NE,FI,FII,FIII&
           ,DFI,DFII,DFIII)
       NTI=NT*FI
       NTII=NT*FII
       NTIII=NT*FIII
       IF (J.EQ.1) THEN
          NUM=NUM+NT*(-DFI+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(-FI+FIII)
       ENDIF
       IF (J.GT.1) THEN
          NUM=NUM+NT*(DFII+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(FII+FIII)
       ENDIF
! FOR THE NEXT STEP IT IS CONVENIENT TO STORE SOME PARAMETERS:
       IF (J.EQ.LINE_ZN(L)) THEN
          NTR=NT
          SELECT CASE (LINE_ION(L))
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
                PRINT*, 'CURRENT VALUE: ', LINE_ION(L)
                STOP
          END SELECT
       ENDIF
    ENDDO
    !
    DN1DNE=(1D0-NUM)/DEN
    !
    DNELECDPG_TEMP=1D0/((DN1DNE*DENO_SA+1D0) * KBOL*TEMP)
    !
    IF (PRESENT(DNE)) DNE=DNELECDPG_TEMP
    !
    DNHYDDPG_TEMP=1D0/((DENO_SA + 1.D0/DN1DNE)*KBOL*TEMP)
    ! THERE EXISTS A VARIABLE THAT ALREADY IS USED FOR THIS! IT IS IN
    ! DERIVVAR MODULE! DNEDPG_TEMP=DNELECDPG_TEMP
    !
    ! We need the derivatives of DNIRDPG_TEMP
    DNIDNH=10.D0**(ABUND(LINE_ZN(L))-12.D0)
    DFIRDNH=DFIR/DN1DNE
    DNIRDNH=DNIDNH*FIR+NTR*DFIRDNH
    DNIRDNE=DN1DNE*10D0**(ABUND(LINE_ZN(L))-12D0)*FIR+NH*DNIDNH*DFIR
    !
    DNIRDPG_TEMP=1.D0/((DENO_SA/DNIRDNH+1.D0/DNIRDNE)*KBOL*TEMP)
    !
  END SUBROUTINE GET_NE_TEMP_SNE
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_TEMP_SNE_ALT(TEMP,PGAS,NE,DNE)
    ! Determines the number of electrons per cm^3 and its derivative
    ! with respect to the pressure at constant temperature: DNELECDPG_TEMP
    ! Input
    REAL(DP), INTENT(IN)                :: TEMP, PGAS
    REAL(DP), INTENT(IN)                :: NE
    ! Output
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                             :: J
    REAL(DP)                            :: NH
    REAL(DP)                            :: NUM, DEN, DN1DNE
    REAL(DP)                            :: NT, NTI, NTII, NTIII
    REAL(DP)                            :: FI, FII, FIII, DFI, DFII, DFIII
    !
    REAL(DP)                              :: DNTI_N, DNTII_N, DNTIII_N
    REAL(DP)                              :: DNTI_D, DNTII_D, DNTIII_D
    !
    REAL(DP)                              :: DNT_N, DNT_D, DNELEC
    !
    ! Estimate number of hydrogen atoms per cm^3
    NH=(PGAS-NE*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
    !
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    FI=0D0
    FII=0D0
    FIII=0D0
    DFI=0D0
    DFII=0D0
    DFIII=0D0
    ! Loop through all elements
    DN1DNE=0D0
    NUM=0D0
    DEN=0D0
    !
    DNT_N=0.D0
    DNT_D=0.D0
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms 
       ! of specie J
       CALL ION_CALC_NETEMP_NEW(J,TEMP,NE,NT,NTI,NTII,NTIII&
           ,DNTI_N,DNTII_N,DNTIII_N,DNTI_D,DNTII_D,DNTIII_D)
       ! In the case of hydrogen (J=1) NTI corresponds to H-
       ! (electron capture)
       IF (J.GT.1) THEN
          DNT_N=DNT_N+DNTII_N+DNTIII_N
          DNT_D=DNT_D+DNTII_D+DNTIII_D
       ENDIF
       IF (J.EQ.1) THEN
          DNT_N=DNT_N-DNTI_N+DNTIII_N
          DNT_D=DNT_D-DNTI_D+DNTIII_D
       ENDIF
    ENDDO
    !
    DNELEC=DNT_N/(1.D0-DNT_D)
    !
    DNELECDPG_TEMP=DNELEC
    IF (PRESENT(DNE)) DNE=DNELECDPG_TEMP
    DNHYDDPG_TEMP=1.D0/DENO_SA*(1.D0/KBOL/TEMP-DNELECDPG_TEMP)
    !
  END SUBROUTINE GET_NE_TEMP_SNE_ALT
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_MUTEM(PG,DENS,TEMP,MU,NHYD,NELEC)
    ! Obtains temperature and molecular weight that
    ! are consistent with a given gas pressure and density
    REAL(DP), INTENT(IN)              :: PG, DENS
    REAL(DP), INTENT(INOUT)             :: TEMP, MU
    REAL(DP), INTENT(INOUT), OPTIONAL   :: NHYD, NELEC
    ! Internal
    REAL(DP)                          :: NH, NE
    REAL(DP)                          :: ERROR
    INTEGER                           :: ITER
    REAL(DP)                          :: DD
    REAL(DP)                          :: PATOM, PELEC
    ! Assume MU = 1 to begin with
    MU=1D0
    TEMP=0D0
    !
    ERROR=1E3
    ITER=1
    !
    DO WHILE (ERROR.GT.1E-6)
       TEMP=DBLE(PG)*MAMU*MU/(DENS*KBOL)
       CALL GET_RHO(TEMP,PG,NH,NE,PATOM,PELEC,DD,MU)
       ERROR=ABS((DENS-DD)/DENS)
       ITER=ITER+1
       IF (ITER.GE.100) THEN
          PRINT*,'No convergence achieved in GET_MUTEM after 100 iterations'
          STOP
       ENDIF
    ENDDO
    IF (PRESENT(NHYD).AND.PRESENT(NELEC)) THEN
       NELEC=NE
       NHYD=NH
    ENDIF
    !
  END SUBROUTINE GET_MUTEM
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_PG(TEMP,PGAS,NE,DNE)
    ! Determines the electron pressure and its derivative with respect to T
    ! at constant gas pressure: DNELECDTEMP_PG
    ! ABOVE SUBROUTINE ASSUMES DNh=0 YET, IT IS NOT
    !
    ! Input
    REAL(DP), INTENT(IN)                  :: TEMP, PGAS
    ! Output
    REAL(DP), INTENT(INOUT)                 :: NE
    REAL(DP), INTENT(INOUT), OPTIONAL       :: DNE
    ! Internal
    INTEGER                               :: J, ITER
    REAL(DP)                              :: NH, ERROR, DERROR, NELEC&
        , NELEC_NEW, NELEC_OLD
    REAL(DP)                              :: NELEC_EST, DNELEC_NEW&
        , DNELEC_OLD, DNELEC_EST
    REAL(DP)                              :: NT, NTI, NTII, NTIII, DNTI&
        , DNTII, DNTIII
    !
    ! Estimate number of electrons per cm^3
    NELEC=0.1*PGAS/(KBOL*TEMP)
    ! Estimate number of hydrogen atoms per cm^3
    NH=(PGAS-NELEC*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
    !
    ERROR=1D0
    DERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    DNELEC_NEW=0D0
    DNELEC_OLD=0D0
    !
    ! Start loop until convergence
    !
    DO WHILE (DERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       DNELEC_OLD = DNELEC_NEW
       !
       NELEC_EST=0D0
       DNELEC_EST=0D0
       ! Initialize arrays
       NT=0D0
       NTI=0D0
       NTII=0D0
       NTIII=0D0
       DNTI=0D0
       DNTII=0D0
       DNTIII=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms
          ! of specie J
          CALL ION_CALC_TEMPPGAS(J,TEMP,PGAS,NELEC_NEW,DNELEC_NEW,NT&
              ,NTI,NTII,NTIII,DNTI,DNTII,DNTIII)
          ! In the case of hydrogen (J=1) NTI corresponds to H-
          ! (electron capture)
          IF (J.GT.1) THEN
             NELEC_EST=NELEC_EST+NTII+NTIII
             DNELEC_EST=DNELEC_EST+DNTII+DNTIII
          ENDIF
          IF (J.EQ.1) THEN
             NELEC_EST=NELEC_EST-NTI+NTIII
             DNELEC_EST=DNELEC_EST-DNTI+DNTIII
          ENDIF
       ENDDO
       !
       IF (NELEC_EST.LT.0.D0) NELEC_EST=0.D0
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       !
       DNELEC_NEW=(DNELEC_EST+DNELEC_OLD)/2D0
       DERROR=ABS((DNELEC_NEW-DNELEC_OLD)/DNELEC_NEW)
       !
       ! Update number of hydrogen atoms per cm^3
       NH=(PGAS-NELEC_NEW*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
       !
       ITER = ITER + 1
       !
    ENDDO
    !
    NE=NELEC_NEW
    !
    ! Above calculation for Ne is NOT the one with the new dne:
    ! we repeat an additional iteration to make them consistent
    DNELEC_EST=0D0
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    DNTI=0D0
    DNTII=0D0
    DNTIII=0D0
    ! Loop through all elements
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms 
       ! of specie J
       CALL ION_CALC_TEMPPGAS(J,TEMP,PGAS,NELEC_NEW,DNELEC_NEW&
           ,NT,NTI,NTII,NTIII,DNTI,DNTII,DNTIII)
       ! In the case of hydrogen (J=1) NTI corresponds to H-
       ! (electron capture)
       IF (J.GT.1) THEN
          DNELEC_EST=DNELEC_EST+DNTII+DNTIII
       ENDIF
       IF (J.EQ.1) THEN
          DNELEC_EST=DNELEC_EST-DNTI+DNTIII
       ENDIF
    ENDDO
    !
    DNELECDTEMP_PGAS=DNELEC_EST
    DNELECDTEMP_PG=DNELEC_EST
    IF (PRESENT(DNE)) DNE=DNELECDTEMP_PGAS
    !
  END SUBROUTINE GET_NE_PG
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_DNEDRHO_TEMP_NIT(TEMP,RHO,NE,DNE)
    ! Determines the number of electrons per cm^3 and its derivative
    ! with respect to the pressure at constant temperature: DNELECDPG_TEMP
    ! Input
    REAL(DP), INTENT(IN)                :: TEMP, RHO
    REAL(DP), INTENT(IN)                :: NE
    ! Output
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                             :: J
    REAL(DP)                            :: NH
    REAL(DP)                            :: NUM, DEN
    REAL(DP)                            :: NT, NTI, NTII, NTIII
    REAL(DP)                            :: FI, FII, FIII
    REAL(DP)                            :: DFI, DFII, DFIII
    REAL(DP)                            :: DENOSUM
    !
    ! Determine NH
    NH=(RHO-NE*MELE)/DENO_SAM
    !
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    FI=0D0
    FII=0D0
    FIII=0D0
    DFI=0D0
    DFII=0D0
    DFIII=0D0
    ! Loop through all elements
    !
    DENOSUM=0D0
    NUM=0D0
    DEN=0D0
    !
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms of
       ! specie J
       CALL ION_CALC_NETEMP(J,TEMP,NE,FI,FII,FIII&
           ,DFI,DFII,DFIII)
       NTI=NT*FI
       NTII=NT*FII
       NTIII=NT*FIII
       !
       IF (J.EQ.1) THEN
          NUM=NUM+NT*(-DFI+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(-FI+FIII)
       ENDIF
       IF (J.GT.1) THEN
          NUM=NUM+NT*(DFII+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(FII+FIII)
       ENDIF
    ENDDO
    !
    DENOSUM=(1D0-NUM)/DEN
    !
    DNELECDRHO_TEMP=1D0/(DENOSUM*DENO_SAM+MELE)
    IF (PRESENT(DNE)) DNE=DNELECDRHO_TEMP
    !
    !DNHYDDRHO_TEMP=-MELE/DENO_SAM*DNELECDRHO_TEMP
    DNHYDDRHO_TEMP=(1D0-DNELECDRHO_TEMP*MELE)/DENO_SAM
    !
  END SUBROUTINE GET_NE_DNEDRHO_TEMP_NIT
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NE_DNEDRHO_TEMP(TEMP,RHO,NE,DNE)
    ! Determines the number of electrons per cm^3 and its derivative
    ! with respect to the pressure at constant temperature: DNELECDPG_TEMP
    ! Input
    REAL(DP), INTENT(IN)                :: TEMP, RHO
!    INTEGER, INTENT(IN)                 :: L
    ! Output
    REAL(DP), INTENT(INOUT)               :: NE
    REAL(DP), INTENT(INOUT), OPTIONAL     :: DNE
    ! Internal
    INTEGER                             :: J, ITER
    REAL(DP)                            :: NH, ERROR, NELEC, NELEC_NEW&
        , NELEC_OLD
    REAL(DP)                            :: NELEC_EST, NUM, DEN
    REAL(DP)                            :: NT, NTI, NTII, NTIII
    REAL(DP)                            :: FI, FII, FIII
    REAL(DP)                            :: DFI, DFII, DFIII
    REAL(DP)                            :: DENOSUM
    !
    ! Determine NH
    NH=RHO/DENO_SAM
    ! Estimate number of electrons per cm^3
    NELEC=0.1*NH*DENO_SA
    !
    ERROR=1D0
    ITER=1
    NELEC_NEW=NELEC
    NELEC_OLD=NELEC
    !
    ! Start loop until convergence
    !
    DO WHILE (ERROR.GT.1E-6.AND.ITER.LE.100)
       NELEC_OLD = NELEC_NEW
       NELEC_EST=0D0
       ! Initialize arrays
       NT=0D0
       NTI=0D0
       NTII=0D0
       NTIII=0D0
       FI=0D0
       FII=0D0
       FIII=0D0
       DFI=0D0
       DFII=0D0
       DFIII=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms
          ! of specie J
          CALL ION_CALC_NETEMP(J,TEMP,NELEC_NEW,FI,FII,FIII&
              ,DFI,DFII,DFIII)
          NTI=NT*FI
          NTII=NT*FII
          NTIII=NT*FIII
          ! In the case of hydrogen (J=1) NTI corresponds to H-
          ! (electron capture)
          IF (J.GT.1) THEN
             NELEC_EST=NELEC_EST+NTII+NTIII
          ENDIF
          IF (J.EQ.1) THEN
             NELEC_EST=NELEC_EST-NTI+NTIII
          ENDIF
       ENDDO
       !
       IF (NELEC_EST.LT.0.D0) NELEC_EST=0.D0
       !
       NELEC_NEW=(NELEC_EST+NELEC_OLD)/2D0
       ERROR=ABS((NELEC_NEW-NELEC_OLD)/NELEC_NEW)
       ! Update number of hydrogen atoms per cm^3
       NH=(RHO-NELEC_NEW*MELE)/DENO_SAM
       !
       ITER = ITER + 1
    ENDDO
    !
    NE=NELEC_NEW
    !
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    FI=0D0
    FII=0D0
    FIII=0D0
    DFI=0D0
    DFII=0D0
    DFIII=0D0
    ! Loop through all elements
    !
    DENOSUM=0D0
    NUM=0D0
    DEN=0D0
    !
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms of
       ! specie J
       CALL ION_CALC_NETEMP(J,TEMP,NE,FI,FII,FIII&
           ,DFI,DFII,DFIII)
       NTI=NT*FI
       NTII=NT*FII
       NTIII=NT*FIII
       !
       IF (J.EQ.1) THEN
          NUM=NUM+NT*(-DFI+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(-FI+FIII)
       ENDIF
       IF (J.GT.1) THEN
          NUM=NUM+NT*(DFII+DFIII)
          DEN=DEN+10D0**(ABUND(J)-12D0)*(FII+FIII)
       ENDIF
    ENDDO
    !
    DENOSUM=(1D0-NUM)/DEN
    !
    DNELECDRHO_TEMP=1D0/(DENOSUM*DENO_SAM+MELE)
    IF (PRESENT(DNE)) DNE=DNELECDRHO_TEMP
    !
    !DNHYDDRHO_TEMP=-MELE/DENO_SAM*DNELECDRHO_TEMP
    DNHYDDRHO_TEMP=(1D0-DNELECDRHO_TEMP*MELE)/DENO_SAM
    !
  END SUBROUTINE GET_NE_DNEDRHO_TEMP
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_NETEMP(I,TEMP,NE,FI,FII,FIII,DFI,DFII,DFIII)
    !
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(DP), INTENT(IN)            :: TEMP
    REAL(DP), INTENT(IN)            :: NE
    ! Output
    REAL(DP), INTENT(INOUT)           :: FI, FII, FIII
    REAL(DP), INTENT(INOUT), OPTIONAL :: DFI, DFII, DFIII
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    !
    REAL(DP)                        :: DENO
    !
    FI=0D0
    FII=0D0
    FIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    DENO=NE**2.D0*ALPHAI*ALPHAII+NE*ALPHAII+1.D0
    FI=NE**2.D0*ALPHAI*ALPHAII/DENO
    FII=NE*ALPHAII/DENO
    FIII=1.D0/DENO
    !
    ! Now derivatives of FI, FII and FIII with respect to the electron density at constant temperature
    !
    IF (PRESENT(DFI).AND.PRESENT(DFII).AND.PRESENT(DFIII)) THEN
       DFI=0D0
       DFII=0D0
       DFIII=0D0
       !
       DFI=FI/NE*2D0-FI*(NE*2D0*ALPHAI*ALPHAII+ALPHAII)/DENO
       DFII=ALPHAII/DENO-FII*(NE*2*ALPHAI*ALPHAII+ALPHAII)/DENO
       DFIII=-FIII*(NE*2*ALPHAI*ALPHAII+ALPHAII)/DENO
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_NETEMP
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_NETEMP_NEW(I,TEMP,NE,N,NI,NII,NIII&
      ,DNINUM,DNIINUM,DNIIINUM,DNIDEN,DNIIDEN,DNIIIDEN)
    !
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(DP), INTENT(IN)            :: TEMP
    REAL(DP), INTENT(IN)            :: NE, N
    ! Output
    REAL(DP), INTENT(INOUT)           :: NI, NII, NIII
    REAL(DP), INTENT(INOUT), OPTIONAL :: DNINUM, DNIINUM, DNIIINUM
    REAL(DP), INTENT(INOUT), OPTIONAL :: DNIDEN, DNIIDEN, DNIIIDEN
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    !
    REAL(DP)                        :: DENOI
    !
    NI=0D0
    NII=0D0
    NIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    NI=N/(1D0+(1D0/(NE*ALPHAI))*(1D0+(1D0/(NE*ALPHAII))))
    NII=N/(NE*ALPHAI+1D0+(1D0/(NE*ALPHAII)))
    NIII=N/(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))
    !
    ! Now derivatives of FI, FII and FIII with respect to the electron density at constant temperature
    !
    IF (PRESENT(DNINUM).AND.PRESENT(DNIIIDEN)) THEN
       !
       DENOI=((ALPHAI*NE+1.D0)*ALPHAII)*NE+1.D0
       !
       DNINUM=NI/N*10D0**(ABUND(I)-12D0)/DENO_SA/KBOL/TEMP
       DNIINUM=NII/N*10D0**(ABUND(I)-12D0)/DENO_SA/KBOL/TEMP
       DNIIINUM=NIII/N*10D0**(ABUND(I)-12D0)/DENO_SA/KBOL/TEMP
       !
       DNIDEN=-NI/N*10D0**(ABUND(I)-12D0)/DENO_SA&
           +N*(2.D0*ALPHAI*ALPHAII*NE)/DENOI&
           -NI*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI
       DNIIDEN=-NII/N*10D0**(ABUND(I)-12D0)/DENO_SA&
           +N*ALPHAII/DENOI&
           -NII*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI
       DNIIIDEN=-NIII/N*10D0**(ABUND(I)-12D0)/DENO_SA&
           +N*0.D0/DENOI&
           -NIII*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_NETEMP_NEW
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_PGTEMP(I,TEMP,NE,DNE,FI,FII,FIII,DFI,DFII,DFIII)
    !
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(DP), INTENT(IN)            :: TEMP
    REAL(DP), INTENT(IN)            :: NE,DNE
    ! Output
    REAL(DP), INTENT(INOUT)           :: FI, FII, FIII
    REAL(DP), INTENT(INOUT), OPTIONAL :: DFI, DFII, DFIII
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    !
    REAL(DP)                        :: DENO, DENOI
    !
    FI=0D0
    FII=0D0
    FIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*TEMP**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    DENO=NE**2.D0*ALPHAI*ALPHAII+NE*ALPHAII+1.D0
    FI=NE**2.D0*ALPHAI*ALPHAII/DENO
    FII=NE*ALPHAII/DENO
    FIII=1.D0/DENO
    !
    ! Now derivatives of FI, FII and FIII with respect to the electron density at constant temperature
    !
    IF (PRESENT(DFI).AND.PRESENT(DFII).AND.PRESENT(DFIII)) THEN
       DFI=0D0
       DFII=0D0
       DFIII=0D0
       !
       DENOI=((ALPHAI*NE+1.D0)*ALPHAII)*NE+1.D0
       !
       DFI=2.D0*ALPHAI*ALPHAII*NE*DNE/DENOI&
           -FI*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI*DNE
       DFII=ALPHAII*DNE/DENOI&
           -FII*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI*DNE
       DFIII=0.D0*DNE/DENOI&
           -FIII*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI*DNE
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_PGTEMP
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_TEMPX(INDEX,TEMP,NE,N,NI,DNIDT_PG,DNIDT_RHO)
    INTEGER,   INTENT(IN)             :: INDEX
    REAL(DP),  INTENT(IN)             :: TEMP
    REAL(DP),  INTENT(IN)             :: NE, N
    REAL(DP),  INTENT(INOUT)            :: NI
    REAL(DP),  INTENT(INOUT), OPTIONAL  :: DNIDT_PG, DNIDT_RHO
    ! Internal
    REAL(DP)                          :: UI, UII, UIII
    REAL(DP)                          :: DUI, DUII, DUIII
    REAL(DP)                          :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                          :: DLAMELEC, DEXPOI, DEXPOII, DALPHAI&
        , DALPHAII
    REAL(DP)                          :: DDENODTEMP_PG, DDENODTEMP_RHO
    REAL(DP)                          :: DNUMEDTEMP_PG, DNUMEDTEMP_RHO
    REAL(DP)                          :: DENO, NUME
    INTEGER                           :: ZN, ION
    !
    ZN=LINE_ZN(INDEX)
    ION=LINE_ION(INDEX)
    NI=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(ZN,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(ZN)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(ZN)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    DLAMELEC =-3D0/(2D0*TEMP)
    DEXPOI=-CSAHA2*XI(ZN)/TEMP**2D0
    DEXPOII=-CSAHA2*XII(ZN)/TEMP**2D0
    !--------------------
    ! Regular derivatives
    !--------------------
    DALPHAI=ALPHAI*(DUI-DUII+DLAMELEC+DEXPOI)
    DALPHAII=ALPHAII*(DUII-DUIII+DLAMELEC+DEXPOII)
    !
    DENO=NE**2D0*ALPHAI*ALPHAII+NE*ALPHAII+1D0
    DDENODTEMP_PG=2D0*NE*DNELECDTEMP_PG*ALPHAI*ALPHAII+NE**2D0&
        *DALPHAI*ALPHAII+NE**2D0*ALPHAI*DALPHAII+ &
        DNELECDTEMP_PG*ALPHAII+NE*DALPHAII
    DDENODTEMP_RHO=2D0*NE*DNELECDTEMP_RHO*ALPHAI*ALPHAII+NE**2D0&
        *DALPHAI*ALPHAII+NE**2D0*ALPHAI*DALPHAII+ &
        DNELECDTEMP_RHO*ALPHAII+NE*DALPHAII
    !
    SELECT CASE(ION)
       CASE(1)
          NUME=NE**2D0*ALPHAI*ALPHAII
          DNUMEDTEMP_PG=2D0*NE*DNELECDTEMP_PG*ALPHAI*ALPHAII+NE**2D0&
               *DALPHAI*ALPHAII+NE**2D0*ALPHAI*DALPHAII
          DNUMEDTEMP_RHO=2D0*NE*DNELECDTEMP_RHO*ALPHAI*ALPHAII+NE**2D0&
               *DALPHAI*ALPHAII+NE**2D0*ALPHAI*DALPHAII
       CASE(2)
          NUME=NE*ALPHAII
          DNUMEDTEMP_PG=DNELECDTEMP_PG*ALPHAII+NE*DALPHAII
          DNUMEDTEMP_RHO=DNELECDTEMP_RHO*ALPHAII+NE*DALPHAII
       CASE(3)
          NUME=1D0
          DNUMEDTEMP_PG=0D0
          DNUMEDTEMP_RHO=0D0
       CASE DEFAULT
         PRINT*, 'IONISATION LEVELS ABOVE 3 ARE NOT ALLOWED'
         STOP
    END SELECT
    !
    NI=N*NUME/DENO
    !
    IF (PRESENT(DNIDT_PG).AND.PRESENT(DNIDT_RHO)) THEN
       DNIDT_PG=0D0
       DNIDT_RHO=0D0
       DNIDT_PG=(NUME/DENO)*10D0**(ABUND(ZN)-12D0)*DNHYDDTEMP_PG&
           +N*(DNUMEDTEMP_PG*DENO-DDENODTEMP_PG*NUME)/DENO**2D0
       DNIDT_RHO=(NUME/DENO)*10D0**(ABUND(ZN)-12D0)*DNHYDDTEMP_RHO&
           +N*(DNUMEDTEMP_RHO*DENO-DDENODTEMP_RHO*NUME)/DENO**2D0
    ENDIF
    !
  END SUBROUTINE ION_CALC_TEMPX
  !
  !------------------------------------------------
  !
  SUBROUTINE SHORT_PRE_ION_CALC_TEMPRHO(TEMP,NE,DER)
    !
    ! Input
    REAL(DP), INTENT(IN)            :: TEMP, NE
    LOGICAL, INTENT(IN)             :: DER
    ! Internal
    INTEGER                         :: I
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII
    REAL(DP)                        :: DLAMELEC, DEXPOI, DEXPOII
    REAL(DP)                        :: DN
    !
    ! Loop through all elements
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3.0D0/2.0D0)      ! This has units of cm^3
    DO I=1,NELEM
      ! Get partition functions
      CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
      ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
      EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
      EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
      ! Check over/under flows
      IF (EXPOI .GT. 0.0D0) EXPOI=DMIN1(EXPOI,650.0D0)
      IF (EXPOI .LT. 0.0D0) EXPOI=DMAX1(EXPOI,-650.0D0)
      IF (EXPOII .GT. 0.0D0) EXPOII=DMIN1(EXPOII,650.0D0)
      IF (EXPOII .LT. 0.0D0) EXPOII=DMAX1(EXPOII,-650.0D0)
      !
      EXPOI = DEXP(EXPOI)
      EXPOII = DEXP(EXPOII)
      !
      VALPHAI(I)=(UI/(2.0D0*UII))*LAMELEC*EXPOI
      VALPHAII(I)=(UII/(2.0D0*UIII))*LAMELEC*EXPOII
      !
      !
      ! Now derivatives with respect to the temperature at constant density
      !
      IF (DER) THEN
        !
        DLAMELEC =-3D0/(2D0*TEMP)
        DEXPOI=-CSAHA2*XI(I)/TEMP**2D0
        DEXPOII=-CSAHA2*XII(I)/TEMP**2D0
        !--------------------
        ! Regular derivatives
        !--------------------
        VDALPHAI(I)=VALPHAI(I)*(DUI-DUII+DLAMELEC+DEXPOI)
        VDALPHAII(I)=VALPHAII(I)*(DUII-DUIII+DLAMELEC+DEXPOII)
        !
      ENDIF
    ENDDO
    !----
    !
  END SUBROUTINE SHORT_PRE_ION_CALC_TEMPRHO
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE SHORT_ION_CALC_TEMPRHO(TEMP,NE,DNE,DER)
    !
    ! Input
    REAL(DP), INTENT(IN)            :: TEMP, NE, DNE
    LOGICAL, INTENT(IN)             :: DER
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                        :: DLAMELEC, DEXPOI, DEXPOII&
        , DALPHAI, DALPHAII
    REAL(DP)                        :: DN
    !
    VNI(:)=VNT(:)/(1.0D0+(1.0D0/(NE*VALPHAI(:)))*(1.0D0+(1.0D0/(NE*VALPHAII(:)))))
    VNII(:)=VNT(:)/(NE*VALPHAI(:)+1.0D0+(1.0D0/(NE*VALPHAII(:))))
    VNIII(:)=VNT(:)/(1.0D0+NE*VALPHAII(:)*(NE*VALPHAI(:)+1.0D0))
    !
    ! Now derivatives with respect to the temperature at constant density
    !
    IF (DER) THEN
       !
       VDNT(:)=-MELE/DENO_SAM*DNE*10D0**(ABUND(:)-12D0)
       !
       VDNI(:)=VDNT(:)*VNI(:)/VNT(:)&
            -VNT(:)*(1D0+(1D0/(NE*VALPHAI(:)))*(1D0+1D0/(NE*VALPHAII(:))))**(-2D0)&
            *(1D0/(NE**2D0*VALPHAI(:)))*((-1D0/VALPHAI(:))*(DNE*VALPHAI(:)+NE*VDALPHAI(:))&
            *(1D0+1D0/(NE*VALPHAII(:)))-(1D0/(NE*VALPHAII(:)**2D0))&
            *(DNE*VALPHAII(:)+NE*VDALPHAII(:)))
       !
       VDNII(:)=VDNT(:)*VNII(:)/VNT(:)&
            -VNT(:)*(NE*VALPHAI(:)+1D0+1D0/(NE*VALPHAII(:)))**(-2D0)&
            *(DNE*VALPHAI(:)+NE*VDALPHAI(:)-(DNE*VALPHAII(:)+NE*VDALPHAII(:))&
            /(NE**2D0*VALPHAII(:)**2D0))
       !
       VDNIII(:)=VDNT(:)*VNIII(:)/VNT(:)&
            -VNT(:)*(1D0+NE*VALPHAII(:)*(NE*VALPHAI(:)+1D0))**(-2D0)&
            *(2D0*NE*DNE*VALPHAI(:)*VALPHAII(:)+NE**2D0*VDALPHAII(:)*VALPHAI(:)&
            +NE**2D0*VALPHAII(:)*VDALPHAI(:)+DNE*VALPHAII(:)+NE*VDALPHAII(:))
       !
    ENDIF
    !----
    !
  END SUBROUTINE SHORT_ION_CALC_TEMPRHO
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_TEMPRHO(I,TEMP,NE,DNE,N,NI,NII,NIII,DNI,DNII,DNIII)
    !
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(DP), INTENT(IN)            :: TEMP, NE, N, DNE
    ! Output
    REAL(DP), INTENT(INOUT)           :: NI, NII, NIII
    REAL(DP), INTENT(INOUT), OPTIONAL :: DNI, DNII, DNIII
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                        :: DLAMELEC, DEXPOI, DEXPOII&
        , DALPHAI, DALPHAII
    REAL(DP)                        :: DN
    !
    NI=0.0D0
    NII=0.0D0
    NIII=0.0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3.0D0/2.0D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0.0D0) EXPOI=DMIN1(EXPOI,650.0D0)
    IF (EXPOI .LT. 0.0D0) EXPOI=DMAX1(EXPOI,-650.0D0)
    IF (EXPOII .GT. 0.0D0) EXPOII=DMIN1(EXPOII,650.0D0)
    IF (EXPOII .LT. 0.0D0) EXPOII=DMAX1(EXPOII,-650.0D0)
    !
    EXPOI = DEXP(EXPOI)
    EXPOII = DEXP(EXPOII)
    !
    ALPHAI=(UI/(2.0D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2.0D0*UIII))*LAMELEC*EXPOII
    !
    NI=N/(1.0D0+(1.0D0/(NE*ALPHAI))*(1.0D0+(1.0D0/(NE*ALPHAII))))
    NII=N/(NE*ALPHAI+1.0D0+(1.0D0/(NE*ALPHAII)))
    NIII=N/(1.0D0+NE*ALPHAII*(NE*ALPHAI+1.0D0))
    !
    ! Now derivatives with respect to the temperature at constant density
    !
    IF (PRESENT(DNI).AND.PRESENT(DNII).AND.PRESENT(DNIII)) THEN
       !
       DN=-MELE/DENO_SAM*DNE*10D0**(ABUND(I)-12D0)
       !
       DNI=0D0
       DNII=0D0
       DNIII=0D0
       !
       DLAMELEC =-3D0/(2D0*TEMP)
       DEXPOI=-CSAHA2*XI(I)/TEMP**2D0
       DEXPOII=-CSAHA2*XII(I)/TEMP**2D0
       !--------------------
       ! Regular derivatives
       !--------------------
       DALPHAI=ALPHAI*(DUI-DUII+DLAMELEC+DEXPOI)
       DALPHAII=ALPHAII*(DUII-DUIII+DLAMELEC+DEXPOII)
!!!!!!!!!!       !
!!!!!!!!!!       DNI=DN*NI/N&
!!!!!!!!!!            -N*(1D0+(1D0/(NE*ALPHAI))*(1D0+1D0/(NE*ALPHAII)))**(-2D0)&
!!!!!!!!!!            *(1D0/(NE**2D0*ALPHAI))*((-1D0/ALPHAI)*(DNE*ALPHAI+NE*DALPHAI)&
!!!!!!!!!!            *(1D0+1D0/(NE*ALPHAII))-(1D0/(NE*ALPHAII**2D0))&
!!!!!!!!!!            *(DNE*ALPHAII+NE*DALPHAII))
!!!!!!!!!!       !
!!!!!!!!!!       DNII=DN*NII/N&
!!!!!!!!!!            -N*(NE*ALPHAI+1D0+1D0/(NE*ALPHAII))**(-2D0)&
!!!!!!!!!!            *(DNE*ALPHAI+NE*DALPHAI-(DNE*ALPHAII+NE*DALPHAII)&
!!!!!!!!!!            /(NE**2D0*ALPHAII**2D0))
!!!!!!!!!!       !
!!!!!!!!!!       DNIII=DN*NIII/N&
!!!!!!!!!!            -N*(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))**(-2D0)&
!!!!!!!!!!            *(2D0*NE*DNE*ALPHAI*ALPHAII+NE**2D0*DALPHAII*ALPHAI&
!!!!!!!!!!            +NE**2D0*ALPHAII*DALPHAI+DNE*ALPHAII+NE*DALPHAII)
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_TEMPRHO
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_TEMPPGAS(I,TEMP,PGAS,NE,DNE,N,NI,NII,NIII,DNI,DNII,DNIII)
    !
    ! Input
    INTEGER,  INTENT(IN)              :: I
    REAL(DP), INTENT(IN)              :: TEMP, PGAS
    REAL(DP), INTENT(IN)              :: NE, N, DNE
    ! Output
    REAL(DP), INTENT(INOUT)           :: NI, NII, NIII
    REAL(DP), INTENT(INOUT), OPTIONAL :: DNI, DNII, DNIII
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                        :: DLAMELEC, DEXPOI, DEXPOII&
        , DALPHAI, DALPHAII
    REAL(DP)                        :: DN
    !
    NI=0D0
    NII=0D0
    NIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    NI=N/(1D0+(1D0/(NE*ALPHAI))*(1D0+(1D0/(NE*ALPHAII))))
    NII=N/(NE*ALPHAI+1D0+(1D0/(NE*ALPHAII)))
    NIII=N/(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))
    !
    DN=-1.D0/DENO_SA*(PGAS/(KBOL*TEMP*TEMP)+DNE)*10D0**(ABUND(I)-12D0)
    !
    ! Now derivatives with respect to the temperature at constant density
    !
    IF (PRESENT(DNI).AND.PRESENT(DNII).AND.PRESENT(DNIII)) THEN
       DNI=0D0
       DNII=0D0
       DNIII=0D0
       !
       DLAMELEC =-3D0/(2D0*TEMP)
       DEXPOI=-CSAHA2*XI(I)/TEMP**2D0
       DEXPOII=-CSAHA2*XII(I)/TEMP**2D0
       !------------------------------------------------------------------------
       ! Derivatives considering that the derivatives of the partition functions
       ! are in fact the derivatives of the natural logarithm of the partition
       ! functions
       !------------------------------------------------------------------------
       DALPHAI=ALPHAI*(DUI-DUII+DLAMELEC+DEXPOI)
       DALPHAII=ALPHAII*(DUII-DUIII+DLAMELEC+DEXPOII)
       ! -----------------------------------------------------------------------
       DNI=DN*NI/N&
            -N*(1D0+(1D0/(NE*ALPHAI))*(1D0+1D0/(NE*ALPHAII)))**(-2D0)&
            *(1D0/(NE**2D0*ALPHAI))*((-1D0/ALPHAI)*(DNE*ALPHAI+NE*DALPHAI)&
            *(1D0+1D0/(NE*ALPHAII))-(1D0/(NE*ALPHAII**2D0))&
            *(DNE*ALPHAII+NE*DALPHAII))
       !
       DNII=DN*NII/N&
            -N*(NE*ALPHAI+1D0+1D0/(NE*ALPHAII))**(-2D0)&
            *(DNE*ALPHAI+NE*DALPHAI-(DNE*ALPHAII+NE*DALPHAII)&
            /(NE**2D0*ALPHAII**2D0))
       !
       DNIII=DN*NIII/N&
            -N*(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))**(-2D0)&
            *(2D0*NE*DNE*ALPHAI*ALPHAII+NE**2D0*DALPHAII*ALPHAI&
            +NE**2D0*ALPHAII*DALPHAI+DNE*ALPHAII+NE*DALPHAII)
       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_TEMPPGAS
  !
  !------------------------------------------------
  !
  SUBROUTINE ION_CALC_TEMPRHO_NIT(I,TEMP,NE,N,NI,NII,NIII&
      ,DNINUM,DNIINUM,DNIIINUM,DNIDEN,DNIIDEN,DNIIIDEN)
    !
    ! WARNING!! NOT SURE IF THIS IS CORRECT!!
    !
    ! This routine gives the thing inside the sums for the numerator and 
    ! denominator.
    ! 
    ! 
    !
    ! Input
    INTEGER,  INTENT(IN)            :: I
    REAL(DP), INTENT(IN)            :: TEMP
    REAL(DP), INTENT(IN)            :: NE, N
    ! Output
    REAL(DP), INTENT(INOUT)           :: NI, NII, NIII
    REAL(DP), INTENT(INOUT), OPTIONAL :: DNINUM, DNIINUM, DNIIINUM
    REAL(DP), INTENT(INOUT), OPTIONAL :: DNIDEN, DNIIDEN, DNIIIDEN
    ! Internal
    REAL(DP)                        :: UI, UII, UIII
    REAL(DP)                        :: DUI, DUII, DUIII
    REAL(DP)                        :: LAMELEC, EXPOI, EXPOII, ALPHAI, ALPHAII
    REAL(DP)                        :: DLAMELEC, DEXPOI, DEXPOII&
        , DALPHAI, DALPHAII
    REAL(DP)                        :: DN, DENOI
    !
    NI=0D0
    NII=0D0
    NIII=0D0
    ! Get partition functions
    CALL PARTITION_FUNCTION(I,TEMP,UI,UII,UIII,DUI,DUII,DUIII)
    ! ALPHAI and ALPHAII factors: these are the factors for the Saha equation
    LAMELEC = CSAHA1*DBLE(TEMP)**(-3D0/2D0)      ! This has units of cm^3
    EXPOI=CSAHA2*XI(I)/TEMP                      ! Dimensionless
    EXPOII=CSAHA2*XII(I)/TEMP                    ! Dimensionless
    ! Check over/under flows
    IF (EXPOI .GT. 0D0) EXPOI=DMIN1(EXPOI,650D0)
    IF (EXPOI .LT. 0D0) EXPOI=DMAX1(EXPOI,-650D0)
    IF (EXPOII .GT. 0D0) EXPOII=DMIN1(EXPOII,650D0)
    IF (EXPOII .LT. 0D0) EXPOII=DMAX1(EXPOII,-650D0)
    !
    EXPOI = EXP(EXPOI)
    EXPOII = EXP(EXPOII)
    !
    ALPHAI=(UI/(2D0*UII))*LAMELEC*EXPOI
    ALPHAII=(UII/(2D0*UIII))*LAMELEC*EXPOII
    !
    NI=N/(1D0+(1D0/(NE*ALPHAI))*(1D0+(1D0/(NE*ALPHAII))))
    NII=N/(NE*ALPHAI+1D0+(1D0/(NE*ALPHAII)))
    NIII=N/(1D0+NE*ALPHAII*(NE*ALPHAI+1D0))
    !
!    DN=-MELE/DENO_SAM*DNE*10D0**(ABUND(I)-12D0)
    DN=-MELE/DENO_SAM*10D0**(ABUND(I)-12D0)
    !
    ! Now derivatives with respect to the temperature at constant density
    !
    IF (PRESENT(DNINUM).AND.PRESENT(DNIIIDEN)) THEN
       DNINUM=0D0
       DNIINUM=0D0
       DNIIINUM=0D0
       DNIDEN=0D0
       DNIIDEN=0D0
       DNIIIDEN=0D0
       !
       DLAMELEC =-3D0/(2D0*TEMP)
       DEXPOI=-CSAHA2*XI(I)/TEMP**2D0
       DEXPOII=-CSAHA2*XII(I)/TEMP**2D0
       !--------------------
       ! Regular derivatives
       !--------------------
       DALPHAI=ALPHAI*(DUI-DUII+DLAMELEC+DEXPOI)
       DALPHAII=ALPHAII*(DUII-DUIII+DLAMELEC+DEXPOII)
       !
       DENOI=((ALPHAI*NE+1.D0)*ALPHAII)*NE+1.D0
       !
       DNINUM=N*NE*NE*(DALPHAI*ALPHAII+ALPHAI*DALPHAII)/DENOI&
           -NI*(DALPHAI*ALPHAII*NE+(ALPHAI*NE+1.D0)*DALPHAII)*NE/DENOI
       DNIINUM=N*NE*DALPHAII/DENOI&
           -NII*(DALPHAI*ALPHAII*NE+(ALPHAI*NE+1.D0)*DALPHAII)*NE/DENOI
       DNIIINUM=N*0.D0/DENOI&
           -NIII*(DALPHAI*ALPHAII*NE+(ALPHAI*NE+1.D0)*DALPHAII)*NE/DENOI
       DNIDEN=DN*NI/N&
           +N*(2.D0*ALPHAI*ALPHAII*NE)/DENOI&
           -NI*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI
       DNIIDEN=DN*NII/N&
           +N*ALPHAII/DENOI&
           -NII*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI
       DNIIIDEN=DN*NIII/N&
           +N*0.D0/DENOI&
           -NIII*(2.D0*ALPHAI*ALPHAII*NE+ALPHAII)/DENOI
!       !
    ENDIF
    !----
    !
  END SUBROUTINE ION_CALC_TEMPRHO_NIT
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NHNEMU(T,P,D,NH,NE,MU)
    !  This routine assumes temperature, gas pressure and density are known
    !  and calculates number of hydrogen atoms per cm^3, number of electrons
    ! per cm^3 and molecular weight.
    REAL(DP), INTENT(IN)     :: T, P, D
    REAL(DP), INTENT(INOUT)    :: NH, NE
    REAL(SP), INTENT(INOUT)    :: MU
    !
    REAL(DP)                 :: NHYD,NELEC,PATOM,PELEC,DENS,MOLECW
    !
    !NH=(DBLE(P)-DBLE(D)*KBOL*DBLE(T)/MELE)&
    !    /(KBOL*DBLE(T)*(DENO_SA-DENO_SAM/MELE))
    !NE=(DBLE(D)-NH*DENO_SAM)/MELE
    !MU=REAL(D*KBOL*T/(REAL(MAMU)*P))
    !
    CALL GET_RHO(T,P,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
    !
    NH=NHYD
    NE=NELEC
    MU=REAL(MOLECW)
    !
    IF (NE.LT.0) THEN
       PRINT*, ''
       PRINT*, 'I AM CHEMICAL EQUILIBRIM-'
       PRINT*, '  IN GET_NHNEMU SUBROUTINE:'
       PRINT*, 'ERROR!! ELECTRON NUMBER IS BELOW ZERO:', NE
       PRINT*, 'STOPPING'
       STOP
    ENDIF
    !
  END SUBROUTINE GET_NHNEMU
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_DNHYDDTEMP(T,P,D,M,NELEC,DNHDT_P,DNHDT_R)
    !
    REAL(SP), INTENT(IN)              :: T, P, D, M
    REAL(DP), INTENT(IN)              :: NELEC
    REAL(DP), INTENT(INOUT), OPTIONAL   :: DNHDT_P,DNHDT_R
    ! Internal
    REAL(DP)                          :: T8, P8, D8, M8
    !
    T8=DBLE(T)
    D8=DBLE(D)
    P8=DBLE(P)
    M8=DBLE(M)
    !
    ! This gives the derivative of ne with respect to PG at constant temperature
    !
    CALL GET_NE_TEMP_SNE_ALT(T8, P8, NELEC)
    !
    ! Now the derivative of ne with respect to T at constant density
    !
    CALL GET_NE_RHO_NIT(T8, D8, NELEC)
    !
    ! Derivative of ne with repect to RHO at constant temperature
    CALL GET_NE_DNEDRHO_TEMP_NIT(T8,D8,NELEC)
    !
    ! Derivative of ne with repect to T at constant gas pressure
    CALL GET_DNE_PG(T8,P8,NELEC)
    !
    ! Derivative of NHYD with repect to T at constant gas pressure
    DNHYDDTEMP_PG=(-1D0/DENO_SA)*(P8/(KBOL*T8**2D0)+DNELECDTEMP_PG)
    ! Derivative of NHYD with repect to T at constant density
    DNHYDDTEMP_RHO=(-MELE/DENO_SAM)*DNELECDTEMP_RHO
    ! Derivative of NHYD with repect to RHO at constant temperature
    !
    IF (PRESENT(DNHDT_P).AND.PRESENT(DNHDT_R)) THEN
       DNHDT_P=DNHYDDTEMP_PG
       DNHDT_R=DNHYDDTEMP_RHO
    ENDIF


!!!!!!!CALL CHECK_DNHYDDTEMP(T,P,D,M,NELEC)



    !
  END SUBROUTINE GET_DNHYDDTEMP






  !
  !------------------------------------------------
  !
  SUBROUTINE GET_DNE_PG(TEMP,PGAS,NE)
    ! Determines the electron pressure and its derivative with respect to T
    ! at constant gas pressure: DNELECDTEMP_PG
    ! ABOVE SUBROUTINE ASSUMES DNh=0 YET, IT IS NOT
    !
    ! Input
    REAL(DP), INTENT(IN)                  :: TEMP, PGAS, NE
    ! Internal
    INTEGER                               :: J, ITER
    REAL(DP)                              :: NH, DERROR
    REAL(DP)                              :: DNELEC_EST 
    REAL(DP)                              :: DNELEC_NEW 
    REAL(DP)                              :: DNELEC_OLD
    REAL(DP)                              :: NT, NTI, NTII, NTIII, DNTI&
        , DNTII, DNTIII
    !
    ! Estimate number of hydrogen atoms per cm^3
    NH=(PGAS-NE*KBOL*TEMP)/(KBOL*TEMP*DENO_SA)
    !
    DNELEC_OLD=0D0
    DNELEC_NEW=0D0
    ITER=0
    DERROR=1.0D20
    DO WHILE (DERROR.GT.1D-6.AND.ITER.LE.100)
       !
       DNELEC_OLD = DNELEC_NEW
       !
       DNELEC_EST=0D0
       ! Initialize arrays
       NT=0D0
       NTI=0D0
       NTII=0D0
       NTIII=0D0
       DNTI=0D0
       DNTII=0D0
       DNTIII=0D0
       ! Loop through all elements
       DO J=1,NELEM
          !
          NT=NH*10D0**(ABUND(J)-12D0)
          ! Calculate neutrals, single ionized and double ionized atoms
          ! of specie J
          CALL ION_CALC_TEMPPGAS(J,TEMP,PGAS,NE,DNELEC_NEW,NT&
              ,NTI,NTII,NTIII,DNTI,DNTII,DNTIII)
          ! In the case of hydrogen (J=1) NTI corresponds to H-
          ! (electron capture)
          IF (J.GT.1) THEN
             DNELEC_EST=DNELEC_EST+DNTII+DNTIII
          ENDIF
          IF (J.EQ.1) THEN
             DNELEC_EST=DNELEC_EST-DNTI+DNTIII
          ENDIF
       ENDDO
       !
       DNELEC_NEW=(DNELEC_EST+DNELEC_OLD)/2D0
       DERROR=ABS((DNELEC_NEW-DNELEC_OLD)/DNELEC_NEW)
       !
       ITER = ITER + 1
       !
    ENDDO
    !
    !
    ! Above calculation for Ne is NOT the one with the new dne:
    ! we repeat an additional iteration to make them consistent

    DNELEC_EST=0.0D0
    NT=0D0
    NTI=0D0
    NTII=0D0
    NTIII=0D0
    DNTI=0D0
    DNTII=0D0
    DNTIII=0D0
    ! Loop through all elements
    DO J=1,NELEM
       !
       NT=NH*10D0**(ABUND(J)-12D0)
       ! Calculate neutrals, single ionized and double ionized atoms 
       ! of specie J
       CALL ION_CALC_TEMPPGAS(J,TEMP,PGAS,NE,DNELEC_NEW&
           ,NT,NTI,NTII,NTIII,DNTI,DNTII,DNTIII)
       ! In the case of hydrogen (J=1) NTI corresponds to H-
       ! (electron capture)
       IF (J.GT.1) THEN
          DNELEC_EST=DNELEC_EST+DNTII+DNTIII
       ENDIF
       IF (J.EQ.1) THEN
          DNELEC_EST=DNELEC_EST-DNTI+DNTIII
       ENDIF
    ENDDO
    !
    DNELECDTEMP_PG=DNELEC_EST
    !
  END SUBROUTINE GET_DNE_PG
  !
  !------------------------------------------------
  !






  !
  SUBROUTINE CHECK_DNHYDDTEMP(T,P,D,M,NE)
    !
    REAL(SP), INTENT(IN)              :: T, P, D, M
    REAL(DP), INTENT(IN)              :: NE
    ! Internal
    REAL(DP)                          :: T8, P8, D8, M8
    REAL(DP)                          :: N1, N2, D1, D2
    REAL(DP) :: NHYD,PATOM,PELEC,DENS,MOLECW,NELEC,DNE
    !
    T8=DBLE(T)
    D8=DBLE(D)
    P8=DBLE(P)
    M8=DBLE(M)
    NELEC=NE
    !
    ! This gives the derivative of ne with respect to PG at constant temperature
    !
    CALL GET_NE_TEMP_SNE_ALT(T8, P8, NELEC)


D1=P8*1.01
    NELEC=NE
CALL GET_RHO(T8,D1,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, T8
N1=NHYD

D2=P8*0.99
    NELEC=NE
CALL GET_RHO(T8,D2,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, T8
N2=NHYD

PRINT*, ' * DNHYDDPG_TEMP * '
PRINT*, '     Analytical:', DNHYDDPG_TEMP
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'



D1=P8*1.01
    NELEC=NE
CALL GET_RHO(T8,D1,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, T8
N1=NELEC
    
D2=P8*0.99
    NELEC=NE
CALL GET_RHO(T8,D2,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, T8
N2=NELEC

PRINT*, ' * DNELECDPG_TEMP * '
PRINT*, '     Analytical:', DNELECDPG_TEMP
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'



    !
    ! Now the derivative of ne with respect to T at constant density
    !
    NELEC=NE
    CALL GET_NE_RHO_NIT(T8, D8, NELEC)

!!!!!!
!!!!!!DNELECDTEMP_RHO=DNELEC
!!!!!!

D1=T8*1.01
    NELEC=NE
CALL GET_PG(D1,D8,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, D8
N1=NELEC
    
D2=T8*0.99
    NELEC=NE
CALL GET_PG(D2,D8,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, D8
N2=NELEC

PRINT*, ' * DNELECDTEMP_RHO * '
PRINT*, '     Analytical:', DNELECDTEMP_RHO
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'


    !
    ! Derivative of ne with repect to RHO at constant temperature
    NELEC=NE
    CALL GET_NE_DNEDRHO_TEMP_NIT(T8,D8,NELEC)


D1=D8*1.01
    NELEC=NE
CALL GET_PG(T8,D1,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, T8
N1=NHYD

D2=D8*0.99
    NELEC=NE
CALL GET_PG(T8,D2,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, T8
N2=NHYD

PRINT*, ' * DNHYDDRHO_TEMP * '
PRINT*, '     Analytical:', DNHYDDRHO_TEMP
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'




D1=D8*1.01
    NELEC=NE
CALL GET_PG(T8,D1,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, T8
N1=NELEC

D2=D8*0.99
    NELEC=NE
CALL GET_PG(T8,D2,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, T8
N2=NELEC

PRINT*, ' * DNELECDRHO_TEMP * '
PRINT*, '     Analytical:', DNELECDRHO_TEMP
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'


    NELEC=NE
    CALL GET_DNE_PG(T8,P8,NELEC)

    !
! .TO BE REMOVED
    ! Derivative of NHYD with repect to T at constant gas pressure
    DNHYDDTEMP_PG=(-1D0/DENO_SA)*(P8/(KBOL*T8**2D0)+DNELECDTEMP_PG)
! TO BE REMOVED.
    ! Derivative of NHYD with repect to T at constant density
    DNHYDDTEMP_RHO=(-MELE/DENO_SAM)*DNELECDTEMP_RHO



D1=T8*1.01
    NELEC=NE
CALL GET_RHO(D1,P8,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, P8
N1=NELEC
    
D2=T8*0.99
    NELEC=NE
CALL GET_RHO(D2,P8,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, P8
N2=NELEC

PRINT*, ' * DNELECDTEMP_PG * '
PRINT*, '     Analytical:', DNELECDTEMP_PG
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'





D1=T8*1.01
    NELEC=NE
CALL GET_RHO(D1,P8,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, P8
N1=NHYD
    
D2=T8*0.99
    NELEC=NE
CALL GET_RHO(D2,P8,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
PRINT*, P8
N2=NHYD

PRINT*, ' * DNHYDDTEMP_PG * '
PRINT*, '     Analytical:', DNHYDDTEMP_PG
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'




D1=T8*1.01
    NELEC=NE
CALL GET_PG(D1,D8,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, D8
N1=NHYD
    
D2=T8*0.99
    NELEC=NE
CALL GET_PG(D2,D8,NHYD,NELEC,PATOM,PELEC,MOLECW)
PRINT*, D8
N2=NHYD

PRINT*, ' * DNHYDDTEMP_RHO * '
PRINT*, '     Analytical:', DNHYDDTEMP_RHO
PRINT*, '      Numerical:', (N1-N2)/(D1-D2)
print*, '-------------------------------------------------'







STOP
STOP
STOP
STOP
STOP
STOP
STOP
STOP



  END SUBROUTINE CHECK_DNHYDDTEMP
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_NHNEMU_G(T,P,D,NH,NE,MU,PEL)
    !
    USE INVERT_PARAM, ONLY: INV_ATMPAR
    USE CODE_MODES, ONLY: MRESPFUNCT
    !
    REAL(SP), INTENT(INOUT)     :: T, P, D, MU, PEL
    REAL(DP), INTENT(INOUT)    :: NH, NE
    !
    REAL(DP)                         :: PATOM, PELEC
    REAL(DP)                         :: T8, D8, P8, MU8
    !
    T8=DBLE(T)
    D8=DBLE(D)
    P8=DBLE(P)
    MU8=DBLE(MU)
    PATOM=0.D0
    PELEC=0.D0
    !

    IF (INV_ATMPAR(3).EQV..FALSE.) THEN
      CALL GET_RHO(T8,P8,NH,NE,PATOM,PELEC,D8,MU8)
      D=REAL(D8)
!PRINT*, T8,P8,D8,D
    ELSE IF (INV_ATMPAR(2).EQV..FALSE.) THEN
      CALL GET_PG(T8,D8,NH,NE,PATOM,PELEC,MU8)
      P=REAL(PATOM+PELEC)
    ELSE IF (INV_ATMPAR(1).EQV..FALSE.) THEN
      CALL GET_MUTEM(P8,D8,T8,MU8,NH,NE)
      T=REAL(T8)
    ELSE
      CALL GET_RHO(T8,P8,NH,NE,PATOM,PELEC,D8,MU8)
      D=REAL(D8)
    ENDIF

    MU=REAL(MU8)
    PEL=REAL(PELEC)
    !
    IF (NE.LT.0) THEN
       PRINT*, ''
       PRINT*, 'I AM CHEMICAL EQUILIBRIM-'
       PRINT*, '  IN GET_NHNEMU SUBROUTINE:'
       PRINT*, 'ERROR!! ELECTRON NUMBER IS BELOW ZERO:', NE, T8,P8,NH
       PRINT*, 'STOPPING'
       STOP
    ENDIF
    !
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      CALL GET_DNHYDDTEMP(T,P,D,MU,NE)
    ENDIF
    !
  END SUBROUTINE GET_NHNEMU_G
  !
  !================================================
  !
END MODULE CHEMICAL_EQUILIBRIUM
!
