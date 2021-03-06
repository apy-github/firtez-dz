!
MODULE RTESOLVER
  !
  !================================================
  !
  ! May 22, 2014
  ! KIS, Freiburg
  !
  USE CONS_PARAM, ONLY: DP, HPLA, LIGHT, KBOL, IMAT
  USE ATM_PARAM, ONLY: TEM, PG, RHO
  USE ABSORPTION_MATRIX
  USE PHYS_PARAM, ONLY: SYN1D, DSYN1D, EVOLG
  USE CODE_MODES, ONLY: MRESPFUNCT, HYDROSTATIC, HYDROSTATICDER
  USE GRID_PARAM, ONLY: ZZ, NZ
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP), DIMENSION(4,4)                  :: EVOL
  REAL(DP), DIMENSION(4,4)                  :: DEVOLDVLOS, DEVOLDT_PG, DEVOLDT_RHO
  REAL(DP), DIMENSION(4,4)                  :: DEVOLDBX, DEVOLDBY, DEVOLDBZ
  REAL(DP), DIMENSION(4,4)                  :: DEVOLDPG_TEMP, DEVOLDRHO_TEMP
  REAL(DP), DIMENSION(:), ALLOCATABLE :: DPGDTK_HE
  !
  PUBLIC :: SET_BOUNDARY
  PUBLIC :: ANALYTICAL_SOLVER
  PUBLIC :: OLD_ANALYTICAL_SOLVER
  PUBLIC :: DPGDTK_HE
  !
  PRIVATE :: EVOL_OPERATOR
  PRIVATE :: SOURCE_FUNCTION
  PRIVATE :: SOURCE_FUNCTION_DER
  !
  PRIVATE :: EVOL
  PRIVATE :: DEVOLDVLOS, DEVOLDT_PG, DEVOLDT_RHO
  PRIVATE :: DEVOLDBX, DEVOLDBY, DEVOLDBZ
  PRIVATE :: DEVOLDPG_TEMP, DEVOLDRHO_TEMP
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! analytical_solver
  ! evol_operator
  ! source_function
  ! source_function_der
  ! set_boundary
  !
  !------------------------------------------------
  !
  SUBROUTINE ANALYTICAL_SOLVER(K,L)
    !
    USE INVERT_PARAM, ONLY: INV_STK, INV_ATMPAR
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)       :: K,L
    REAL(DP), DIMENSION(4)       :: SF
    REAL(DP), DIMENSION(4,1)     :: SFM, IM, IN
    REAL(DP)                     :: DSF
    ! For derivatives
    REAL(DP), DIMENSION(4,1)     :: RTE_DT_PG, RTE_DT_RHO, RTE_DPG_T, RTE_DRHO_T
    REAL(DP), DIMENSION(4,1)     :: RTE_DBX, RTE_DBY, RTE_DBZ
    REAL(DP), DIMENSION(4,1)     :: DSFM
    !
    INTEGER :: PCNT
    !
    CALL EVOL_OPERATOR(K,L)
    !
    CALL SOURCE_FUNCTION(K,L,SF)
    SFM(:,1)=SF(:)
    !
    ! If LTE
    IN(:,1)=SYN1D(:,L)
    IM=MATMUL(IMAT-EVOL,SFM)+MATMUL(EVOL,IN)
    SYN1D(:,L)=IM(:,1)
    !----------------------------------------
    ! If Non-LTE
    !IN(:,1)=STOKES(I,J,K-1,L,:)
    !IM=REAL(MATMUL(IMAT-EVOL,SFM)+MATMUL(EVOL,IN))
    !STOKES(I,J,K,L,:)=IM(:,1)
    !----------------------------------------
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      !
      ! Derivatives of IM:
      !
      ! Derivative of the Source function:
      CALL SOURCE_FUNCTION_DER(K,L,DSF)
      DSFM(:,:)=0.0D0
      DSFM(1,1)=DSF
      !
      ! dT at constant PG
      RTE_DT_PG=MATMUL(-DEVOLDT_PG,SFM)+MATMUL(IMAT-EVOL,DSFM)&
          +MATMUL(DEVOLDT_PG,IN)
      PCNT=0
      IF (INV_ATMPAR(1).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DT_PG(:,1)
      ENDIF
      ! dT at constant RHO (auxiliar)
      RTE_DT_RHO=MATMUL(-DEVOLDT_RHO,SFM)+MATMUL(IMAT-EVOL,DSFM)&
          +MATMUL(DEVOLDT_RHO,IN)
      !
      ! dP at constant T
      RTE_DPG_T=MATMUL(-DEVOLDPG_TEMP,SFM)+MATMUL(DEVOLDPG_TEMP,IN)
      IF (INV_ATMPAR(2).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DPG_T(:,1)
      ENDIF
      !
      ! dRHO at constant T
      RTE_DRHO_T=MATMUL(-DEVOLDRHO_TEMP,SFM)+MATMUL(DEVOLDRHO_TEMP,IN)
      IF (INV_ATMPAR(3).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DRHO_T(:,1)
      ENDIF
      ! dBx
      RTE_DBX=MATMUL(-DEVOLDBX,SFM)+MATMUL(DEVOLDBX,IN)
      IF (INV_ATMPAR(4).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DBX(:,1)
      ENDIF
      !
      ! dBy
      RTE_DBY=MATMUL(-DEVOLDBY,SFM)+MATMUL(DEVOLDBY,IN)
      IF (INV_ATMPAR(5).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DBY(:,1)
      ENDIF
      !
      ! dBz
      RTE_DBZ=MATMUL(-DEVOLDBZ,SFM)+MATMUL(DEVOLDBZ,IN)
      IF (INV_ATMPAR(6).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DBZ(:,1)
      ENDIF
      !
      ! dVLOS
      IM=MATMUL(-DEVOLDVLOS,SFM)+MATMUL(DEVOLDVLOS,IN)
      IF (INV_ATMPAR(7).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=IM(:,1)
      ENDIF
      !
      ! dP0
      IF (INV_ATMPAR(8).EQV..TRUE.) THEN
        PCNT=PCNT+1
        DSYN1D(:,PCNT,L,K)=RTE_DPG_T(:,1)
      ENDIF
      !
      !
      EVOLG(:,:,L,K)=EVOL(:,:)
      ! The response function has to be calculated with respect to the
      ! outgoing Stokes spectra, not outgoing each layer.
    ENDIF
    !
  END SUBROUTINE ANALYTICAL_SOLVER
  !
  !------------------------------------------------
  !
  SUBROUTINE OLD_ANALYTICAL_SOLVER(K,L)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)       :: K,L
    REAL(DP), DIMENSION(4)       :: SF
    REAL(DP), DIMENSION(4,1)     :: SFM, IM, STKIN
    REAL(DP)                     :: DSF
    ! For derivatives
    REAL(DP), DIMENSION(4,1)     :: RTE_DT_PG, RTE_DT_RHO, RTE_DPG_T, RTE_DRHO_T
    REAL(DP), DIMENSION(4,1)     :: RTE_DBX, RTE_DBY, RTE_DBZ
    REAL(DP), DIMENSION(4,1)     :: DSFM
    !
    REAL(DP), DIMENSION(4,4,NZ)     :: EVOLZ
    INTEGER :: KB
!
    CALL EVOL_OPERATOR(K,L)
    !
    CALL SOURCE_FUNCTION(K,L,SF)
    SFM(:,1)=SF(:)
    !
    ! If LTE
    STKIN(:,1)=SYN1D(:,L)
    IM=MATMUL(IMAT-EVOL,SFM)+MATMUL(EVOL,STKIN)
    SYN1D(:,L)=IM(:,1)
    !----------------------------------------
    ! If Non-LTE
    !STKIN(:,1)=STOKES(I,J,K-1,L,:)
    !IM=REAL(MATMUL(IMAT-EVOL,SFM)+MATMUL(EVOL,STKIN))
    !STOKES(I,J,K,L,:)=IM(:,1)
    !----------------------------------------
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      !
      EVOLG(:,:,L,K)=EVOL(:,:)
      ! Derivatives of IM:
      !
      ! Derivative of the Source function:
      CALL SOURCE_FUNCTION_DER(K,L,DSF)
      DSFM(:,:)=0.0D0
      DSFM(1,1)=DSF
      !

      ! First, dp_t
      ! It is so because if H. Eq., it is needed for dt_p
      !
      ! dP at constant T
      RTE_DPG_T=MATMUL(-DEVOLDPG_TEMP,SFM)+MATMUL(DEVOLDPG_TEMP,STKIN)
      DSYN1D(:,2,L,K)=RTE_DPG_T(:,1)

      !
      ! dT at constant PG
      ! 1) Local term:
      RTE_DT_PG=MATMUL(-DEVOLDT_PG,SFM)+MATMUL(IMAT-EVOL,DSFM)&
          +MATMUL(DEVOLDT_PG,STKIN)
      ! 2) Global term through HE?:
      IF (HYDROSTATICDER.EQV..TRUE.) THEN
        IF (HYDROSTATIC.EQV..TRUE.) THEN
            ! We never enter here if we are in k=1 so we do not need to...
            ! ...avoid K<=1:
            STKIN(:,:)=0
            EVOLZ(:,:,:)=EVOLG(:,:,L,:)
            DO KB=1,K,1
              IF (EVOLZ(1,1,KB).GT.1.0D-5) THEN
                STKIN(:,1)=DSYN1D(:,2,L,KB)*DPGDTK_HE(KB)+MATMUL(EVOLZ(:,:,KB), STKIN(:,1))
              ENDIF
            ENDDO
            RTE_DT_PG(:,1)=RTE_DT_PG(:,1)+STKIN(:,1)
            STKIN(:,1)=SYN1D(:,L)
        ENDIF
      ENDIF ! HYDROSTATICDER
      DSYN1D(:,1,L,K)=RTE_DT_PG(:,1)
      ! dT at constant RHO
      RTE_DT_RHO=MATMUL(-DEVOLDT_RHO,SFM)+MATMUL(IMAT-EVOL,DSFM)&
          +MATMUL(DEVOLDT_RHO,STKIN)
      !
      ! dRHO at constant T
      RTE_DRHO_T=MATMUL(-DEVOLDRHO_TEMP,SFM)+MATMUL(DEVOLDRHO_TEMP,STKIN)
      DSYN1D(:,3,L,K)=RTE_DRHO_T(:,1)
      ! dBx
      RTE_DBX=MATMUL(-DEVOLDBX,SFM)+MATMUL(DEVOLDBX,STKIN)
      DSYN1D(:,4,L,K)=RTE_DBX(:,1)
      !
      ! dBy
      RTE_DBY=MATMUL(-DEVOLDBY,SFM)+MATMUL(DEVOLDBY,STKIN)
      DSYN1D(:,5,L,K)=RTE_DBY(:,1)
      !
      ! dBz
      RTE_DBZ=MATMUL(-DEVOLDBZ,SFM)+MATMUL(DEVOLDBZ,STKIN)
      DSYN1D(:,6,L,K)=RTE_DBZ(:,1)
      !
      ! dVLOS
      IM=MATMUL(-DEVOLDVLOS,SFM)+MATMUL(DEVOLDVLOS,STKIN)
      DSYN1D(:,7,L,K)=IM(:,1)
      !
      ! dP0
      DSYN1D(:,8,L,K)=DSYN1D(:,8,L,K)*0.E0+1.E0
      !
      !
      ! The response function has to be calculated with respect to the
      ! outgoing Stokes spectra, not outgoing each layer.
    ENDIF
    !
  END SUBROUTINE OLD_ANALYTICAL_SOLVER
  !
  !------------------------------------------------
  !
  SUBROUTINE EVOL_OPERATOR(K,L)
    !
    IMPLICIT NONE
    !
    INTEGER,   INTENT(IN)           :: K,L
    !
    REAL(DP)                        :: ETA_I, ETA_Q, ETA_U, ETA_V
    REAL(DP)                        :: RHO_Q, RHO_U, RHO_V
    REAL(DP)                        :: ETA2, RHO2, ETARHO, SIGMA, THETA, RPART, EXPETAI
    REAL(DP)                        :: LAMBDA1, LAMBDA2, COSHL1, SINHL1, COSL2, SINL2
    REAL(DP), DIMENSION(4,4)        :: M1, M2, M3, M4
    LOGICAL                         :: NAN, CHECKNAN
    !
    ! Dummy variables for derivative
    REAL(DP)                        :: DIFFETA2RHO2, DAUX1, DAUX2, DAUX3, DELTAZ
    REAL(DP)                        :: DTHETADETAQ, DTHETADETAU, DTHETADETAV, DTHETADRHOQ, DTHETADRHOU, DTHETADRHOV
    REAL(DP)                        :: DLAMBDA1DETAQ, DLAMBDA1DETAU, DLAMBDA1DETAV&
        , DLAMBDA1DRHOQ, DLAMBDA1DRHOU, DLAMBDA1DRHOV
    REAL(DP)                        :: DLAMBDA2DETAQ, DLAMBDA2DETAU, DLAMBDA2DETAV&
        , DLAMBDA2DRHOQ, DLAMBDA2DRHOU, DLAMBDA2DRHOV
    ! M1 derivatives are 0
    REAL(DP), DIMENSION(4,4)        :: DM2DETAQ, DM2DETAU, DM2DETAV, DM2DRHOQ, DM2DRHOU, DM2DRHOV
    REAL(DP), DIMENSION(4,4)        :: DM3DETAQ, DM3DETAU, DM3DETAV, DM3DRHOQ, DM3DRHOU, DM3DRHOV
    REAL(DP), DIMENSION(4,4)        :: DM4DETAQ, DM4DETAU, DM4DETAV, DM4DRHOQ, DM4DRHOU, DM4DRHOV
    !
    REAL(DP), DIMENSION(4,4)        :: DEVOLDETAI, DEVOLDETAQ, DEVOLDETAU, DEVOLDETAV&
        , DEVOLDRHOQ, DEVOLDRHOU, DEVOLDRHOV
    REAL(DP), DIMENSION(4,4)        :: TERM1, TERM2, TERM3, TERM4, TERM5, TERM6, TERM7, TERM8
    !
    !
    ! NEW: September 21, 2017
    !
    LOGICAL                         :: EVOL_IMAT
    !
    ! END NEW: September 21, 2017
    !
    !
    M1(:,:)=0.0D0
    M2(:,:)=0.0D0
    M3(:,:)=0.0D0
    M4(:,:)=0.0D0
    SIGMA=0.
    !
    ETA_I=DBLE(ETAI(L))
    ETA_Q=DBLE(ETAQ(L))
    ETA_U=DBLE(ETAU(L))
    ETA_V=DBLE(ETAV(L))
    RHO_Q=DBLE(RHOQ(L))
    RHO_U=DBLE(RHOU(L))
    RHO_V=DBLE(RHOV(L))

    EVOL_IMAT = .FALSE.
    IF (ETA_I.EQ.(ETA_I+1.D0)) EVOL_IMAT = .TRUE.
    IF (ETA_Q.EQ.(ETA_Q+1.D0)) EVOL_IMAT = .TRUE.
    IF (ETA_U.EQ.(ETA_U+1.D0)) EVOL_IMAT = .TRUE.
    IF (ETA_V.EQ.(ETA_V+1.D0)) EVOL_IMAT = .TRUE.
    IF (RHO_Q.EQ.(RHO_Q+1.D0)) EVOL_IMAT = .TRUE.
    IF (RHO_U.EQ.(RHO_U+1.D0)) EVOL_IMAT = .TRUE.
    IF (RHO_V.EQ.(RHO_V+1.D0)) EVOL_IMAT = .TRUE.
    IF (EVOL_IMAT.EQV..TRUE.) THEN
      ETA_I=1.D15
      ETA_Q=1.D15
      ETA_U=1.D15
      ETA_V=1.D15
      RHO_Q=1.D15
      RHO_U=1.D15
      RHO_V=1.D15
    ENDIF
    ! Checking for NaN
    NAN=.FALSE.
    NAN=CHECKNAN(DBLE(ETA_I))
    NAN=CHECKNAN(DBLE(ETA_Q))
    NAN=CHECKNAN(DBLE(ETA_U))
    NAN=CHECKNAN(DBLE(ETA_V))
    NAN=CHECKNAN(DBLE(RHO_Q))
    NAN=CHECKNAN(DBLE(RHO_U))
    NAN=CHECKNAN(DBLE(RHO_V))
    IF (NAN.EQV..TRUE.) THEN
       PRINT*,'NaN detected in EVOL_OPERATOR'
       PRINT*,ETA_I,ETA_Q,ETA_U,ETA_V,RHO_Q,RHO_U,RHO_V
       PRINT*,TEM(K),RHO(K),PG(K)
       PRINT*,K,L
       STOP
    ENDIF
    !
    ! This part follows Landi Degl'Innocenti & Landi Degl'Innocenti
    ! Solar Physics, 1985, 239-250
    !
    ETA2=ETA_Q**2.+ETA_U**2.+ETA_V**2.
    RHO2=RHO_Q**2.+RHO_U**2.+RHO_V**2.
    ETARHO=ETA_Q*RHO_Q+ETA_U*RHO_U+ETA_V*RHO_V
    RPART=SQRT(((ETA2-RHO2)**2.)/4.+ETARHO**2.)+(ETA2-RHO2)/2.
    LAMBDA1=SQRT(RPART)
    !
    RPART=SQRT(((ETA2-RHO2)**2.)/4.+ETARHO**2.)-(ETA2-RHO2)/2.
    LAMBDA2=SQRT(RPART)
    IF (ETARHO.GT.0) SIGMA=1.
    IF (ETARHO.LT.0) SIGMA=-1.
    THETA=2.*SQRT(((ETA2-RHO2)**2.)/4.+ETARHO**2.)
    !
    ! M1 matrix is the indentiy matrix
    M1=IMAT
    ! M2 matrix
    M2(2,1)=LAMBDA2*ETA_Q-SIGMA*LAMBDA1*RHO_Q
    M2(3,1)=LAMBDA2*ETA_U-SIGMA*LAMBDA1*RHO_U
    M2(4,1)=LAMBDA2*ETA_V-SIGMA*LAMBDA1*RHO_V
    !M2(1,2)=REAL(LAMBDA2*ETA_Q-SIGMA*LAMBDA1*RHO_Q)
    M2(1,2)=M2(2,1)
    M2(3,2)=SIGMA*LAMBDA1*ETA_V+LAMBDA2*RHO_V
    M2(4,2)=-SIGMA*LAMBDA1*ETA_U-LAMBDA2*RHO_U
    !M2(1,3)=REAL(LAMBDA2*ETA_U-SIGMA*LAMBDA1*RHO_U)
    !M2(2,3)=REAL(-SIGMA*LAMBDA1*ETA_V-LAMBDA2*RHO_V)
    M2(1,3)=M2(3,1)
    M2(2,3)=-M2(3,2)
    M2(4,3)=SIGMA*LAMBDA1*ETA_Q+LAMBDA2*RHO_Q
    !M2(1,4)=REAL(LAMBDA2*ETA_V-SIGMA*LAMBDA1*RHO_V)
    !M2(2,4)=REAL(SIGMA*LAMBDA1*ETA_U+LAMBDA2*RHO_U)
    !M2(3,4)=REAL(-SIGMA*LAMBDA1*ETA_Q-LAMBDA2*RHO_Q)
    M2(1,4)=M2(4,1)
    M2(2,4)=-M2(4,2)
    M2(3,4)=-M2(4,3)
    M2(:,:)=M2(:,:)/THETA
    ! M3 matrix
    M3(2,1)=LAMBDA1*ETA_Q+SIGMA*LAMBDA2*RHO_Q
    M3(3,1)=LAMBDA1*ETA_U+SIGMA*LAMBDA2*RHO_U
    M3(4,1)=LAMBDA1*ETA_V+SIGMA*LAMBDA2*RHO_V
    !M3(1,2)=REAL(LAMBDA1*ETA_Q+SIGMA*LAMBDA2*RHO_Q)
    M3(1,2)=M3(2,1)
    M3(3,2)=-SIGMA*LAMBDA2*ETA_V+LAMBDA1*RHO_V
    M3(4,2)=SIGMA*LAMBDA2*ETA_U-LAMBDA1*RHO_U
    !M3(1,3)=REAL(LAMBDA1*ETA_U+SIGMA*LAMBDA2*RHO_U)
    !M3(2,3)=REAL(SIGMA*LAMBDA2*ETA_V-LAMBDA1*RHO_V)
    M3(1,3)=M3(3,1)
    M3(2,3)=-M3(3,2)
    M3(4,3)=-SIGMA*LAMBDA2*ETA_Q+LAMBDA1*RHO_Q
    !M3(1,4)=REAL(LAMBDA1*ETA_V+SIGMA*LAMBDA2*RHO_V)
    !M3(2,4)=REAL(-SIGMA*LAMBDA2*ETA_U+LAMBDA1*RHO_U)
    !M3(3,4)=REAL(SIGMA*LAMBDA2*ETA_Q-LAMBDA1*RHO_Q)
    M3(1,4)=M3(4,1)
    M3(2,4)=-M3(4,2)
    M3(3,4)=-M3(4,3)
    M3(:,:)=M3(:,:)/THETA
    ! M4 matrix
    M4(1,1)=(ETA2+RHO2)/2.
    M4(2,1)=ETA_V*RHO_U-ETA_U*RHO_V
    M4(3,1)=ETA_Q*RHO_V-ETA_V*RHO_Q
    M4(4,1)=ETA_U*RHO_Q-ETA_Q*RHO_U
    !M4(1,2)=REAL(ETA_U*RHO_V-ETA_V*RHO_U)
    M4(1,2)=-M4(2,1)
    !M4(2,2)=REAL(ETA_Q**2.+RHO_Q**2.-(ETA2+RHO2)/2.)
    M4(2,2)=ETA_Q**2.+RHO_Q**2.-M4(1,1)
    M4(3,2)=ETA_Q*ETA_U+RHO_Q*RHO_U
    M4(4,2)=ETA_V*ETA_Q+RHO_V*RHO_Q
    !M4(1,3)=REAL(ETA_V*RHO_Q-ETA_Q*RHO_V)
    !M4(2,3)=REAL(ETA_Q*ETA_U+RHO_Q*RHO_U)
    M4(1,3)=-M4(3,1)
    M4(2,3)=M4(3,2)
    !M4(3,3)=REAL(ETA_U**2.+RHO_U**2.-(ETA2+RHO2)/2.)
    M4(3,3)=ETA_U**2.+RHO_U**2.-M4(1,1)
    M4(4,3)=ETA_U*ETA_V+RHO_U*RHO_V
    !M4(1,4)=REAL(ETA_Q*RHO_U-ETA_U*RHO_Q)
    !M4(2,4)=REAL(ETA_V*ETA_Q+RHO_V*RHO_Q)
    !M4(3,4)=REAL(ETA_U*ETA_V+RHO_U*RHO_V)
    M4(1,4)=-M4(4,1)
    M4(2,4)=M4(4,2)
    M4(3,4)=M4(4,3)
    !M4(4,4)=REAL(ETA_V**2.+RHO_V**2.-(ETA2+RHO2)/2.)
    M4(4,4)=ETA_V**2.+RHO_V**2.-M4(1,1)
    M4(:,:)=2.*M4(:,:)/THETA
    !
    IF (K.EQ.NZ) THEN
      DELTAZ=1D5*ABS(ZZ(K)-ZZ(K-1))
    ELSE
      DELTAZ=1D5*ABS(ZZ(K+1)-ZZ(K))
    ENDIF
    !
    ! RF P0
!.    DSTOKES(K,L,:,8)=1.E0!DELTAZ
    !DSYN(K,L,:,8)=1.E0!DELTAZ
    ! END RF P0
    !
    IF (THETA.EQ.0) THEN
      M2(:,:)=0.0D0
      M3(:,:)=0.0D0
      M4(:,:)=0.0D0
    ENDIF
    !
    ! Multiplicative factors
    COSHL1=COSH(LAMBDA1*DELTAZ)
    SINHL1=SINH(LAMBDA1*DELTAZ)
    ! FOR GIVEN COMBINATIONS OF T, PG, AND RHO SINHL1 AND COSHL1 MIGHT BECOME
    ! EXTREMELY HUGE, INFINITY FOR PRACTICAL PURPOSES. EVENTUALLY, THIS GIVES RISE 
    ! TO NaN. TO AVOID SUCH SITUATION, WE LIMIT THE LARGEST VALUES OF BOTH COSHL1
    ! AND SINHL1. THIS HAS NO UNEXPECTED INFLUENCE ON THE CODE PERFORMANCE BECAUSE
    ! FOR LARGE ENOUGH COSHL1 AND SINHL1, EVOL BECOMES 0
    IF (COSHL1.EQ.(COSHL1+1.D0)) COSHL1=1.D15
    IF (SINHL1.EQ.(SINHL1+1.D0)) SINHL1=1.D15
    !WRITE(*,*) COSHL1, SINHL1
    !
    COSL2=COS(LAMBDA2*DELTAZ)
    SINL2=SIN(LAMBDA2*DELTAZ)
    EXPETAI=EXP(-ETA_I*DELTAZ)
    !
    ! Evolution operator
    EVOL=EXPETAI*(0.5D0*(COSHL1+COSL2)*M1-SINL2*M2-SINHL1*M3+0.5D0*(COSHL1-COSL2)*M4)
    IF (EVOL_IMAT.EQV..TRUE.) EVOL=IMAT
    !
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      !    
      ! Evolution operator derivatives:
      !
      !
      ! Some useful variables: (see documentation)
      DIFFETA2RHO2=ETA2-RHO2
      DAUX1=SQRT(0.25D0*DIFFETA2RHO2**2+ETARHO**2)
      DAUX2=0.5D0/SQRT(DAUX1+0.5D0*DIFFETA2RHO2)
      IF (DAUX2.EQ.(DAUX2+1.D0)) THEN
         DAUX2=0.D0
      ENDIF
      DAUX3=0.5D0/SQRT(DAUX1-0.5D0*DIFFETA2RHO2)
      NAN=.FALSE.
      NAN=CHECKNAN(DBLE(DAUX3))
      IF (NAN.EQV..TRUE.) DAUX3=0.D0
      ! it is convenient to change it: (see documentation)
      IF (DAUX1 .EQ. 0.D0) THEN
        DAUX1=0.D0
        DAUX3=0.D0
      ELSE
        DAUX1=1.D0/DAUX1
      ENDIF
      ! Since ETARHO is not used anymore in the original implementation, we multiplied it by 2
      ETARHO=2.*ETARHO
      !
      ! Derivatives of THETA with respect ETAs and RHOs:
      !
      DTHETADETAQ=DAUX1*(DIFFETA2RHO2*ETA_Q+ETARHO*RHO_Q)
      DTHETADETAU=DAUX1*(DIFFETA2RHO2*ETA_U+ETARHO*RHO_U)
      DTHETADETAV=DAUX1*(DIFFETA2RHO2*ETA_V+ETARHO*RHO_V)
      DTHETADRHOQ=DAUX1*(-DIFFETA2RHO2*RHO_Q+ETARHO*ETA_Q)
      DTHETADRHOU=DAUX1*(-DIFFETA2RHO2*RHO_U+ETARHO*ETA_U)
      DTHETADRHOV=DAUX1*(-DIFFETA2RHO2*RHO_V+ETARHO*ETA_V)
      !
      ! Derivatives of LAMBDA1 with respect ETAs and RHOs:
      !
      DLAMBDA1DETAQ=DAUX2*(0.5D0*DTHETADETAQ+ETA_Q)
      DLAMBDA1DETAU=DAUX2*(0.5D0*DTHETADETAU+ETA_U)
      DLAMBDA1DETAV=DAUX2*(0.5D0*DTHETADETAV+ETA_V)
      DLAMBDA1DRHOQ=DAUX2*(0.5D0*DTHETADRHOQ-RHO_Q)
      DLAMBDA1DRHOU=DAUX2*(0.5D0*DTHETADRHOU-RHO_U)
      DLAMBDA1DRHOV=DAUX2*(0.5D0*DTHETADRHOV-RHO_V)
      !
      ! Derivatives of LAMBDA2 with respect ETAs and RHOs:
      !
      DLAMBDA2DETAQ=DAUX3*(0.5D0*DTHETADETAQ-ETA_Q)
      DLAMBDA2DETAU=DAUX3*(0.5D0*DTHETADETAU-ETA_U)
      DLAMBDA2DETAV=DAUX3*(0.5D0*DTHETADETAV-ETA_V)
      DLAMBDA2DRHOQ=DAUX3*(0.5D0*DTHETADRHOQ+RHO_Q)
      DLAMBDA2DRHOU=DAUX3*(0.5D0*DTHETADRHOU+RHO_U)
      DLAMBDA2DRHOV=DAUX3*(0.5D0*DTHETADRHOV+RHO_V)
      !
      ! Derivatives of M2 to M4 with respect ETAs and RHOs:A
      !
      ! Zeroed matrixes
      !
      DM2DETAQ(:,:)=0.0
      DM2DETAU(:,:)=0.0
      DM2DETAV(:,:)=0.0
      DM2DRHOQ(:,:)=0.0
      DM2DRHOU(:,:)=0.0
      DM2DRHOV(:,:)=0.0
      DM3DETAQ(:,:)=0.0
      DM3DETAU(:,:)=0.0
      DM3DETAV(:,:)=0.0
      DM3DRHOQ(:,:)=0.0
      DM3DRHOU(:,:)=0.0
      DM3DRHOV(:,:)=0.0
      DM4DETAQ(:,:)=0.0
      DM4DETAU(:,:)=0.0
      DM4DETAV(:,:)=0.0
      DM4DRHOQ(:,:)=0.0
      DM4DRHOU(:,:)=0.0
      DM4DRHOV(:,:)=0.0
      IF (THETA.NE.0.D0) THEN
        !
        ! DM2
        !
        ! DM2dETAQ
        DM2DETAQ(2,1)=DLAMBDA2DETAQ*ETA_Q+LAMBDA2-SIGMA*DLAMBDA1DETAQ*RHO_Q
        DM2DETAQ(3,1)=DLAMBDA2DETAQ*ETA_U-SIGMA*DLAMBDA1DETAQ*RHO_U
        DM2DETAQ(4,1)=DLAMBDA2DETAQ*ETA_V-SIGMA*DLAMBDA1DETAQ*RHO_V
        DM2DETAQ(1,2)=DM2DETAQ(2,1)
        DM2DETAQ(3,2)=SIGMA*DLAMBDA1DETAQ*ETA_V+DLAMBDA2DETAQ*RHO_V
        DM2DETAQ(4,2)=-SIGMA*DLAMBDA1DETAQ*ETA_U-DLAMBDA2DETAQ*RHO_U
        DM2DETAQ(1,3)=DM2DETAQ(3,1)
        DM2DETAQ(2,3)=-DM2DETAQ(3,2)
        DM2DETAQ(4,3)=SIGMA*(DLAMBDA1DETAQ*ETA_Q+LAMBDA1)+DLAMBDA2DETAQ*RHO_Q
        DM2DETAQ(1,4)=DM2DETAQ(4,1)
        DM2DETAQ(2,4)=-DM2DETAQ(4,2)
        DM2DETAQ(3,4)=-DM2DETAQ(4,3)
        DM2DETAQ(:,:)=(-DTHETADETAQ*M2+DM2DETAQ(:,:))/THETA
        !
        ! DM2dETAU
        DM2DETAU(2,1)=DLAMBDA2DETAU*ETA_Q-SIGMA*DLAMBDA1DETAU*RHO_Q
        DM2DETAU(3,1)=DLAMBDA2DETAU*ETA_U+LAMBDA2-SIGMA*DLAMBDA1DETAU*RHO_U
        DM2DETAU(4,1)=DLAMBDA2DETAU*ETA_V-SIGMA*DLAMBDA1DETAU*RHO_V
        DM2DETAU(1,2)=DM2DETAU(2,1)
        DM2DETAU(3,2)=SIGMA*DLAMBDA1DETAU*ETA_V+DLAMBDA2DETAU*RHO_V
        DM2DETAU(4,2)=-SIGMA*(DLAMBDA1DETAU*ETA_U+LAMBDA1)-DLAMBDA2DETAU*RHO_U
        DM2DETAU(1,3)=DM2DETAU(3,1)
        DM2DETAU(2,3)=-DM2DETAU(3,2)
        DM2DETAU(4,3)=SIGMA*DLAMBDA1DETAU*ETA_Q+DLAMBDA2DETAU*RHO_Q
        DM2DETAU(1,4)=DM2DETAU(4,1)
        DM2DETAU(2,4)=-DM2DETAU(4,2)
        DM2DETAU(3,4)=-DM2DETAU(4,3)
        DM2DETAU(:,:)=(-DTHETADETAU*M2+DM2DETAU(:,:))/THETA
        !
        ! DM2dETAV
        DM2DETAV(2,1)=DLAMBDA2DETAV*ETA_Q-SIGMA*DLAMBDA1DETAV*RHO_Q
        DM2DETAV(3,1)=DLAMBDA2DETAV*ETA_U-SIGMA*DLAMBDA1DETAV*RHO_U
        DM2DETAV(4,1)=DLAMBDA2DETAV*ETA_V+LAMBDA2-SIGMA*DLAMBDA1DETAV*RHO_V
        DM2DETAV(1,2)=DM2DETAV(2,1)
        DM2DETAV(3,2)=SIGMA*(DLAMBDA1DETAV*ETA_V+LAMBDA1)+DLAMBDA2DETAV*RHO_V
        DM2DETAV(4,2)=-(SIGMA*DLAMBDA1DETAV*ETA_U+DLAMBDA2DETAV*RHO_U)
        DM2DETAV(1,3)=DM2DETAV(3,1)
        DM2DETAV(2,3)=-DM2DETAV(3,2)
        DM2DETAV(4,3)=SIGMA*DLAMBDA1DETAV*ETA_Q+DLAMBDA2DETAV*RHO_Q
        DM2DETAV(1,4)=DM2DETAV(4,1)
        DM2DETAV(2,4)=-DM2DETAV(4,2)
        DM2DETAV(3,4)=-DM2DETAV(4,3)
        DM2DETAV(:,:)=(-DTHETADETAV*M2+DM2DETAV(:,:))/THETA
        !
        ! DM2dRHOQ
        DM2DRHOQ(2,1)=DLAMBDA2DRHOQ*ETA_Q-SIGMA*(DLAMBDA1DRHOQ*RHO_Q+LAMBDA1)
        DM2DRHOQ(3,1)=DLAMBDA2DRHOQ*ETA_U-SIGMA*DLAMBDA1DRHOQ*RHO_U
        DM2DRHOQ(4,1)=DLAMBDA2DRHOQ*ETA_V-SIGMA*DLAMBDA1DRHOQ*RHO_V
        DM2DRHOQ(1,2)=DM2DRHOQ(2,1)
        DM2DRHOQ(3,2)=SIGMA*DLAMBDA1DRHOQ*ETA_V+DLAMBDA2DRHOQ*RHO_V
        DM2DRHOQ(4,2)=-(SIGMA*DLAMBDA1DRHOQ*ETA_U+DLAMBDA2DRHOQ*RHO_U)
        DM2DRHOQ(1,3)=DM2DRHOQ(3,1)
        DM2DRHOQ(2,3)=-DM2DRHOQ(3,2)
        DM2DRHOQ(4,3)=SIGMA*DLAMBDA1DRHOQ*ETA_Q+DLAMBDA2DRHOQ*RHO_Q+LAMBDA2
        DM2DRHOQ(1,4)=DM2DRHOQ(4,1)
        DM2DRHOQ(2,4)=-DM2DRHOQ(4,2)
        DM2DRHOQ(3,4)=-DM2DRHOQ(4,3)
        DM2DRHOQ(:,:)=(-DTHETADRHOQ*M2+DM2DRHOQ(:,:))/THETA
        !
        ! DM2dRHOU
        DM2DRHOU(2,1)=DLAMBDA2DRHOU*ETA_Q-SIGMA*DLAMBDA1DRHOU*RHO_Q
        DM2DRHOU(3,1)=DLAMBDA2DRHOU*ETA_U-SIGMA*(DLAMBDA1DRHOU*RHO_U+LAMBDA1)
        DM2DRHOU(4,1)=DLAMBDA2DRHOU*ETA_V-SIGMA*DLAMBDA1DRHOU*RHO_V
        DM2DRHOU(1,2)=DM2DRHOU(2,1)
        DM2DRHOU(3,2)=SIGMA*DLAMBDA1DRHOU*ETA_V+DLAMBDA2DRHOU*RHO_V
        DM2DRHOU(4,2)=-(SIGMA*DLAMBDA1DRHOU*ETA_U+(DLAMBDA2DRHOU*RHO_U+LAMBDA2))
        DM2DRHOU(1,3)=DM2DRHOU(3,1)
        DM2DRHOU(2,3)=-DM2DRHOU(3,2)
        DM2DRHOU(4,3)=SIGMA*DLAMBDA1DRHOU*ETA_Q+DLAMBDA2DRHOU*RHO_Q
        DM2DRHOU(1,4)=DM2DRHOU(4,1)
        DM2DRHOU(2,4)=-DM2DRHOU(4,2)
        DM2DRHOU(3,4)=-DM2DRHOU(4,3)
        DM2DRHOU(:,:)=(-DTHETADRHOU*M2+DM2DRHOU(:,:))/THETA
        !
        ! DM2dRHOV
        DM2DRHOV(2,1)=DLAMBDA2DRHOV*ETA_Q-SIGMA*DLAMBDA1DRHOV*RHO_Q
        DM2DRHOV(3,1)=DLAMBDA2DRHOV*ETA_U-SIGMA*DLAMBDA1DRHOV*RHO_U
        DM2DRHOV(4,1)=DLAMBDA2DRHOV*ETA_V-SIGMA*(DLAMBDA1DRHOV*RHO_V+LAMBDA1)
        DM2DRHOV(1,2)=DM2DRHOV(2,1)
        DM2DRHOV(3,2)=SIGMA*DLAMBDA1DRHOV*ETA_V+(DLAMBDA2DRHOV*RHO_V+LAMBDA2)
        DM2DRHOV(4,2)=-(SIGMA*DLAMBDA1DRHOV*ETA_U+DLAMBDA2DRHOV*RHO_U)
        DM2DRHOV(1,3)=DM2DRHOV(3,1)
        DM2DRHOV(2,3)=-DM2DRHOV(3,2)
        DM2DRHOV(4,3)=SIGMA*DLAMBDA1DRHOV*ETA_Q+DLAMBDA2DRHOV*RHO_Q
        DM2DRHOV(1,4)=DM2DRHOV(4,1)
        DM2DRHOV(2,4)=-DM2DRHOV(4,2)
        DM2DRHOV(3,4)=-DM2DRHOV(4,3)
        DM2DRHOV(:,:)=(-DTHETADRHOV*M2+DM2DRHOV(:,:))/THETA
        !
        ! DM3
        !
        ! DM3dETAQ
        DM3DETAQ(2,1)=DLAMBDA1DETAQ*ETA_Q+LAMBDA1+SIGMA*DLAMBDA2DETAQ*RHO_Q
        DM3DETAQ(3,1)=DLAMBDA1DETAQ*ETA_U+SIGMA*DLAMBDA2DETAQ*RHO_U
        DM3DETAQ(4,1)=DLAMBDA1DETAQ*ETA_V+SIGMA*DLAMBDA2DETAQ*RHO_V
        DM3DETAQ(1,2)=DM3DETAQ(2,1)
        DM3DETAQ(3,2)=-SIGMA*DLAMBDA2DETAQ*ETA_V+DLAMBDA1DETAQ*RHO_V
        DM3DETAQ(4,2)=SIGMA*DLAMBDA2DETAQ*ETA_U-DLAMBDA1DETAQ*RHO_U
        DM3DETAQ(1,3)=DM3DETAQ(3,1)
        DM3DETAQ(2,3)=-DM3DETAQ(3,2)
        DM3DETAQ(4,3)=-SIGMA*(DLAMBDA2DETAQ*ETA_Q+LAMBDA2)+DLAMBDA1DETAQ*RHO_Q
        DM3DETAQ(1,4)=DM3DETAQ(4,1)
        DM3DETAQ(2,4)=-DM3DETAQ(4,2)
        DM3DETAQ(3,4)=-DM3DETAQ(4,3)
        DM3DETAQ(:,:)=(-DTHETADETAQ*M3+DM3DETAQ(:,:))/THETA
        !
        ! DM3dETAU
        DM3DETAU(2,1)=DLAMBDA1DETAU*ETA_Q+SIGMA*DLAMBDA2DETAU*RHO_Q
        DM3DETAU(3,1)=DLAMBDA1DETAU*ETA_U+LAMBDA1+SIGMA*DLAMBDA2DETAU*RHO_U
        DM3DETAU(4,1)=DLAMBDA1DETAU*ETA_V+SIGMA*DLAMBDA2DETAU*RHO_V
        DM3DETAU(1,2)=DM3DETAU(2,1)
        DM3DETAU(3,2)=-(SIGMA*DLAMBDA2DETAU*ETA_V-DLAMBDA1DETAU*RHO_V)
        DM3DETAU(4,2)=SIGMA*(DLAMBDA2DETAU*ETA_U+LAMBDA2)-DLAMBDA1DETAU*RHO_U
        DM3DETAU(1,3)=DM3DETAU(3,1)
        DM3DETAU(2,3)=-DM3DETAU(3,2)
        DM3DETAU(4,3)=-SIGMA*DLAMBDA2DETAU*ETA_Q+DLAMBDA1DETAU*RHO_Q
        DM3DETAU(1,4)=DM3DETAU(4,1)
        DM3DETAU(2,4)=-DM3DETAU(4,2)
        DM3DETAU(3,4)=-DM3DETAU(4,3)
        DM3DETAU(:,:)=(-DTHETADETAU*M3+DM3DETAU(:,:))/THETA
        !
        ! DM3dETAV
        DM3DETAV(2,1)=DLAMBDA1DETAV*ETA_Q+SIGMA*DLAMBDA2DETAV*RHO_Q
        DM3DETAV(3,1)=DLAMBDA1DETAV*ETA_U+SIGMA*DLAMBDA2DETAV*RHO_U
        DM3DETAV(4,1)=DLAMBDA1DETAV*ETA_V+LAMBDA1+SIGMA*DLAMBDA2DETAV*RHO_V
        DM3DETAV(1,2)=DM3DETAV(2,1)
        DM3DETAV(3,2)=-SIGMA*(DLAMBDA2DETAV*ETA_V+LAMBDA2)+DLAMBDA1DETAV*RHO_V
        DM3DETAV(4,2)=SIGMA*DLAMBDA2DETAV*ETA_U-DLAMBDA1DETAV*RHO_U
        DM3DETAV(1,3)=DM3DETAV(3,1)
        DM3DETAV(2,3)=-DM3DETAV(3,2)
        DM3DETAV(4,3)=-SIGMA*DLAMBDA2DETAV*ETA_Q+DLAMBDA1DETAV*RHO_Q
        DM3DETAV(1,4)=DM3DETAV(4,1)
        DM3DETAV(2,4)=-DM3DETAV(4,2)
        DM3DETAV(3,4)=-DM3DETAV(4,3)
        DM3DETAV(:,:)=(-DTHETADETAV*M3+DM3DETAV(:,:))/THETA
        !
        ! DM3dRHOQ
        DM3DRHOQ(2,1)=DLAMBDA1DRHOQ*ETA_Q+SIGMA*(DLAMBDA2DRHOQ*RHO_Q+LAMBDA2)
        DM3DRHOQ(3,1)=DLAMBDA1DRHOQ*ETA_U+SIGMA*DLAMBDA2DRHOQ*RHO_U
        DM3DRHOQ(4,1)=DLAMBDA1DRHOQ*ETA_V+SIGMA*DLAMBDA2DRHOQ*RHO_V
        DM3DRHOQ(1,2)=DM3DRHOQ(2,1)
        DM3DRHOQ(3,2)=-SIGMA*DLAMBDA2DRHOQ*ETA_V+DLAMBDA1DRHOQ*RHO_V
        DM3DRHOQ(4,2)=SIGMA*DLAMBDA2DRHOQ*ETA_U-DLAMBDA1DRHOQ*RHO_U
        DM3DRHOQ(1,3)=DM3DRHOQ(3,1)
        DM3DRHOQ(2,3)=-DM3DRHOQ(3,2)
        DM3DRHOQ(4,3)=-SIGMA*DLAMBDA2DRHOQ*ETA_Q+DLAMBDA1DRHOQ*RHO_Q+LAMBDA1
        DM3DRHOQ(1,4)=DM3DRHOQ(4,1)
        DM3DRHOQ(2,4)=-DM3DRHOQ(4,2)
        DM3DRHOQ(3,4)=-DM3DRHOQ(4,3)
        DM3DRHOQ(:,:)=(-DTHETADRHOQ*M3+DM3DRHOQ(:,:))/THETA
        !
        ! DM3dRHOU
        DM3DRHOU(2,1)=DLAMBDA1DRHOU*ETA_Q+SIGMA*DLAMBDA2DRHOU*RHO_Q
        DM3DRHOU(3,1)=DLAMBDA1DRHOU*ETA_U+SIGMA*(DLAMBDA2DRHOU*RHO_U+LAMBDA2)
        DM3DRHOU(4,1)=DLAMBDA1DRHOU*ETA_V+SIGMA*DLAMBDA2DRHOU*RHO_V
        DM3DRHOU(1,2)=DM3DRHOU(2,1)
        DM3DRHOU(3,2)=-SIGMA*DLAMBDA2DRHOU*ETA_V+DLAMBDA1DRHOU*RHO_V
        DM3DRHOU(4,2)=SIGMA*DLAMBDA2DRHOU*ETA_U-(DLAMBDA1DRHOU*RHO_U+LAMBDA1)
        DM3DRHOU(1,3)=DM3DRHOU(3,1)
        DM3DRHOU(2,3)=-DM3DRHOU(3,2)
        DM3DRHOU(4,3)=-SIGMA*DLAMBDA2DRHOU*ETA_Q+DLAMBDA1DRHOU*RHO_Q
        DM3DRHOU(1,4)=DM3DRHOU(4,1)
        DM3DRHOU(2,4)=-DM3DRHOU(4,2)
        DM3DRHOU(3,4)=-DM3DRHOU(4,3)
        DM3DRHOU(:,:)=(-DTHETADRHOU*M3+DM3DRHOU(:,:))/THETA
        !
        ! DM3dRHOV
        DM3DRHOV(2,1)=DLAMBDA1DRHOV*ETA_Q+SIGMA*DLAMBDA2DRHOV*RHO_Q
        DM3DRHOV(3,1)=DLAMBDA1DRHOV*ETA_U+SIGMA*DLAMBDA2DRHOV*RHO_U
        DM3DRHOV(4,1)=DLAMBDA1DRHOV*ETA_V+SIGMA*(DLAMBDA2DRHOV*RHO_V+LAMBDA2)
        DM3DRHOV(1,2)=DM3DRHOV(2,1)
        DM3DRHOV(3,2)=-SIGMA*DLAMBDA2DRHOV*ETA_V+DLAMBDA1DRHOV*RHO_V+LAMBDA1
        DM3DRHOV(4,2)=SIGMA*DLAMBDA2DRHOV*ETA_U-DLAMBDA1DRHOV*RHO_U
        DM3DRHOV(1,3)=DM3DRHOV(3,1)
        DM3DRHOV(2,3)=-DM3DRHOV(3,2)
        DM3DRHOV(4,3)=-SIGMA*DLAMBDA2DRHOV*ETA_Q+DLAMBDA1DRHOV*RHO_Q
        DM3DRHOV(1,4)=DM3DRHOV(4,1)
        DM3DRHOV(2,4)=-DM3DRHOV(4,2)
        DM3DRHOV(3,4)=-DM3DRHOV(4,3)
        DM3DRHOV(:,:)=(-DTHETADRHOV*M3+DM3DRHOV(:,:))/THETA
        !
        ! DM4
        !
        ! DM4dETAQ
        DM4DETAQ(1,1)=ETA_Q
        DM4DETAQ(3,1)=RHO_V
        DM4DETAQ(4,1)=-RHO_U
        DM4DETAQ(2,2)=ETA_Q
        DM4DETAQ(3,2)=ETA_U
        DM4DETAQ(4,2)=ETA_V
        DM4DETAQ(1,3)=-RHO_V
        DM4DETAQ(2,3)=ETA_U
        DM4DETAQ(3,3)=-ETA_Q
        DM4DETAQ(1,4)=RHO_U
        DM4DETAQ(2,4)=ETA_V
        DM4DETAQ(4,4)=-ETA_Q
        DM4DETAQ(:,:)=(-DTHETADETAQ*M4+2.*DM4DETAQ(:,:))/THETA
        !
        ! DM4dETAU
        DM4DETAU(1,1)=ETA_U
        DM4DETAU(2,1)=-RHO_V
        DM4DETAU(4,1)=RHO_Q
        DM4DETAU(1,2)=RHO_V
        DM4DETAU(2,2)=-ETA_U
        DM4DETAU(3,2)=ETA_Q
        DM4DETAU(2,3)=ETA_Q
        DM4DETAU(3,3)=ETA_U
        DM4DETAU(4,3)=ETA_V
        DM4DETAU(1,4)=-RHO_Q
        DM4DETAU(3,4)=ETA_V
        DM4DETAU(4,4)=-ETA_U
        DM4DETAU(:,:)=(-DTHETADETAU*M4+2.*DM4DETAU(:,:))/THETA
        !
        ! DM4dETAV
        DM4DETAV(1,1)=ETA_V
        DM4DETAV(2,1)=RHO_U
        DM4DETAV(3,1)=-RHO_Q
        DM4DETAV(1,2)=-RHO_U
        DM4DETAV(2,2)=-ETA_V
        DM4DETAV(4,2)=ETA_Q
        DM4DETAV(1,3)=RHO_Q
        DM4DETAV(3,3)=-ETA_V
        DM4DETAV(4,3)=ETA_U
        DM4DETAV(2,4)=ETA_Q
        DM4DETAV(3,4)=ETA_U
        DM4DETAV(4,4)=ETA_V
        DM4DETAV(:,:)=(-DTHETADETAV*M4+2.*DM4DETAV(:,:))/THETA
        !
        ! DM4dRHOQ
        DM4DRHOQ(1,1)=RHO_Q
        DM4DRHOQ(3,1)=-ETA_V
        DM4DRHOQ(4,1)=ETA_U
        DM4DRHOQ(2,2)=RHO_Q
        DM4DRHOQ(3,2)=RHO_U
        DM4DRHOQ(4,2)=RHO_V
        DM4DRHOQ(1,3)=ETA_V
        DM4DRHOQ(2,3)=RHO_U
        DM4DRHOQ(3,3)=-RHO_Q
        DM4DRHOQ(1,4)=-ETA_U
        DM4DRHOQ(2,4)=RHO_V
        DM4DRHOQ(4,4)=-RHO_Q
        DM4DRHOQ(:,:)=(-DTHETADRHOQ*M4+2.*DM4DRHOQ(:,:))/THETA
        !
        ! DM4dRHOU
        DM4DRHOU(1,1)=RHO_U
        DM4DRHOU(2,1)=ETA_V
        DM4DRHOU(4,1)=-ETA_Q
        DM4DRHOU(1,2)=-ETA_V
        DM4DRHOU(2,2)=-RHO_U
        DM4DRHOU(3,2)=RHO_Q
        DM4DRHOU(2,3)=RHO_Q
        DM4DRHOU(3,3)=RHO_U
        DM4DRHOU(4,3)=RHO_V
        DM4DRHOU(1,4)=ETA_Q
        DM4DRHOU(3,4)=RHO_V
        DM4DRHOU(4,4)=-RHO_U
        DM4DRHOU(:,:)=(-DTHETADRHOU*M4+2.*DM4DRHOU(:,:))/THETA
        !
        ! DM4dRHOV
        DM4DRHOV(1,1)=RHO_V
        DM4DRHOV(2,1)=-ETA_U
        DM4DRHOV(3,1)=ETA_Q
        DM4DRHOV(1,2)=ETA_U
        DM4DRHOV(2,2)=-RHO_V
        DM4DRHOV(4,2)=RHO_Q
        DM4DRHOV(1,3)=-ETA_Q
        DM4DRHOV(3,3)=-RHO_V
        DM4DRHOV(4,3)=RHO_U
        DM4DRHOV(2,4)=RHO_Q
        DM4DRHOV(3,4)=RHO_U
        DM4DRHOV(4,4)=RHO_V
        DM4DRHOV(:,:)=(-DTHETADRHOV*M4+2.*DM4DRHOV(:,:))/THETA
        !
      ENDIF
      !
      ! Derivatives of Stokes-attenuation operator with respect ETAs and RHOs:
      !
      ! dEVOLdETAI
      DEVOLDETAI=-DELTAZ*EVOL
      !
      ! dEVOLdETAQ
      TERM1(:,:)=0.5D0*(SINHL1*DLAMBDA1DETAQ-SINL2*DLAMBDA2DETAQ)*DELTAZ*M1(:,:)
      TERM2(:,:)=0.0D0
      TERM3(:,:)=-COSL2*DELTAZ*DLAMBDA2DETAQ*M2(:,:)
      TERM4(:,:)=-SINL2*DM2DETAQ(:,:)
      TERM5(:,:)=-COSHL1*DELTAZ*DLAMBDA1DETAQ*M3(:,:)
      TERM6(:,:)=-SINHL1*DM3DETAQ(:,:)
      TERM7(:,:)=0.5D0*(SINHL1*DLAMBDA1DETAQ+SINL2*DLAMBDA2DETAQ)*DELTAZ*M4(:,:)
      TERM8(:,:)=0.5D0*(COSHL1-COSL2)*DM4DETAQ(:,:)
      DEVOLDETAQ(:,:)=EXPETAI*(TERM1(:,:)+TERM2(:,:)+TERM3(:,:)+TERM4(:,:)+TERM5(:,:)+TERM6(:,:)+TERM7(:,:)+TERM8(:,:))
      !
      ! dEVOLdETAU
      TERM1(:,:)=0.5D0*(SINHL1*DLAMBDA1DETAU-SINL2*DLAMBDA2DETAU)*DELTAZ*M1(:,:)
      TERM2(:,:)=0.0D0
      TERM3(:,:)=-COSL2*DELTAZ*DLAMBDA2DETAU*M2(:,:)
      TERM4(:,:)=-SINL2*DM2DETAU(:,:)
      TERM5(:,:)=-COSHL1*DELTAZ*DLAMBDA1DETAU*M3(:,:)
      TERM6(:,:)=-SINHL1*DM3DETAU(:,:)
      TERM7(:,:)=0.5D0*(SINHL1*DLAMBDA1DETAU+SINL2*DLAMBDA2DETAU)*DELTAZ*M4(:,:)
      TERM8(:,:)=0.5D0*(COSHL1-COSL2)*DM4DETAU(:,:)
      DEVOLDETAU(:,:)=EXPETAI*(TERM1(:,:)+TERM2(:,:)+TERM3(:,:)+TERM4(:,:)+TERM5(:,:)+TERM6(:,:)+TERM7(:,:)+TERM8(:,:))
      !
      ! dEVOLdETAV
      TERM1(:,:)=0.5D0*(SINHL1*DLAMBDA1DETAV-SINL2*DLAMBDA2DETAV)*DELTAZ*M1(:,:)
      TERM2(:,:)=0.0D0
      TERM3(:,:)=-COSL2*DELTAZ*DLAMBDA2DETAV*M2(:,:)
      TERM4(:,:)=-SINL2*DM2DETAV(:,:)
      TERM5(:,:)=-COSHL1*DELTAZ*DLAMBDA1DETAV*M3(:,:)
      TERM6(:,:)=-SINHL1*DM3DETAV(:,:)
      TERM7(:,:)=0.5D0*(SINHL1*DLAMBDA1DETAV+SINL2*DLAMBDA2DETAV)*DELTAZ*M4(:,:)
      TERM8(:,:)=0.5D0*(COSHL1-COSL2)*DM4DETAV(:,:)
      DEVOLDETAV(:,:)=EXPETAI*(TERM1(:,:)+TERM2(:,:)+TERM3(:,:)+TERM4(:,:)+TERM5(:,:)+TERM6(:,:)+TERM7(:,:)+TERM8(:,:))
      !
      ! dEVOLdRHOQ
      TERM1(:,:)=0.5D0*(SINHL1*DLAMBDA1DRHOQ-SINL2*DLAMBDA2DRHOQ)*DELTAZ*M1(:,:)
      TERM2(:,:)=0.0D0
      TERM3(:,:)=-COSL2*DELTAZ*DLAMBDA2DRHOQ*M2(:,:)
      TERM4(:,:)=-SINL2*DM2DRHOQ(:,:)
      TERM5(:,:)=-COSHL1*DELTAZ*DLAMBDA1DRHOQ*M3(:,:)
      TERM6(:,:)=-SINHL1*DM3DRHOQ(:,:)
      TERM7(:,:)=0.5D0*(SINHL1*DLAMBDA1DRHOQ+SINL2*DLAMBDA2DRHOQ)*DELTAZ*M4(:,:)
      TERM8(:,:)=0.5D0*(COSHL1-COSL2)*DM4DRHOQ(:,:)
      DEVOLDRHOQ(:,:)=EXPETAI*(TERM1(:,:)+TERM2(:,:)+TERM3(:,:)+TERM4(:,:)+TERM5(:,:)+TERM6(:,:)+TERM7(:,:)+TERM8(:,:))
      !
      ! dEVOLdRHOU
      TERM1(:,:)=0.5D0*(SINHL1*DLAMBDA1DRHOU-SINL2*DLAMBDA2DRHOU)*DELTAZ*M1(:,:)
      TERM2(:,:)=0.0D0
      TERM3(:,:)=-COSL2*DELTAZ*DLAMBDA2DRHOU*M2(:,:)
      TERM4(:,:)=-SINL2*DM2DRHOU(:,:)
      TERM5(:,:)=-COSHL1*DELTAZ*DLAMBDA1DRHOU*M3(:,:)
      TERM6(:,:)=-SINHL1*DM3DRHOU(:,:)
      TERM7(:,:)=0.5D0*(SINHL1*DLAMBDA1DRHOU+SINL2*DLAMBDA2DRHOU)*DELTAZ*M4(:,:)
      TERM8(:,:)=0.5D0*(COSHL1-COSL2)*DM4DRHOU(:,:)
      DEVOLDRHOU(:,:)=EXPETAI*(TERM1(:,:)+TERM2(:,:)+TERM3(:,:)+TERM4(:,:)+TERM5(:,:)+TERM6(:,:)+TERM7(:,:)+TERM8(:,:))
      !
      ! dEVOLdRHOV
      TERM1(:,:)=0.5D0*(SINHL1*DLAMBDA1DRHOV-SINL2*DLAMBDA2DRHOV)*DELTAZ*M1(:,:)
      TERM2(:,:)=0.0D0
      TERM3(:,:)=-COSL2*DELTAZ*DLAMBDA2DRHOV*M2(:,:)
      TERM4(:,:)=-SINL2*DM2DRHOV(:,:)
      TERM5(:,:)=-COSHL1*DELTAZ*DLAMBDA1DRHOV*M3(:,:)
      TERM6(:,:)=-SINHL1*DM3DRHOV(:,:)
      TERM7(:,:)=0.5D0*(SINHL1*DLAMBDA1DRHOV+SINL2*DLAMBDA2DRHOV)*DELTAZ*M4(:,:)
      TERM8(:,:)=0.5D0*(COSHL1-COSL2)*DM4DRHOV(:,:)
      DEVOLDRHOV(:,:)=EXPETAI*(TERM1(:,:)+TERM2(:,:)+TERM3(:,:)+TERM4(:,:)+TERM5(:,:)+TERM6(:,:)+TERM7(:,:)+TERM8(:,:))
      !
      ! Now we have EVOL derivatives with respect to the various ETAs and RHOs
      ! Next step is to get the derivatives with respect to the physical parameters:
      ! T, b, gamma, phi, vlos
      !
      ! Zeroed matrixes:
      !
      DEVOLDBX(:,:)=0.0D0
      DEVOLDBY(:,:)=0.0D0
      DEVOLDBZ(:,:)=0.0D0
      DEVOLDVLOS(:,:)=0.0D0
      DEVOLDT_PG(:,:)=0.0D0
      DEVOLDT_RHO(:,:)=0.0D0
      DEVOLDPG_TEMP(:,:)=0.0D0
      DEVOLDRHO_TEMP(:,:)=0.0D0
      !
      IF (EVOL_IMAT.EQV..FALSE.) THEN
        !
        ! EVOL derivative with respect to Bx
        DEVOLDBX(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DBX(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DBX(L)&
            +DEVOLDETAU(:,:)*DETAU_DBX(L)&
            +DEVOLDETAV(:,:)*DETAV_DBX(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DBX(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DBX(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DBX(L)
        !
        ! EVOL derivative with respect to By
        DEVOLDBY(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DBY(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DBY(L)&
            +DEVOLDETAU(:,:)*DETAU_DBY(L)&
            +DEVOLDETAV(:,:)*DETAV_DBY(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DBY(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DBY(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DBY(L)
        !
        ! EVOL derivative with respect to Bz
        DEVOLDBZ(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DBZ(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DBZ(L)&
            +DEVOLDETAU(:,:)*DETAU_DBZ(L)&
            +DEVOLDETAV(:,:)*DETAV_DBZ(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DBZ(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DBZ(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DBZ(L)
        !
        ! EVOL derivative with respect to VLOS
        DEVOLDVLOS(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DVLOS(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DVLOS(L)&
            +DEVOLDETAU(:,:)*DETAU_DVLOS(L)&
            +DEVOLDETAV(:,:)*DETAV_DVLOS(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DVLOS(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DVLOS(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DVLOS(L)
        !
        ! EVOL derivative with respect to T at constant PG
        DEVOLDT_PG(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DTPG(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DTPG(L)&
            +DEVOLDETAU(:,:)*DETAU_DTPG(L)&
            +DEVOLDETAV(:,:)*DETAV_DTPG(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DTPG(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DTPG(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DTPG(L)
        !
        ! EVOL derivative with respect to T at constant RHO
        DEVOLDT_RHO(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DTRHO(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DTRHO(L)&
            +DEVOLDETAU(:,:)*DETAU_DTRHO(L)&
            +DEVOLDETAV(:,:)*DETAV_DTRHO(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DTRHO(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DTRHO(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DTRHO(L)
        !
        ! EVOL derivative with respect to PG at constant T
        DEVOLDPG_TEMP(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DPGT(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DPGT(L)&
            +DEVOLDETAU(:,:)*DETAU_DPGT(L)&
            +DEVOLDETAV(:,:)*DETAV_DPGT(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DPGT(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DPGT(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DPGT(L)
        !
        ! EVOL derivative with respect to rho at constant T
        DEVOLDRHO_TEMP(:,:)=&
            DEVOLDETAI(:,:)*DETAI_DRHOT(L)&
            +DEVOLDETAQ(:,:)*DETAQ_DRHOT(L)&
            +DEVOLDETAU(:,:)*DETAU_DRHOT(L)&
            +DEVOLDETAV(:,:)*DETAV_DRHOT(L)&
            +DEVOLDRHOQ(:,:)*DRHOQ_DRHOT(L)&
            +DEVOLDRHOU(:,:)*DRHOU_DRHOT(L)&
            +DEVOLDRHOV(:,:)*DRHOV_DRHOT(L)
        !
        ! END NEW: September 21, 2017
        !
      ENDIF
    ENDIF
    !
  END SUBROUTINE EVOL_OPERATOR
  !
  !------------------------------------------------
  !
   SUBROUTINE SOURCE_FUNCTION(K,L,SF)
    !
    USE FORWARD_PARAM, ONLY: NUML, PIXEL_END, PIXEL_INI &
        , WAVE, LINE_L0, IND_LINE, POPL, POPU, LINE_NUM
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)                  :: K,L
    REAL(DP), INTENT(OUT), DIMENSION(4)     :: SF
    INTEGER                                 :: M
    REAL(DP)                                :: WAVELENGTH=0.D0
    !
    LINESELECT: DO M=1,NUML
       IF (L.GE.PIXEL_INI(M).AND.L.LE.PIXEL_END(M)) THEN
          WAVELENGTH=LINE_L0(IND_LINE(M))+WAVE(L)/1.0D3
          EXIT LINESELECT
       ENDIF
    ENDDO LINESELECT
    ! cm !
    WAVELENGTH=WAVELENGTH*1D-8
    ! We assume black body radiation (LTE) + fixed NLTE departure coefficients  + 
    ! + no polarization for the source function
!print*, shape(popu)
    SF(:)=0.0D0
    SF(1)=2.0D0*HPLA*LIGHT**2.0D0/(WAVELENGTH**5.0D0)&
        *(1.0D0/((POPL(M,K)/POPU(M,K))*EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(K)))-1.0D0))
    !
  END SUBROUTINE SOURCE_FUNCTION
  !
  !------------------------------------------------
  !
  SUBROUTINE SOURCE_FUNCTION_DER(K,L,DSF)
    ! Calculates the derivative of the source function with respect to the temperature
    !
    USE FORWARD_PARAM, ONLY: NUML, PIXEL_END, PIXEL_INI &
        , WAVE, LINE_L0, IND_LINE, POPL, POPU, LINE_NUM
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)                  :: K,L
    REAL(DP), INTENT(OUT)                   :: DSF
    INTEGER                                 :: M
    REAL(DP)                                :: WAVELENGTH=0.D0
    REAL(DP), PARAMETER                     :: FACTOR = 8.5685641D-6 ! h^2*c^3/K
    !
    LINESELECT: DO M=1,NUML
       IF (L.GE.PIXEL_INI(M).AND.L.LE.PIXEL_END(M)) THEN
          WAVELENGTH=LINE_L0(IND_LINE(M))+WAVE(L)/1.0D3
          EXIT LINESELECT
       ENDIF
    ENDDO LINESELECT
    ! cm !
    WAVELENGTH=WAVELENGTH*1D-8
    !
    DSF = REAL((2.0D0*FACTOR/(WAVELENGTH**6.0D0*TEM(K)**2.0D0))&
        *((POPL(M,K)/POPU(M,K))*EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(K)))&
        /((POPL(M,K)/POPU(M,K))*EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(K)))-1.0D0)**2.0D0))
    !
  END SUBROUTINE SOURCE_FUNCTION_DER
  !
! Obsolete after including nlte departure coefficients:
! Obsolete after including nlte departure coefficients:  SUBROUTINE SOURCE_FUNCTION(K,L,SF)
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    USE FORWARD_PARAM, ONLY: NUML, PIXEL_END, PIXEL_INI &
! Obsolete after including nlte departure coefficients:        , WAVE, LINE_L0, IND_LINE
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    IMPLICIT NONE
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    INTEGER,    INTENT(IN)                  :: K,L
! Obsolete after including nlte departure coefficients:    REAL(DP), INTENT(OUT), DIMENSION(4)     :: SF
! Obsolete after including nlte departure coefficients:    INTEGER                                 :: M
! Obsolete after including nlte departure coefficients:    REAL(DP)                                :: WAVELENGTH=0.D0
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    DO M=1,NUML
! Obsolete after including nlte departure coefficients:       !IF (L.GE.PIXEL_INI(M).AND.M.LE.PIXEL_END(M)) WAVELENGTH&
! Obsolete after including nlte departure coefficients:       !    =LINE_L0(IND_LINE(M))+WAVE(L)/1.0D3
! Obsolete after including nlte departure coefficients:       IF (L.GE.PIXEL_INI(M).AND.L.LE.PIXEL_END(M)) THEN
! Obsolete after including nlte departure coefficients:         WAVELENGTH = LINE_L0(IND_LINE(M))+WAVE(L)/1.0D3
! Obsolete after including nlte departure coefficients:         EXIT
! Obsolete after including nlte departure coefficients:       ENDIF
! Obsolete after including nlte departure coefficients:    ENDDO
! Obsolete after including nlte departure coefficients:    ! cm !
! Obsolete after including nlte departure coefficients:    WAVELENGTH=WAVELENGTH*1D-8
! Obsolete after including nlte departure coefficients:    ! We assume black body radiation (LTE) and no polarization for the 
! Obsolete after including nlte departure coefficients:    ! source function
! Obsolete after including nlte departure coefficients:    SF(:)=0.0D0
! Obsolete after including nlte departure coefficients:    SF(1)=2.0D0*HPLA*LIGHT**2.0D0/(WAVELENGTH**5.0D0)&
! Obsolete after including nlte departure coefficients:        *(1.0D0/(EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(K)))-1.0D0))
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:  END SUBROUTINE SOURCE_FUNCTION
! Obsolete after including nlte departure coefficients:  !
! Obsolete after including nlte departure coefficients:  !------------------------------------------------
! Obsolete after including nlte departure coefficients:  !
! Obsolete after including nlte departure coefficients:  SUBROUTINE SOURCE_FUNCTION_DER(K,L,DSF)
! Obsolete after including nlte departure coefficients:    ! Calculates the derivative of the source function with respect to the temperature
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    USE FORWARD_PARAM, ONLY: NUML, PIXEL_END, PIXEL_INI &
! Obsolete after including nlte departure coefficients:        , WAVE, LINE_L0, IND_LINE
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    IMPLICIT NONE
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    INTEGER,    INTENT(IN)                  :: K,L
! Obsolete after including nlte departure coefficients:    REAL(DP), INTENT(OUT)                   :: DSF
! Obsolete after including nlte departure coefficients:    INTEGER                                 :: M
! Obsolete after including nlte departure coefficients:    REAL(DP)                                :: WAVELENGTH=0.D0
! Obsolete after including nlte departure coefficients:    REAL(DP), PARAMETER                     :: FACTOR = 8.5685641D-6 ! h^2*c^3/K
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    DO M=1,NUML
! Obsolete after including nlte departure coefficients:       !IF (L.GE.PIXEL_INI(M).AND.M.LE.PIXEL_END(M)) WAVELENGTH&
! Obsolete after including nlte departure coefficients:       !    =LINE_L0(IND_LINE(M))+WAVE(L)/1.0D3
! Obsolete after including nlte departure coefficients:       IF (L.GE.PIXEL_INI(M).AND.L.LE.PIXEL_END(M)) THEN
! Obsolete after including nlte departure coefficients:         WAVELENGTH = LINE_L0(IND_LINE(M))+WAVE(L)/1.0D3
! Obsolete after including nlte departure coefficients:         EXIT
! Obsolete after including nlte departure coefficients:       ENDIF
! Obsolete after including nlte departure coefficients:    ENDDO
! Obsolete after including nlte departure coefficients:    ! cm !
! Obsolete after including nlte departure coefficients:    WAVELENGTH=WAVELENGTH*1D-8
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:    DSF = REAL((2.0D0*FACTOR/(WAVELENGTH**6.0D0*TEM(K)**2.0D0))&
! Obsolete after including nlte departure coefficients:        *(EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(K)))&
! Obsolete after including nlte departure coefficients:        /(EXP(HPLA*LIGHT/(WAVELENGTH*KBOL*TEM(K)))-1.0D0)**2.0D0))
! Obsolete after including nlte departure coefficients:    !
! Obsolete after including nlte departure coefficients:  END SUBROUTINE SOURCE_FUNCTION_DER
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_BOUNDARY(K,L)!,SF,DSF)
    !
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN)                  :: K,L
    REAL(DP),   DIMENSION(4)                :: SF
!    REAL(SP),   DIMENSION(4,DER_ARGS)      :: DSF
    ! At the lower boundary we assume the emergent radiation is the 
    ! source function
    CALL SOURCE_FUNCTION(K,L,SF)
    SYN1D(:,L)=SF
    DSYN1D(:,:,L,K)=0.0D0
    !
  END SUBROUTINE SET_BOUNDARY
  !
  !================================================
  !
END MODULE RTESOLVER
!
