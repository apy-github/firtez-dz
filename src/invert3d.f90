!
MODULE INVERSION
  !
  !================================================
  !
  USE CODE_MODES, ONLY: TWODSPATIALSPD, COUPLED &
      , MSYNTHESIS, MRESPFUNCT, NAMEMODEL, NAMEPROFILE &
      , MGETTAU, MINVERSION, MGETHEQ, VREGULARIZATION &
      , MPROGRES
  USE INVERT_PARAM, ONLY: INV_ATMPAR, MAXITER, NFREQ &
      , CURIC, IMASK, AM_I_DONE, INV_MAS, INV_STK &
      , NSTKINV, STOP_CRIT, IFREEP, JACOB, DELTA &
      , ASSIST_B, ASSIST_T, ASSIST_P, ASSIST_V, MAXSTEPS &
      , SAMPLED_MOD, PENALTY, PEN_RES, PEN_HSS

  USE FORWARD, ONLY: FORWARD3D, UPDATE_BEST_TO_CURRENT, UPDATE_CURRENT_TO_BEST
  USE user_mpi, ONLY: mpi__myrank, mpi__ierror &
      , MPI_INTEGER, MPI_COMM_WORLD, MPI_LOGICAL
  USE MISC, ONLY: WRITE_MODEL, WRITE_PROFILES
  USE EXTEND_2D, ONLY: EXTEND_MODEL, SETUP_MASK
  USE GET_DMODEL, ONLY: GET_DMODEL3DC, GET_DMODEL3DS &
      , IFREE_VAR_SPACE, SET_WEIGHTS
  USE GRID_PARAM, ONLY: NZ,NY,NX, ZZ, YY, XX
  USE CONS_PARAM, ONLY: HPLA, LIGHT, KBOL, DP, SP
  USE PHYS_PARAM
  !
  USE FORWARD_PARAM, ONLY: NUMW, PIXEL_INI, PIXEL_END &
      , INDEX, WAVE, LINE_L0, HYDRO_TOP, IND_LINE
  !
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_DP,ALLOCATE_2D_DP &
      ,ALLOCATE_2D_SP,ALLOCATE_3D_DP,ALLOCATE_4D_SP
  !
  IMPLICIT NONE
  !
  LOGICAL :: WSET
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: INVERT3D
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! invert3d
  ! initialize_cycle
  ! smooth_model
  ! set_tau
  ! store_cycle
  ! solve_dmodel
  ! make_decision
  ! check_chi2_xy
  !
  !------------------------------------------------
  !
  SUBROUTINE INVERT3D()
    !
    INTEGER                   :: I
    LOGICAL                   :: STOP_ITERATION!, CRFS
    INTEGER                   :: CRFS
    INTEGER                   :: IITER
    INTEGER                   :: CNT
    !
    IF (mpi__myrank.EQ.0) PRINT*, 'INVERSION3D:'
    !
    WSET=.FALSE.
    IITER=1
    IF (VREGULARIZATION) THEN
      IF (MPROGRES.EQV..FALSE.) IITER=MAXITER
    ENDIF
    ! If we are inverting, maxiter is greater or equal 1
    CNT = 0
    DO I=IITER,MAXITER
      !
      CALL INITIALIZE_CYCLE(I,STOP_ITERATION,CRFS)
      !
      DO WHILE (STOP_ITERATION)
        !
        CALL MODELLING(CRFS)
        !
        ! We now calculate the perturbation to be applied...
        ! ...consistent with the previous modelizing:
        CALL SOLVE_DMODEL()
        !
        CALL MAKE_DECISION(STOP_ITERATION)
    CNT = CNT + 1
    !IF (CNT.GT.1) EXIT
        !
      ENDDO
      !
      ! To be moved to another subroutine:
      !
      CALL END_CYCLE(I)
      !
    ENDDO ! Inversion cycles
    !
    IF ((MSYNTHESIS.EQV..TRUE.) &
        .OR.(MGETTAU.EQV..TRUE.) &
        .OR.(MGETHEQ.EQV..TRUE.)) THEN
      ! First, we solve the forward model once.
      ! this is done always!!!! no matter if we are...
      ! ...synthesising, getting tau, or inverting
      CURIC=MAXITER
      IF (MRESPFUNCT.EQV..TRUE.) THEN
        CALL IFREE_VAR_SPACE()
        IF (mpi__myrank.eq.0) THEN
          ! IFREEP already updated. Various COEF, already updated.
          ! We do now allocate DSYN and BEST_DSYN
          CALL ALLOCATE_4D_SP(DSYN, NFREQ, IFREEP, NY, NX, 'ALLOCATE DSYN(I)')
        ELSE
          CALL ALLOCATE_2D_SP(EQVDSYN, NFREQ, IFREEP, 'ALLOCATE EQVDSYN(I)')
        ENDIF
      ENDIF
      CALL FORWARD3D()
    ENDIF
    !
  END SUBROUTINE INVERT3D
  !
  !------------------------------------------------
  !
  SUBROUTINE MODELLING (VCRF)
    !
    !USE FORWARD_PARAM, ONLY: N_FWD_MODELLING
    !
    INTEGER, INTENT(INOUT) :: VCRF
    !
    !IF (LOGIC.EQV..FALSE.) THEN
    IF (MOD(VCRF,1).EQ.0) THEN
      CALL FORWARD3D()
      !LOGIC=.FALSE.
    ELSE
      MRESPFUNCT=.FALSE.
      CALL FORWARD3D()
      IF (mpi__myrank.EQ.0) DSYN(:,:,:,:)=BEST_DSYN(:,:,:,:)
      MRESPFUNCT=.TRUE.
      !LOGIC=.TRUE.
    ENDIF
    VCRF=VCRF+1
    !
  END SUBROUTINE MODELLING
!><O><  SUBROUTINE MODELLING (LOGIC)
!><O><    !
!><O><    USE FORWARD_PARAM, ONLY: N_FWD_MODELLING
!><O><    USE INVERT_PARAM, ONLY: MAXSTEPS
!><O><    !
!><O><    LOGICAL, INTENT(INOUT) :: LOGIC
!><O><    !
!><O><    IF (LOGIC.EQV..FALSE.) THEN
!><O><      MRESPFUNCT=.FALSE.
!><O><      CALL FORWARD3D()
!><O><      IF (mpi__myrank.EQ.0) DSYN(:,:,:,:)=BEST_DSYN(:,:,:,:)
!><O><      MRESPFUNCT=.TRUE.
!><O><      LOGIC=.TRUE.
!><O><    ELSE
!><O><      CALL FORWARD3D()
!><O><    !><  IF (mpi__myrank.EQ.0) THEN
!><O><    !><    IF (SUM(INV_MAS(1,:,:))/(1.D0*NX*NY).LT.MAXSTEPS/2) LOGIC=.FALSE.
!><O><    !><  ENDIF
!><O><    !><  CALL MPI_BCAST(LOGIC,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi__ierror)
!><O><      LOGIC=.FALSE.
!><O><    ENDIF
!><O><    !
!><O><  END SUBROUTINE MODELLING
  !                                                               
  !------------------------------------------------
  ! 
  SUBROUTINE INITIALIZE_CYCLE(I,STOP_IT,CRF)
    !
    USE COUPLED_INVERSION, ONLY: INIT_COUPLED_INVERSION_DYNAMIC_VARS
    USE INVERT_PARAM, ONLY: INU
    USE EXTEND_2D, ONLY: SMOOTHING, SMOOTHING_B
    USE CODE_MODES, ONLY: MSMOOTHING
    !
    INTEGER, INTENT(IN) :: I
    INTEGER, INTENT(INOUT) :: CRF
    LOGICAL, INTENT(INOUT) :: STOP_IT!, CRF
    !
    !INTEGER                :: STEP, XINIT, YINIT
    !
    ! Set CRF:
    !CRF=.TRUE.
    CRF=0
    ! Set CURIC:
    CURIC=I
    !
    ! Update Cycle dimension dependent variables:
    CALL IFREE_VAR_SPACE()
    !
    INU=NFREQ!-IFREEP
    IF (INU.LT.0) THEN
      IF (mpi__myrank.EQ.1) PRINT* &
          , 'Warning: Number of free parameters larger than number of observations!'
      INU=ABS(INU)
    ENDIF
    !
    IF (mpi__myrank.eq.0) THEN
      !
      ! IFREEP already updated. Various COEF, already updated.
      ! We do now allocate DSYN and BEST_DSYN
      !
      CALL ALLOCATE_4D_SP(DSYN, NFREQ, IFREEP, NY, NX, 'ALLOCATE DSYN(I)')
      CALL ALLOCATE_4D_SP(BEST_DSYN, NFREQ, IFREEP, NY, NX &
          , 'ALLOCATE BEST_DSYN(I)')
      !
      PRINT*, 'NDSVD: ', NSTKINV*NUMW, ' ; MFIT: ', IFREEP
    ELSE ! Master. Slaves:
      CALL ALLOCATE_2D_SP(EQVDSYN, NFREQ, IFREEP, 'ALLOCATE EQVDSYN(I)')
      CALL ALLOCATE_2D_DP(JACOB,NFREQ,IFREEP,'JACOB IN GET_DMODEL')
      CALL ALLOCATE_1D_DP(DELTA,IFREEP,'DELTA IN GET_DMODEL')
      CALL ALLOCATE_1D_DP(SAMPLED_MOD,IFREEP,'SAMPLED_MOD IN GET_DMODEL')
      CALL ALLOCATE_1D_DP(PENALTY,IFREEP,'PENALTY IN GET_DMODEL')
      CALL ALLOCATE_1D_DP(PEN_RES,IFREEP,'PEN_RES IN GET_DMODEL')
      CALL ALLOCATE_2D_DP(PEN_HSS,IFREEP,IFREEP,'PEN_HSS IN GET_DMODEL')
    ENDIF ! Slaves.
    !
    IF (WSET.EQV..FALSE.) THEN
      CALL SET_WEIGHTS()
      WSET=.TRUE.
    ENDIF
    !
    STOP_IT=.TRUE.
    !
    AM_I_DONE(:,:)=0
    !
    IF (mpi__myrank.eq.0) THEN
      IMASK=REAL(1-AM_I_DONE)
      ! 1- STEPS_WITHOUT_IMPROVEMENT
      INV_MAS(1,:,:)=0.0D0
      ! 2- CURCHI:
      INV_MAS(2,:,:)=1.0D29
      ! 3- PRECHI:
      INV_MAS(3,:,:)=1.0D29
      ! 4- PRECYC
      INV_MAS(4,:,:)=0.0D0
      ! 9- LAMBDA
      INV_MAS(9,:,:)=1.0D+1
    ENDIF
    !
    IF (TWODSPATIALSPD.EQV..TRUE.) THEN
      ! IF WE SET THIS OPTION, NOT ALL THE PIXELS ARE INVERTED AT...
      ! ... EVERY CYCLE BUT ONLY A BUNCH OF THEM ARE.
      !
      CALL SETUP_MASK(CURIC,MAXITER)
      !
    ENDIF
    !
    IF (COUPLED.EQV..TRUE.) THEN
      CALL INIT_COUPLED_INVERSION_DYNAMIC_VARS()
    ENDIF
    !
    CALL SET_TAU(I)
    !
    CALL UPDATE_CURRENT_TO_BEST()
MSMOOTHING=.TRUE.
!MSMOOTHING=.FALSE.
IF (MSMOOTHING.EQV..TRUE.) THEN
  CALL SMOOTHING_B()
  CALL FORWARD3D()
ENDIF

    !
    !CALL STORE_CYCLE(I, 'pit_')
    !
  END SUBROUTINE INITIALIZE_CYCLE
  !
  !------------------------------------------------
  !
  SUBROUTINE END_CYCLE(I)
    !
    USE COUPLED_INVERSION, ONLY: END_COUPLED_INVERSION_DYNAMIC_VARS
    USE HEIGHT_HANDLER, ONLY: P0_COEFS, PG_COEFS, RH_COEFS, TM_COEFS &
        , BX_COEFS, BY_COEFS, BZ_COEFS, VZ_COEFS
    !
    INTEGER, INTENT(IN)      :: I
    !
    AM_I_DONE(:,:)=0
    !
    !PRINT*, 'Update best to current'
    CALL UPDATE_BEST_TO_CURRENT()
    !
    !CALL STORE_CYCLE(I, 'it_')
    !
    ! Here, we already have the final model. As a last step, we calculate the...
    ! ... errors:
    !CALL GET_ERRORS3DS()
    MRESPFUNCT=.FALSE.
!!    IF (mpi__myrank.EQ.0) PRINT*, 'BB', SUM(ABS(SYN3D)), 'BB'
    CALL FORWARD3D()
!!    IF (mpi__myrank.EQ.0) PRINT*, 'BB', SUM(ABS(SYN3D)), 'BB'
    MRESPFUNCT=.TRUE.
    !
    !CALL UPDATE_BEST_TO_CURRENT()
    !
    !CALL STORE_CYCLE(I, 'it_')
    !
    IF (ALLOCATED(P0_COEFS)) DEALLOCATE(P0_COEFS)
    IF (ALLOCATED(PG_COEFS)) DEALLOCATE(PG_COEFS)
    IF (ALLOCATED(RH_COEFS)) DEALLOCATE(RH_COEFS)
    IF (ALLOCATED(TM_COEFS)) DEALLOCATE(TM_COEFS)
    IF (ALLOCATED(BX_COEFS)) DEALLOCATE(BX_COEFS)
    IF (ALLOCATED(BY_COEFS)) DEALLOCATE(BY_COEFS)
    IF (ALLOCATED(BZ_COEFS)) DEALLOCATE(BZ_COEFS)
    IF (ALLOCATED(VZ_COEFS)) DEALLOCATE(VZ_COEFS)
    !
    IF (mpi__myrank.eq.0) THEN
      ! IFREEP already updated. Various COEF, already updated.
      ! We do now allocate DSYN and BEST_DSYN
      DEALLOCATE(DSYN)
      DEALLOCATE(BEST_DSYN)
      !
      ! Update current model to the best one:
      !MODEL2D_RCV(:,:,1)=MODEL2D_RCV(:,:,2)
      !
    ELSE
      DEALLOCATE(EQVDSYN)
      DEALLOCATE(JACOB)
      DEALLOCATE(DELTA)
      DEALLOCATE(SAMPLED_MOD)
      DEALLOCATE(PENALTY)
      DEALLOCATE(PEN_RES)
      DEALLOCATE(PEN_HSS)
    ENDIF

    IF (TWODSPATIALSPD.EQV..TRUE.) THEN
      ! IF WE SET THIS OPTION, NOT ALL THE PIXELS ARE INVERTED AT...
      ! ... EVERY CYCLE BUT ONLY A BUNCH OF THEM ARE.
      ! 11 -> NOT INVERTED
      ! 0  -> INVERTED
      !
      CALL EXTEND_MODEL(CURIC,MAXITER)
      !
      IF (CURIC.NE.MAXITER) THEN
        AM_I_DONE(:,:)=0
        CALL FORWARD3D()
      ENDIF
      !XINIT=MAXVAL((/STEP/2,1/))
      !YINIT=MAXVAL((/STEP/2,1/))
      !
      !AM_I_DONE(XINIT:NX:STEP,YINIT:NY:STEP)=0
      !
    ENDIF
    !
    IF (COUPLED.EQV..TRUE.) THEN
      CALL END_COUPLED_INVERSION_DYNAMIC_VARS()
    ENDIF
    !
  END SUBROUTINE END_CYCLE
  !
  !------------------------------------------------
  !
  SUBROUTINE SMOOTH_MODEL()
    !
    USE PHYS_PARAM
    USE FORWARD_PARAM, ONLY: N_FWD_MODELLING
    !
    USE USER_FFTW3, ONLY: FFTK2D, FFTY2D, GPLAN2D &
        , FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
    !
    INTEGER      :: J,K
    REAL(DP), DIMENSION(NZ,NY,NX) :: ARRAY
    !
    IF (COUPLED.EQV..TRUE.) THEN
      !
      !IF ((CURIC.EQ.1).AND.(MOD(N_FWD_MODELLING+1,3).EQ.0)) THEN
      IF (MOD(N_FWD_MODELLING+1,3).EQ.0) THEN
      !IF (CURIC.EQ.1) THEN
PRINT*, '*** Smoothing ***'
        !
        DO J=1,8
          !
          IF (INV_ATMPAR(J).EQV..FALSE.) CYCLE
          !
          SELECT CASE (J)
            CASE(1)
              ARRAY=DBLE(TEM3D)
            CASE(2)
              ARRAY=DBLE(PG3D)
            CASE(3)
              ARRAY=DBLE(RHO3D)
            CASE(4)
              ARRAY=DBLE(BX3D)
            CASE(5)
              ARRAY=DBLE(BY3D)
            CASE(6)
              ARRAY=DBLE(BZ3D)
            CASE(7)
              ARRAY=DBLE(VZ3D)
            CASE(8)
              ARRAY=DBLE(PG3D)
          END SELECT
          !
          DO K=1,NZ
            ! DFT synthetic data:
            CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),ARRAY(K,:,:),FFTY2D)
            !
            ! Multiply transformed arrays:
            FFTY2D=FFTY2D*FFTK2D
            ! DFT back the product:
            CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,ARRAY(K,:,:))
          ENDDO
          !
          SELECT CASE (J)
            CASE(1)
              TEM3D=REAL(ARRAY)
            CASE(2)
              PG3D=REAL(ARRAY)
            CASE(3)
              RHO3D=REAL(ARRAY)
            CASE(4)
              BX3D=REAL(ARRAY)
            CASE(5)
              BY3D=REAL(ARRAY)
            CASE(6)
              BZ3D=REAL(ARRAY)
            CASE(7)
              VZ3D=REAL(ARRAY)
            CASE(8)
              PG3D=REAL(ARRAY)
          END SELECT
          !
        ENDDO
        !
      ENDIF
      !
    ENDIF

  END SUBROUTINE SMOOTH_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_TAU(I)
    !
    USE PHYS_PARAM
    !
    INTEGER, INTENT(IN)  :: I
    !
    INTEGER      :: J,K
    LOGICAL      :: BU_INV, BU_SYN, BU_HYD, BU_RFC
    !
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: TAU1TEM
    REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: VECTORS
    REAL(DP)                 :: LK_HC, ATN2
    INTEGER                  :: ITZ, WQOFF, WUOFF, WVOFF
    REAL(DP), DIMENSION(1,1)  :: KOL
    !
    IF ((INV_ATMPAR(2).EQV..FALSE.).AND.(INV_ATMPAR(8).EQV..FALSE.)) THEN
IF (mpi__myrank.EQ.0) PRINT*, 'CHECK TAU'
      BU_SYN=MSYNTHESIS
      BU_INV=MINVERSION
      BU_HYD=HYDRO_TOP
      BU_RFC=MRESPFUNCT
      MSYNTHESIS=.FALSE.
      MINVERSION=.FALSE.
      MRESPFUNCT=.FALSE.
      MGETTAU=.TRUE.
      IF (I.EQ.1) THEN
      !
      !
      !
      !
      ! Assist B:
        IF (ASSIST_B.EQV..TRUE.) THEN
IF (mpi__myrank.EQ.0) PRINT*, '  Assist B:'
          IF (mpi__myrank.EQ.0) THEN
            IF (ALL(INV_STK(2:4)).EQV..TRUE.) THEN
              CALL ALLOCATE_3D_DP(TAU1TEM, NY, NX, 4, 'A')
              ITZ=PIXEL_END(1)-PIXEL_INI(1)+1
              CALL ALLOCATE_2D_DP(VECTORS, ITZ, 3, 'B')
              VECTORS(:,1)=SUM(SUM(OBS3D(PIXEL_INI(1):PIXEL_END(1),:,:) &
                  ,DIM=2),DIM=2)                                                                 
              VECTORS(:,1)=VECTORS(:,1)/SUM(ABS(VECTORS(:,1)))
              ! First derivative
              VECTORS(2:ITZ-1,2)=(VECTORS(3:ITZ,1)-VECTORS(1:ITZ-2,1))/2.0D0
              VECTORS(1,2)=VECTORS(2,2)
              VECTORS(ITZ,2)=VECTORS(ITZ-1,2)
              VECTORS(:,2)=VECTORS(:,2)/SUM(ABS(VECTORS(:,2)))
              ! Second derivative
              VECTORS(2:ITZ-1,3)=(VECTORS(3:ITZ,2)-VECTORS(1:ITZ-2,2))/2.0D0
              VECTORS(1,3)=VECTORS(2,3)
              VECTORS(ITZ,3)=VECTORS(ITZ-1,3)
              VECTORS(:,3)=VECTORS(:,3)/SUM(ABS(VECTORS(:,3)))
              !
              WQOFF=NUMW*1
              WUOFF=NUMW*2
              WVOFF=NUMW*3
              !
              DO J=1,NX
                DO K=1,NY
  ! StkV:
                  TAU1TEM(K,J,3)=SUM(OBS3D(PIXEL_INI(1)+3*NUMW:PIXEL_END(1) &
                      +3*NUMW,K,J)*VECTORS(:,2))
  ! StkQ:
                  TAU1TEM(K,J,1)=SUM(OBS3D(PIXEL_INI(1)+1*NUMW:PIXEL_END(1) &
                      +1*NUMW,K,J)*VECTORS(:,3))
  ! StkU:
                  TAU1TEM(K,J,2)=SUM(OBS3D(PIXEL_INI(1)+2*NUMW:PIXEL_END(1) &
                      +2*NUMW,K,J)*VECTORS(:,3))
                ENDDO
              ENDDO
              ! guess Bz accordingly
              !TAU1TEM(:,:,3)=TAU1TEM(:,:,3)/MAX(ABS(TAU1TEM(:,:,3)))
              TAU1TEM(:,:,3)=TAU1TEM(:,:,3)/MAXVAL(ABS(TAU1TEM(:,:,3)))
              TAU1TEM(:,:,4)=DSQRT(TAU1TEM(:,:,1)**2+TAU1TEM(:,:,2)**2)
              TAU1TEM(:,:,4)=TAU1TEM(:,:,4)/MAXVAL(TAU1TEM(:,:,4))
              DO J=1,NX
                DO K=1,NY
                  ATN2=ATAN2(TAU1TEM(K,J,2),TAU1TEM(K,J,1))/2.0D0
                  !ATN2=ATAN2(TAU1TEM(K,J,2),TAU1TEM(K,J,1))/2.0D0
                  BX3D(:,K,J)=-REAL(SIN(ATN2))*1000.E0*TAU1TEM(K,J,4)
                  BY3D(:,K,J)=REAL(COS(ATN2))*1000.E0*TAU1TEM(K,J,4)
                  BZ3D(:,K,J)=REAL(TAU1TEM(K,J,3))*(-1000.E0)
                ENDDO
              ENDDO
              !
              DEALLOCATE(TAU1TEM)
              DEALLOCATE(VECTORS)
              !
              BEST_BZ3D(:,:,:)=BZ3D(:,:,:)
              BEST_BY3D(:,:,:)=BY3D(:,:,:)
              BEST_BX3D(:,:,:)=BX3D(:,:,:)
              !
            ENDIF ! Invert polarized?
            !
          ENDIF ! Master
        ENDIF ! Assist B
        !
        !
        !
        !
        ! Assist Vlos:
        IF (ASSIST_V.EQV..TRUE.) THEN
IF (mpi__myrank.EQ.0) PRINT*, '  Assist Vlos:'
          IF (mpi__myrank.EQ.0) THEN
            CALL ALLOCATE_3D_DP(TAU1TEM, NY, NX, 1, 'A')
            ITZ=PIXEL_END(1)-PIXEL_INI(1)+1
            CALL ALLOCATE_2D_DP(VECTORS, ITZ, 2, 'B')
            VECTORS(:,1)=SUM(SUM(OBS3D(PIXEL_INI(1):PIXEL_END(1),:,:) &
                ,DIM=2),DIM=2)                                                                 
            VECTORS(:,1)=VECTORS(:,1)/SUM(ABS(VECTORS(:,1)))
            ! First derivative
            VECTORS(2:ITZ-1,2)=(VECTORS(3:ITZ,1)-VECTORS(1:ITZ-2,1))/2.0D0
            VECTORS(1,2)=VECTORS(2,2)
            VECTORS(ITZ,2)=VECTORS(ITZ-1,2)
            VECTORS(:,2)=VECTORS(:,2)/SUM(ABS(VECTORS(:,2)))
            ! They are almost ortogonal, so we approximate the coefficient by ...                                                                   
            ! ... the dot product
            !
            KOL=MATMUL(TRANSPOSE(VECTORS(:,2:2)), VECTORS(:,1:1))
            DO J=1,NX
              DO K=1,NY
! Di:
                TAU1TEM(K,J,1)=SUM(OBS3D(PIXEL_INI(1):PIXEL_END(1),K,J) &
                    *VECTORS(:,2))-SUM(OBS3D(PIXEL_INI(1):PIXEL_END(1),K,J) &
                    *VECTORS(:,1))*KOL(1,1)
              ENDDO
            ENDDO
            ! guess vz accordingly
            TAU1TEM(:,:,1)=TAU1TEM(:,:,1)/MAXVAL(ABS(TAU1TEM(:,:,1)))
            DO J=1,NX
              DO K=1,NY
                VZ3D(:,K,J)=REAL(TAU1TEM(K,J,1))*(-2.0E5)
              ENDDO
            ENDDO
            !
            DEALLOCATE(TAU1TEM)
            DEALLOCATE(VECTORS)
            !
            BEST_VZ3D(:,:,:)=VZ3D(:,:,:)
            !
          ENDIF ! Master
        ENDIF ! Assist Vlos
        !
        !
        !
        !
        !
        ! Assist Temperature:
        IF (ASSIST_T.EQV..TRUE.) THEN
IF (mpi__myrank.EQ.0) PRINT*, '  Assist T:'
          IF (mpi__myrank.EQ.0) THEN
            CALL ALLOCATE_3D_DP(TAU1TEM, NY, NX, 2, 'A')
            TAU1TEM(:,:,1)=SUM(DBLE(OBS3D(PIXEL_INI(1):PIXEL_END(1),:,:)) &
                ,DIM=1)/DBLE(PIXEL_END(1)-PIXEL_INI(1)+1)
            LK_HC = (LINE_L0(IND_LINE(1))*1.0D-8 * KBOL)/ (HPLA * LIGHT)
            TAU1TEM(:,:,2)=LK_HC * DLOG(1.0D0+1.0D0/TAU1TEM(:,:,1) &
                *(DEXP(1.0D0/6.39D3/LK_HC)-1.0D0))
            TAU1TEM(:,:,2)=1.0D0/TAU1TEM(:,:,2)
            !
            ITZ=CEILING(REAL(NZ)/3.)
            TAU1TEM(:,:,2)=TAU1TEM(:,:,2)-TEM3D(ITZ,:,:)
            TAU3D5(:,:,:)=1.
            TAU3D5(ITZ,:,:)=0.
            DO J=1,NX
              DO K=1,NY
                TEM3D(:,K,J)=TEM3D(:,K,J)+TAU1TEM(K,J,2)
              ENDDO
            ENDDO
            !
            BEST_TEM3D(:,:,:)=TEM3D(:,:,:)
            !
            DEALLOCATE(TAU1TEM)
          ENDIF ! Master
        ENDIF ! Assist T
        !
        !
        !
        !
        !
        ! Assist Gas pressure:
        IF (ASSIST_P.EQV..TRUE.) THEN
IF (mpi__myrank.EQ.0) PRINT*, '  Assist P:'
          IF (mpi__myrank.EQ.0) THEN
            PRINT*, 'SET_TAU'
            ITZ=CEILING(REAL(NZ)/3.)
            TAU3D5(:,:,:)=1.0E0
            TAU3D5(ITZ,:,:)=0.0E0
            DO J=1,NX
              DO K=1,NY
            !    PG3D(:,K,J)=PG3D(:,K,J)*1.333E5 &
            !      /PG3D(MINLOC(ABS(TAU3D5(:,K,J)), 1), K, J)
                PG3D(:,K,J)=PG3D(:,K,J)*1.333E5 &
                  /PG3D(ITZ, K, J)
              ENDDO
            ENDDO
            !
            BEST_PG3D(:,:,:)=PG3D(:,:,:)
            !
          ENDIF ! Master
        ENDIF ! Assist Pg
        !
        !
        !
        IF ( (ASSIST_T.EQV..TRUE.) .OR. (ASSIST_P.EQV..TRUE.) ) THEN
          !
          HYDRO_TOP=.FALSE.!BU_HYD
          CALL FORWARD3D()
          !
        ENDIF
        !
        HYDRO_TOP=BU_HYD
  ! Get profiles
  !PRINT*, 'Get profiles after assistance:'
        MSYNTHESIS=.TRUE.
        CALL FORWARD3D()
      ENDIF
!
      HYDRO_TOP=BU_HYD
      MSYNTHESIS=BU_SYN
      MINVERSION=BU_INV
      MRESPFUNCT=BU_RFC
    ENDIF
    !
  END SUBROUTINE SET_TAU
  !
  !------------------------------------------------
  !
  SUBROUTINE STORE_CYCLE(IT, TEXT)
    !
    INTEGER, INTENT(IN)     :: IT
    CHARACTER(*)            :: TEXT
    !
    CHARACTER(LEN=1)         :: ITCNT1
    CHARACTER(LEN=2)         :: ITCNT2
    CHARACTER(LEN=3)         :: ITCNT3
    INTEGER                  :: I, OFF
    !
    IF (mpi__myrank.eq.0) THEN

      IF (IT.LT.10) THEN
        WRITE (ITCNT1, "(I1)") IT
      ELSE IF (IT.LT.100) THEN
        WRITE (ITCNT2, "(I2)") IT
      ELSE
        WRITE (ITCNT3, "(I3)") IT
      ENDIF
      PRINT*, 'SAVING IT FINAL MODEL AND PROFILES:'

      IF (ALLOCATED(SYN)) DEALLOCATE(SYN)
      CALL ALLOCATE_4D_SP(SYN,4,NUMW,NY,NX,'SYN')
      OFF=0
      DO I=1,4
        IF (INV_STK(I).EQV..TRUE.) THEN
          SYN(I,:,:,:)=SYN3D(OFF+1:OFF+NUMW,:,:)
          OFF=OFF+NUMW
        ENDIF
      ENDDO
      IF (IT.LT.10) THEN
        CALL WRITE_PROFILES(TRIM(TEXT)//TRIM(ITCNT1)//'_'//TRIM(NAMEPROFILE) &
            , SHAPE(SYN), INDEX, WAVE, SYN, 3000, 3)
        DEALLOCATE(SYN)

        CALL WRITE_MODEL(TRIM(TEXT)//TRIM(ITCNT1)//'_'//TRIM(NAMEMODEL) &
            , NX, NY, NZ&
            , TEM3D, PG3D, RHO3D, BX3D, BY3D, BZ3D, VZ3D&
            , PEL3D, MW3D, TAU3D5, XX, YY, ZZ, 3000, 3)
      ELSE IF (IT.LT.100) THEN
        CALL WRITE_PROFILES(TRIM(TEXT)//TRIM(ITCNT2)//'_'//TRIM(NAMEPROFILE) &
            , SHAPE(SYN), INDEX, WAVE, SYN, 3000, 3)
        DEALLOCATE(SYN)

        CALL WRITE_MODEL(TRIM(TEXT)//TRIM(ITCNT2)//'_'//TRIM(NAMEMODEL) &
            , NX, NY, NZ&
            , TEM3D, PG3D, RHO3D, BX3D, BY3D, BZ3D, VZ3D&
            , PEL3D, MW3D, TAU3D5, XX, YY, ZZ, 3000, 3)
      ELSE
        CALL WRITE_PROFILES(TRIM(TEXT)//TRIM(ITCNT3)//'_'//TRIM(NAMEPROFILE) &
            , SHAPE(SYN), INDEX, WAVE, SYN, 3000, 3)
        DEALLOCATE(SYN)

        CALL WRITE_MODEL(TRIM(TEXT)//TRIM(ITCNT3)//'_'//TRIM(NAMEMODEL) &
            , NX, NY, NZ&
            , TEM3D, PG3D, RHO3D, BX3D, BY3D, BZ3D, VZ3D&
            , PEL3D, MW3D, TAU3D5, XX, YY, ZZ, 3000, 3)
      ENDIF
      !
    ENDIF
    !
    !
  END SUBROUTINE STORE_CYCLE
  !
  !------------------------------------------------
  !
  SUBROUTINE SOLVE_DMODEL()
    !
!><    REAL(DP)     :: AVGCUR, AVGPRE
    !
    IF (.NOT.COUPLED) THEN
!PRINT*, 'Not coupled!'
      CALL GET_DMODEL3DS()
    ELSE
      CALL GET_DMODEL3DC()
    ENDIF
    !
!><    IF (mpi__myrank.eq.0) THEN
!><      AVGCUR=SUM(INV_MAS(2,:,:)*IMASK)/REAL(NX*NY-SUM(AM_I_DONE))
!><      AVGPRE=SUM(INV_MAS(3,:,:)*IMASK)/REAL(NX*NY-SUM(AM_I_DONE))
!><      !IF ((AVGCUR.GT.9.9E28).AND.(AVGPRE.GT.9.9E28)) &
!><      !IF ((INV_MAS(2,1).GT.9.9E28).AND.(INV_MAS(3,1).GT.9.9E28)) &
!><      !    !DMODEL=DMODEL*0.E0
!><    ENDIF
    !
  END SUBROUTINE SOLVE_DMODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE MAKE_DECISION(STOP_IT)
    !
    LOGICAL, INTENT(INOUT) :: STOP_IT
    !
    IF (mpi__myrank.eq.0) THEN
      WRITE(*,*) '<STEPS>=', SUM(INV_MAS(1,:,:))/(1.D0*NX*NY)&
          ,';',MINVAL(INV_MAS(1,:,:)),';', MAXVAL(INV_MAS(1,:,:))
      WRITE(*,*) '<CHI>=', SUM(INV_MAS(3,:,:)*IMASK)/(SUM(IMASK))
      WRITE(*,*) '<LAMBDA>=', SUM(INV_MAS(9,:,:))/(1.D0*NX*NY)
      !
      ! Check for bad pixels. If any, replace it by the surrounding atmosphere
      CALL CHECK_CHI2_XY()
      !
      STOP_CRIT=0
      IF (NX*NY-SUM(AM_I_DONE).EQ.0) THEN
        STOP_CRIT=1
      ENDIF

      IF (SUM(INV_MAS(9,:,:)*(1.-AM_I_DONE))/SUM(1.0D0-AM_I_DONE).GE.1.D4) THEN
        STOP_CRIT=1
      ENDIF
!PRINT*, ' >>>>>>>> ', CURIC, MAXITER, MAXSTEPS, MAXSTEPS-(MAXITER-CURIC)
      IF (SUM(INV_MAS(1,:,:))/(1.D0*NX*NY).GE.MAX(1,(MAXSTEPS-(MAXITER-CURIC)))) THEN
      !IF (SUM(INV_MAS(1,:,:))/(1.D0*NX*NY).GE.MAXSTEPS) THEN
        STOP_CRIT=1
      ENDIF

    ENDIF ! Master
    !
    ! BROADCAST DECISION FROM MASTER TO SLAVES:
    CALL MPI_BCAST(STOP_CRIT,1,MPI_INTEGER,0,MPI_COMM_WORLD &
                ,mpi__ierror)
    ! BROADCAST MASK FROM MASTER TO SLAVES:
    CALL MPI_BCAST(AM_I_DONE,NX*NY,MPI_INTEGER,0,MPI_COMM_WORLD &
                ,mpi__ierror)
    !
    IF (STOP_CRIT.GT.0) STOP_IT=.FALSE.
    !
  END SUBROUTINE MAKE_DECISION
  !
  !------------------------------------------------
  !
  SUBROUTINE CHECK_CHI2_XY()
    !
    IF (mpi__myrank.eq.0) THEN
      !
  !    IF ( (NX.GT.1) .OR. (NY.GT.1) ) THEN

  !      

  !    ENDIF
      !
    ENDIF ! Master
    !
  END SUBROUTINE CHECK_CHI2_XY
  !
  !================================================
  !
END MODULE INVERSION
!
