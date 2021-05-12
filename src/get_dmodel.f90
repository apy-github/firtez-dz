!
MODULE GET_DMODEL
  !
  !================================================
  !
  USE USER_MPI
  !
  USE CODE_MODES, ONLY: VREGULARIZATION, INVERSIONTYPE, MVERBOSE
  USE CONS_PARAM, ONLY: SP, DP
  USE GRID_PARAM, ONLY: ZZ, NZ, NY, NX
  USE PHYS_PARAM
  USE FORWARD_PARAM, ONLY: ATM_ARGS, NUMW
  USE INVERT_PARAM, ONLY: INV_MAS, INV_SLA, NSLAB_PER_FREEV &
      , INV_ATMPAR, MAXSTEPS, INV_STK, NFREQ, NFREEP, AM_I_DONE &
      , MSGSIZE, MAXITER, ISIGMAP, WSTK, NSTKINV, ISIGMA, CURIC &
      , AUTOWEIGHT, NFREEV, YGUESS, YTOFIT, CYCPOW, NSLB_MAX &
      , JACOB, DELTA, ATM_FACTOR, MPERT, IFREEP, INU, RSVDTOL &
      , NJEVALS, TOFFSET, SAMPLED_MOD, PENALTY, PEN_RES &
      , REG_TYPE_FREEV, PEN_HSS, PEN_TYP, PEN_ALP
  !
  USE ATM_PARAM, ONLY: SPLIT_MODEL
  !
  ! Include COEFS:
  USE HEIGHT_HANDLER
  !
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_3D_DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: PEN_FACTOR
  REAL(DP), DIMENSION(8), PARAMETER :: PARAM_NORMS=(/6000.0d0, 1.0d5, 1.0d-5, 1.0d3, 1.0d3, 1.0d3, 1.0d6,6.0d2/)
  !                                                 tem     ,  pgas,    rho,    bx,    by,    bz,  vlos,   p0
  CHARACTER*4, DIMENSION(8), PARAMETER :: PARAM_LABELS=(/'Tem ', 'PGas', 'Rho ', 'Bx  ' &
      , 'By  ', 'Bz  ', 'Vlos','P0  '/)
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: GET_DMODEL3DC
  PUBLIC :: GET_DMODEL3DS
  PUBLIC :: IFREE_VAR_SPACE
  PUBLIC :: SET_WEIGHTS
  !PUBLIC :: PERT_TEMP
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! get_gdmodel
  ! update_chi2
  ! update_stokes_model_rf
  ! update_smrf_coupled
  ! update_smrf_noncoupled
  ! update_smrf_1d
  ! update_chi2_coupled
  ! update_chi2_noncoupled
  ! get_chi2_1d
  ! split_current_model
  ! join_both_models
  ! join_best_model
  ! split_new_model
  ! esplit_model3d
  ! setup_invmas
  ! set_isigmap
  ! ifree_var_space
  ! start_lm
  ! get_positions
  ! get_int_coefs
  ! update_coefs
  ! normalize_jacob
  ! new_solve_perturbation
  ! solve_perturbation
  ! get_chi
  ! expand_vector
  ! contract_vector
  ! limit_perturbation
  ! apply_perturbation
  ! update_lambda
  ! calculate_chi2
  ! new_calculate_chi2
  ! get_pert_1d
  ! get_error_1d
  ! new_solve_errors
  ! get_dmodel3ds
  ! get_errors3ds
  ! get_dmodel3dc
  ! update_stokes_model
  ! from_da_to_model
  ! new_sort_delta_to_atmosphere
  ! new_rearange_da_model
  ! set_weights
  ! set_weights3d
  ! setup_invmas_weight
  ! set_weights1d
  ! stddev
  !
  !------------------------------------------------
  !


  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE SPLIT_CURRENT_MODEL(J,I)
    !
    IMPLICIT NONE
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
       CASE DEFAULT
         CYCLE
       ENDSELECT
!IF (MAXVAL(MODEL2D_SND(:,K,1)).NE.MAXVAL(MODEL2D_SND(:,K,1))) THEN
IF (MAXVAL(ABS(MODEL2D_SND(:,K,1))).NE.MAXVAL(ABS(MODEL2D_SND(:,K,1)))) THEN
  PRINT*, 'Modified model NAN: ', MODEL2D_SND(:,K,1) &
      , ' ; K= ', K, ' ; J= ', J, ' ; I= ', I
ENDIF
       ARR (:,J,I) = MODEL2D_SND(:,K,1)

       NULLIFY(ARR)
    ENDDO
    ZZ(:)=MODEL2D_SND(:,13,1)
    !
    !
    DO K=1,ATM_ARGS-1
       SELECT CASE (K)
       CASE(1)
          ARR => BEST_TEM3D
       CASE(2)
          ARR => BEST_PG3D
       CASE(3)
          ARR => BEST_RHO3D
       CASE(4)
          ARR => BEST_PEL3D
       CASE(5)
          ARR => BEST_MW3D
       CASE(6)
          ARR => BEST_BX3D
       CASE(7)
          ARR => BEST_BY3D
       CASE(8)
          ARR => BEST_BZ3D
       CASE(9)
          ARR => BEST_VX3D
       CASE(10)
          ARR => BEST_VY3D
       CASE(11)
          ARR => BEST_VZ3D
       CASE(12)
          ARR => BEST_TAU3D5
       CASE DEFAULT
         CYCLE
       ENDSELECT
       ARR (:,J,I) = MODEL2D_SND(:,K,2)
IF (MAXVAL(ABS(MODEL2D_SND(:,K,2))).NE.MAXVAL(ABS(MODEL2D_SND(:,K,2)))) THEN
  PRINT*, 'Modified model NAN: ', MODEL2D_SND(:,K,2) &
      , ' ; K= ', K, ' ; J= ', J, ' ; I= ', I
ENDIF

       NULLIFY(ARR)
    ENDDO
    ZZ(:)=MODEL2D_SND(:,13,2)
    !
    !PRINT*, 'Uno: ', SUM(MODEL2D_SND(:,:,1),1), I, J
    !PRINT*, 'Dos: ', SUM(MODEL2D_SND(:,:,2),1), I, J

    !
  END SUBROUTINE SPLIT_CURRENT_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE JOIN_BOTH_MODELS(J,I)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: I, J
    INTEGER               :: K
    REAL(SP), POINTER     :: ARR(:,:,:)
    !
    ! Current model
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
       MODEL2D_SND(:,K,1) = ARR (:,J,I)
    !
       NULLIFY(ARR)
    ENDDO
    MODEL2D_SND(:,13,1) = ZZ(:)
    !
    ! Best model:
    !
    DO K=1,ATM_ARGS-1
       SELECT CASE (K)
       CASE(1)
          ARR => BEST_TEM3D
       CASE(2)
          ARR => BEST_PG3D
       CASE(3)
          ARR => BEST_RHO3D
       CASE(4)
          ARR => BEST_PEL3D
       CASE(5)
          ARR => BEST_MW3D
       CASE(6)
          ARR => BEST_BX3D
       CASE(7)
          ARR => BEST_BY3D
       CASE(8)
          ARR => BEST_BZ3D
       CASE(9)
          ARR => BEST_VX3D
       CASE(10)
          ARR => BEST_VY3D
       CASE(11)
          ARR => BEST_VZ3D
       CASE(12)
          ARR => BEST_TAU3D5
       ENDSELECT
       MODEL2D_SND(:,K,2) = ARR (:,J,I)
    !
       NULLIFY(ARR)
    ENDDO
    MODEL2D_SND(:,13,2) = ZZ(:)
    !
  END SUBROUTINE JOIN_BOTH_MODELS
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE ESPLIT_MODEL3D(I,J)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)   :: I, J
    INTEGER               :: K
    REAL(SP), POINTER     :: ARR(:,:,:)
    !
    DO K=1,ATM_ARGS-1
       SELECT CASE (K)
       CASE(1)
          ARR => ETEM3D
       CASE(2)
          ARR => EPG3D
       CASE(3)
          ARR => ERHO3D
       !CASE(4)
       !   ARR => EPEL3D
       !CASE(5)
       !   ARR => EMW3D
       CASE(6)
          ARR => EBX3D
       CASE(7)
          ARR => EBY3D
       CASE(8)
          ARR => EBZ3D
       CASE(9)
          ARR => EVX3D
       CASE(10)
          ARR => EVY3D
       CASE(11)
          ARR => EVZ3D
       !CASE(12)
       !   ARR => ETAU3D5
       CASE DEFAULT
         CYCLE
       ENDSELECT
       ARR (:,J,I) = MODEL1D_SND(:,K)
    !
       NULLIFY(ARR)
    ENDDO
    !ZZ(:)=MODEL1D_SND(:,13)
    !
  END SUBROUTINE ESPLIT_MODEL3D
  !
  !------------------------------------------------
  !
  SUBROUTINE SETUP_INVMAS(JJ,II)
    !
    INTEGER, INTENT(IN)   :: II, JJ
    !
    INTEGER               :: OFFSET
    !
    ! Observed
    OFFSET=10
    INV_MAS(OFFSET:MSGSIZE,JJ,II)=0.D0
    !
    ! Difference between observed and synthetic profiles to be sent:
    !
    INV_MAS(OFFSET:MSGSIZE,JJ,II)=OBS3D(:,JJ,II)-SYN3D(:,JJ,II)
    !
  END SUBROUTINE SETUP_INVMAS
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_ISIGMAP(INW)
  !
    INTEGER, INTENT(IN)   :: INW
  !
    INTEGER               :: IPAR, I
  !
  !
  IPAR=0
  ISIGMAP(:)=0.D0
  !
  DO I=1,4
    IF (INV_STK(I).EQV..TRUE.) THEN
      ISIGMAP(IPAR*INW+1:(IPAR+1)*INW)=WSTK(I)**2/ISIGMA**2
      IPAR=IPAR+1
    ENDIF
  ENDDO
  !
  END SUBROUTINE SET_ISIGMAP
  !
  !------------------------------------------------
  !
  INTEGER FUNCTION GET_NJEVALS(APAR)
    !
    INTEGER, INTENT(IN) :: APAR
    !
    IF (APAR.EQ.1) THEN
      ! If I am Temperature:
      IF (CURIC.EQ.0) THEN
        GET_NJEVALS=NZ
      ELSE IF (CURIC.EQ.MAXITER) THEN
        GET_NJEVALS=NZ
      ELSE
        GET_NJEVALS=MIN(NZ,INT(CYCPOW**(CURIC+TOFFSET-1)))
        !
      ENDIF
    ELSE
      ! If I am not Temperature:
      IF (CURIC.EQ.0) THEN
        GET_NJEVALS=NZ
      ELSE IF (CURIC.EQ.1) THEN
        GET_NJEVALS=1
      ELSE IF (CURIC.EQ.MAXITER) THEN
        GET_NJEVALS=NZ
      ELSE
        GET_NJEVALS=MIN(NZ,INT(CYCPOW**(CURIC-1)))
        !
      ENDIF
    ENDIF

    !!RETURN, GET_NJEVALS

  END FUNCTION
  !
  !------------------------------------------------
  !
  SUBROUTINE IFREE_VAR_SPACE()
    !
    INTEGER                  :: I
    INTEGER                  :: OFFSET
    INTEGER                  :: CNT_FREEV
    INTEGER                  :: INSLABS
    !
    OFFSET=0
    !
    MPERT=0.15D0
    IF (RSVDTOL.EQ.0) RSVDTOL=1.0D-4
    !
    IFREEP=OFFSET
    CNT_FREEV=1
    DO I=1,SIZE(INV_ATMPAR)
      !
      NJEVALS=GET_NJEVALS(I)
      !
      ! Warning, testing the possibility of having T with more...
      ! ...freedom from the beginning than the others.
!PRINT*, I, INV_ATMPAR(I), NSLB_MAX(I), NJEVALS, NZ, INSLABS
      CALL UPDATE_COEFS(I,INV_ATMPAR(I),NSLB_MAX(I),NJEVALS,NZ &
          ,INSLABS)

      IF(INV_ATMPAR(I).EQV..TRUE.) THEN
        NSLAB_PER_FREEV(CNT_FREEV)=INSLABS
        CNT_FREEV=CNT_FREEV+1
      ENDIF

    ENDDO
    !
  END SUBROUTINE IFREE_VAR_SPACE

  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_POSITIONS(INZ, INJ, IPOS)
    !
    INTEGER, INTENT(IN)              :: INZ, INJ
    DOUBLE PRECISION, INTENT(INOUT)  :: IPOS(INJ)
    !
    DOUBLE PRECISION                 :: STEP
    INTEGER                          :: I
    !
    STEP=DBLE(INZ-1)/DBLE(INJ-1)
    !
    DO I=1,INJ
      IPOS(I)=(DBLE(I)-1.0d0)*STEP+1.D0
    ENDDO
    !
  END SUBROUTINE GET_POSITIONS
  SUBROUTINE REFERENCE_GET_POSITIONS(INZ, INJ, IPOS)
    !
    INTEGER, INTENT(IN)              :: INZ, INJ
    DOUBLE PRECISION, INTENT(INOUT)  :: IPOS(INJ)
    !
    DOUBLE PRECISION                 :: STEP
    INTEGER                          :: I
    !
    !
    STEP=DBLE(INZ-1)/DBLE(INJ)
    !
    DO I=1,INJ
      IPOS(I)=(DBLE(I)-0.5D0)*STEP+1.D0
    ENDDO
    !
  END SUBROUTINE REFERENCE_GET_POSITIONS
  !
  !------------------------------------------------
  !
  SUBROUTINE REFERENCE_GET_INT_COEFS(INZ,INJ,ICOEFS)
    !
    INTEGER, INTENT(IN)                                    :: INZ, INJ
    !
    REAL(DP),DIMENSION(:,:), POINTER, INTENT(INOUT)  :: ICOEFS
    !
    REAL(DP), DIMENSION(INJ)                       :: RINTERVALS
    REAL(DP), DIMENSION(INJ+4)                     :: BRINTERVALS
    REAL(DP), DIMENSION(INZ, INJ+2)                :: RCOEFS
    REAL(DP), DIMENSION(INZ, INJ)                  :: OCOEFS
    !
    INTEGER                                                :: I, J
    !
!    PRINT*, '**GET_INT_COEFS**'
    !
    RINTERVALS(:)=0.D0
    BRINTERVALS(:)=0.D0
    RCOEFS(:,:)=0.D0
    OCOEFS(:,:)=0.D0
    !
    IF (INJ.LT.INZ) THEN
      !
      CALL REFERENCE_GET_POSITIONS(INZ, INJ, RINTERVALS)
      !
      BRINTERVALS(1)=-1.D0
      BRINTERVALS(2)=0.D0
      BRINTERVALS(3:3+INJ-1)=RINTERVALS(:)
      BRINTERVALS(INJ+3)=DBLE(INZ)+1.D0
      BRINTERVALS(INJ+4)=DBLE(INZ)+2.D0
      !
      RCOEFS(:,:)=0.D0
      !
      DO I=2,INJ+2+1
        DO J=1,INZ
          IF ((J.GE.BRINTERVALS(I-1)).AND.(J.LE.BRINTERVALS(I+1))) THEN
            IF (J.LT.BRINTERVALS(I)) THEN
              ! INTEGRATION BELOW
              RCOEFS(J,I-1)=(DBLE(J)-BRINTERVALS(I-1)) &
                  /(BRINTERVALS(I)-BRINTERVALS(I-1))
            ELSE
              ! INTEGRATION ABOVE
              RCOEFS(J,I-1)=(BRINTERVALS(I+1)-DBLE(J)) &
                  /(BRINTERVALS(I+1)-BRINTERVALS(I))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !
    ELSE IF (INZ.EQ.INJ) THEN
      ! IF NZ
      DO I=1,INZ
        RCOEFS(I,1+I)=1.D0
      ENDDO
    ELSE
      PRINT*, 'NJ CANNOT BE LARGER THAN NZ!'
      STOP
    ENDIF
    !
    DO I=1,INJ+2
      IF (INJ.EQ.1) THEN
        J=1
      ELSE
        IF (I.LE.2) THEN
          J=1
        ELSE IF (I.GE.(INJ+2-1)) THEN
          J=INJ+2-2
        ELSE
          J=I-1
        ENDIF
      ENDIF
      !
      OCOEFS(:,J)=OCOEFS(:,J)+RCOEFS(:,I)
    ENDDO
    !
    ICOEFS=OCOEFS
    !
  END SUBROUTINE REFERENCE_GET_INT_COEFS
  SUBROUTINE GET_INT_COEFS(INZ,INJ,ICOEFS)
    !
    INTEGER, INTENT(IN)                                    :: INZ, INJ
    !
    REAL(DP),DIMENSION(:,:), POINTER, INTENT(INOUT)  :: ICOEFS
    !
    REAL(DP), DIMENSION(INJ)                       :: RINTERVALS
    REAL(DP), DIMENSION(INJ+4)                     :: BRINTERVALS
    REAL(DP), DIMENSION(INZ, INJ+2)                :: RCOEFS
    REAL(DP), DIMENSION(INZ, INJ)                  :: OCOEFS
    !
    INTEGER :: I, J
    REAL(DP) :: JI
    REAL(DP) :: JM
    REAL(DP) :: JF
    !
!    PRINT*, '**GET_INT_COEFS**'
    !
    RINTERVALS(:)=0.D0
    BRINTERVALS(:)=0.D0
    RCOEFS(:,:)=0.D0
    OCOEFS(:,:)=0.D0
    !
    IF (INJ.EQ.1) THEN
      OCOEFS(:,1) = 1.0d0!DBLE(INJ) / DBLE(INZ)
    ELSE IF (INJ.LT.INZ) THEN
      !
      CALL GET_POSITIONS(INZ, INJ, RINTERVALS)
      !
      DO J=1,INJ
        JI = RINTERVALS(MAX(J-1,1))
        JF = RINTERVALS(MIN(J+1,INJ))
        JM = RINTERVALS(J)
!if (inj.eq.1)        print*, J, JI, JF
        DO I=1,INZ
          IF (I.LT.JI) CYCLE
          IF (I.GT.JF) CYCLE
!if (inj.eq.1) PRINT*, 'I=', I, 'JI=', JI, ' ; JM=', JM, ' ; JF=', JF
          IF (I.EQ.JM) THEN
            OCOEFS(I,J) = 1.0d0
!PRINT*, OCOEFS(I,J)
          ELSE IF (I.LT.JM) THEN
            OCOEFS(I,J) = (1.0d0-0.0d0) / (JM-JI) * (DBLE(I)-JI)
!PRINT*, OCOEFS(I,J)
          ELSE IF (I.GT.JM) THEN
            OCOEFS(I,J) = (1.0d0-0.0d0) / (JF-JM) * (JF-DBLE(I))
!PRINT*, OCOEFS(I,J)
          ENDIF
        ENDDO
!if (inj.eq.1)         PRINT*, OCOEFS(:,J)
!PRINT*, ''
      ENDDO
      !
    ELSE IF (INZ.EQ.INJ) THEN
      ! IF NZ
      DO I=1,INZ
        OCOEFS(I,I)=1.D0
      ENDDO
    ELSE
      PRINT*, 'NJ CANNOT BE LARGER THAN NZ!'
      STOP
    ENDIF
!!!!!    !
!!!!!    DO I=1,INJ+2
!!!!!      IF (INJ.EQ.1) THEN
!!!!!        J=1
!!!!!      ELSE
!!!!!        IF (I.LE.2) THEN
!!!!!          J=1
!!!!!        ELSE IF (I.GE.(INJ+2-1)) THEN
!!!!!          J=INJ+2-2
!!!!!        ELSE
!!!!!          J=I-1
!!!!!        ENDIF
!!!!!      ENDIF
!!!!!      !
!!!!!      OCOEFS(:,J)=OCOEFS(:,J)+RCOEFS(:,I)
!!!!!    ENDDO
    !
    ICOEFS=OCOEFS
    !
  END SUBROUTINE GET_INT_COEFS
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_COEFS(ATMIND, LOGIC, VMAX, INJ, INZ, NSLAB_PAR)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_2D_DP
    !
    INTEGER, INTENT(IN)                     :: ATMIND
    LOGICAL, INTENT(IN)                     :: LOGIC
    INTEGER, INTENT(IN)                     :: VMAX, INJ, INZ
    INTEGER, INTENT(INOUT)                  :: NSLAB_PAR
    !
    INTEGER                                 :: INJEVALS
    REAL(DP)                                :: STEP
    !
    REAL(DP),DIMENSION(:,:), POINTER  :: PTR_COEF
    !
!    PRINT*, '**UPDATE_COEFS**'
    STEP=0.D0
    !
!PRINT*, 'ATMIND, LOGIC', ATMIND, LOGIC
    IF (LOGIC.EQV..TRUE.) THEN
      INJEVALS=MIN(INJ,VMAX)
      STEP=DBLE(INZ)/DBLE(INJEVALS)
      ! COEFS:
      SELECT CASE (ATMIND)
      CASE(1)
       ! INJEVALS=MIN(INJ*CYCPOW,VMAX)
       ! STEP=DBLE(INZ)/DBLE(INJEVALS)
         CALL ALLOCATE_2D_DP(TM_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>TM_COEFS
      CASE(2)
         CALL ALLOCATE_2D_DP(PG_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>PG_COEFS
      CASE(3)
         CALL ALLOCATE_2D_DP(RH_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>RH_COEFS
      CASE(4)
         CALL ALLOCATE_2D_DP(BX_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>BX_COEFS
      CASE(5)
         CALL ALLOCATE_2D_DP(BY_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>BY_COEFS
      CASE(6)
         CALL ALLOCATE_2D_DP(BZ_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>BZ_COEFS
      CASE(7)
         CALL ALLOCATE_2D_DP(VZ_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>VZ_COEFS
      CASE(8)
         CALL ALLOCATE_2D_DP(P0_COEFS, INZ, INJEVALS &
          , 'ICOEFS IN SETUP_LM_LV1')
          PTR_COEF=>P0_COEFS
      ENDSELECT  
      !
      CALL GET_INT_COEFS(INZ, INJEVALS, PTR_COEF)
      !
IF (mpi__myrank.eq.0) PRINT*, IFREEP, INJEVALS
      IFREEP=IFREEP+INJEVALS
      NSLAB_PAR=INJEVALS
      !
      NULLIFY(PTR_COEF)
      !
    ENDIF
!IF (mpi__myrank.eq.0) PRINT*, ' [[ end UPDATE_COEFS ]]', IFREEP
!    STOP
!    STOP
!    STOP
    !
  END SUBROUTINE UPDATE_COEFS
  !
  !------------------------------------------------
  !
  SUBROUTINE NORMALIZE_JACOB()
    !
    INTEGER                           :: ITPAR
    INTEGER                           :: ITVAR
    INTEGER                           :: JINIT
    INTEGER                           :: JDELT
    REAL(DP)                          :: NORM
    REAL(DP)                          :: FACTOR
    !REAL(DP),DIMENSION(IFREEP,IFREEP) :: PRODMAT
    !
    INTEGER                           :: J
    !
    !
    ITVAR=0
    !
    ! Apply the square root of the noise/weights
    !
    DO ITPAR=1,IFREEP
      JACOB(:,ITPAR)=JACOB(:,ITPAR) * DSQRT(ISIGMAP)/DSQRT(INU)
    ENDDO
    !
    ! Calculate and apply the normalization factors:
    !
    JINIT=1
    DO ITPAR=1,SIZE(INV_ATMPAR)
      IF (INV_ATMPAR(ITPAR).EQV..TRUE.) THEN

        NORM = PARAM_NORMS(ITPAR)

        ITVAR=ITVAR+1
        JDELT=NSLAB_PER_FREEV(ITVAR)-1

        !PRINT*, ' I=', ITPAR, ' ; Norm=', NORM
        ATM_FACTOR(ITPAR)=1.D0*NORM

        JACOB(:,JINIT:JINIT+JDELT) &
            = JACOB(:,JINIT:JINIT+JDELT) * ATM_FACTOR(ITPAR)

        JINIT=JINIT+JDELT+1

      ENDIF ! Do I invert that parameter?
    ENDDO ! Inv. parameters


    !
  END SUBROUTINE NORMALIZE_JACOB
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE NEW_SOLVE_PERTURBATION()
    !
    !USE LM, ONLY: PGET_PERTURBATION
    USE LM, ONLY: GET_PERTURBATION
    !
    INTEGER       :: I, CNT
    REAL(DP)      :: FACTOR, AMINAN
    !
    FACTOR=1.0D0
    !
    IF (.NOT.VREGULARIZATION) THEN
      PEN_RES(:)=0.0D0
      PEN_HSS(:,:)=0.0D0
    ELSE
      PEN_RES=PEN_RES*PEN_FACTOR
      PEN_HSS=PEN_HSS*PEN_FACTOR
    ENDIF
    !
    CALL GET_PERTURBATION(IFREEP, NFREQ&
        , JACOB, INV_SLA(10:MSGSIZE) &
        , INV_SLA(9), DELTA &
        , PEN_RES, PEN_HSS &
        , RSVDTOL)
    !
    ! Check and filter out NaN
    AMINAN=SUM(DELTA)
    IF (AMINAN.NE.AMINAN) THEN
      PRINT*, 'Warning! I am ', mpi__myrank &
          , ' and my perturbation has NaN!'
      PRINT*, '    -> Filtering them out'
      CNT=0
      DO I=1,SIZE(DELTA)
        AMINAN=DELTA(I)
        IF (AMINAN.NE.AMINAN) THEN
          CNT=CNT+1
          DELTA(I)=0.D0
        ENDIF
      ENDDO
      PRINT*, '    -> ', CNT, 'Filtered!'
      STOP
    ENDIF
    !
    !
    ! Once used, I clean the jacobian:
    ! Remember that EQVDSYN stores the best equiv. response function.
    JACOB(:,:)=0.D0
    !
  END SUBROUTINE NEW_SOLVE_PERTURBATION
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE EXPAND_VECTOR(NJ, IVECTOR, COEFS, OVECTOR)
    !
    INTEGER, INTENT(IN)                        :: NJ
    REAL(DP), DIMENSION(NZ,NJ), INTENT(IN)     :: COEFS
    REAL(DP), DIMENSION(NJ), INTENT(IN)        :: IVECTOR
    REAL(DP), DIMENSION(NZ), INTENT(INOUT)     :: OVECTOR
    !
    INTEGER                                    :: K, J
    !
    OVECTOR(:)=0
    !
    DO K=1,NZ
      !
      DO J=1,NJ
        OVECTOR(K)=OVECTOR(K)+IVECTOR(J)*COEFS(K,J)
      ENDDO
    ENDDO
    !
  END SUBROUTINE EXPAND_VECTOR
  !
  !------------------------------------------------
  !
  SUBROUTINE CONTRACT_VECTOR(NJ, IVECTOR, COEFS, OVECTOR)
    !
    INTEGER, INTENT(IN)                        :: NJ
    REAL(DP), DIMENSION(NZ,NJ), INTENT(IN)     :: COEFS
    REAL(SP), DIMENSION(NZ), INTENT(IN)        :: IVECTOR
    REAL(DP), DIMENSION(NJ), INTENT(INOUT)     :: OVECTOR
    !
    INTEGER                                    :: K, J
    !
    OVECTOR(:)=0
    !
    DO J=1,NJ
      DO K=1,NZ
      !
        OVECTOR(J)=OVECTOR(J)+IVECTOR(K)*COEFS(K,J)
      ENDDO
    ENDDO
    OVECTOR(:)=OVECTOR(:)/SUM(COEFS,DIM=1)
    !
  END SUBROUTINE CONTRACT_VECTOR
  !
  !------------------------------------------------
  !
  SUBROUTINE LIMIT_PERTURBATION(PTR_MODEL2D,PTR_ITDELTA,PTR_ATMFACTOR,ERRORS)
    !
    !REAL(SP)                                :: IPGTOP
    !REAL(DP)                                :: ITDELTA
    !
    LOGICAL, INTENT(IN), OPTIONAL           :: ERRORS
    !
    INTEGER                                 :: IPAR
    INTEGER                                 :: VPAR, JPAR
    REAL(DP)                                :: ABSFACTOR
    !
    REAL(SP), POINTER, DIMENSION(:,:)       :: PTR_MODEL2D
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ITDELTA
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ATMFACTOR
    !
    IF (PRESENT(ERRORS)) THEN
    !><  PRINT*, '***ERRORS***'
      CALL REARANGE_DA_MODEL(PTR_MODEL2D, PTR_ITDELTA, PTR_ATMFACTOR,.TRUE.)
    ELSE
      CALL REARANGE_DA_MODEL(PTR_MODEL2D, PTR_ITDELTA, PTR_ATMFACTOR)
    ENDIF
    !
    !
  END SUBROUTINE LIMIT_PERTURBATION
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_LAMBDA(CHICUR, CHIPRE, LAMBDA, STEPS, PRECYC, IFACTOR)                                                            
    !
    REAL(DP), INTENT(INOUT)     :: CHICUR, CHIPRE, LAMBDA, STEPS, PRECYC                                                              
    REAL(DP), INTENT(IN), OPTIONAL :: IFACTOR
    !
    REAL(DP)                       :: FACTOR, LUPDATE, LMIN
    !
    LMIN=1.0D-2
    !
    FACTOR=5.D-2
    IF (PRESENT(IFACTOR)) FACTOR=DBLE(IFACTOR)
!
    IF ((CHIPRE-CHICUR).GT.0) THEN
      ! CHI IMPROVES:
      ! WE ONLY UPDATE LAMBDA IF THE PREVIOUS WAS AN IMPROVEMENT:                                                                  
      IF ((CHIPRE-CHICUR).GT.FACTOR*CHIPRE) THEN
        STEPS=MAXVAL((/0.D0, STEPS-1.D0/))
      ELSE
        STEPS=STEPS+1.0D0
      ENDIF
      LUPDATE=DSQRT(10.D0)
      CHIPRE=CHICUR

      LAMBDA=LAMBDA/LUPDATE
    ELSE
      IF (ABS(PRECYC-1.D0).LT.1.D-2) THEN
        ! CHI DOES NOT IMPROVE BUT IN THE PREVIOUS STEP, WAS IMPROVING!                                                               
        PRECYC=0.D0
! WE CONTINUE DECREASING LAMBDA TO CHECK THAT WE LOST THE THING
!        LAMBDA=LAMBDA/10.D0!*2.D-1
! WE COME BACK TO THE PREVIOUS LAMBDA TO CHECK THAT WE LOST THE THING                                                                 
        LAMBDA=LAMBDA*3.0D0!*2.D-1
      ELSE
        LAMBDA=LAMBDA*10.0D0
      ENDIF
      STEPS=STEPS+1.D0
    ENDIF
    LAMBDA=MAX(LMIN, LAMBDA)
    !
  END SUBROUTINE UPDATE_LAMBDA
  !
  !------------------------------------------------
  !
  SUBROUTINE NEW_CALCULATE_CHI2()
    !
    ! 1) Use OBS and current SYN to estimate X2:
    INV_SLA(2)=SUM((ROBS1D(:)-RSYN1D(:))**2*ISIGMAP)/INU
    !
    ! 1b) Add reg term?
    IF (VREGULARIZATION) THEN
      PEN_FACTOR=1.0d0
      !IF (SUM(PEN_RES*PEN_RES).LT.(INV_SLA(2)*1.0D-2)) THEN
      !  PEN_FACTOR=DSQRT(INV_SLA(2)/SUM(PEN_RES*PEN_RES)*5.0D-1)
      !ENDIF
      !
      IF (MVERBOSE.GT.0) THEN
        IF (MVERBOSE.LT.3) THEN
            PRINT*, ' Chi: ', INV_SLA(2), ' ; Reg. Chi: ', SUM(PEN_RES*PEN_RES)! &
      !          , ' ; MReg. Chi: ', SUM(PEN_RES*PEN_FACTOR*PEN_FACTOR*PEN_RES)
        ELSE
            PRINT*, ' Chi: ', INV_SLA(2), ' ; Reg. Chi: ', SUM(PEN_RES*PEN_RES) &
                , ' ; MReg. Chi: ', SUM(PEN_RES*PEN_FACTOR*PEN_FACTOR*PEN_RES) &
                , ' ; P(factor)=', PEN_FACTOR
        ENDIF
      ENDIF
      !
      INV_SLA(2)=INV_SLA(2)+SUM(PEN_RES*PEN_FACTOR*PEN_FACTOR*PEN_RES)
      !
    ENDIF
    !
    ! Once we know the curren Chi2, we can compare with the previous one ...
    ! ... and so, update, if needed, the model and RF and update lambda ...
    ! ... of this pixel:
    !
    IF (INV_SLA(2).LT.INV_SLA(3)) THEN
      !
      ! Update best model. Remember that MODEL is a 3D array ...
      ! ... with NZxATM_ARGSx2, with the current model in the ...
      ! ... first position of the third index and the best one ...
      ! ... in the second position of the last index:
      !
      MODEL2D_RCV(:,:,2)=MODEL2D_RCV(:,:,1)
      !
      ! Update best equivalent response function. Remember that ...
      ! ... JACOB stores the current equiv. response function and ...
      ! ... EQVDSYN carries the best equiv. response function.
      !
      EQVDSYN(:,:)=REAL(JACOB(:,:))
      !
      ! Update best synthetic Stokes profiles:
      BEST_RSYN1D(:)=RSYN1D(:)
      !
    ELSE
      JACOB(:,:)=DBLE(EQVDSYN(:,:))
      MODEL2D_RCV(:,:,1)=MODEL2D_RCV(:,:,2)
      !
      ! If using regularization, update penalties:
      CALL GET_IT_PENALTIES()
      !
    ENDIF
    !
    ! Now, update lambda according to the improvement of the fit.
    CALL UPDATE_LAMBDA(INV_SLA(2), INV_SLA(3) &
        , INV_SLA(9), INV_SLA(1), INV_SLA(4) &
        , 0.01D0)
        !, 0.01D0)
!(1.0D0-DBLE(CURIC)/DBLE(MAXITER)*0.0D0)*1.D-2)
    !
    ! Now, we calculate the difference between OBS and BEST_SYN to...
    ! ... be used for the delta calculation:
    INV_SLA(10:MSGSIZE)=(ROBS1D(:)-BEST_RSYN1D(:))*DSQRT(ISIGMAP)/DSQRT(INU)
!PRINT*, 'diff', (ROBS1D(:)-BEST_RSYN1D(:))
!PRINT*, 'sigma, nu', DSQRT(ISIGMAP), DSQRT(INU)
    !
    !
  END SUBROUTINE NEW_CALCULATE_CHI2
  !
  !------------------------------------------------
  !


  !
  !
  !
  SUBROUTINE GET_IT_PENALTIES()
    !
    INTEGER :: JINIT
    INTEGER :: ITPAR
    INTEGER :: ITVAR
    INTEGER :: JI
    INTEGER :: JD
    INTEGER :: JF
    INTEGER :: K
    INTEGER :: CNT
    INTEGER :: NPEN
    REAL(DP) :: NORM
    REAL(DP) :: ALP
    !
    REAL(DP),DIMENSION(:,:), POINTER  :: PTR_COEF
    REAL(SP),DIMENSION(:), POINTER    :: PTR_MODEL
    REAL(SP),DIMENSION(:), POINTER    :: PTR_ZMODEL
    !
    REAL(DP),DIMENSION(:,:), ALLOCATABLE  :: IT_LACOB
    REAL(DP),DIMENSION(:), ALLOCATABLE  :: ITZ
    !
    REAL(DP) :: GDZ
    REAL(DP) :: DZ
    !
    !
    PENALTY(:)=1.0d0
    !
    JI=1
    ITVAR=0
    DO ITPAR=1,SIZE(INV_ATMPAR)
      IF (INV_ATMPAR(ITPAR).EQV..TRUE.) THEN
        ITVAR=ITVAR+1

        !PTR_ZMODEL => MODEL2D_RCV(:,13,1)
        PTR_ZMODEL => MODEL2D_RCV(:,12,1) ! Use tau instead
        SELECT CASE(ITPAR)
          CASE(1)
             PTR_COEF => TM_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,1,1)
          CASE(2)
             PTR_COEF => PG_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,2,1)
          CASE(3)
             PTR_COEF => RH_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,3,1)
          CASE(4)
             PTR_COEF => BX_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,6,1)
          CASE(5)
             PTR_COEF => BY_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,7,1)
          CASE(6)
             PTR_COEF => BZ_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,8,1)
          CASE(7)
             PTR_COEF => VZ_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,11,1)
          CASE(8)
             PTR_COEF => P0_COEFS(:,:)
             PTR_MODEL => MODEL2D_RCV(:,2,1)
        ENDSELECT  
        NORM = PARAM_NORMS(ITPAR)
        ALP=DSQRT(DBLE(PEN_ALP(ITPAR))/DBLE(IFREEP))

        JD = NSLAB_PER_FREEV(ITVAR)
        JF=JI+JD-1

        IF (MVERBOSE.GT.6) THEN
          SAMPLED_MOD(JI:JF) = 0
          CALL CONTRACT_VECTOR(JD, MODEL2D_RCV(:,12,1), PTR_COEF, SAMPLED_MOD(JI:JF))
          CNT=0
          print*, ' ***Tau: *** '!, JI, JF, ' ; ', SAMPLED_MOD(JI:JF)
          DO K=JI,JF
            IF (SAMPLED_MOD(K).LT.(-4)) CYCLE
            IF (SAMPLED_MOD(K).GT.0.5) CYCLE
            CNT=CNT+1
            PRINT*, SAMPLED_MOD(K)
          ENDDO
          PRINT*, ' Number of nodes inside ltauE(0.5,-4)=', CNT, ' out of a maximum of ', JF-JI+1 &
              , ' nodes for ', PARAM_LABELS(ITPAR)
          SAMPLED_MOD(JI:JF) = 0
        ENDIF

        CALL CONTRACT_VECTOR(JD, PTR_MODEL, PTR_COEF, SAMPLED_MOD(JI:JF))
        SAMPLED_MOD(JI:JF)=SAMPLED_MOD(JI:JF)/NORM

        NPEN=JD
        ALLOCATE(IT_LACOB(NPEN,JD))
        PENALTY(JI:JF)=1.0d0! / dble(JD)

        IT_LACOB(:,:)=0.0D0
        PEN_RES(JI:JF)=0.0D0
        SELECT CASE(PEN_TYP(ITPAR))
          CASE(0)
            PENALTY(JI:JF)=1.0d0 / dble(JD)
            DO K=1,NPEN
              IT_LACOB(K,K)=1.0D0*PENALTY(K+JI-1)
            ENDDO
            ! Difference to 0:
!            PENALTY(JI:JF)=PENALTY(JI:JF)*(SAMPLED_MOD(JI:JF)-0.0D0)
          CASE(1)
            ! Difference to next node:
            IF (JD.GT.1) THEN
              ! Second, fill the lacobian for this physical parameter:
              DO K=2,JD
                IT_LACOB(K,K)=1.0D0*PENALTY(K+JI-1)
                IT_LACOB(K,K-1)=-1.0D0*PENALTY(K+JI-1)
              ENDDO              ! First, estimate the penalty for each node:

            ENDIF
 
        ENDSELECT

        ! Now, estimate the residual term of the regularization term:
        IT_LACOB=IT_LACOB*ALP
        !
        !
        ! Testing:
        PENALTY(JI:JF) = MATMUL(IT_LACOB, SAMPLED_MOD(JI:JF))
        ! Testing.
        ! 
        ! 
        PEN_RES(JI:JF)=MATMUL(TRANSPOSE(IT_LACOB), PENALTY(JI:JF))
        PEN_HSS(JI:JF,JI:JF)=MATMUL(TRANSPOSE(IT_LACOB), IT_LACOB)
!PRINT*, ':::::', shape(IT_LACOB), ';', shape(PEN_HSS(JI:JF,JI:JF))
!!!!        PRINT*, ' Alpha; JI, JF'
!!!!        PRINT*, ALP, JI, JF
!!!!!        PRINT*, ' Model:'
!!!!!        PRINT*, SAMPLED_MOD(JI:JF)
!!!!        PRINT*, ' Penalty:', PENALTY
!!!!!        PRINT*, ' Penalty:'
!!!!!        PRINT*, SUM(PENALTY(JI:JF)*PENALTY(JI:JF))
!!!!!        PRINT*, ' Lesidue:'
!!!!!        PRINT*, PEN_RES
!!!!        PRINT*, ' Lacobian:'
!!!!        DO K=1,NPEN
!!!!          PRINT*, IT_LACOB(K,:)
!!!!        ENDDO
!!!!        PRINT*, ' Lessian:'
!!!!        DO K=1,SIZE(SAMPLED_MOD)
!!!!          PRINT*, PEN_HSS(max(1,K-3):min(K+3,SIZE(SAMPLED_MOD)),K)
!!!!        ENDDO
        DEALLOCATE(IT_LACOB)
!!!!!        DEALLOCATE(ITZ)

        JI=JI+JD

        NULLIFY(PTR_COEF)
        NULLIFY(PTR_MODEL)
        NULLIFY(PTR_ZMODEL)
      ENDIF ! Do I invert that parameter?
    ENDDO ! Inv. parameters
!!!!    PRINT*, ' sampled mod:', SAMPLED_MOD
!print*, ' Check if this is working!!'
!stop
    !
  END SUBROUTINE GET_IT_PENALTIES
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_PERT_1D()
    !
USE ATM_PARAM

    REAL(SP), POINTER, DIMENSION(:,:)       :: MODEL2D
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ITDELTA
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ATMFACTOR
    !
    !

! First, estimate the penalties for this iteration:

    CALL GET_IT_PENALTIES()


    ! First, we estimate the chi2:
    CALL NEW_CALCULATE_CHI2()
    ! From here on, both models stored in MODEL2D_RCV are the same
    !
    ! Now that we have already update the model and the equivalent ...
    ! ... response function (if needed), we split the model into ...
    ! ... the various atmospheric physical parameters:
    CALL SPLIT_MODEL(MODEL2D_RCV(:,:,2))

    !
    ! The equivalent response functions are multiplied by the square root...
    ! ... of sigma and also, normalized by the various atmospheric factors:
    !
    CALL NORMALIZE_JACOB()
    !
    CALL NEW_SOLVE_PERTURBATION()
    !
    MODEL2D=>MODEL2D_RCV(:,:,1)
    PTR_ITDELTA=>DELTA(:)
    PTR_ATMFACTOR=>ATM_FACTOR(:)
    !
    CALL LIMIT_PERTURBATION(MODEL2D,PTR_ITDELTA,PTR_ATMFACTOR)
    !
    NULLIFY(MODEL2D)
    NULLIFY(PTR_ITDELTA)
    NULLIFY(PTR_ATMFACTOR)
    !
    !
  END SUBROUTINE GET_PERT_1D
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_ERROR_1D(MODEL)
    !
    REAL(SP), INTENT(INOUT), TARGET :: MODEL(NZ,ATM_ARGS)
    !
    REAL(SP), POINTER, DIMENSION(:,:)       :: MODEL2D
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ITDELTA
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ATMFACTOR
    !
    ! Should be removed in the near future
    !CALL START_LM()
    !
    ! First, we estimate the chi2:
    !
    CALL NEW_CALCULATE_CHI2()
    !
    ! Now  we split the model into ...
    ! ... the various atmospheric physical parameters:
    !
    ! The equivalent response functions are multiplied by the square root...
    ! ... of sigma and also, normalized by the various atmospheric factors:
    !
    CALL NORMALIZE_JACOB()
    !
    CALL NEW_SOLVE_ERRORS(INV_SLA(2)/INU)
    !
    MODEL(:,:)=0.E0
    MODEL2D=>MODEL
    PTR_ITDELTA=>DELTA(:)
    PTR_ATMFACTOR=>ATM_FACTOR(:)
    !
    CALL LIMIT_PERTURBATION(MODEL2D,PTR_ITDELTA,PTR_ATMFACTOR,.TRUE.)
    !
!    PRINT*, SUM(MODEL,1)
    NULLIFY(MODEL2D)
    NULLIFY(PTR_ITDELTA)
    NULLIFY(PTR_ATMFACTOR)
!    PRINT*, SUM(MODEL,1)
    !
  END SUBROUTINE GET_ERROR_1D
  !
  !------------------------------------------------
  !
  SUBROUTINE NEW_SOLVE_ERRORS(NUME)
    !
    USE LM
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_2D_DP
    !
    REAL(DP),INTENT(IN)                       :: NUME
    !
    REAL(DP), ALLOCATABLE, DIMENSION(:,:)     :: HESS, COV
    REAL(SP)                                  :: STD
    INTEGER                                   :: I
    !
!    PRINT*, 'RCHI2: ', NUME, ISIGMA
    !
    CALL ALLOCATE_2D_DP(HESS, IFREEP, IFREEP, 'HESS')
    CALL ALLOCATE_2D_DP(COV, IFREEP, IFREEP, 'IHESSCOV')
    CALL HESSIAN(IFREEP, NFREQ, TRANSPOSE(JACOB), HESS)
    !
    CALL GET_INVERSE_COV(IFREEP, HESS, COV)
    !
    DO I=1,IFREEP
      !><DELTA(I)=DSQRT(NUME*COV(I,I))
      IF (HESS(I,I).GT. 1.D-32) THEN
        DELTA(I)=DSQRT(NUME/HESS(I,I))
      ELSE
        DELTA(I)=1.D32
      ENDIF
!      PRINT*, 'I: ', I, 'DELTA(I): ', DELTA(I), DSQRT(NUME*COV(I,I)), DSQRT(NUME/HESS(I,I))!, NUME &
!          , ';', 1.D0/HESS(I,I), DSQRT(IHESSCOV(I,I))
    ENDDO
    !
    DEALLOCATE(HESS)
    DEALLOCATE(COV)
    !
  END SUBROUTINE NEW_SOLVE_ERRORS
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_DMODEL3DS()
  !
    !
    INTEGER              :: I, J, L, M
    INTEGER              :: RJ, RI
    INTEGER              :: max_work
    INTEGER              :: NMAX_FWD
    !
    IF (mpi__myrank .EQ. 0) PRINT*, 'GET_DMODEL(START): ', NX*NY, SUM(AM_I_DONE)
    !IF (mpi__myrank .EQ. 0) DMODEL(:,:,:)=0.D0
    ! WORK DIVISION:
    ! INITIALISE X AND Y INDEXES
    I=1
    J=1
    ! LOOP THE THROUGH THE WHOLE DATA IN STEPS OF WORKING SLAVES (MASTER DOES NOT WORK)
    NMAX_FWD=NX*NY-SUM(AM_I_DONE)
    DO L=1,NMAX_FWD,mpi__size-1
      !
      ! CHECK IF WE ARE IN THE LAST X-Y ITERATION (IF IT IS THE LAST, IT MIGHT BE THAT
      ! NOT ALL THE SLAVES HAS SOMETHING TO DO
      ! max_work IS THE VARIABLE THAT ENCODES THIS INFO
      max_work=mpi__size-1
      IF ( (NMAX_FWD-L+1).LT.mpi__size) THEN
         max_work=NMAX_FWD-L+1
      ENDIF
      !
      ! ONLY THE PROCS THAT HAS TO WORK IN THIS L (X-Y) ITERATION MUST DO SOMETHING:
      IF (mpi__myrank.LE.max_work) THEN
        !
        ! MASTER SENDS I,J ATMOSPHERE MODELS TO SLAVES AND WAIT FOR THE STOKES AND DERIVATIVES:
        IF (mpi__myrank .EQ. 0) THEN
          !
          ! LOOP THROUGH THE NUMBER OF SLAVES WORKING
          RJ=J
          RI=I
!WRITE(*,*) 'Dmodel: ', I, ' of ', NX, ' and ', J, ' of ', NY

          M=0
          DO WHILE (M.LT.max_work)
            IF (MOD(J,NY+1).EQ.0) THEN
              I=I+1
              J=1
            ENDIF
            J=J+1
            IF (AM_I_DONE(J-1,I).GT.0.5) CYCLE
!PRINT*, 'DMODEL: ', I, J-1
            ! NOW, IF WE PASS THIS CLAUSE, WE MUST ACCESS J-1 IN THE FOLLOWING!
            ! M STARTS AT 0, NOT 1, SO ADD +1
            !
            ! Send Stokes profiles (both observed and synthetic) as well ...
            ! ... as some inversion parameters:
            !
            !CALL SETUP_INVMAS(J-1,I)
!            CALL MPI_SEND(INV_MAS(:,J-1,I),MSGSIZE &
!                ,MPI_DOUBLE_PRECISION,M+1,1,MPI_COMM_WORLD &
!                ,mpi__ierror)
            ! 
            CALL MPI_SEND(OBS3D(:,J-1,I),NFREQ &
                ,MPI_DOUBLE_PRECISION,M+1,52,MPI_COMM_WORLD &
                ,mpi__ierror)
            ! 
            CALL MPI_SEND(SYN3D(:,J-1,I),NFREQ &
                ,MPI_DOUBLE_PRECISION,M+1,53,MPI_COMM_WORLD &
                ,mpi__ierror)
            ! 
            CALL MPI_SEND(BEST_SYN(:,J-1,I),NFREQ &
                ,MPI_DOUBLE_PRECISION,M+1,54,MPI_COMM_WORLD &
                ,mpi__ierror)
            ! 
            ! Send current response functions:
            CALL MPI_SEND(INV_MAS(1:9,J-1,I),9 &
                ,MPI_DOUBLE_PRECISION,M+1,1,MPI_COMM_WORLD &
                ,mpi__ierror)
            ! 
            ! Send current response functions:
            ! Send current response functions:
            !
            CALL MPI_SEND(DSYN(:,:,J-1,I),NFREQ*IFREEP&
                ,MPI_REAL,M+1,2,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Send best response functions:
            !
            CALL MPI_SEND(BEST_DSYN(:,:,J-1,I),NFREQ*IFREEP&
                ,MPI_REAL,M+1,3,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Send SIGMA:
            !
            CALL MPI_SEND(ISIGMAP3D(:,J-1,I),NFREQ,MPI_DOUBLE_PRECISION &
                ,M+1,4,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Finally, we also need to send both, current atmosphere model:
            ! Encode a I, J set of thermodynamic variables in a 2d array to be send
            !
            CALL JOIN_BOTH_MODELS(J-1,I)
            !
            ! Send atmospheric parameters:
            CALL MPI_SEND(MODEL2D_SND,2*NZ*ATM_ARGS,MPI_REAL,M+1,5,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! UPDATE M
            M=M+1
            !
          ENDDO
          ! RESET i AND j TO RECOVER SLAVE OUTPUT:
          J=RJ
          I=RI
          ! THIS WEIRD THING IS JUST NEEDED AS LONG AS WE AVOID ASYNCHRONICITY
          ! WHENEVER WE PROCEED WITH THAT PROPERTY, AN ASYNCHRONOUS SEND, AVOIDS
          ! BLOCKING COMMUNICATIONS AND SO, IT IS EASIER TO IMPLEMENT.
          ! YET, WE DO NOT USE THIS THING FOR SOME REASON I AM NOT CONVINCED OF.
          M=0
          DO WHILE (M.LT.max_work)
            IF (MOD(J,NY+1).EQ.0) THEN
              I=I+1
              J=1
            ENDIF
            J=J+1
            IF (AM_I_DONE(J-1,I).GT.0.5) CYCLE
            !
            ! NOW, IF WE PASS THIS CLAUSE, WE MUST ACCESS J-1 IN THE FOLLOWING!
            ! NOW, M+1 IS THE SENDER, NOT M (IT STARTS WITH 0, MASTER)
            ! RECEIVE MODEL PERTURBATION:
            CALL MPI_RECV(MODEL2D_SND,2*NZ*ATM_ARGS,MPI_REAL,M+1,6 &
                ,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            !
            CALL SPLIT_CURRENT_MODEL(J-1,I)
            !
            CALL MPI_RECV(INV_MAS(1:9,J-1,I),9,MPI_DOUBLE_PRECISION &
                ,M+1,7,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            !
            ! Best equiv. response function:
            !
            CALL MPI_RECV(BEST_DSYN(:,:,J-1,I),NFREQ*NFREEP&
                ,MPI_REAL,M+1,8,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Best Stokes spectra:
            !
            CALL MPI_RECV(BEST_SYN(:,J-1,I),NFREQ&
                ,MPI_DOUBLE_PRECISION,M+1,9,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            !
            IF (INV_MAS(1,J-1,I).GT.MAXSTEPS) AM_I_DONE(J-1,I)=1
            !
            M=M+1
            !
          ENDDO
          ! FINALLY, MASTER CHECKS THAT WE HAVE FILLED THE STOKES ARRAY, DERIVATIVES AND SOME ATMOSPHERIC PARAMETERS:
        ELSE
          ! Slaves:
          IF (mpi__myrank.LE.max_work) THEN
            !
            ! Receive from master
            !
            ! Observed Stokes profiles:
            CALL MPI_RECV(ROBS1D(:),NFREQ,MPI_DOUBLE_PRECISION,0,52&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Current synthetic profiles:
            CALL MPI_RECV(RSYN1D(:),NFREQ,MPI_DOUBLE_PRECISION,0,53&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Best synthetic profiles:
            CALL MPI_RECV(BEST_RSYN1D(:),NFREQ,MPI_DOUBLE_PRECISION,0,54&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
!><            CALL MPI_RECV(INV_SLA(:),MSGSIZE,MPI_DOUBLE_PRECISION,0,1&
!><                ,MPI_COMM_WORLD,MPI__STATUS &
!><                ,mpi__ierror)
            CALL MPI_RECV(INV_SLA(1:9),9,MPI_DOUBLE_PRECISION,0,1&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            !
            ! Current equiv. response function:
            !
            CALL MPI_RECV(EQVDSYN(:,:),NFREQ*NFREEP&
                ,MPI_REAL,0,2,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            JACOB(:,:)=DBLE(EQVDSYN(:,:))
            !
            !
            ! Best equiv. response function:
            !
            CALL MPI_RECV(EQVDSYN(:,:),NFREQ*NFREEP&
                ,MPI_REAL,0,3,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Receive SIGMA:
            !
            CALL MPI_RECV(ISIGMAP(:),NFREQ,MPI_DOUBLE_PRECISION &
                ,0,4,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            CALL MPI_RECV(MODEL2D_RCV,2*NZ*ATM_ARGS &
                ,MPI_REAL,0,5,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Do my duties:
            !
            ! Solve LM in 1D
            CALL GET_PERT_1D()
            !PRINT*, 'Propuesta: ', MODEL2D_RCV(:,1,1)
            !PRINT*, 'Salvado: ', MODEL2D_RCV(:,1,2)
            !
            ! Send products back to master
            !
            CALL MPI_SEND(MODEL2D_RCV,2*NZ*ATM_ARGS,MPI_REAL,0,6,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Send Inversion auxiliars
            !
            CALL MPI_SEND(INV_SLA(1:9),9,MPI_DOUBLE_PRECISION,0,7&
                ,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            !
            ! Send best response functions:
            !
            CALL MPI_SEND(EQVDSYN(:,:),NFREQ*IFREEP&
                ,MPI_REAL,0,8,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            !
            ! Send best Stokes spectra:
            !
            CALL MPI_SEND(BEST_RSYN1D(:),NFREQ&
                ,MPI_DOUBLE_PRECISION,0,9,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
  IF (mpi__myrank .EQ. 0) THEN
!IF (MAXVAL(ABS(MODEL))) THEN
    !WRITE(*,*) SUM(DMODEL)/REAL(NX*NY*NZ), 'PERTURBATION AVERAGE'
    PRINT*, 'GET_DMODEL(END): ', NX*NY, SUM(AM_I_DONE)
  ENDIF
  !   
  END SUBROUTINE GET_DMODEL3DS
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_ERRORS3DS()
  !
    !
    INTEGER              :: I, J, L, M
    INTEGER              :: RJ, RI
    INTEGER              :: max_work
    !
    ! WORK DIVISION:
    ! INITIALISE X AND Y INDEXES
    I=1
    J=0
    ! LOOP THE THROUGH THE WHOLE DATA IN STEPS OF WORKING SLAVES (MASTER DOES NOT WORK)
    DO L=1,NX*NY,mpi__size-1
      !
      ! CHECK IF WE ARE IN THE LAST X-Y ITERATION (IF IT IS THE LAST, IT MIGHT BE THAT
      ! NOT ALL THE SLAVES HAS SOMETHING TO DO
      ! max_work IS THE VARIABLE THAT ENCODES THIS INFO
      max_work=mpi__size-1
      IF ( (NX*NY-L+1).LT.mpi__size) THEN
         max_work=NX*NY-L+1
      ENDIF
      !
      ! ONLY THE PROCS THAT HAS TO WORK IN THIS L (X-Y) ITERATION MUST DO SOMETHING:
      IF (mpi__myrank.LE.max_work) THEN
        !
        ! MASTER SENDS I,J ATMOSPHERE MODELS TO SLAVES AND WAIT FOR THE STOKES AND DERIVATIVES:
        IF (mpi__myrank .EQ. 0) THEN
          !
          ! LOOP THROUGH THE NUMBER OF SLAVES WORKING
          RJ=J
          RI=I
          DO M=1,max_work
            J=J+1
            ! Send inversion stuff:
            !
            ! Diff. between the observed and best synthetic spectra
            CALL SETUP_INVMAS(J,I)
            CALL MPI_SEND(INV_MAS(:,J,I),MSGSIZE,MPI_DOUBLE_PRECISION,M,1,MPI_COMM_WORLD &
                ,mpi__ierror)
            ! 
            !
            ! Best equiv. response function:
            !
!PRINT*, 'GET_ERROR3DS: ', I, J, SUM(ABS(BEST_DSYN(:,:,J,I)))
            CALL MPI_SEND(BEST_DSYN(:,:,J,I),NFREQ*IFREEP&
                ,MPI_REAL,M,2,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Finally, we also need to send the best atmosphere model:
            ! encode a I, J set of thermodynamic variables in a 3d array to be send
            CALL JOIN_BOTH_MODELS(J,I)
            ! Send atmospheric parameters:
            CALL MPI_SEND(MODEL2D_SND(:,:,2),NZ*ATM_ARGS,MPI_REAL,M,3,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            IF (MOD(J,NY).EQ.0) THEN
              I=I+1
              J=0
            ENDIF
          ENDDO
          ! RESET i AND j TO RECOVER SLAVE OUTPUT:
          J=RJ
          I=RI
          ! THIS WEIRD THING IS JUST NEEDED AS LONG AS WE AVOID ASYNCHRONICITY
          ! WHENEVER WE PROCEED WITH THAT PROPERTY, AN ASYNCHRONOUS SEND, AVOIDS
          ! BLOCKING COMMUNICATIONS AND SO, IT IS EASIER TO IMPLEMENT.
          ! YET, WE DO NOT USE THIS THING FOR SOME REASON I AM NOT CONVINCED OF.
          DO M=1,max_work
            J=J+1
            ! RECEIVE MODEL PERTURBATION:
            CALL MPI_RECV(MODEL1D_SND,NZ*ATM_ARGS,MPI_REAL,M,4 &
                ,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            !
            CALL ESPLIT_MODEL3D(I,J)
            !
            IF (MOD(J,NY).EQ.0) THEN
              I=I+1
              J=0
            ENDIF
          ENDDO
          ! FINALLY, MASTER CHECKS THAT WE HAVE FILLED THE STOKES ARRAY, DERIVATIVES AND SOME ATMOSPHERIC PARAMETERS:
        ELSE
          !SLAVES:
          ! RECEIVE MODEL
          IF (mpi__myrank.LE.max_work) THEN
            !
            ! Receive from master
            !
            CALL MPI_RECV(INV_SLA(:),MSGSIZE,MPI_DOUBLE_PRECISION,0,1&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Current equiv. response function:
            !
            CALL MPI_RECV(EQVDSYN(:,:),NFREQ*NFREEP&
                ,MPI_REAL,0,2,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            JACOB(:,:)=DBLE(EQVDSYN(:,:))
            !
            ! Best model up to now
            !
            CALL MPI_RECV(MODEL1D_RCV,NZ*ATM_ARGS,MPI_REAL,0,3,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Calculate 1D errors:
            !
            CALL GET_ERROR_1D(MODEL1D_RCV)
            !
            ! Send the products: perturbation of the model:
            !
!PRINT*, 'sERRORS: ', SUM(MODEL1D_RCV,1)
            CALL MPI_SEND(MODEL1D_RCV,NZ*ATM_ARGS,MPI_REAL,0,4,MPI_COMM_WORLD &
                ,mpi__ierror)
!PRINT*, 'sERRORS: ', SUM(MODEL1D_RCV,1)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
  !   
  END SUBROUTINE GET_ERRORS3DS
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_DMODEL3DC()
    !
    PRINT*, 'Moved to another place'
    !   
  END SUBROUTINE GET_DMODEL3DC
  !
  !------------------------------------------------
  !
  SUBROUTINE SORT_DELTA_TO_ATMOSPHERE(PTR_DELTA,FACTOR,IPAR,VOFFSET,JOFFSET &
      ,INF,SUP,NJ,COEFS,ABSOLUTE)
    !
    REAL(DP), INTENT(IN), DIMENSION(:), POINTER  :: PTR_DELTA
    REAL(DP), INTENT(IN)                         :: FACTOR, SUP, INF
    REAL(DP), INTENT(IN)                         :: ABSOLUTE
    REAL(SP), DIMENSION(:), INTENT(INOUT),POINTER       :: IPAR
    INTEGER, INTENT(IN)                          :: NJ
    REAL(DP), DIMENSION(:,:), INTENT(IN), POINTER  :: COEFS
    !
    INTEGER, INTENT(INOUT)                     :: VOFFSET, JOFFSET
    !
    REAL(SP), DIMENSION(NZ)                    :: OPAR
    REAL(DP), DIMENSION(NJ)                    :: ITDELTA
    REAL(DP), DIMENSION(NJ)                    :: PERTURBATION
    REAL(DP), DIMENSION(NJ)                    :: AVGPAR
    !
    INTEGER                                    :: K, J, COUNTER
    !REAL(DP)                                   :: PERTURBATION, MVAL, SVAL
    REAL(DP), DIMENSION(NZ)                    :: NZPERT
    LOGICAL                                    :: LGCPRNT
    !
!    PRINT*, '**SORT_DELTA_TO_ATMOSPHERE**'
    !
    COUNTER=1
    !
    NZPERT(:)=0
    !
!
! 2020:
    ITDELTA(:)=PTR_DELTA(JOFFSET+1:JOFFSET+NJ)
    ITDELTA(:)=ITDELTA(:)*FACTOR
!PRINT*, ITDELTA
! 2020.
    !
    !
    ! Either if an absolute or relative allowed perturbation criteria is...
    ! ... supplied, delta is not allowed to take larger than the larger...
    ! ... value of this criterium:
    CALL CONTRACT_VECTOR(NJ, IPAR, COEFS, AVGPAR)
    !
    DO J=1,NJ
      IF (ABSOLUTE.GT.0) THEN
        PERTURBATION(J)=MAXVAL((/ABSOLUTE,(ABS(MPERT*AVGPAR(J)))/))
      ELSE
        PERTURBATION(J)=ABS(MPERT*AVGPAR(J))
      ENDIF
      IF (AVGPAR(J).EQ.0) PERTURBATION(J)=ABSOLUTE
    ENDDO
    !
    !
    DO J=1,NJ
      IF (ABS(ITDELTA(J)).GT.PERTURBATION(J)) THEN
        ITDELTA(J)=ITDELTA(J)/ABS(ITDELTA(J))*PERTURBATION(J)
      ENDIF
    ENDDO
    PTR_DELTA(JOFFSET+1:JOFFSET+NJ)=ITDELTA(:)
    !
    !
      !DO J=1,NJ
      !  IF (AVGPAR(J)+ITDELTA(J).LT.INF) THEN
      !    ITDELTA(J) = INF - AVGPAR(J)
      !  ENDIF
      !  IF (AVGPAR(J)+ITDELTA(J).GT.SUP) THEN
      !    ITDELTA(J) = SUP - AVGPAR(J)
      !  ENDIF
      !ENDDO
    CALL EXPAND_VECTOR(NJ, AVGPAR+ITDELTA, COEFS, NZPERT)

    IF (INVERSIONTYPE.EQ.1) THEN
      OPAR(:)=REAL(NZPERT)
      DO K=1,NZ
        IF (OPAR(K).LT.INF) THEN
          OPAR(K)=INF
        ENDIF
        IF (OPAR(K).GT.SUP) THEN
          OPAR(K)=SUP
        ENDIF
      ENDDO
      IPAR(:)=OPAR(:)-IPAR(:)
    ENDIF
    !
    ! SIR/NICOLE like perts:
    IF (INVERSIONTYPE.EQ.0) THEN
      NZPERT = NZPERT - IPAR
      CALL CONTRACT_VECTOR(NJ, REAL(NZPERT), COEFS, ITDELTA)
      CALL EXPAND_VECTOR(NJ, ITDELTA, COEFS, NZPERT)
    PTR_DELTA(JOFFSET+1:JOFFSET+NJ)=ITDELTA(:)
!
      OPAR(:)=IPAR(:)+REAL(NZPERT)
      !
      DO K=1,NZ
        IF (OPAR(K).LT.INF) THEN
          NZPERT(K)=INF-IPAR(K)
        ENDIF
        IF (OPAR(K).GT.SUP) THEN
          NZPERT(K)=SUP-IPAR(K)
        ENDIF
      ENDDO
      !
      IPAR(:)=REAL(NZPERT(:))
    !ELSE
    ENDIF
    !! STiC like inversion:
    !  CALL EXPAND_VECTOR(NJ, ITDELTA+AVGPAR, COEFS, NZPERT)
    !  OPAR(:)=REAL(NZPERT)
    !  DO K=1,NZ
    !    IF (OPAR(K).LT.INF) THEN
    !      OPAR(K)=INF
    !    ENDIF
    !    IF (OPAR(K).GT.SUP) THEN
    !      OPAR(K)=SUP
    !    ENDIF
    !  ENDDO
    !  IPAR(:)=OPAR(:)
    !ENDIF
    !
    ! UPDATE VALUE
    VOFFSET=VOFFSET+1
    JOFFSET=JOFFSET+NJ
    !
  END SUBROUTINE SORT_DELTA_TO_ATMOSPHERE
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE REARANGE_DA_MODEL(PTR_MODEL2D, PTR_DELTA, PTR_FACTORS,ERRORS)
    !
    !USE COUPLED_PARAM, ONLY: RCV_DA1D, RCV_NFACT1D
    !
    REAL(SP), INTENT(INOUT), POINTER, DIMENSION(:,:) :: PTR_MODEL2D
    REAL(DP), INTENT(IN), POINTER, DIMENSION(:) :: PTR_DELTA
    REAL(DP), INTENT(IN), POINTER, DIMENSION(:) :: PTR_FACTORS
    LOGICAL, INTENT(IN), OPTIONAL               :: ERRORS
    !
    INTEGER                      :: VPAR, JPAR, I, CNT_IPARV, K
    !REAL(DP)                     :: ABSFACTOR
    REAL(SP), ALLOCATABLE, DIMENSION(:) :: BUMOD
    REAL(SP) :: TLLIM, TULIM
    !
    REAL(DP), DIMENSION(8)          :: LOWLIM
    REAL(DP), DIMENSION(8)          :: UPPLIM
    INTEGER, DIMENSION(8)           :: ATMPOS
    REAL(DP), DIMENSION(8)          :: ABSFACTOR
    !
    REAL(DP),DIMENSION(:,:), POINTER  :: PTR_COEF
    REAL(SP),DIMENSION(:), POINTER    :: PTR_MODEL
    !
!    PRINT*, '**REARANGE_DA_MODEL**'
    !
    !
    ATMPOS(1)=1
    ATMPOS(2)=2
    ATMPOS(3)=3
    ATMPOS(4)=6
    ATMPOS(5)=7
    ATMPOS(6)=8
    ATMPOS(7)=11
    ATMPOS(8)=2
    !
    IF (PRESENT(ERRORS)) THEN
      LOWLIM(:)=0.0D0
      !
      UPPLIM(:)=1.0D20
      !
      ABSFACTOR(:)=1.0D30
      !
    ELSE
      LOWLIM(1)=2500.0D0
      LOWLIM(2)=1.0D-2
      LOWLIM(3)=1.0D-10
      LOWLIM(4)=-5.0D+4
      LOWLIM(5)=-5.0D+4
      LOWLIM(6)=-5.0D+4
      LOWLIM(7)=-1.0D+7
      LOWLIM(8)=1.0D-8
      !
      UPPLIM(1)=1.8D+4
      UPPLIM(2)=1.0D+12
      UPPLIM(3)=1.0D+10
      UPPLIM(4)=5.0D+4
      UPPLIM(5)=5.0D+4
      UPPLIM(6)=5.0D+4
      UPPLIM(7)=1.0D+7
      UPPLIM(8)=1.0D12
      !
      IF (CURIC.EQ.1) THEN
        ABSFACTOR(1)=-100.!-1.0D4
        ABSFACTOR(2)=-1.0D6
        ABSFACTOR(3)=-1.0D2
        ABSFACTOR(4)=987.654D0
        ABSFACTOR(5)=987.654D0
        ABSFACTOR(6)=987.654D0
        ABSFACTOR(7)=98765.4321D0
        ABSFACTOR(8)=-1.0D0
      ELSE
        ABSFACTOR(1)=-100.0!-1.0D5
        ABSFACTOR(2)=-1.0D7
        ABSFACTOR(3)=-1.0D2
        ABSFACTOR(4)=98.7654D0
        ABSFACTOR(5)=98.7654D0
        ABSFACTOR(6)=98.7654D0
        ABSFACTOR(7)=9876.54321D0
        ABSFACTOR(8)=-1.0D0
      ENDIF
    ENDIF ! Errors
    !
    VPAR=0
    JPAR=0
    !
    CNT_IPARV=0
    !
    DO I=1,8
      IF (INV_ATMPAR(I).EQV..TRUE.) THEN
        !
        CNT_IPARV=CNT_IPARV+1
        !
        ! Assign pointer:
        SELECT CASE (I)
        CASE(1)
           PTR_COEF => TM_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,1)
        CASE(2)
           PTR_COEF => PG_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,2)
        CASE(3)
           PTR_COEF => RH_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,3)
        CASE(4)
           PTR_COEF => BX_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,6)
        CASE(5)
           PTR_COEF => BY_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,7)
        CASE(6)
           PTR_COEF => BZ_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,8)
        CASE(7)
           PTR_COEF => VZ_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,11)
        CASE(8)
           PTR_COEF => P0_COEFS(:,:)
           PTR_MODEL => PTR_MODEL2D(:,2)
        ENDSELECT  
        !
        ALLOCATE(BUMOD(SIZE(PTR_COEF,1)))
        BUMOD(:)=PTR_MODEL(:)
        !
        CALL SORT_DELTA_TO_ATMOSPHERE(PTR_DELTA,PTR_FACTORS(I) &
            ,PTR_MODEL,VPAR,JPAR,LOWLIM(I),UPPLIM(I) &
            , SIZE(PTR_COEF,2), PTR_COEF, PARAM_NORMS(I)*MPERT)!ABSFACTOR(I))
        !
        ! If I am dealing with T:
        IF (I.EQ.1) THEN
          !TLLIM=-1.5E0
          TLLIM=-4.5E0
          TLLIM=-11.5E0
          TULIM=0.5E0
          !TULIM=7.5E0
          DO K=SIZE(BUMOD),1,-1
            ! Atmosphere above sensitivity area:
            IF (PTR_MODEL2D(K,12).LT.TLLIM) THEN
              PTR_MODEL(K)=EXP(-(ABS(PTR_MODEL2D(K,12)-TLLIM)/0.5)**2)*PTR_MODEL(K)
            ELSE IF (PTR_MODEL2D(K,12).GT.TULIM) THEN
              PTR_MODEL(K)=PTR_MODEL(K+1)
            ENDIF
          ENDDO
        ENDIF 
        PTR_MODEL(:)=PTR_MODEL(:)+BUMOD(:)

IF (.FALSE.) THEN
        IF (INVERSIONTYPE.EQ.0) THEN
          !! If I am dealing with T:
          !IF (I.EQ.1) THEN
          !  TLLIM=-1.5E0
          !  TLLIM=-4.5E0
          !  TLLIM=-11.5E0
          !  TULIM=0.5E0
          !  DO K=SIZE(BUMOD),1,-1
          !    ! Atmosphere above sensitivity area:
          !    IF (PTR_MODEL2D(K,12).LT.TLLIM) THEN
          !      PTR_MODEL(K)=EXP(-(ABS(PTR_MODEL2D(K,12)-TLLIM)/0.5)**2)*PTR_MODEL(K)
          !    ELSE IF (PTR_MODEL2D(K,12).GT.TULIM) THEN
          !      PTR_MODEL(K)=PTR_MODEL(K+1)
          !    ENDIF
          !  ENDDO
          !ENDIF 
          PTR_MODEL(:)=PTR_MODEL(:)+BUMOD(:)
        ENDIF
ENDIF
        !
        DEALLOCATE(BUMOD)
        !
        NULLIFY(PTR_COEF)
        NULLIFY(PTR_MODEL)
        !
      ENDIF
    ENDDO
    !
  END SUBROUTINE REARANGE_DA_MODEL
  !
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_WEIGHTS()
    !
    USE MISC
    !
    INTEGER               :: IPAR, I
    !
    IF (AUTOWEIGHT.EQV..TRUE.) THEN
      CALL SET_WEIGHTS3D()
    ELSE
      IF (mpi__myrank.EQ.0) THEN 
        !
        IPAR=0
        ISIGMAP3D(:,:,:)=0.D0
PRINT*, 'Weights: ', WSTK, 'NUMW: ', NUMW
        DO I=1,4
          IF (INV_STK(I).EQV..TRUE.) THEN
            ISIGMAP3D(IPAR*NUMW+1:(IPAR+1)*NUMW,:,:)=WSTK(I)**2/ISIGMA**2
            IPAR=IPAR+1
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
  END SUBROUTINE SET_WEIGHTS
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_WEIGHTS3D()
    !
    INTEGER              :: I, J, L, M
    INTEGER              :: RJ, RI
    INTEGER              :: max_work
    INTEGER              :: NMAX_FWD
    !
    IF (mpi__myrank .EQ. 0) PRINT*, 'GET_AUTOWEIGHTS(START): ', NX*NY
    ! WORK DIVISION:
    ! INITIALISE X AND Y INDEXES
    I=1
    J=1
    ! LOOP THE THROUGH THE WHOLE DATA IN STEPS OF WORKING SLAVES (MASTER DOES NOT WORK)
    NMAX_FWD=NX*NY
    DO L=1,NMAX_FWD,mpi__size-1
      !
      ! CHECK IF WE ARE IN THE LAST X-Y ITERATION (IF IT IS THE LAST, IT MIGHT BE THAT
      ! NOT ALL THE SLAVES HAS SOMETHING TO DO
      ! max_work IS THE VARIABLE THAT ENCODES THIS INFO
      max_work=mpi__size-1
      IF ( (NMAX_FWD-L+1).LT.mpi__size) THEN
         max_work=NMAX_FWD-L+1
      ENDIF
      !
      ! ONLY THE PROCS THAT HAS TO WORK IN THIS L (X-Y) ITERATION MUST DO SOMETHING:
      IF (mpi__myrank.LE.max_work) THEN
        !
        ! MASTER SENDS I,J ATMOSPHERE MODELS TO SLAVES AND WAIT FOR THE STOKES AND DERIVATIVES:
        IF (mpi__myrank .EQ. 0) THEN
          !
          ! LOOP THROUGH THE NUMBER OF SLAVES WORKING
          RJ=J
          RI=I

          M=0
          DO WHILE (M.LT.max_work)
            IF (MOD(J,NY+1).EQ.0) THEN
              I=I+1
              J=1
            ENDIF
            J=J+1
            ! NOW, IF WE PASS THIS CLAUSE, WE MUST ACCESS J-1 IN THE FOLLOWING!
            ! M STARTS AT 0, NOT 1, SO ADD +1
            !
            ! Send Stokes profiles (observed)
            CALL SETUP_INVMAS_WEIGHT(J-1,I)
            CALL MPI_SEND(INV_MAS(:,J-1,I),MSGSIZE &
                ,MPI_DOUBLE_PRECISION,M+1,1,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! UPDATE M
            M=M+1
            !
          ENDDO
          ! RESET i AND j TO RECOVER SLAVE OUTPUT:
          J=RJ
          I=RI
          ! THIS WEIRD THING IS JUST NEEDED AS LONG AS WE AVOID ASYNCHRONICITY
          ! WHENEVER WE PROCEED WITH THAT PROPERTY, AN ASYNCHRONOUS SEND, AVOIDS
          ! BLOCKING COMMUNICATIONS AND SO, IT IS EASIER TO IMPLEMENT.
          ! YET, WE DO NOT USE THIS THING FOR SOME REASON I AM NOT CONVINCED OF.
          M=0
          DO WHILE (M.LT.max_work)
            IF (MOD(J,NY+1).EQ.0) THEN
              I=I+1
              J=1
            ENDIF
            J=J+1
            !
            ! NOW, IF WE PASS THIS CLAUSE, WE MUST ACCESS J-1 IN THE FOLLOWING!
            ! NOW, M+1 IS THE SENDER, NOT M (IT STARTS WITH 0, MASTER)
            !
            CALL MPI_RECV(INV_MAS(:,J-1,I),MSGSIZE,MPI_DOUBLE_PRECISION &
                ,M+1,2,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            ISIGMAP3D(:,J-1,I)=INV_MAS(10:MSGSIZE,J-1,I)
            !
            M=M+1
            !
          ENDDO
          ! FINALLY, MASTER CHECKS THAT WE HAVE FILLED THE STOKES ARRAY, DERIVATIVES AND SOME ATMOSPHERIC PARAMETERS:
        ELSE
          ! Slaves:
          IF (mpi__myrank.LE.max_work) THEN
            !
            ! Receive from master
            !
            ! It is sometimes useful to prevent Slaves from exhaustively ask for
            ! the message:
            !
            CALL MPI_RECV(INV_SLA(:),MSGSIZE,MPI_DOUBLE_PRECISION,0,1&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Do my duties:
            !
            ! Solve LM in 1D
            CALL SET_WEIGHTS1D()
            !
            ! Set isigmap:
            WSTK(:)=INV_SLA(5:8)
            CALL SET_ISIGMAP(NUMW)
            INV_SLA(10:MSGSIZE)=ISIGMAP
            !
            ! Send products back to master
            !
            CALL MPI_SEND(INV_SLA(:),MSGSIZE,MPI_DOUBLE_PRECISION,0,2&
                ,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
    IF (mpi__myrank .EQ. 0) THEN
      PRINT*, 'SET_WEIGHTS3D(END): ', NX*NY, SUM(AM_I_DONE)
    ENDIF
  !   
  END SUBROUTINE SET_WEIGHTS3D
  !
  !------------------------------------------------
  !
  SUBROUTINE SETUP_INVMAS_WEIGHT(JJ,II)
    !
    INTEGER, INTENT(IN)   :: II, JJ
    !
    INTEGER               :: OFFSET, I
    !
    OFFSET=10
    INV_MAS(OFFSET:MSGSIZE,JJ,II)=0.D0
    !
    ! Observed profiles to be sent:
    !
    INV_MAS(OFFSET:MSGSIZE,JJ,II)=OBS3D(:,JJ,II)
    !
  END SUBROUTINE SETUP_INVMAS_WEIGHT
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_WEIGHTS1D()
    !
    INTEGER    :: I
    INTEGER    :: OFFSET
    REAL(DP)   :: STD
    REAL(DP)   :: MINW
    !
    OFFSET=10
    MINW=1.D60
    DO I=1,SIZE(INV_STK)
      IF (INV_STK(I).EQV..TRUE.) THEN
        CALL STDDEV(NUMW,INV_SLA(OFFSET:OFFSET+NUMW-1),STD)
        INV_SLA(4+I)=DSQRT(1.D0/STD**2)
        OFFSET=OFFSET+NUMW
        IF (INV_SLA(4+I).LT.MINW) MINW=INV_SLA(4+I)
      ELSE
        INV_SLA(4+I)=0.D0
      ENDIF
    ENDDO
    !
    !PRINT*, INV_SLA(5:8)
    MINW=MINW*10.0D0
    DO I=1,SIZE(INV_STK)
      IF (INV_SLA(4+I).GT.MINW) INV_SLA(4+I)=MINW
    ENDDO
    !
    MINW=DSQRT(SUM(INV_SLA(5:8)**2))/DSQRT(DBLE(NSTKINV))
    !
    DO I=1,SIZE(INV_STK)
      INV_SLA(4+I)=INV_SLA(4+I)/MINW
    ENDDO
    PRINT*, DSQRT(SUM(INV_SLA(5:8)**2))/DSQRT(DBLE(NSTKINV))
    !
  END SUBROUTINE SET_WEIGHTS1D
  !
  !------------------------------------------------
  !
  SUBROUTINE STDDEV(NN, ARRAY, VAL)                   
    !                                                 
    INTEGER, INTENT(IN)                    :: NN      
    REAL(DP), INTENT(IN), DIMENSION(NN)    :: ARRAY   
    REAL(DP), INTENT(INOUT)                :: VAL     
    !                                                 
    REAL(DP)                               :: MN      
    !                                                 
    INTEGER                                :: I       
    !                                                 
    ! PRINT*, '**STDDEV**'                             
    !                                                 
    MN=SUM(ARRAY/DBLE(NN))                            
    !                                                 
    VAL=0.D0                                          
    DO I=1,NN                                         
      VAL=VAL+(ARRAY(I)-MN)**2                        
    ENDDO                                             
    !                                                 
    VAL=VAL/(DBLE(NN)-1.D0)                           
    !                                                 
    VAL=DSQRT(VAL)                                     
    !                                                                    
  END SUBROUTINE STDDEV                                                  
  !
  !================================================
  !
!!!!SUBROUTINE PERT_TEMP(NJ,NZ,ARR,PER,I,PSIZE)
!!!!
!!!!  INTEGER, INTENT(IN) :: NJ,NZ
!!!!  REAL(SP), DIMENSION(NZ), INTENT(INOUT) :: ARR
!!!!  REAL(DP), INTENT(IN) :: PER
!!!!  INTEGER, INTENT(IN) :: I
!!!!  REAL(DP), INTENT(INOUT) :: PSIZE
!!!!
!!!!  REAL(DP), DIMENSION(NJ) :: AVGPAR
!!!!  REAL(DP), DIMENSION(NZ) :: NZPERT
!!!!  INTEGER :: J
!!!!
!!!!!  PRINT*, ' In pert_temp: ', NJ, NZ, SHAPE(ARR), PER, I
!!!!  PRINT*, SUM(ARR)
!!!!
!!!!!  PRINT*, SHAPE(TM_COEFS)
!!!!  CALL CONTRACT_VECTOR(NJ, ARR, TM_COEFS, AVGPAR)
!!!!
!!!!!  PRINT*, AVGPAR
!!!!
!!!!  DO J=1,NJ
!!!!    IF (J.EQ.I) THEN
!!!!      AVGPAR(J)=PER*AVGPAR(J)
!!!!      PSIZE=AVGPAR(J)
!!!!    ELSE
!!!!      AVGPAR(J)=0.0D0
!!!!    ENDIF
!!!!  ENDDO
!!!!
!!!!  !PRINT*, AVGPAR
!!!!  CALL EXPAND_VECTOR(NJ, AVGPAR, TM_COEFS, NZPERT)
!!!!  !PRINT*, NZPERT
!!!!
!!!!  ARR(:)=ARR(:)+REAL(NZPERT(:))
!!!!
!!!!!  CALL CONTRACT_VECTOR(NJ, ARR, TM_COEFS, AVGPAR)
!!!!!
!!!!!  PRINT*, AVGPAR
!!!!
!!!!
!!!!END SUBROUTINE PERT_TEMP


END MODULE GET_DMODEL
!

