!
MODULE GET_DMODEL
  !
  !================================================
  !
  USE USER_MPI
  !
  USE CONS_PARAM, ONLY: SP, DP
  USE GRID_PARAM, ONLY: ZZ, NZ, NY, NX
  USE PHYS_PARAM
  USE FORWARD_PARAM, ONLY: ATM_ARGS, NUMW
  USE INVERT_PARAM, ONLY: INV_MAS, INV_SLA, NSLAB_PER_FREEV &
      , INV_ATMPAR, MAXSTEPS, INV_STK, NFREQ, NFREEP, AM_I_DONE &
      , MSGSIZE, MAXITER, ISIGMAP, WSTK, NSTKINV, ISIGMA, CURIC &
      , AUTOWEIGHT, NFREEV, YGUESS, YTOFIT, CYCPOW, NSLB_MAX &
      , JACOB, DELTA, ATM_FACTOR, MPERT, IFREEP, INU, SVDTOL &
      , NJEVALS, TOFFSET
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
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: GET_DMODEL3DC
  PUBLIC :: GET_DMODEL3DS
  PUBLIC :: IFREE_VAR_SPACE
  PUBLIC :: SET_WEIGHTS
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
  SUBROUTINE GET_GDMODEL()
    !
    ! 1) Update X2:
    CALL UPDATE_CHI2()
    !
    ! 2) Update Model atmosphere and response functions:
    CALL UPDATE_STOKES_MODEL_RF()
    !
    ! 3) :
    !
    ! 4) :
    !
  END SUBROUTINE GET_GDMODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_CHI2()
    !
    USE CODE_MODES, ONLY: COUPLED
    !
    ! This subroutine updates INV_MAS(2,:,:)...
    ! ... and INV_MAS(10:MSFSIZE,:,:)
    ! If coupled inversion:
    IF (COUPLED .EQV. .TRUE.) THEN
      CALL UPDATE_CHI2_COUPLED()
    ELSE ! Otherwise:
      CALL UPDATE_CHI2_NONCOUPLED()
    ENDIF
    !
    ! Check bad pixels and change them before proceeding?
    ! I might change 
    !
  END SUBROUTINE UPDATE_CHI2
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_STOKES_MODEL_RF()
    !
    USE CODE_MODES, ONLY: COUPLED
    !
    ! This subroutine updates INV_MAS(2,:,:)...
    ! ... and INV_MAS(10:MSFSIZE,:,:)
    ! If coupled inversion:
    IF (COUPLED .EQV. .TRUE.) THEN
      CALL UPDATE_SMRF_COUPLED()
    ELSE ! Otherwise:
      CALL UPDATE_SMRF_NONCOUPLED()
    ENDIF
    !
    !
    ! Update AM_I_DONE!!!
    ! AQI
    ! Check bad pixels and change them before proceeding?
    ! I might change 
    !
  END SUBROUTINE UPDATE_STOKES_MODEL_RF
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_SMRF_COUPLED()
    !
    USE COUPLED_PARAM, ONLY: COU_CHISQ,COU_OCHISQ &
        ,COU_SVDLAMBDA,COU_STEPS,COU_PRECYC
    !
    !LOGICAL, INTENT(INOUT)    :: IMP, FAC
    LOGICAL                   :: IMP, FAC
    !
    REAL(DP)                  :: COU_CUR, COU_PRE
    REAL(DP)                  :: COU_LAMBDA, COU_STEP
    !
    INTEGER                   :: I, J
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      COU_CUR=SUM(INV_MAS(2,:,:))/DBLE(SIZE(INV_MAS(2,:,:)))*1.0001D0
      COU_PRE=SUM(INV_MAS(3,:,:))/DBLE(SIZE(INV_MAS(3,:,:)))
      COU_LAMBDA=COU_SVDLAMBDA(1,1)
      COU_STEP=DBLE(COU_STEPS(1,1))/10.0D0
      !
      PRINT*, 'Update data'
      PRINT*, ' PrevX2: ', COU_PRE, ' ; CurX2: ', COU_CUR &
          , 'A: ', COU_LAMBDA, 'S: ', COU_STEP
      !
      IMP=.FALSE.
      FAC=.FALSE.
      !
      IMP=(COU_CUR.LT.COU_PRE)
      !
PRINT*, ' Before: ', COU_CUR, COU_PRE, COU_LAMBDA, COU_STEP, COU_PRECYC
      CALL UPDATE_LAMBDA(COU_CUR, COU_PRE, COU_LAMBDA, COU_STEP, COU_PRECYC, 0.01D0)
PRINT*, ' After: ', COU_CUR, COU_PRE, COU_LAMBDA, COU_STEP, COU_PRECYC
      !
      IF (IMP.EQV..TRUE.) THEN
        !
        ! Model:
        !
        BEST_TEM3D=TEM3D
        BEST_BX3D=BX3D
        BEST_BY3D=BY3D
        BEST_BZ3D=BZ3D
        BEST_VX3D=VX3D
        BEST_VY3D=VY3D
        BEST_VZ3D=VZ3D
        BEST_PG3D=PG3D
        BEST_RHO3D=RHO3D
        BEST_PEL3D=PEL3D
        BEST_MW3D=MW3D
        BEST_TAU3D5=TAU3D5
        !
        ! RFs:
        BEST_DSYN=DSYN
        !
        ! Stokes:
        BEST_SYN=SYN3D
        !
      ENDIF  ! Chi2 improvement
      !
      !COU_CHISQ(:,:)=COU_CUR
      INV_MAS(2,:,:)=COU_CUR
      !COU_OCHISQ(:,:)=COU_PRE
      INV_MAS(3,:,:)=COU_PRE
      !
      ! To change AQI
      COU_SVDLAMBDA(:,:)=COU_LAMBDA
      COU_STEPS(:,:)=NINT(COU_STEP*10.0D0)
      !
      PRINT*, ' PostX2: ', COU_PRE, ' Steps: ', COU_STEP&
          , ' Fac: ', FAC, ' IMP: ', IMP, 'A: ', COU_LAMBDA
            !
    ENDIF
    !
  END SUBROUTINE UPDATE_SMRF_COUPLED
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_SMRF_NONCOUPLED()
    !
    INTEGER              :: I, J, L, M
    INTEGER              :: RJ, RI
    INTEGER              :: max_work
    INTEGER              :: NMAX_FWD
    !
    IF (mpi__myrank .EQ. 0) PRINT*, 'GET_DMODEL(START): ', NX*NY, SUM(AM_I_DONE)
    ! Work division:
    ! Initialise x and y indexes
    I=1
    J=1
    ! Loop the through the whole data in steps of working slaves (master does not work)
    NMAX_FWD=NX*NY-SUM(AM_I_DONE)
    DO L=1,NMAX_FWD,mpi__size-1
      !
      ! Check if we are in the last x-y iteration (if it is the last, it might be that...
      ! ...not all the slaves has something to do...
      ! ...max_work IS THE VARIABLE THAT ENCODES THIS INFO.
      max_work=mpi__size-1
      IF ( (NMAX_FWD-L+1).LT.mpi__size) THEN
         max_work=NMAX_FWD-L+1
      ENDIF
      !
      ! Only the procs that has to work in this l (x-y) iteration must do something:
      IF (mpi__myrank.LE.max_work) THEN
        !
        ! Master sends i,j atmosphere models to slaves and wait for the stokes and derivatives:
        IF (mpi__myrank .EQ. 0) THEN
          !
          ! Loop through the number of slaves working
          RJ=J
          RI=I
          M=0
          DO WHILE (M.LT.max_work)
            IF (MOD(J,NY+1).EQ.0) THEN
              I=I+1
              J=1
            ENDIF
            J=J+1
            IF (AM_I_DONE(J-1,I).GT.0.5) CYCLE
            ! Now, if we pass this clause, we must access j-1 in the following!
            ! M starts at 0, not 1, so add +1
            !
            CALL MPI_SEND(INV_MAS(2:3,J-1,I),2 &
                ,MPI_DOUBLE_PRECISION,M+1,1,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
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
            ! Now, if we pass this clause, we must access j-1 in the following!
            ! Now, m+1 is the sender, not m (it starts with 0, master)
            !
            ! Receive model perturbation:
            CALL MPI_RECV(MODEL2D_SND,2*NZ*ATM_ARGS,MPI_REAL,M+1,5 &
                ,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            !
            CALL SPLIT_CURRENT_MODEL(J-1,I)
            !
            ! Best equiv. response function:
            !
            CALL MPI_RECV(BEST_DSYN(:,:,J-1,I),NFREQ*NFREEP&
                ,MPI_REAL,M+1,7,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            M=M+1
            !
          ENDDO
        ELSE ! Master. Slaves:
          ! Slaves:
          IF (mpi__myrank.LE.max_work) THEN
            !
            ! Receive from master
            !
            CALL MPI_RECV(INV_SLA(2:3),2,MPI_DOUBLE_PRECISION,0,1&
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
            !
            ! Best equiv. response function:
            !
            CALL MPI_RECV(EQVDSYN(:,:),NFREQ*NFREEP&
                ,MPI_REAL,0,3,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Receive both models:
            !
            CALL MPI_RECV(MODEL2D_RCV,2*NZ*ATM_ARGS &
                ,MPI_REAL,0,5,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Do my duties:
            CALL UPDATE_SMRF_1D()
            !
            ! Send products back to master
            !
            CALL MPI_SEND(MODEL2D_RCV,2*NZ*ATM_ARGS,MPI_REAL,0,5,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Send best response functions:
            !
            CALL MPI_SEND(EQVDSYN(:,:),NFREQ*IFREEP&
                ,MPI_REAL,0,7,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
    !
    IF (mpi__myrank .EQ. 0) THEN
      PRINT*, 'GET_DMODEL(END): ', NX*NY, SUM(AM_I_DONE)
    ENDIF
    !
  END SUBROUTINE UPDATE_SMRF_NONCOUPLED
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_SMRF_1D()
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
    ENDIF
    !
  END SUBROUTINE UPDATE_SMRF_1D
  





























  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_CHI2_COUPLED()
    !
    USE INVERT_PARAM, ONLY:NSTKINV,NSLAB_PER_FREEV,AM_I_DONE
    !
    USE USER_FFTW3, ONLY: FFTK2D, FFTY2D, GPLAN2D &
        , FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_3D_DP
    !
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: OWNCHI, OWNDIFF
    INTEGER :: K
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      CALL ALLOCATE_3D_DP(OWNCHI,NSTKINV*NUMW,NY,NX,'OWNCHI')
      CALL ALLOCATE_3D_DP(OWNDIFF,NSTKINV*NUMW,NY,NX,'OWNCHI')
      !
      DO K=1,NSTKINV*NUMW
        !
        ! DFT synthetic data:
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),SYN3D(K,:,:),FFTY2D)
        !
        ! Multiply transformed arrays:
        FFTY2D=FFTY2D*FFTK2D
        !
        ! DFT back the product:
        ! During testing, we overwrite OBS3D to compare between us and MvN
        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,OWNDIFF(K,:,:))
        !
        ! DFT synthetic data:
        OWNCHI(K,:,:)=DSQRT(ISIGMAP3D(K,:,:))
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),OWNCHI(K,:,:),FFTY2D)
        !
        ! Multiply transformed arrays:
        FFTY2D=FFTY2D*FFTK2D
        !
        ! DFT back the product:
        ! During testing, we overwrite OBS3D to compare between us and MvN
        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,OWNCHI(K,:,:))
      ENDDO
      !
      OWNDIFF=(OBS3D-OWNDIFF)*OWNCHI
      INV_MAS(10:MSGSIZE,:,:)=(OBS3D-OWNDIFF)*OWNCHI
      !
      ! X2:
      OWNCHI=OWNDIFF**2
      !COU_CHISQ=SUM(OWNCHI,DIM=1)
      INV_MAS(2,:,:)=SUM(OWNCHI,DIM=1)
      !
    ENDIF ! Master
    !
  END SUBROUTINE UPDATE_CHI2_COUPLED
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_CHI2_NONCOUPLED()
    !
    INTEGER              :: I, J, L, M
    INTEGER              :: RJ, RI
    INTEGER              :: max_work
    INTEGER              :: NMAX_FWD
    !
    IF (mpi__myrank .EQ. 0) PRINT*, 'UPDATE_CHI2_NOCOUPLED(START): ', NX*NY, SUM(AM_I_DONE)
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
            CALL SETUP_INVMAS(J-1,I)
            CALL MPI_SEND(INV_MAS(:,J-1,I),MSGSIZE &
                ,MPI_DOUBLE_PRECISION,M+1,1,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Send SIGMA:
            !
            CALL MPI_SEND(ISIGMAP3D(:,J-1,I),NFREQ,MPI_DOUBLE_PRECISION &
                ,M+1,4,MPI_COMM_WORLD &
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
            !
            CALL MPI_RECV(INV_MAS(:,J-1,I),MSGSIZE,MPI_DOUBLE_PRECISION &
                ,M+1,6,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            !
            IF (INV_MAS(1,J-1,I).GT.MAXSTEPS) AM_I_DONE(J-1,I)=1
            !
            M=M+1
            !
          ENDDO
          !
        ELSE
          ! Slaves:
          IF (mpi__myrank.LE.max_work) THEN
            !
            ! Receive from master
            !
            CALL MPI_RECV(INV_SLA(:),MSGSIZE,MPI_DOUBLE_PRECISION,0,1&
                ,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Receive SIGMA:
            !
            CALL MPI_RECV(ISIGMAP(:),NFREQ,MPI_DOUBLE_PRECISION &
                ,0,4,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! Do (O-S)/N
            ! Do SUM((O-S)**2/N**2)
            CALL GET_CHI2_1D()
            !
            ! Send result back to master
            !
            CALL MPI_SEND(INV_SLA(:),MSGSIZE,MPI_DOUBLE_PRECISION,0,6&
                ,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
    !
    IF (mpi__myrank .EQ. 0) THEN
      PRINT*, 'UPDATE_CHI2_NOCOUPLED(END): ', NX*NY, SUM(AM_I_DONE)
    ENDIF
    !
  END SUBROUTINE UPDATE_CHI2_NONCOUPLED
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_CHI2_1D()
    !
    INV_SLA(10:MSGSIZE)=INV_SLA(10:MSGSIZE) * DSQRT(ISIGMAP(:))
    INV_SLA(2)=SUM((INV_SLA(10:MSGSIZE))**2)
    !
  END SUBROUTINE GET_CHI2_1D
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
  SUBROUTINE JOIN_BEST_MODEL(J,I)
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
       MODEL1D_SND(:,K) = ARR (:,J,I)
    !
       NULLIFY(ARR)
    ENDDO
    MODEL1D_SND(:,13) = ZZ(:)
    !
  END SUBROUTINE JOIN_BEST_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE SPLIT_NEW_MODEL(J,I)
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
IF (MAXVAL(MODEL1D_SND(:,K)).NE.MAXVAL(MODEL1D_SND(:,K))) THEN
  PRINT*, 'Modified model NAN: ', MODEL1D_SND(:,K) &
      , ' ; K= ', K, ' ; J= ', J, ' ; I= ', I
ENDIF
       ARR (:,J,I) = MODEL1D_SND(:,K)

       NULLIFY(ARR)
    ENDDO
    ZZ(:)=MODEL1D_SND(:,13)
    !
  END SUBROUTINE SPLIT_NEW_MODEL
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
  SUBROUTINE IFREE_VAR_SPACE()
    !
    INTEGER                  :: I
    INTEGER                  :: OFFSET
    INTEGER                  :: CNT_FREEV
    INTEGER                  :: INSLABS
    !
    OFFSET=0
    !
    !MPERT=(90.D0-DBLE(CURIC)/DBLE(MAXITER)*8.0D1)*1.D-2
    !MPERT=10.D0*DEXP(DLOG(1.0D0/1.0D1)*(DBLE(CURIC)/DBLE(MAXITER)))*1.D-2
    !MPERT=50.D0*DEXP(DLOG(1.0D1/5.0D1)*(DBLE(CURIC)/DBLE(MAXITER)))*1.D-2
    MPERT=0.5D0
    SVDTOL=4
    IF (CURIC.EQ.1) THEN
      NJEVALS=1
    ELSE IF (CURIC.EQ.MAXITER) THEN
      NJEVALS=NZ
    ELSE
      NJEVALS=INT(CYCPOW**(CURIC-1))
!
!
IF (mpi__myrank.EQ.1) PRINT*, 'NJEVALS, CURIC, MAXITER, CYCPOW: ', NJEVALS, CURIC, MAXITER, CYCPOW
IF (NJEVALS.GT.NZ) NJEVALS=NZ
!
    ENDIF
IF (mpi__myrank.EQ.1) PRINT*, 'NJEVALS: ', NJEVALS
    !
    IFREEP=OFFSET
    CNT_FREEV=1
    DO I=1,SIZE(INV_ATMPAR)
      !
      ! Warning, testing the possibility of having T with more...
      ! ...freedom from the beginning than the others.

      IF (I.EQ.1) THEN
IF (mpi__myrank.EQ.1) PRINT*, 'I am ', I, ' in cycle: ', CURIC, ' and I am using ', NJEVALS, ' perturbation'
        CALL UPDATE_COEFS(I,INV_ATMPAR(I),NSLB_MAX(I),NJEVALS,NZ &
            ,INSLABS)
      ELSE
        IF (CURIC.GT.TOFFSET) THEN
          CALL UPDATE_COEFS(I,INV_ATMPAR(I),NSLB_MAX(I) &
               ,INT(CYCPOW**(CURIC-TOFFSET-1)),NZ &
              ,INSLABS)
        ELSE
IF (mpi__myrank.EQ.1) PRINT*, 'I am ', I, ' in cycle: ', CURIC &
    , ' and I am using just 1 perturbation to allow T more freedom' &
    , TOFFSET
          CALL UPDATE_COEFS(I,INV_ATMPAR(I),NSLB_MAX(I),1,NZ &
              ,INSLABS)
        ENDIF
      ENDIF
! Warning end.
      IF(INV_ATMPAR(I)) THEN
        NSLAB_PER_FREEV(CNT_FREEV)=INSLABS
        CNT_FREEV=CNT_FREEV+1
      ENDIF
    ENDDO
IF (mpi__myrank.EQ.1) PRINT*, NSLAB_PER_FREEV
    !
  END SUBROUTINE IFREE_VAR_SPACE
  !
  !------------------------------------------------
  !
  SUBROUTINE START_LM()
    !
    !
    !
  END SUBROUTINE START_LM
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
    !
    STEP=DBLE(INZ-1)/DBLE(INJ)
    !
    DO I=1,INJ
      IPOS(I)=(DBLE(I)-0.5D0)*STEP+1.D0
    ENDDO
    !
  END SUBROUTINE GET_POSITIONS
  !
  !------------------------------------------------
  !
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
      CALL GET_POSITIONS(INZ, INJ, RINTERVALS)
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
      IFREEP=IFREEP+INJEVALS
      NSLAB_PAR=INJEVALS
      !
      NULLIFY(PTR_COEF)
      !
    ENDIF
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
    REAL(DP),DIMENSION(IFREEP,IFREEP) :: PRODMAT
    !
    !
    ITVAR=0
    !
    ! Apply the square root of the noise/weights
    !
    DO ITPAR=1,IFREEP
      JACOB(:,ITPAR)=JACOB(:,ITPAR) * DSQRT(ISIGMAP)
    ENDDO
    !
    ! Calculate and apply the normalization factors:
    !
    JINIT=1
    DO ITPAR=1,SIZE(INV_ATMPAR)
      IF (INV_ATMPAR(ITPAR).EQV..TRUE.) THEN

        ITVAR=ITVAR+1
        JDELT=NSLAB_PER_FREEV(ITVAR)-1


! Method 2:
!        PRODMAT(:,:)=MATMUL(TRANSPOSE(JACOB),JACOB)
!        NORM=MAXVAL(ABS(PRODMAT(JINIT:JINIT+JDELT,JINIT:JINIT+JDELT)))
!        NORM=DSQRT(NORM)
! Method 1:
        NORM=SUM(ABS(JACOB(:,JINIT:JINIT+JDELT)))


        ATM_FACTOR(ITPAR)=1.D0/MAXVAL((/NORM,1.0D-18/))

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
  SUBROUTINE NEW_SOLVE_PERTURBATION()
    !
    USE LM, ONLY: GET_PERTURBATION
    !
    INTEGER       :: I, CNT
    REAL(DP)      :: FACTOR, AMINAN
    !
    FACTOR=1.0D0
    !
    CALL GET_PERTURBATION(IFREEP, NFREQ&
        , JACOB, INV_SLA(10:MSGSIZE) &
        , INV_SLA(9), DELTA, SVDTOL)
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
  SUBROUTINE SOLVE_PERTURBATION()
    !
    USE LM, ONLY: GET_PERTURBATION
    !
    INTEGER       :: I, CNT
    REAL(DP)      :: FACTOR, AMINAN
    !
    FACTOR=1.0D0
    !
!IF (mpi__myrank.EQ.1) PRINT*, 'NFREQ: ', NFREQ, 'INV(9): ', INV_SLA(9)
    CALL GET_PERTURBATION(IFREEP, NFREQ&
        , JACOB*FACTOR, INV_SLA(10:9+NFREQ) &
        * SQRT(ISIGMAP)*FACTOR, INV_SLA(9), DELTA, SVDTOL)
    !
    DELTA(:)=DELTA(:)/FACTOR/FACTOR
    !
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
  END SUBROUTINE SOLVE_PERTURBATION
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_CHI(ITCHI)
    !
    REAL(DP), INTENT(INOUT)   :: ITCHI
    !
    ITCHI=SUM((YTOFIT-YGUESS)**2*ISIGMAP)/INU
    !
  END SUBROUTINE GET_CHI
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
      CALL NEW_REARANGE_DA_MODEL(PTR_MODEL2D, PTR_ITDELTA, PTR_ATMFACTOR,.TRUE.)
    ELSE
      CALL NEW_REARANGE_DA_MODEL(PTR_MODEL2D, PTR_ITDELTA, PTR_ATMFACTOR)
    ENDIF
    !
    !
  END SUBROUTINE LIMIT_PERTURBATION
  !
  !------------------------------------------------
  !
  SUBROUTINE APPLY_PERTURBATION()
    !
    !USE PHYS_PARAM, ONLY: MODEL1D_RCV
    USE COUPLED_PARAM, ONLY: RCV_DA1D, RCV_NFACT1D
    !USE INVERT_PARAM, ONLY: INV_ATMPAR
    !
    ! TBUC LOGICAL, INTENT(IN), OPTIONAL           :: ERRORS
    !
    REAL(SP), POINTER, DIMENSION(:,:)       :: MODEL2D
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ITDELTA
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ATMFACTOR
    !
    REAL(DP), TARGET, DIMENSION(8)          :: FULL_FACTOR
    INTEGER                                 :: I, CNT
    !
    ! Coupled thing has a reduced version of the factors:
    CNT=0
    FULL_FACTOR(:)=0.D0
    DO I=1,SIZE(INV_ATMPAR)
      IF (INV_ATMPAR(I).EQV..TRUE.) THEN
        CNT=CNT+1
        FULL_FACTOR(I)=RCV_NFACT1D(CNT)
      ENDIF
    ENDDO
    !
    MODEL2D=>MODEL1D_RCV
    PTR_ITDELTA=>RCV_DA1D(:)
    PTR_ATMFACTOR=>FULL_FACTOR(:)!!!RCV_NFACT1D(:)
    !
      CALL NEW_REARANGE_DA_MODEL(MODEL2D, PTR_ITDELTA, PTR_ATMFACTOR)
    !
    NULLIFY(MODEL2D)
    NULLIFY(PTR_ITDELTA)
    NULLIFY(PTR_ATMFACTOR)
    !
  END SUBROUTINE APPLY_PERTURBATION
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_LAMBDA(CHICUR, CHIPRE, LAMBDA, STEPS, PRECYC, IFACTOR)                                                            
    !
    REAL(DP), INTENT(INOUT)     :: CHICUR, CHIPRE, LAMBDA, STEPS, PRECYC                                                              
    REAL(DP), INTENT(IN), OPTIONAL :: IFACTOR
    !
    REAL(DP)                       :: FACTOR, LUPDATE
    !
    FACTOR=5.D-2
    IF (PRESENT(IFACTOR)) FACTOR=DBLE(IFACTOR)
!
    IF ((CHIPRE-CHICUR).GT.0) THEN
      ! CHI IMPROVES:
      ! WE ONLY UPDATE LAMBDA IF THE PREVIOUS WAS AN IMPROVEMENT:                                                                  
      IF ((CHIPRE-CHICUR).GT.FACTOR*CHIPRE) THEN
        STEPS=MAXVAL((/0.D0, STEPS-1.D0/))
        CHIPRE=CHICUR
        LUPDATE=10.D0
      ELSE
        STEPS=STEPS+0.5D0
        LUPDATE=1.0D0/2.0D0
      ENDIF
      CHIPRE=CHICUR

      LAMBDA=LAMBDA/LUPDATE
    ELSE
      IF (LAMBDA.LT.1.D-4) LAMBDA=1.D-4
      IF (ABS(PRECYC-1.D0).LT.1.D-2) THEN
        ! CHI DOES NOT IMPROVE BUT IN THE PREVIOUS STEP, WAS IMPROVING!                                                               
        PRECYC=0.D0
! WE CONTINUE DECREASING LAMBDA TO CHECK THAT WE LOST THE THING
!        LAMBDA=LAMBDA/10.D0!*2.D-1
! WE COME BACK TO THE PREVIOUS LAMBDA TO CHECK THAT WE LOST THE THING                                                                 
        LAMBDA=LAMBDA*9.0D0!*2.D-1
      ELSE
          LAMBDA=LAMBDA*25.0D0
      ENDIF
      STEPS=STEPS+1.D0
    ENDIF
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
    INV_SLA(10:MSGSIZE)=(ROBS1D(:)-BEST_RSYN1D(:))*DSQRT(ISIGMAP)
    !
    !
  END SUBROUTINE NEW_CALCULATE_CHI2
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_PERT_1D()
    !
    REAL(SP), POINTER, DIMENSION(:,:)       :: MODEL2D
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ITDELTA
    REAL(DP), POINTER, DIMENSION(:)         :: PTR_ATMFACTOR
    !
    ! First, we estimate the chi2:
    CALL NEW_CALCULATE_CHI2()
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
    CALL START_LM()
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
    USE COUPLED_PARAM
    USE COUPLED_INVERSION
    !USE INVERT_PARAM, ONLY:NSTKINV,NSLAB_PER_FREEV,AM_I_DONE
    !
    USE USER_FFTW3, ONLY: FFTK2D, FFTY2D, GPLAN2D &
        , FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
    !
    USE MISC, ONLY: WRITE_PROFILES
    USE coupled_matrix_inversion, ONLY: BUILD_ALPHA &
        ,BUILD_BETA,BUILD_HINV,CALCULATE_DELTA
    USE TIME_PARAM
    !
    USE MISC
    !
    INTEGER     :: I, J, K, S
    INTEGER     :: FFST, DFFST
    !
    INTEGER     :: MFT
    INTEGER     :: CNT
    REAL(DP)    :: SVDEPS
    REAL(DP)    :: INVNORM
    !
    INTEGER     :: NDFIT
    INTEGER     :: NFREEPAR
    !
    INTEGER     :: OFF
    REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: OWNCHI, OWNDIFF
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: ABSAVG
    REAL(DP) :: LAMBDAP
    !
    INTEGER :: MSTOP, SSTOP
    !
    LOGICAL :: DOINTR
    !
    REAL(DP)    :: MVAL
    REAL(DP), DIMENSION(IFREEP, IFREEP) :: PRODMAT
    !
    IF (mpi__myrank.EQ.0) WRITE(*,*) 'I AM get_dmodel3d BUT I AM NOT IMPLEMENTED YET'
    !
    MFT=IFREEP
    !
IF (mpi__myrank.EQ.0) PRINT*, 'INIT_EVERY_CYC_EVERY_IT'
    CALL INIT_EVERY_CYC_EVERY_IT()

    !Now, only master has duties:
    IF (mpi__myrank.EQ.0) THEN
      !
      ! Our alternative implementation:
      !
      !
      ! 1- Build alpha: J^T,J
      !
      CALL ALLOCATE_3D_DP(OWNCHI,NSTKINV*NUMW,NY,NX,'OWNCHI')
      CALL ALLOCATE_3D_DP(OWNDIFF,NSTKINV*NUMW,NY,NX,'OWNCHI')
      !
      DO K=1,NSTKINV*NUMW
        !
        ! DFT synthetic data:
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),SYN3D(K,:,:),FFTY2D)
        !
        ! Multiply transformed arrays:
        FFTY2D=FFTY2D*FFTK2D
        !
        ! DFT back the product:
        ! During testing, we overwrite OBS3D to compare between us and MvN
        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,OWNDIFF(K,:,:))
        !
        ! DFT synthetic data:
        OWNCHI(K,:,:)=DSQRT(ISIGMAP3D(K,:,:))
        !OWNCHI(K,:,:)=ISIGMAP3D(K,:,:)
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),OWNCHI(K,:,:),FFTY2D)
        !
        ! Multiply transformed arrays:
        FFTY2D=FFTY2D*FFTK2D
        !
        ! DFT back the product:
        ! During testing, we overwrite OBS3D to compare between us and MvN
        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,OWNCHI(K,:,:))
      ENDDO
      OWNDIFF=(OBS3D-OWNDIFF)*OWNCHI
      OWNCHI=OWNDIFF**2
      !
      !
      !
      ! We first check if we have improved the fitting:
      COU_CHISQ=SUM(OWNCHI,DIM=1)
      CALL UPDATE_STOKES_MODEL()
!
! Test!!! Beta uses best synthetic profiles, not the last one:
!
      DO K=1,NSTKINV*NUMW
        !
        ! DFT synthetic data:
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),BEST_SYN(K,:,:),FFTY2D)
        !
        ! Multiply transformed arrays:
        FFTY2D=FFTY2D*FFTK2D
        !
        ! DFT back the product:
        ! During testing, we overwrite OBS3D to compare between us and MvN
        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,OWNDIFF(K,:,:))
        !
        ! DFT synthetic data:
        OWNCHI(K,:,:)=DSQRT(ISIGMAP3D(K,:,:))
        !OWNCHI(K,:,:)=ISIGMAP3D(K,:,:)
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),OWNCHI(K,:,:),FFTY2D)
        !
        ! Multiply transformed arrays:
        FFTY2D=FFTY2D*FFTK2D
        !
        ! DFT back the product:
        ! During testing, we overwrite OBS3D to compare between us and MvN
        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,OWNCHI(K,:,:))
      ENDDO
      OWNDIFF=(OBS3D-OWNDIFF)*OWNCHI
      !
! Test!!! Beta uses best synthetic profiles, not the last one.
      ! And we build alpha:
      !
      ! Get alpha scaling factors:
      !CALL GET_COU_FACTORS()
      COU_ALPHA(:,:,:,:)=BEST_DSYN(:,:,:,:)
      DO K=1,IFREEP
        COU_ALPHA(:,K,:,:)=COU_ALPHA(:,K,:,:)*DSQRT(ISIGMAP3D(:,:,:))
      ENDDO
      PRINT*, SHAPE(PRODMAT), ' -> ' &
          , SHAPE(MATMUL(TRANSPOSE(COU_ALPHA(:,:,1,1)), COU_ALPHA(:,:,1,1)))
      DO I=1,NX
        DO J=1,NY
          FFST=1
          PRODMAT(:,:)=MATMUL(TRANSPOSE(COU_ALPHA(:,:,J,I)), COU_ALPHA(:,:,J,I))
          DO K=1,NFREEV
            !
            DFFST=NSLAB_PER_FREEV(K)-1
            PRODMAT(1:1,1:1)=DSQRT(ABS(PRODMAT(FFST:FFST+DFFST,FFST:FFST+DFFST)))
            INVNORM=PRODMAT(1,1)
            !INVNORM=SUM(ABS(COU_ALPHA(:,FFST:FFST+DFFST,J,I)))!/DBLE(NSTKINV*NUMW)*1.0D0
            !INVNORM=MAXVAL(DABS(COU_ALPHA(:,FFST:FFST+DFFST,J,I)))!/DBLE(NSTKINV*NUMW)
            INVNORM=DMAX1(1.D-10,INVNORM)
IF ((I.EQ.3).AND.(J.EQ.4)) PRINT*, 'FFST: ', FFST, 'DFFST: ', DFFST &
    , 'NSLAB_PER_FREEV(K): ', NSLAB_PER_FREEV(K), 'K: ', K &
    , 'FFST+(DFFST+1): ', FFST+(DFFST+1), INVNORM &
    , SHAPE(COU_ALPHA(:,FFST:FFST+DFFST,J,I)), SHAPE(COU_ALPHA)
            COU_NFACT(K,J,I)=1.0D0/INVNORM!SUM(ABS(COU_ALPHA(:,FFST:FFST+DFFST,J,I)))
            IF (COU_NFACT(K,J,I).NE.COU_NFACT(K,J,I)) COU_NFACT(K,J,I)=0.D0
            COU_ALPHA(:,FFST:FFST+DFFST,J,I)=COU_ALPHA(:,FFST:FFST+DFFST,J,I) &
                *COU_NFACT(K,J,I)
            FFST=FFST+(DFFST+1)
          ENDDO
        ENDDO
      ENDDO
!STOP

!
! For testing purposes only:
CALL WRITE_BIN('cou_alpha.bin', SIZE(SHAPE(COU_ALPHA)), SHAPE(COU_ALPHA) &
    ,SIZE(COU_ALPHA,KIND=8),REAL(COU_ALPHA), 3000, 3)
CALL WRITE_BIN('cou_beta.bin', SIZE(SHAPE(OWNDIFF)), SHAPE(OWNDIFF) &
    ,SIZE(OWNCHI,KIND=8),REAL(OWNCHI), 3000, 3)
! End purposes only.

      !
      !
      DOINTR=.TRUE.
      DO WHILE (DOINTR.EQV..TRUE.)
      CALL BUILD_ALPHA(SHAPE(COU_ALPHA),COU_ALPHA)
      !
      ! 2- Build beta: J^T,D
      CALL BUILD_BETA(SHAPE(COU_ALPHA),COU_ALPHA,SHAPE(OWNDIFF),OWNDIFF)
      !
      !
      !
      ! 3- Build alpha inverse:
      SVDEPS=10.0D0**(-4)!(-NINT(INV_MAS(6,1,1)))
      LAMBDAP=SUM(COU_SVDLAMBDA)/DBLE(SIZE(COU_SVDLAMBDA))
      CALL BUILD_HINV(SHAPE(COU_PSF),COU_PSF,(/COU_NPY/2+1,COU_NPX/2+1/) &
          ,COU_BLCKSZ,SVDEPS,LAMBDAP)
      !
      !
      !
      ! 4- Calculate delta (iterative):
        CALL CALCULATE_DELTA(COU_DA,DOINTR)
        IF (DOINTR.EQV..TRUE.) THEN
          COU_SVDLAMBDA=COU_SVDLAMBDA*10.0D0
          COU_STEPS(:,:)=COU_STEPS(:,:)+2
          IF (COU_STEPS(1,1).GT.50) DOINTR=.FALSE.
        ENDIF
      ENDDO
!><      DO K=1,IFREEP
!><        MVAL=SUM(COU_DA(K,:,:))/DBLE(NY)/DBLE(NX)
!><        !COU_DA(K,:,:)=MVAL+(COU_DA(K,:,:)-MVAL)/10.0D0
!><        !
!><        ! DFT synthetic data:
!><        OWNCHI(K,:,:)=COU_DA(K,:,:)-MVAL
!><        CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),OWNCHI(K,:,:),FFTY2D)
!><        !
!><        ! Multiply transformed arrays:
!><        FFTY2D=FFTY2D*FFTK2D
!><!        FFTY2D(NY/2-5:,:)=0.0D0
!><!        FFTY2D(:,NX/2+1-5:NX/2+1+5)=0.0D0
!><        !
!><        ! DFT back the product:
!><        ! During testing, we overwrite OBS3D to compare between us and MvN
!><        CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,COU_DA(K,:,:))
!><        COU_DA(K,:,:)=(COU_DA(K,:,:)*4.0D0/5.0D0+OWNCHI(K,:,:)*1.0D0/5.0D0 &
!><            )+MVAL
!><        !
!><      ENDDO
      IF (ALLOCATED(OWNCHI)) DEALLOCATE(OWNCHI)
      IF (ALLOCATED(OWNDIFF)) DEALLOCATE(OWNDIFF)
      !
      !
    ENDIF ! Master stuff


IF (mpi__myrank.EQ.0) THEN
PRINT*, 'Da(NAN): ', MINVAL(COU_DA(:,:,:)), MAXVAL(COU_DA(:,:,:))
PRINT*, 'Is this right (2)?: ', SUM(SUM(ABS(COU_DA),3),2), SHAPE(COU_DA)
  MSTOP=1
  !
  ! Avoid slaves consuming all the CPU:
  DO I=1,mpi__size-1
    CALL MPI_SEND(MSTOP,1,MPI_INTEGER,I,101,MPI_COMM_WORLD &
        ,mpi__ierror)
  ENDDO
ELSE !Master. Slaves:
  CALL NICE_WAITING(0,101)
  CALL MPI_RECV(SSTOP,1,MPI_INTEGER,0,101,MPI_COMM_WORLD &
      ,mpi__status, mpi__ierror)
  
ENDIF !Slaves.
    !
    ! step5: Transform from michiel's output to this code parameters:
    ! All the pixels are modified because we want to check if this...
    ! ... new estimation does improve the fit
    CALL FROM_DA_TO_MODEL()
    !
    ! Finally, check how many iterations have not suceeded. If above 10...
    ! ... continue to the next cycle:
    IF (mpi__myrank.EQ.0) THEN
      IF (REAL(COU_STEPS(1,1))/10.0D0.GT.MAXSTEPS) AM_I_DONE(:,:)=1
      PRINT*, 'Steps left: ', NX*NY-SUM(AM_I_DONE) &
          , ' S=',REAL(COU_STEPS(1,1))/10.
    ENDIF
    CALL MPI_BCAST(AM_I_DONE,NX*NY,MPI_INTEGER,0,MPI_COMM_WORLD &
                ,mpi__ierror)
IF (mpi__myrank.EQ.0) THEN
  PRINT*, 'COU_SVDLAMBDA: ', MINVAL(COU_SVDLAMBDA), MAXVAL(COU_SVDLAMBDA)
  PRINT*, 'COU_STEPS: ', MINVAL(COU_STEPS), MAXVAL(COU_STEPS)
ENDIF
  !   
  END SUBROUTINE GET_DMODEL3DC
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_STOKES_MODEL()
    !
    USE COUPLED_PARAM, ONLY: COU_CHISQ,COU_OCHISQ &
        ,COU_SVDLAMBDA,COU_STEPS,COU_PRECYC
    !
    !LOGICAL, INTENT(INOUT)    :: IMP, FAC
    LOGICAL                   :: IMP, FAC
    !
    REAL(DP)                  :: COU_CUR, COU_PRE
    REAL(DP)                  :: COU_LAMBDA, COU_STEP
    !
    INTEGER                   :: I, J
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      COU_CUR=SUM(COU_CHISQ)/DBLE(SIZE(COU_CHISQ))*1.0001D0
      COU_PRE=SUM(COU_OCHISQ)/DBLE(SIZE(COU_CHISQ))
      COU_LAMBDA=COU_SVDLAMBDA(1,1)
      COU_STEP=DBLE(COU_STEPS(1,1))/10.0D0
      !
      PRINT*, 'Update data'
      PRINT*, ' PrevX2: ', COU_PRE, ' ; CurX2: ', COU_CUR &
          , 'A: ', COU_LAMBDA, 'S: ', COU_STEP
      !
      IMP=.FALSE.
      FAC=.FALSE.
      !
      IMP=(COU_CUR.LT.COU_PRE)
      !
PRINT*, ' Before: ', COU_CUR, COU_PRE, COU_LAMBDA, COU_STEP, COU_PRECYC
      CALL UPDATE_LAMBDA(COU_CUR, COU_PRE, COU_LAMBDA, COU_STEP, COU_PRECYC, 0.01D0)
PRINT*, ' After: ', COU_CUR, COU_PRE, COU_LAMBDA, COU_STEP, COU_PRECYC
      !
      IF (IMP.EQV..TRUE.) THEN
        !
        ! Model:
PRINT*, 'Updating!'
        BEST_TEM3D=TEM3D
        BEST_BX3D=BX3D
        BEST_BY3D=BY3D
        BEST_BZ3D=BZ3D
        BEST_VX3D=VX3D
        BEST_VY3D=VY3D
        BEST_VZ3D=VZ3D
        BEST_PG3D=PG3D
        BEST_RHO3D=RHO3D
        BEST_PEL3D=PEL3D
        BEST_MW3D=MW3D
        BEST_TAU3D5=TAU3D5
        !
        ! Stokes:
        BEST_DSYN=DSYN
        ! Stokes:
        BEST_SYN=SYN3D
        !
        ! This is done in UPDATE_LAMBDA! COU_OCHISQ=COU_CHISQ
      ENDIF  ! Chi2 improvement
      !
      COU_CHISQ(:,:)=COU_CUR
      COU_OCHISQ(:,:)=COU_PRE
      COU_SVDLAMBDA(:,:)=COU_LAMBDA
      COU_STEPS(:,:)=NINT(COU_STEP*10.0D0)
      PRINT*, ' PostX2: ', COU_PRE, ' Steps: ', COU_STEP&
          , ' Fac: ', FAC, ' IMP: ', IMP, 'A: ', COU_LAMBDA
            !
    ENDIF
    !
  END SUBROUTINE UPDATE_STOKES_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE FROM_DA_TO_MODEL()
    !
    USE COUPLED_PARAM, ONLY: COU_DA, COU_NFACT, RCV_DA1D, RCV_NFACT1D
    !
    INTEGER              :: I, J, L, M
    INTEGER              :: RJ, RI
    INTEGER              :: max_work
    INTEGER              :: NMAX_FWD
    !
!    LOGICAL        :: READY
    !
    IF (mpi__myrank .EQ. 0) PRINT*, 'FROM_DA_TO_MODEL(START): '
    !
    !
    ! WORK DIVISION:
    ! INITIALISE X AND Y INDEXES
    I=1
    J=1
    ! Loop the through the whole data in steps of working slaves (master does not work)
    NMAX_FWD=NX*NY!!!-SUM(AM_I_DONE)
    !
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
        ! Master sends I,J atmosphere models to slaves and wait for the stokes and derivatives:
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
            ! Now, we must access j-1 in the following!
            ! m starts at 0, not 1, so add +1
            !
            ! Send inversion stuff
            ! 
            ! We need to send the current atmosphere model:
            ! Encode a i, j set of thermodynamic variables in a 2d array to be send
            CALL JOIN_BEST_MODEL(J-1,I)
            ! Send atmospheric parameters:
            CALL MPI_SEND(MODEL1D_SND,NZ*ATM_ARGS,MPI_REAL,M+1,3,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! We also send the perturbation calculated by Michiel's routines:
            ! Its dimension is the number of free parameters, i.e. the number...
            ! ... of parameters inverted times the number of slabs considered.
            CALL MPI_SEND(COU_DA(:,J-1,I),IFREEP,MPI_DOUBLE_PRECISION,M+1,4,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
            ! Finally, we need the jacobian normalization factors used:
            ! Its dimension is the number of physical parameters inverted
            CALL MPI_SEND(COU_NFACT(:,J-1,I),NFREEV,MPI_DOUBLE_PRECISION &
                ,M+1,5,MPI_COMM_WORLD &
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
            ! Now, we must access j-1 in the following!
            ! Now, m+1 is the sender, not m (it starts with 0, master).
            ! Receive new model:
            CALL MPI_RECV(MODEL1D_SND,NZ*ATM_ARGS,MPI_REAL,M+1,8 &
                ,MPI_COMM_WORLD, MPI__STATUS &
                ,mpi__ierror)
            !
            ! Fill 3D model accordingly:
            CALL SPLIT_NEW_MODEL(J-1,I)
            !
            M=M+1
            !
          ENDDO
          !
        ELSE
          !SLAVES:
          IF (mpi__myrank.LE.max_work) THEN
            !
            ! Avoid slaves from waiting in agressive mode:
            IF (L.EQ.1) THEN
              CALL NICE_WAITING(0,3)
            ENDIF
            !
            ! Receive model
            CALL MPI_RECV(MODEL1D_RCV,NZ*ATM_ARGS,MPI_REAL,0,3,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            ! Receive da
            CALL MPI_RECV(RCV_DA1D,IFREEP,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            ! Receive factor
            CALL MPI_RECV(RCV_NFACT1D,NFREEV,MPI_DOUBLE_PRECISION &
                ,0,5,MPI_COMM_WORLD,MPI__STATUS &
                ,mpi__ierror)
            !
            ! 
            ! Arrange Da to M
            CALL APPLY_PERTURBATION()
            !
            ! Send model to master:
            CALL MPI_SEND(MODEL1D_RCV,NZ*ATM_ARGS,MPI_REAL,0,8,MPI_COMM_WORLD &
                ,mpi__ierror)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
    !
    IF (mpi__myrank .EQ. 0) PRINT*, 'FROM_DA_TO_MODEL(end). '
    !
  END SUBROUTINE FROM_DA_TO_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE NEW_SORT_DELTA_TO_ATMOSPHERE(PTR_DELTA,FACTOR,IPAR,VOFFSET,JOFFSET &
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
    ITDELTA(:)=PTR_DELTA(JOFFSET+1:JOFFSET+NJ)*FACTOR
    !
    !
    ! Either if an absolute or relative allowed perturbation criteria is...
    ! ... supplied, delta is not allowed to take larger than the larger...
    ! ... value of this criterium:
    CALL CONTRACT_VECTOR(NJ, IPAR, COEFS, AVGPAR)
    !
    IF (ABSOLUTE.GT.0) THEN
      DO J=1,NJ
        PERTURBATION(J)=MAXVAL((/ABSOLUTE,(ABS(MPERT*AVGPAR(J)))/))
      ENDDO
    ELSE
      DO J=1,NJ
        PERTURBATION(J)=ABS(MPERT*AVGPAR(J))
      ENDDO
    ENDIF
    !
    !
    DO J=1,NJ
      IF (ABS(ITDELTA(J)).GT.PERTURBATION(J)) THEN
        ITDELTA(J)=ITDELTA(J)/ABS(ITDELTA(J))*PERTURBATION(J)
      ENDIF
    ENDDO
    !
    !
    !
    CALL EXPAND_VECTOR(NJ, ITDELTA, COEFS, NZPERT)
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
    IPAR(:)=IPAR(:)+REAL(NZPERT(:))
    !
    ! UPDATE VALUE
    VOFFSET=VOFFSET+1
    JOFFSET=JOFFSET+NJ
    !
  END SUBROUTINE NEW_SORT_DELTA_TO_ATMOSPHERE
  !
  !------------------------------------------------
  !
  SUBROUTINE NEW_REARANGE_DA_MODEL(PTR_MODEL2D, PTR_DELTA, PTR_FACTORS,ERRORS)
    !
    !USE COUPLED_PARAM, ONLY: RCV_DA1D, RCV_NFACT1D
    !
    REAL(SP), INTENT(INOUT), POINTER, DIMENSION(:,:) :: PTR_MODEL2D
    REAL(DP), INTENT(IN), POINTER, DIMENSION(:) :: PTR_DELTA
    REAL(DP), INTENT(IN), POINTER, DIMENSION(:) :: PTR_FACTORS
    LOGICAL, INTENT(IN), OPTIONAL               :: ERRORS
    !
    INTEGER                      :: VPAR, JPAR, I, CNT_IPARV
    !REAL(DP)                     :: ABSFACTOR
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
!><      PRINT*, '+++ERRORS+++'
!><      PRINT*, SUM(PTR_MODEL2D,1)
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
      LOWLIM(4)=-5.0D+5
      LOWLIM(5)=-5.0D+5
      LOWLIM(6)=-5.0D+5
      LOWLIM(7)=-1.0D+9
      LOWLIM(8)=1.0D-8
      !
      UPPLIM(1)=1.8D+4
      UPPLIM(2)=1.0D+12
      UPPLIM(3)=1.0D+10
      UPPLIM(4)=5.0D+5
      UPPLIM(5)=5.0D+5
      UPPLIM(6)=5.0D+5
      UPPLIM(7)=1.0D+9
      UPPLIM(8)=1.0D12
      !
      IF (CURIC.EQ.1) THEN
        ABSFACTOR(1)=-1.0D4
        ABSFACTOR(2)=-1.0D6
        ABSFACTOR(3)=-1.0D2
        ABSFACTOR(4)=1.0D4
        ABSFACTOR(5)=1.0D4
        ABSFACTOR(6)=1.0D4
        ABSFACTOR(7)=1.0D6
        ABSFACTOR(8)=-1.0D0
      ELSE
        ABSFACTOR(1)=-1.0D5
        ABSFACTOR(2)=-1.0D7
        ABSFACTOR(3)=-1.0D2
        ABSFACTOR(4)=-1.0D2
        ABSFACTOR(5)=-1.0D2
        ABSFACTOR(6)=-1.0D2
        ABSFACTOR(7)=1.0D4!1.0D4
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
        CALL NEW_SORT_DELTA_TO_ATMOSPHERE(PTR_DELTA,PTR_FACTORS(I) &
            ,PTR_MODEL,VPAR,JPAR,LOWLIM(I),UPPLIM(I) &
            , SIZE(PTR_COEF,2), PTR_COEF, ABSFACTOR(I))
        !
        NULLIFY(PTR_COEF)
        NULLIFY(PTR_MODEL)
        !
      ENDIF
    ENDDO
    !
  END SUBROUTINE NEW_REARANGE_DA_MODEL
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
END MODULE GET_DMODEL
!

