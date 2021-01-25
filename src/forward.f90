!
MODULE FORWARD
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP
  USE TIME_PARAM
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER              :: SIZE1,SIZE2
  INTEGER, DIMENSION(:,:), ALLOCATABLE, ASYNCHRONOUS :: BUFFER
  integer :: nbyt1, nbyt2, bpos1, bpos2
  !
  PUBLIC :: FORWARD3D
  PUBLIC :: UPDATE_BEST_TO_CURRENT
  PUBLIC :: UPDATE_Current_to_best
  !
  PRIVATE :: JOIN_MODEL2
  PRIVATE :: DSYN_PROPAGATION
  PRIVATE :: CALCULATE_EQUIVALENT_DSYN
  PRIVATE :: CALCULATE_EQUIVALENT_SYN1D
  PRIVATE :: FORWARD1D
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! forward3d
  ! update_best_to_current
  ! join_model2
  ! dsyn_propagation
  ! calculate_equivalent_dsyn
  ! calculate_equivalent_syn1d
  ! forward1d
  ! split_model3d
  !
  !------------------------------------------------
  !
  SUBROUTINE FORWARD3D()
    !
    USE user_mpi
    !
    USE GRID_PARAM, ONLY: NX, NY, NZ
    USE PHYS_PARAM, ONLY: MODEL1D_SND, SYN3D, DSYN, SYN5D &
        , MODEL1D_RCV, RSYN1D, EQVDSYN, FS_SYN1D, TAU3DLIN
    USE INVERT_PARAM, ONLY: NFREQ, AM_I_DONE, IFREEP
    USE FORWARD_PARAM, ONLY: ATM_ARGS, N_FWD_MODELLING, FULL_STOKES, NUMW &
        , TAULIN
    USE CODE_MODES, ONLY: MINVERSION, MRESPFUNCT, MSYNTHESIS, MTAULIN
    !
    INTEGER              :: I, J, L, M, RJ,RI, max_work
    INTEGER              :: IT_INDX
    INTEGER              :: NMAX_FWD



    nbyt1 = NZ*ATM_ARGS*4
    nbyt2 = NZ*ATM_ARGS*4
    IF ((MINVERSION.EQV..TRUE.).OR.(MSYNTHESIS.EQV..TRUE.)) THEN
      nbyt2 = nbyt2 + NFREQ * 8
    ENDIF
    IF (MRESPFUNCT.EQV..TRUE.) THEN
       nbyt2 = nbyt2 + NFREQ*IFREEP*4
    ENDIF

    size1=nbyt1/4
    size2=nbyt2/4

    IF (mpi__myrank .EQ. 0) THEN
      ALLOCATE(BUFFER(MAX(size1,size2),mpi__size-1) )
      IF (allocated(mpi__vrequest)) deallocate(mpi__vrequest)
      allocate(mpi__vrequest(mpi__size-1))
    ELSE
      ALLOCATE(BUFFER(MAX(size1,size2),1) )
    ENDIF











    !
    ! WORK DIVISION:
    ! INITIALISE X AND Y INDEXES
    I=1
    J=1
    ! LOOP THE THROUGH THE WHOLE DATA IN STEPS OF WORKING SLAVES (MASTER DOES NOT WORK)
    IF (mpi__myrank .EQ. 0) THEN
      PRINT*, 'Forward (start): ', NX*NY, SUM(AM_I_DONE)
      CALL SYSTEM_CLOCK(COUNT=TSTART,COUNT_RATE=TRATE,COUNT_MAX=TMAX)
      TIME1C=DBLE(TSTART)
    ENDIF
    !
call keep_workers_waiting()
    !
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
        ! MASTER SENDS I,J ATMOSPHERE MODELS TO SLAVES AND WAIT FOR...
        ! ...THE STOKES AND DERIVATIVES:
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
            IF (AM_I_DONE(J-1,I).GT.0.5) CYCLE
            ! NOW, IF WE PASS THIS CLAUSE, WE MUST ACCESS J-1 IN THE FOLLOWING!
            !
            ! Encode a I, J set of thermodynamic variables in a 2d array to be send
            !!!!!!!!CALL JOIN_MODEL2(I,J-1)
            CALL JOIN_MODEL2_NBCK(I,J-1,m+1)

!!!!!!            bpos1=0
!!!!!!            call mpi_pack(MODEL1D_SND, NZ*ATM_ARGS,MPI_REAL, buffer(1:size1,m+1) &
!!!!!!                , nbyt1, bpos1, MPI_COMM_WORLD, mpi__ierror)



            ! Send send atmospheric parameters:
          !!!!  CALL MPI_SEND(MODEL1D_SND,NZ*ATM_ARGS,MPI_REAL,M+1,1 &
          !!!!      ,MPI_COMM_WORLD,mpi__ierror)
            CALL MPI_ISEND(buffer(1:size1,m+1),size1,MPI_integer,M+1,21 &
                ,MPI_COMM_WORLD,mpi__vrequest(m+1),mpi__ierror)
            ! Update M
            M=M+1
          ENDDO
          !
          ! EVERYTHING SENT
          M=0
          DO WHILE (M.LT.max_work)
            CALL MPI_WAIT(mpi__vrequest(m+1), mpi__status, mpi__ierror)
            M=M+1
          ENDDO
          ! NOW, RECEIVE RESULTS:
          J=RJ
          I=RI
          !
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
!!!!!!!!!!!!!!            !
!!!!!!!!!!!!!!            ! 
!!!!!!!!!!!!!!            ! Send back full_stokes if needed:
!!!!!!!!!!!!!!            IF (FULL_STOKES.EQV..TRUE.) THEN
!!!!!!!!!!!!!!              CALL MPI_RECV(SYN5D(:,:,:,J-1,I),4*NUMW*NZ,MPI_REAL,M+1,5 &
!!!!!!!!!!!!!!                  ,MPI_COMM_WORLD,mpi__status,mpi__ierror)
!!!!!!!!!!!!!!            ENDIF ! Full Stokes.
!!!!!!!!!!!!!!            !
!!!!!!!!!!!!!!            IF (MTAULIN.EQV..TRUE.) THEN
!!!!!!!!!!!!!!              CALL MPI_RECV(TAU3DLIN(:,:,J-1,I),NUMW*NZ,MPI_REAL,M+1,6 &
!!!!!!!!!!!!!!                  ,MPI_COMM_WORLD,mpi__status,mpi__ierror)
!!!!!!!!!!!!!!            ENDIF


!!!!!!!!!!!              CALL MPI_RECV(BUFFER(1:SIZE2,m+1),SIZE2,MPI_INTEGER,M+1,20&
!!!!!!!!!!!                  ,MPI_COMM_WORLD, MPI__STATUS,mpi__ierror)
              CALL MPI_IRECV(BUFFER(1:SIZE2,M+1),SIZE2,MPI_INTEGER,M+1,20&
                  ,MPI_COMM_WORLD, mpi__vrequest(m+1),mpi__ierror)


            M=M+1
          ENDDO
          !
          ! EVERYTHING RECEIVED, NOW MOVE IT FROM TMP BUFFER TO TARGET STORAGE:
          J=RJ
          I=RI
          !
          M=0
          DO WHILE (M.LT.max_work)
            IF (MOD(J,NY+1).EQ.0) THEN
              I=I+1
              J=1
            ENDIF
            J=J+1
            IF (AM_I_DONE(J-1,I).GT.0.5) CYCLE


            CALL MPI_WAIT(mpi__vrequest(m+1), mpi__status, mpi__ierror)

              bpos2=0

              call mpi_unpack(BUFFER(:,M+1), nbyt2, bpos2, MODEL1D_SND &
                  ,NZ*ATM_ARGS, MPI_REAL, MPI_COMM_WORLD, mpi__ierror)
              IF ((MINVERSION.EQV..TRUE.).OR.(MSYNTHESIS.EQV..TRUE.)) THEN
                call mpi_unpack(BUFFER(:,M+1), nbyt2, bpos2, SYN3D(:,J-1,I)&
                    ,NFREQ,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi__ierror)
              ENDIF

              IF (MRESPFUNCT.EQV..TRUE.) THEN
                call mpi_unpack(BUFFER(:,M+1), nbyt2, bpos2, DSYN(:,:,J-1,I)&
                    ,NFREQ*IFREEP,MPI_REAL, MPI_COMM_WORLD, mpi__ierror)
              ENDIF

            
            CALL SPLIT_MODEL3D(I,J-1)
            
            !
            M=M+1
            !
          ENDDO
          ! DONE WITH THE MASTER
        ELSE
          !SLAVES:
          ! RECEIVE MODEL
          IF (mpi__myrank.LE.max_work) THEN
            !
         !!!!!   CALL MPI_RECV(MODEL1D_RCV,NZ*ATM_ARGS,MPI_REAL,0,1,MPI_COMM_WORLD &
         !!!!!       ,mpi__status, mpi__ierror)
            CALL MPI_RECV(BUFFER(1:SIZE1,1),SIZE1,MPI_INTEGER,0,21,MPI_COMM_WORLD &
                ,mpi__status, mpi__ierror)
            bpos1=0
            call mpi_unpack(BUFFER(1:SIZE1,1), nbyt1, bpos1, MODEL1D_RCV &
                ,NZ*ATM_ARGS,MPI_REAL, MPI_COMM_WORLD, mpi__ierror)
            !
            ! Solve the rte in 1d
            !
            CALL FORWARD1D(MODEL1D_RCV)
            !
            ! Send the products: stokes and derivatives to master:
            bpos2=0
            call mpi_pack(MODEL1D_RCV, NZ*ATM_ARGS,MPI_REAL, buffer(:,1) &
                , nbyt2, bpos2, MPI_COMM_WORLD, mpi__ierror)
            IF ((MINVERSION.EQV..TRUE.).OR.(MSYNTHESIS.EQV..TRUE.)) THEN
              call mpi_pack(RSYN1D,NFREQ,MPI_DOUBLE_PRECISION, buffer(:,1) &
                  , nbyt2, bpos2, MPI_COMM_WORLD, mpi__ierror)
            ENDIF

            IF (MRESPFUNCT.EQV..TRUE.) THEN
              call mpi_pack(EQVDSYN,NFREQ*IFREEP,MPI_REAL, buffer(:,1) &
                  , nbyt2, bpos2, MPI_COMM_WORLD, mpi__ierror)
            ENDIF
!!!!!!!!!!!!!!!!!            !
!!!!!!!!!!!!!!!!!            ! 
!!!!!!!!!!!!!!!!!            ! Send back full_stokes if needed:
!!!!!!!!!!!!!!!!!            IF (FULL_STOKES.EQV..TRUE.) THEN
!!!!!!!!!!!!!!!!!              CALL MPI_SEND(FS_SYN1D,4*NUMW*NZ,MPI_REAL,0,5 &
!!!!!!!!!!!!!!!!!                  ,MPI_COMM_WORLD,mpi__ierror)
!!!!!!!!!!!!!!!!!            ENDIF ! Full Stokes.
!!!!!!!!!!!!!!!!!            !
!!!!!!!!!!!!!!!!!            IF (MTAULIN.EQV..TRUE.) THEN
!!!!!!!!!!!!!!!!!              CALL MPI_SEND(TAULIN(:,:),NUMW*NZ,MPI_REAL,0,6 &
!!!!!!!!!!!!!!!!!                  ,MPI_COMM_WORLD,mpi__ierror)
!!!!!!!!!!!!!!!!!            ENDIF

              CALL MPI_SEND(BUFFER(1:SIZE2,1),SIZE2,MPI_INTEGER,0,20&
                  ,MPI_COMM_WORLD, mpi__ierror)
            !
          ENDIF
          !
        ENDIF
      ENDIF
    ENDDO
    IF (mpi__myrank .EQ. 0) THEN
      PRINT*, 'Forward (end): ', NX*NY, SUM(AM_I_DONE)
      N_FWD_MODELLING = N_FWD_MODELLING + 1
      CALL SYSTEM_CLOCK(COUNT=TSTART,COUNT_RATE=TRATE,COUNT_MAX=TMAX)
      TIME2C=DBLE(TSTART)
      PRINT*, 'Forward 3D: ', (TIME2C-TIME1C)/1000.
      !!!!    PRINT*, ' [>>>>>>>>> ', SUM(SYN3D), ' <<<<<<<<<<<<] '
      !!!!      IF (MRESPFUNCT.EQV..TRUE.) THEN
      !!!!    PRINT*, ' [>>>>>>>>> ', SUM(DSYN), ' <<<<<<<<<<<<] '
      !!!!      ENDIF
    ENDIF

    DEALLOCATE(BUFFER)
    !   
  END SUBROUTINE FORWARD3D
  !
  !------------------------------------------------
  !
  SUBROUTINE UPDATE_BEST_TO_CURRENT()
    !
    USE PHYS_PARAM
    USE user_mpi, ONLY: mpi__myrank
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      TEM3D=BEST_TEM3D
      PG3D=BEST_PG3D
      RHO3D=BEST_RHO3D
      PEL3D=BEST_PEL3D
      MW3D=BEST_MW3D
      BX3D=BEST_BX3D
      BY3D=BEST_BY3D
      BZ3D=BEST_BZ3D
      VZ3D=BEST_VZ3D
      TAU3D5=BEST_TAU3D5
      !
      SYN3D(:,:,:)=BEST_SYN(:,:,:)
      !
      DSYN(:,:,:,:)=BEST_DSYN(:,:,:,:)
      !
    ENDIF
    !
  END SUBROUTINE UPDATE_BEST_TO_CURRENT
  SUBROUTINE UPDATE_CURRENT_TO_BEST()
    !
    USE PHYS_PARAM
    USE user_mpi, ONLY: mpi__myrank
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      BEST_TEM3D=TEM3D
      BEST_PG3D=PG3D
      BEST_RHO3D=RHO3D
      BEST_PEL3D=PEL3D
      BEST_MW3D=MW3D
      BEST_BX3D=BX3D
      BEST_BY3D=BY3D
      BEST_BZ3D=BZ3D
      BEST_VZ3D=VZ3D
      BEST_TAU3D5=TAU3D5
      !
      BEST_SYN(:,:,:)=SYN3D(:,:,:)
      !
      BEST_DSYN(:,:,:,:)=DSYN(:,:,:,:)
      !
    ENDIF
    !
  END SUBROUTINE UPDATE_CURRENT_TO_BEST
  !
  !------------------------------------------------
  !
  SUBROUTINE JOIN_MODEL2(I,J)
    !
    USE GRID_PARAM, ONLY: ZZ
    USE PHYS_PARAM
    USE FORWARD_PARAM, ONLY: ATM_ARGS
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
  END SUBROUTINE JOIN_MODEL2
  !
  !------------------------------------------------
  !
  SUBROUTINE JOIN_MODEL2_NBCK(I,J,m)
    !
    USE GRID_PARAM, ONLY: ZZ, NZ
    USE PHYS_PARAM
    USE FORWARD_PARAM, ONLY: ATM_ARGS
    USE user_mpi
    !
    INTEGER, INTENT(IN)   :: I, J, m
    INTEGER               :: K
    REAL(SP), POINTER     :: ARR(:,:,:)
    !
    bpos1=0
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
       !MODEL1D_SND(:,K) = ARR (:,J,I)
!print*, K, ';', SUM(ARR(:,J,I))
       call mpi_pack(ARR(:,J,I), NZ,MPI_real, buffer(1:size1,m) &
           , nbyt1, bpos1, MPI_COMM_WORLD, mpi__ierror)

    !
       NULLIFY(ARR)
    ENDDO
    !MODEL1D_SND(:,13) = ZZ(:)
    call mpi_pack(REAL(zz(:)), NZ,MPI_REAL, buffer(1:size1,m) &
        , nbyt1, bpos1, MPI_COMM_WORLD, mpi__ierror)
!PRINT*, 13, ';', SUM(zz)
    !
  END SUBROUTINE JOIN_MODEL2_NBCK
  !
  !------------------------------------------------
  !
  SUBROUTINE DSYN_PROPAGATION(PG,KINIT)
    !
    USE GRID_PARAM, ONLY: NZ
    USE FORWARD_PARAM, ONLY: NUMW, DER_ARGS, GAIN1D
    USE PHYS_PARAM, ONLY: DSYN1D, EVOLG
    USE INVERT_PARAM, ONLY: INV_ATMPAR, NFREEV
    !
    REAL(SP), INTENT(IN), DIMENSION(NZ) :: PG
    INTEGER, INTENT(IN), OPTIONAL :: KINIT
    !
    INTEGER  :: K, L, P
    !
    REAL(DP), DIMENSION(4,4,NUMW) :: ACCMAT
    REAL(DP), DIMENSION(4,4)      :: IT_EVOL
    !
    INTEGER :: KMIN
    !
    IF (PRESENT(KINIT)) THEN
      KMIN=KINIT
    ELSE
      KMIN=1
    ENDIF

    !
    ACCMAT(:,:,:)=0.0D0
    DO K=1,4
      ACCMAT(K,K,:)=1.0D0
    ENDDO
    !
    DO K=NZ-1,KMIN,-1
      DO L=1,NUMW
        IT_EVOL(:,:)=EVOLG(:,:,L,K+1)
        ACCMAT(:,:,L)=MATMUL(ACCMAT(:,:,L),IT_EVOL)
        DO P=1,NFREEV
          DSYN1D(:,P,L,K)=MATMUL(ACCMAT(:,:,L), DSYN1D(:,P,L,K))*GAIN1D(:,L)
        ENDDO
      ENDDO
      IF (INV_ATMPAR(8).EQV..TRUE.) THEN
        DSYN1D(:,NFREEV,:,K)=PG(K)/PG(NZ)*DSYN1D(:,NFREEV,:,K)
      ENDIF
    ENDDO
    !
  END SUBROUTINE DSYN_PROPAGATION
  !
  !------------------------------------------------
  !
  SUBROUTINE OLD_DSYN_PROPAGATION(PG)
    !
    USE GRID_PARAM, ONLY: NZ
    USE FORWARD_PARAM, ONLY: NUMW, DER_ARGS
    USE PHYS_PARAM, ONLY: DSYN1D, EVOLG
    USE INVERT_PARAM, ONLY: INV_ATMPAR
    !
    REAL(SP), INTENT(IN), DIMENSION(NZ) :: PG
    !
    INTEGER  :: K, L, P
    !
    REAL(DP), DIMENSION(4,4,NUMW) :: ACCMAT
    REAL(DP), DIMENSION(4,4)      :: IT_EVOL
    !
    ACCMAT(:,:,:)=0.E0
    DO K=1,4
      ACCMAT(K,K,:)=1.E0
    ENDDO
    !
    DO K=NZ-1,1,-1
      DO L=1,NUMW
        IT_EVOL(:,:)=EVOLG(:,:,L,K+1)
        ACCMAT(:,:,L)=MATMUL(ACCMAT(:,:,L),IT_EVOL)
        DO P=1,DER_ARGS
          IF (P.NE.2) THEN
            IF (INV_ATMPAR(P).EQV..FALSE.) CYCLE
          ELSE
            IF ((INV_ATMPAR(P).EQV..FALSE.).AND.(INV_ATMPAR(8).EQV..FALSE.)) CYCLE
          ENDIF
          DSYN1D(:,P,L,K)=MATMUL(ACCMAT(:,:,L), DSYN1D(:,P,L,K))
        ENDDO
      ENDDO
      IF (INV_ATMPAR(8).EQV..TRUE.) THEN
        DSYN1D(:,8,:,K)=DSYN1D(:,8,:,K)*PG(K)/PG(NZ)*DSYN1D(:,2,:,K)
      ENDIF
    ENDDO
    !
  END SUBROUTINE OLD_DSYN_PROPAGATION
  !
  !------------------------------------------------
  !
  SUBROUTINE CALCULATE_EQUIVALENT_DSYN()
    !
    USE GRID_PARAM, ONLY: NZ
    USE HEIGHT_HANDLER, ONLY: TM_COEFS, PG_COEFS, P0_COEFS, RH_COEFS, BX_COEFS &
        , BY_COEFS, BZ_COEFS, VZ_COEFS
    USE PHYS_PARAM, ONLY: DSYN1D, EQVDSYN
    USE INVERT_PARAM, ONLY: INV_ATMPAR, INV_STK, NSLAB_PER_FREEV
    USE FORWARD_PARAM, ONLY: NUMW, NUML, PIXEL_INI, PIXEL_END, GAIN1D, LSF_VALID
    USE CODE_MODES, ONLY: MLSF
    !
    USE USER_FFTW3
    !
    INTEGER   :: IP, II, ISCNT, NJ, IPCNT, DFFST, JCNT
    REAL(DP),DIMENSION(:,:), POINTER  :: PTR_COEF
    !
    INTEGER   :: JJ, WI, WF, JK
    REAL(DP)  :: DPSYN(NUMW)
    REAL(DP)  :: DUM(NUMW,NZ)
    REAL(DP)  :: GDUM(NUMW)
    !
    !IF (MRESPFUNCT.EQV..TRUE.) THEN
      EQVDSYN(:,:)=0.E0
      IPCNT=0
      JCNT=0
      DO IP=1,SIZE(INV_ATMPAR)
        !
        IF (INV_ATMPAR(IP).EQV..TRUE.) THEN
          IPCNT=IPCNT+1
          SELECT CASE(IP)
            CASE(1)
              PTR_COEF=>TM_COEFS
            CASE(2)
              PTR_COEF=>PG_COEFS
            CASE(3)
              PTR_COEF=>RH_COEFS
            CASE(4)
              PTR_COEF=>BX_COEFS
            CASE(5)
              PTR_COEF=>BY_COEFS
            CASE(6)
              PTR_COEF=>BZ_COEFS
            CASE(7)
              PTR_COEF=>VZ_COEFS
            CASE(8)
              PTR_COEF=>P0_COEFS
          ENDSELECT
          !
          DFFST=NSLAB_PER_FREEV(IPCNT)-1
          NJ=SIZE(PTR_COEF,2)
          ISCNT=0
          !
          DO II=1,SIZE(INV_STK)
            !
            IF (INV_STK(II).EQV..TRUE.) THEN
              !
              DUM(:,:)=DSYN1D(II,IP,:,:)
              GDUM(:)=GAIN1D(II,:)
              !!!!!!!!!!DUM(:,:)=DSYN1D(II,IPCNT,:,:)
              !
              EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+1:JCNT+NJ) &
                  = MATMUL(DUM,PTR_COEF(:,:))

              DO JJ=1,NJ
                EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+JJ) &
                    = EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+JJ) * GDUM(:)
              ENDDO

              !
              IF (MLSF.EQV..TRUE.) THEN
                DO JK=1,NJ
                  FFTY1D(:)=0.0D0
                  DPSYN(:)=DBLE(EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+JK))
                  ! Loop through all the given spectral regions:
                  DO JJ=1,NUML
                    !
            IF (.NOT.LSF_VALID(JJ)) CYCLE
                    !
                    WI=PIXEL_INI(JJ)
                    WF=PIXEL_END(JJ)
                    !
                    ! Direct transformation:
                    CALL FFTW_EXECUTE_DFT_R2C(GPLAN1D(1,JJ),DPSYN(WI:WF) &
                        , FFTY1D(FFTINI(JJ):FFTEND(JJ)))
!                  ENDDO
                  !
                  ! Once DFT are stored for all the spectral regions ...
                  ! ... perform the transformed multiplication:
                  FFTY1D(FFTINI(JJ):FFTEND(JJ))=FFTY1D(FFTINI(JJ):FFTEND(JJ))*FFTK1D(FFTINI(JJ):FFTEND(JJ))
                  ! And now, come back to the real space:
                  ! Loop through all the given spectral regions:
!                  DO JJ=1,NUML
                    !
                    WI=PIXEL_INI(JJ)
                    WF=PIXEL_END(JJ)
                    !
                    ! Direct transformation:
                    CALL FFTW_EXECUTE_DFT_C2R(GPLAN1D(2,JJ) &
                        ,FFTY1D(FFTINI(JJ):FFTEND(JJ)),DPSYN(WI:WF))
                    DPSYN(WI:WF)=DPSYN(WI:WF)!/DBLE(WF-WI+1)
                  ENDDO
                  EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+JK)=REAL(DPSYN)
                  !
                ENDDO ! Equivalent heights
                !
              ENDIF ! Convolution with LSF
!              EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+1:JCNT+NJ) &
!                  = MATMUL(DSYN1D(II,IP,:,:),PTR_COEF(:,:))
              !
              ! Update counter
              ISCNT=ISCNT+1
              !
            ENDIF ! Invert Stokes
          ENDDO !Stokes
          JCNT=JCNT+(DFFST+1)
          NULLIFY(PTR_COEF)
        ENDIF !Invert atm. param.
      ENDDO ! Atm. param.
    !ENDIF
    !
  END SUBROUTINE CALCULATE_EQUIVALENT_DSYN
  !
  !------------------------------------------------
  !
  SUBROUTINE OLD_CALCULATE_EQUIVALENT_DSYN()
    !
    USE HEIGHT_HANDLER, ONLY: TM_COEFS, PG_COEFS, P0_COEFS, RH_COEFS, BX_COEFS &
        , BY_COEFS, BZ_COEFS, VZ_COEFS
    USE PHYS_PARAM, ONLY: DSYN1D, EQVDSYN
    USE INVERT_PARAM, ONLY: INV_ATMPAR, INV_STK, NSLAB_PER_FREEV
    USE FORWARD_PARAM, ONLY: NUMW, NUML, PIXEL_INI, PIXEL_END
    USE CODE_MODES, ONLY: MLSF
    !
    USE USER_FFTW3
    !
    INTEGER   :: IP, II, ISCNT, NJ, IPCNT, DFFST, JCNT
    REAL(DP),DIMENSION(:,:), POINTER  :: PTR_COEF
    !
    INTEGER   :: JJ, WI, WF, JK
    REAL(DP)  :: DPSYN(NUMW)
    !
    !IF (MRESPFUNCT.EQV..TRUE.) THEN
      EQVDSYN(:,:)=0.E0
      IPCNT=0
      JCNT=0
      DO IP=1,SIZE(INV_ATMPAR)
        !
        IF (INV_ATMPAR(IP).EQV..TRUE.) THEN
          IPCNT=IPCNT+1
          SELECT CASE(IP)
            CASE(1)
              PTR_COEF=>TM_COEFS
            CASE(2)
              PTR_COEF=>PG_COEFS
            CASE(3)
              PTR_COEF=>RH_COEFS
            CASE(4)
              PTR_COEF=>BX_COEFS
            CASE(5)
              PTR_COEF=>BY_COEFS
            CASE(6)
              PTR_COEF=>BZ_COEFS
            CASE(7)
              PTR_COEF=>VZ_COEFS
            CASE(8)
              PTR_COEF=>P0_COEFS
          ENDSELECT
          !
          DFFST=NSLAB_PER_FREEV(IPCNT)-1
          NJ=SIZE(PTR_COEF,2)
          ISCNT=0
          !
          DO II=1,SIZE(INV_STK)
            !
            IF (INV_STK(II).EQV..TRUE.) THEN
              !
              EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+1:JCNT+NJ) &
                  = MATMUL(DSYN1D(II,IP,:,:),PTR_COEF(:,:))
              !
              IF (MLSF.EQV..TRUE.) THEN
                DO JK=1,NJ
                  FFTY1D(:)=0.0D0
                  DPSYN(:)=DBLE(EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+JK))
                  ! Loop through all the given spectral regions:
                  DO JJ=1,NUML
                    !
                    WI=PIXEL_INI(JJ)
                    WF=PIXEL_END(JJ)
                    !
                    ! Direct transformation:
                    CALL FFTW_EXECUTE_DFT_R2C(GPLAN1D(1,JJ),DPSYN(WI:WF) &
                        , FFTY1D(FFTINI(JJ):FFTEND(JJ)))
                  ENDDO
                  !
                  ! Once DFT are stored for all the spectral regions ...
                  ! ... perform the transformed multiplication:
                  FFTY1D=FFTY1D*FFTK1D
                  ! And now, come back to the real space:
                  ! Loop through all the given spectral regions:
                  DO JJ=1,NUML
                    !
                    WI=PIXEL_INI(JJ)
                    WF=PIXEL_END(JJ)
                    !
                    ! Direct transformation:
                    CALL FFTW_EXECUTE_DFT_C2R(GPLAN1D(2,JJ) &
                        ,FFTY1D(FFTINI(JJ):FFTEND(JJ)),DPSYN(WI:WF))
                    DPSYN(WI:WF)=DPSYN(WI:WF)!/DBLE(WF-WI+1)
                  ENDDO
                  EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+JK)=REAL(DPSYN)
                  !
                ENDDO ! Equivalent heights
                !
              ENDIF ! Convolution with LSF
!              EQVDSYN(ISCNT*NUMW+1:(ISCNT+1)*NUMW,JCNT+1:JCNT+NJ) &
!                  = MATMUL(DSYN1D(II,IP,:,:),PTR_COEF(:,:))
              !
              ! Update counter
              ISCNT=ISCNT+1
              !
            ENDIF ! Invert Stokes
          ENDDO !Stokes
          JCNT=JCNT+(DFFST+1)
          NULLIFY(PTR_COEF)
        ENDIF !Invert atm. param.
      ENDDO ! Atm. param.
    !ENDIF
    !
  END SUBROUTINE OLD_CALCULATE_EQUIVALENT_DSYN
  !
  !------------------------------------------------
  !
  SUBROUTINE CALCULATE_EQUIVALENT_SYN1D()
    !
    USE PHYS_PARAM, ONLY: SYN1D, RSYN1D
    USE INVERT_PARAM, ONLY: INV_STK
    USE FORWARD_PARAM, ONLY: NUMW, NUML, PIXEL_INI, PIXEL_END, LSF_VALID, GAIN1D
    USE CODE_MODES, ONLY: MLSF
    !
    USE USER_FFTW3
    !
    INTEGER   :: IP, II, ISCNT, NJ, IPCNT, DFFST, JCNT
    !
    INTEGER   :: JJ, WI, WF
    REAL(DP)  :: DPSYN(NUMW)
    !
    ISCNT=0
    !
    SYN1D(:,:)=SYN1D(:,:)*GAIN1D(:,:)
    DO II=1,SIZE(INV_STK)
      IF (INV_STK(II).EQV..TRUE.) THEN
        !
        ! Here, we have to be careful if a LSF is given:
!PRINT*, 'Previous: ', SYN1D(II,:)
        IF (MLSF.EQV..TRUE.) THEN
          FFTY1D(:)=0.0D0
          DPSYN(:)=DBLE(SYN1D(II,:))
          ! Loop through all the given spectral regions:
          DO JJ=1,NUML
            !
            IF (.NOT.LSF_VALID(JJ)) CYCLE
            !
            WI=PIXEL_INI(JJ)
            WF=PIXEL_END(JJ)
            !
            ! Direct transformation:
            CALL FFTW_EXECUTE_DFT_R2C(GPLAN1D(1,JJ),DPSYN(WI:WF) &
                , FFTY1D(FFTINI(JJ):FFTEND(JJ)))
            !
            ! Once DFT are stored for all the spectral regions ...
            ! ... perform the transformed multiplication:
            FFTY1D(FFTINI(JJ):FFTEND(JJ))=FFTY1D(FFTINI(JJ):FFTEND(JJ)) &
                *FFTK1D(FFTINI(JJ):FFTEND(JJ))
            ! And now, come back to the real space:
            ! Loop through all the given spectral regions:
            !
            WI=PIXEL_INI(JJ)
            WF=PIXEL_END(JJ)
            !
            ! Direct transformation:
            CALL FFTW_EXECUTE_DFT_C2R(GPLAN1D(2,JJ) &
                ,FFTY1D(FFTINI(JJ):FFTEND(JJ)),DPSYN(WI:WF))
          ENDDO
          SYN1D(II,:)=REAL(DPSYN)
        ENDIF
        !
        ! Once convolution is done (if required)
        RSYN1D(ISCNT*NUMW+1:(ISCNT+1)*NUMW) &
            = SYN1D(II,:)
        ISCNT=ISCNT+1
      ENDIF ! Invert Stokes
    ENDDO !Stokes
    !
  END SUBROUTINE CALCULATE_EQUIVALENT_SYN1D
  !
  !------------------------------------------------
  !
  SUBROUTINE FORWARD1D(MODEL)
    !
    USE GRID_PARAM, ONLY: NZ, ZZ
    USE RTESOLVER, ONLY: SET_BOUNDARY, ANALYTICAL_SOLVER, OLD_ANALYTICAL_SOLVER, DPGDTK_HE
    USE LINE_OPACITY, ONLY: OPAC_LINE
    USE CONT_OPACITY, ONLY: OPAC_CONTINUUM
    USE DAMPING, ONLY: GET_DAMPING
    USE ABSORPTION_MATRIX, ONLY: GET_ABS_MAT, ETAI
    USE ATM_PARAM
    USE PHYS_PARAM
    USE HYDROSTATIC_TOOLS, ONLY: RK4_INTEGRATION_LOG &
        , RK4_INTEGRATION_LOG_LTAUM05
    USE CODE_MODES
    USE CHEMICAL_EQUILIBRIUM, ONLY: GET_NHNEMU_G, GET_DNHYDDTEMP
    USE FORWARD_PARAM
    USE CONS_PARAM, ONLY: KBOL, MAMU
    USE INVERT_PARAM, ONLY: NFREEV
    !
    REAL(SP), INTENT(INOUT), TARGET  :: MODEL(NZ,ATM_ARGS)
    !
    INTEGER                          :: K, L, B
    REAL(DP)                         :: NHYD, NELEC, DAMP
    REAL(DP),DIMENSION(NZ)           :: VNHYD, VNELEC
    !
    REAL(DP) :: KC5THRESHOLD
    INTEGER :: KTHRESHOLD
    !
    KC5THRESHOLD=2.0D10
    !
    CALL SPLIT_MODEL(MODEL)
    !
    KTHRESHOLD=1
    !
    ! FIRST CHECK IF WE NEED TO IMPOSE HE:
    IF (HYDROSTATIC.EQV..TRUE.) THEN
      !
      ALLOCATE(DPGDTK_HE(NZ))
      DPGDTK_HE(:)=0.0D0
      !


!        CALL RK4_INTEGRATION_LOG(NZ,TEM,PG,RHO,PEL,MW,ZZ)
      !CALL RK4_INTEGRATION_LOG_BOTTOM(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ)
      IF (HYDRO_TOP.EQV..TRUE.) THEN
        CALL RK4_INTEGRATION_LOG(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ)
      ELSE
        CALL RK4_INTEGRATION_LOG_LTAUM05(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ,TAU5)
      ENDIF
      VNELEC=PEL(:)/KBOL/TEM(:)
    ENDIF
    !
    ! WE NOW CHECK WHETHER HE IS THE ONLY THING TO BE DONE OR NOT
    IF (MGETHEQ.EQV..FALSE.) THEN
      !
      DO K=1,NZ
         !
         IF (HYDROSTATIC.EQV..FALSE.) THEN
           CALL GET_NHNEMU_G(TEM(K),PG(K),RHO(K),NHYD,NELEC,MW(K),PEL(K))
         ELSE
           NHYD=VNHYD(K)
           NELEC=VNELEC(K)
         ENDIF
         !
         IF (MRESPFUNCT.EQV..TRUE.) THEN
           CALL GET_DNHYDDTEMP(TEM(K),PG(K),RHO(K),MW(K),NELEC)
           CALL GET_HYD_DERIVATIVES(VNHYD,VNELEC,K)
         ENDIF
         !
         ! Continuum opacity at 5000 A (only for conversion to TAU5)
         CALL OPAC_CONTINUUM(TEM(K),5000.0D0,NELEC,NHYD,KC5(K))
         KC5(K)=KC5(K)/NHYD*(1.0D0/MAMU)/MW(K)
         !
! Skip too deep layers, as their opacity is so large that they become the boundary condition:
IF (KC5(K).LT.KC5THRESHOLD) THEN
         !
         IF ((MSYNTHESIS.EQV..TRUE.).OR.(MINVERSION.EQV..TRUE.)) THEN
           ! Clean KLIN & KC
           KLIN(:,:)=0.0D0
           KC(:)=0.0D0
           !
           ! JUST FOR SPECTRAL LINE UNDER CONSIDERATION
           DO L=1,NUML
              !
              CALL OPAC_CONTINUUM(TEM(K),LINE_L0(BLENDSID(L,1)),NELEC&
                   ,NHYD,KC(L))
              !
              DO B=1,BLENDSMAX+1
                !
                IF (BLENDSID(L,B).GT.0.1) THEN
                  !
                  CALL OPAC_LINE(BLENDSID(L,B),TEM(K),NELEC,NHYD,KLIN(L,B))
                  !
                  CALL GET_DAMPING(LINE_ZN(BLENDSID(L,B)), LINE_L0(BLENDSID(L,B))&
                       *1D-8, LINE_ION(BLENDSID(L,B)), EPLOW(BLENDSID(L,B)), NHYD&
                       , TEM(K), SIGMA(BLENDSID(L,B)), ALPHA(BLENDSID(L,B)), DAMP)
                  !
                  CALL GET_ABS_MAT(K,L,B,BLENDSID(L,B),DAMP,BLENDSDIFF(L,B))
                  !
                ENDIF
              ENDDO
              !
           ENDDO
ELSE
KTHRESHOLD=K ! It stores the largest height that we can neglect
ENDIF

           ! FOR EVERY WAVELENGTH
           DO L=1,NUMW
              IF (MTAULIN.EQV..TRUE.) THEN
                KLINTAU(L,K)=ETAI(L)/NHYD*(1.0D0/MAMU)/MW(K)
              ENDIF
              !
              IF ( (K.EQ.1) .OR. (KC5(K).GE.KC5THRESHOLD) ) THEN
                 CALL SET_BOUNDARY(K,L)
              ELSE
                 !CALL ANALYTICAL_SOLVER(K,L)
                 CALL OLD_ANALYTICAL_SOLVER(K,L)
              ENDIF
              !
           ENDDO
           ! 
           ! DO I WANT TO STORE THE WHOLE STOKES VECTOR (i.e. THE WHOLE HEIGHT GRID?)?
           IF (FULL_STOKES.EQV..TRUE.) THEN
             FS_SYN1D(:,:,K)=REAL(SYN1D(:,:))*GAIN1D(:,:)
           ENDIF
           !
         ENDIF
      ENDDO
      !
      ! NORMALIZATION:
      IF ((MSYNTHESIS.EQV..TRUE.).OR.(MINVERSION.EQV..TRUE.)) THEN
        !
        IF (MRESPFUNCT.EQV..TRUE.) THEN
          CALL OLD_DSYN_PROPAGATION(PG)
          CALL CALCULATE_EQUIVALENT_DSYN()
        ENDIF
        ! END NEW TEST
        CALL CALCULATE_EQUIVALENT_SYN1D()
        !
      ENDIF
      ! Conversion to optical depth (not really needed but for test purposes)
      TAU5(NZ)=1E-9
      DO K=NZ-1,1,-1
         TAU5(K)=TAU5(K+1)&
             + (((REAL(KC5(K))*RHO(K)+REAL(KC5(K+1))*RHO(K+1))/2.E0)&
             * ((ZZ(K+1)-ZZ(K))*1000.E0*100.E0))
      ENDDO
      TAU5(:)=LOG10(TAU5(:))
      !
      IF (MTAULIN.EQV..TRUE.) THEN
        TAULIN(:,NZ)=1E-9
        DO K=NZ-1,1,-1
           TAULIN(:,K)=TAULIN(:,K+1)&
               + (((REAL(KLINTAU(:,K))*RHO(K)+REAL(KLINTAU(:,K+1))*RHO(K+1))/2.E0)&
               * ((ZZ(K+1)-ZZ(K))*1000.E0*100.E0))
        ENDDO
        TAULIN(:,:)=LOG10(TAULIN(:,:))
      ENDIF
      !
    ENDIF
    !
    !
    !
    IF (HYDROSTATIC.EQV..TRUE.) THEN
      DEALLOCATE(DPGDTK_HE)
    ENDIF
    !
    !
    !
    !
  END SUBROUTINE FORWARD1D
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_HYD_DERIVATIVES(HH, EE, K)
    !
    USE GRID_PARAM, ONLY: NZ, ZZ
    USE CODE_MODES, ONLY: HYDROSTATIC
    USE ATM_PARAM, ONLY: TEM,PG,RHO,MW
    USE RTESOLVER, ONLY: DPGDTK_HE
    USE CONS_PARAM, ONLY: KBOL, MAMU, GRAV
    !
    REAL(DP), DIMENSION(NZ), INTENT(IN) :: HH, EE
    INTEGER, INTENT(IN) :: K
    !
    INTEGER :: I
    !
    LOGICAL :: CENT_DER
    !
    REAL(DP) :: CTE!NHYD1,NELEC1,PATOM1,PELEC1,DENS1,MOLECW1
    !

    IF (HYDROSTATIC.EQV..TRUE.) THEN
      DPGDTK_HE(:) = 0.0D0
      IF (K.GT.1) THEN
        CTE = MAMU / KBOL * GRAV * (ZZ(K)-ZZ(K-1)) * 1.0D3 * 1.0D2
        DO I=1,K-1
          DPGDTK_HE(I) = (PG(I) / PG(K)) * (- PG(K) * MW(K) * CTE / TEM(K) / TEM(K))
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE GET_HYD_DERIVATIVES
  !
  !------------------------------------------------
  !
  SUBROUTINE OLD_FORWARD1D(MODEL)
    !
    USE GRID_PARAM, ONLY: NZ, ZZ
    USE RTESOLVER, ONLY: SET_BOUNDARY, ANALYTICAL_SOLVER
    USE LINE_OPACITY, ONLY: OPAC_LINE
    USE CONT_OPACITY, ONLY: OPAC_CONTINUUM
    USE DAMPING, ONLY: GET_DAMPING
    USE ABSORPTION_MATRIX, ONLY: GET_ABS_MAT, ETAI
    USE ATM_PARAM
    USE PHYS_PARAM
    USE HYDROSTATIC_TOOLS, ONLY: RK4_INTEGRATION_LOG &
        , RK4_INTEGRATION_LOG_LTAUM05
    USE CODE_MODES
    USE CHEMICAL_EQUILIBRIUM, ONLY: GET_NHNEMU_G, GET_DNHYDDTEMP
    USE FORWARD_PARAM
    USE CONS_PARAM, ONLY: KBOL, MAMU
    !
    REAL(SP), INTENT(INOUT), TARGET  :: MODEL(NZ,ATM_ARGS)
    !
    INTEGER                          :: K, L, B
    REAL(DP)                         :: NHYD, NELEC, DAMP
    REAL(DP),DIMENSION(NZ)           :: VNHYD, VNELEC
    !
    CALL SPLIT_MODEL(MODEL)
    !
    !
    ! FIRST CHECK IF WE NEED TO IMPOSE HE:
    IF (HYDROSTATIC.EQV..TRUE.) THEN
!        CALL RK4_INTEGRATION_LOG(NZ,TEM,PG,RHO,PEL,MW,ZZ)
      !CALL RK4_INTEGRATION_LOG_BOTTOM(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ)
      IF (HYDRO_TOP.EQV..TRUE.) THEN
        CALL RK4_INTEGRATION_LOG(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ)
      ELSE
        CALL RK4_INTEGRATION_LOG_LTAUM05(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ,TAU5)
      ENDIF
      VNELEC=PEL(:)/KBOL/TEM(:)
    ENDIF
    !
    ! WE NOW CHECK WHETHER HE IS THE ONLY THING TO BE DONE OR NOT
    IF (MGETHEQ.EQV..FALSE.) THEN
      !
      DO K=1,NZ
         !
         IF (HYDROSTATIC.EQV..FALSE.) THEN
           CALL GET_NHNEMU_G(TEM(K),PG(K),RHO(K),NHYD,NELEC,MW(K),PEL(K))
         ELSE
           NHYD=VNHYD(K)
           NELEC=VNELEC(K)
         ENDIF
         !
         IF (MRESPFUNCT.EQV..TRUE.) THEN
           CALL GET_DNHYDDTEMP(TEM(K),PG(K),RHO(K),MW(K),NELEC)
         ENDIF
         !
         ! Continuum opacity at 5000 A (only for conversion to TAU5)
         CALL OPAC_CONTINUUM(TEM(K),5000.0D0,NELEC,NHYD,KC5(K))
         KC5(K)=KC5(K)/NHYD*(1.0D0/MAMU)/MW(K)
         !
         IF ((MSYNTHESIS.EQV..TRUE.).OR.(MINVERSION.EQV..TRUE.)) THEN
           ! Clean KLIN & KC
           KLIN(:,:)=0.0D0
           KC(:)=0.0D0
           !
           ! JUST FOR SPECTRAL LINE UNDER CONSIDERATION
           DO L=1,NUML
              !
              CALL OPAC_CONTINUUM(TEM(K),LINE_L0(BLENDSID(L,1)),NELEC&
                   ,NHYD,KC(L))
              !
              DO B=1,BLENDSMAX+1
                !
                IF (BLENDSID(L,B).GT.0.1) THEN
                  !
                  CALL OPAC_LINE(BLENDSID(L,B),TEM(K),NELEC,NHYD,KLIN(L,B))
                  !
                  CALL GET_DAMPING(LINE_ZN(BLENDSID(L,B)), LINE_L0(BLENDSID(L,B))&
                       *1D-8, LINE_ION(BLENDSID(L,B)), EPLOW(BLENDSID(L,B)), NHYD&
                       , TEM(K), SIGMA(BLENDSID(L,B)), ALPHA(BLENDSID(L,B)), DAMP)
                  !
                  CALL GET_ABS_MAT(K,L,B,BLENDSID(L,B),DAMP,BLENDSDIFF(L,B))
                  !
                ENDIF
              ENDDO
              !
           ENDDO
           ! FOR EVERY WAVELENGTH
           DO L=1,NUMW
              IF (MTAULIN.EQV..TRUE.) THEN
                KLINTAU(L,K)=ETAI(L)/NHYD*(1.0D0/MAMU)/MW(K)
              ENDIF
              !
              IF (K.EQ.1) THEN
                 CALL SET_BOUNDARY(K,L)
              ELSE
                 CALL ANALYTICAL_SOLVER(K,L)
              ENDIF
              !
           ENDDO
           ! 
           ! DO I WANT TO STORE THE WHOLE STOKES VECTOR (i.e. THE WHOLE HEIGHT GRID?)?
           IF (FULL_STOKES.EQV..TRUE.) THEN
             FS_SYN1D(:,:,K)=REAL(SYN1D(:,:))
           ENDIF
           !
         ENDIF
      ENDDO
      !
      ! NORMALIZATION:
      IF ((MSYNTHESIS.EQV..TRUE.).OR.(MINVERSION.EQV..TRUE.)) THEN
        SYN1D(:,:)=SYN1D*GAIN1D
        !
        ! NEW TEST:
        IF (FULL_STOKES.EQV..TRUE.) THEN
          DO K=1,NZ
            FS_SYN1D(:,:,K)=FS_SYN1D(:,:,K)*GAIN1D(:,:)
          ENDDO
        ENDIF
        !
        IF (MRESPFUNCT.EQV..TRUE.) THEN
          CALL DSYN_PROPAGATION(PG)
          !
          DO K=1,NZ
            DO L=1,DER_ARGS
              DSYN1D(:,L,:,K)=DSYN1D(:,L,:,K)*GAIN1D(:,:)
            ENDDO
          ENDDO
          !
          CALL CALCULATE_EQUIVALENT_DSYN()
        ENDIF
        ! END NEW TEST
        CALL CALCULATE_EQUIVALENT_SYN1D()
        !
      ENDIF
      ! Conversion to optical depth (not really needed but for test purposes)
      TAU5(NZ)=1E-9
      DO K=NZ-1,1,-1
         TAU5(K)=TAU5(K+1)&
             + (((REAL(KC5(K))*RHO(K)+REAL(KC5(K+1))*RHO(K+1))/2.E0)&
             * ((ZZ(K+1)-ZZ(K))*1000.E0*100.E0))
      ENDDO
      TAU5(:)=LOG10(TAU5(:))
      !
      IF (MTAULIN.EQV..TRUE.) THEN
        TAULIN(:,NZ)=1E-9
        DO K=NZ-1,1,-1
           TAULIN(:,K)=TAULIN(:,K+1)&
               + (((REAL(KLINTAU(:,K))*RHO(K)+REAL(KLINTAU(:,K+1))*RHO(K+1))/2.E0)&
               * ((ZZ(K+1)-ZZ(K))*1000.E0*100.E0))
        ENDDO
        TAULIN(:,:)=LOG10(TAULIN(:,:))
      ENDIF
      !
    ENDIF
    !
  END SUBROUTINE OLD_FORWARD1D
  !
  !------------------------------------------------
  !
  SUBROUTINE SPLIT_MODEL3D(I,J)
    !
    USE CONS_PARAM, ONLY: SP
    USE GRID_PARAM, ONLY: ZZ
    USE PHYS_PARAM
    USE FORWARD_PARAM, ONLY: ATM_ARGS
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
       ENDSELECT
       ARR (:,J,I) = MODEL1D_SND(:,K)
       NULLIFY(ARR)
    ENDDO
    ZZ(:)=MODEL1D_SND(:,13)
    !
  END SUBROUTINE SPLIT_MODEL3D
  !
  !================================================
  !
END MODULE FORWARD
!

