!
MODULE coupled_matrix_inversion
  !
  USE CONS_PARAM, ONLY: SP, DP
  USE ALLOCATE_UTILS
  USE MISC, ONLY: WRITE_BIN
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:)  :: ALP
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: BETA
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: HINV
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: HPREINV
  !REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)  :: WINV
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: WINV
  !
  INTEGER :: NPAR,NY,NX
  INTEGER :: NB,NTR,INVSIZE
  !
  REAL(DP) :: SVDTOL
  REAL(DP) :: LAMBDAP
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: PSF
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)  :: IPSF
  INTEGER,DIMENSION(2)                   :: CTP
  INTEGER :: NPY,NPX
  !
  REAL(DP) :: SIZEFAC
  !
  PUBLIC :: BUILD_ALPHA, BUILD_BETA
  PUBLIC :: BUILD_HINV, CALCULATE_DELTA
  !
  PRIVATE :: ALP, BETA, HINV, WINV, HPREINV
  PRIVATE :: NPAR,NY,NX
  PRIVATE :: NB,NTR,INVSIZE
  PRIVATE :: SVDTOL, LAMBDAP
  PRIVATE :: CTP, PSF, NPY, NPX
  PRIVATE :: SIZEFAC
  !PRIVATE :: SCTP, SPSF, NSPY, NSPX
  !
  PRIVATE :: INV_MAT, ALPHA_DOT_DELTA, HINV_DOT_BETA, DEALLOCATE_ALL
  !
  CONTAINS
  !
  SUBROUTINE BUILD_ALPHA(SH,JAC)
    !
    INTEGER,DIMENSION(4),INTENT(IN) :: SH
    REAL(DP), INTENT(IN), DIMENSION(SH(1),SH(2),SH(3),SH(4))  :: JAC
    !
    INTEGER :: I,J
    !
    CALL ALLOCATE_4D_DP(ALP,SH(2),SH(2),SH(3),SH(4),'C')
    NPAR=SH(2)
    NY=SH(3)
    NX=SH(4)
    SIZEFAC=1.0D0/DBLE(NY)/DBLE(NX)/DBLE(NPAR)
    !
    DO I=1,NX
      DO J=1,NY
        ALP(:,:,J,I)=MATMUL(TRANSPOSE(JAC(:,:,J,I)), JAC(:,:,J,I))
      ENDDO
    ENDDO
    !
    CALL WRITE_BIN('alpha.bin',SIZE(SHAPE(ALP)), SHAPE(ALP) &
        , SIZE(ALP,KIND=8),REAL(ALP),3000,3)
    !
  END SUBROUTINE BUILD_ALPHA
  !
  !
  SUBROUTINE BUILD_BETA(SH4,JAC,SH3,DAR)
    !
    INTEGER,DIMENSION(4),INTENT(IN) :: SH4
    REAL(DP), INTENT(IN), DIMENSION(SH4(1),SH4(2),SH4(3),SH4(4))  :: JAC
    INTEGER,DIMENSION(3),INTENT(IN) :: SH3
    REAL(DP), INTENT(IN), DIMENSION(SH3(1),SH3(2),SH3(3))  :: DAR
    !
    INTEGER :: I,J
    !
    CALL ALLOCATE_3D_DP(BETA,NPAR,NY,NX,'C')
    DO I=1,NX
      DO J=1,NY
        BETA(:,J,I)=MATMUL(TRANSPOSE(JAC(:,:,J,I)), DAR(:,J,I))
      ENDDO
    ENDDO
    !
  !><  CALL WRITE_BIN('ojac.bin',SIZE(SHAPE(JAC)), SHAPE(JAC) &
  !><      , SIZE(JAC,KIND=8),REAL(JAC),3000,3)
  !><  CALL WRITE_BIN('odar.bin',SIZE(SHAPE(DAR)), SHAPE(DAR) &
  !><      , SIZE(DAR,KIND=8),REAL(DAR),3000,3)
  !><  CALL WRITE_BIN('obeta.bin',SIZE(SHAPE(BETA)), SHAPE(BETA) &
  !><      , SIZE(BETA,KIND=8),REAL(BETA),3000,3)
    !
  END SUBROUTINE BUILD_BETA
  !
  !
  SUBROUTINE BUILD_HINV(SHP,EPSF2D,ECTP,ENB,ESVDTOL,ELAMBDAP)!,ESCTP)
    !
    ! Alpha:
    ! PSF:
    INTEGER,DIMENSION(2),INTENT(IN) :: SHP
    REAL(DP), INTENT(IN), DIMENSION(SHP(1),SHP(2))  :: EPSF2D
    INTEGER,DIMENSION(2),INTENT(IN) :: ECTP
    ! Hinvs:
    INTEGER,INTENT(IN) :: ENB
    REAL(DP),INTENT(IN) :: ESVDTOL
    REAL(DP),INTENT(IN) :: ELAMBDAP
    !INTEGER,INTENT(IN) :: ESCTP
    !
    INTEGER :: ITN, IT2, IT3, IP2, IP3
    INTEGER :: ITN2, ITN3, ITP2, ITP3
    INTEGER :: ACC2, ACC3, MACC2, MACC3
    !
    INTEGER :: I
    !
    !
    !
    !
!><MV TO SUBROUTINE    REAL(DP),ALLOCATABLE, DIMENSION(:,:)                 :: M
!><MV TO SUBROUTINE!
!><MV TO SUBROUTINE    REAL(DP),ALLOCATABLE, DIMENSION(:,:)                 :: CM
!><MV TO SUBROUTINE    INTEGER                                              :: INFO, LWORK
!><MV TO SUBROUTINE    INTEGER                                              :: IH
!><MV TO SUBROUTINE    REAL(DP),ALLOCATABLE, DIMENSION(:,:)                 :: SVDU, SVDVY, SVDIS
!><MV TO SUBROUTINE    REAL(DP),ALLOCATABLE, DIMENSION(:)                   :: SVDS
!><MV TO SUBROUTINE    REAL(DP),ALLOCATABLE, DIMENSION(:)                   :: WORK
!><MV TO SUBROUTINE    INTEGER                                              :: CNT
    !
    !
    !
    ! Set internal parameters:
    NB=MINVAL((/NX*NY,ENB/))
    SVDTOL=ESVDTOL
    LAMBDAP=ELAMBDAP
    IF (LAMBDAP.LT.1.0D0) LAMBDAP=1.0D0
!PRINT*, 'PSF:'
    CTP(:)=ECTP(:)
    NPY=SIZE(EPSF2D,1)
    NPX=SIZE(EPSF2D,2)
    CALL ALLOCATE_2D_DP(PSF,NPY,NPX,'Internal PSF')
    CALL ALLOCATE_2D_DP(IPSF,NPY,NPX,'Internal PSF')
    IPSF(:,:)=EPSF2D(:,:)
    PSF(:,:)=EPSF2D(:,:)
    PRINT*, SUM(PSF), SUM(IPSF), SUM(EPSF2D)
    PRINT*, CTP
    !PSF(CTP(1)-1:CTP(1)+1,CTP(2)-1:CTP(2)+1) &
    !    =EPSF2D(CTP(1)-1:CTP(1)+1,CTP(2)-1:CTP(2)+1)
    !
    ! Set for the first time INVSIZE and NTR
    ! 
    ! Npar times batch size
    INVSIZE=NB*NPAR
    ! Ntr:
    NTR = NX*NY/NB
    IF ( MOD(NX*NY, NB) .NE. 0) THEN
      NTR=NTR+1
    ENDIF
    !
!PRINT*, 'Allocating:'
    !
    CALL ALLOCATE_3D_DP(HINV,INVSIZE,INVSIZE,NTR &
        ,' Allocating inverse matrix')
    CALL ALLOCATE_3D_DP(HPREINV,INVSIZE,INVSIZE,NTR &
        ,' Allocating inverse matrix')
    !CALL ALLOCATE_3D_DP(WINV,INVSIZE,INVSIZE,NTR,' Allocating inverse matrix')
    CALL ALLOCATE_2D_DP(WINV,NY,NX,' Allocating inverse matrix')
    !
!PRINT*, 'Loops:'
!><!$OMP PARALLEL DEFAULT(PRIVATE), NUM_THREADS(41), &
!><!$OMP& SHARED(HINV,LAMBDAP,SVDTOL,INVSIZE,PSF,ALP,NPAR,NB,NY,NX,CTP,NPX,NPY,NTR,WINV)
!$OMP PARALLEL DEFAULT(PRIVATE), NUM_THREADS(41), &
!$OMP& FIRSTPRIVATE(LAMBDAP,SVDTOL,INVSIZE,NPAR,NB,NY,NX,CTP,NPX,NPY,NTR), &
!$OMP& SHARED(HINV,PSF,ALP,WINV,IPSF,HPREINV)
!$OMP DO SCHEDULE(DYNAMIC)
    DO ITN=1,NTR
!PRINT*, ITN, NTR
      !
      ITN2=(ITN-1)*NB
      ITN3=(ITN-1)*NB
      ! Real 2D identity
      DO IT3=1,NX
        DO IT2=1,NY
          ! The y coordinate of the matrix only depends on the pixel considered:
          ACC3 = (IT2-1) + (IT3-1) * NY
          IF ( (ACC3 .LT. ITN3) .OR. (ACC3+1 .GT. ITN3+INVSIZE/NPAR)) CYCLE
          ! For each x, y of the image, we consider the PSF effect:
          DO IP3=1,NPX
            DO IP2=1,NPY
              ! x,y identity of the pixel considered:
              ITP2 = (IT2-1) + ((IP2-1)-(CTP(1)-1))
              ITP3 = (IT3-1) + ((IP3-1)-(CTP(2)-1))
              ! Check if PSF (x,y) considered is actually out of area:
              IF ( (ITP2 .LT. 0) .OR. (ITP2+1 .GT. NY) ) CYCLE
              IF ( (ITP3 .LT. 0) .OR. (ITP3+1 .GT. NX) ) CYCLE
              ! Identify (x,y) in aprime:
              ACC2 = ITP2 + ITP3 * NY
              IF ( (ACC2 .LT. ITN2) .OR. (ACC2+1 .GT. ITN2+INVSIZE/NPAR) ) CYCLE
              !
              MACC2 = MOD(ACC2,NB) * NPAR
              MACC3 = MOD(ACC3,NB) * NPAR
              !
              HINV(MACC2+1:MACC2+NPAR,MACC3+1:MACC3+NPAR,ITN) &
                  = ALP(:,:,ITP2+1,ITP3+1) * IPSF(IP2,IP3)
              !WINV(MACC2+1:MACC2+NPAR,MACC3+1:MACC3+NPAR,ITN) &
              !    = PSF(IP2,IP3) * 1.0D0
              WINV(ITP2+1,ITP3+1)=WINV(ITP2+1,ITP3+1) &
                  +IPSF(IP2,IP3)
!PRINT*, PSF(IP2,IP3)
              !
            ENDDO !ITP2
          ENDDO !ITP3
          ! PSF END
        ENDDO !IT2
      ENDDO !IT3
!PRINT*, 'Modify lambdas:'
      !
      ! Modify diagonal with lambda:
      HPREINV(:,:,ITN)=HINV(I,I,ITN)
      DO I=1,INVSIZE
PRINT*, 'ITN: ', ITN, 'I: ', I, ' ; ', HINV(I,I,ITN)
        HINV(I,I,ITN) = HINV(I,I,ITN) + HINV(I,I,ITN) * LAMBDAP
      ENDDO
      !
      CALL INV_MAT(HINV(:,:,ITN))
!>< MV TO SUBROUTINE:!PRINT*, 'Invert matrix:'
!>< MV TO SUBROUTINE:    IF (SUM(ABS(HINV(:,:,ITN))).LT.1.0D-12) CYCLE
!>< MV TO SUBROUTINE:    ALLOCATE(SVDU(INVSIZE,INVSIZE))
!>< MV TO SUBROUTINE:    ALLOCATE(SVDVY(INVSIZE,INVSIZE))
!>< MV TO SUBROUTINE:    ALLOCATE(SVDIS(INVSIZE,INVSIZE))
!>< MV TO SUBROUTINE:    ALLOCATE(M(INVSIZE,INVSIZE))
!>< MV TO SUBROUTINE:    ALLOCATE(CM(INVSIZE,INVSIZE))
!>< MV TO SUBROUTINE:    ALLOCATE(SVDS(INVSIZE))
!>< MV TO SUBROUTINE:    SVDU(:,:)=0.D0
!>< MV TO SUBROUTINE:    SVDVY(:,:)=0.D0
!>< MV TO SUBROUTINE:    SVDIS(:,:)=0.D0
!>< MV TO SUBROUTINE:    CM(:,:)=0.D0
!>< MV TO SUBROUTINE:    M(:,:)=0.D0
!>< MV TO SUBROUTINE:    SVDS(:)=0.D0
!>< MV TO SUBROUTINE:    !
!>< MV TO SUBROUTINE:    ! SVD estimation of workspace
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(1):'
!>< MV TO SUBROUTINE:    ALLOCATE(WORK(1))
!>< MV TO SUBROUTINE:    WORK(:)=0
!>< MV TO SUBROUTINE:    LWORK=-1
!>< MV TO SUBROUTINE:    CM(:,:)=HINV(:,:,ITN)
!>< MV TO SUBROUTINE:    CALL DGESVD('a','a',INVSIZE,INVSIZE, CM, INVSIZE, SVDS, SVDU, INVSIZE, SVDVY, INVSIZE&
!>< MV TO SUBROUTINE:        , WORK,LWORK,INFO)
!>< MV TO SUBROUTINE:    LWORK=INT(WORK(1))
!>< MV TO SUBROUTINE:    DEALLOCATE(WORK)
!>< MV TO SUBROUTINE:    !
!>< MV TO SUBROUTINE:    ! SVD
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(2):'
!>< MV TO SUBROUTINE:    ALLOCATE(WORK(LWORK))
!>< MV TO SUBROUTINE:    WORK(:)=0
!>< MV TO SUBROUTINE:    CM(:,:)=HINV(:,:,ITN)
!>< MV TO SUBROUTINE:    CALL DGESVD('a','a',INVSIZE,INVSIZE, CM, INVSIZE, SVDS, SVDU, INVSIZE, SVDVY, INVSIZE&
!>< MV TO SUBROUTINE:        , WORK,LWORK,INFO)
!>< MV TO SUBROUTINE:    DEALLOCATE(WORK)
!>< MV TO SUBROUTINE:    !
!>< MV TO SUBROUTINE:!><IS IT RIGHT?PRINT*, ITN, SUM(ABS(SVDS)), SUM(ABS(SVDU)), SUM(ABS(SVDVY)) &
!>< MV TO SUBROUTINE:!><IS IT RIGHT?    , INVSIZE, SUM(ABS(CM)), SUM(ABS(HINV(:,:,ITN)))
!>< MV TO SUBROUTINE:PRINT*, SUM(ABS(HINV(:,:,ITN)))
!>< MV TO SUBROUTINE:!><NR><    CM(:,:)=HINV(:,:,ITN)
!>< MV TO SUBROUTINE:!><NR><    CALL SVDCMP2(CM,INVSIZE,INVSIZE,INVSIZE,INVSIZE,SVDS,SVDVY)
!>< MV TO SUBROUTINE:!><NR><    SVDU(:,:)=CM(:,:)
!>< MV TO SUBROUTINE:!><NR><    SVDVY(:,:)=TRANSPOSE(SVDVY)
!>< MV TO SUBROUTINE:!><NR><    CM(:,:)=HINV(:,:,ITN)
!>< MV TO SUBROUTINE:!><PRINT*, ITN, SUM(ABS(SVDS)), SUM(ABS(SVDU)), SUM(ABS(SVDVY)) &
!>< MV TO SUBROUTINE:!><    , INVSIZE, SUM(ABS(CM)), SUM(ABS(HINV(:,:,ITN)))
!>< MV TO SUBROUTINE:!><STOP
!>< MV TO SUBROUTINE:    ! Inverse matrix
!>< MV TO SUBROUTINE:    !
!>< MV TO SUBROUTINE:    ! Discard the less relevant elements
!>< MV TO SUBROUTINE:    CNT=0
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(3):'
!>< MV TO SUBROUTINE:    DO IH=1,INVSIZE
!>< MV TO SUBROUTINE:      IF (ABS(SVDS(IH)) .LT. (SVDTOL * MAXVAL(ABS(SVDS)))) THEN
!>< MV TO SUBROUTINE:        SVDIS(IH,IH)=0.D0
!>< MV TO SUBROUTINE:      ELSE
!>< MV TO SUBROUTINE:        CNT=CNT+1
!>< MV TO SUBROUTINE:        SVDIS(IH,IH)=1.D0/SVDS(IH)
!>< MV TO SUBROUTINE:      ENDIF
!>< MV TO SUBROUTINE:    ENDDO
!>< MV TO SUBROUTINE:    !
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(4):'
!>< MV TO SUBROUTINE:!><    M(:,:)=MATMUL(TRANSPOSE(SVDVY), MATMUL(SVDIS,TRANSPOSE(SVDU)))
!>< MV TO SUBROUTINE:!><PRINT*, 'SVD(5):'
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(4):'
!>< MV TO SUBROUTINE:    M(:,:)=TRANSPOSE(SVDU)
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(4.1):'
!>< MV TO SUBROUTINE:    M(:,:)=MATMUL(SVDIS,M(:,:))
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(4.2):'
!>< MV TO SUBROUTINE:    M(:,:)=MATMUL(TRANSPOSE(SVDVY),M(:,:))
!>< MV TO SUBROUTINE:!PRINT*, 'SVD(5):'
!>< MV TO SUBROUTINE:    HINV(:,:,ITN)=M(:,:)
!>< MV TO SUBROUTINE:PRINT*, ITN, SUM(ABS(HINV(:,:,ITN)))
!>< MV TO SUBROUTINE:!><PRINT*, ITN, SUM(ABS(SVDS)), SUM(ABS(SVDU)), SUM(ABS(SVDVY)), SUM(ABS(HINV(:,:,ITN))) &
!>< MV TO SUBROUTINE:!><    , LWORK, INFO, INVSIZE, SUM(ABS(CM)), SUM(ABS(M)), SUM(ABS(WINV(:,:,ITN)))
!>< MV TO SUBROUTINE:    DEALLOCATE(SVDU)
!>< MV TO SUBROUTINE:    DEALLOCATE(SVDVY)
!>< MV TO SUBROUTINE:    DEALLOCATE(SVDIS)
!>< MV TO SUBROUTINE:    DEALLOCATE(M)
!>< MV TO SUBROUTINE:    DEALLOCATE(CM)
!>< MV TO SUBROUTINE:    DEALLOCATE(SVDS)
!>< MV TO SUBROUTINE:!   

      !
    ENDDO ! ITN
!$OMP END DO
!$OMP END PARALLEL 
    !
    CALL WRITE_BIN('cou_hinv.bin',SIZE(SHAPE(HINV)), SHAPE(HINV) &
        , SIZE(HINV,KIND=8),REAL(HINV),3000,3)
    CALL WRITE_BIN('cou_hpreinv.bin',SIZE(SHAPE(HPREINV)), SHAPE(HPREINV) &
        , SIZE(HPREINV,KIND=8),REAL(HPREINV),3000,3)
!    CALL WRITE_BIN('cou_winv.bin',SIZE(SHAPE(WINV)), SHAPE(WINV) &
!        , SIZE(WINV,KIND=8),REAL(WINV),3000,3)
!STOP
    !STOP
    !
  END SUBROUTINE BUILD_HINV
  !
  SUBROUTINE INV_MAT(M)
    !
    REAL(DP), INTENT(INOUT), DIMENSION(INVSIZE,INVSIZE)  :: M
!
    REAL(DP), DIMENSION(INVSIZE,INVSIZE)                 :: CM
    INTEGER                                              :: INFO, LWORK
    INTEGER                                              :: I
    REAL(DP), DIMENSION(INVSIZE,INVSIZE)                 :: SVDU, SVDVY, SVDIS
    REAL(DP), DIMENSION(INVSIZE)                         :: SVDS
    REAL(DP), ALLOCATABLE, DIMENSION(:)                  :: WORK
    INTEGER                                              :: CNT
    !
!!PRINT*, 'SVD(0):'
    SVDU(:,:)=0.D0
    SVDVY(:,:)=0.D0
    SVDIS(:,:)=0.D0
    CM(:,:)=0.D0
    SVDS(:)=0.D0
    !
    ! SVD estimation of workspace
!PRINT*, 'SVD(1):'
    ALLOCATE(WORK(1))
    WORK(:)=0
    LWORK=-1
    CM(:,:)=M(:,:)
    CALL DGESVD('a','a',INVSIZE,INVSIZE, CM, INVSIZE, SVDS, SVDU, INVSIZE, SVDVY, INVSIZE&
        , WORK,LWORK,INFO)
    LWORK=INT(WORK(1))
    DEALLOCATE(WORK)
    !
    ! SVD
!PRINT*, 'SVD(2):'
    ALLOCATE(WORK(LWORK))
    WORK(:)=0
    CM(:,:)=M(:,:)
    CALL DGESVD('a','a',INVSIZE,INVSIZE, CM, INVSIZE, SVDS, SVDU, INVSIZE, SVDVY, INVSIZE&
        , WORK,LWORK,INFO)
    DEALLOCATE(WORK)
    ! Inverse matrix
    !
!><NR><    CM(:,:)=M(:,:)
!><NR><    CALL SVDCMP2(CM,INVSIZE,INVSIZE,INVSIZE,INVSIZE,SVDS,SVDVY)
!><NR><    SVDU(:,:)=CM(:,:)
!><NR><    SVDVY(:,:)=TRANSPOSE(SVDVY)
!><NR><    ! Discard the less relevant elements
    CNT=0
!PRINT*, 'SVD(3):'
    DO I=1,INVSIZE
      IF (ABS(SVDS(I)) .LT. (SVDTOL * MAXVAL(ABS(SVDS)))) THEN
        SVDIS(I,I)=0.D0
      ELSE
        CNT=CNT+1
        SVDIS(I,I)=1.D0/SVDS(I)
      ENDIF
    ENDDO
    !
!PRINT*, 'SVD(4):'
!><    M(:,:)=MATMUL(TRANSPOSE(SVDVY), MATMUL(SVDIS,TRANSPOSE(SVDU)))
!><PRINT*, 'SVD(5):'
!PRINT*, 'SVD(4):'
    M(:,:)=TRANSPOSE(SVDU)
!PRINT*, 'SVD(4.1):'
    M(:,:)=MATMUL(SVDIS,M(:,:))
!PRINT*, 'SVD(4.2):'
    M(:,:)=MATMUL(TRANSPOSE(SVDVY),M(:,:))
!PRINT*, 'SVD(5):'

    !
  END SUBROUTINE INV_MAT

  !
  SUBROUTINE ALPHA_DOT_DELTA(SH3,IVEC,OVEC)
    !
    ! Delta:
    INTEGER,DIMENSION(3),INTENT(IN) :: SH3
    REAL(DP), INTENT(IN), DIMENSION(SH3(1),SH3(2),SH3(3))  :: IVEC
    !
    REAL(DP), INTENT(INOUT), DIMENSION(SH3(1),SH3(2),SH3(3))  :: OVEC
    !
    REAL(DP), DIMENSION(NPAR,NPAR)       :: MAT
    REAL(DP), DIMENSION(1,NPAR)          :: PROD1, PROD2
    REAL(DP), DIMENSION(1,1)             :: PROD
    INTEGER :: IT2, IT3, IP2, IP3
    INTEGER :: ITP2, ITP3
    !
    INTEGER :: I
    !
    !
    OVEC(:,:,:)=0.0D0
    ! Real 2D identity
!$OMP PARALLEL DEFAULT(PRIVATE), NUM_THREADS(41), &
!$OMP& SHARED(OVEC,IVEC,NX,NY,NPX,NPY,CTP,ALP,NPAR,LAMBDAP,PSF)
!$OMP DO SCHEDULE(GUIDED)
    DO IT3=1,NX
      DO IT2=1,NY
        ! For each x, y of the image, we consider the PSF effect:
        ! Test smaller fwd: DO IP3=1,NPX
        ! Test smaller fwd:   DO IP2=1,NPY
        DO IP3=NPX/2+1-5,NPX/2+1+5
          DO IP2=NPY/2+1-5,NPY/2+1+5
            ! x,y identity of the pixel considered:
            ITP2 = (IT2-1) + ((IP2-1)-(CTP(1)-1))
            ITP3 = (IT3-1) + ((IP3-1)-(CTP(2)-1))
            ! Check if PSF (x,y) considered is actually out of area:
            IF ( (ITP2 .LT. 0) .OR. (ITP2+1 .GT. NY) ) CYCLE
            IF ( (ITP3 .LT. 0) .OR. (ITP3+1 .GT. NX) ) CYCLE
            !
            MAT(:,:) = alp(:,:,ITP2+1,ITP3+1)
            !
            IF ( (IT2-1 .EQ. ITP2) .AND. (IT3-1 .EQ. ITP3) ) THEN
              ! Modify diagonal with lambda:
              DO I=1,NPAR
                MAT(I,I) = MAT(I,I) + MAT(I,I) * LAMBDAP
              ENDDO
              !
            ENDIF
      !
            DO I=1,NPAR
              PROD1(1,:)=MAT(:,I) * PSF(IP2,IP3)
              PROD2(1,:)=IVEC(:,ITP2+1,ITP3+1)
              PROD=MATMUL(PROD1, TRANSPOSE(PROD2))
              OVEC(I,IT2,IT3) = OVEC(I,IT2,IT3) + PROD(1,1)
            ENDDO
            !
          ENDDO !ITP2
        ENDDO !ITP3
        ! PSF END
      ENDDO !IT2
    ENDDO !IT3
!$OMP END DO
!$OMP END PARALLEL 
    !
    !
  END SUBROUTINE ALPHA_DOT_DELTA
  !
!>< OBS ><  !
!>< OBS ><  SUBROUTINE HINV_DOT_BETA(SHV,IVEC,OVEC)
!>< OBS ><    !
!>< OBS ><    INTEGER,DIMENSION(3),INTENT(IN) :: SHV
!>< OBS ><    REAL(DP), INTENT(IN), DIMENSION(SHV(1),SHV(2),SHV(3))  :: IVEC
!>< OBS ><    !
!>< OBS ><    REAL(DP), INTENT(INOUT), DIMENSION(SHV(1),SHV(2),SHV(3))  :: OVEC
!>< OBS ><    !
!>< OBS ><    INTEGER :: ITN
!>< OBS ><    INTEGER :: ITPX,ITPY,ITIPX,ITIPY
!>< OBS ><    INTEGER :: CNT, I
!>< OBS ><    !
!>< OBS ><    REAL(DP), DIMENSION(INVSIZE)  :: RES
!>< OBS ><    REAL(DP), DIMENSION(INVSIZE,1)  :: VEC
!>< OBS ><    REAL(DP), DIMENSION(1,1)  :: PROD
!>< OBS ><    !
!>< OBS ><    ! Remember some dimensions:
!>< OBS ><    !NB=SH3(1)/SHV(1)
!>< OBS ><    !NTR=(SHV(3)*SHV(2))/NB
!>< OBS ><    !IF ( MOD(SHV(3)*SHV(2), NB) .NE. 0) THEN
!>< OBS ><    !  NTR=NTR+1
!>< OBS ><    !ENDIF
!>< OBS ><    !
!>< OBS ><    ! Initialize some counters:
!>< OBS ><    ITPX=-1
!>< OBS ><    ITPY=-1
!>< OBS ><    ITIPX=-1
!>< OBS ><    ITIPY=-1
!>< OBS ><    !
!>< OBS ><    ! Proceed with the product:
!>< OBS ><    DO ITN=1,NTR
!>< OBS ><      !
!>< OBS ><      VEC(:,:)=0.0D0
!>< OBS ><      RES(:)=0.0D0
!>< OBS ><      !
!>< OBS ><      ! Fill VEC:
!>< OBS ><      CNT=0
!>< OBS ><      DO WHILE (CNT.LT.(INVSIZE/NPAR) )
!>< OBS ><        ITPY=ITPY+1
!>< OBS ><        IF ( MOD(ITPY,SHV(2)) .EQ. 0) THEN
!>< OBS ><          ITPX=ITPX+1
!>< OBS ><          ITPY=0
!>< OBS ><        ENDIF
!>< OBS ><        IF ( ((ITPY+1).GT.SHV(2)).OR.((ITPX+1) .GT. SHV(3)) ) EXIT
!>< OBS ><        !
!>< OBS ><        ! (NPAR encodes the number of free parameters)
!>< OBS ><        VEC(CNT*NPAR+1:(CNT+1)*NPAR,1)=IVEC(:,ITPY+1,ITPX+1)
!>< OBS ><        CNT=CNT+1
!>< OBS ><      ENDDO ! End building VEC
!>< OBS ><      !
!>< OBS ><      ! Make matrix product:
!>< OBS ><      DO I=1,CNT*NPAR
!>< OBS ><        !PRINT*, SHAPE(HINV(1:CNT*NPAR,I:I,ITN))
!>< OBS ><        !PRINT*, SHAPE(VEC(1:CNT*NPAR,:))
!>< OBS ><        PROD=MATMUL(TRANSPOSE(HINV(1:CNT*NPAR,I:I,ITN)), VEC(1:CNT*NPAR,:))
!>< OBS ><        RES(I)=PROD(1,1)
!>< OBS ><      ENDDO
!>< OBS ><      !
!>< OBS ><      ! Fill OVEC:
!>< OBS ><      CNT=0
!>< OBS ><      DO WHILE (CNT.LT.(INVSIZE/NPAR) )
!>< OBS ><        ITIPY=ITIPY+1
!>< OBS ><        IF ( MOD(ITIPY,SHV(2)) .EQ. 0) THEN
!>< OBS ><          ITIPX=ITIPX+1
!>< OBS ><          ITIPY=0
!>< OBS ><        ENDIF
!>< OBS ><        IF ( ((ITIPY+1).GT.SHV(2)).OR.((ITIPX+1) .GT. SHV(3)) ) EXIT
!>< OBS ><        !
!>< OBS ><        ! (NPAR encodes the number of free parameters)
!>< OBS ><        OVEC(:,ITIPY+1,ITIPX+1)=RES(CNT*NPAR+1:(CNT+1)*NPAR)
!>< OBS ><        CNT=CNT+1
!>< OBS ><      ENDDO ! End building VEC
!>< OBS ><      !
!>< OBS ><    ENDDO ! Number of inverse matrixes
!>< OBS ><    !
!>< OBS ><  END SUBROUTINE HINV_DOT_BETA
!>< OBS ><  !
  !
  SUBROUTINE HINV_DOT_BETA(SHV,IVEC,OVEC)
    !
    INTEGER,DIMENSION(3),INTENT(IN) :: SHV
    REAL(DP), INTENT(IN), DIMENSION(SHV(1),SHV(2),SHV(3))  :: IVEC
    !
    REAL(DP), INTENT(INOUT), DIMENSION(SHV(1),SHV(2),SHV(3))  :: OVEC
    !
    INTEGER :: ITN
    INTEGER :: ITPX,ITPY
    INTEGER :: CNT, I
    !
    REAL(DP), DIMENSION(INVSIZE)  :: RES
    REAL(DP), DIMENSION(INVSIZE,1)  :: VEC
    REAL(DP), DIMENSION(1,1)  :: PROD
    !
    LOGICAL :: SIB
    !
    OVEC(:,:,:)=0.0D0
    !GCNT=0
    !
    ! Proceed with the product:
!$OMP PARALLEL DEFAULT(PRIVATE), NUM_THREADS(41), &
!$OMP& SHARED(NTR,NX,NY,INVSIZE,NPAR,IVEC,OVEC,HINV,WINV)
!$OMP DO SCHEDULE(GUIDED)
    DO ITN=1,NTR
      !
      VEC(:,:)=0.0D0
      RES(:)=0.0D0
      !
      ! Fill VEC:
      CNT=0
      SIB=.FALSE.
!PRINT*,ITN,NTR,NX,NY,CNT,SIB,SUM(ABS(VEC)),SUM(ABS(RES)),NPAR,INVSIZE&
    !,SUM(ABS(IVEC)),SUM(ABS(OVEC)), SUM(ABS(HINV)), SUM(ABS(WINV))
      DO ITPX=1,NX
        DO ITPY=1,NY
          !
          IF ((((ITPX-1)*NY+ITPY)*NPAR).LE.(INVSIZE*(ITN-1))) CYCLE
          !
          VEC(CNT*NPAR+1:(CNT+1)*NPAR,1)=IVEC(:,ITPY,ITPX)
!PRINT*, ITPY, ITPX
          CNT=CNT+1
          !GCNT=GCNT+1
          !
          IF ((CNT*NPAR+1).GT.INVSIZE) SIB=.TRUE.
          IF (SIB.EQV..TRUE.) EXIT
          !
        ENDDO ! End NY
        IF (SIB.EQV..TRUE.) EXIT
      ENDDO ! End NX
      ! End building VEC
      !
      ! Make matrix product:
      DO I=1,CNT*NPAR
        !PROD=MATMUL(TRANSPOSE(HINV(1:CNT*NPAR,I:I,ITN)), VEC(1:CNT*NPAR,:))
        PROD=MATMUL(TRANSPOSE(&
            HINV(1:CNT*NPAR,I:I,ITN)), VEC(1:CNT*NPAR,:))
            !HINV(1:CNT*NPAR+1,I:I,ITN)*WINV(1:CNT*NPAR+1,I:I,ITN)), VEC(1:CNT*NPAR+1,:))
        RES(I)=PROD(1,1)!*SUM(WINV(1:CNT*NPAR,I,ITN)/NPAR)!**2
      ENDDO
      !
      ! Fill OVEC:
      CNT=0
      SIB=.FALSE.
      DO ITPX=1,NX
        DO ITPY=1,NY
          !
          IF ((((ITPX-1)*NY+ITPY)*NPAR).LE.(INVSIZE*(ITN-1))) CYCLE
          !
          OVEC(:,ITPY,ITPX)=RES(CNT*NPAR+1:(CNT+1)*NPAR)
          CNT=CNT+1
          !
          IF ((CNT*NPAR+1).GT.INVSIZE) SIB=.TRUE.
          IF (SIB.EQV..TRUE.) EXIT
          !
        ENDDO ! End NY
        IF (SIB.EQV..TRUE.) EXIT
      ENDDO ! End NX
!PRINT*,ITN,NTR,NX,NY,CNT,SIB,SUM(ABS(VEC)),SUM(ABS(RES)),NPAR,INVSIZE&
    !,SUM(ABS(IVEC)),SUM(ABS(OVEC)), SUM(ABS(HINV)), SUM(ABS(WINV))
!PRINT*, ''
      ! End fill OVEC:
      !
    ENDDO ! Number of inverse matrixes
!$OMP END DO
!$OMP END PARALLEL 
!PRINT*, GCNT, 'A'
    !
  END SUBROUTINE HINV_DOT_BETA
  !
  SUBROUTINE CALCULATE_DELTA(DELTA,FAILED)
    !
    REAL(DP), INTENT(INOUT),ALLOCATABLE, DIMENSION(:,:,:)  :: DELTA
    LOGICAL, INTENT(INOUT)   :: FAILED
    !
    REAL(DP), DIMENSION(NPAR,NY,NX)  :: PBETAI
    REAL(DP), DIMENSION(NPAR,NY,NX)  :: BETAI
    REAL(DP), DIMENSION(NPAR,NY,NX)  :: DELTAI
    REAL(DP), DIMENSION(NPAR,NY,NX)  :: DIFFI
    REAL(DP), DIMENSION(NPAR,NY,NX)  :: CORRI
    !
    INTEGER :: CNT, I
    REAL(DP) :: SQRCORR, SQRVAR, PSQRCORR, DSQRCORR
    !
    ! Problem to be solved:
    !
    ! alpha (*) delta = beta
    !
    ! Knowns: alpha, beta
    ! Unknown: delta
    !
    ! We would like to calculate:
    ! h=(alpha)^-1
    ! Not possible, instead we use a block-diagonal approach of h:
    ! hinv~~h
    !
    ! Method: (M. van Noort 2012)
    !
    ! delta_i = hinv (*) beta
    !
    ! Iterate:
    !
    !    1- alpha (*) delta_i = beta_i
    !    2- diff_i = beta - beta_i
    !    3- corr_i = hinv (*) diff_i
    !    4- delta_i+1 = delta_i +/- corr_i
    !
PRINT*, ' '
PRINT*, ' '
PRINT*, ' '
PRINT*, ' CALCULATE_DELTA: '
PRINT*, 'BETA: ', MAXVAL(BETA), SUM(BETA)
PRINT*, 'HINV: ', MAXVAL(HINV), SUM(HINV)
PRINT*, 'SVDTOL: ', SVDTOL
PRINT*, 'LAMBDAP: ', LAMBDAP
    IF (ALLOCATED(DELTA)) DEALLOCATE(DELTA)
    CALL ALLOCATE_3D_DP(DELTA,NPAR,NY,NX, 'DELTA')
    !
    CALL HINV_DOT_BETA(SHAPE(BETA),BETA,DELTAI)
    CALL WRITE_BIN('deltai.bin',SIZE(SHAPE(DELTAI)), SHAPE(DELTAI) &
        , SIZE(DELTAI,KIND=8),REAL(DELTAI),3000,3)
    !
    CNT=0
    PBETAI(:,:,:)=1.0D29
    SQRCORR=1.0D0
    SQRVAR=1.0D0
    PSQRCORR=1.0D29
    DSQRCORR=1.0D0
    DO WHILE ( (CNT.LT.1000) .AND. (SQRCORR.GT.1.0D-6) &
        .AND. (DSQRCORR.GT.0.0D0) )
        !.AND. (SQRVAR.GT.1.0D-16) .AND. (DSQRCORR.GT.0.0D0) )
      !
      ! step 1:
      CALL ALPHA_DOT_DELTA(SHAPE(DELTAI),DELTAI,BETAI)
      !
      ! step 2:
      DIFFI(:,:,:)=BETA(:,:,:)-BETAI(:,:,:)
      !
      ! step 3:
      CALL HINV_DOT_BETA(SHAPE(DIFFI),DIFFI,CORRI)
      !
      ! step 4:
      DO I=1,NPAR
        CORRI(I,:,:)=CORRI(I,:,:)*WINV(:,:)
      ENDDO
      DELTAI(:,:,:)=DELTAI(:,:,:)+CORRI(:,:,:)
      !
      ! Update conditional clauses:
      CNT=CNT+1
      SQRCORR=SUM((DIFFI/(ABS(BETA)+1.0D-3))**2)*SIZEFAC
      SQRVAR=SUM(((BETAI-PBETAI)/(ABS(BETA)+1.0D-3))**2)*SIZEFAC
      DSQRCORR=PSQRCORR-SQRCORR
      !
      IF (MOD(CNT,100).EQ.0) THEN
        PRINT*, 'IT: ', CNT, SQRCORR, SQRVAR, DSQRCORR
      ENDIF
      PBETAI(:,:,:)=BETAI(:,:,:)
      PSQRCORR=SQRCORR
      !
    ENDDO
    PRINT*, CNT, SQRCORR, SQRVAR, DSQRCORR
    !
    DELTA(:,:,:)=0.0D0
    FAILED=.TRUE.
    IF (SQRCORR.LT.1.0D-6) THEN
      DELTA(:,:,:)=DELTAI(:,:,:)
      FAILED=.FALSE.
    ENDIF
!    IF (SQRVAR.LT.1.0D-16) THEN
!      DELTA(:,:,:)=DELTAI(:,:,:)/4.0D0
!      FAILED=.FALSE.
!    ENDIF
    !
!><    CALL WRITE_BIN('beta.bin',SIZE(SHAPE(BETA)), SHAPE(BETA) &
!><        , SIZE(BETA,KIND=8),REAL(BETA),3000,3)
!><    CALL WRITE_BIN('betai.bin',SIZE(SHAPE(BETAI)), SHAPE(BETAI) &
!><        , SIZE(BETAI,KIND=8),REAL(BETAI),3000,3)
    CALL WRITE_BIN('delta.bin',SIZE(SHAPE(DELTA)), SHAPE(DELTA) &
        , SIZE(DELTA,KIND=8),REAL(DELTA),3000,3)
    ! Remove all the ALLOCATED variables:
    CALL DEALLOCATE_ALL()
PRINT*, ' '
PRINT*, ' '
PRINT*, ' '
    !
  END SUBROUTINE CALCULATE_DELTA
  !
  SUBROUTINE DEALLOCATE_ALL()
    !
    IF (ALLOCATED(ALP)) DEALLOCATE(ALP)
    IF (ALLOCATED(BETA)) DEALLOCATE(BETA)
    IF (ALLOCATED(HINV)) DEALLOCATE(HINV)
    IF (ALLOCATED(WINV)) DEALLOCATE(WINV)
    IF (ALLOCATED(PSF)) DEALLOCATE(PSF)
    IF (ALLOCATED(IPSF)) DEALLOCATE(IPSF)
    !IF (ALLOCATED(TTT)) DEALLOCATE(TTT)
    !
  END SUBROUTINE DEALLOCATE_ALL
  !
END MODULE coupled_matrix_inversion
!
