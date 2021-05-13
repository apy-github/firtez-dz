!
MODULE LM
  !
  !================================================
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  !
  PUBLIC PGET_PERTURBATION
  PUBLIC GET_PERTURBATION
  PUBLIC HESSIAN
  PUBLIC GET_INVERSE_COV
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! hessian
  ! get_term1
  ! get_inverse
  ! get_perturbation
  ! get_inverse_cov
  !
  !------------------------------------------------
  !
  SUBROUTINE HESSIAN(NP, NF, JACOBIAN, HSS)
    !
    INTEGER                                           :: NP, NF
    !
    REAL*8, INTENT(IN), DIMENSION(:,:)      :: JACOBIAN(NP,NF)
    REAL*8, INTENT(INOUT), DIMENSION(:,:)   :: HSS(NP,NP)
    !
    HSS(:,:)=MATMUL(JACOBIAN, TRANSPOSE(JACOBIAN))
    !
  END SUBROUTINE HESSIAN
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_TERM1(NP, HSS, LAMBDAP, TERM1)
    !
    INTEGER                                         :: NP, I
    !
    REAL*8                                :: LAMBDAP
    !
    REAL*8, INTENT(IN), DIMENSION(:,:)    :: HSS(NP,NP)
    REAL*8, INTENT(INOUT), DIMENSION(:,:) :: TERM1(NP,NP)
    !
    TERM1(:,:)=HSS(:,:)
    DO I=1,NP
      TERM1(I,I)=TERM1(I,I)+HSS(I,I)*lambdaP
    ENDDO
    !
  END SUBROUTINE GET_TERM1
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_INVERSE(NP, M, IM, SVDTOL)
    !
    INTEGER      :: NP, I
    !
    REAL*8, INTENT(IN), DIMENSION(:,:)  :: M(NP,NP)
    REAL*8, INTENT(INOUT), DIMENSION(:,:)  :: IM(NP,NP)
    REAL*8, INTENT(IN), OPTIONAL          :: SVDTOL
    !
    REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: CM
    INTEGER                               :: INFO, LWORK
    REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SVDU, SVDVY, SVDIS
    REAL*8, ALLOCATABLE, DIMENSION(:)   :: SVDS, WORK
    INTEGER                               :: ISVDTOL
    INTEGER                                       :: CNT
    REAL*8 :: SVDTHRES
    !
    ALLOCATE(SVDU(NP,NP))
    SVDU(:,:)=0.D0
    ALLOCATE(SVDVY(NP,NP))
    SVDVY(:,:)=0.D0
    ALLOCATE(SVDIS(NP,NP))
    SVDIS(:,:)=0.D0
    ALLOCATE(CM(NP,NP))
    CM(:,:)=0.D0
    ALLOCATE(SVDS(NP))
    SVDS(:)=0.D0
    !
    ! SVD estimation of workspace
    ALLOCATE(WORK(1))
    WORK(:)=0
    LWORK=-1
    CM(:,:)=M(:,:)
    CALL DGESVD('a','a',NP,NP, CM, NP, SVDS, SVDU, NP, SVDVY, NP&
        , WORK,LWORK,INFO)
    LWORK=INT(WORK(1))
    DEALLOCATE(WORK)
    !
    ! SVD
    ALLOCATE(WORK(LWORK))
    WORK(:)=0
    CM(:,:)=M(:,:)
    CALL DGESVD('a','a',NP,NP, CM, NP, SVDS, SVDU, NP, SVDVY, NP&
        , WORK,LWORK,INFO)
    DEALLOCATE(WORK)
    ! Inverse matrix
    !
    ! 1- Discard the less relevant elements:
    SVDTHRES = MAXVAL(ABS(SVDS))
    IF (PRESENT(SVDTOL)) THEN
      SVDTHRES = SVDTHRES * SVDTOL
    ELSE
      ! Default
      SVDTHRES = SVDTHRES * 10.D-4
    END IF
    !
    CNT=0
    DO I=1,NP
      IF (ABS(SVDS(I)) .LT. SVDTHRES) THEN
        SVDIS(I,I)=0.D0
      ELSE
        CNT=CNT+1
        SVDIS(I,I)=1.D0/SVDS(I)
      ENDIF
    ENDDO
    !
    IM(:,:)=MATMUL(TRANSPOSE(SVDVY), MATMUL(SVDIS,TRANSPOSE(SVDU)))
    !
    DEALLOCATE(SVDU)
    DEALLOCATE(SVDVY)
    DEALLOCATE(SVDS)
    DEALLOCATE(SVDIS)
    !
  END SUBROUTINE GET_INVERSE
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_PERTURBATION(NP, NF, IJACOBIAN, DIFF, LAMBDAP, DELTA, REG_PEN, REG_HSS, DTYPE, SVDTOL)
    !
    INTEGER, INTENT(IN) :: NP, NF
    !
    REAL*8, INTENT(IN), DIMENSION(:,:)      :: IJACOBIAN(NF,NP)
    REAL*8, INTENT(IN)                      :: LAMBDAP
    REAL*8, INTENT(IN), DIMENSION(:)        :: DIFF(NF)
    !
    REAL*8, INTENT(IN), DIMENSION(:)        :: REG_PEN(NP)
    REAL*8, INTENT(IN), DIMENSION(:,:)      :: REG_HSS(NP,NP)
    !
    REAL*8, INTENT(INOUT), DIMENSION(:)     :: DELTA(NP)
    INTEGER, INTENT(IN)                     :: DTYPE
    REAL*8, INTENT(IN), OPTIONAL            :: SVDTOL
    !
    REAL*8, DIMENSION(:,:)                  :: JACOBIAN(NP,NF)
    REAL*8, DIMENSION(:,:)                  :: HSS(NP,NP)
    REAL*8, DIMENSION(:,:)                  :: TERM1(NP,NP)
    REAL*8, DIMENSION(:,:)                  :: ITERM1(NP,NP)
    REAL*8, DIMENSION(:,:)                  :: TERM2B(NP,NF)
    REAL*8, DIMENSION(:,:)                  :: TERM2(NP)
    !
    INTEGER :: I
    !
    JACOBIAN(:,:)=TRANSPOSE(IJACOBIAN(:,:))
    ! Get hessian:
    CALL HESSIAN(NP, NF, JACOBIAN, HSS)
    !
    ! Get Term 2:
    TERM2(:)=MATMUL(JACOBIAN(:,:),DIFF(:))
    !
    ! If reg. subtract reg term:
    TERM2(:)=TERM2(:)-REG_PEN(:)
    !
    ! If reg. add both hessians
    HSS(:,:)=HSS(:,:)+REG_HSS(:,:)
    !
    !
    ! Get Term 1:
    !
    CALL GET_TERM1(NP, HSS, LAMBDAP, TERM1)
    !
    ! Solve linear system:
    SELECT CASE (DTYPE)
      CASE(2)
        CALL NEW_CALCULATE_DELTA_UBICGSTAB(NP,DELTA,TERM1,TERM2)
      CASE(1)
        !
        ! Get Term 1 inverse:
        IF (PRESENT(SVDTOL)) THEN
          CALL GET_INVERSE(NP, TERM1, ITERM1, SVDTOL)
        ELSE
          CALL GET_INVERSE(NP, TERM1, ITERM1)
        ENDIF
        !
        ! Get Delta:
        !
        DELTA(:)=MATMUL(ITERM1(:,:), TERM2(:))
      CASE DEFAULT
        PRINT*, ""
        PRINT*, ""
        PRINT*, " ERROR: LM does not support DTYPE=", DTYPE
        PRINT*, ""
        PRINT*, ""
        STOP
    ENDSELECT
    !
  END SUBROUTINE GET_PERTURBATION
  !
  SUBROUTINE GET_INVERSE_COV(NP, M, IM_COV)
    !
    REAL*8, INTENT(IN), DIMENSION(:,:)    :: M(NP,NP)
    REAL*8, INTENT(INOUT), DIMENSION(:,:) :: IM_COV(NP,NP)
    !
    INTEGER                                         :: NP, I
    !
    REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: CM
    INTEGER                                         :: INFO, LWORK
    REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: SVDU, SVDVY, SVDIS2
    REAL*8, ALLOCATABLE, DIMENSION(:)     :: SVDS, WORK
    !
    ALLOCATE(SVDU(NP,NP))
    SVDU(:,:)=0.D0
    ALLOCATE(SVDVY(NP,NP))
    SVDVY(:,:)=0.D0
    ALLOCATE(SVDIS2(NP,NP))
    SVDIS2(:,:)=0.D0
    ALLOCATE(CM(NP,NP))
    CM(:,:)=0.D0
    ALLOCATE(SVDS(NP))
    SVDS(:)=0.D0
    !
    ! SVD estimation of workspace
    ALLOCATE(WORK(1))
    WORK(:)=0
    LWORK=-1
    CM(:,:)=M(:,:)
    CALL DGESVD('a','a',NP,NP, CM, NP, SVDS, SVDU, NP, SVDVY, NP&
        , WORK,LWORK,INFO)
    LWORK=INT(WORK(1))
    DEALLOCATE(WORK)
    !
    ! SVD
    ALLOCATE(WORK(LWORK))
    WORK(:)=0
    CM(:,:)=M(:,:)
    CALL DGESVD('a','a',NP,NP, CM, NP, SVDS, SVDU, NP, SVDVY, NP&
        , WORK,LWORK,INFO)
    DEALLOCATE(WORK)
    ! Inverse matrix covariance:
    !
    ! Build the diagonal matrix:
    DO I=1,NP
      SVDIS2(I,I)=1.0D0/SVDS(I)!/SVDS(I)
      IF (SVDS(I).LT.1.D-24) SVDIS2(I,I)=1.D24
      !PRINT*, SVDIS2(I,I)
    ENDDO
    !><DO I=1,NP
    !><  SVDIS2(I,I)=1.D0/SVDS(I)
    !><  IF (SVDS(I).LT.1.D-12) SVDIS2(I,I)=1.D12
    !><  PRINT*, SVDIS2(I,I)
    !><ENDDO
    !
    ! CALCULATE COVARIANCE MATRIX
!TEST      IM_COV(:,:)=MATMUL(SVDU,MATMUL(SVDIS2,TRANSPOSE(SVDU)))/(DBLE(NP)-1.D0)
    !IM_COV(:,:)=MATMUL(TRANSPOSE(SVDVY), MATMUL(SVDIS2,TRANSPOSE(SVDU)))
    IM_COV(:,:)=MATMUL(SVDU, MATMUL(SVDIS2,TRANSPOSE(SVDU)))
    !
    DEALLOCATE(SVDU)
    DEALLOCATE(SVDVY)
    DEALLOCATE(SVDS)
    DEALLOCATE(SVDIS2)
    !
  END SUBROUTINE GET_INVERSE_COV
  !
  !------------------------------------------------
  !
  SUBROUTINE PGET_PERTURBATION(NP, NF, IJACOBIAN &
      , DIFF, LAMBDAP, DELT, REG_PEN, REG_HSS, DTYPE, SVDTOL)
    !
    INTEGER                                           :: NP, NF
    !
    REAL*8, INTENT(IN), DIMENSION(:,:)      :: IJACOBIAN(NF,NP)
    REAL*8, INTENT(IN)                      :: LAMBDAP
    REAL*8, INTENT(IN), DIMENSION(:)        :: DIFF(NF)
    !
    REAL*8, INTENT(INOUT), DIMENSION(:)     :: DELT(NP)
    !
    REAL*8, INTENT(IN), DIMENSION(:)        :: REG_PEN(NP)
    REAL*8, INTENT(IN), DIMENSION(:,:)      :: REG_HSS(NP,NP)
    !
    INTEGER, INTENT(IN)                     :: DTYPE
    REAL*8, INTENT(IN), OPTIONAL                     :: SVDTOL
    !
    REAL*8, DIMENSION(:,:), ALLOCATABLE     :: JACOBIAN2
    REAL*8, DIMENSION(:,:), ALLOCATABLE     :: LESSIAN2
    !
    REAL*8, DIMENSION(:), ALLOCATABLE       :: RV2
    REAL*8, DIMENSION(:), ALLOCATABLE       :: DELT2
    INTEGER, DIMENSION(:)     :: NZE(NP)
    INTEGER :: I, J, NPNZ, NPNZ2
    !
    ! First, count number of non-zero elements:
    NZE(:)=0
    DO I=1,NP
      IF (SUM(ABS(IJACOBIAN(:,I))).GT.1.0D-19) THEN
        NZE(I)=1
      ENDIF
    ENDDO
    !
    NPNZ=SUM(NZE)
    !
    ALLOCATE(JACOBIAN2(NF,NPNZ))
    ALLOCATE(DELT2(NPNZ))
    !
    JACOBIAN2(:,:)=0.0D0
    DELT2(:)=0.0D0
    !
    ALLOCATE(LESSIAN2(NPNZ,NPNZ))
    ALLOCATE(RV2(NPNZ))
    !
    LESSIAN2(:,:)=0.0D0
    RV2(:)=0.0D0
    !
    NPNZ=0
    DO I=1,NP
      IF (NZE(I).GT.0) THEN
        NPNZ=NPNZ+1
        JACOBIAN2(:,NPNZ)=IJACOBIAN(:,I)
      ENDIF
    ENDDO
    !
PRINT*, ' [', 'Using ', NPNZ, ' out of ', NP,']'
    !
    NPNZ=0
    DO I=1,NP
      IF (NZE(I).GT.0) THEN
        NPNZ=NPNZ+1
        RV2(NPNZ)=REG_PEN(I)
      ENDIF
    ENDDO
    NPNZ=0
    DO I=1,NP
      IF (NZE(I).GT.0) THEN
        NPNZ=NPNZ+1
        NPNZ2=0
        DO J=1,NP
          IF (NZE(J).GT.0) THEN
            NPNZ2=NPNZ2+1
            LESSIAN2(NPNZ,NPNZ2)=REG_HSS(I,J)
          ENDIF
        ENDDO
      ENDIF
    ENDDO


print*, shape(lessian2), ';', shape(rv2), ';', shape(delt2)
!stop
    CALL GET_PERTURBATION(NPNZ, NF, JACOBIAN2, DIFF &
        , LAMBDAP, DELT2, RV2, LESSIAN2, DTYPE, SVDTOL)
    !
    NPNZ=0
    DO I=1,NP
      IF (NZE(I).GT.0) THEN
        NPNZ=NPNZ+1
        DELT(I)=DELT2(NPNZ)
      ENDIF
    ENDDO
    !
    DEALLOCATE(JACOBIAN2)
    DEALLOCATE(DELT2)
    DEALLOCATE(LESSIAN2)
    DEALLOCATE(RV2)
    !
  END SUBROUTINE PGET_PERTURBATION
  !
  !------------------------------------------------
  !
  SUBROUTINE NEW_CALCULATE_DELTA_UBICGSTAB(NPAR,DELT,ALP,BET)
    !
    ! Unpreconditioned BiCGSTAB
    !   Implemented following NPAC Technical Report SCCS 691:
    !      Conjugate Gradient Algorithms in Fortran 90 and High...
    !      ...Performance Fortran
    !      K.A.Hawick, K.Dincer, G.Robinson, G.C.Fox
    !
    INTEGER, INTENT(IN)   :: NPAR
    REAL*8, INTENT(INOUT), DIMENSION(NPAR)  :: DELT
    REAL*8, INTENT(IN), DIMENSION(NPAR,NPAR)  :: ALP
    REAL*8, INTENT(IN), DIMENSION(NPAR)  :: BET
    !
    REAL*8, DIMENSION(NPAR) :: XV, PV, RV, QV, CHECK!, XV0
    REAL*8 :: CRHO, CALPHA, CBETA, CRHO0
    LOGICAL :: FAILED
    !
    INTEGER :: ITER
    REAL*8 :: CONDIT
    REAL*8 :: DCONDIT
    REAL*8 :: PCONDIT
    !
    DELT(:)=0.0D0
    !
    XV(:) = 0.0D0
    !
    PV(:) = BET(:)
    RV(:) = MATMUL(ALP(:,:),XV(:))
    RV(:) = BET(:)-RV(:)
    !
    CRHO=SUM(RV(:)*RV(:))
    !
    QV(:)=MATMUL(ALP(:,:),PV(:))
    !
    CALPHA=SUM(PV(:)*QV(:))
    CALPHA = CRHO / CALPHA
    XV = XV + CALPHA * PV
    RV = RV - CALPHA * QV
    !
    PCONDIT=1.0D29
    DCONDIT=-1.0D29
    !
    DO ITER=1,100000
      ! Update:
      CRHO0 = CRHO
      !
      CRHO=SUM(RV(:)*RV(:))
      CBETA = CRHO / CRHO0
      !
      PV = RV + CBETA * PV
      QV(:)=MATMUL(ALP(:,:),PV(:))
      !
      CALPHA=SUM(PV(:)*QV(:))
      CALPHA = CRHO / CALPHA
      !
      XV = XV + CALPHA * PV
      RV = RV - CALPHA * QV
      !
      CHECK(:)=MATMUL(ALP(:,:),XV(:))
      CONDIT=SUM(((CHECK-BET)/(DABS(BET)+1.0D-5))**2)
      CONDIT=CONDIT/DBLE(NPAR)
      IF ((CONDIT.LT.1.0D-3).OR.(DCONDIT.GT.1.0D30)) THEN
!        PRINT*, ' UBICGSTAB method:'
!        PRINT*, 'NITER= ', ITER, ' ; CONDIT= ', CONDIT
!        PRINT*, '--------------------------------------'
        DELT(:)=XV(:)
        FAILED=.FALSE.
        EXIT
      ENDIF
      !
      DCONDIT=CONDIT-PCONDIT
      PCONDIT=CONDIT
      !
    ENDDO
!    PRINT*, 'NITER= ', ITER, ' ; CONDIT= ', CONDIT
!    PRINT*, '**************************************'
    IF (FAILED.EQV..TRUE.) DELT(:)=XV(:)/10.0D0
    !
  END SUBROUTINE NEW_CALCULATE_DELTA_UBICGSTAB
  !
  !
END MODULE LM
!

