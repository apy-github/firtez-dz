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
!!!!!!!!!!!  PUBLIC PGET_PERTURBATION
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
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: JACOBIAN(NP,NF)
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:)   :: HSS(NP,NP)
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
    DOUBLE PRECISION                                :: LAMBDAP
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)    :: HSS(NP,NP)
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:) :: TERM1(NP,NP)
!    DOUBLE PRECISION, DIMENSION(:)    :: TP(NP)
    !
    TERM1(:,:)=HSS(:,:)
    DO I=1,NP
!      PRINT*, ' H(I)=', HSS(I,I)
!      TP(I)=TERM1(I,I)
      TERM1(I,I)=TERM1(I,I)+HSS(I,I)*lambdaP
    ENDDO
!PRINT*, TP
    !
  END SUBROUTINE GET_TERM1
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_INVERSE(NP, M, IM, SVDTOL)
    !
    INTEGER      :: NP, I
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)  :: M(NP,NP)
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:)  :: IM(NP,NP)
    INTEGER, INTENT(IN), OPTIONAL                  :: SVDTOL
    !
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: CM
    INTEGER                               :: INFO, LWORK
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SVDU, SVDVY, SVDIS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: SVDS, WORK
    INTEGER                               :: ISVDTOL
    INTEGER                                       :: CNT
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
    IF (PRESENT(SVDTOL)) THEN
      ISVDTOL=SVDTOL
    ELSE
      ISVDTOL=4
    END IF
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
    CNT=0
    DO I=1,NP
      IF (ABS(SVDS(I)) .LT. (10.D0**(-ISVDTOL) * MAXVAL(ABS(SVDS)))) THEN
        SVDIS(I,I)=0.D0
      ELSE
        CNT=CNT+1
        SVDIS(I,I)=1.D0/SVDS(I)
      ENDIF
    ENDDO
    !
    IM(:,:)=MATMUL(TRANSPOSE(SVDVY), MATMUL(SVDIS,TRANSPOSE(SVDU)))

!    DO I=1,NP
!      PRINT*, 'IHESS(I,I)=', IM(I,I)
!    ENDDO
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
  SUBROUTINE GET_PERTURBATION(NP, NF, IJACOBIAN, DIFF, LAMBDAP, DELTA, REG_PEN, REG_HSS, SVDTOL)
    !
    INTEGER                                           :: NP, NF
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: IJACOBIAN(NF,NP)
    DOUBLE PRECISION, INTENT(IN)                      :: LAMBDAP
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:)        :: DIFF(NF)
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:)        :: REG_PEN(NP)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: REG_HSS(NP,NP)
    !
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:)     :: DELTA(NP)
    INTEGER, INTENT(IN), OPTIONAL                     :: SVDTOL
    !
    DOUBLE PRECISION, DIMENSION(:,:)                  :: JACOBIAN(NP,NF)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: HSS(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM1(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: ITERM1(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM2B(NP,NF)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM2(NP)
    !
    JACOBIAN(:,:)=TRANSPOSE(IJACOBIAN(:,:))
    ! Get hessian:
    CALL HESSIAN(NP, NF, JACOBIAN, HSS)
    !
    ! Get Term 1:
    !
    ! If reg. add both hessians
    HSS(:,:)=HSS(:,:)+REG_HSS(:,:)
    !
    CALL GET_TERM1(NP, HSS, LAMBDAP, TERM1)
    !
    ! Get Term 2:
    TERM2(:)=MATMUL(JACOBIAN(:,:),DIFF(:))
    !
    ! If reg. subtract reg term:
    TERM2(:)=TERM2(:)-REG_PEN(:)
    !
    !
    ! Get Term 1 inverse:
!    CALL NEW_CALCULATE_DELTA_UBICGSTAB(NP,DELTA,TERM1,TERM2)
    !
    IF (PRESENT(SVDTOL)) THEN
      CALL GET_INVERSE(NP, TERM1, ITERM1, SVDTOL)
    ELSE
      CALL GET_INVERSE(NP, TERM1, ITERM1)
    ENDIF
    !
    ! Get Delta:
    !
    DELTA(:)=MATMUL(ITERM1(:,:), TERM2(:))
!PRINT*, DELTA
!PRINT*, TERM2
!PRINT*, MATMUL(TERM1(:,:), DELTA)
    !
  END SUBROUTINE GET_PERTURBATION
  !
  !------------------------------------------------
  !
  SUBROUTINE HESSIAN_REG(NP, NF, JACOBIAN, LACOBIAN, HSS)
    !
    INTEGER                                           :: NP, NF, I
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: JACOBIAN(NP,NF)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: LACOBIAN(NP,NP)
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:)   :: HSS(NP,NP)
    !
    HSS(:,:)=MATMUL(JACOBIAN, TRANSPOSE(JACOBIAN))
    HSS(:,:)=HSS(:,:)+MATMUL(LACOBIAN, TRANSPOSE(LACOBIAN))
    !
  END SUBROUTINE HESSIAN_REG
  !
  !------------------------------------------------
  !

  SUBROUTINE GET_PERTURBATION_REG(NP, NF, IJACOBIAN, LACOBIAN &
      , DIFF, RV,LAMBDAP, DELT, SVDTOL)
    !
    INTEGER                                           :: NP, NF
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: IJACOBIAN(NF,NP)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: LACOBIAN(NP,NP)
    DOUBLE PRECISION, INTENT(IN)                      :: LAMBDAP
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:)        :: DIFF(NF)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:)        :: RV(NP)
    !
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:)     :: DELT(NP)
    INTEGER, INTENT(IN), OPTIONAL                     :: SVDTOL
    !
    DOUBLE PRECISION, DIMENSION(:,:)                  :: JACOBIAN(NP,NF)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: HSS(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM1(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: ITERM1(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM2(NP)
    !
    JACOBIAN(:,:)=TRANSPOSE(IJACOBIAN(:,:))
    ! Get hessian:
    CALL HESSIAN_REG(NP, NF, JACOBIAN, LACOBIAN, HSS)
    !
    ! Get Term 1:
    CALL GET_TERM1(NP, HSS, LAMBDAP, TERM1)
    !
    ! Get Term 2:
    TERM2(:)=MATMUL(JACOBIAN,DIFF)
    TERM2(:)=TERM2(:)+MATMUL(TRANSPOSE(LACOBIAN),RV)
    !
    !
    ! Get Term 1 inverse:
    !
    IF (PRESENT(SVDTOL)) THEN
      CALL GET_INVERSE(NP, TERM1, ITERM1, SVDTOL)
    ELSE
      CALL GET_INVERSE(NP, TERM1, ITERM1)
    ENDIF
    !
    ! Get Delta:
    !
    DELT(:)=MATMUL(ITERM1, TERM2)
    !
  END SUBROUTINE GET_PERTURBATION_REG
  !
  !------------------------------------------------
  !
  SUBROUTINE GET_INVERSE_COV(NP, M, IM_COV)
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)    :: M(NP,NP)
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:,:) :: IM_COV(NP,NP)
    !
    INTEGER                                         :: NP, I
    !
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: CM
    INTEGER                                         :: INFO, LWORK
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: SVDU, SVDVY, SVDIS2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: SVDS, WORK
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
  !SUBROUTINE PGET_PERTURBATION_REG(NP, NF, IJACOBIAN, LACOBIAN &
  !    , DIFF, RV, LAMBDAP, DELT, SVDTOL)
!!!!!!!!!!!!!  SUBROUTINE PGET_PERTURBATION(NP, NF, IJACOBIAN &
!!!!!!!!!!!!!      , DIFF, LAMBDAP, DELT, SVDTOL)
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    INTEGER                                           :: NP, NF
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: IJACOBIAN(NF,NP)
!!!!!!!!!!!!!    DOUBLE PRECISION, INTENT(IN)                      :: LAMBDAP
!!!!!!!!!!!!!    DOUBLE PRECISION, INTENT(IN), DIMENSION(:)        :: DIFF(NF)
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:)     :: DELT(NP)
!!!!!!!!!!!!!    INTEGER, INTENT(IN), OPTIONAL                     :: SVDTOL
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE     :: JACOBIAN2
!!!!!!!!!!!!!    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE     :: LACOBIAN2
!!!!!!!!!!!!!    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE       :: RV2
!!!!!!!!!!!!!    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE       :: DELT2
!!!!!!!!!!!!!    INTEGER, DIMENSION(:)     :: NZE(NP)
!!!!!!!!!!!!!    INTEGER :: I, J, NPNZ
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    ! First, count number of non-zero elements:
!!!!!!!!!!!!!    NZE(:)=0
!!!!!!!!!!!!!    DO I=1,NP
!!!!!!!!!!!!!      IF (SUM(ABS(IJACOBIAN(:,I))).GT.1.0D-19) THEN
!!!!!!!!!!!!!        NZE(I)=1
!!!!!!!!!!!!!      ENDIF
!!!!!!!!!!!!!    ENDDO
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    NPNZ=SUM(NZE)
!!!!!!!!!!!!!    ALLOCATE(JACOBIAN2(NF,NPNZ))
!!!!!!!!!!!!!    ALLOCATE(DELT2(NPNZ))
!!!!!!!!!!!!!    JACOBIAN2(:,:)=0.0D0
!!!!!!!!!!!!!    DELT2(:)=0.0D0
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    NPNZ=0
!!!!!!!!!!!!!    DO I=1,NP
!!!!!!!!!!!!!      IF (NZE(I).GT.0) THEN
!!!!!!!!!!!!!        NPNZ=NPNZ+1
!!!!!!!!!!!!!        JACOBIAN2(:,NPNZ)=IJACOBIAN(:,I)
!!!!!!!!!!!!!      ENDIF
!!!!!!!!!!!!!    ENDDO
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!!PRINT*, ' [', 'Using ', NPNZ, ' out of ', NP,']'
!!!!!!!!!!!!!    CALL GET_PERTURBATION(NPNZ, NF, JACOBIAN2, DIFF &
!!!!!!!!!!!!!        , LAMBDAP, DELT2, SVDTOL)
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    NPNZ=0
!!!!!!!!!!!!!    DO I=1,NP
!!!!!!!!!!!!!      IF (NZE(I).GT.0) THEN
!!!!!!!!!!!!!        NPNZ=NPNZ+1
!!!!!!!!!!!!!        DELT(I)=DELT2(NPNZ)
!!!!!!!!!!!!!      ENDIF
!!!!!!!!!!!!!    ENDDO
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!    DEALLOCATE(JACOBIAN2)
!!!!!!!!!!!!!    DEALLOCATE(DELT2)
!!!!!!!!!!!!!    !
!!!!!!!!!!!!!  END SUBROUTINE PGET_PERTURBATION
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
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NPAR)  :: DELT
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NPAR,NPAR)  :: ALP
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NPAR)  :: BET
    !
    DOUBLE PRECISION, DIMENSION(NPAR) :: XV, PV, RV, QV, CHECK!, XV0
    DOUBLE PRECISION :: CRHO, CALPHA, CBETA, CRHO0
    LOGICAL :: FAILED
    !
    INTEGER :: ITER
    DOUBLE PRECISION :: CONDIT
    DOUBLE PRECISION :: DCONDIT
    DOUBLE PRECISION :: PCONDIT
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
!><  !================================================
!><  !
!><  SUBROUTINE NEW_CALCULATE_DELTA_UBICGSTAB(NPAR,DELT,ALP,BET,FAILED)
!><    !
!><    ! Unpreconditioned BiCGSTAB
!><    !   Implemented following NPAC Technical Report SCCS 691:
!><    !      Conjugate Gradient Algorithms in Fortran 90 and High...
!><    !      ...Performance Fortran
!><    !      K.A.Hawick, K.Dincer, G.Robinson, G.C.Fox
!><    !
!><    INTEGER, INTENT(IN)   :: NPAR
!><    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NPAR)  :: DELT
!><    LOGICAL, INTENT(INOUT)   :: FAILED
!><    DOUBLE PRECISION, INTENT(IN), DIMENSION(NPAR,NPAR)  :: ALP
!><    DOUBLE PRECISION, INTENT(IN), DIMENSION(NPAR)  :: BET
!><    !
!><    DOUBLE PRECISION, DIMENSION(NPAR) :: XV, PV, RV, QV, XV0, CHECK
!><    DOUBLE PRECISION :: CRHO, CALPHA, CBETA, CRHO0
!><    !
!><    INTEGER :: ITER
!><    DOUBLE PRECISION :: CONDIT
!><    DOUBLE PRECISION :: DCONDIT
!><    DOUBLE PRECISION :: PCONDIT
!><    !
!><    DELT(:)=0.0D0
!><    !
!><    XV(:) = 0.0D0
!><    !
!><    PV(:) = BET(:)
!><    RV(:) = BET(:)
!><    !
!><    CRHO=SUM(RV(:)*RV(:))
!><    !
!><    QV(:)=MATMUL(ALP(:,:),PV(:))
!><    !
!><    CALPHA=SUM(PV(:)*QV(:))
!><    CALPHA = CRHO / CALPHA
!><    XV = XV + CALPHA * PV
!><    RV = RV - CALPHA * QV
!><    !
!><    PCONDIT=1.0D29
!><    !
!><    DO ITER=1,100
!><      ! Update:
!><      XV0(:) = XV(:)
!><      CRHO0 = CRHO
!><      !
!><      CRHO=SUM(RV(:)*RV(:))
!><PRINT*, 'CRHO0: ', CRHO0
!><      CBETA = CRHO / CRHO0
!><      !
!><      PV = RV + CBETA * PV
!><      QV(:)=MATMUL(ALP(:,:),PV(:))
!><      !
!><      CALPHA=SUM(PV(:)*QV(:))
!><PRINT*, 'CALPHA: ', CALPHA
!><      CALPHA = CRHO / CALPHA
!><      !
!><      XV = XV + CALPHA * PV
!><      RV = RV - CALPHA * QV
!><      !
!><      CHECK(:)=MATMUL(ALP(:,:),XV(:))
!><      CONDIT=SUM(((CHECK-BET)/(DABS(BET)+1.0D-5))**2)
!><      CONDIT=CONDIT/DBLE(NPAR)
!><      !  PRINT*, ' UBICGSTAB method:'
!><      !  PRINT*, 'NITER= ', ITER, ' ; CONDIT= ', CONDIT &
!><      !      , ' ; DCONDIT= ', DCONDIT, ' ; PCONDIT= ' &
!><      !      , PCONDIT
!><      !  PRINT*, '--------------------------------------'
!><      IF ((CONDIT.LT.1.0D-3).OR.(DCONDIT.GT.0.0D0)) THEN
!><        PRINT*, ' UBICGSTAB method:'
!><        PRINT*, 'NITER= ', ITER, ' ; CONDIT= ', CONDIT
!><        PRINT*, '--------------------------------------'
!><        DELT(:)=XV(:)
!><        FAILED=.FALSE.
!><        EXIT
!><      ENDIF
!><      !
!><      DCONDIT=CONDIT-PCONDIT
!><      PCONDIT=CONDIT
!><      !
!><    ENDDO
!><    IF (FAILED.EQV..TRUE.) DELT(:)=XV(:)/10.0D0
!><    !
!><  END SUBROUTINE NEW_CALCULATE_DELTA_UBICGSTAB
!><  !  !
!><  !================================================
  !
END MODULE LM
!

