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
  PUBLIC GET_PERTURBATION
  PUBLIC HESSIAN
  PUBLIC GET_INVERSE_COV
  !
  PRIVATE GET_TERM1
  PRIVATE GET_INVERSE
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
  SUBROUTINE GET_PERTURBATION(NP, NF, IJACOBIAN, DIFF, LAMBDAP, DELTA, SVDTOL)
    !
    INTEGER                                           :: NP, NF
    !
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:)      :: IJACOBIAN(NF,NP)
    DOUBLE PRECISION, INTENT(IN)                      :: LAMBDAP
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:)        :: DIFF(NF)
    !
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:)     :: DELTA(NP)
    INTEGER, INTENT(IN), OPTIONAL                     :: SVDTOL
    !
    DOUBLE PRECISION, DIMENSION(:,:)                  :: JACOBIAN(NP,NF)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: HSS(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM1(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: ITERM1(NP,NP)
    DOUBLE PRECISION, DIMENSION(:,:)                  :: TERM2B(NP,NF)
    !
    JACOBIAN(:,:)=TRANSPOSE(IJACOBIAN(:,:))
    ! Get hessian:
    CALL HESSIAN(NP, NF, JACOBIAN, HSS)
    !
    ! Get Term 1:
    CALL GET_TERM1(NP, HSS, LAMBDAP, TERM1)
    !
    ! Get Term 2:
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
    TERM2B(:,:)=MATMUL(ITERM1, JACOBIAN)
    DELTA(:)=MATMUL(TERM2B, DIFF)
    !
  END SUBROUTINE GET_PERTURBATION
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
  !================================================
  !
END MODULE LM
!
