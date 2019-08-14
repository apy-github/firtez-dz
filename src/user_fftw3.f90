!
MODULE USER_FFTW3
  !
  !================================================
  !
  USE, INTRINSIC :: ISO_C_BINDING
  USE CONS_PARAM, ONLY: DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INCLUDE 'fftw3.f03'
  !
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:) :: FFTK1D, FFTY1D
  TYPE(C_PTR), ALLOCATABLE, DIMENSION(:,:) :: GPLAN1D
  INTEGER, ALLOCATABLE, DIMENSION(:) :: FFTINI, FFTEND
  !
  ! 2D convolution (spatial)
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: FFTK2D, FFTY2D
  TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: GPLAN2D
  !
  ! Smoothing 2D
  COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: FFTS2D, FFTSY2D
  TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: GPLANS2D
  ! Accesibility:
  !   Variables:
  PUBLIC   :: FFTK1D, FFTY1D
  PUBLIC   :: GPLAN1D
  PUBLIC   :: FFTINI, FFTEND
  !
  PUBLIC   :: FFTK2D, FFTY2D
  PUBLIC   :: GPLAN2D
  !
  PUBLIC   :: FFTS2D, FFTSY2D
  PUBLIC   :: GPLANS2D
  !   Functions:
  PUBLIC :: FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
  PUBLIC   :: MAKE_1D_PLAN, MAKE_1D_IPLAN
  !
  PUBLIC   :: MAKE_2D_PLAN
  PUBLIC   :: MAKE_2D_IPLAN
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! make_1d_plan
  ! make_1d_iplan
  ! perform_fftf_conv1d
  ! make_2d_plan
  ! make_2d_iplan
  ! perform_fftf_conv2d
  ! remove_plan
  !
  !------------------------------------------------
  !
  SUBROUTINE MAKE_1D_PLAN(IPLAN, SY, Y, SFY, FFTY)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: IPLAN
    INTEGER, INTENT(IN) :: SY, SFY
    REAL(DP), INTENT(INOUT), DIMENSION(SY) :: Y
    COMPLEX(DP), INTENT(INOUT), DIMENSION(SFY) :: FFTY
    !
    iplan = fftw_plan_dft_r2c_1d(SY,Y,FFTY,FFTW_ESTIMATE)
    !
  END SUBROUTINE MAKE_1D_PLAN
  !
  !------------------------------------------------
  !
  SUBROUTINE MAKE_1D_IPLAN(IPLAN,SFY,FFTY,SY,Y)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: IPLAN
    INTEGER, INTENT(IN) :: SY, SFY
    REAL(DP), INTENT(INOUT), DIMENSION(SY) :: Y
    COMPLEX(DP), INTENT(INOUT), DIMENSION(SFY) :: FFTY
    !
    iplan = fftw_plan_dft_c2r_1d(SY,FFTY,Y,FFTW_ESTIMATE)
    !
  END SUBROUTINE MAKE_1D_IPLAN
  !
  !------------------------------------------------
  !
  SUBROUTINE PERFORM_FFTF_CONV1D(FPLAN, IPLAN, SY, Y, K, SFY, FFTY, FFTK, RES)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: FPLAN, IPLAN
    INTEGER, INTENT(IN) :: SY, SFY
    REAL(DP), INTENT(INOUT), DIMENSION(SY) :: Y, K, RES
    COMPLEX(DP), INTENT(INOUT), DIMENSION(SFY) :: FFTY, FFTK
    COMPLEX(DP), DIMENSION(SFY)                :: FFTP
    !
    ! Forward discrete transformation for both, function and kernel:
    CALL FFTW_EXECUTE_DFT_R2C(FPLAN,Y,FFTY) 
    CALL FFTW_EXECUTE_DFT_R2C(FPLAN,K,FFTK) 
    !
    ! Multiply the transformed products
    FFTP=FFTK*FFTY
    ! Back transform the product
    CALL FFTW_EXECUTE_DFT_C2R(IPLAN,FFTP,RES) 
    !
    ! Returned tranformation is multiplied by the size of the array
    RES=RES/DBLE(SY)
    !
  END SUBROUTINE PERFORM_FFTF_CONV1D
  !
  !------------------------------------------------
  !
  SUBROUTINE MAKE_2D_PLAN(IPLAN, SHY, Y, SHFY, FFTY)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: IPLAN
    INTEGER, INTENT(IN), DIMENSION(2) :: SHY, SHFY
    REAL(DP), INTENT(INOUT), DIMENSION(SHY(1), SHY(2)) :: Y
    COMPLEX(DP), INTENT(INOUT), DIMENSION(SHFY(1),SHFY(2)) :: FFTY
    !
    iplan = fftw_plan_dft_r2c_2d(SHY(2),SHY(1),Y,FFTY,FFTW_ESTIMATE)
    !
  END SUBROUTINE MAKE_2D_PLAN
  !
  !------------------------------------------------
  !
  SUBROUTINE MAKE_2D_IPLAN(IPLAN,SHFY,FFTY,SHY,Y)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: IPLAN
    INTEGER, INTENT(IN), DIMENSION(2) :: SHY, SHFY
    REAL(DP), INTENT(INOUT), DIMENSION(SHY(1),SHY(2)) :: Y
    COMPLEX(DP), INTENT(INOUT), DIMENSION(SHFY(1),SHFY(2)) :: FFTY
    !
    iplan = fftw_plan_dft_c2r_2d(SHY(2),SHY(1),FFTY,Y,FFTW_ESTIMATE)
    !
  END SUBROUTINE MAKE_2D_IPLAN
  !
  !------------------------------------------------
  !
  SUBROUTINE PERFORM_FFTF_CONV2D(FPLAN, IPLAN, SHY, Y, K, SHFY, FFTY, FFTK, RES)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: FPLAN, IPLAN
    INTEGER, INTENT(IN), DIMENSION(2) :: SHY, SHFY
    REAL(DP), INTENT(INOUT), DIMENSION(SHY(1), SHY(2)) :: Y, K, RES
    COMPLEX(DP), INTENT(INOUT), DIMENSION(SHFY(1), SHFY(2)) :: FFTY, FFTK
    COMPLEX(DP), DIMENSION(SHFY(1),SHFY(2))                :: FFTP
    !
    ! Forward discrete transformation for both, function and kernel:
    CALL FFTW_EXECUTE_DFT_R2C(FPLAN,Y,FFTY) 
    CALL FFTW_EXECUTE_DFT_R2C(FPLAN,K,FFTK) 
    !
    ! Multiply the transformed products
    FFTP=FFTK*FFTY
    ! Back transform the product
    CALL FFTW_EXECUTE_DFT_C2R(IPLAN,FFTP,RES) 
    !
    ! Returned tranformation is multiplied by the size of the array
    RES=RES/PRODUCT(DBLE(SHY))
    !
  END SUBROUTINE PERFORM_FFTF_CONV2D
  !
  !------------------------------------------------
  !
  SUBROUTINE REMOVE_PLAN(IPLAN)
    IMPLICIT NONE
    TYPE(C_PTR), INTENT(INOUT) :: IPLAN
    CALL fftw_destroy_plan(iplan)
  END SUBROUTINE REMOVE_PLAN
  !
  !================================================
  !
END MODULE USER_FFTW3
!
