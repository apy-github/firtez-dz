!
MODULE ALLOCATE_UTILS
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: ALLOCATE_1D_IP
  PUBLIC :: ALLOCATE_1D_DP
  PUBLIC :: ALLOCATE_2D_IP
  PUBLIC :: ALLOCATE_2D_DP
  PUBLIC :: ALLOCATE_3D_DP
  PUBLIC :: ALLOCATE_4D_DP
  PUBLIC :: ALLOCATE_1D_SP
  PUBLIC :: ALLOCATE_L1D_SP
  PUBLIC :: ALLOCATE_2D_SP
  PUBLIC :: ALLOCATE_3D_SP
  PUBLIC :: ALLOCATE_4D_SP
  PUBLIC :: ALLOCATE_5D_SP
  PUBLIC :: ALLOCATE_6D_SP
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! allocate_1d_ip
  ! allocate_1d_dp
  ! allocate_2d_ip
  ! allocate_2d_dp
  ! allocate_3d_dp
  ! allocate_4d_dp
  ! allocate_1d_sp
  ! allocate_l1d_sp
  ! allocate_2d_sp
  ! allocate_3d_sp
  ! allocate_4d_sp
  ! allocate_5d_sp
  ! allocate_6d_sp
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_1D_IP(ARRAY1D,N1,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1
    INTEGER, ALLOCATABLE, INTENT(INOUT)     :: ARRAY1D(:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY1D)) DEALLOCATE(ARRAY1D)
    !
    ALLOCATE(ARRAY1D(N1), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY1D(:)=0
    !
  END SUBROUTINE ALLOCATE_1D_IP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_1D_DP(ARRAY1D,N1,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1
    REAL(DP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY1D(:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY1D)) DEALLOCATE(ARRAY1D)
    !
    ALLOCATE(ARRAY1D(N1), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY1D(:)=0.D0
    !
  END SUBROUTINE ALLOCATE_1D_DP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_2D_IP(ARRAY2D,N1,N2,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1, N2
    INTEGER, ALLOCATABLE, INTENT(INOUT)     :: ARRAY2D(:,:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY2D)) DEALLOCATE(ARRAY2D)
    !
    ALLOCATE(ARRAY2D(N1,N2), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY2D(:,:)=0
    !
  END SUBROUTINE ALLOCATE_2D_IP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_2D_DP(ARRAY2D,N1,N2,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1, N2
    REAL(DP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY2D(:,:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY2D)) DEALLOCATE(ARRAY2D)
    !
    ALLOCATE(ARRAY2D(N1, N2), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY2D(:,:)=0.D0
    !
  END SUBROUTINE ALLOCATE_2D_DP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_3D_DP(ARRAY3D,N1,N2,N3,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1, N2, N3
    REAL(DP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY3D(:,:,:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY3D)) DEALLOCATE(ARRAY3D)
    !
    ALLOCATE(ARRAY3D(N1, N2, N3), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY3D(:,:,:)=0.D0
    !
  END SUBROUTINE ALLOCATE_3D_DP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_4D_DP(ARRAY4D,N1,N2, N3,N4,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1, N2, N3, N4
    REAL(DP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY4D(:,:,:,:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY4D)) DEALLOCATE(ARRAY4D)
    !
    ALLOCATE(ARRAY4D(N1, N2, N3, N4), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY4D(:,:,:,:)=0.D0
    !
  END SUBROUTINE ALLOCATE_4D_DP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_1D_SP(ARRAY1D,N1,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1
    REAL(SP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY1D(:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY1D)) DEALLOCATE(ARRAY1D)
    !
    ALLOCATE(ARRAY1D(N1), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY1D(:)=0.E0
    !
  END SUBROUTINE ALLOCATE_1D_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_L1D_SP(ARRAY1D,N1,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER(KIND=8), INTENT(IN)             :: N1
    REAL(SP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY1D(:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY1D)) DEALLOCATE(ARRAY1D)
    !
    ALLOCATE(ARRAY1D(N1), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY1D(:)=0.E0
    !
  END SUBROUTINE ALLOCATE_L1D_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_2D_SP(ARRAY2D,N1,N2,TEXT)
    !
    CHARACTER(*)                            :: TEXT
    INTEGER, INTENT(IN)                     :: N1, N2
    REAL(SP), ALLOCATABLE, INTENT(INOUT)    :: ARRAY2D(:,:)
    INTEGER                                 :: IERR
    !
    IF (ALLOCATED(ARRAY2D)) DEALLOCATE(ARRAY2D)
    !
    ALLOCATE(ARRAY2D(N1, N2), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY2D(:,:)=0.E0
    !
  END SUBROUTINE ALLOCATE_2D_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_3D_SP(ARRAY3D,N1,N2,N3,TEXT)
    !
    CHARACTER(*)                         :: TEXT
    INTEGER, INTENT(IN)                  :: N1, N2, N3
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: ARRAY3D(:,:,:)
    INTEGER                              :: IERR
    !
    IF (ALLOCATED(ARRAY3D)) DEALLOCATE(ARRAY3D)
    !
    ALLOCATE(ARRAY3D(N1, N2, N3), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY3D(:,:,:)=0.E0
    !
  END SUBROUTINE ALLOCATE_3D_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_4D_SP(ARRAY4D,N1,N2,N3,N4,TEXT)
    !
    CHARACTER(*)                         :: TEXT
    INTEGER, INTENT(IN)                  :: N1, N2, N3, N4
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: ARRAY4D(:,:,:,:)
    INTEGER                              :: IERR
    !
    IF (ALLOCATED(ARRAY4D)) DEALLOCATE(ARRAY4D)
    !
    ALLOCATE(ARRAY4D(N1,N2,N3,N4), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.',IERR,N1,N2,N3,N4
       STOP
    ENDIF
    !
    ARRAY4D(:,:,:,:)=0.E0
    !
  END SUBROUTINE ALLOCATE_4D_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_5D_SP(ARRAY5D,N1,N2,N3,N4,N5,TEXT)
    !
    CHARACTER(*)                         :: TEXT
    INTEGER, INTENT(IN)                  :: N1, N2, N3, N4, N5
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: ARRAY5D(:,:,:,:,:)
    INTEGER                              :: IERR
    !
    IF (ALLOCATED(ARRAY5D)) DEALLOCATE(ARRAY5D)
    !
    ALLOCATE(ARRAY5D(N1, N2, N3, N4, N5), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.'
       STOP
    ENDIF
    !
    ARRAY5D(:,:,:,:,:)=0.E0
    !
  END SUBROUTINE ALLOCATE_5D_SP
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_6D_SP(ARRAY6D,N1,N2,N3,N4,N5,N6,TEXT)
    !
    CHARACTER(*)                         :: TEXT
    INTEGER, INTENT(IN)                  :: N1, N2, N3, N4, N5, N6
    REAL(SP), ALLOCATABLE, INTENT(INOUT) :: ARRAY6D(:,:,:,:,:,:)
    INTEGER                              :: IERR
    !
    IF (ALLOCATED(ARRAY6D)) DEALLOCATE(ARRAY6D)
    !
    ALLOCATE(ARRAY6D(N1, N2, N3, N4, N5, N6), STAT=IERR)
    IF (IERR.NE.0) THEN
       PRINT*,'Not enought memory to allocate ', TEXT, ' array. STOP.', IERR
       STOP
    ENDIF
    !
    ARRAY6D(:,:,:,:,:,:)=0.E0
    !
  END SUBROUTINE ALLOCATE_6D_SP
  !
  !================================================
  !
END MODULE ALLOCATE_UTILS
!
