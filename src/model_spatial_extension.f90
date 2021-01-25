!
MODULE EXTEND_2D
  !
  !================================================
  !
  USE user_mpi, ONLY: mpi__myrank
  USE CONS_PARAM, ONLY: SP
  USE GRID_PARAM, ONLY: NX, NY
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_2D_SP
  USE FORWARD_PARAM, ONLY: ATM_ARGS
  USE INVERT_PARAM, ONLY: AM_I_DONE, IMASK, CYCPOW
  USE PHYS_PARAM
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: EXTEND_MODEL
  PUBLIC :: SETUP_MASK
  PUBLIC :: SMOOTHING
  PUBLIC :: SMOOTHING_B
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! setup_mask
  ! extend_model
  !
  !------------------------------------------------
  !
  SUBROUTINE SETUP_MASK(IT, NMAX)
    !
    INTEGER, INTENT(IN)         :: IT, NMAX
    !
    INTEGER            :: XSTEP, YSTEP, XINIT, YINIT
    !
IF (IT.NE.NMAX) THEN
    !XSTEP=MINVAL((/NMAX+1-IT, NX/))
    !YSTEP=MINVAL((/NMAX+1-IT, NY/))
    XSTEP=MINVAL((/CYCPOW**(NMAX-IT), NX/))
    YSTEP=MINVAL((/CYCPOW**(NMAX-IT), NY/))
    !
    XINIT=MAXVAL((/XSTEP/2,1/))
    YINIT=MAXVAL((/YSTEP/2,1/))
    !
    AM_I_DONE(:,:)=1
    AM_I_DONE(YINIT:NY:YSTEP,XINIT:NX:XSTEP)=0
    IMASK=REAL(1-AM_I_DONE)
    !
IF (mpi__myrank.EQ.0) THEN
PRINT*, NMAX
PRINT*, IT
PRINT*, XSTEP
PRINT*, YSTEP
PRINT*, XINIT
PRINT*, YINIT
PRINT*, SIZE(AM_I_DONE), SUM(AM_I_DONE)
ENDIF
ELSE
  AM_I_DONE(:,:)=0
  IMASK(:,:)=1.E0
ENDIF
    !
  END SUBROUTINE SETUP_MASK
  !
  !------------------------------------------------
  !
  SUBROUTINE EXTEND_MODEL(IT, NMAX)
    !
    INTEGER, INTENT(IN)         :: IT, NMAX
    !
    INTEGER            :: J, K, M, L, I
    INTEGER            :: XSTEP, YSTEP, XINIT, YINIT, XNUM, YNUM
    !
    REAL(SP), DIMENSION(:,:), ALLOCATABLE :: XWEI, YWEI
    REAL(SP), DIMENSION(:,:), ALLOCATABLE :: WEI1, WEI2, WEI3, WEI4
    !
    INTEGER                               :: XINI_SQR, XFIN_SQR, YINI_SQR, YFIN_SQR
    REAL(SP), POINTER                     :: ARR(:,:,:)
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      IF (IT.NE.NMAX) THEN
        !
        XSTEP=MINVAL((/CYCPOW**(NMAX-IT), NX/))
        YSTEP=MINVAL((/CYCPOW**(NMAX-IT), NY/))
        !
        XINIT=MAXVAL((/XSTEP/2,1/))
        YINIT=MAXVAL((/YSTEP/2,1/))
        !
        !DO J=XINIT,NX,XSTEP
        !  DO K=YINIT,NY,YSTEP
        !    TEM3D(:,K,J)=1000.D0*(J/REAL(NX)*REAL(NY)+K)
        !  ENDDO
        !ENDDO
        !
        CALL ALLOCATE_2D_SP(XWEI, XSTEP+1, YSTEP+1 &
            , 'XWEI, model_spatial_extension.f90')
        CALL ALLOCATE_2D_SP(YWEI, XSTEP+1, YSTEP+1 &
            , 'YWEI, model_spatial_extension.f90')
        !
        CALL ALLOCATE_2D_SP(WEI1, XSTEP+1, YSTEP+1 &
            , 'WEI1, model_spatial_extension.f90')
        CALL ALLOCATE_2D_SP(WEI2, XSTEP+1, YSTEP+1 &
            , 'WEI2, model_spatial_extension.f90')
        CALL ALLOCATE_2D_SP(WEI3, XSTEP+1, YSTEP+1 &
            , 'WEI3, model_spatial_extension.f90')
        CALL ALLOCATE_2D_SP(WEI4, XSTEP+1, YSTEP+1 &
            , 'WEI4, model_spatial_extension.f90')
        !
        DO J=1,YSTEP+1
          YWEI(J,:)=ABS((REAL(YSTEP-(J-1)))/REAL(YSTEP))
        ENDDO
        DO J=1,XSTEP+1
          XWEI(:,J)=ABS(REAL(XSTEP-(J-1))/REAL(XSTEP))
        ENDDO
        !
        WEI1=XWEI*YWEI
        WEI2=XWEI*(1.-YWEI)
        WEI3=(1.-XWEI)*YWEI
        WEI4=(1.-XWEI)*(1.-YWEI)
        !
        XNUM=CEILING(REAL(NX-XINIT)/REAL(XSTEP))
        YNUM=CEILING(REAL(NY-YINIT)/REAL(YSTEP))
        !
        !
        DO I=1,ATM_ARGS-1
          !
          SELECT CASE (I)
          CASE(1)
             ARR => TEM3D
          CASE(2)
             ARR => PG3D
          CASE(3)
             ARR => RHO3D
          CASE(6)
            ARR => BX3D
          CASE(7)
            ARR => BY3D
          CASE(8)
            ARR => BZ3D
          CASE(11)
            ARR => VZ3D
          CASE(12)
            ARR => TAU3D5
          CASE DEFAULT
            CYCLE
          ENDSELECT
          !
          DO J=1,XNUM-1
            DO K=1,YNUM-1
              ! DEFINE SQUARE:
              XINI_SQR=(J-1)*XSTEP+XINIT
              XFIN_SQR=XINI_SQR+XSTEP
              YINI_SQR=(K-1)*YSTEP+YINIT
              YFIN_SQR=YINI_SQR+YSTEP
              ! TREAT EACH SQUARE:
              DO L=1,XSTEP+1
                DO M=1,YSTEP+1
                  !
                  ARR(:, YINI_SQR+(M-1), XINI_SQR+(L-1)) &
                      =WEI1(M,L)*ARR(:,YINI_SQR,XINI_SQR) &
                      +WEI2(M,L)*ARR(:,YFIN_SQR,XINI_SQR) &
                      +WEI3(M,L)*ARR(:,YINI_SQR,XFIN_SQR) &
                      +WEI4(M,L)*ARR(:,YFIN_SQR,XFIN_SQR)
                  !
                ENDDO
              ENDDO
              !
            ENDDO
          ENDDO
          !
          ! EDGES:
          ! X
          IF (XINIT.NE.1) THEN
            DO J=1,XINIT
              ARR(:,:,J)=ARR(:,:,XINIT)
            ENDDO
          ENDIF
          IF (XINIT+(XNUM-1)*XSTEP.NE.NX) THEN
            DO J=XINIT+(XNUM-1)*XSTEP,NX
              ARR(:,:,J)=ARR(:,:,XINIT+(XNUM-1)*XSTEP)
            ENDDO
          ENDIF
          ! Y
          IF (YINIT.NE.1) THEN
            DO J=1,YINIT
              ARR(:,J,:)=ARR(:,YINIT,:)
            ENDDO
          ENDIF
          IF (YINIT+(YNUM-1)*YSTEP.NE.NY) THEN
            DO J=YINIT+(YNUM-1)*YSTEP,NY
              ARR(:,J,:)=ARR(:,YINIT+(YNUM-1)*YSTEP,:)
            ENDDO
          ENDIF
          !
          NULLIFY(ARR)
          !
        ENDDO
        !
      ENDIF
      !
    ENDIF
    !
  END SUBROUTINE EXTEND_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE SMOOTHING()
    !
    USE user_mpi, ONLY: mpi__ierror, MPI_COMM_WORLD
    USE CONS_PARAM, ONLY: DP
    USE USER_FFTW3, ONLY: FFTS2D, GPLANS2D, FFTSY2D, MAKE_2D_PLAN &
        , MAKE_2D_IPLAN, FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
    USE GRID_PARAM, ONLY: NY, NX, NZ
    USE PHYS_PARAM
    USE INVERT_PARAM, ONLY: INV_ATMPAR, CURIC, MAXITER
    !
    INTEGER      :: J,K
    REAL(DP), DIMENSION(NZ,NY,NX) :: ARRAY
    !
    INTEGER                    :: I
    INTEGER                    :: NY2, NX2
    REAL(DP), DIMENSION(NY,NX) :: X2D, Y2D, SKERNEL
    REAL(DP)                   :: ESIG, ISIG
    !
    NY2=NY/2+1
    NX2=NX/2+1
    !
    X2D(:,:)=0.0D0
    Y2D(:,:)=0.0D0
    DO I=1,NY
      Y2D(I,:)=DBLE(I-NY2)
    ENDDO
    DO I=1,NX
      X2D(:,I)=DBLE(I-NX2)
    ENDDO
    !
    ESIG=3.0D0-2.0D0*DBLE(CURIC)/DBLE(MAXITER)!2.0D0
PRINT*, 'Esig: ', ESIG
    ISIG=0.9D0
    !
    SKERNEL(:,:)=0.0D0
    SKERNEL(:,:)=DEXP(-Y2D**2/(2.0D0*ESIG**2)) &
        *DEXP(-X2D**2/(2.0D0*ESIG**2))
!#    SKERNEL(:,:)=SKERNEL(:,:) &
!#        -(DEXP(-Y2D**2/(2.0D0*ISIG**2)) &
!#        *DEXP(-X2D**2/(2.0D0*ISIG**2)))
    !
    ! Shift:
    SKERNEL=CSHIFT(SKERNEL,SHIFT=NY2-1,DIM=1)
    SKERNEL=CSHIFT(SKERNEL,SHIFT=NX2-1,DIM=2)
    !
    ! Normalize:
    SKERNEL=SKERNEL/SUM(SKERNEL)
    ! Second normalization: (this is due to the factorization of fftw3)
    SKERNEL=SKERNEL/DBLE(NX)/DBLE(NY)
    !
    ALLOCATE(GPLANS2D(2))
    ALLOCATE(FFTS2D(NY2,NX))
    ALLOCATE(FFTSY2D(NY2,NX))
    FFTS2D(:,:)=0.0D0
    FFTSY2D(:,:)=0.0D0
    !
    ! Define a forward plan:
    CALL MAKE_2D_PLAN(GPLANS2D(1), SHAPE(SKERNEL), SKERNEL, SHAPE(FFTS2D), FFTS2D)
    ! Define an inverse plan:
    CALL MAKE_2D_IPLAN(GPLANS2D(2), SHAPE(FFTS2D), FFTS2D, SHAPE(SKERNEL), SKERNEL)
    !
    ! Store kernel transformation:
    CALL FFTW_EXECUTE_DFT_R2C(GPLANS2D(1),SKERNEL,FFTS2D)
    !
    !
    IF (mpi__myrank.EQ.0) THEN
PRINT*, '*** Smoothing ***'
      !
      DO J=1,8
        !
        IF (INV_ATMPAR(J).EQV..FALSE.) CYCLE
        !
        SELECT CASE (J)
          CASE(1)
            ARRAY=DBLE(BEST_TEM3D)
          CASE(2)
            ARRAY=DBLE(BEST_PG3D)
          CASE(3)
            ARRAY=DBLE(BEST_RHO3D)
          CASE(4)
            ARRAY=DBLE(BEST_BX3D)
          CASE(5)
            ARRAY=DBLE(BEST_BY3D)
          CASE(6)
            ARRAY=DBLE(BEST_BZ3D)
          CASE(7)
            ARRAY=DBLE(BEST_VZ3D)
          CASE(8)
            ARRAY=DBLE(BEST_PG3D)
        END SELECT
        !
        DO K=1,NZ
          ! DFT synthetic data:
          CALL FFTW_EXECUTE_DFT_R2C(GPLANS2D(1),ARRAY(K,:,:),FFTSY2D)
          !
          ! Multiply transformed arrays:
          FFTSY2D=FFTSY2D*FFTS2D
          ! DFT back the product:
          CALL FFTW_EXECUTE_DFT_C2R(GPLANS2D(2),FFTSY2D,ARRAY(K,:,:))
        ENDDO
        !
        SELECT CASE (J)
          CASE(1)
            TEM3D=REAL(ARRAY)
            BEST_TEM3D=REAL(ARRAY)
          CASE(2)
            PG3D=REAL(ARRAY)
            BEST_PG3D=REAL(ARRAY)
          CASE(3)
            RHO3D=REAL(ARRAY)
            BEST_RHO3D=REAL(ARRAY)
          CASE(4)
            BX3D=REAL(ARRAY)
            BEST_BX3D=REAL(ARRAY)
          CASE(5)
            BY3D=REAL(ARRAY)
            BEST_BY3D=REAL(ARRAY)
          CASE(6)
            BZ3D=REAL(ARRAY)
            BEST_BZ3D=REAL(ARRAY)
          CASE(7)
            VZ3D=REAL(ARRAY)
            BEST_VZ3D=REAL(ARRAY)
          CASE(8)
            PG3D=REAL(ARRAY)
            BEST_PG3D=REAL(ARRAY)
        END SELECT
        !
      ENDDO
      !
    ENDIF
    !
    DEALLOCATE(GPLANS2D)
    DEALLOCATE(FFTS2D)
    DEALLOCATE(FFTSY2D)
    !
  END SUBROUTINE SMOOTHING
  !
  !------------------------------------------------
  !
  SUBROUTINE SMOOTHING_B()
    !
    USE user_mpi, ONLY: mpi__ierror, MPI_COMM_WORLD
    USE CONS_PARAM, ONLY: DP
    USE USER_FFTW3, ONLY: FFTS2D, GPLANS2D, FFTSY2D, MAKE_2D_PLAN &
        , MAKE_2D_IPLAN, FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
    USE GRID_PARAM, ONLY: NY, NX, NZ
    USE PHYS_PARAM
    USE INVERT_PARAM, ONLY: INV_ATMPAR, CURIC, MAXITER
    !
    INTEGER      :: J,K
    REAL(DP), DIMENSION(NZ,NY,NX) :: ARRAY
    !
    INTEGER                    :: I
    INTEGER                    :: NY2, NX2
    REAL(DP), DIMENSION(NY,NX) :: X2D, Y2D, SKERNEL
    REAL(DP)                   :: ESIG, ISIG
    !
    NY2=NY/2+1
    NX2=NX/2+1
    !
    X2D(:,:)=0.0D0
    Y2D(:,:)=0.0D0
    DO I=1,NY
      Y2D(I,:)=DBLE(I-NY2)
    ENDDO
    DO I=1,NX
      X2D(:,I)=DBLE(I-NX2)
    ENDDO
    !
    ESIG=3.0D0-2.0D0*DBLE(CURIC)/DBLE(MAXITER)!2.0D0
!PRINT*, 'Esig: ', ESIG
    ISIG=0.9D0
    !
    SKERNEL(:,:)=0.0D0
    SKERNEL(:,:)=DEXP(-Y2D**2/(2.0D0*ESIG**2)) &
        *DEXP(-X2D**2/(2.0D0*ESIG**2))
!#    SKERNEL(:,:)=SKERNEL(:,:) &
!#        -(DEXP(-Y2D**2/(2.0D0*ISIG**2)) &
!#        *DEXP(-X2D**2/(2.0D0*ISIG**2)))
    !
    ! Shift:
    SKERNEL=CSHIFT(SKERNEL,SHIFT=NY2-1,DIM=1)
    SKERNEL=CSHIFT(SKERNEL,SHIFT=NX2-1,DIM=2)
    !
    ! Normalize:
    SKERNEL=SKERNEL/SUM(SKERNEL)
    ! Second normalization: (this is due to the factorization of fftw3)
    SKERNEL=SKERNEL/DBLE(NX)/DBLE(NY)
    !
    ALLOCATE(GPLANS2D(2))
    ALLOCATE(FFTS2D(NY2,NX))
    ALLOCATE(FFTSY2D(NY2,NX))
    FFTS2D(:,:)=0.0D0
    FFTSY2D(:,:)=0.0D0
    !
    ! Define a forward plan:
    CALL MAKE_2D_PLAN(GPLANS2D(1), SHAPE(SKERNEL), SKERNEL, SHAPE(FFTS2D), FFTS2D)
    ! Define an inverse plan:
    CALL MAKE_2D_IPLAN(GPLANS2D(2), SHAPE(FFTS2D), FFTS2D, SHAPE(SKERNEL), SKERNEL)
    !
    ! Store kernel transformation:
    CALL FFTW_EXECUTE_DFT_R2C(GPLANS2D(1),SKERNEL,FFTS2D)
    !
    !
    IF (mpi__myrank.EQ.0) THEN
PRINT*, '*** Smoothing ***'
      !
      DO J=1,8
        !
        IF (INV_ATMPAR(J).EQV..FALSE.) CYCLE
        !
        SELECT CASE (J)
          CASE(1)
            ARRAY=DBLE(BEST_TEM3D)
          CASE(2)
            ARRAY=DBLE(BEST_PG3D)
          CASE(3)
            ARRAY=DBLE(BEST_RHO3D)
          CASE(4)
            ARRAY=DBLE(BEST_BX3D)
          CASE(5)
            ARRAY=DBLE(BEST_BY3D)
          CASE(6)
            ARRAY=DBLE(BEST_BZ3D)
          CASE(7)
            ARRAY=DBLE(BEST_VZ3D)
          CASE(8)
            ARRAY=DBLE(BEST_PG3D)
        END SELECT
        !
        IF (NX.GT.4) THEN
          ARRAY(:,:,2:NX-1)=(ARRAY(:,:,3:NX)+ARRAY(:,:,1:NX-2)+ARRAY(:,:,2:NX-1)*2.0D0)/4.0E0
        ENDIF
        IF (NY.GT.4) THEN
          ARRAY(:,2:NY-1,:)=(ARRAY(:,3:NY,:)+ARRAY(:,1:NY-2,:)+ARRAY(:,2:NY-1,:)*2.0D0)/4.0E0
        ENDIF
        IF (NZ.GT.4) THEN
          ARRAY(2:NZ-1,:,:)=(ARRAY(3:NZ,:,:)+ARRAY(1:NZ-2,:,:)+ARRAY(2:NZ-1,:,:)*2.0D0)/4.0E0
        ENDIF

!!!        ARRAY(:,:,2:NX-1)=(ARRAY(:,:,3:NX)+ARRAY(:,:,1:NX-2))/2.0E0
!!!        ARRAY(:,2:NY-1,:)=(ARRAY(:,3:NY,:)+ARRAY(:,1:NY-2,:))/2.0E0
!!!        ARRAY(2:NZ-1,:,:)=(ARRAY(3:NZ,:,:)+ARRAY(1:NZ-2,:,:))/2.0E0
        !
        SELECT CASE (J)
          CASE(1)
            TEM3D=REAL(ARRAY)
            BEST_TEM3D=REAL(ARRAY)
          CASE(2)
            PG3D=REAL(ARRAY)
            BEST_PG3D=REAL(ARRAY)
          CASE(3)
            RHO3D=REAL(ARRAY)
            BEST_RHO3D=REAL(ARRAY)
          CASE(4)
            BX3D=REAL(ARRAY)
            BEST_BX3D=REAL(ARRAY)
          CASE(5)
            BY3D=REAL(ARRAY)
            BEST_BY3D=REAL(ARRAY)
          CASE(6)
            BZ3D=REAL(ARRAY)
            BEST_BZ3D=REAL(ARRAY)
          CASE(7)
            VZ3D=REAL(ARRAY)
            BEST_VZ3D=REAL(ARRAY)
          CASE(8)
            PG3D=REAL(ARRAY)
            BEST_PG3D=REAL(ARRAY)
        END SELECT
        !
      ENDDO
      !
    ENDIF
    !
    DEALLOCATE(GPLANS2D)
    DEALLOCATE(FFTS2D)
    DEALLOCATE(FFTSY2D)
    !
  END SUBROUTINE SMOOTHING_B
  !
  !================================================
  !
END MODULE EXTEND_2D
!
  SUBROUTINE SMOOTH_MODEL()
    !
  END SUBROUTINE SMOOTH_MODEL
