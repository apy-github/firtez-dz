!
MODULE COUPLED_INVERSION
  !
  USE COUPLED_PARAM
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    SUBROUTINE INIT_COUPLED_INVERSION_STATIC_VARS()
      !
      USE INVERT_PARAM, ONLY: NSTKINV, NFREEV
      USE GRID_PARAM, ONLY: NX,NY
      USE FORWARD_PARAM, ONLY: NUMW
      USE ALLOCATE_UTILS, ONLY: ALLOCATE_2D_DP,ALLOCATE_3D_DP &
          , ALLOCATE_2D_IP, ALLOCATE_1D_IP, ALLOCATE_1D_SP, ALLOCATE_1D_DP
      USE user_mpi, ONLY: mpi__myrank
      !
      INTEGER :: NDSVD
      INTEGER :: NTWL
      INTEGER :: NWL
      !
      NDSVD=NSTKINV*NUMW
      NTWL=NUMW
      NWL=NUMW
      !
      IF (mpi__myrank.EQ.0) THEN
        !
        ! SON INDEPENDIENTES DE ITERACION
        CALL ALLOCATE_3D_DP(COU_BETA,NDSVD,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_BETA')
        CALL ALLOCATE_3D_DP(COU_OBS,NWL,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_OBS')
        CALL ALLOCATE_3D_DP(COU_YMOD,NWL,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_YMOD')
        CALL ALLOCATE_3D_DP(COU_SIG,NWL,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_SIG')
        CALL ALLOCATE_3D_DP(COU_NFACT,NFREEV,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_NFACT')
        !
        ! SON INDEPENDIENTES DE ITERACION
        CALL ALLOCATE_2D_DP(COU_CH,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_CH')
        CALL ALLOCATE_2D_DP(COU_CHISQ,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_CHISQ')
        CALL ALLOCATE_2D_DP(COU_OCHISQ,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_OCHISQ')
        !CALL ALLOCATE_2D_IP(COU_ID,NY,NX &
        !    ,'INIT_COUPLED_INVERSION_VARS, COU_ID')
        !CALL ALLOCATE_2D_IP(COU_IWL1,NY,NX &
        !    ,'INIT_COUPLED_INVERSION_VARS, COU_IWL1')
        !CALL ALLOCATE_2D_IP(COU_IWLN,NY,NX &
        !    ,'INIT_COUPLED_INVERSION_VARS, COU_IWLN')
        !
        ! SON INDEPENDIENTES DE ITERACION
        CALL ALLOCATE_2D_DP(COU_SVDLAMBDA,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_SVDLAMBDA')
        CALL ALLOCATE_2D_IP(COU_STEPS,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_STEPS')
        CALL ALLOCATE_1D_IP(COU_LISTA,1 &
            ,'INIT_COUPLED_INVERSION_VARS, COU_LISTA')
        CALL ALLOCATE_1D_IP(COU_IPGRP,1 &
            ,'INIT_COUPLED_INVERSION_VARS, COU_IPGRP')
        !
        ! COU_ID STORES THE NUMBER OF ELEMENTS IN NDSVD ALREADY FILLED...
        ! ... IT MUST BE 0 AT THE BEGINING AND MUST BE CHANGED BY beta2dc
        !
        !COU_IWL1(:,:)=1
        !COU_IWLN(:,:)=NWL
        !
        COU_NPGRP=1
        !
        ! COU_SVDLAMBDA stores the lambda parameter of the lm for each x,y
        !COU_SVDLAMBDA(:,:)=1.D0
        !COU_SVDLAMBDA=1.0D0
        !COU_STEPS=0.0D0
        !
        ! COU_OCHISQ stores the lambda parameter of the lm for each x,y
        !COU_OCHISQ(:,:)=1.D30
      ELSE
        CALL ALLOCATE_1D_DP(RCV_NFACT1D,NFREEV,'COU: RCV_DA1D')
      ENDIF
      !
    END SUBROUTINE INIT_COUPLED_INVERSION_STATIC_VARS
    !
    !
    SUBROUTINE INIT_COUPLED_INVERSION_DYNAMIC_VARS()
      !
      USE INVERT_PARAM, ONLY: NSTKINV, IFREEP
      USE GRID_PARAM, ONLY: NX,NY
      USE FORWARD_PARAM, ONLY: NUMW
      USE ALLOCATE_UTILS, ONLY: ALLOCATE_3D_DP,ALLOCATE_4D_DP,ALLOCATE_1D_DP
      USE user_mpi, ONLY: mpi__myrank
      !
      INTEGER :: NDSVD
      INTEGER :: MFIT
      INTEGER :: NTWL
      INTEGER :: NWL
      !
!PRINT*, 'I am ', mpi__myrank &
!    , ' and I am starting INIT_COUPLED_INVERSION_DYNAMIC_VARS'
      NDSVD=NSTKINV*NUMW
      MFIT=IFREEP
      NTWL=NUMW
      NWL=NUMW
      !
      IF (mpi__myrank.EQ.0) THEN
        !
        ! SON DEPENDIENTES DE ITERACION
        CALL ALLOCATE_4D_DP(COU_ALPHA,NDSVD,MFIT,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_ALPHA')
        CALL ALLOCATE_4D_DP(COU_DYDA,NTWL,MFIT,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_DYDA')
        CALL ALLOCATE_3D_DP(COU_DA,MFIT,NY,NX &
            ,'INIT_COUPLED_INVERSION_VARS, COU_DA')
        !
        ! As we are starting a new cycle for a larger number of slabs...
        ! ... set these two to starting point:
        COU_OCHISQ(:,:)=1.D30
        COU_CHISQ(:,:)=1.D30   ! Esto... no deberia ser 0?
        COU_CH(:,:)=1.D30      ! Esto... no deberia ser 0?
        COU_SVDLAMBDA(:,:)=1.D3
        COU_STEPS(:,:)=0
        !COU_SVDLAMBDA=1.0D0
        !COU_STEPS=0.0D0
        COU_PRECYC=0.0D0
        !
      ELSE
        CALL ALLOCATE_1D_DP(RCV_DA1D,IFREEP,'COU: RCV_DA1D')
      ENDIF
!PRINT*, 'I am ', mpi__myrank &
!    , ' and I am ending INIT_COUPLED_INVERSION_DYNAMIC_VARS'
      !
    END SUBROUTINE INIT_COUPLED_INVERSION_DYNAMIC_VARS
    !
    SUBROUTINE END_COUPLED_INVERSION_STATIC_VARS()
      !
      USE user_mpi, ONLY: mpi__myrank
      !
      IF (mpi__myrank.EQ.0) THEN
        !
        DEALLOCATE(COU_BETA)
        DEALLOCATE(COU_OBS)
        DEALLOCATE(COU_YMOD)
        DEALLOCATE(COU_SIG)
        DEALLOCATE(COU_NFACT)
        !
        DEALLOCATE(COU_CH)
        DEALLOCATE(COU_CHISQ)
        DEALLOCATE(COU_OCHISQ)
        !DEALLOCATE(COU_ID)
        !DEALLOCATE(COU_IWL1)
        !DEALLOCATE(COU_IWLN)
        !
        DEALLOCATE(COU_SVDLAMBDA)
        DEALLOCATE(COU_LISTA)
        DEALLOCATE(COU_IPGRP)
        DEALLOCATE(COU_STEPS)
        !
      ELSE
        DEALLOCATE(RCV_NFACT1D)
      ENDIF
       !
    END SUBROUTINE END_COUPLED_INVERSION_STATIC_VARS
    !
    SUBROUTINE INIT_EVERY_CYC_EVERY_IT()
      !
      USE user_mpi, ONLY: mpi__myrank
      USE INVERT_PARAM, ONLY: AM_I_DONE
      !
      IF (mpi__myrank.EQ.0) THEN
        !
        ! 
        !COU_ID(:,:)=1
        COU_BETA(:,:,:)=0.0D0
        COU_ALPHA(:,:,:,:)=0.0D0
        !COU_DYDA(:,:,:,:)=0.0D0
        COU_DA(:,:,:)=0.0D0
        COU_CHISQ(:,:)=0.0D0
        COU_CH(:,:)=0.0D0
        AM_I_DONE(:,:)=0
        !
      ENDIF
      !
    END SUBROUTINE INIT_EVERY_CYC_EVERY_IT
    !
    !
    SUBROUTINE END_COUPLED_INVERSION_DYNAMIC_VARS()
      !
      USE user_mpi, ONLY: mpi__myrank
      !
!PRINT*, 'I am ', mpi__myrank &
!    , ' and I am starting END_COUPLED_INVERSION_DYNAMIC_VARS'
      IF (mpi__myrank.EQ.0) THEN
        !
        ! Depend on the cycle number of slabs considered:
        DEALLOCATE(COU_ALPHA)
        !DEALLOCATE(COU_DYDA)
        DEALLOCATE(COU_DA)
        !
      ELSE
        DEALLOCATE(RCV_DA1D)
      ENDIF
      !
!PRINT*, 'I am ', mpi__myrank &
!    , ' and I am ending END_COUPLED_INVERSION_DYNAMIC_VARS'
    END SUBROUTINE END_COUPLED_INVERSION_DYNAMIC_VARS
    !
END MODULE COUPLED_INVERSION
!
