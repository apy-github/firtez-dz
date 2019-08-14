!
MODULE MAGSPLIT
  !
  !================================================
  !
  ! J M Borrero
  ! Jan 18, 2012
  !
  USE CONS_PARAM, ONLY: DP
  USE FORWARD_PARAM, ONLY: JL, JU, LINE_L0, SL, SU, LL, LU
  USE CODE_MODES, ONLY: INPUTFILE
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_DP, ALLOCATE_1D_IP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP),       ALLOCATABLE     :: LAMBDA_SC(:)         ! Lambda Split Component
  REAL(DP),       ALLOCATABLE     :: STRENGTH_SC(:)       ! Strength of the Split Component
  INTEGER,        ALLOCATABLE     :: TYPE_SC(:)           ! Type of the Split Component
  INTEGER                         :: NPI, NSR, NSB, NTOT
  !
  PUBLIC :: LAMBDA_SC
  PUBLIC :: STRENGTH_SC
  PUBLIC :: TYPE_SC
  PUBLIC :: NTOT
  PUBLIC :: ZEEMAN_COMP
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! zeeman_comp
  !
  !------------------------------------------------
  !
  SUBROUTINE ZEEMAN_COMP(I)
    !
    INTEGER                 :: I
    !
    REAL(DP)                :: GLOW, GUPP, GEFF
    REAL(DP)                :: SLOW, SUPP, LLOW, LUPP, JLOW, JUPP
    REAL(DP)                :: DELTALB
    REAL(DP)                :: SUM_PI, SUM_SB, SUM_SR
    INTEGER                 :: J, K, COUNT, NMJL, NMJU, DELTAM, DELTAJ
    REAL(DP), ALLOCATABLE   :: ML(:), MU(:)
    !
    ! Number of transitions for each type
    NPI = INT(2*MIN(JL(I),JU(I))+1)    ! Pi-components
    NSR = INT(JL(I)+JU(I))             ! Sigma-red components
    NSB = INT(JL(I)+JU(I))             ! Sigma-blue components 
    NTOT = NPI + NSR + NSB
    DELTAJ = INT(JU(I)-JL(I))
    ! Selection rule in quantum number J
    IF (ABS(DELTAJ).GT.1) THEN
       PRINT*,'Selection rules imply that DELTAJ=-1,0,+1. However,'
       PRINT*,'in '//TRIM(INPUTFILE)//' the spectral line', LINE_L0(I)
       PRINT*,'does not verify such rule. Please correct the electronic configuration. STOP.'
       STOP
    ENDIF
    !
    ! We use double precision because this will be employed to determine wavelength shifts
    !
    SLOW=DBLE(SL(I))
    SUPP=DBLE(SU(I))
    LLOW=DBLE(LL(I))
    LUPP=DBLE(LU(I))
    JLOW=DBLE(JL(I))
    JUPP=DBLE(JU(I))
    !
    ! Lande factor lower level
    IF (JLOW.NE.0) THEN
       GLOW=1.5D0+(SLOW*(SLOW+1.0D0)-LLOW*(LLOW+1.0D0))/(2.0D0*JLOW*(JLOW+1.0D0))
    ELSE
       GLOW=0.0D0
    ENDIF
    ! Lande factor upper level
    IF (JUPP.NE.0) THEN
       GUPP=1.5D0+(SUPP*(SUPP+1.0D0)-LUPP*(LUPP+1.0D0))/(2.0D0*JUPP*(JUPP+1.0D0))
    ELSE
       GUPP=0.0D0
    ENDIF
    !
    ! Effective Lande factor
    GEFF=0.5D0*(GLOW+GUPP)+0.25D0*(GLOW-GUPP)*(JLOW*(JLOW+1.0D0)-JUPP*(JUPP+1.0D0))
    ! Allocate array for wavelength shifts
    !ALLOCATE( LAMBDA_SC(NTOT), STRENGTH_SC(NTOT), TYPE_SC(NTOT))
    CALL ALLOCATE_1D_DP(LAMBDA_SC, NTOT, 'ALLOCATE lambda_sc IN magsplit.f90')
    CALL ALLOCATE_1D_DP(STRENGTH_SC, NTOT, 'ALLOCATE strength_sc IN magsplit.f90')
    CALL ALLOCATE_1D_IP(TYPE_SC, NTOT, 'ALLOCATE type_sc IN magsplit.f90')
    !
! TO BE DEALLOCATED: LAMBDA_SC
! TO BE DEALLOCATED: STRENGTH_SC
! TO BE DEALLOCATED: TYPE_SC
    !
    ! Allocate array for number of M-sublevels that split out of one J level
    NMJL=INT(2*JL(I)+1)
    NMJU=INT(2*JU(I)+1)
    !ALLOCATE (ML(NMJL), MU(NMJU))
    CALL ALLOCATE_1D_DP(ML, NMJL, 'ALLOCATE ml IN magsplit.f90')
    CALL ALLOCATE_1D_DP(MU, NMJU, 'ALLOCATE mu IN magsplit.f90')
    !
    ! Create M quantum numbers
    DO J=1,NMJL
       ML(J)=DBLE(J-1.0D0-JLOW)
    ENDDO
    DO J=1,NMJU
       MU(J)=DBLE(J-1.0D0-JUPP)
    ENDDO
    ! Start Loop to determine allowed transitions and their wavelength shift and strength
    COUNT = 0
    DO J=1,NMJL
       DO K=1,NMJU
          DELTAM=INT(MU(K)-ML(J))
          IF (ABS(DELTAM).LE.1) THEN
             COUNT=COUNT+1
             ! Units are mA/Gauss
             DELTALB=REAL(4.66685D-10*(ML(J)*GLOW-MU(K)*GUPP)*LINE_L0(I)**2)
             TYPE_SC(COUNT)=DELTAM
             LAMBDA_SC(COUNT)=DELTALB
             SELECT CASE (DELTAJ)
                !
             CASE(-1)
                SELECT CASE (DELTAM)
                CASE(-1)
                   STRENGTH_SC(COUNT) = (JLOW+MU(K))*(JUPP+MU(K)+2.0D0)
                CASE(0)
                   STRENGTH_SC(COUNT) = 2.0D0*(JLOW**2.0D0-MU(K)**2.0D0)
                CASE(1)
                   STRENGTH_SC(COUNT) = (JLOW-MU(K))*(JUPP-MU(K)+2.0D0)
                END SELECT
                !
             CASE(0)
                SELECT CASE (DELTAM)
                CASE(-1)
                   STRENGTH_SC(COUNT) = (JUPP-MU(K))*(JUPP+MU(K)+1.0D0)
                CASE(0)
                   STRENGTH_SC(COUNT) = 2.0D0*MU(K)**2.0D0
                CASE(1)
                   STRENGTH_SC(COUNT) = (JUPP+MU(K))*(JUPP-MU(K)+1.0D0)
                END SELECT
             CASE(1)
                SELECT CASE (DELTAM)
                CASE(-1)
                   STRENGTH_SC(COUNT) = (JUPP-MU(K))*(JLOW-MU(K))
                CASE(0)
                   STRENGTH_SC(COUNT) = 2.0D0*(JUPP**2.0D0-MU(K)**2.0D0)
                CASE(1)
                   STRENGTH_SC(COUNT) = (JUPP+MU(K))*(JLOW+MU(K))
                END SELECT
             END SELECT
          ENDIF
       ENDDO
    ENDDO
    !
    DEALLOCATE( ML, MU )
    !
    !
    ! Normalize the strength of the Zeeman pattern individually for PI, SIGMA_BLUE and SIGMA_RED components
    ! 
    SUM_PI=0.0D0
    SUM_SR=0.0D0
    SUM_SB=0.0D0
    DO J=1,NTOT
       SELECT CASE (TYPE_SC(J))
       CASE(-1)
          SUM_SR=SUM_SR+STRENGTH_SC(J)
       CASE(0)
          SUM_PI=SUM_PI+STRENGTH_SC(J)
       CASE(1)
          SUM_SB=SUM_SB+STRENGTH_SC(J)
       END SELECT
    ENDDO
    !
    DO J=1,NTOT
       SELECT CASE (TYPE_SC(J))
       CASE(-1)
          STRENGTH_SC(J)=STRENGTH_SC(J)/SUM_SR
       CASE(0)
          STRENGTH_SC(J)=STRENGTH_SC(J)/SUM_PI
       CASE(1)
          STRENGTH_SC(J)=STRENGTH_SC(J)/SUM_SB
       END SELECT
    ENDDO
    !
  END SUBROUTINE ZEEMAN_COMP
  !
  !================================================
  !
END MODULE MAGSPLIT
!
