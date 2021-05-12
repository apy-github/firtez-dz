!
MODULE HYDROSTATIC_TOOLS
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP, GRAV, MAMU, KBOL
  USE CHEMICAL_EQUILIBRIUM, ONLY: GET_RHO
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER       :: K
  REAL(DP)      :: DZ
  REAL(DP)      :: NHYD,NELEC,PATOM,PELEC,DENS,MOLECW
  !
  REAL(DP)      ::  T0
  REAL(DP)      ::  T1
  REAL(DP)      ::  U0
  REAL(DP)      ::  U1
  !
  PUBLIC :: RK4_INTEGRATION_LOG_BOTTOM
  PUBLIC :: RK4_INTEGRATION_LOG
  PUBLIC :: RK4_INTEGRATION_LOG_LTAUM05
  PRIVATE :: RK4_SOLVE_LOG
  PRIVATE :: RK4_INTEGRATION_LOG_MED
  !
  PRIVATE :: K
  PRIVATE :: DZ
  PRIVATE :: NHYD,NELEC,PATOM,PELEC,DENS,MOLECW
  !
  PRIVATE ::  T0
  PRIVATE ::  T1
  PRIVATE ::  U0
  PRIVATE ::  U1
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! rk4_solve_log
  ! rk4_integration_log_bottom
  ! rk4_integration_log
  ! rk4_integration_log_ltaum05
  ! rk4_integration_log_med
  ! rk4_log
  !
  !------------------------------------------------
  !
  SUBROUTINE RK4_SOLVE_LOG(ZI, PI, GI, RCTE, MH, TEM, NEL, DPDZ)
    !
    REAL(DP), INTENT(IN)    :: ZI
    REAL(DP), INTENT(IN)    :: PI
    REAL(DP), INTENT(IN)    :: GI
    REAL(DP), INTENT(IN)    :: RCTE
    REAL(DP), INTENT(IN)    :: MH
    REAL(DP), INTENT(IN)    :: TEM
    !
    REAL(DP), INTENT(INOUT) :: DPDZ, NEL
    REAL(DP)                :: INHYD,IPATOM,IPELEC,IDENS,IMOLECW

    ! Calculate the gas density provided the equation of state:
    CALL GET_RHO(TEM,PI &
        ,INHYD,NEL,IPATOM,IPELEC,IDENS,IMOLECW)
    !
    DPDZ=-IMOLECW*MH*GI/(RCTE*TEM)
    !
  END SUBROUTINE RK4_SOLVE_LOG
  !
  !------------------------------------------------
  !
  SUBROUTINE RK4_INTEGRATION_LOG_BOTTOM(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ)
    !
    INTEGER, INTENT(IN)                     :: NZ
    REAL(SP), INTENT(IN), DIMENSION(NZ)     :: TEM
    REAL(SP), INTENT(INOUT), DIMENSION(NZ)  :: PG, RHO, PEL, MW
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: ZZ
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: VNHYD
    !
    REAL(DP), DIMENSION(NZ)                 :: PG8, RHO8, PEL8, MW8, T8
    !
    PG8(:)=0.D0
    RHO8(:)=0.D0
    PEL8(:)=0.D0
    MW8(:)=0.D0
    T8(:)=0.D0
    DZ=0.D0
    NHYD=0.D0
    NELEC=0.D0
    PATOM=0.D0
    PELEC=0.D0
    DENS=0.D0
    MOLECW=0.D0
    !
    T0=0.D0
    T1=0.D0
    U0=0.D0
    U1=0.D0
    !
    T8(:)=DBLE(TEM)
    !
    PG8(1)=DLOG(DBLE(PG(1)))
    RHO8(1)=DBLE(RHO(1))
    PEL8(1)=DBLE(PEL(1))
    MW8(1)=DBLE(MW(1))
    !
    T0 = ZZ(1)*1000.*100.
    U0 = PG8(1)
    !
    CALL GET_RHO(T8(1),DEXP(U0)&
        ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
    RHO8(1)=DENS
    PEL8(1)=PELEC
    MW8(1)=MOLECW
    VNHYD(1)=NHYD
    !
    DO K=2,NZ
      !
      DZ=DBLE(ZZ(K)-ZZ(K-1))*1000.0D0*100.0D0
      !
      T1=T0+DZ
      !
      CALL RK4_LOG(T0, U0, DZ, T8(K-1) &
          , T8(K), NELEC, RK4_SOLVE_LOG, U1)
      !
      CALL GET_RHO(T8(K),DEXP(U1)&
          ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
      !
      RHO8(K)=DENS
      PEL8(K)=PELEC
      MW8(K)=MOLECW
      PG8(K)=U1
      VNHYD(K)=NHYD
      U0=U1
      T0=T1
      !
    ENDDO
    PG(:)=REAL(DEXP(PG8(:)))
    RHO(:)=REAL(RHO8(:))
    PEL(:)=REAL(PEL8(:))
    MW(:)=REAL(MW8(:))
    !
  END SUBROUTINE RK4_INTEGRATION_LOG_BOTTOM
  !
  !------------------------------------------------
  !
  SUBROUTINE RK4_INTEGRATION_LOG(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ)
    !
    INTEGER, INTENT(IN)                     :: NZ
    REAL(SP), INTENT(IN), DIMENSION(NZ)     :: TEM
    REAL(SP), INTENT(INOUT), DIMENSION(NZ)  :: PG, RHO, PEL, MW
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: ZZ
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: VNHYD
    !
    REAL(DP), DIMENSION(NZ)                 :: PG8, RHO8, PEL8, MW8, T8
    !
    PG8(:)=0.D0
    RHO8(:)=0.D0
    PEL8(:)=0.D0
    MW8(:)=0.D0
    T8(:)=0.D0
    DZ=0.D0
    NHYD=0.D0
    NELEC=0.D0
    PATOM=0.D0
    PELEC=0.D0
    DENS=0.D0
    MOLECW=0.D0
    !
    T0=0.D0
    T1=0.D0
    U0=0.D0
    U1=0.D0
    !
    T8(:)=DBLE(TEM)
    !
    PG8(NZ)=DLOG(DBLE(PG(NZ)))
    RHO8(NZ)=DBLE(RHO(NZ))
    PEL8(NZ)=DBLE(PEL(NZ))
    NELEC=PEL8(NZ)/KBOL/T8(NZ)
    MW8(NZ)=DBLE(MW(NZ))
    !
    T0 = ZZ(NZ)*1000.0D0*100.0D0
    U0 = PG8(NZ)
    !
    CALL GET_RHO(T8(NZ),DEXP(U0)&
        ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
    RHO8(NZ)=DENS
    PEL8(NZ)=PELEC
    MW8(NZ)=MOLECW
    VNHYD(NZ)=NHYD
    !
    DO K=NZ-1,1,-1
      !
      DZ=DBLE(ZZ(K)-ZZ(K+1))*1000.D0*100.D0
      !
      T1=T0+DZ
      !
      CALL RK4_LOG(T0, U0, DZ, T8(K+1) &
          , T8(K), NELEC, RK4_SOLVE_LOG, U1)
      !
      CALL GET_RHO(T8(K),DEXP(U1)&
          ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
      !
      RHO8(K)=DENS
      PEL8(K)=PELEC
      MW8(K)=MOLECW
      PG8(K)=U1
      VNHYD(K)=NHYD
      U0=U1
      T0=T1
      !
    ENDDO
    PG(:)=REAL(DEXP(PG8(:)))
    RHO(:)=REAL(RHO8(:))
    PEL(:)=REAL(PEL8(:))
    MW(:)=REAL(MW8(:))
    !
  END SUBROUTINE RK4_INTEGRATION_LOG
  !
  !------------------------------------------------
  !
  SUBROUTINE RK4_INTEGRATION_LOG_LTAUM05(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ,TAU5)
    !
    INTEGER, INTENT(IN)                     :: NZ
    REAL(SP), INTENT(IN), DIMENSION(NZ)     :: TEM, TAU5
    REAL(SP), INTENT(INOUT), DIMENSION(NZ)  :: PG, RHO, PEL, MW
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: ZZ
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: VNHYD
    !
    INTEGER                                 :: IZ
    !
    IZ=MINLOC(ABS(TAU5(:)),1)
    !
    CALL RK4_INTEGRATION_LOG_MED(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ,IZ)
    !
  END SUBROUTINE RK4_INTEGRATION_LOG_LTAUM05
  !
  !------------------------------------------------
  !
  SUBROUTINE RK4_INTEGRATION_LOG_MED(NZ,TEM,PG,RHO,PEL,MW,VNHYD,ZZ,IZ)
    !
    !
    INTEGER, INTENT(IN)                     :: NZ, IZ
    REAL(SP), INTENT(IN), DIMENSION(NZ)     :: TEM
    REAL(SP), INTENT(INOUT), DIMENSION(NZ)  :: PG, RHO, PEL, MW
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: ZZ
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)  :: VNHYD
    !
    REAL(DP), DIMENSION(NZ)                 :: PG8, RHO8, PEL8, MW8, T8
    !
    PG8(:)=0.D0
    RHO8(:)=0.D0
    PEL8(:)=0.D0
    MW8(:)=0.D0
    T8(:)=0.D0
    DZ=0.D0
    NHYD=0.D0
    NELEC=0.D0
    PATOM=0.D0
    PELEC=0.D0
    DENS=0.D0
    MOLECW=0.D0
    !
    T0=0.D0
    T1=0.D0
    U0=0.D0
    U1=0.D0
    !
    T8(:)=DBLE(TEM)
    !
    PG8(IZ)=DLOG(DBLE(PG(IZ)))
    RHO8(IZ)=DBLE(RHO(IZ))
    PEL8(IZ)=DBLE(PEL(IZ))
    MW8(IZ)=DBLE(MW(IZ))
    NELEC=PEL8(IZ)/KBOL/T8(IZ)
    !
    T0 = ZZ(IZ)*1000.D0*100.D0
    U0 = PG8(IZ)
    !
    CALL GET_RHO(T8(IZ),DEXP(U0)&
        ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
    RHO8(IZ)=DENS
    PEL8(IZ)=PELEC
    MW8(IZ)=MOLECW
    VNHYD(IZ)=NHYD
    !
    ! WE FIRST GO UPWARDS
    DO K=IZ+1,NZ,1
      !
      DZ=DBLE(ZZ(K)-ZZ(K-1))*1000.D0*100.D0
      !
      T1=T0+DZ
      !
      CALL RK4_LOG(T0, U0, DZ, T8(K-1) &
          , T8(K), NELEC, RK4_SOLVE_LOG, U1)
      !
      CALL GET_RHO(T8(K),DEXP(U1)&
          ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
      !
      RHO8(K)=DENS
      PEL8(K)=PELEC
      MW8(K)=MOLECW
      PG8(K)=U1
      VNHYD(K)=NHYD
      U0=U1
      T0=T1
      !
    ENDDO
    !
    ! WE THEN GO DOWNWARDS:
    T0 = ZZ(IZ)*1000.0D0*100.0D0
    U0 = PG8(IZ)
    NELEC=PEL8(IZ)/KBOL/T8(IZ)
    !
    DO K=IZ-1,1,-1
      !
      DZ=DBLE(ZZ(K)-ZZ(K+1))*1000.D0*100.D0
      !
      T1=T0+DZ
      !
      CALL RK4_LOG(T0, U0, DZ, T8(K+1) &
          , T8(K), NELEC, RK4_SOLVE_LOG, U1)
      !
      CALL GET_RHO(T8(K),DEXP(U1)&
          ,NHYD,NELEC,PATOM,PELEC,DENS,MOLECW)
      !
      RHO8(K)=DENS
      PEL8(K)=PELEC
      MW8(K)=MOLECW
      PG8(K)=U1
      VNHYD(K)=NHYD
      U0=U1
      T0=T1
      !
    ENDDO
    !
    PG(:)=REAL(DEXP(PG8(:)))
    RHO(:)=REAL(RHO8(:))
    PEL(:)=REAL(PEL8(:))
    MW(:)=REAL(MW8(:))
    !
  END SUBROUTINE RK4_INTEGRATION_LOG_MED
  !
  !-----------------------------------------------
  !
  SUBROUTINE RK4_LOG(IT0,IY0,IDT,TM0,TM1,NELEC,ESUBROUTINE,IYO)
    !
    !Following wikipedia:
    !
    REAL(DP), INTENT(IN)                :: IT0
    REAL(DP), INTENT(IN)                :: IY0
    REAL(DP), INTENT(IN)                :: IDT
    REAL(DP), INTENT(IN)                :: TM0
    REAL(DP), INTENT(IN)                :: TM1
    REAL(DP), INTENT(INOUT)             :: NELEC
    REAL(DP), INTENT(INOUT)             :: IYO
    !
    EXTERNAL :: ESUBROUTINE
    !
    REAL(DP)     :: K1, K2, K3, K4
    REAL(DP)     :: IT1, IT2, IT3
    REAL(DP)     :: IY1, IY2, IY3
    REAL(DP)     :: DPDZ
    REAL(DP)     :: TEM
    !
    ! For a RK4, we need four evaluations of the function:
    TEM = TM0
    CALL ESUBROUTINE(IT0,IY0,GRAV,KBOL,MAMU,TEM,NELEC,DPDZ)
    K1 = DPDZ * IDT
    !
    IT1=IT0+IDT/2.0D0
    IY1=IY0+K1/2.0D0
    TEM = TM0 / IDT * (IT0+IDT-IT1) + TM1 / IDT * (IT1-IT0)
    CALL ESUBROUTINE(IT1,IY1,GRAV,KBOL,MAMU,TEM,NELEC,DPDZ)
    K2 = DPDZ * IDT
    !
    IT2=IT0+IDT/2.0D0
    IY2=IY0+K2/2.0D0
    !TEM = TM0 / IDT * (IT0+IDT-IT1) + TM1 / IDT * (IT1-IT0)
    ! It is the same point so it is not necessary to repeat it
    CALL ESUBROUTINE(IT2,IY2,GRAV,KBOL,MAMU,TEM,NELEC,DPDZ)
    K3 = DPDZ * IDT
    !
    IT3=IT0+IDT
    IY3=IY0+K3
    TEM = TM1
    CALL ESUBROUTINE(IT3,IY3,GRAV,KBOL,MAMU,TEM,NELEC,DPDZ)
    K4 = DPDZ * IDT
    !
    IYO = IY0 + 1.0D0 / 6.0D0 * ( K1 + 2.0D0 * K2 + 2.0D0 * K3 + K4 )
    !
  END SUBROUTINE RK4_LOG
  !
  !================================================
  !
END MODULE 
!
