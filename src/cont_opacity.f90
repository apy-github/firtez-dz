!
MODULE CONT_OPACITY
  !
  !================================================
  !
  ! J M Borrero
  ! June 14th, 2017
  ! KIS Freiburg
  ! Adapted from routines by Andres Asensio at the IAC.
  ! These routines also calculate the derivatives at constant density and gas 
  ! pressure and to work with the number of electrons and hydrogen atoms per cm^3
  !
  USE CONS_PARAM, ONLY: DP, SP, KBOL, LIGHT, HPLA, MELE &
      , QELE, EVOLT, DPI
  USE DERIVVAR
  USE CHEMICAL_EQUILIBRIUM, ONLY: ION_CALC_TEMPRHO
  USE CODE_MODES, ONLY: MRESPFUNCT
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP)               :: HMINUS_FF, DHMINUS_FF_RHO, DHMINUS_FF_PG
  REAL(DP)               :: HMINUS_BF, DHMINUS_BF_RHO, DHMINUS_BF_PG
  REAL(DP)               :: THOMSON, DTHOMSON_PG, DTHOMSON_RHO
  REAL(DP)               :: RAYLEIGH_H
  REAL(DP)               :: HYD, DHYD_PG, DHYD_RHO
  ! 
  REAL(DP)               :: DHMINUS_FFDPG_TEMP, DHMINUS_BFDPG_TEMP, DTHOMSONDPG_TEMP
  REAL(DP)               :: DHMINUS_FFDRHO_TEMP, DHMINUS_BFDRHO_TEMP, DTHOMSONDRHO_TEMP
  ! 
  PUBLIC :: OPAC_CONTINUUM
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! opac_hminus_ff
  ! opac_hminus_bf
  ! opac_thomson
  ! opac_rayleigh_h
  ! opac_hydrogen
  ! opac_continuum
  !
  !------------------------------------------------
  !
  SUBROUTINE OPAC_HMINUS_FF(T,NE,LAMBDA0,OPAC,DOPACDT_RHO,DOPACDT_PG)
    ! This routine calculates the continuum opacity per hydrogen atom due to free-free 
    ! transitions in negative hydrogen atom. Also calculates derivatives with respect to 
    ! temperature at constant density and gas pressure.
    ! Andres Asension cites:  John 1989 A&A 193, 189
    ! T      : temperature in K (1400 < T < 100080)
    ! LAMBDA0: wavelength in A (greater than 1880 A)
    ! NE     : number of electrons per cm^3
    IMPLICIT NONE
    ! Input/Output
    REAL(DP),   INTENT(IN)                  :: T
    REAL(DP),   INTENT(IN)                  :: NE, LAMBDA0
    REAL(DP),   INTENT(OUT), OPTIONAL       :: OPAC, DOPACDT_RHO, DOPACDT_PG
    ! Internal
    REAL(DP)    , DIMENSION(6)              :: A1, B1, C1, D1, E1, F1, COM1
    REAL(DP)    , DIMENSION(4)              :: A2, B2, C2, D2, E2, F2, COM2
    REAL(DP)                                :: THETA, DTHETA, LAMBDA0MIC, PART1, PART2
    ! Reset values
    HMINUS_FF = 0D0
    DHMINUS_FF_RHO = 0D0
    DHMINUS_FF_PG = 0D0
    ! Contants for polynomials
    A1 = (/0D0,2483.346D0,-3449.889D0,2200.04D0,-696.271D0,88.283D0/)
    B1 = (/0D0,285.827D0,-1158.382D0,2427.719D0,-1841.4D0,444.517D0/)
    C1 = (/0D0,-2054.291D0,8746.523D0,-13651.105D0,8624.97D0,-1863.864D0/)
    D1 = (/0D0,2827.776D0,-11485.632D0,16755.524D0,-10051.53D0,2095.288D0/)
    E1 = (/0D0,-1341.537D0,5303.609D0,-7510.494D0,4400.067D0,-901.788D0/)
    F1 = (/0D0,208.952D0,-812.939D0,1132.738D0,-655.02D0,132.985D0/)
    A2 = (/518.1021D0,473.2636D0,-482.2089D0,115.5291D0/)
    B2 = (/-734.8666D0,1443.4137D0,-737.1616D0,169.6374D0/)
    C2 = (/1021.1775D0,-1977.3395D0,1096.8827D0,-245.649D0/)
    D2 = (/-479.0721D0,922.3575D0,-521.1341D0,114.243D0/)
    E2 = (/93.1373D0,-178.9275D0,101.7963D0,-21.9972D0/)
    F2 = (/-6.4285D0,12.36D0,-7.0571D0,1.5097D0/)
    !
    IF (LAMBDA0.LT.1800D0) THEN 
       PRINT*,'Error in HMINUS_FF. Wavelength must be larger than 1800 A'
       STOP
    ENDIF
    ! Wavelength to microns
    LAMBDA0MIC = LAMBDA0/1D4
    THETA = 5040D0/DBLE(T)
    DTHETA=-THETA/DBLE(T)
    !
    IF (LAMBDA0MIC.LT.0.3645D0) THEN
       COM2 = A2*LAMBDA0MIC**2D0 + B2 + C2/LAMBDA0MIC + D2/LAMBDA0MIC**2D0 + E2/LAMBDA0MIC**3D0 + F2/LAMBDA0MIC**4D0
       PART1 = COM2(1)*THETA + COM2(2)*THETA**1.5D0 + COM2(3)*THETA**2D0 + COM2(4)*THETA**2.5D0
       PART2 = (COM2(1) + 1.5D0*COM2(2)*THETA**0.5D0 + 2D0*COM2(3)*THETA + 2.5D0*COM2(4)*THETA**1.5D0)*DTHETA
    ELSE
       COM1 = A1*LAMBDA0MIC**2D0 + B1 + C1/LAMBDA0MIC + D1/LAMBDA0MIC**2D0 + E1/LAMBDA0MIC**3D0 + F1/LAMBDA0MIC**4D0
       PART1 = COM1(1)*THETA + COM1(2)*THETA**1.5D0 + COM1(3)*THETA**2D0 + COM1(4)*THETA**2.5D0 + &
            COM1(5)*THETA**3D0 + COM1(6)*THETA**3.5D0
       PART2= (COM1(1) + 1.5D0*COM1(2)*THETA**0.5D0 + 2D0*COM1(3)*THETA + &
            2.5D0*COM1(4)*THETA**1.5D0 + 3D0*COM1(5)*THETA**2D0 + 3.5D0*COM1(6)*THETA**2.5D0)*DTHETA
    ENDIF
    !
    HMINUS_FF = 1D-29*PART1*KBOL*NE*DBLE(T)
    !
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      DHMINUS_FF_RHO = HMINUS_FF*(DNELECDTEMP_RHO/NE+1D0/DBLE(T)+(PART2/PART1))
      DHMINUS_FF_PG =  HMINUS_FF*(DNELECDTEMP_PG/NE+1D0/DBLE(T)+(PART2/PART1))
      ! 
      DHMINUS_FFDPG_TEMP=HMINUS_FF/NE*DNELECDPG_TEMP
      DHMINUS_FFDRHO_TEMP=HMINUS_FF/NE*DNELECDRHO_TEMP
      ! 
      IF (PRESENT(OPAC).AND.PRESENT(DOPACDT_RHO).AND.PRESENT(DOPACDT_PG)) THEN
         OPAC =  HMINUS_FF
         DOPACDT_RHO = DHMINUS_FF_RHO 
         DOPACDT_PG = DHMINUS_FF_PG
      ENDIF
    ENDIF
    !
  END SUBROUTINE OPAC_HMINUS_FF
  !
  !------------------------------------------------
  !
  SUBROUTINE OPAC_HMINUS_BF(T,NE,LAMBDA0,OPAC,DOPACDT_RHO,DOPACDT_PG)
    ! This routine calculates the continuum opacity per hydrogen atom due to bound-free 
    ! transitions in negative hydrogen atom. Also calculates derivatives with respect to 
    ! temperature at constant density and gas pressure.
    ! Andres Asensio cites:  John 1989 A&A 193, 189
    ! T      : temperature in K (1400 < T < 100080)
    ! LAMBDA0: wavelength in A (greater than 1880 A)
    ! NE     : number of electrons per cm^3
    IMPLICIT NONE
    ! Input/Output
    REAL(DP),   INTENT(IN)                  :: T
    REAL(DP),   INTENT(IN)                  :: NE, LAMBDA0
    REAL(DP),   INTENT(OUT), OPTIONAL       :: OPAC, DOPACDT_RHO, DOPACDT_PG
    ! Internal
    REAL(DP),   PARAMETER                   :: LAMBDAP = 1.6419D0
    REAL(DP),   PARAMETER                   :: CTE = 0.75D-18 
    REAL(DP),   PARAMETER                   :: ALPHA = (HPLA*LIGHT/KBOL)*1D4
    REAL(DP),   DIMENSION(6)                :: CC
    REAL(DP)                                :: LAMBDA0MIC, COM, SIGMA, PART, DPART
    ! Reset values
    HMINUS_BF = 0D0
    DHMINUS_BF_RHO = 0D0
    DHMINUS_BF_PG = 0D0
    !
    CC = (/152.519D0,49.534D0,-118.858D0,92.536D0,-34.194D0,4.982D0/)
    ! Wavelength to microns
    LAMBDA0MIC = LAMBDA0/1D4
    !
    IF (LAMBDA0MIC.LT.LAMBDAP) THEN
       COM=1D0/LAMBDA0MIC-1D0/LAMBDAP
       SIGMA = CC(1) + CC(2)*COM**0.5D0 + CC(3)*COM + CC(4)*COM**1.5D0 + CC(5)*COM**2D0 + CC(6)*COM**2.5D0
       SIGMA = CTE*SIGMA*LAMBDA0MIC**3D0*COM**1.5D0
       PART = T**(-2.5D0)*DEXP(ALPHA/(T*LAMBDAP))*(1D0-DEXP(-ALPHA/(T*LAMBDA0MIC)))
       DPART = (PART/T)*(-2.5D0+ALPHA*COM/T-(ALPHA/(LAMBDA0MIC*T))*(1D0/(1D0-DEXP(-ALPHA/(LAMBDA0MIC*T)))))
       HMINUS_BF = PART*SIGMA*NE*KBOL*T
       IF (MRESPFUNCT.EQV..TRUE.) THEN
         DHMINUS_BF_RHO = HMINUS_BF*(DPART/PART+1D0/T+DNELECDTEMP_RHO/NE)
         DHMINUS_BF_PG = HMINUS_BF*(DPART/PART+1D0/T+DNELECDTEMP_PG/NE)
         ! 
         DHMINUS_BFDPG_TEMP=HMINUS_BF/NE*DNELECDPG_TEMP
         DHMINUS_BFDRHO_TEMP=HMINUS_BF/NE*DNELECDRHO_TEMP
       ENDIF
       ! 
    ELSE
       HMINUS_BF=0D0
       DHMINUS_BF_RHO = 0D0 
       DHMINUS_BF_PG = 0D0
       DHMINUS_BFDPG_TEMP=0D0
       DHMINUS_BFDRHO_TEMP=0D0
    ENDIF
    !
    IF (PRESENT(OPAC).AND.PRESENT(DOPACDT_RHO).AND.PRESENT(DOPACDT_PG)) THEN
       OPAC =  HMINUS_BF
       DOPACDT_RHO = DHMINUS_BF_RHO 
       DOPACDT_PG = DHMINUS_BF_PG
    ENDIF
    !
  END SUBROUTINE OPAC_HMINUS_BF
  !
  !------------------------------------------------
  !
  !----------------------------------------------
  !----------------------------------------------
  SUBROUTINE OPAC_THOMSON(NE,NH,OPAC,DOPACDT_RHO,DOPACDT_PG)
    ! Calculates the continuum absorption coefficient per hydrogen atom 
    ! due to Thomson scattering (scattering with free electrons)
    ! Also calculates the partial derivatives with respect to the
    ! temperature at constant density and gas pressure
    IMPLICIT NONE
    ! Input/Output
    REAL(DP), INTENT(IN)            :: NE, NH
    REAL(DP), INTENT(OUT), OPTIONAL :: OPAC, DOPACDT_RHO, DOPACDT_PG
    ! Internal
    REAL(DP), PARAMETER             :: CT = (8D0*DPI/3D0)*((QELE/LIGHT)**4D0)/MELE**2D0
    ! Reset values
    THOMSON = 0D0
    DTHOMSON_PG = 0D0
    DTHOMSON_RHO = 0D0
    !
    THOMSON = CT*(NE/NH)
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      DTHOMSON_PG = THOMSON*(DNELECDTEMP_PG/NE-DNHYDDTEMP_PG/NH)
      DTHOMSON_RHO = THOMSON*(DNELECDTEMP_RHO/NE-DNHYDDTEMP_RHO/NH)
      DTHOMSONDPG_TEMP=THOMSON/NE*DNELECDPG_TEMP-THOMSON/NH*DNHYDDPG_TEMP
      DTHOMSONDRHO_TEMP=THOMSON/NE*DNELECDRHO_TEMP-THOMSON/NH*DNHYDDRHO_TEMP
      !
      IF (PRESENT(OPAC).AND.PRESENT(DOPACDT_RHO).AND.PRESENT(DOPACDT_PG)) THEN
         OPAC =  THOMSON
         DOPACDT_RHO = DTHOMSON_RHO 
         DOPACDT_PG = DTHOMSON_PG
      ENDIF
    ENDIF
    !
  END SUBROUTINE OPAC_THOMSON
  !
  !------------------------------------------------
  !
  SUBROUTINE OPAC_RAYLEIGH_H(LAMBDA0)
    ! Continuum absorption coefficient per hydrogen atom
    ! due to Rayleigh scattering by hydrogen
    ! Comments by Andres Asensio: Landi Degl'Innocenti 1976, A&ASS, 25, 379 
    ! (he quotes Dalgarno (1962) but the polynomial has a mistake)
    ! Derivatives with respect to T are zero
    IMPLICIT NONE
    ! Input/Output
    REAL(DP), INTENT(IN)          :: LAMBDA0
    ! Internal
    REAL(DP), DIMENSION(3)        :: CC
    ! Reset values
    RAYLEIGH_H = 0D0
    !
    CC = (/5.799D-13,1.422D-6,2.784D0/)
    !
    RAYLEIGH_H = (CC(1)+(CC(2)+CC(3)/LAMBDA0**2D0) / LAMBDA0**2D0) / LAMBDA0**4D0
    !
  END SUBROUTINE OPAC_RAYLEIGH_H
  !
  !------------------------------------------------
  !
  SUBROUTINE OPAC_HYDROGEN(T,LAMBDA0,OPAC,DOPACDT_RHO,DOPACDT_PG)
    ! Conrinuum absorption coefficient per hydrogen atom
    ! due to bound-bound and bound-free transitions
    ! by neutral hydrogen
    ! Comments by Andres Asensio: Landi Degl'Innocenti 1976, A&ASS, 25, 379
    ! These equations assume that the partition function for hydrogen is 2
    ! and this is only valid for T < 12000 K. Above this, we have only
    ! an approximation
    IMPLICIT NONE
    ! Input/Output 1E-8*13.595D0*EVOLT/(HPLA*LIGHT)
    REAL(DP), INTENT(IN)            :: T
    REAL(DP), INTENT(IN)            :: LAMBDA0
    REAL(DP), INTENT(OUT), OPTIONAL :: OPAC, DOPACDT_RHO, DOPACDT_PG
    ! Internal
    REAL(DP), PARAMETER             :: IOPOT = 13.595 ! eV
    REAL(DP), PARAMETER             :: C1 = 1E-8*IOPOT*EVOLT/(HPLA*LIGHT)
    ! C1 = 1.09651067903121578E-003
    REAL(DP), PARAMETER             :: C2 = (1E-8)**3D0*64D0*DPI**4D0*MELE*QELE**10D0/(3D0*SQRT(3D0)*LIGHT**4D0*HPLA**6D0)
    ! C2 = 1.04490915480325020E-026
    REAL(DP), PARAMETER             :: C3 = HPLA*LIGHT/(KBOL*1E-8)
    ! C3 = 143886436.07434255
    REAL(DP), PARAMETER             :: C4 = IOPOT*EVOLT/KBOL
    ! C4 = 157773.01682218799
    INTEGER                         :: N0, I
    REAL(DP)                        :: SUM, DSUM
    ! Reset values
    HYD = 0D0
    DHYD_PG = 0D0
    DHYD_RHO = 0D0
    ! Andres' implementation
    N0=1+FLOOR(SQRT(C1*LAMBDA0))
    SUM=0D0
    DSUM=0D0
    IF (N0.LE.8) THEN
       SUM = DEXP(C4/(T*DBLE(N0)**2D0))*DBLE(N0)**(-3D0)
       DSUM = -C4*DEXP(C4/(T*DBLE(N0)**2D0))*DBLE(N0)**(-5D0)*T**(-2D0)
       DO I=N0+1,8
          SUM = SUM+DEXP(C4/(T*DBLE(I)**2D0))*DBLE(I)**(-3D0)
          DSUM = DSUM - C4*DEXP(C4/(T*DBLE(I)**2D0))*DBLE(I)**(-5D0)*T**(-2D0)
       ENDDO
       SUM = SUM + (0.117D0 + DEXP(C4/(T*9D0**2D0)))*(T/(2D0*C4))
       DSUM = DSUM + (1D0/(2D0*C4))*(0.117D0+DEXP(C4/(9D0**2D0*T)))- &
            (1D0/(2D0*9D0**2D0*T))*DEXP(C4/(T*9D0**2D0))
    ELSE
       SUM = (0.117D0 + DEXP(C4/(T*DBLE(N0)**2D0)))*(T/(2D0*C4))
       DSUM = (1D0/(2D0*C4))*(0.117D0+DEXP(C4/(9D0**2D0*T)))- &
            (1D0/(2D0*9D0**2D0*T))*DEXP(C4/(T*9D0**2D0))
    ENDIF
    HYD = C2*SUM*(1D0-DEXP(-C3/(T*LAMBDA0)))*DEXP(-C4/T)*LAMBDA0**3D0
    DHYD_RHO=HYD*(DSUM/SUM+C4/(T**2D0))-(C2*C3*LAMBDA0**2D0/T**2D0)*DEXP(-C4/T)*DEXP(-C3/(LAMBDA0*T))
    ! All dependences are with T, so the derivative at constant gas pressure and density are the same
    DHYD_PG=DHYD_RHO
    IF (PRESENT(OPAC).AND.PRESENT(DOPACDT_RHO).AND.PRESENT(DOPACDT_PG)) THEN
       OPAC =  HYD
       DOPACDT_RHO = DHYD_RHO
       DOPACDT_PG = DHYD_PG
    ENDIF
    !
  END SUBROUTINE OPAC_HYDROGEN
  !
  !------------------------------------------------
  !
  SUBROUTINE OPAC_CONTINUUM(T,LAMBDA0,NE,NH,KC,DKCDT_RHO,DKCDT_PG,DKCDPG_TEM)
    ! This routine calculates the continuum absorption coefficient and
    ! its derivatives with respect to temperature at constant density
    ! and constant gas pressure (separately).
    IMPLICIT NONE
    ! Input/Output
    REAL(SP), INTENT(IN)              :: T
    REAL(DP), INTENT(IN)              :: LAMBDA0, NE, NH
    REAL(DP), INTENT(OUT)             :: KC
    ! Internal
    REAL(DP)                          :: NI,NII,NIII,FACTOR
    ! Optional
    REAL(DP), INTENT(OUT), OPTIONAL   :: DKCDT_RHO, DKCDT_PG, DKCDPG_TEM
    !
    REAL(DP)                          :: T8
    !
    T8=DBLE(T)
    !
    ! Determine number of negative hydrogen ions per cm^3
    CALL ION_CALC_TEMPRHO(1,T8,NE,0D0,NH,NI,NII,NIII)
    !
    CALL OPAC_HMINUS_FF(T8,NE,LAMBDA0)
    CALL OPAC_HMINUS_BF(T8,NE,LAMBDA0)
    CALL OPAC_THOMSON(NE,NH)
    !
    CALL OPAC_RAYLEIGH_H(LAMBDA0)
    CALL OPAC_HYDROGEN(T8,LAMBDA0)
    !
    FACTOR=HMINUS_FF+HMINUS_BF+THOMSON+RAYLEIGH_H+HYD
    ! All opacities are per hydrogen atom, so the total is:
    KC=NH*FACTOR
    !
    ! Now derivatives
    IF (MRESPFUNCT.EQV..TRUE.) THEN
      DKCDTEMP_RHO=FACTOR*DNHYDDTEMP_RHO+NH*(DHMINUS_FF_RHO&
          +DHMINUS_BF_RHO+DTHOMSON_RHO+DHYD_RHO)
      DKCDTEMP_PG=FACTOR*DNHYDDTEMP_PG+NH*(DHMINUS_FF_PG&
          +DHMINUS_BF_PG+DTHOMSON_PG+DHYD_PG)
      DKCDPG_TEMP=FACTOR*DNHYDDPG_TEMP+NH*(DHMINUS_FFDPG_TEMP&
          +DHMINUS_BFDPG_TEMP+DTHOMSONDPG_TEMP)
      DKCDRHO_TEMP=FACTOR*DNHYDDRHO_TEMP+NH*(DHMINUS_FFDRHO_TEMP&
          +DHMINUS_BFDRHO_TEMP+DTHOMSONDRHO_TEMP)
      !
      ! Optional output
      IF (PRESENT(DKCDT_RHO).AND.PRESENT(DKCDT_PG).AND.PRESENT(DKCDPG_TEM)) THEN
         DKCDT_RHO  = DKCDTEMP_RHO
         DKCDT_PG   = DKCDTEMP_PG
         DKCDPG_TEM = DKCDPG_TEMP
      ENDIF
      !
    ENDIF
    !
  END SUBROUTINE OPAC_CONTINUUM
  !
  !================================================
  !
END MODULE CONT_OPACITY
!
