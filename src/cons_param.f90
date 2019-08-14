!
MODULE CONS_PARAM
  !
  ! J M Borrero
  ! Jan 7, 2007
  ! HAO-NCAR for HMI-Stanford
  !
  !================================================
  !
  IMPLICIT NONE
  !
  INTEGER, PUBLIC, PARAMETER  :: DP = KIND(1.D0)
  INTEGER, PUBLIC, PARAMETER  :: SP = KIND(1.)
  REAL(DP), PUBLIC, PARAMETER :: DPI = 3.141592653589793238462643D0
  REAL(SP), PUBLIC, PARAMETER :: SPI = 3.141592653589793238462643E0
  REAL(DP), PUBLIC, PARAMETER :: DD2R = DPI/180D0                                   ! deg -> rad (double precision)
  REAL(SP), PUBLIC, PARAMETER :: SD2R = SPI/180.0E0                                 ! deg -> rad (single precision)
  REAL(DP), PUBLIC, PARAMETER :: LIGHT = 2.99792458D+10                             ! Speed of light (cm/s)
  REAL(DP), PUBLIC, PARAMETER :: KBOL = 1.38056503D-16                              ! Boltzmann constant in (erg/K = cm^2 g/s^2 K)
  REAL(DP), PUBLIC, PARAMETER :: MHYD = 1.67262158D-24                              ! Proton mass (g)
  REAL(DP), PUBLIC, PARAMETER :: MELE = 9.10938215D-28                              ! Electron mass (g)
  REAL(DP), PUBLIC, PARAMETER :: QELE = 4.80320425D-10                              ! Electron charge (erg*cm)^(1/2) = statCoulomb
  REAL(DP), PUBLIC, PARAMETER :: EVOLT = 1.60217656D-12                             ! Electronvolt (ergs = cm^2 g/s^2)
  REAL(DP), PUBLIC, PARAMETER :: MAMU = 1.66053892D-24                              ! Atomic mass unit (g): 1/12 weight of carbon atom
  REAL(DP), PUBLIC, PARAMETER :: HPLA = 6.62606957D-27                              ! Planck's constant (erg s = cm^2 g/s)
  REAL(DP), PUBLIC, PARAMETER :: CSAHA1 = HPLA**3.0D0/(2.0D0*DPI*MELE*KBOL)**(3.0D0/2.0D0)  ! Constant # 1 in Saha Equation
  REAL(DP), PUBLIC, PARAMETER :: CSAHA2 = EVOLT/KBOL                                ! Constant # 2 in Saha Equation
  REAL(DP), PUBLIC, PARAMETER :: RBOHR = HPLA**2.0D0/(4.0D0*DPI**2.0D0*MELE*QELE**2.0D0)    ! Bohr's radius (cm)
  REAL(DP), PUBLIC, PARAMETER :: GRAV = 2.7414D4                                    ! Solar surface gravity (cm/s^2)
  REAL(DP), PUBLIC :: IMAT(4,4)                                          ! Identity matrix
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  !------------------------------------------------
  !
  !
  !================================================
  !
END MODULE CONS_PARAM
!
