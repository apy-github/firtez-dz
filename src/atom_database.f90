!
MODULE ATOM_DATABASE
  !
  !================================================
  !
  ! J M Borrero
  ! KIS. March 25, 2013
  !
  ! Basically cut/paste from atmdatb.f in SIR
  ! which is basically the same as provided by Axel Wittmann
  !
  USE CONS_PARAM, ONLY: SP, DP, MAMU
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER,  PARAMETER                   :: NELEM = 92
  REAL(DP), DIMENSION(NELEM)            :: XI, XII, MATOM, ABUND
  REAL(DP)                              :: DENO_SA, DENO_SAM
  INTEGER,  DIMENSION(NELEM)            :: ZATOM
  !
  PUBLIC :: XII
  PUBLIC :: XI
  PUBLIC :: MATOM
  PUBLIC :: ABUND
  PUBLIC :: NELEM
  PUBLIC :: DENO_SA
  PUBLIC :: DENO_SAM
  !
  PUBLIC :: ATOM_INIT
  PUBLIC :: PARTITION_FUNCTION
  PUBLIC :: REFRAX
  PRIVATE :: ZATOM
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! atom_init
  ! subroutine
  ! refrax
  !
  !------------------------------------------------
  !
  SUBROUTINE ATOM_INIT
    !
    INTEGER                             :: I
    !
    ! Atomic number
    DO I=1,NELEM
       ZATOM(I)=I
    ENDDO
    ! First ionization potential
    ! In the case of hydrogen refers to binding energy of H- ion
    XI(:)=(/0.754D0,24.58D0,5.39D0,9.32D0,8.298D0,11.256D0,14.529D0,13.614D0 &
        ,17.418D0,21.559D0,5.138D0,7.644D0,5.984D0,8.149D0,10.474D0,10.357D0 &
        ,13.012D0,15.755D0,4.339D0,6.111D0,6.538D0,6.825D0,6.738D0,6.763D0 &
        ,7.432D0,7.896D0,7.863D0,7.633D0,7.724D0,9.391D0,5.997D0,7.88D0 &
        ,9.81D0,9.75D0,11.840D0,13.996D0,4.176D0,5.692D0,6.377D0,6.838D0 &
        ,6.881D0,7.10D0,7.28D0,7.36D0,7.46D0,8.33D0,7.574D0,8.991D0,5.785D0 &
        ,7.34D0,8.64D0,9.01D0,10.454D0,12.127D0,3.893D0,5.210D0,5.577D0 &
        ,5.466D0,5.422D0,5.489D0,5.554D0,5.631D0,5.666D0,6.141D0,5.852D0 &
        ,5.927D0,6.018D0,6.101D0,6.184D0,6.254D0,5.426D0,6.650D0,7.879D0 &
        ,7.980D0,7.870D0,8.70D0,9.10D0,9.00D0,9.22D0,10.43D0,6.105D0 &
        ,7.415D0,7.287D0,8.43D0,9.30D0,10.745D0,4.D0,5.276D0,6.9D0,6.0D0 &
        ,6.0D0,6.0D0/)
    ! Second ionization potential
    ! In the case of hydrogen refers in fact ot first ionization potential
    XII(:)=(/13.595,54.403,75.62,18.21,25.15,24.376,29.59,35.11,34.98,41.07, &
         7.290,15.03,18.823,16.34,19.72,23.405,23.798,27.62,31.81,11.868, &
         12.891,13.63,14.205,16.493,15.636,16.178,17.052,18.15,20.286,17.96, &
         20.509,15.93,18.63,21.50,21.60,24.565,27.50,11.027,12.233,13.13, &
         14.316,16.15,15.26,16.76,18.07,19.42,21.48,16.904,18.86,14.63, &
         16.50,18.60,19.09,21.20,25.10,10.001,11.060,10.850,10.550,10.730, &
         10.899,11.069,11.241,12.090,11.519,11.670,11.800,11.930,12.050, &
         12.184,13.900,14.900,16.2,17.7,16.60,17.00,20.00,18.56,20.50,18.75, &
         20.42,15.03,16.68,19.,20.,20.,22.,10.144,12.1,12.,12.,12./)
!    ! Abundances THEVENIN 1989
    ABUND(:)=(/12.00,11.00,1.00,1.15,2.60,8.69,7.99,8.91,4.56,8.00,6.28, &
         7.53,6.43,7.50,5.45,7.21,5.50,6.58,5.05,6.36,2.99,4.88,3.91,5.61, &
         5.47,7.46,4.85,6.18,4.24,4.60,2.88,3.57,2.39,3.35,2.63,3.21,2.60, &
         2.93,2.18,2.46,1.46,2.10,0.00,1.78,1.10,1.69,0.94,1.86,1.66,2.00, &
         1.00,2.25,1.51,2.19,1.12,2.18,1.07,1.58,0.76,1.40,0.00,0.88,0.48, &
         1.13,0.20,1.07,0.26,0.93,0.00,1.08,0.76,0.88,-0.09,0.98,0.26,1.45, &
         1.36,1.80,1.13,1.27,0.90,1.90,0.71,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0, &
         0.02,-8.0,-0.47/)
     ! Abundances Vigeesh
!    ABUND(:)=(/12.00,10.93,1.10,1.40,2.79,8.41,7.80,8.66,4.48,8.08,6.33,7.58&
!        ,6.47,7.55,5.45,7.33,5.28,6.40,5.12,6.36,3.17,5.02,4.00,5.67,5.39,7.50&
!        ,4.92,6.25,4.21,4.60,2.88,3.41,2.37,3.41,2.63,3.31,2.60,2.97,2.24,2.60&
!        ,1.42,1.92,-8.00,1.84,1.12,1.69,1.24,1.77,0.82,2.00,1.00,2.24,1.51&
!        ,2.17,1.13,2.13,1.17,1.58,0.71,1.50,-8.00,1.01,0.51,1.12,0.35,1.14&
!        ,0.51,0.93,0.15,1.08,0.06,0.88,-0.13,0.69,0.28,1.45,1.35,1.80,0.88&
!        ,1.13,0.83,1.95,0.71,-8.00,-8.00,-8.00,-8.00,-8.00,-8.00,0.09&
!        ,-8.00,-0.50/)
    ! Atomic mass
    MATOM(:)=(/1.008,4.003,6.939,9.012,10.811,12.011,14.007,16.,18.998, &
         20.183,22.99,24.312,26.982,28.086,30.974,32.064,35.453,39.948, &
         39.102,40.08,44.956,47.90,50.942,51.996,54.938,55.847,58.933, &
         58.71,63.54,65.37,69.72,72.59,74.92,78.96,79.91,83.80,85.47,87.62, &
         88.905,91.22,92.906,95.94, 99.00,101.07,102.9,106.4,107.87,112.40, &
         114.82,118.69,121.75,127.6,126.9,131.3,132.9,137.34,138.91,140.12, &
         140.91,144.24,147.00,150.35,151.96,157.25,158.92,162.50,164.93, &
         167.26,168.93,173.04,174.97,178.49,180.95,183.85,186.2,190.2,192.2, &
         195.09,196.97,200.59,204.37,207.19,208.98,210.,211.,222.,223.,226.1, &
         227.1,232.04,231.,238.03/)
    ! Calculate SA & SAM
    DENO_SA=0D0
    DENO_SAM=0D0
    DO I=1,NELEM
       DENO_SA=DENO_SA+10D0**(ABUND(I)-12D0)
       DENO_SAM=DENO_SAM+MAMU*MATOM(I)*10D0**(ABUND(I)-12D0)
    ENDDO
    !
  END SUBROUTINE ATOM_INIT
  !
  !------------------------------------------------
  !
  PURE SUBROUTINE PARTITION_FUNCTION(INDEX,T,UI,UII,UIII,DUI,DUII,DUIII)
    !
    INTEGER,    INTENT(IN)                :: INDEX
    REAL(DP),   INTENT(IN)                :: T
    REAL(DP),   INTENT(OUT)               :: UI, UII, UIII
    REAL(DP),   INTENT(OUT)               :: DUI, DUII, DUIII
    REAL(DP)                              :: X, Y, DX
    !
    UI=0.
    UII=0.
    UIII=0.
    DUI=0.
    DUII=0.
    DUIII=0.
    !
    X=DLOG(5040./T)
    Y=1E-3*T
    DX=-1./T
    !
    SELECT CASE (INDEX)
       CASE(1)
          ! Hydrogen
          UI=1.    ! H-
          UII=2.   ! H
          IF (T .GT. 1.3E4) THEN
             UII=1.51+3.8E-5*T
             DUII=3.8E-5/UII
          ENDIF
          IF (T .GT. 1.62E4) THEN
             UII=11.41+T*(-1.1428E-3+T*3.52E-8)
             DUII=(-1.1428E-3+2.*T*3.52E-8)/UII
          ENDIF
          UIII=1.  ! H+
       CASE(2)
          ! Helium
          UI=1.
          IF (T .GT. 3E4) THEN
             UI=14.8+T*(-9.4103E-4+T*1.6095E-8)
             DUI=(-9.4103E-4+2*T*1.6095E-8)/UI
          ENDIF
          UII=2.
          UIII=1.
       CASE(3)
          ! Lithium
          UI=2.081-Y*(6.8926E-2-Y*1.4081E-2)
          DUI=1.E-3*(-6.8926E-2+2*Y*1.4081E-2)/UI
          IF (T .GT. 6E3) THEN
             UI=3.4864+T*(-7.3292E-4+T*8.5586E-8)
             DUI=(-7.3292E-4+2*T*8.5586E-8)/UI
          ENDIF
          UII=1.
          UIII=2.
       CASE(4)
          ! Berylium
          UI=MAX(1.,.631+7.032E-5*T)
          IF (UI .NE. 1.) DUI=7.032E-5/UI
          UII=2.
          UIII=1.
       CASE(5)
          ! Boron
          UI=5.9351+1.0438E-2*Y
          DUI=1.E-3*1.0438e-2/UI
          UII=1.
          UIII=2.
       CASE(6)
          ! Carbon
          UI=8.6985+Y*(2.0485E-2+Y*(1.7629E-2-3.9091E-4*Y))
          DUI=1.E-3*((2.0485E-2+2*Y*1.7629E-2-3*Y*Y*3.9091E-4)/UI)
          IF (T .GT. 1.2E4) THEN
             UI=13.97+T*(-1.3907E-3+T*9.0844E-8)
             DUI=(-1.3907E-3+2*T*9.0844E-8)/UI
          END IF
          UII=5.838+1.6833E-5*T
          DUII=1.6833E-5/UII
          IF (T .GT. 2.4E4) THEN
             UII=10.989+T*(-6.9347E-4+T*2.0861E-8)
             DUII=(-6.9347E-4+2*T*2.0861E-8)/UII
          ENDIF
          UIII=1.
          IF (T .GT. 1.95E4) THEN
             UIII=-0.555+8E-5*T
             DUIII=8E-5/UIII
          ENDIF
       CASE(7)
          ! Nitrogen
          UI=3.9914+Y*(1.7491E-2-Y*(1.0148E-2-Y*1.7138E-3))
          DUI=1.E-3*((1.7491E-2-2*Y*1.0148E-2+3*Y*Y*1.7138E-3)/UI)
          IF (T .GT. 8800.) THEN
             UI=2.171+2.54E-4*T
             DUI=2.54E-4/UI
          ENDIF
          IF (T .GT. 1.8E4) THEN
             UI=11.396+T*(-1.7139E-3+T*8.633E-8)
             DUI=(-1.7139E-3+2*T*8.633e-8)/UI
          ENDIF
          UII=8.060+1.420E-4*T
          DUII=1.420E-4/UII
          IF (T .GT. 3.3E4) THEN
             UII=26.793+T*(-1.8931E-3+T*4.4612E-8)
             DUII=(-1.8931E-3+2*T*4.4612E-8)/UII
          ENDIF
          UIII=5.9835+T*(-2.6651E-5+T*1.8228E-9)
          DUIII=(-2.6651E-5+2*T*1.8228E-9)/UIII
          IF (T.LT.7310.5) THEN
             UIII=5.89
             DUIII=0.
          ENDIF
       CASE(8)
          ! Oxygen
          UI=8.29+1.10E-4*T
          DUI=(1.10E-4)/UI
          IF (T .GT. 1.9E4) THEN
             UI=66.81+T*(-6.019E-3+T*1.657E-7)
             DUI=(-6.019E-3+2*T*1.657E-7)/UI
          ENDIF
          UII=MAX(4.,3.51+8.E-5*T)
          IF (UII .ne. 4.) DUII=8.E-5/UII
          IF (T .GT. 3.64E4) THEN
             UII=68.7+T*(-4.216E-3+T*6.885E-8)
             DUII=(-4.216E-3+2*T*6.885E-8)/UII
          ENDIF
          UIII=7.865+1.1348E-4*T
          DUIII=1.1348E-4/UIII
       CASE(9)
          ! Fluorine
          UI=4.5832+Y*(.77683+Y*(-.20884+Y*(2.6771E-2-1.3035E-3*Y)))
          DUI=.77683-2*Y*.20884+3*Y*Y*2.6771E-2-4*Y*Y*Y*1.3035E-3
          DUI=1.E-3*(DUI/UI)
          IF (T .GT. 8750.) THEN
             UI=5.9
             DUI=0.
          ENDIF
          IF (T .GT. 2E4) THEN
             UI=15.16+T*(-9.229E-4+T*2.312E-8)
             DUI=(-9.229E-4+2*T*2.312E-8)/UI
          ENDIF
          UII=8.15+8.9E-5*T
          DUII=8.9E-5/UII
          UIII=2.315+1.38E-4*T
          DUIII=1.38E-4/UIII
       CASE(10)
          ! Neon
          UI=1.
          IF (T .GT. 2.69E4) THEN
             UI=26.3+T*(-2.113E-3+T*4.359E-8)
             DUI=(-2.113E-3+2*T*4.359E-8)/UI
          ENDIF
          UII=5.4+4.E-5*T
          DUII=4.E-5/UII
          UIII=7.973+7.956E-5*T
          DUIII=7.956E-5/UIII
       CASE(11)
          ! Sodium
          UI=MAX(2.,1.72+9.3E-5*T)
          IF (UI .NE. 2.) DUI=9.3E-5/UI
          IF (T .GT. 5400.) THEN
             UI=-0.83+5.66E-4*T
             DUI=5.66E-4/UI
          ENDIF
          IF (T .GT. 8.5E3) THEN
             UI=4.5568+T*(-1.2415E-3+T*1.3861E-7)
             DUI=(-1.2415E-3+2*T*1.3861E-7)/UI
          ENDIF
          UII=1.000
          UIII=5.69+5.69E-6*T
          DUIII=5.69E-6/UIII
       CASE(12)
          ! Magnesium
          UI=1.+EXP(-4.027262-X*(6.173172+X*(2.889176+X*(2.393895+.784131*X))))
          DUI=-DX*(6.173172+X*(2*2.889176+X*(3*2.393895+4*X*.784131)))
          DUI=DUI*(UI-1.)/UI
          IF (T .GT. 8E3) THEN
             UI=2.757+T*(-7.8909E-4+T*7.4531E-8)
             DUI=(-7.8909E-4+2*T*7.4531E-8)/UI
          ENDIF
          UII=2.+EXP(-7.721172-X*(7.600678+X*(1.966097+.212417*X)))
          DUII=-DX*(7.0600678+X*(2*1.966097+3*X*.212417))
          DUII=DUII*(UII-2.)/UII
          IF (T .GT. 2E4) THEN
             UII=7.1041+T*(-1.0817E-3+T*4.7841E-8)
             DUII=(-1.0817E-3+2*T*4.7841E-8)/UII
          ENDIF
          UIII=1.0000
       CASE(13)
          ! Aluminium
          UI=5.2955+Y*(.27833-Y*(4.7529E-2-Y*3.0199E-3))
          DUI=1.E-3*(.27833-Y*(2*4.7529E-2-3*Y*3.0199E-3))/UI
          UII=MAX(1.,.725+3.245E-5*T)
          IF (UII .NE. 1.) DUII=3.245E-5/UII
          IF (T .GT. 2.24E4) THEN
             UII=61.06+T*(-5.987E-3+T*1.485E-7)
             DUII=(-5.987E-3+2*T*1.485E-7)/UII
          ENDIF
          UIII=MAX(2.,1.976+3.43E-6*T)
          IF (UIII .NE. 2.) DUIII=3.43E-6/UIII
          IF (T .GT. 1.814E4) THEN
             UIII=3.522+T*(-1.59E-4+T*4.382E-9)
             DUIII=(-1.59E-4+2*T*4.382E-9)/UIII
          ENDIF
       CASE(14)
          ! Silicon
          UI=6.7868+Y*(.86319+Y*(-.11622+Y*(.013109-6.2013E-4*Y)))
          DUI=1.E-3*(.86319+Y*(-2*.11622+Y*(3*.013109-4*Y*6.2013E-4)))/UI
          IF (T .GT. 1.04E4) THEN
             UI=86.01+T*(-1.465E-2+T*7.282E-7)
             DUI=(-1.465E-2+2*T*7.282E-7)/UI
          ENDIF
          UII=5.470+4.E-5*T
          DUII=4.E-5/UII
          IF (T .GT. 1.8E4) THEN
             UII=26.44+T*(-2.22E-3+T*6.188E-8)
             DUII=(-2.22E-3+2*T*6.188E-8)/UII
          ENDIF
          UIII=MAX(1.,.911+1.1E-5*T)
          IF (UIII .NE. 1.) DUIII=1.1E-5/UIII
          IF (T .GT. 3.33E4) THEN
             UIII=19.14+T*(-1.408E-3+T*2.617E-8)
             DUIII=(-1.408E-3+2*T*2.617E-8)/UIII
          ENDIF
       CASE(15)
          ! Phosphorus
          UI=4.2251+Y*(-.22476+Y*(.057306-Y*1.0381E-3))
          DUI=1.E-3*(-.22476+Y*(2*.057306-3*Y*1.0381E-3))/UI
          IF (T .GT. 6.E3) THEN
             UI=1.56+5.2E-4*T
             DUI=5.2E-4/UI
          ENDIF
          UII=4.4151+Y*(2.2494+Y*(-.55371+Y*(.071913-Y*3.5156E-3)))
          DUII=1.E-3*(2.2494+Y*(-2*.55371+Y*(3*.071913-4*Y*3.5156E-3)))/UII
          IF (T.GT.7250.) THEN
             UII=4.62+5.38E-4*T
             DUII=5.38E-4/UII
          ENDIF
          UIII=5.595+3.4E-5*T
          DUIII=3.4E-5/UIII
       CASE(16)
          ! Sulfur
          UI=7.5+2.15E-4*T
          DUI=2.15E-4/UI
          IF (T .GT. 1.16E4) THEN
             UI=38.76+T*(-4.906E-3+T*2.125E-7)
             DUI=(-4.906E-3+2*T*2.125E-7)/UI
          ENDIF
          UII=2.845+2.43E-4*T
          DUII=2.43E-4/UII
          IF (T .GT. 1.05E4) THEN
             UII=6.406+T*(-1.68E-4+T*1.323E-8)
             DUII=(-1.68E-4+T*2*1.323E-8)/UII
          ENDIF
          UIII=7.38+1.88E-4*T
          DUIII=1.88E-4/UIII
       CASE(17)
          ! Chlorine
          UI=5.2+6.E-5*T
          DUI=6.E-5/UI
          IF (T .GT. 1.84E4) THEN
             UI=-81.6+4.8E-3*T
             DUI=4.8E-3/UI
          ENDIF
          UII=7.0+2.43E-4*T
          DUII=2.43E-4/UII
          UIII=2.2+2.62E-4*T
          DUIII=2.62E-4/UIII
       CASE(18)
          ! Argon
          UI=1.000
          UII=5.20+3.8E-5*T
          DUII=3.8E-5/UII
          UIII=7.474+1.554E-4*T
          DUIII=1.554E-4/UIII
       CASE(19)
          ! Potassium
          UI=1.9909+Y*(.023169-Y*(.017432-Y*4.0938E-3))
          DUI=1.E-3*(.023169-Y*(2*.017432-Y*3*4.0938E-3))/UI
          IF (T .GT. 5800.) THEN
             UI=-9.93+2.124E-3*T
             DUI=2.124E-3/UI
          ENDIF
          UII=1.000
          UIII=5.304+1.93E-5*T
          DUIII=1.93E-5/UIII
       CASE(20)
          ! Calcium
          UI=1.+EXP(-1.731273-X*(5.004556+X*(1.645456+X*(1.326861+.508553*X))))
          DUI=-DX*(5.004556+X*(2*1.645456+X*(3*1.326861+4*.508553*X)))
          DUI=DUI*(UI-1.)/UI
          UII=2.+EXP(-1.582112-X*(3.996089+X*(1.890737+.539672*X)))
          DUII=-DX*(3.996089+X*(2*1.890737+X*3*.539672))
          DUII=DUII*(UII-2.)/UII
          UIII=1.000
       CASE(21)
          ! Scandium
          UI=4.+EXP(2.071563+X*(-1.2392+X*(1.173504+.517796*X)))
          DUI=DX*(-1.2392+X*(2*1.173504+3*X*.517796))
          DUI=DUI*(UI-4.)/UI
          UII=3.+EXP(2.988362+X*(-.596238+.054658*X))
          DUII=DX*(-.596238+2*X*.054658)
          DUII=DUII*(UII-3.)/UII
          UIII=10.
       CASE(22)
          ! Titanium
          UI=5.+EXP(3.200453+X*(-1.227798+X*(.799613+.278963*X)))
          DUI=DX*(-1.227798+X*(2*.799613+3*X*.278963))
          DUI=DUI*(UI-5.)/UI
          IF (T .LT. 5.5E3) THEN
             UI=16.37+T*(-2.838E-4+T*5.819E-7)
             DUI=(-2.838E-4+2*T*5.819E-7)/UI
          ENDIF
          UII=4.+EXP(3.94529+X*(-.551431+.115693*X))
          DUII=DX*(-.551431+2.*X*.115693)
          DUII=DUII*(UII-4.)/UII
          UIII=16.4+8.5E-4*T
          DUIII=8.5E-4/UIII
       CASE(23)
          ! Vanadium
          UI=4.+EXP(3.769611+X*(-.906352+X*(.724694+.1622*X)))
          DUI=DX*(-.906352+X*(2*.724694+3.*X*.1622))
          DUI=DUI*(UI-4.)/UI
          UII=1.+EXP(3.755917+X*(-.757371+.21043*X))
          DUII=DX*(-.757371+X*2.*.21043)
          DUII=DUII*(UII-1.)/UII
          UIII=-18.+1.03E-2*T
          DUIII=1.03E-2/UIII
          IF (T.LT.2.25E3) THEN
             UIII=2.4E-3*T
             DUIII=2.4E-3/UIII
          ENDIF
       CASE(24)
          ! Chromium
          UI=7.+EXP(1.225042+X*(-2.923459+X*(.154709+.09527*X)))
          DUI=DX*(-2.923459+X*(2*.154709+X*3*.09527))
          DUI=DUI*(UI-7.)/UI
          UII=6.+EXP(.128752-X*(4.143973+X*(1.096548+.230073*X)))
          DUII=-DX*(4.143973+X*(2*1.096548+3.*X*.230073))
          DUII=DUII*(UII-6.)/UII
          UIII=10.4+2.1E-3*T
          DUIII=2.1E-3/UIII
       CASE(25)
          ! Manganese
          UI=6.+EXP(-.86963-X*(5.531252+X*(2.13632+X*(1.061055+.265557*X))))
          DUI=-DX*(5.531252+X*(2*2.13632+X*(3*1.061055+4*.265557*X)))
          DUI=DUI*(UI-6.)/UI
          UII=7.+EXP(-.282961-X*(3.77279+X*(.814675+.159822*X)))
          DUII=-DX*(3.77279+X*(2*.814675+3*X*.159822))
          DUII=DUII*(UII-7.)/UII
          UIII=10.
       CASE(26)
          ! Iron
          UI=9.+EXP(2.930047+X*(-.979745+X*(.76027+.118218*X)))
          DUI=DX*(-.979745+X*(2*.76027+3*X*.118218))
          DUI=DUI*(UI-9.)/UI
          IF (T .LT. 4E3) THEN
             UI=15.85+T*(1.306E-3+T*2.04E-7)
             DUI=(1.306E-3+2*T*2.04E-7)/UI
          ENDIF
          IF (T .GT. 9E3) THEN
             UI=39.149+T*(-9.5922E-3+T*1.2477E-6)
             DUI=(-9.5922E-3+2*T*1.2477E-6)/UI
          ENDIF
          UII=10.+EXP(3.501597+X*(-.612094+.280982*X))
          DUII=DX*(-.612094+2*X*.280982)
          DUII=DUII*(UII-10.)/UII
          IF (T .GT. 1.8E4) THEN
             UII=68.356+T*(-6.1104E-3+T*5.1567E-7)
             DUII=(-6.1104E-3+2*T*5.1567E-7)/UII
          ENDIF
          UIII=17.336+T*(5.5048E-4+T*5.7514E-8)
          DUIII=(5.5048E-4+2*T*5.7514E-8)/UIII
       CASE(27)
          ! Cobalt
          UI=8.65+4.9E-3*T
          DUI=4.9e-3/UI
          UII=11.2+3.58E-3*T
          DUII=3.58e-3/UII
          UIII=15.0+1.42E-3*T
          DUIII=1.42E-3/UIII
       CASE(28)
          ! Nickel
          UI=9.+EXP(3.084552+X*(-.401323+X*(.077498-.278468*X)))
          DUI=DX*(-.401323+X*(2*.077498-3.*X*.278468))
          DUI=DUI*(UI-9.)/UI
          UII=6.+EXP(1.593047-X*(1.528966+.115654*X))
          DUII=-DX*(1.528966+2.*X*.115654)
          DUII=DUII*(UII-6.)/UII
          UIII=13.3+6.9E-4*T
          DUIII=6.9E-4/UIII
       CASE(29)
          ! Copper
          UI=MAX(2.,1.50+1.51E-4*T)
          IF (UI.NE.2.)DUI=1.51E-4/UI
          IF (T.GT.6250.) THEN
             UI=-.3+4.58E-4*T
             DUI=4.58E-4/UI
          ENDIF
          UII=MAX(1.,.22+1.49E-4*T)
          IF (UII.NE.1.)DUII=1.49E-4/UII
          UIII=8.025+9.4E-5*T
          DUIII=9.4E-5/UIII
       CASE(30)
          ! Zinc
          UI=MAX(1.,.632+5.11E-5*T)
          IF (UI.NE.1.)DUI=5.11E-5/UI
          UII=2.00
          UIII=1.00
       CASE(31)
          ! Gallium
          UI=1.7931+Y*(1.9338+Y*(-.4643+Y*(.054876-Y*2.5054E-3)))
          DUI=1.E-3*(1.9338+Y*(-2.*.4643+Y*(3*.054876-Y*4*2.5054E-3)))/UI
          IF (T.GT.6.E3) THEN
             UI=4.18+2.03E-4*T
             DUI=2.03E-4/UI
          ENDIF
          UII=1.0
          UIII=2.0
       CASE(32)
          ! Germanium
          UI=6.12+4.08E-4*T
          DUI=4.08E-4/UI
          UII=3.445+1.78E-4*T
          DUII=1.78E-4/UII
          UIII=1.1
       CASE(33)
          ! Arsenic
          UI=2.65+3.65E-4*T
          DUI=3.65E-4/UI
          UII=-.25384+Y*(2.284+Y*(-.33383+Y*(.030408-Y*1.1609E-3)))
          DUII=1.E-3*(2.284+Y*(-2*.33383+Y*(3*.030408-4*Y*1.1609E-3)))/UII
          IF (T.GT.1.2E4) THEN
             UII=8.
             DUII=0.
          ENDIF
          UIII=8.
       CASE(34)
          ! Selenium
          UI=6.34+1.71E-4*T
          DUI=1.71E-4/UI
          UII=4.1786+Y*(-.15392+3.2053E-2*Y)
          DUII=1.E-3*(-.15392+3.2053E-2*Y*2)/UII
          UIII=8.
       CASE(35)
          ! Bromine
          UI=4.12+1.12E-4*T
          DUI=1.12E-4/UI
          UII=5.22+3.08E-4*T
          DUII=3.08E-4/UII
          UIII=2.3+2.86E-4*T
          DUIII=2.86E-4/UIII
       CASE(36)
          ! Krypton
          UI=1.00
          UII=4.11+7.4E-5*T
          DUII=7.4E-5/UII
          UIII=5.35+2.23E-4*T
          DUIII=2.23E-4/UIII
       CASE(37)
          ! Rubidium
          UI=MAX(2.,1.38+1.94E-4*T)
          IF (UI.NE.2.) DUI=1.94E-4/UI
          IF (T.GT.6250.) THEN
             UI=-14.9+2.79E-3*T
             DUI=2.79E-3/UI
          ENDIF
          UII=1.000
          UIII=4.207+4.85E-5*T
          DUIII=4.85E-5/UIII
       CASE(38)
          ! Strontium
          UI=.87127+Y*(.20148+Y*(-.10746+Y*(.021424-Y*1.0231E-3)))
          DUI=1.E-3*(.20148+Y*(-.10746*2.+Y*(.021424*3.-Y*1.0231E-3*4.)))
          DUI=DUI/UI
          IF (T.GT.6500.) THEN
             UI=-6.12+1.224E-3*T
             DUI=1.224E-3/UI
          ENDIF
          UII=MAX(2.,.84+2.6E-4*t)
          IF (UII.NE.2.) DUII=2.6E-4/UII
          UIII=1.0
          !
          ! From now on we ignore derivatives
          !
       CASE(39)
          ! Yttrium
          UI=.2+2.58E-3*T
          DUI=2.58E-3/UI
          UII=7.15+1.855E-3*T
          DUII=1.855E-3/UII
          UIII=9.71+9.9E-5*T
          DUIII=9.9E-5/UIII
       CASE(40)
          ! Zirconium
          UI=76.31+T*(-1.866E-2+T*2.199E-6)
          DUI=(-1.866E-2+T*2.199E-6*2.)/UI
          IF (T.LT.6236.) UI=6.8+T*(2.806E-3+T*5.386E-7)
          IF (T.LT.6236.) DUI=(2.806E-3+T*5.386E-7*2.)/UI
          UII=4.+EXP(3.721329-.906502*X)
          DUII=-DX*.906502
          DUII=DUII*(UII-4.)/UII
          UIII=12.3+1.385E-3*T
          DUIII=1.385E-3/UIII
       CASE(41)
          ! Niobium
          UI=MAX(1.,-19.+1.43E-2*T)
          IF (UI.NE.1.) DUI=1.43E-2/UI
          UII=-4.+1.015E-2*T
          DUII=1.015E-2/UII
          UIII=25.
       CASE(42)
          ! Molybdenum
          UI=MAX(7.,2.1+1.5E-3*T)
          IF (UI.NE.7.) DUI=1.5E-3/UI
          IF (T.GT.7.E3) UI=-38.1+7.28E-3*T
          IF (T.GT.7.E3) DUI=7.28E-3/UI
          UII=1.25+1.17E-3*T
          DUII=1.17E-3/UII
          IF (T.GT.6900.) UII=-28.5+5.48E-3*T
          IF (T.GT.6900.) DUII=5.48E-3/UII
          UIII=24.04+1.464E-4*T
          DUIII=1.464E-4/UIII
       CASE(43)
          ! Technetium
          UI=4.439+Y*(.30648+Y*(1.6525+Y*(-.4078+Y*(.048401-Y*2.1538E-3))))
          IF (T.GT.6.E3) UI=24.
          IF (T.LE.6.E3) DUI=1.E-3*(.30648+Y*(1.6525*2.+Y* &
              (-.4078*3.+Y*(.048401*4.-Y*2.1538E-3*5))))/UI
          UII=8.1096+Y*(-2.963+Y*(2.369+Y*(-.502+Y*(.049656-Y*1.9087E-3))))
          IF (T.GT.6.E3) UII=17.
          IF (T.LE.6.E3) DUII=1.E-3*(-2.963+Y*(2.369*2.+Y* &
              (-.502*3.+Y*(.049656*4.-Y*1.9087E-3*5.))))/UII
          UIII=220.
       CASE(44)
          !Ruthenium
          UI=-3.+7.17E-3*T
          DUI=7.17E-3/UI
          UII=3.+4.26E-3*T
          DUII=4.26E-3/UII
          UIII=22.
       CASE(45)
          ! Rhodium
          UI=6.9164+Y*(3.8468+Y*(.043125-Y*(8.7907E-3-Y*5.9589E-4)))
          DUI=1.E-3*(3.8468+Y*(.043125*2.-Y*(8.7907E-3*3.-Y*5.9589E-4*4.)))/UI
          UII=7.2902+Y*(1.7476+Y*(-.038257+Y*(2.014E-3+Y*2.1218E-4)))
          DUII=1.E-3*(1.7476+Y*(-.038257*2.+Y*(2.014E-3*3.+Y*2.1218E-4*4.)))/UII
          UIII=30.
       CASE(46)
          ! Palladium
          UI=MAX(1.,-1.75+9.86E-4*T)
          IF(UI.NE.1.)DUI=9.86E-4/UI
          UII=5.60+3.62E-4*T
          DUII=3.62E-4/UII
          UIII=20.
       CASE(47)
          ! Silver
          UI=MAX(2.,1.537+7.88E-5*T)
          IF(UI.NE.2.)DUI=7.88E-5/UI
          UII=MAX(1.,0.73+3.4E-5*T)
          IF(UII.NE.1.)DUII=3.4E-5/UII
          UIII=6.773+1.248E-4*T
          DUIII=1.248E-4/UIII
       CASE(48)
          ! Cadmium
          UI=MAX(1.,.43+7.6E-5*T)
          IF(UI.NE.1.)DUI=7.6E-5/UI
          UII=2.
          UIII=1.
       CASE(49)
          ! Indium
          UI=2.16+3.92E-4*T
          DUI=3.92E-4/UI
          UII=1.
          UIII=2.
       CASE(50)
          ! Tin
          UI=2.14+6.16E-4*T
          DUI=6.16E-4/UI
          UII=2.06+2.27E-4*T
          DUII=2.27E-4/UII
          UIII=1.05
       CASE(51)
          ! Antinomy
          UI=2.34+4.86E-4*T
          DUI=4.86E-4/UI
          UII=.69+5.36E-4*T
          DUII=5.36E-4/UII
          UIII=3.5
       CASE(52)
          ! Tellurium
          UI=3.948+4.56E-4*T
          DUI=4.56E-4/UI
          UII=4.2555+Y*(-.25894+Y*(.06939-Y*2.4271E-3))
          DUII=1.E-3*(-.25894+Y*(.06939*2.-Y*2.4271E-3*3.))/UII
          IF (T.GT.1.2E4) UII=7.
          UIII=5.
       CASE(53)
          ! Iodine
          UI=MAX(4.,3.8+9.5E-5*T)
          IF(UI.NE.4.)DUI=9.5E-5/UI
          UII=4.12+3.E-4*T
          DUII=3.E-4/UII
          UIII=7.
       CASE(54)
          ! Xenon
          UI=1.
          UII=3.75+6.876E-5*T
          DUII=6.876E-5/UII
          UIII=4.121+2.323E-4*T
          DUIII=2.323E-4/UIII
       CASE(55)
          ! Cesium
          UI=MAX(2.,1.56+1.67E-4*T)
          IF(UI.NE.2.)DUI=1.67E-4/UI
          IF (T.GT.4850.) UI=-2.680+1.04E-3*T
          IF (T.GT.4850.) DUI=1.04E-3/UI
          UII=1.
          UIII=3.769+4.971E-5*T
          DUIII=4.971E-5/UIII
       CASE(56)
          ! Barium
          UI=MAX(1.,-1.8+9.85E-4*T)
          IF (UI.NE.1.) DUI=9.85E-4/UI
          IF (T.GT.6850.) THEN
             UI=-16.2+3.08E-3*T
             DUI=3.08E-3/UI
          ENDIF
          UII=1.11+5.94E-4*T
          DUII=5.94E-4/UII
          UIII=1.
       CASE(57)
          ! Lanthanum
          UI=15.42+9.5E-4*T
          DUI=9.5E-4/UI
          IF (T.GT.5060.) UI=1.+3.8E-3*T
          IF (T.GT.5060.) DUI=3.8E-3/UI
          UII=13.2+3.56E-3*T
          DUII=3.56E-3/UII
          UIII=12.
       CASE(58)
          ! Cerium
          UI=9.+EXP(5.202903+X*(-1.98399+X*(.119673+.179675*X)))
          DUI=DX*(-1.98399+X*(.119673*2.+.179675*3.*X))
          DUI=DUI*(UI-9.)/UI
          UII=8.+EXP(5.634882-X*(1.459196+X*(.310515+.052221*X)))
          DUII=-DX*(1.459196+X*(.310515*2.+.052221*3.*X))
          DUII=DUII*(UII-8.)/UII
          UIII=9.+EXP(3.629123-X*(1.340945+X*(.372409+X*(.03186-.014676*X))))
          DUIII=-DX*(1.340945+X*(.372409*2.+X*(.03186*3.-.014676*4.*X)))
          DUIII=DUIII*(UIII-9.)/UIII
       CASE(59)
          ! Praseodymium
          UII=9.+EXP(4.32396-X*(1.191467+X*(.149498+.028999*X)))
          DUII=-DX*(1.191467+X*(.149498*2.+.028999*3.*X))
          DUII=DUII*(UII-9.)/UII
          UI=UII
          DUI=DUII
          UIII=10.+EXP(3.206855+X*(-1.614554+X*(.489574+.277916*X)))
          DUIII=DX*(-1.614554+X*(.489574*2.+.277916*3.*X))
          DUIII=DUIII*(UIII-10.)/UIII
       CASE(60)
          ! Neodymium
          UI=9.+EXP(4.456882+X*(-2.779176+X*(.082258+X*(.50666+.127326*X))))
          DUI=DX*(-2.779176+X*(.082258*2.+X*(.50666*3.+.127326*4.*X)))
          DUI=DUI*(UI-9.)/UI
          UII=8.+EXP(4.689643+X*(-2.039946+X*(.17193+X*(.26392+.038225*X))))
          DUII=DX*(-2.039946+X*(.17193*2.+X*(.26392*3.+.038225*4.*X)))
          DUII=DUII*(UII-8.)/UII
          UIII=UII
          DUIII=DUII
       CASE(61)
          ! Promethium
          UI=20.
          UII=25.
          UIII=100.
       CASE(62)
          ! Samarium
          UI=1.+EXP(3.549595+X*(-1.851549+X*(.9964+.566263*X)))
          DUI=DX*(-1.851549+X*(.9964*2.+.566263*3.*X))
          DUI=DUI*(UI-1.)/UI
          UII=2.+EXP(4.052404+X*(-1.418222+X*(.358695+.161944*X)))
          DUII=DX*(-1.418222+X*(.358695*2.+.161944*3.*X))
          DUII=DUII*(UII-2.)/UII
          UIII=1.+EXP(3.222807-X*(.699473+X*(-.056205+X*(.533833+.251011*X))))
          DUIII=-DX*(.699473+X*(-.056205*2.+X*(.533833*3.+.251011*4.*X)))
          DUIII=DUIII*(UIII-1.)/UIII
       CASE(63)
          ! Europium
          UI=8.+EXP(1.024374-X*(4.533653+X*(1.540805+X*(.827789+.286737*X))))
          DUI=-DX*(4.533653+X*(1.540805*2.+X*(.827789*3.+.286737*4.*X)))
          DUI=DUI*(UI-8.)/UI
          UII=9.+EXP(1.92776+X*(-1.50646+X*(.379584+.05684*X)))
          DUII=DX*(-1.50646+X*(.379584*2.+.05684*3.*X))
          DUII=DUII*(UII-9.)/UII
          UIII=8.
       CASE(64)
          ! Gadolinium
          UI=5.+EXP(4.009587+X*(-1.583513+X*(.800411+.388845*X)))
          DUI=DX*(-1.583513+X*(.800411*2.+.388845*3.*X))
          DUI=DUI*(UI-5.)/UI
          UII=6.+EXP(4.362107-X*(1.208124+X*(-.074813+X*(.076453+.055475*X))))
          DUII=-DX*(1.208124+X*(-.074813*2.+X*(.076453*3.+.055475*4.*X)))
          DUII=DUII*(UII-6.)/UII
          UIII=5.+EXP(3.412951-X*(.50271+X*(.042489-4.017E-3*X)))
          DUIII=-DX*(.50271+X*(.042489*2.-4.017E-3*3.*X))
          DUIII=DUIII*(UIII-5.)/UIII
       CASE(65)
          ! Terbium
          UI=16.+EXP(4.791661+X*(-1.249355+X*(.570094+.240203*X)))
          DUI=DX*(-1.249355+X*(.570094*2.+.240203*3.*X))
          DUI=DUI*(UI-16.)/UI
          UII=15.+EXP(4.472549-X*(.295965+X*(5.88E-3+.131631*X)))
          DUII=-DX*(.295965+X*(5.88E-3*2.+.131631*3.*X))
          DUII=DUII*(UII-15.)/UII
          UIII=UII
          DUIII=DUII
       CASE(66)
          ! Dysprosium
          UI=17.+EXP(3.029646-X*(3.121036+X*(.086671-.216214*X)))
          DUI=-DX*(3.121036+X*(.086671*2.-.216214*3.*X))
          DUI=DUI*(UI-17.)/UI
          UII=18.+EXP(3.465323-X*(1.27062+X*(-.382265+X*(.431447+.303575*X))))
          DUII=-DX*(1.27062+X*(-.382265*2.+X*(.431447*3.+.303575*4.*X)))
          DUII=DUII*(UII-18.)/UII
          UIII=UII
          DUIII=DUII
       CASE(67)
          ! Holmium
          UIII=16.+EXP(1.610084-X*(2.373926+X*(.133139-.071196*X)))
          DUIII=-DX*(2.373926+X*(.133139*2.-.071196*3.*X))
          DUIII=DUIII*(UIII-16.)/UIII
          UI=UIII
          DUI=DUIII
          UII=UIII
          DUII=DUIII
       CASE(68)
          ! Erbium
          UI=13.+EXP(2.895648-X*(2.968603+X*(.561515+X*(.215267+.095813*X))))
          DUI=-DX*(2.968603+X*(.561515*2.+X*(.215267*3.+.095813*4.*X)))
          DUI=DUI*(UI-13.)/UI
          UII=14.+EXP(3.202542-X*(.852209+X*(-.226622+X*(.343738+.186042*X))))
          DUII=-DX*(.852209+X*(-.226622*2.+X*(.343738*3.+.186042*4.*X)))
          DUII=DUII*(UII-14.)/UII
          UIII=UII
          DUIII=DUII
       CASE(69)
          ! Thulium
          UI=8.+EXP(1.021172-X*(4.94757+X*(1.081603+.034811*X)))
          DUI=-DX*(4.94757+X*(1.081603*2.+.034811*3.*X))
          DUI=DUI*(UI-8.)/UI
          UII=9.+EXP(2.173152+X*(-1.295327+X*(1.940395+.813303*X)))
          DUII=DX*(-1.295327+X*(1.940395*2.+.813303*3.*X))
          DUII=DUII*(UII-9.)/UII
          UIII=8.+EXP(-.567398+X*(-3.383369+X*(.799911+.554397*X)))
          DUIII=DX*(-3.383369+X*(.799911*2.+.554397*3.*X))
          DUIII=DUIII*(UIII-8.)/UIII
       CASE(70)
          ! Ytterbium
          UI=1.+EXP(-2.350549-X*(6.688837+X*(1.93869+.269237*X)))
          DUI=-DX*(6.688837+X*(1.93869*2.+.269237*3.*X))
          DUI=DUI*(UI-1.)/UI
          UII=2.+EXP(-3.047465-X*(7.390444+X*(2.355267+.44757*X)))
          DUII=-DX*(7.390444+X*(2.355267*2.+.44757*3.*X))
          DUII=DUII*(UII-2.)/UII
          UIII=1.+EXP(-6.192056-X*(10.560552+X*(4.579385+.940171*X)))          
          DUIII=-DX*(10.560552+X*(4.579385*2.+.940171*3.*X))
          DUIII=DUIII*(UIII-1.)/UIII
       CASE(71)
          ! Lutetium
          UI=4.+EXP(1.537094+X*(-1.140264+X*(.608536+.193362*X)))
          DUI=DX*(-1.140264+X*(.608536*2.+.193362*3.*X))
          DUI=DUI*(UI-4.)/UI
          UII=MAX(1.,0.66+1.52E-4*T)
          IF(UII.NE.1.)DUII=1.52E-4/UII
          IF (T.GT.5250.) UII=-1.09+4.86E-4*T
          IF (T.GT.5250.) DUII=4.86E-4/UII
          UIII=5.
       CASE(72)
          ! Hafnium
          UI=4.1758+Y*(.407+Y*(.57862-Y*(.072887-Y*3.6848E-3)))
          DUI=1.E-3*(.407+Y*(.57862*2.-Y*(.072887*3.-Y*3.6848E-3*4.)))/UI
          UII=-2.979+3.095E-3*T
          DUII=3.095E-3/UII
          UIII=30.
! NEW: apy 20170605 BEGINS
       CASE(73)
          ! Tantalum
          UI=3.0679+Y*(.81776+Y*(.34936+Y*(7.4861E-3+Y*3.0739E-4)))
          DUI=1.E-3*(.81776+Y*(.34936*2.+Y*(7.4861E-3*3.+Y*3.0739E-4*4.)))/UI
          UII=1.6834+Y*(2.0103+Y*(.56443-Y*(.031036-Y*8.9565E-4)))
          DUII=1.E-3*(2.0103+Y*(.56443*2.-Y*(.031036*3.-Y*8.9565E-4*4.)))/UII
          UIII=15.
       CASE(74)
          ! Tungsten
          UI=.3951+Y*(-.25057+Y*(1.4433+Y*(-.34373+Y*(.041924-Y*1.84E-3))))
          IF(T.GT.1.2E4) UI=23.
          IF(T.LE.1.2E4) DUI=1.E-3*(-.25057+Y*(1.4433*2.+Y* &
              (-.34373*3.+Y*(.041924*4.-Y*1.84E-3*5.))))/UI
          UII=1.055+Y*(1.0396+Y*(.3303-Y*(8.4971E-3-Y*5.5794E-4)))
          DUII=1.E-3*(1.0396+Y*(.3303*2.-Y*(8.4971E-3*3.-Y*5.5794E-4*4.)))/UII
          UIII=20.
       CASE(75)
          ! Rhenium
          UI=5.5671+Y*(.72721+Y*(-.42096+Y*(.09075-Y*3.9331E-3)))
          IF(T.GT.1.2E4) UI=29.
          IF(T.LE.1.2E4) DUI=1.E-3*(.72721+Y*(-.42096*2.+Y* &
              (.09075*3.-Y*3.9331E-3*4.)))/UI
          UII=6.5699+Y*(.59999+Y*(-.28532+Y*(.050724-Y*1.8544E-3)))
          IF(T.GT.1.2E4) UII=22.
          IF(T.LE.1.2E4) DUII=1.E-3*(.59999+Y*(-.28532*2.+Y* &
              (.050724*3.-Y*1.8544E-3*4.)))/UII
          UIII=20.
       CASE(76)
          ! Osmium
          UI=8.6643+Y*(-.32516+Y*(.68181-Y*(.044252-Y*1.9975E-3)))
          DUI=1.E-3*(-.32516+Y*(.68181*2.-Y*(.044252*3.-Y*1.9975E-3*4.)))/UI
          UII=9.7086+Y*(-.3814+Y*(.65292-Y*(.064984-Y*2.8792E-3)))
          DUII=1.E-3*(-.3814+Y*(.65292*2.-Y*(.064984*3.-Y*2.8792E-3*4.)))/UII
          UIII=10.
       CASE(77)
          ! Iridium
          UI=11.07+Y*(-2.412+Y*(1.9388+Y*(-.34389+Y*(.033511-1.3376E-3*Y))))
          IF(T.GT.1.2E4) UI=30.
          IF(T.LE.1.2E4)DUI=1.E-3*(-2.412+Y*(1.9388*2.+Y*(-.34389*3.+Y* &
              (.033511*4.-1.3376E-3*5.*Y))))/UI
          UII=15.
          UIII=20.
       CASE(78)
          ! Plutonium
          UI=16.4+1.27E-3*T
          DUI=1.27E-3/UI
          UII=6.5712+Y*(-1.0363+Y*(.57234-Y*(.061219-2.6878E-3*Y)))
          DUII=1.E-3*(-1.0363+Y*(.57234*2.-Y*(.061219*3.-Y*4*2.6878E-3)))
          DUII=DUII/UII
          UIII=15.
       CASE(79)
          ! Gold
          UI=1.24+2.79E-4*T
          DUI=2.79E-4/UI
          UII=1.0546+Y*(-.040809+Y*(2.8439E-3+Y*1.6586E-3))
          DUII=1.E-3*(-.040809+Y*(2.8439E-3*2+Y*1.6586E-3*3))/UII
          UIII=7.
       CASE(80)
          ! Mercury
          UI=1.0
          UII=2.0
          UIII=AMAX1(1.,.669+3.976E-5*T)
          IF(UIII.NE.1.)DUIII=3.976E-5/UIII
       CASE(81)
          ! Thallium
          UI=AMAX1(2.,0.63+3.35E-4*T)
          IF(UI.NE.2.)DUI=3.35E-4/UI
          UII=1.0
          UIII=2.0
       CASE(82)
          ! Lead
          UI=AMAX1(1.,0.42+2.35E-4*T)
          IF(UI.NE.1.)DUI=2.35E-4/UI
          IF(T.GT.6125.) UI=-1.2+5.E-4*T
          IF(T.GT.6125.)DUI=5.E-4/UI
          UII=AMAX1(2.,1.72+7.9E-5*T)
          IF(UII.NE.2.)DUII=7.9E-5/UII
          UIII=1.0
       CASE(83)
          ! Bismuth
          UI=2.78+2.87E-4*T
          DUI=2.87e-4/UI
          UII=AMAX1(1.,.37+1.41E-4*T)
          IF(UII.NE.1.)DUII=1.41E-4/UII
          UIII=2.5 ! APPROXIMATELY apy: O_O
       CASE(84)
          ! Polonium
          UI=5.
          UII=5.
          UIII=4.
       CASE(85)
          ! Astatine
          UI=4.
          UII=6.
          UIII=6.
       CASE(86)
          ! Radon
          UI=1.
          UII=4.
          UIII=6.
       CASE(87)
          ! Francium
          UI=2.
          UII=1.
          UIII=4.5
       CASE(88)
          ! Radium
          UI=1.
          UII=2.
          UIII=1.
       CASE(89)
          ! Actinium
          UI=6.
          UII=3.
          UIII=7.
       CASE(90)
          ! Thorium
          UI=8.
          UII=8.
          UIII=8.
       CASE(91)
          ! Protactinium
          UI=50.
          UII=50.
          UIII=50.
       CASE(92)
          ! Uranium
          UI=25.
          UII=25.
          UIII=25.
          !
    END SELECT
    !
  END SUBROUTINE PARTITION_FUNCTION
  !
  !------------------------------------------------
  !
  SUBROUTINE REFRAX (W,RINDEX)
    ! Based of refrax.f from SIR distribution (Ruiz Cobo & Del Toro Iniesta 1992)
    ! Inout/Output
    REAL(DP), INTENT(IN)           :: W            ! Wavelength in microns
    REAL(DP), INTENT(INOUT)        :: RINDEX
    ! Internal: tabulated values
    REAL(DP), PARAMETER            :: T = 15.0D0   ! temperature in Celcius
    REAL(DP), PARAMETER            :: P = 760.0D0  ! 1 atm pressure in mmHg
    REAL(DP), PARAMETER            :: H = 0.0D0    ! Humidity
    ! Internal
    REAL(DP)                       :: S, R, Y, WVSP, WVP, ALPHA
    !
    S=1.0D0/W**2.0D0
    R=6.4328D-5+2.94981D-2/(146.0D0-S)+2.554D-4/(41.0D0-S)
    Y=1.0D0/(273.15D0+T)
    ALPHA=3.67D-3+3.3D-4*DEXP((0.19D0-W)/0.12D0)
    IF (T.LE.0.0D0) THEN
       ! Valid only when T < 0 !!
       WVSP=10.0D0**(77.4021323D0+Y*(-9.6982D4+Y*(5.50733046D7&
            +Y*(-1.70804596D10+Y*(2.96751446D12-Y&
            *(2.73866936D14-1.04883576D16*Y))))))
    ELSE
       ! Valid only when T > 0 !!
       WVSP=10.0D0**(-55.1132754D0+Y*(1.19551822D5+Y*(-9.70108858D7&
            +Y*(4.12856913D10+Y*(-9.87931835D12+Y&
            *(1.25868645D15-6.66885272D16*Y))))))
    ENDIF
    !
    WVP=1.0D-2*H*WVSP
    RINDEX=1.0D0+1.31579D-3*R*P/(1.0D0+ALPHA*(T-1.0D0)/(1.0D0+15.0D0*ALPHA))&
        -5.49D-8*WVP/(1.0D0+ALPHA*T)
    !
  END SUBROUTINE REFRAX
  !
  !================================================
  !
END MODULE ATOM_DATABASE
!
