MODULE SPLINES
  !
  ! J M Borrero
  ! July 23, 2013
  ! KIS, Freiburg
  ! Based on Numerical Recipes
  !
  USE CONS_PARAM
  !
  CONTAINS
  !
  ! splie2
  ! splin2
  ! spline
  ! splint
  !
  ! -----------------
  ! Subroutine splie2
  ! -----------------
  SUBROUTINE SPLIE2(X2A,YA,M,N,Y2A)
    IMPLICIT NONE
    INTEGER,       INTENT(IN)    :: M,N
    REAL(DP),      INTENT(IN)    :: X2A(N), YA(M,N)
    REAL(DP),      INTENT(INOUT) :: Y2A(M,N)
    !
    REAL(DP),          PARAMETER :: EPS=1.E30
    INTEGER,           PARAMETER :: NN=1000
    REAL(DP)                     :: YTMP(NN),Y2TMP(NN)
    INTEGER                      :: J,K
    !
!INCLUDE CODE HERE
    !
    RETURN
    !
  END SUBROUTINE SPLIE2
  ! -----------------
  ! Subroutine splin2
  ! -----------------
  SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
    IMPLICIT NONE
    INTEGER,        INTENT(IN)  :: M,N
    REAL(DP),       INTENT(IN)  :: X1A(M),X2A(N),YA(M,N),Y2A(M,N),X1,X2 
    REAL(DP),       INTENT(INOUT) :: Y
    !
    REAL(DP),         PARAMETER :: EPS=1.E30
    INTEGER,          PARAMETER :: NN=1000
    REAL(DP)                    :: YTMP(NN),Y2TMP(NN),YYTMP(NN)
    INTEGER                     :: J,K
    !
!INCLUDE CODE HERE
    !
  END SUBROUTINE SPLIN2
  ! -----------------
  ! Subroutine spline
  ! -----------------
  SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
    IMPLICIT NONE
    INTEGER,        INTENT(IN)  :: N
    REAL(DP),       INTENT(IN)  :: X(N),Y(N), YP1, YPN
    REAL(DP),       INTENT(INOUT) :: Y2(N)
    !
    REAL(DP)                    :: U(N), P, QN, UN, SIG
    INTEGER                     :: I,K
    !
!INCLUDE CODE HERE
    !
  END SUBROUTINE SPLINE
  !
  ! Subroutine splint
  !
  SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
    IMPLICIT NONE
    INTEGER,                 INTENT(IN)  :: N
    REAL(DP),                INTENT(IN)  :: XA(N),YA(N),Y2A(N),X
    REAL(DP),                INTENT(INOUT) :: Y
    !
    REAL(DP)                             :: A,B,H
    INTEGER                              :: K,KLO, KHI
    !
!INCLUDE CODE HERE
    !
  END SUBROUTINE SPLINT
  !
END MODULE SPLINES
