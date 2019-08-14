!
MODULE LOG
  !
  !================================================
  !
  ! J M Borrero 
  ! Jan 11, 2009
  ! HAO-NCAR & MPS-Lindau for HMI-Stanford
  !
  USE CONS_PARAM, ONLY: DP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER,   PARAMETER    :: LOGUNIT=108
  !
  PUBLIC :: LOG_INIT
  PUBLIC :: LOGW_CHAR_REAL
  PUBLIC :: LOGW_CHAR
  PUBLIC :: LOGW_CHAR_CHAR
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! log_init
  ! log_close
  ! logw_char
  ! logw_char_real
  ! logw_char_int
  ! logw_char_char
  ! logw_inv_details
  ! logw_summary
  !
  !------------------------------------------------
  !
  SUBROUTINE LOG_INIT
    OPEN(UNIT=LOGUNIT,FILE='info.log',STATUS='UNKNOWN')
    CALL LOGW_CHAR('TEST version 1.00: Started !!')
    CALL LOGW_CHAR('==========================================================')
    CALL LOGW_CHAR('This is a Free Program originally developed by J M Borrero')
    CALL LOGW_CHAR('----------------------------------------------------------')
    !CALL LOGW_CHAR_INT(22,'Number of Fits files: ',NFIT)
    !CALL LOGW_CHAR_INT(26,'Number of Pixel per slit: ',NPIX)
    !CALL LOGW_CHAR_CHAR(17,'Input Directory: ',LEN_TRIM(WORKDIR),TRIM(WORKDIR))
    !CALL LOGW_CHAR_CHAR(18,'Output Directory: ',LEN_TRIM(OUTDIR),TRIM(OUTDIR))
    CALL LOGW_CHAR(' ')
    !CALL LOGW_CHAR_CHAR(22,'Input Fits file list: ',LEN_TRIM(WORKDIR)+LEN_TRIM(LIST) &
    !     ,TRIM(WORKDIR)//TRIM(LIST))
    !CALL LOGW_CHAR_CHAR(21,'Output File results: ',LEN_TRIM(OUTDIR)+LEN_TRIM(OUTFIL) &
    !     ,TRIM(OUTDIR)//TRIM(OUTFIL))
    !IF (QUICKLOOK.EQ..TRUE.) THEN
    !   CALL LOGW_CHAR(60,'----------------------------------------------------------')
    !   CALL LOGW_CHAR(41,'WARNING: This is a QuickLook inversion !!')
    !   CALL LOGW_CHAR_REAL(43,'Minimum Polarization to perform inversion: ',TREPOL)
    !   CALL LOGW_CHAR(60,'----------------------------------------------------------')
    !ENDIF
  END SUBROUTINE LOG_INIT
  !
  !------------------------------------------------
  !
  SUBROUTINE LOG_CLOSE
    CLOSE(UNIT=LOGUNIT)
  END SUBROUTINE LOG_CLOSE
  !
  !------------------------------------------------
  !
  SUBROUTINE LOGW_CHAR(MESS)
    CHARACTER(*)           :: MESS
    WRITE(UNIT=LOGUNIT,FMT='(A)') MESS
  END SUBROUTINE LOGW_CHAR
  !
  !------------------------------------------------
  !
  SUBROUTINE LOGW_CHAR_REAL(MESS,REA)
    CHARACTER(*)           :: MESS
    REAL(DP)                 :: REA
    WRITE(UNIT=LOGUNIT,FMT='(A,1x,F12.4)') MESS,REA
  END SUBROUTINE LOGW_CHAR_REAL
  !
  !------------------------------------------------
  !
  SUBROUTINE LOGW_CHAR_INT (MESS,INT)
    CHARACTER(*)           :: MESS
    INTEGER                  :: INT
    WRITE(UNIT=LOGUNIT,FMT='(A,1x,I8)') MESS,INT
  END SUBROUTINE LOGW_CHAR_INT
  !
  !------------------------------------------------
  !
  SUBROUTINE LOGW_CHAR_CHAR(MESS1,MESS2)
    CHARACTER(*)           :: MESS1
    CHARACTER(*)           :: MESS2
    WRITE(UNIT=LOGUNIT,FMT='(A,1x,A)') MESS1,MESS2
  END SUBROUTINE LOGW_CHAR_CHAR
  !
  !------------------------------------------------
  !
  SUBROUTINE LOGW_INV_DETAILS
    INTEGER              :: I, NF
    CHARACTER*11         :: PAR
    NF=0
    CALL LOGW_CHAR('----------------------------------------------------------')
    DO I=1,12
          SELECT CASE (I)
          CASE(1)
             PAR='Absopt.Coef'
          CASE(2)
             PAR='Inclination'
          CASE(3)
             PAR='Azimuth.Ang'
          CASE(4)
             PAR='Damping.Par'
          CASE(5)
             PAR='Doppl.Width'
          CASE(6)
             PAR='Magnetic Fi'
          CASE(7)
             PAR='LOS Velocit'
          CASE(8)
             PAR='Con.SourceF'
          CASE(9)
             PAR='SourceF.Gra'
          CASE(10)
             PAR='Mag.Fil.Fac'
          CASE(11)
             PAR='Macroturb.V'
          CASE(12)
             PAR='NonMag.VLOS'
          END SELECT
       !IF (FREE(I).EQ..TRUE.) THEN
       !   CALL LOGW_CHAR_CHAR(11,PAR,6,' : YES')
       !   NF=NF+1
       !ENDIF
       !IF (FREE(I).EQ..FALSE.) CALL LOGW_CHAR_CHAR(11,PAR,5,' : NO')
    ENDDO
    CALL LOGW_CHAR_INT('Total Num of Free Parameters: ',NF)
    CALL LOGW_CHAR('----------------------------------------------------------')

  END SUBROUTINE LOGW_INV_DETAILS
  !
  !------------------------------------------------
  !
  SUBROUTINE LOGW_SUMMARY(TIME)!,N10,N25)!,MEANC)
    !INTEGER         :: N10, N25
    !REAL(DP)        :: MEANC, MINC!, MAXC
    REAL(4)         :: TIME
    CALL LOGW_CHAR('----------------------------------------------------------')
    CALL LOGW_CHAR_REAL('Time [sec] employed: ',DBLE(TIME))
    !CALL LOGW_CHAR_REAL(27,'Inversion Speed [pix/sec]: ',DBLE(NFIT*NPIX)/DBLE(TIME))
    !IF (QUICKLOOK.EQ..FALSE.) THEN
    !   CALL LOGW_CHAR(60,'----------------------------------------------------------')
    !   CALL LOGW_CHAR_INT(32,'Number of pixels with Chi2 > 10 : ',N10)
    !   CALL LOGW_CHAR_INT(32,'Number of pixels with Chi2 > 25 : ',N25)
    !   CALL LOGW_CHAR_REAL(20,'Mean Chi2 achieved: ',MEANC)
    !   CALL LOGW_CHAR(1,'')
    !ENDIF
    CALL LOGW_CHAR('==========================================================')
    CALL LOGW_CHAR('TEST version 1.00: Finished !!')
  ENDSUBROUTINE LOGW_SUMMARY
  !
  !================================================
  !
END MODULE LOG
!
