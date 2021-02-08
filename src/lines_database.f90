!
MODULE LINES_DATABASE
  !
  !================================================
  !
  USE user_mpi, ONLY: mpi__myrank
  USE CONS_PARAM, ONLY: SP, DP
  USE FORWARD_PARAM
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_1D_DP, ALLOCATE_2D_IP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  REAL(DP), ALLOCATABLE, PRIVATE   :: ELOWLIM(:), EUPPLIM(:)
  !
  CHARACTER*800, ALLOCATABLE, DIMENSION(:)    :: LDBASE_READ_LINES
  !
  PUBLIC :: INIT_LINEDATABASE_VARS
  !PUBLIC :: READ_LINES_DATABASE
  PUBLIC :: NREAD_LINES_DATABASE
  PUBLIC :: FINI_LINEDATABASE_VARS
  PUBLIC :: NCOLUMNS
  !
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! init_linedatabase_vars
  ! read_linesdatabase_row
  ! fini_linedatabase_vars
  ! linesdatabase_get13
  ! linesdatabase_get8
  ! read_lines_database
  ! elecconf
  ! optical_elec_tran
  !
  !------------------------------------------------
  !
  SUBROUTINE INIT_LINEDATABASE_VARS(NL)
    !
    INTEGER, INTENT(IN)                           :: NL
    !
    CALL ALLOCATE_1D_IP(LINE_POS, NL, 'ALLOCATE line_pos IN lines_database.f90')
    CALL ALLOCATE_1D_IP(LINE_NUM, NL, 'ALLOCATE line_num IN lines_database.f90')
    CALL ALLOCATE_1D_DP(LINE_L0, NL, 'ALLOCATE line_l0 IN lines_database.f90')
    CALL ALLOCATE_1D_IP(LINE_ZN, NL, 'ALLOCATE line_zn IN lines_database.f90')
    CALL ALLOCATE_1D_IP(LINE_ION, NL, 'ALLOCATE line_ion IN lines_database.f90')
    CALL ALLOCATE_1D_DP(LOGGF, NL, 'ALLOCATE loggf IN lines_database.f90')
    CALL ALLOCATE_1D_DP(EPLOW, NL, 'ALLOCATE eplow IN lines_database.f90')
    CALL ALLOCATE_1D_DP(ALPHA, NL, 'ALLOCATE alpha IN lines_database.f90')
    CALL ALLOCATE_1D_DP(SIGMA, NL, 'ALLOCATE sigma IN lines_database.f90')
    CALL ALLOCATE_1D_DP(ELOWLIM, NL, 'ALLOCATE elowlim IN lines_database.f90')
    CALL ALLOCATE_1D_DP(EUPPLIM, NL, 'ALLOCATE eupplim IN lines_database.f90')
    CALL ALLOCATE_1D_DP(SL, NL, 'ALLOCATE sl IN lines_database.f90')
    CALL ALLOCATE_1D_DP(SU, NL, 'ALLOCATE su IN lines_database.f90')
    CALL ALLOCATE_1D_DP(LL, NL, 'ALLOCATE ll IN lines_database.f90')
    CALL ALLOCATE_1D_DP(LU, NL, 'ALLOCATE lu IN lines_database.f90')
    CALL ALLOCATE_1D_DP(JL, NL, 'ALLOCATE jl IN lines_database.f90')
    CALL ALLOCATE_1D_DP(JU, NL, 'ALLOCATE ju IN lines_database.f90')
    CALL ALLOCATE_2D_IP(OETRANSITION, NL, 2 &
        , 'ALLOCATE oetransition IN lines_database.f90')
    !
  END SUBROUTINE INIT_LINEDATABASE_VARS
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_LINESDATABASE_ROW(READLINE, I, NCOL, LOW, UPP, SHELL)
    !
    USE user_mpi
    !
    CHARACTER*800, INTENT(IN)                   :: READLINE
    !
    INTEGER, INTENT(IN)                         :: NCOL, I
    !
    CHARACTER*2, INTENT(INOUT)                  :: SHELL
    CHARACTER*7, INTENT(INOUT)                  :: LOW, UPP
    !
    ! Check number of columns found and act accordingly
    IF (NCOL.EQ.13) THEN
      !
      !WRITE(*,*) '13 columns format'
      CALL LINESDATABASE_GET13(READLINE,LINE_NUM(I), LINE_ZN(I), LINE_ION(I) &
          , LINE_L0(I), LOW, UPP, LOGGF(I), EPLOW(I), ALPHA(I), SIGMA(I) &
          , ELOWLIM(I), EUPPLIM(I), SHELL)
!IF (mpi__myrank.eq.0) PRINT*, 'Here: ', LINE_NUM(I), LINE_ZN(I), LINE_ION(I) &
!          , LINE_L0(I), LOW, UPP, LOGGF(I), EPLOW(I), ALPHA(I), SIGMA(I) &
!          , ELOWLIM(I), EUPPLIM(I), SHELL

    ELSE IF (NCOL.EQ.8) THEN
      !WRITE(*,*) '8 columns format'
      CALL LINESDATABASE_GET8(READLINE,LINE_NUM(I), LINE_ZN(I), LINE_ION(I) &
          , LINE_L0(I), LOW, UPP &
          , LOGGF(I), EPLOW(I))
    ELSE
      WRITE(*,*) 'YOU SHOULD NOT BE HERE!, I AM READ_LINESDATABASE_ROW SUBROUTINE IN lines_database.f90'
      WRITE(*,*) 'YOU SHOULD NOT BE HERE!, I AM READ_LINESDATABASE_ROW SUBROUTINE IN lines_database.f90'
     
    ENDIF
    !
  END SUBROUTINE READ_LINESDATABASE_ROW
  !
  !------------------------------------------------
  !
  SUBROUTINE FINI_LINEDATABASE_VARS()
    !
    DEALLOCATE(LINE_POS)
    DEALLOCATE(LINE_NUM)
    DEALLOCATE(LINE_L0)
    DEALLOCATE(LINE_ZN)
    DEALLOCATE(LINE_ION)
    DEALLOCATE(LOGGF)
    DEALLOCATE(EPLOW)
    DEALLOCATE(ALPHA)
    DEALLOCATE(SIGMA)
    DEALLOCATE(ELOWLIM)
    DEALLOCATE(EUPPLIM)
    DEALLOCATE(SL)
    DEALLOCATE(SU)
    DEALLOCATE(LL)
    DEALLOCATE(LU)
    DEALLOCATE(JL)
    DEALLOCATE(JU)
    DEALLOCATE(OETRANSITION)
    !
  END SUBROUTINE FINI_LINEDATABASE_VARS
  !
  !------------------------------------------------
  !
  INTEGER FUNCTION NCOLUMNS(INLINE)
     !
     IMPLICIT NONE
     !
     CHARACTER(*), INTENT(IN)       :: INLINE
     INTEGER                         :: I, NCOL
     !
     LOGICAL                         :: PRECHAR, SPECIALCHAR
     !
     NCOL=0
     !
     PRECHAR=.FALSE.
     SPECIALCHAR=.FALSE.
     !
     DO I=1,LEN(INLINE)

       IF (INLINE(I:I).NE.' ') THEN
         ! AM I SPECIAL?
         IF (INLINE(I:I).EQ.',') THEN
           SPECIALCHAR=.TRUE.
           CYCLE
         ENDIF
         !WAS THE PREVIOUS THING A NON EMPTY SPACE?
         IF ((PRECHAR.EQV..FALSE.).AND.(SPECIALCHAR.EQV..FALSE.)) THEN
!PRINT*, I
             !PRINT*, 'LINEA:', INLINE(1:I)
             NCOL=NCOL+1
             !PRINT*, NCOL
! CAUTION, THIS WAY, WE ARE GETTING WHERE A NEW ENTRY OF THE LINE BEGINS,...
! ...NOT DIRECTLY IDENTIFYING THE COLUMNS THEMSELVES
         ELSE
           SPECIALCHAR=.FALSE.
         ENDIF
         ! UPDATE PRECHAR:
         PRECHAR=.TRUE.
       ELSE
         PRECHAR=.FALSE.
       ENDIF

     ENDDO
     NCOLUMNS=NCOL
     !
  END FUNCTION NCOLUMNS
  !
  !------------------------------------------------
  !
  SUBROUTINE LINESDATABASE_GET13(READLINE,INUM,IZN,IION,IL0,LOW,UPP,ILOGGF,IEPLOW&
      ,IALPHA,ISIGMA,IELOWLIM,IEUPPLIM,SHELL)
    !
    CHARACTER(*), INTENT(IN)                    :: READLINE
    !
    CHARACTER*2, INTENT(INOUT)                  :: SHELL
    CHARACTER*7, INTENT(INOUT)                  :: LOW, UPP
    !
    INTEGER                                     :: INUM, IZN, IION
    !
    REAL(DP)                                    :: IL0, ILOGGF, IEPLOW, IALPHA, ISIGMA
    REAL(DP)                                    :: IELOWLIM,IEUPPLIM
    !
    INTEGER                                     :: I, PRE, START, INCOL, J
    INTEGER                                     :: ITERM
    CHARACTER*7                                 :: ATERM
    REAL(DP)                                    :: DTERM
    REAL(SP)                                    :: RTERM
    LOGICAL                                     :: ISREAL
    !
    ! Column order:
    !    Col 1: Line database ID
    !    Col 2: Atomic number
    !    Col 3: Ionisation stage
    !    Col 4: Transition central wavelength in angstrom
    !    Col 5: Lower level electronic term
    !    Col 6: Upper level electronic term
    !    Col 7: gf decimal logarithm
    !    Col 8: Lower level energy in eV
    !    Col 9: ALPHA (CHANGE)
    !    Col 10: SIGMA (CHANGE)
    !    Col 11: SHELL
    !    Col 12: LOWER ENERGY ABE (CHANGE)
    !    Col 13: UPPER ENERGY ABE (CHANGE)
    INCOL=0
    PRE=0
    START=1
    DO I=1,LEN(READLINE)
      IF(READLINE(I:I).EQ.' ') THEN
        IF (PRE.EQ.1) THEN
          ! Give appropiate format depending on the column:
          IF ((INCOL.EQ.1).OR.(INCOL.EQ.2).OR.(INCOL.EQ.3)) THEN
            READ(READLINE(START:I-1), '(I7)') ITERM
          ELSE IF ((INCOL.EQ.5).OR.(INCOL.EQ.6).OR.(INCOL.EQ.11)) THEN
            READ(READLINE(START:I-1), '(A)') ATERM
          ELSE
            ! Check if the number is given as real or integer:
            ISREAL=.FALSE.
            DO J=1,LEN(READLINE(START:I-1))
              IF(READLINE(START+J-1:START+J-1).EQ.'.') THEN
                ISREAL=.TRUE.
                EXIT
              ENDIF
            ENDDO
!IF (mpi__myrank.EQ.0) PRINT*, LEN(READLINE(START:I-1)), J, READLINE(START:START+J) &
    !,READLINE(START:I-1)
            IF (ISREAL.EQV..TRUE.) THEN
              READ(READLINE(START:I-1), '(f18.8)') RTERM
              DTERM=DBLE(RTERM)
              !IF (LEN(READLINE(START:I-1)).EQ.9) THEN
                !IF (J.EQ.5) THEN
                  !READ(READLINE(START:I-1), '(f9.4)') RTERM
!IF (mpi__myrank.EQ.0) PRINT*, RTERM
                !ENDIF
              !ENDIF
!IF (mpi__myrank.EQ.0) PRINT*, RTERM
              DTERM=DBLE(RTERM)
!IF (mpi__myrank.EQ.0) PRINT*, RTERM, DTERM
            ELSE
              READ(READLINE(START:I-1), '(i14)') ITERM
              DTERM=REAL(ITERM)
            ENDIF
!IF (mpi__myrank.EQ.0) PRINT*, RTERM, DTERM
            DTERM=DBLE(NINT(DTERM*1.0D6,KIND=8))*1.0D-6
!IF (mpi__myrank.EQ.0) PRINT*, RTERM, DTERM
          ENDIF
          !PRINT*, READLINE(START:I-1), INCOL
          ! ASSIGN TERM TO EACH COLUMN:
          SELECT CASE(INCOL)
             CASE(1)
                INUM=ITERM
             CASE(2)
                IZN=ITERM
             CASE(3)
                IION=ITERM
             CASE(4)
!IF (mpi__myrank.EQ.0) PRINT*, 'L0: ', DTERM, READLINE(START:I-1)
                IL0=DBLE(DTERM)
             CASE(5)
                LOW=ATERM
             CASE(6)
                UPP=ATERM
             CASE(7)
!IF (mpi__myrank.EQ.0) PRINT*, 'LOGGF: ', DTERM, READLINE(START:I-1)
                ILOGGF=DBLE(DTERM)
             CASE(8)
                IEPLOW=DBLE(DTERM)
             CASE(9)
                IALPHA=DBLE(DTERM)
             CASE(10)
                ISIGMA=DBLE(DTERM)
             CASE(11)
                SHELL=ATERM(1:3)
             CASE(12)
                IELOWLIM=DBLE(DTERM)
             CASE(13)
                IEUPPLIM=DBLE(DTERM)
          END SELECT
          !END ASSIGNMENT
        ENDIF
        PRE=0
        CYCLE
      ENDIF
      !
      IF (PRE.EQ.0) THEN
        INCOL=INCOL+1
        PRE=1
        START=I
      ENDIF
    ENDDO
    !
  END SUBROUTINE LINESDATABASE_GET13
  !
  !------------------------------------------------
  !
  SUBROUTINE LINESDATABASE_GET8(READLINE,INUM,IZN,IION,IL0,LOW,UPP,ILOGGF,IEPLOW)
    !
    CHARACTER(*), INTENT(IN)                    :: READLINE
    !
    CHARACTER*7, INTENT(INOUT)                  :: LOW, UPP
    !
    INTEGER                                     :: INUM, IZN, IION
    !
    REAL(DP)                                    :: IL0, ILOGGF, IEPLOW
    !
    INTEGER                                     :: I, PRE, START, INCOL
    INTEGER                                     :: ITERM
    CHARACTER*7                                 :: ATERM
    REAL(DP)                                    :: DTERM
    !
    ! Column order:
    !    Col 1: Line database ID
    !    Col 2: Atomic number
    !    Col 3: Ionisation stage
    !    Col 4: Transition central wavelength in angstrom
    !    Col 5: Lower level electronic term
    !    Col 6: Upper level electronic term
    !    Col 7: gf decimal logarithm
    !    Col 8: Lower level energy in eV
    INCOL=0
    PRE=0
    START=0
    DO I=1,LEN(READLINE)
      IF(READLINE(I:I).EQ.' ') THEN
        IF (PRE.EQ.1) THEN
          ! Give appropiate format depending on the column:
          IF ((INCOL.EQ.1).OR.(INCOL.EQ.2).OR.(INCOL.EQ.3)) THEN
            READ(READLINE(START:I), '(I7)') ITERM
          ELSE IF ((INCOL.EQ.5).OR.(INCOL.EQ.6)) THEN
            READ(READLINE(START:I), '(A)') ATERM
          ELSE
            READ(READLINE(START:I), '(e18.8)') DTERM
            DTERM=DBLE(NINT(DTERM*1.D4))*1.D-4
          ENDIF
          PRINT*, READLINE(START:I), INCOL
          ! ASSIGN TERM TO EACH COLUMN:
          SELECT CASE(INCOL)
             CASE(1)
                INUM=ITERM
             CASE(2)
                IZN=ITERM
             CASE(3)
                IION=ITERM
             CASE(4)
                IL0=DBLE(DTERM)
             CASE(5)
                LOW=ATERM
             CASE(6)
                UPP=ATERM
             CASE(7)
                ILOGGF=DBLE(DTERM)
             CASE(8)
                IEPLOW=DBLE(DTERM)
          END SELECT
          !END ASSIGNMENT
        ENDIF
        PRE=0
        CYCLE
      ENDIF
      !
      IF (PRE.EQ.0) THEN
        INCOL=INCOL+1
        PRE=1
        START=I
      ENDIF
    ENDDO
    !
  END SUBROUTINE LINESDATABASE_GET8
  !
  !------------------------------------------------
  !
  SUBROUTINE NREAD_LINES_DATABASE
    !
    USE LOG
    USE DAMPING, ONLY: GET_NEFF, GET_COLLISIONAL_PARAM
    USE CODE_MODES, ONLY: LINEPATH
    USE user_mpi, ONLY: mpi__ierror, mpi_character, mpi_comm_world, mpi_integer
    !
    IMPLICIT NONE
    INTEGER                                     :: I, NL, IERR, FLAG
    CHARACTER*2                                 :: SHELL
    CHARACTER*7                                 :: LOW, UPP
    CHARACTER*100                               :: FIN
    REAL(DP)                                    :: NLOW_EFF, NUPP_EFF
    !
    INTEGER                                     :: NCOL
    CHARACTER*800                               :: READLINE
    !
    ! Does lines_database exist ?
    !
    IF (mpi__myrank.eq.0) THEN
      !
      OPEN(UNIT=1, FILE=TRIM(LINEPATH)//'lines_database.dat' &
          , STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
      IF (IERR.NE.0) THEN
  PRINT*, TRIM(LINEPATH)
         PRINT*,'The file containing database about spectral lines (lines_database.dat)'
         PRINT*,'could not be fonud in the directory of the source code. STOP'
         CLOSE(UNIT=1)
         STOP
      ENDIF
      !
      ! Determine number of lines inside
      !
      NL=0
      IERR=0 ! reset ierr
      DO WHILE (IERR.EQ.0)
         READ(UNIT=1, FMT='(A)', IOSTAT=IERR) READLINE
         IF (IERR.EQ.0) NL=NL+1
         IF (IERR.NE.0) EXIT
      ENDDO
      CLOSE(1)
      NUML_DATABASE=NL-1

      CALL MPI_BCAST(NUML_DATABASE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__ierror)
      ALLOCATE(LDBASE_READ_LINES(NUML_DATABASE))

      ! Read the file:
      OPEN(UNIT=1, FILE=TRIM(LINEPATH)//'lines_database.dat' &
          , STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
      DO I=1,NUML_DATABASE
         !
         FLAG = 0
         IERR = 0
         !
         READ(UNIT=1, FMT='(A)', IOSTAT=IERR) LDBASE_READ_LINES(I)
      ENDDO

      READ(UNIT=1,FMT='(A)', IOSTAT=IERR) FIN
      IF (TRIM(FIN).NE.'END') THEN
         PRINT*, TRIM(FIN)
         PRINT*,'File lines_database.dat could not be properly read. Check format. STOP'
         STOP
      ENDIF
      CLOSE(1)
      !
      ! Broadcast info:
      ! 1 byte per character -> use sizeof
      CALL MPI_BCAST(LDBASE_READ_LINES, SIZEOF(LDBASE_READ_LINES), MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi__ierror)

    ELSE
      ! Workers:
      CALL MPI_BCAST(NUML_DATABASE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__ierror)
      ALLOCATE(LDBASE_READ_LINES(NUML_DATABASE))
      !
      ! 1 byte per character -> use sizeof
      CALL MPI_BCAST(LDBASE_READ_LINES, SIZEOF(LDBASE_READ_LINES), MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi__ierror)

    ENDIF
    !
    ! Read atomic data
    ! 
    CALL INIT_LINEDATABASE_VARS(NUML_DATABASE)
    !
!    OPEN(UNIT=1, FILE=TRIM(LINEPATH)//'lines_database.dat' &
!        , STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
    DO I=1,NUML_DATABASE
       !
       FLAG = 0
       IERR = 0
       !
       !READ(UNIT=1, FMT='(A)', IOSTAT=IERR) READLINE
       READLINE=LDBASE_READ_LINES(I)
       NCOL=NCOLUMNS(READLINE)
       !PRINT*, I, NCOL, TRIM(READLINE)
       IF((NCOL.NE.13).AND.(NCOL.NE.8)) THEN
         PRINT*, I, NCOL, READLINE
         PRINT*,'WRONG lines_database.dat FORMAT!!!'
         PRINT*,'STOPPING!!!'
         STOP
       ENDIF
       CALL READ_LINESDATABASE_ROW(READLINE, I, NCOL, LOW, UPP, SHELL)
       !
       LINE_POS(I)=I
       !
       ! Determine SLJ quantum numbers from electronic configuration
       !
       CALL ELECCONF(LOW,SL(I),LL(I),JL(I))
       CALL ELECCONF(UPP,SU(I),LU(I),JU(I))
       !
       ! Determine the transition type for the optical electron
       !
       CALL OPTICAL_ELEC_TRAN(SHELL,OETRANSITION(I,:),FLAG)
       !----------------------------------------
       ! We determine the collisional parameters
       !----------------------------------------
       ! Case 1: parameters provided in lines_database.dat
       IF (ALPHA(I) .NE.0 .AND. SIGMA(I) .NE. 0) THEN
          IF (mpi__myrank.eq.0) THEN
            CALL LOGW_CHAR_REAL('Collisional parameters provided by user for ',LINE_L0(I))
            CALL LOGW_CHAR_REAL('Collisional cross section is ',SIGMA(I))
            CALL LOGW_CHAR_REAL('Temperature parameter is ',ALPHA(I))
          ENDIF
       ENDIF
       ! Case 2: parameters not provided by lines_database.dat
       ! but transition of the optical electron is available
       ! so we attempt an interpolation from Anstee, Barklem and O'Mara tables
       IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) .AND. (FLAG .EQ.0)) THEN
       ! Case 2.1: ALPHA and SIGMA parameters not provided but ELIMS are provided:
          IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) &
              .AND. ( ELOWLIM(I) .NE. 0D0 .AND. EUPPLIM(I) .NE. 0D0) ) THEN
            CALL GET_NEFF(LINE_ZN(I),EPLOW(I),LINE_L0(I),LINE_ION(I) &
                ,NLOW_EFF,NUPP_EFF,ELOWLIM(I),EUPPLIM(I))
          ELSE
       ! Case 2.2: interpolate without ELIM information, n* are then overestimated
             ! Determine effective quantum n numbers
             CALL GET_NEFF(LINE_ZN(I),EPLOW(I),LINE_L0(I),LINE_ION(I),NLOW_EFF,NUPP_EFF)
          ENDIF

          ! Interpolate from tables
          CALL GET_COLLISIONAL_PARAM(NLOW_EFF,NUPP_EFF,OETRANSITION(I,:) &
              ,SIGMA(I),ALPHA(I),IERR)
          ! Check for errors:

          IF (mpi__myrank.eq.0) THEN
            IF (IERR .NE. 0) CALL LOGW_CHAR_REAL( &
                'Collisional parameters not available for ',LINE_L0(I))
            IF (IERR .EQ. 1) CALL LOGW_CHAR_CHAR( &
                SHELL,'transition does not have tabulated tables.')
            IF (IERR .EQ. 2) CALL LOGW_CHAR_REAL( &
                'Lower effective quantum number is outside tables ',NLOW_EFF)
            IF (IERR .EQ. 3) CALL LOGW_CHAR_REAL( &
                'Upper effective quantum number is outside tables ',NUPP_EFF)
            IF (IERR .EQ. 0) THEN
               CALL LOGW_CHAR_REAL('Collisional parameters interpolated for ',LINE_L0(I))
               CALL LOGW_CHAR_REAL('Collisional cross section is ',SIGMA(I))
               CALL LOGW_CHAR_REAL('Temperature parameter is ',ALPHA(I))
            ENDIF
          ENDIF
       ENDIF
       ! Case 3: parameters not provided by lines_database.dat
       ! and transition of the optical electron is not available
       IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) .AND. (FLAG .NE.0)) THEN
          IF (mpi__myrank.eq.0) THEN
            CALL LOGW_CHAR_REAL('Collisional parameters not available for ',LINE_L0(I))
            CALL LOGW_CHAR('For this line we will use Van der Waals broadening')
          ENDIF
       ENDIF
    ENDDO
    !
    IF (mpi__myrank.eq.0) THEN
      CALL LOGW_CHAR('File containing database of spectral lines (lines_database.dat) has been correctly read.')
    ENDIF
    DEALLOCATE(LDBASE_READ_LINES)

    !-------------------------------------
  END SUBROUTINE NREAD_LINES_DATABASE
  !
  !------------------------------------------------
  !
!!!!  SUBROUTINE READ_LINES_DATABASE
!!!!    !
!!!!    USE LOG
!!!!    USE DAMPING, ONLY: GET_NEFF, GET_COLLISIONAL_PARAM
!!!!    USE CODE_MODES, ONLY: LINEPATH
!!!!    !
!!!!    IMPLICIT NONE
!!!!    INTEGER                                     :: I, NL, IERR, FLAG
!!!!    CHARACTER*2                                 :: SHELL
!!!!    CHARACTER*7                                 :: LOW, UPP
!!!!    CHARACTER*100                               :: FIN
!!!!    REAL(DP)                                    :: NLOW_EFF, NUPP_EFF
!!!!    !
!!!!    INTEGER                                     :: NCOL
!!!!    CHARACTER*800                               :: READLINE
!!!!    !
!!!!    ! Does lines_database exist ?
!!!!    !
!!!!    OPEN(UNIT=1, FILE=TRIM(LINEPATH)//'lines_database.dat' &
!!!!        , STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!!!!    IF (IERR.NE.0) THEN
!!!!PRINT*, TRIM(LINEPATH)
!!!!       PRINT*,'The file containing database about spectral lines (lines_database.dat)'
!!!!       PRINT*,'could not be fonud in the directory of the source code. STOP'
!!!!       CLOSE(UNIT=1)
!!!!       STOP
!!!!    ENDIF
!!!!    !
!!!!    ! Determine number of lines inside
!!!!    !
!!!!    NL=0
!!!!    IERR=0 ! reset ierr
!!!!    DO WHILE (IERR.EQ.0)
!!!!       READ(UNIT=1, FMT='(A)', IOSTAT=IERR) READLINE
!!!!       IF (IERR.EQ.0) NL=NL+1
!!!!       IF (IERR.NE.0) EXIT
!!!!    ENDDO
!!!!    CLOSE(1)
!!!!    NUML_DATABASE=NL-1
!!!!    !
!!!!    ! Read atomic data
!!!!    ! 
!!!!    CALL INIT_LINEDATABASE_VARS(NL)
!!!!    !
!!!!    OPEN(UNIT=1, FILE=TRIM(LINEPATH)//'lines_database.dat' &
!!!!        , STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!!!!    DO I=1,NL-1
!!!!       !
!!!!       FLAG = 0
!!!!       IERR = 0
!!!!       !
!!!!       READ(UNIT=1, FMT='(A)', IOSTAT=IERR) READLINE
!!!!       NCOL=NCOLUMNS(READLINE)
!!!!       !PRINT*, I, NCOL
!!!!       IF((NCOL.NE.13).AND.(NCOL.NE.8)) THEN
!!!!         PRINT*, I, NCOL, READLINE
!!!!         PRINT*,'WRONG lines_database.dat FORMAT!!!'
!!!!         PRINT*,'STOPPING!!!'
!!!!         STOP
!!!!       ENDIF
!!!!       CALL READ_LINESDATABASE_ROW(READLINE, I, NCOL, LOW, UPP, SHELL)
!!!!       !
!!!!       LINE_POS(I)=I
!!!!       !
!!!!       ! Determine SLJ quantum numbers from electronic configuration
!!!!       !
!!!!       CALL ELECCONF(LOW,SL(I),LL(I),JL(I))
!!!!       CALL ELECCONF(UPP,SU(I),LU(I),JU(I))
!!!!       !
!!!!       ! Determine the transition type for the optical electron
!!!!       !
!!!!       CALL OPTICAL_ELEC_TRAN(SHELL,OETRANSITION(I,:),FLAG)
!!!!       !----------------------------------------
!!!!       ! We determine the collisional parameters
!!!!       !----------------------------------------
!!!!       ! Case 1: parameters provided in lines_database.dat
!!!!       IF (ALPHA(I) .NE.0 .AND. SIGMA(I) .NE. 0) THEN
!!!!          CALL LOGW_CHAR_REAL('Collisional parameters provided by user for ',LINE_L0(I))
!!!!          CALL LOGW_CHAR_REAL('Collisional cross section is ',SIGMA(I))
!!!!          CALL LOGW_CHAR_REAL('Temperature parameter is ',ALPHA(I))
!!!!       ENDIF
!!!!       ! Case 2: parameters not provided by lines_database.dat
!!!!       ! but transition of the optical electron is available
!!!!       ! so we attempt an interpolation from Anstee, Barklem and O'Mara tables
!!!!       IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) .AND. (FLAG .EQ.0)) THEN
!!!!       ! Case 2.1: ALPHA and SIGMA parameters not provided but ELIMS are provided:
!!!!          IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) &
!!!!              .AND. ( ELOWLIM(I) .NE. 0D0 .AND. EUPPLIM(I) .NE. 0D0) ) THEN
!!!!            CALL GET_NEFF(LINE_ZN(I),EPLOW(I),LINE_L0(I),LINE_ION(I) &
!!!!                ,NLOW_EFF,NUPP_EFF,ELOWLIM(I),EUPPLIM(I))
!!!!          ELSE
!!!!       ! Case 2.2: interpolate without ELIM information, n* are then overestimated
!!!!             ! Determine effective quantum n numbers
!!!!             CALL GET_NEFF(LINE_ZN(I),EPLOW(I),LINE_L0(I),LINE_ION(I),NLOW_EFF,NUPP_EFF)
!!!!          ENDIF
!!!!
!!!!          ! Interpolate from tables
!!!!          CALL GET_COLLISIONAL_PARAM(NLOW_EFF,NUPP_EFF,OETRANSITION(I,:) &
!!!!              ,SIGMA(I),ALPHA(I),IERR)
!!!!          ! Check for errors:
!!!!          IF (IERR .NE. 0) CALL LOGW_CHAR_REAL( &
!!!!              'Collisional parameters not available for ',LINE_L0(I))
!!!!          IF (IERR .EQ. 1) CALL LOGW_CHAR_CHAR( &
!!!!              SHELL,'transition does not have tabulated tables.')
!!!!          IF (IERR .EQ. 2) CALL LOGW_CHAR_REAL( &
!!!!              'Lower effective quantum number is outside tables ',NLOW_EFF)
!!!!          IF (IERR .EQ. 3) CALL LOGW_CHAR_REAL( &
!!!!              'Upper effective quantum number is outside tables ',NUPP_EFF)
!!!!          IF (IERR .EQ. 0) THEN
!!!!             CALL LOGW_CHAR_REAL('Collisional parameters interpolated for ',LINE_L0(I))
!!!!             CALL LOGW_CHAR_REAL('Collisional cross section is ',SIGMA(I))
!!!!             CALL LOGW_CHAR_REAL('Temperature parameter is ',ALPHA(I))
!!!!          ENDIF
!!!!       ENDIF
!!!!       ! Case 3: parameters not provided by lines_database.dat
!!!!       ! and transition of the optical electron is not available
!!!!       IF ( (ALPHA(I) .EQ. 0D0 .OR. SIGMA(I).EQ.0D0) .AND. (FLAG .NE.0)) THEN
!!!!          CALL LOGW_CHAR_REAL('Collisional parameters not available for ',LINE_L0(I))
!!!!          CALL LOGW_CHAR('For this line we will use Van der Waals broadening')
!!!!       ENDIF
!!!!    ENDDO
!!!!    !
!!!!    READ(UNIT=1,FMT='(A)', IOSTAT=IERR) FIN
!!!!    IF (TRIM(FIN).NE.'END') THEN
!!!!       PRINT*,'File lines_database.dat could not be properly read. Check format. STOP'
!!!!       STOP
!!!!    ENDIF
!!!!    CLOSE(1)
!!!!    CALL LOGW_CHAR('File containing database of spectral lines (lines_database.dat) has been correctly read.')
!!!!    !-------------------------------------
!!!!  END SUBROUTINE READ_LINES_DATABASE
  !
  !------------------------------------------------
  !
  SUBROUTINE ELECCONF(ELC,SOUT,LOUT,JOUT)
    !
    IMPLICIT NONE
    !
    CHARACTER*7                     :: ELC
    CHARACTER*1                     :: L
    INTEGER                         :: IERR
    REAL(8)                         :: SOUT, LOUT, JOUT
    REAL(4)                         :: S, J
    !
    READ(ELC(1:3),'(F3.1)') S
    READ(ELC(4:4),'(A1)') L
    READ(ELC(5:7),'(F3.1)') J
    JOUT=DBLE(J)
    SOUT=(DBLE(S)-1D0)/2D0
    IERR=1
    !
    SELECT CASE (L)
    CASE('S')
       LOUT=0D0
       IERR=0
    CASE('P')
       LOUT=1D0
       IERR=0
    CASE('D')
       LOUT=2D0
       IERR=0
    CASE('F')
       LOUT=3D0
       IERR=0
    CASE('G')
       LOUT=4D0
       IERR=0
    CASE('H')
       LOUT=5D0
       IERR=0
    CASE('I')
       LOUT=6D0
       IERR=0
    END SELECT
    !
    IF (IERR.NE.0) THEN
       PRINT*,'Unknown electronic configuration '//TRIM(L)//' in lines_database.dat. STOP'
       STOP
    ENDIF
    !
  ENDSUBROUTINE ELECCONF
  !
  !------------------------------------------------
  !
  SUBROUTINE OPTICAL_ELEC_TRAN(SHELL,TRANSITION,FLAG)
    ! This routine determines the transition type for the optical electron
    ! that is, the last electron in the electronic configuration. In general
    ! it is different from the transition that gives raise to the spectral line
    !
    IMPLICIT NONE
    !
    CHARACTER*2,               INTENT(IN)     :: SHELL
    INTEGER,                   INTENT(OUT)    :: TRANSITION(2)
    INTEGER,                   INTENT(INOUT)  :: FLAG
    CHARACTER*1                               :: ORBITAL
    INTEGER                                   :: I
    !
    DO I=1,2
       ORBITAL=SHELL(I:I)
       FLAG = 1
       SELECT CASE (ORBITAL)
          CASE('s')
             TRANSITION(I)=0
             FLAG=0
          CASE('S')
             TRANSITION(I)=0
             FLAG=0
          CASE('p')
             TRANSITION(I)=1
             FLAG=0
          CASE('P')
             TRANSITION(I)=1
             FLAG=0
          CASE('d')
             TRANSITION(I)=2
             FLAG=0
          CASE('D')
             TRANSITION(I)=2
             FLAG=0
          CASE('f')
             TRANSITION(I)=3
             FLAG=0
          CASE('F')
             TRANSITION(I)=3
             FLAG=0
       END SELECT
    ENDDO
    ! If flag = 1 we will not be able to determine the collisional damping with ABO theory
  END SUBROUTINE OPTICAL_ELEC_TRAN
  !
  !================================================
  !
END MODULE LINES_DATABASE
!
