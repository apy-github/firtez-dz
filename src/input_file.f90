!
MODULE INPUT_FILE
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP
  USE user_mpi, ONLY: mpi__myrank
  !
  IMPLICIT NONE
  !
  CHARACTER*800    :: INPUT_READ_LINES(1000)
  INTEGER          :: INPUT_READ_NUMBER
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: READ_INPUT_FILE, CLOSE_READ_INPUT_FILE, READ_MINPUT_FILE
  !
  !PRIVATE :: READ_BOX
  !PRIVATE :: READ_LINES
!  PRIVATE :: SPLIT_READ_LINE_SPECTRAL
!  PRIVATE :: SPLIT_READ_LINE_BLENDS_NUMBER
!  PRIVATE :: SPLIT_READ_LINE_BLENDS
!  PRIVATE :: SET_VAR_INVERSION
!  PRIVATE :: SET_STK_INVERSION
!  PRIVATE :: SPLIT_READ_LINE_WEIGHTS
!  PRIVATE :: SPLIT_READ_LINE_PSF
!  PRIVATE :: SET_CONTINUUM_NORMALIZATION
!  PRIVATE :: SET_COUPLED_INVERSION
!  PRIVATE :: SPLIT_LSF_LINE
!  PRIVATE :: SET_LSF
!  PRIVATE :: SET_MISC
!  PRIVATE :: UPDATE_FREE_VAR
!  PRIVATE :: FREE_VAR_SPACE
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! read_input_file
  ! read_box
  ! read_lines
  ! split_read_line_spectral
  ! set_var_inversion
  ! set_stk_inversion
  ! split_read_line_weights
  ! split_read_line_psf
  ! set_continuum_normalization
  ! set_coupled_inversion
  ! split_lsf_line
  ! set_lsf
  ! set_misc
  ! close_read_input_file
  ! update_free_var
  ! free_var_space
  !
  !________________________________________________
  !
  SUBROUTINE READ_MINPUT_FILE()
    !
    USE CODE_MODES
    USE user_mpi, ONLY: mpi__ierror, MPI_CHARACTER, MPI_COMM_WORLD, MPI_INTEGER
    !
    CHARACTER*800                              :: LINE
    INTEGER                                    :: IERR, I, J, LINES_READ
    INTEGER                                    :: NUM_LINES_OBS
    !
    INTEGER                                    :: CNT
    !
    ! Does input file exist ?
    !
    INPUT_READ_LINES(:)=' '
    IF (mpi__myrank.eq.0) THEN

      CALL GETARG(1,INPUTFILE)
      OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR &
          , FORM='FORMATTED')
      IF (IERR.NE.0) THEN
        IF (mpi__myrank.EQ.0) THEN
          PRINT*,'The file containing information about the inversion box (spectral lines, sizes, etc.)'//TRIM(INPUTFILE)
          PRINT*,'could not be fonud in the directory of the source code. STOP'
          CLOSE(UNIT=1)
          STOP
        ENDIF
      ENDIF
      !
      ! Read the whole input file:
      !
      INPUT_READ_NUMBER=0
      DO WHILE (INPUT_READ_NUMBER.LT.1000)
        READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
        IF (IERR.LT.0) EXIT
        LINE=TRIM(LINE)
        ! Avoid empty lines:
        IF (LEN_TRIM(LINE).EQ.0) CYCLE
        ! Avoid comments:
        IF (LINE(1:1).EQ.'#') CYCLE
        ! Remove comments from valid lines:
        CNT=SCAN(LINE,'#')
        IF (CNT.NE.0) LINE=LINE(1:CNT-1)
        !
        ! Check if we are reaching "END:"
        ! It is actually not needed anymore, but if used
        ! ... we do not read anything else:
        IF (TRIM(LINE).EQ.'END:') EXIT
        ! Now, we store it
        INPUT_READ_NUMBER=INPUT_READ_NUMBER+1
        INPUT_READ_LINES(INPUT_READ_NUMBER)=LINE
        !
      ENDDO
      !
      CLOSE(1)
      !PRINT*, mpi__myrank, SIZEOF(INPUT_READ_LINES)
      CALL MPI_BCAST(INPUT_READ_LINES, SIZEOF(INPUT_READ_LINES), MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi__ierror)
      CALL MPI_BCAST(INPUT_READ_NUMBER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__ierror)
      !
    ELSE
      !
      !PRINT*, mpi__myrank, SIZEOF(INPUT_READ_LINES)
      CALL MPI_BCAST(INPUT_READ_LINES, SIZEOF(INPUT_READ_LINES), MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi__ierror)
      CALL MPI_BCAST(INPUT_READ_NUMBER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__ierror)
      !
    ENDIF
    !
    ! Once read, we first check paths:
    !
    !-----------------------------
    DATAPATH='./'
    OUTPPATH='./'
    MODLPATH='./'
    LINEPATH='./'
    !-----------------------------
    CNT = 0
    !
    DO WHILE (CNT.LE.INPUT_READ_NUMBER)
      !
      CNT = CNT + 1
      !
      IF (TRIM(INPUT_READ_LINES(CNT)).EQ.'MISC SETUP:') THEN
        CALL NSET_MISC(CNT)
        EXIT
      ENDIF
    ENDDO
    !
  END SUBROUTINE READ_MINPUT_FILE
  !
  !________________________________________________
  !
  SUBROUTINE NSET_MISC(START_LINE)
    !
    USE user_mpi
    USE INVERT_PARAM, ONLY: ASSIST_T, ASSIST_P, ASSIST_B, ASSIST_V &
        , RSVDTOL, RSIGMAP, ISIGMA
    USE FORWARD_PARAM, ONLY: FULL_STOKES
    USE CODE_MODES, ONLY: INPUTFILE, DATAPATH, OUTPPATH, LINEPATH &
        , MODLPATH, MVERBOSE
    !
    INTEGER, INTENT(IN)              :: START_LINE
    CHARACTER*800                    :: LINE
    INTEGER                          :: J
    LOGICAL                          :: CONDIT
    REAL(SP)                         :: DUMSP
    INTEGER                          :: DUMIP
    LOGICAL                          :: DUMSVDL, DUMSIGL
    INTEGER                          :: CNT
    LOGICAL                          :: EXISTS
    !
    !
    CONDIT=.TRUE.
    CNT=START_LINE+1
    !
    DO WHILE (CNT.LE.INPUT_READ_NUMBER)
      !
      DUMSVDL=.FALSE.
      DUMSIGL=.FALSE.
      !
      LINE=INPUT_READ_LINES(CNT)
      !
      ! Check if the section has finished:
      !
      DO J=1,800
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          EXIT
        ENDIF
      ENDDO
      !
      CNT=CNT+1
      !
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'TEMP', 0, ASSIST_T, DUMSP, DUMIP)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PGAS', 0, ASSIST_P, DUMSP, DUMIP)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BVEC', 0, ASSIST_B, DUMSP, DUMIP)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'VLOS', 0, ASSIST_V, DUMSP, DUMIP)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'FULL_STOKES', 0, FULL_STOKES, DUMSP, DUMIP)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'SVDTOL', 0, DUMSVDL, DUMSP, DUMIP)
      IF (DUMSVDL.EQV..TRUE.) RSVDTOL=DBLE(DUMSP)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'NOISE', 0, DUMSIGL, DUMSP, DUMIP)
      IF (DUMSIGL.EQV..TRUE.) ISIGMA=DBLE(DUMSP)
      CALL SPLIT_READ_LINE_STRING(LINE, 'DATAPATH', DATAPATH)
      CALL SPLIT_READ_LINE_STRING(LINE, 'MODPATH', MODLPATH)
      CALL SPLIT_READ_LINE_STRING(LINE, 'OUTPATH', OUTPPATH)
      CALL SPLIT_READ_LINE_STRING(LINE, 'SLDPATH', LINEPATH)
      !
    ENDDO
    !
    ! Can I write into outpath?
    !
    IF (mpi__myrank.EQ.0) THEN
      IF (TRIM(OUTPPATH).NE.'./') THEN
        INQUIRE(FILE=TRIM(OUTPPATH)//"/./", EXIST=EXISTS)
        IF (EXISTS.EQV..FALSE.) THEN
          PRINT*, ""
          PRINT*, "Error!"
          PRINT*, "I cannot write to output directory: ", TRIM(OUTPPATH)
          PRINT*, ""
          CALL MPI_ABORT(MPI_COMM_WORLD)
        ENDIF ! Writing error
      ENDIF ! Path != default
    ENDIF ! Master
    !
  END SUBROUTINE NSET_MISC
  !
  !________________________________________________
  !
  SUBROUTINE NREAD_BOX(START_LINE)
    !
    USE CODE_MODES, ONLY: INPUTFILE
    USE GRID_PARAM, ONLY: NX, DX, NY, DY, NZ, DZ
    !
    INTEGER, INTENT(INOUT)     :: START_LINE
    INTEGER                    :: I, IERR, IOS
    !
    DO I=1,3
      !
      START_LINE=START_LINE+1
      !
      SELECT CASE (I)
        CASE(1)
          !
          ! Read new line: NX
          !
          CALL SPLIT_READ_BOX_LINE(TRIM(INPUT_READ_LINES(START_LINE)), NX, DX)
        CASE(2)
          !
          ! Read new line: NY
          !
          CALL SPLIT_READ_BOX_LINE(TRIM(INPUT_READ_LINES(START_LINE)), NY, DY)
        CASE(3)
          !
          ! Read new line: NZ
          !
          CALL SPLIT_READ_BOX_LINE(TRIM(INPUT_READ_LINES(START_LINE)), NZ, DZ)
      ENDSELECT
   ENDDO
    !
  END SUBROUTINE NREAD_BOX
  !________________________________________________
  !
  SUBROUTINE NREAD_LINES(START_LINE)
    !
    USE CODE_MODES, ONLY: INPUTFILE
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_1D_DP &
        , ALLOCATE_2D_IP, ALLOCATE_2D_DP
    USE FORWARD_PARAM
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    INTEGER,   INTENT(INOUT)              :: START_LINE
    CHARACTER*800                         :: LINE
    INTEGER                               :: IERR, I, K, J
    !
    INTEGER                               :: NUMC
    INTEGER                               :: BLENDSNUM
    INTEGER                               :: NUM_LINES_OBS
    LOGICAL                               :: END_FIELD
    CHARACTER*800                         :: LINES(100)
    !
    ! Determine number of lines inside and allocate arrays
    !
    END_FIELD=.FALSE.
    IERR=0
    NUM_LINES_OBS=0
    
    DO WHILE (START_LINE.LE.INPUT_READ_NUMBER)
       !
       START_LINE=START_LINE+1
       LINE=INPUT_READ_LINES(START_LINE)
!PRINT*, TRIM(LINE)
       !
       DO J=1,LEN(LINE)
         IF (':'.EQ.LINE(J:J)) THEN
           END_FIELD=.TRUE.
           START_LINE=START_LINE-1
           EXIT
         ENDIF
       ENDDO
       IF (END_FIELD.EQV..TRUE.) EXIT
       NUM_LINES_OBS=NUM_LINES_OBS+1
       LINES(NUM_LINES_OBS)=LINE
    ENDDO
    CLOSE(1)

    !
    NUML=NUM_LINES_OBS
    CALL ALLOCATE_1D_IP(IND_LINE,NUML,'ALLOCATE ind_line IN initial_duties.f90')
    CALL ALLOCATE_1D_DP(WAVE_INI,NUML,'ALLOCATE wave_ini IN initial_duties.f90')
    CALL ALLOCATE_1D_DP(DELTA_LAMBDA,NUML,'ALLOCATE delta_lambda IN initial_duties.f90')
    CALL ALLOCATE_1D_IP(NUMWAVE,NUML,'ALLOCATE numwave IN initial_duties.f90')
    !
    !
    ! Start a variable for the maximum number of wavelengths:
    DO I=1,NUML
       LINE=LINES(I)
       NUMC=NCOLUMNS(LINE)

       CALL SPLIT_READ_LINE_SPECTRAL(LINE &
           , IND_LINE(I), NUMWAVE(I), WAVE_INI(I), DELTA_LAMBDA(I))
       DO K=1,NUML_DATABASE
         IF (LINE_NUM(K) .EQ. IND_LINE(I)) THEN
           IND_LINE(I)=K
           EXIT
         ENDIF
       ENDDO
       CALL SPLIT_READ_LINE_BLENDS_NUMBER(LINE, BLENDSNUM)
       IF (BLENDSNUM.GT.BLENDSMAX) BLENDSMAX=BLENDSNUM
    ENDDO
    !
    ! Now, repeat but storing spectral line info:
    ! Read until last line
    !
    ! Now, we already know the maximum number of blends to be considered:
    CALL ALLOCATE_2D_IP(BLENDSID,NUML,BLENDSMAX+1,'BLENDSID IN READ_LINES.')
    CALL ALLOCATE_2D_DP(BLENDSDIFF,NUML,BLENDSMAX+1,'BLENDSDIFF IN READ_LINES.')
    !
    !
    DO I=1,NUML
    ! First column, reference line id
      BLENDSID(I,1)=IND_LINE(I)
      !
      LINE=LINES(I)
      !LINE=INPUT_READ_LINES(START_LINE)
      !
      CALL SPLIT_READ_LINE_BLENDS(LINE, BLENDSMAX, BLENDSID(I,2:))
      DO J=1,BLENDSMAX
        IF (BLENDSID(I,J+1).GT.0) THEN
          DO K=1,NUML_DATABASE
            IF (LINE_NUM(K) .EQ. BLENDSID(I,J+1)) THEN
              BLENDSID(I,J+1)=K
              EXIT
            ENDIF
          ENDDO
          IF (K.GT.NUML_DATABASE) THEN
            PRINT*, ''
            PRINT*, ' ERROR!'
            PRINT*, ''
            PRINT*, '   BLEND ID', BLENDSID(I,J+1), ' IS NOT IN lines_database.dat'
            PRINT*, ''
            STOP
          ENDIF
! DIFFERENCE OF WAVELENGTH IN ANGSTROM:
          BLENDSDIFF(I,J+1)=LINE_L0(BLENDSID(I,J+1))-LINE_L0(IND_LINE(I))
! IN mA:
          BLENDSDIFF(I,J+1)=1.0D3*BLENDSDIFF(I,J+1)
        ENDIF
      ENDDO
      !print*, '??', BLENDSID(I,1), ' ; ', TRIM(LINE)
    ENDDO
!PRINT*, "NREAD_LINES:", IND_LINE,NUML
    !
    NUMW=SUM(NUMWAVE)
!DO I=1,NUML
!  PRINT*, BLENDSID(I,:)
!ENDDO
!STOP
    !
    !LINES_READ=LINES_READ+NUM_LINES_OBS
    !
  END SUBROUTINE NREAD_LINES

  !________________________________________________
  !
  SUBROUTINE READ_INPUT_FILE()
    !
    USE CODE_MODES
    USE GRID_PARAM, ONLY: NZ
    USE INVERT_PARAM
! Obsolete nlte dep. coef.:    USE FORWARD_PARAM, ONLY: LINE_L0, IND_LINE, LINE_POS, IND_LINE, NUML &
! Obsolete nlte dep. coef.:        , NUML_DATABASE, FULL_STOKES, HYDRO_TOP
    USE FORWARD_PARAM, ONLY: LINE_L0, IND_LINE, LINE_POS, IND_LINE, NUML &
           , NUML_DATABASE, FULL_STOKES, HYDRO_TOP, POPL, POPU 
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_2D_DP
    !USE LOG, ONLY: LOGW_CHAR_REAL
    !
    CHARACTER*800                              :: LINE
    CHARACTER(LEN=1)                           :: CHARMODE
    CHARACTER(LEN=1)                           :: CHARREAD
    INTEGER                                    :: IERR, I, J, LINES_READ
    INTEGER,    DIMENSION(NUML_DATABASE)       :: DIF
    !
    LOGICAL                                    :: FILEPROFILE, FILEMODEL
    LOGICAL                                    :: WORKINGMODE
    LOGICAL                                    :: TMPDERS
    INTEGER                                    :: CNT
    !
    !-----------------------------
    !-----------------------------
    BOX = .FALSE.
    LINES = .FALSE.
    SIR = .FALSE.
    SPINOR = .FALSE.
    MURAM = .FALSE.
    WORKINGMODE = .FALSE.
    HYDROSTATIC = .FALSE.
    HYDROSTATICDER = .FALSE.
    TMPDERS=.TRUE.
    !-----------------------------
    !-----------------------------
    FILEPROFILE = .FALSE.
    FILEMODEL = .FALSE.
    MSYNTHESIS = .FALSE.
    MINVERSION = .FALSE.
    MGETTAU = .FALSE.
    MGETHEQ = .FALSE.
    MGETERRORS = .FALSE.
    MRESPFUNCT = .FALSE.
    SAVE_RFS = .FALSE.
    MLSF=.FALSE.
    MVERBOSE=0
    MTAULIN=.FALSE.
    !
    INVERSIONTYPE=0  ! By default, SIR/NICOLE inversion nodes
    !
    VREGULARIZATION=.FALSE.
    PEN_TYP(:)=0
    PEN_ALP(:)=0.0D0
    !-----------------------------
    INV_ATMPAR(:) = .FALSE.
    !-----------------------------
    NSLB_MAX(:)=-1
    !-----------------------------
    INV_STK(:) = .FALSE.
    !
    AUTOWEIGHT=.FALSE.
    !-----------------------------
    ISIGMA=4.E-3
    !-----------------------------
    HSRA_NORMALIZATION = .TRUE.
    !-----------------------------
    CHARMODE='A'
    CHARREAD='A'
    !-----------------------------
    COUPLED=.FALSE.
    !-----------------------------
    TWODSPATIALSPD=.FALSE.
    !TWODSPATIALSPD=.TRUE.
    !-----------------------------
    ASSIST_T=.FALSE.
    ASSIST_P=.FALSE.
    ASSIST_B=.FALSE.
    ASSIST_V=.FALSE.
    !-----------------------------
    FULL_STOKES=.FALSE.
    !-----------------------------
    CNT = 0
    !
    LINE=''
    !
    ! First read BOX and if it is not present, take NX, NY, NZ from model header
    DO WHILE (CNT.LE.INPUT_READ_NUMBER)
       !
       CNT = CNT + 1
       LINE=INPUT_READ_LINES(CNT)
       !
       SELECT CASE (TRIM(LINE))
       CASE('BOX:')
          BOX = .TRUE.
          CALL NREAD_BOX(CNT)
       CASE('FILEPROFILE:')
          FILEPROFILE = .TRUE.
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          NAMEPROFILE=TRIM(LINE)
       CASE('FILEMODEL:')
          FILEMODEL = .TRUE.
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          NAMEMODEL=TRIM(LINE)
       CASE('END:')
          EXIT
       ENDSELECT
    ENDDO
    !
    CALL BOX_AVAILABLE(BOX)
    !
    CNT = 0
    !
    ! Now, we can proceed as always:
    LINE=''
    DO WHILE (CNT.LE.INPUT_READ_NUMBER)
       !
       CNT = CNT + 1
       LINE=INPUT_READ_LINES(CNT)
       !
       SELECT CASE (TRIM(LINE))
       CASE('LINES:')
          LINES = .TRUE.
          CALL NREAD_LINES(CNT)
          ! NZ was already stablished & NUML is determined ny NREAD_LINES
          ! So, now we can allocate the arrays for departure coefficients
          ! We need them even if we will not do NLTE
          CALL ALLOCATE_2D_DP(POPL,NUML,NZ,'POPL')
          CALL ALLOCATE_2D_DP(POPU,NUML,NZ,'POPU')
          ! Ensure LTE to begin with
          POPL(:,:)=1D0
          POPU(:,:)=1D0
       CASE('MODE:')
          WORKINGMODE = .TRUE.
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          CHARMODE=TRIM(LINE)
       CASE('HYDROSTATIC:')
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          CHARREAD=TRIM(LINE)
          IF ((CHARREAD.EQ.'y').OR.(CHARREAD.EQ.'Y')) THEN
            HYDROSTATIC = .TRUE.
            HYDRO_TOP=.TRUE.
            !HYDRO_TOP=.FALSE.
          ENDIF
       CASE('HYDROSTATIC DERIVATIVES:')
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          CHARREAD=TRIM(LINE)
          IF ((CHARREAD.EQ.'n').OR.(CHARREAD.EQ.'N')) THEN
            TMPDERS = .FALSE.
          ELSE
            TMPDERS = .TRUE.
          ENDIF
       CASE('STOKES SETUP:')
          CALL NSET_STK_INVERSION(CNT)
       CASE('INVERSION SETUP:')
          CALL NSET_VAR_INVERSION(CNT)
       CASE('CONTINUUM NORMALIZATION:')
          CALL NSET_CONTINUUM_NORMALIZATION(CNT)
       CASE('COUPLED INVERSION:')
          CALL NSET_COUPLED_INVERSION(CNT)
       CASE('VERTICAL REGULARIZATION TYPE:')
          CALL NSET_VREGULARIZATION_TYPE(CNT)
       CASE('VERTICAL REGULARIZATION NORM:')
          CALL NSET_VREGULARIZATION_NORM(CNT)
       CASE('WRITE RESPONSE FUNCTION:')
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          CHARREAD=TRIM(LINE)
          IF ((CHARREAD.EQ.'y').OR.(CHARREAD.EQ.'Y')) THEN
            SAVE_RFS = .TRUE.
            MRESPFUNCT = .TRUE.
          ENDIF
       CASE('INVERSION NODES:')
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          CHARREAD=TRIM(LINE)
          IF ((CHARREAD.EQ.'e').OR.(CHARREAD.EQ.'E')) THEN
            INVERSIONTYPE=1
          ENDIF
       CASE('LINE OPTICAL DEPTH:')
          CNT = CNT + 1
          LINE=INPUT_READ_LINES(CNT)
          CHARREAD=TRIM(LINE)
          IF ((CHARREAD.EQ.'y').OR.(CHARREAD.EQ.'Y')) THEN
            MTAULIN = .TRUE.
          ENDIF
       CASE('LINE SPREAD FUNCTION:')
          CALL NSET_LSF(CNT)
       CASE('MISC SETUP:')
          CALL NSET_MISC(CNT)
       CASE('NLTE DEPARTURE COEFFICIENTS:')
         CALL NSET_NLTE(CNT)
       CASE('END:')
          EXIT
       ENDSELECT
    ENDDO
    !
    ! Check some stuff:
    !
    !
    IF (HYDROSTATIC.EQV..TRUE.) THEN
      HYDROSTATICDER=.TRUE.
      IF (TMPDERS.EQV..FALSE.) HYDROSTATICDER=.FALSE.
    ENDIF
    !
    IF (FILEMODEL.EQV..FALSE.) THEN
IF (mpi__myrank.EQ.0) THEN
       PRINT*, ''
       PRINT*, 'INPUT FILE MUST! CONTAIN NO EMPTY...'
       PRINT*, '  "FILEMODEL:" FIELD.'
       PRINT*, ''
       STOP
ENDIF
    ENDIF
    !
    IF (WORKINGMODE.EQV..FALSE.) THEN
IF (mpi__myrank.EQ.0) THEN
       PRINT*, ''
       PRINT*, 'INPUT FILE MUST! CONTAIN NO EMPTY...'
       PRINT*, '                       "MODE:" FIELD'
       PRINT*, ''
       STOP
ENDIF
    ENDIF
    !
    IF ((CHARMODE.EQ.'s').OR.(CHARMODE.EQ.'S')) THEN
      MSYNTHESIS = .TRUE.
      !
      ! Check some definition of some parameters:
      IF (NSTKINV.EQ.0) THEN
        NSTKINV=4 
        INV_STK(:)=.TRUE.
      ENDIF ! definition parameters.
      !
      ! Check some variables for RF:
      IF (MRESPFUNCT.EQV..TRUE.) THEN
        DO I=1,SIZE(INV_ATMPAR)
          IF (INV_ATMPAR(I).EQV..TRUE.) NSLB_MAX(I)=NZ
        ENDDO
        CALL FREE_VAR_SPACE(NZ)
      ENDIF ! RF variables.
    ELSE IF ((CHARMODE.EQ.'i').OR.(CHARMODE.EQ.'I')) THEN
      MRESPFUNCT = .TRUE.
      MINVERSION = .TRUE.
      MGETERRORS = .TRUE.
      ! IN INVERSION MODE, WE HAVE TO CHECK WHAT STOKES HAVE TO BE INVERTED:
      IF (ANY(INV_STK).EQV..FALSE.) THEN
        WRITE(*,*) ''
        WRITE(*,*) '"STOKES SETUP:" THIS FIELD IS MANDATORY IN INVERSION MODE...'
        WRITE(*,*) '                                                     STKI...'
        WRITE(*,*) '                                                     STKQ...'
        WRITE(*,*) '                                                     STKU...'
        WRITE(*,*) '                                                     STKV.'
        WRITE(*,*) ''
        STOP
      ENDIF
      CALL FREE_VAR_SPACE(NZ)
    ELSE IF ((CHARMODE.EQ.'t').OR.(CHARMODE.EQ.'T')) THEN
      MGETTAU= .TRUE.
    ELSE IF ((CHARMODE.EQ.'h').OR.(CHARMODE.EQ.'H')) THEN
      MGETHEQ= .TRUE.
    ELSE
       WRITE(*,*) ''
       WRITE(*,*) '"MODE:" FIELD IN INPUT FILE MUST! CONTAIN...'
       WRITE(*,*) '                      SYNTHESIS/synthesis...'
       WRITE(*,*) '                      INVERSION/inversion...'
       WRITE(*,*) '                                    tau/TAU.'
       WRITE(*,*) ''
       STOP
    ENDIF
    !
!><    ! CHECK IF INVERSION PARAMETERS ARE CORRECT:
!><    !   WE DO ALLOW CARTESIAN OR SPHERICAL INVERSION PARAMETERS, NO MIXING IS ALLOWED
!><    IF (&
!><        ( (INV_BXNN.EQV..TRUE.) .OR. (INV_BYNN.EQV..TRUE.) .OR. (INV_BZNN.EQV..TRUE.) ) &
!><       .AND.&
!><        ( (INV_BNNN.EQV..TRUE.) .OR. (INV_GAMN.EQV..TRUE.) .OR. (INV_PHIN.EQV..TRUE.) )&
!><       ) THEN
!><      WRITE(*,*) ''
!><      WRITE(*,*) 'YOU CAN NOT MIX MAGNETIC FIELDS IN SPHERICAL AND CARTESIAN...'
!><      WRITE(*,*) '                                            REFERENCE FRAMES.'
!><      WRITE(*,*) ''
!><      STOP
!><    ENDIF
    !
    IF ((MINVERSION.EQV..TRUE.).AND.(&
        (ANY(INV_ATMPAR).EQV..FALSE.))) THEN 
       WRITE(*,*) ''
       WRITE(*,*) 'IN INVERSION MODE, "INVERSION SETUP: " FIELD IN INPUT FILE...'
       WRITE(*,*) ' MUST! CONTAIN...'
       WRITE(*,*) '                                       TEM...'
       WRITE(*,*) '                                        P0...'
       WRITE(*,*) '                                      PGAS...'
       WRITE(*,*) '                                       RHO...'
       WRITE(*,*) '                                        BX...'
       WRITE(*,*) '                                        BY...'
       WRITE(*,*) '                                        BZ...'
       WRITE(*,*) '                                        VLOS.'
       STOP
    ENDIF
    !
    IF ( (MINVERSION .EQV. .TRUE.).OR.(MSYNTHESIS .EQV. .TRUE.) ) THEN
      ! Check that there are spectral lines in the input file
      IF (LINES.EQV..FALSE.) THEN
IF (mpi__myrank.EQ.0) THEN
         PRINT*,'Missing spectral lines (LINES:) in the input file.'
         PRINT*,'Check '//TRIM(INPUTFILE)
         PRINT*,'STOP'
         STOP
ENDIF
      ENDIF
      ! Check that indexes of spectral lines exist in database
      DO I=1,NUML
         DIF(:)=1E3
         DO J=1,NUML_DATABASE
            !DIF(J)=ABS(LINE_NUM(J)-IND_LINE(I))
            DIF(J)=ABS(LINE_POS(J)-IND_LINE(I))
         ENDDO
         IF (MINVAL(DIF).NE.0) THEN
IF (mpi__myrank.EQ.0) THEN
PRINT*, IND_LINE(I)
            PRINT*,'Spectral line index in '//TRIM(INPUTFILE)//' cannot'
            PRINT*,'be found in databse (lines_database.dat). STOP'
            STOP
ENDIF
         ENDIF
      ENDDO
      IF (FILEPROFILE.EQV..FALSE.) THEN
IF (mpi__myrank.EQ.0) THEN
         PRINT*, ''
         PRINT*, 'INPUT FILE MUST! CONTAIN NO EMPTY...'
         PRINT*, '..."FILEPROFILE:" in Synthesis/inversion mode.'
         PRINT*, ''
         STOP
ENDIF
      ENDIF
    ELSE
      CALL SETUP_FAKE_LINES()
    ENDIF
  !
  END SUBROUTINE READ_INPUT_FILE
  !
  !________________________________________________
  !
  SUBROUTINE BOX_AVAILABLE(LOGIC)
    !
    USE CODE_MODES, ONLY: NAMEMODEL, MODLPATH
    USE MISC, ONLY: READ_SP1D
    USE GRID_PARAM, ONLY: NX, NY, NZ
    !
    LOGICAL, INTENT(IN) :: LOGIC
    !
    REAL(SP), DIMENSION(9) :: TOREAD
    !
    IF (LOGIC.EQV..FALSE.) THEN
      ! If we do not supply BOX info, then we have to get nz, ny, nz from...
      ! ...model hearder:
      CALL READ_SP1D(TRIM(MODLPATH)//TRIM(NAMEMODEL), SHAPE(TOREAD), TOREAD)
      NX=NINT(TOREAD(7))
      NY=NINT(TOREAD(8))
      NZ=NINT(TOREAD(9))
      !
    ENDIF
    !
  END SUBROUTINE BOX_AVAILABLE
  !
  !________________________________________________
  !
  SUBROUTINE SETUP_FAKE_LINES ()
    !
    USE FORWARD_PARAM
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_1D_DP &
        , ALLOCATE_2D_IP, ALLOCATE_2D_DP
    !
    ! Careful, if we are neither in synthesis nor in inversion mode...
    ! ... we initialize the various line variables required:
    !
    NUML=1
    CALL ALLOCATE_1D_IP(IND_LINE,NUML,'ALLOCATE ind_line')
    CALL ALLOCATE_1D_DP(WAVE_INI,NUML,'ALLOCATE wave_ini')
    CALL ALLOCATE_1D_DP(DELTA_LAMBDA,NUML,'ALLOCATE delta_lambda')
    CALL ALLOCATE_1D_IP(NUMWAVE,NUML,'ALLOCATE numwave')
    !
    IND_LINE(1)=1
    NUMWAVE(1)=1
    WAVE_INI(1)=0.D0
    DELTA_LAMBDA(1)=0.0D0
    !
    BLENDSMAX=0
    !
    ! Now, we already know the maximum number of blends to be considered:
    CALL ALLOCATE_2D_IP(BLENDSID,NUML,BLENDSMAX+1,'BLENDSID.')
    CALL ALLOCATE_2D_DP(BLENDSDIFF,NUML,BLENDSMAX+1,'BLENDSDIFF.')
    !
    ! First column, reference line id
    BLENDSID(1,1)=IND_LINE(1)
    !APY
    NUMW=SUM(NUMWAVE)
    !END APY
    !
  END SUBROUTINE SETUP_FAKE_LINES
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_BOX_LINE(LINE, NN, DD)
    !
    CHARACTER(*), INTENT(IN) :: LINE
    INTEGER, INTENT(INOUT) :: NN
    REAL(DP), INTENT(INOUT) :: DD
    !
    REAL(SP)                    :: DVALUE
    INTEGER                     :: IVALUE
    !
    INTEGER                     :: I, NCOL, START, PRE, ISREAL, J
    LOGICAL                     :: ILOGIC
    !
    NCOL=1
    START=1
    PRE=0
    !
    DO I=1,LEN(LINE)
! If we have an empty space, then we might need to interpret the thing:
      IF ((LINE(I:I).EQ.' ') .OR.((I.EQ.LEN(LINE)).AND.(LINE(I:I).NE.' ')))THEN
! We need to do something if there are previous non empty characters (pre=1)
        IF ((PRE.EQ.1).OR.((I.EQ.LEN(LINE)).AND.(LINE(I:I).NE.' '))) THEN
          !READ(LINE(START:I), '(A40)') STERM
          ! ASSIGN TERM TO EACH COLUMN:
          SELECT CASE(NCOL)
            CASE(1)
! Standard value for number of grid points in spatia. direction:
              READ(LINE(START:I), '(i14)') NN
            CASE(2)
              ! First, we need to check if the number ...
              ! ...is given as integer or otherwise:
              ISREAL=0
              DO J=1,LEN(LINE(START:I))
                IF(LINE(START+J:START+J).EQ.'.') ISREAL=1
              ENDDO
              !
              IF (ISREAL.NE.1) THEN
                READ(LINE(START:I), '(i14)') IVALUE
                DD=DBLE(IVALUE)
              ELSE
                READ(LINE(START:I), '(e14.5)') DVALUE
                DD=DBLE(DVALUE)
              ENDIF
          END SELECT
          !End assignment
        ENDIF
        PRE=0
        START=I
        NCOL=NCOL+1
        CYCLE
      ENDIF
! If it is not an empty space and the previous character was an empty...
! ...character, then add a column:
      IF (PRE.EQ.0) THEN
        PRE=1
      ENDIF

    ENDDO
    !
  END SUBROUTINE SPLIT_READ_BOX_LINE
  !
  !________________________________________________
  !
!!!!!!!!!!!  SUBROUTINE READ_BOX(LINES_READ)
!!!!!!!!!!!    !
!!!!!!!!!!!    USE CODE_MODES, ONLY: INPUTFILE
!!!!!!!!!!!    USE GRID_PARAM, ONLY: NX, DX, NY, DY, NZ, DZ
!!!!!!!!!!!    !
!!!!!!!!!!!    INTEGER, INTENT(INOUT)           :: LINES_READ
!!!!!!!!!!!    INTEGER,           PARAMETER     :: NARGS = 2
!!!!!!!!!!!    INTEGER                          :: I, IERR, IOS
!!!!!!!!!!!    INTEGER,       DIMENSION(NARGS)  :: IND_INI,IND_END
!!!!!!!!!!!    CHARACTER*800                    :: LINE
!!!!!!!!!!!    ! Open
!!!!!!!!!!!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!!!!!!!!!!!    ! Read until last line
!!!!!!!!!!!    DO I=1,LINES_READ
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    ! Read new line: NX
!!!!!!!!!!!    READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    CALL SPLIT_READ_BOX_LINE(LINE, NX, DX)
!!!!!!!!!!!    LINES_READ = LINES_READ + 1
!!!!!!!!!!!    IF (IERR.NE.0) THEN
!!!!!!!!!!!IF (mpi__myrank.EQ.0) THEN
!!!!!!!!!!!       PRINT*,'Error while reading '//TRIM(INPUTFILE)//'. STOP.'
!!!!!!!!!!!       PRINT*,'Error in line #',LINES_READ
!!!!!!!!!!!       STOP
!!!!!!!!!!!ENDIF
!!!!!!!!!!!    ENDIF
!!!!!!!!!!!    ! Read new line: NY
!!!!!!!!!!!    READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    CALL SPLIT_READ_BOX_LINE(LINE, NY, DY)
!!!!!!!!!!!    LINES_READ = LINES_READ + 1
!!!!!!!!!!!    IF (IERR.NE.0) THEN
!!!!!!!!!!!IF (mpi__myrank.EQ.0) THEN
!!!!!!!!!!!       PRINT*,'Error while reading '//TRIM(INPUTFILE)//'. STOP.'
!!!!!!!!!!!       PRINT*,'Error in line #',LINES_READ
!!!!!!!!!!!       STOP
!!!!!!!!!!!ENDIF
!!!!!!!!!!!    ENDIF
!!!!!!!!!!!    ! Read new line: NZ
!!!!!!!!!!!    READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    CALL SPLIT_READ_BOX_LINE(LINE, NZ, DZ)
!!!!!!!!!!!    LINES_READ = LINES_READ + 1
!!!!!!!!!!!    IF (IERR.NE.0) THEN
!!!!!!!!!!!IF (mpi__myrank.EQ.0) THEN
!!!!!!!!!!!       PRINT*,'Error while reading '//TRIM(INPUTFILE)//'. STOP.'
!!!!!!!!!!!       PRINT*,'Error in line #',LINES_READ
!!!!!!!!!!!       STOP
!!!!!!!!!!!ENDIF
!!!!!!!!!!!    ENDIF
!!!!!!!!!!!    !
!!!!!!!!!!!    CLOSE(1)
!!!!!!!!!!!    !
!!!!!!!!!!!  END SUBROUTINE READ_BOX
  !
  !________________________________________________
  !
!!!!!!!!!!!  SUBROUTINE READ_LINES(LINES_READ, NUM_LINES_OBS)
!!!!!!!!!!!    !
!!!!!!!!!!!    USE CODE_MODES, ONLY: INPUTFILE
!!!!!!!!!!!    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_1D_DP &
!!!!!!!!!!!        , ALLOCATE_2D_IP, ALLOCATE_2D_DP
!!!!!!!!!!!    USE FORWARD_PARAM
!!!!!!!!!!!    USE LINES_DATABASE, ONLY: NCOLUMNS
!!!!!!!!!!!    !
!!!!!!!!!!!    INTEGER,   INTENT(INOUT)              :: LINES_READ
!!!!!!!!!!!    INTEGER,   INTENT(OUT)                :: NUM_LINES_OBS
!!!!!!!!!!!    INTEGER,   PARAMETER                  :: NARGS = 4
!!!!!!!!!!!    INTEGER                               :: IERR, I, K, J
!!!!!!!!!!!    CHARACTER*800                         :: LINE
!!!!!!!!!!!    !
!!!!!!!!!!!    INTEGER                               :: NUMC
!!!!!!!!!!!    INTEGER                               :: BLENDSNUM
!!!!!!!!!!!    LOGICAL                               :: END_FIELD
!!!!!!!!!!!    ! Open
!!!!!!!!!!!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR&
!!!!!!!!!!!        , FORM='FORMATTED')
!!!!!!!!!!!    ! Read until last line
!!!!!!!!!!!    DO I=1,LINES_READ
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    !
!!!!!!!!!!!    ! Determine number of lines inside and allocate arrays
!!!!!!!!!!!    !
!!!!!!!!!!!    NUM_LINES_OBS = 0
!!!!!!!!!!!    END_FIELD=.FALSE.
!!!!!!!!!!!    DO WHILE (IERR.EQ.0)
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!       IF (IERR.NE.0) EXIT
!!!!!!!!!!!       DO J=1,LEN(LINE)
!!!!!!!!!!!         IF (':'.EQ.LINE(J:J)) THEN
!!!!!!!!!!!           END_FIELD=.TRUE.
!!!!!!!!!!!           EXIT
!!!!!!!!!!!         ENDIF
!!!!!!!!!!!       ENDDO
!!!!!!!!!!!       IF (END_FIELD.EQV..TRUE.) EXIT
!!!!!!!!!!!       NUM_LINES_OBS=NUM_LINES_OBS+1
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    CLOSE(1)
!!!!!!!!!!!    !
!!!!!!!!!!!    NUML=NUM_LINES_OBS
!!!!!!!!!!!    CALL ALLOCATE_1D_IP(IND_LINE,NUML,'ALLOCATE ind_line IN initial_duties.f90')
!!!!!!!!!!!    CALL ALLOCATE_1D_DP(WAVE_INI,NUML,'ALLOCATE wave_ini IN initial_duties.f90')
!!!!!!!!!!!    CALL ALLOCATE_1D_DP(DELTA_LAMBDA,NUML,'ALLOCATE delta_lambda IN initial_duties.f90')
!!!!!!!!!!!    CALL ALLOCATE_1D_IP(NUMWAVE,NUML,'ALLOCATE numwave IN initial_duties.f90')
!!!!!!!!!!!    !
!!!!!!!!!!!    ! Read data for observed spectral lines
!!!!!!!!!!!    !
!!!!!!!!!!!    ! Read until last line
!!!!!!!!!!!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR&
!!!!!!!!!!!        , FORM='FORMATTED')
!!!!!!!!!!!    DO I=1,LINES_READ
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    ! Read info about spectral lines
!!!!!!!!!!!    !
!!!!!!!!!!!    !
!!!!!!!!!!!    ! START A VARIABLE FOR THE MAXIMUM NUMBER OF WAVELENGTHS:
!!!!!!!!!!!    DO I=1,NUML
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!       NUMC=NCOLUMNS(LINE)
!!!!!!!!!!!       CALL SPLIT_READ_LINE_SPECTRAL(LINE &
!!!!!!!!!!!           , IND_LINE(I), NUMWAVE(I), WAVE_INI(I), DELTA_LAMBDA(I))
!!!!!!!!!!!       DO K=1,NUML_DATABASE
!!!!!!!!!!!         IF (LINE_NUM(K) .EQ. IND_LINE(I)) THEN
!!!!!!!!!!!           IND_LINE(I)=K
!!!!!!!!!!!           EXIT
!!!!!!!!!!!         ENDIF
!!!!!!!!!!!       ENDDO
!!!!!!!!!!!       CALL SPLIT_READ_LINE_BLENDS_NUMBER(LINE, BLENDSNUM)
!!!!!!!!!!!       IF (BLENDSNUM.GT.BLENDSMAX) BLENDSMAX=BLENDSNUM
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    CLOSE(1)
!!!!!!!!!!!    !
!!!!!!!!!!!    ! Now, repeat but storing spectral line info:
!!!!!!!!!!!    ! Read until last line
!!!!!!!!!!!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR&
!!!!!!!!!!!        , FORM='FORMATTED')
!!!!!!!!!!!    DO I=1,LINES_READ
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    !
!!!!!!!!!!!    ! NOW, WE ALREADY KNOW THE MAXIMUM NUMBER OF BLENDS TO BE CONSIDERED:
!!!!!!!!!!!    CALL ALLOCATE_2D_IP(BLENDSID,NUML,BLENDSMAX+1,'BLENDSID IN READ_LINES.')
!!!!!!!!!!!    CALL ALLOCATE_2D_DP(BLENDSDIFF,NUML,BLENDSMAX+1,'BLENDSDIFF IN READ_LINES.')
!!!!!!!!!!!    !
!!!!!!!!!!!    ! First column, reference line id
!!!!!!!!!!!    DO I=1,NUML
!!!!!!!!!!!      BLENDSID(I,1)=IND_LINE(I)
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    !
!!!!!!!!!!!    DO I=1,NUML
!!!!!!!!!!!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!!!!!!!!!!!       !
!!!!!!!!!!!       CALL SPLIT_READ_LINE_BLENDS(LINE, BLENDSMAX, BLENDSID(I,2:))
!!!!!!!!!!!       DO J=1,BLENDSMAX
!!!!!!!!!!!         IF (BLENDSID(I,J+1).GT.0) THEN
!!!!!!!!!!!           DO K=1,NUML_DATABASE
!!!!!!!!!!!             IF (LINE_NUM(K) .EQ. BLENDSID(I,J+1)) THEN
!!!!!!!!!!!               BLENDSID(I,J+1)=K
!!!!!!!!!!!               EXIT
!!!!!!!!!!!             ENDIF
!!!!!!!!!!!           ENDDO
!!!!!!!!!!!! DIFFERENCE OF WAVELENGTH IN ANGSTROM:
!!!!!!!!!!!           BLENDSDIFF(I,J+1)=LINE_L0(BLENDSID(I,J+1))-LINE_L0(IND_LINE(I))
!!!!!!!!!!!! IN mA:
!!!!!!!!!!!           BLENDSDIFF(I,J+1)=1.D3*BLENDSDIFF(I,J+1)
!!!!!!!!!!!         ENDIF
!!!!!!!!!!!       ENDDO
!!!!!!!!!!!    ENDDO
!!!!!!!!!!!    CLOSE(1)
!!!!!!!!!!!    !APY
!!!!!!!!!!!    NUMW=SUM(NUMWAVE)
!!!!!!!!!!!    !END APY
!!!!!!!!!!!    LINES_READ=LINES_READ+NUM_LINES_OBS
!!!!!!!!!!!  END SUBROUTINE READ_LINES
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_LINE_SPECTRAL(INLINE, SLID, INW, ISTART, ISTEP)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    !
    CHARACTER(*), INTENT(IN)       :: INLINE
    INTEGER, INTENT(INOUT)         :: SLID, INW
    REAL(DP), INTENT(INOUT)        :: ISTART, ISTEP
    !
    INTEGER                         :: I, J, NCOL, ENDENTRY
    !
    LOGICAL                         :: PRECHAR, SPECIALCHAR
    !
    INTEGER, ALLOCATABLE            :: INIT(:), FINI(:)
    !
    NCOL=NCOLUMNS(INLINE)
    !
    ! WHERE DO COLUMNS BEGIN?:
    CALL ALLOCATE_1D_IP(INIT,NCOL,'INIT ARRAY IN READ_LINES.')
    CALL ALLOCATE_1D_IP(FINI,NCOL,'FINI ARRAY IN READ_LINES.')
    !
    NCOL=0
    !
    PRECHAR=.FALSE.
    SPECIALCHAR=.FALSE.
    !
    DO I=1,LEN(INLINE)
      IF (INLINE(I:I).NE.' ') THEN
        ! Am i special?
        IF (INLINE(I:I).EQ.',') THEN
          SPECIALCHAR=.TRUE.
          CYCLE
        ENDIF
        !Was the previous thing a non empty space?
        IF ((PRECHAR.EQV..FALSE.).AND.(SPECIALCHAR.EQV..FALSE.)) THEN
            NCOL=NCOL+1
            INIT(NCOL)=I
            IF (NCOL.GT.1) FINI(NCOL-1)=I-1
        ELSE
          SPECIALCHAR=.FALSE.
        ENDIF
        ! Update prechar:
        PRECHAR=.TRUE.
      ELSE
        PRECHAR=.FALSE.
      ENDIF
    ENDDO
    FINI(NCOL)=LEN(INLINE)
    !
    ! NOW, WE ASSIGN THE VARIOUS THINGS:
    DO I=1,NCOL
!      PRINT*, INLINE(INIT(I):FINI(I)), I
      SELECT CASE(I)
         CASE(1)
           ! Here, we must be careful, only the main line (in case some ...
!            ... blends are present) has to be read:
           ENDENTRY=FINI(I)
           DO J=INIT(I),FINI(I)
             IF(INLINE(J:J).EQ.',') THEN
               ENDENTRY=J-1
               EXIT
             ENDIF
           ENDDO
!           PRINT*, 'AM I RIGHT?: ', ENDENTRY, INIT(I), FINI(I),&
!               INLINE(INIT(I):ENDENTRY)
           READ(INLINE(INIT(I):ENDENTRY), '(i10)') SLID
         CASE(2)
           READ(INLINE(INIT(I):FINI(I)), '(i10)') INW
         CASE(3)
           READ(INLINE(INIT(I):FINI(I)), '(d12.5)') ISTART
         CASE(4)
           READ(INLINE(INIT(I):FINI(I)), '(d12.5)') ISTEP
      END SELECT

    ENDDO
!PRINT*, 'INIT: ', INIT
!PRINT*, 'FINI: ', FINI
!PRINT*, INW, ISTART, ISTEP
    !
    DEALLOCATE(INIT)
    DEALLOCATE(FINI)

    !
  END SUBROUTINE SPLIT_READ_LINE_SPECTRAL
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_LINE_BLENDS_NUMBER(INLINE, NCOMMA)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    !
    CHARACTER(*), INTENT(IN)       :: INLINE
    INTEGER, INTENT(INOUT)         :: NCOMMA
    !
    INTEGER                         :: I, NCOL
    !
    LOGICAL                         :: PRECHAR, SPECIALCHAR
    !
    INTEGER, ALLOCATABLE            :: INIT(:), FINI(:)
    !
    NCOL=NCOLUMNS(INLINE)
    !
    ! WHERE DO COLUMNS BEGIN?:
    CALL ALLOCATE_1D_IP(INIT,NCOL,'INIT ARRAY IN READ_LINES.')
    CALL ALLOCATE_1D_IP(FINI,NCOL,'FINI ARRAY IN READ_LINES.')
    !
    NCOL=0
    NCOMMA=0
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
            NCOL=NCOL+1
            INIT(NCOL)=I
            IF (NCOL.GT.1) FINI(NCOL-1)=I-1
        ELSE
          SPECIALCHAR=.FALSE.
        ENDIF
        ! UPDATE PRECHAR:
        PRECHAR=.TRUE.
      ELSE
        PRECHAR=.FALSE.
      ENDIF
    ENDDO
    FINI(NCOL)=LEN(INLINE)
    !
    ! NOW, WE ARE ONLY INTERESTED IN THE FIRST ROW, AND IN PARTICULAR...
    ! ... WE ARE ONLY INTERESTED IN THE NUMBER OF COMMAS PRESENT:
    DO I=INIT(1),FINI(1)
      IF(INLINE(I:I).EQ.',') NCOMMA=NCOMMA+1
    ENDDO
    !
    DEALLOCATE(INIT)
    DEALLOCATE(FINI)
    !
  END SUBROUTINE SPLIT_READ_LINE_BLENDS_NUMBER
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_LINE_BLENDS(INLINE, BMAX, BID)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    !
    CHARACTER(*), INTENT(IN)       :: INLINE
    INTEGER, INTENT(INOUT)         :: BMAX
    INTEGER, INTENT(INOUT), DIMENSION(BMAX) :: BID
    !
    INTEGER                         :: I, NCOL
    INTEGER                         :: PRE, CNT
    !
    LOGICAL                         :: PRECHAR, SPECIALCHAR
    !
    INTEGER, ALLOCATABLE            :: INIT(:), FINI(:)
    !
    NCOL=NCOLUMNS(INLINE)
    !
    ! WHERE DO COLUMNS BEGIN?:
    CALL ALLOCATE_1D_IP(INIT,NCOL,'INIT ARRAY IN READ_LINES.')
    CALL ALLOCATE_1D_IP(FINI,NCOL,'FINI ARRAY IN READ_LINES.')
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
            NCOL=NCOL+1
            INIT(NCOL)=I
            IF (NCOL.GT.1) FINI(NCOL-1)=I-1
        ELSE
          SPECIALCHAR=.FALSE.
        ENDIF
        ! UPDATE PRECHAR:
        PRECHAR=.TRUE.
      ELSE
        PRECHAR=.FALSE.
      ENDIF
    ENDDO
    FINI(NCOL)=LEN(INLINE)
    !
    ! NOW, WE ARE ONLY INTERESTED IN THE FIRST ROW, AND IN PARTICULAR...
    ! ... WE ARE WANT TO STORE THE IDs OF THE VARIOUS BLENDS:
    PRE=0
    CNT=1
    BID(:)=-1
    DO I=INIT(1),FINI(1)
      IF(INLINE(I:I).EQ.',') THEN
        IF (PRE.EQ.0) THEN
          PRE=I
          CYCLE
        ELSE
          READ(INLINE(PRE+1:I-1), '(i10)') BID(CNT)
          CNT=CNT+1
          PRE=I
        ENDIF
      ENDIF
    ENDDO
    !
    IF (PRE.NE.0) READ(INLINE(PRE+1:FINI(1)), '(i10)') BID(CNT)
    !
    !PRINT*, BID(:), ';', CNT
    !
    DEALLOCATE(INIT)
    DEALLOCATE(FINI)
    !
  END SUBROUTINE SPLIT_READ_LINE_BLENDS
  !
  !________________________________________________
  !
  SUBROUTINE NSET_VAR_INVERSION(START_LINE)
    !
    USE CODE_MODES, ONLY: INPUTFILE
    USE GRID_PARAM, ONLY: NZ
    USE INVERT_PARAM
    !
    INTEGER, INTENT(INOUT)           ::START_LINE 
    !INTEGER                          :: NUM_LINES_OBS
    INTEGER                          :: IERR
    INTEGER                          :: I, J
    LOGICAL                          :: CONDIT, POWLOG, MITERLOG
    CHARACTER*800                    :: LINE
    !
    REAL(SP)                         :: WEIGHT
    INTEGER                          :: IVALUE
    !
    ! Initialize CYCPOW to default value if not provided:
    CYCPOW=2
    POWLOG=.FALSE.
    ! Initialize MAXSTEPS to default value if not provided:
    MAXSTEPS=5
    MITERLOG=.FALSE.
    !
    !NUM_LINES_OBS=0
    CONDIT=.TRUE.
    !DO WHILE (CONDIT.EQV..TRUE.)
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)

      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)

      ! Check if the section has finished:
      DO J=1,LEN(LINE)
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      IF (CONDIT.EQV..FALSE.) EXIT
      !NUM_LINES_OBS=NUM_LINES_OBS+1
      !
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'TEM', NZ, INV_ATMPAR(1) &
         , WEIGHT,NSLB_MAX(1))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'P0', NZ, INV_ATMPAR(8) &
         , WEIGHT,NSLB_MAX(8))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PGAS', NZ, INV_ATMPAR(2) &
         , WEIGHT,NSLB_MAX(2))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'RHO', NZ, INV_ATMPAR(3) &
         , WEIGHT,NSLB_MAX(3))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BX', NZ, INV_ATMPAR(4) &
         , WEIGHT,NSLB_MAX(4))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BY', NZ, INV_ATMPAR(5) &
         , WEIGHT,NSLB_MAX(5))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BZ', NZ, INV_ATMPAR(6) &
         , WEIGHT,NSLB_MAX(6))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'VLOS', NZ, INV_ATMPAR(7) &
         , WEIGHT,NSLB_MAX(7))
      !
      ! New:
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'POW', 1, POWLOG, WEIGHT, IVALUE)
      IF (POWLOG.EQV..TRUE.) CYCPOW=IVALUE
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'MAXITER', 1, MITERLOG, WEIGHT, IVALUE)
      IF (MITERLOG.EQV..TRUE.) MAXSTEPS=IVALUE
      ! End New.
      !
    ENDDO
    !
    !
IF (mpi__myrank.EQ.0) THEN
!PRINT*, 'POW: ', CYCPOW
!PRINT*, 'MAXSTEPS: ', MAXSTEPS
ENDIF
    !
    !
    !
  END SUBROUTINE NSET_VAR_INVERSION
  !
  !________________________________________________
  !
  SUBROUTINE NSET_VREGULARIZATION_TYPE(START_LINE)
    !
    USE CODE_MODES, ONLY: VREGULARIZATION, INPUTFILE
    USE GRID_PARAM, ONLY: NZ
    USE INVERT_PARAM
    !
    INTEGER, INTENT(INOUT)           :: START_LINE
    INTEGER                          :: IERR
    INTEGER                          :: I, J
    LOGICAL                          :: CONDIT, DUMLOG
    CHARACTER*800                    :: LINE
    !
    REAL(SP)                         :: DUMR
    INTEGER                          :: IVALUE
    !
    ! Switch on vertical regularization:
    VREGULARIZATION=.TRUE.
    !
    CONDIT=.TRUE.
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)

      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)

      ! Check if the section has finished:
      DO J=1,LEN(LINE)
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      IF (CONDIT.EQV..FALSE.) EXIT
      !NUM_LINES_OBS=NUM_LINES_OBS+1
      !
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'TEM', 1, DUMLOG &
         , DUMR,PEN_TYP(1))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'P0', 0, DUMLOG &
         , DUMR,PEN_TYP(8))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PGAS', 1, DUMLOG &
         , DUMR,PEN_TYP(2))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'RHO', 1, DUMLOG &
         , DUMR,PEN_TYP(3))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BX', 0, DUMLOG &
         , DUMR,PEN_TYP(4))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BY', 0, DUMLOG &
         , DUMR,PEN_TYP(5))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BZ', 1, DUMLOG &
         , DUMR,PEN_TYP(6))
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'VLOS', 1, DUMLOG &
         , DUMR,PEN_TYP(7))
      !
    ENDDO
    !CLOSE(1)
    !
    !
IF (mpi__myrank.EQ.0) THEN
PRINT*, 'Regularization: ', PEN_TYP
ENDIF
    !
    !
    !LINES_READ=LINES_READ+NUM_LINES_OBS
    !
  END SUBROUTINE NSET_VREGULARIZATION_TYPE
  !
  !________________________________________________
  !
  SUBROUTINE NSET_VREGULARIZATION_NORM(START_LINE)
    !
    USE CODE_MODES, ONLY: VREGULARIZATION
    USE GRID_PARAM, ONLY: NZ
    USE INVERT_PARAM
    !
    INTEGER, INTENT(INOUT)           :: START_LINE
    INTEGER                          :: IERR
    INTEGER                          :: I, J
    LOGICAL                          :: CONDIT, DUMLOG
    CHARACTER*800                    :: LINE
    !
    INTEGER                         :: DUMI
    !
    ! Switch on vertical regularization:
    VREGULARIZATION=.TRUE.
    !
    CONDIT=.TRUE.
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)

      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)

      ! Check if the section has finished:
      DO J=1,LEN(LINE)
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      IF (CONDIT.EQV..FALSE.) EXIT
      !NUM_LINES_OBS=NUM_LINES_OBS+1
      !
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'TEM', 1, DUMLOG &
         , PEN_ALP(1), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'P0', 0, DUMLOG &
         , PEN_ALP(8), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PGAS', 1, DUMLOG &
         , PEN_ALP(2), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'RHO', 1, DUMLOG &
         , PEN_ALP(3), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BX', 0, DUMLOG &
         , PEN_ALP(4), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BY', 0, DUMLOG &
         , PEN_ALP(5), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BZ', 1, DUMLOG &
         , PEN_ALP(6), DUMI)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'VLOS', 1, DUMLOG &
         , PEN_ALP(7), DUMI)
      !
    ENDDO
    !CLOSE(1)
    !
    !
IF (mpi__myrank.EQ.0) THEN
PRINT*, 'Regularization: ', PEN_ALP
ENDIF
    !
    !
    !LINES_READ=LINES_READ+NUM_LINES_OBS
    !
  END SUBROUTINE NSET_VREGULARIZATION_NORM
  !
  !________________________________________________
  !
!  SUBROUTINE SET_VAR_INVERSION(LINES_READ)
!    !
!    USE CODE_MODES, ONLY: INPUTFILE
!    USE GRID_PARAM, ONLY: NZ
!    USE INVERT_PARAM
!    !
!    INTEGER, INTENT(INOUT)           :: LINES_READ
!    INTEGER                          :: NUM_LINES_OBS
!    INTEGER                          :: IERR
!    INTEGER                          :: I, J
!    LOGICAL                          :: CONDIT, POWLOG, MITERLOG
!    CHARACTER*800                    :: LINE
!    !
!    REAL(SP)                         :: WEIGHT
!    INTEGER                          :: IVALUE
!    !
!    ! Open
!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!    ! Read until last line
!    DO I=1,LINES_READ
!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!    ENDDO
!    !
!    ! Initialize CYCPOW to default value if not provided:
!    CYCPOW=2
!    POWLOG=.FALSE.
!    ! Initialize MAXSTEPS to default value if not provided:
!    MAXSTEPS=5
!    MITERLOG=.FALSE.
!    !
!    NUM_LINES_OBS=0
!    CONDIT=.TRUE.
!    DO WHILE (CONDIT.EQV..TRUE.)
!      READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!      IF (IERR.LT.0) EXIT
!      ! CHECK IF THE SECTION HAS FINISHED:
!      DO J=1,LEN(LINE)
!        IF (LINE(J:J).EQ.':') THEN
!          CONDIT=.FALSE.
!          EXIT
!        ENDIF
!      ENDDO
!      IF (CONDIT.EQV..FALSE.) EXIT
!      NUM_LINES_OBS=NUM_LINES_OBS+1
!      !
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'TEM', NZ, INV_ATMPAR(1) &
!         , WEIGHT,NSLB_MAX(1))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'P0', NZ, INV_ATMPAR(8) &
!         , WEIGHT,NSLB_MAX(8))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PGAS', NZ, INV_ATMPAR(2) &
!         , WEIGHT,NSLB_MAX(2))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'RHO', NZ, INV_ATMPAR(3) &
!         , WEIGHT,NSLB_MAX(3))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BX', NZ, INV_ATMPAR(4) &
!         , WEIGHT,NSLB_MAX(4))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BY', NZ, INV_ATMPAR(5) &
!         , WEIGHT,NSLB_MAX(5))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BZ', NZ, INV_ATMPAR(6) &
!         , WEIGHT,NSLB_MAX(6))
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'VLOS', NZ, INV_ATMPAR(7) &
!         , WEIGHT,NSLB_MAX(7))
!      !
!      ! New:
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'POW', 1, POWLOG, WEIGHT, IVALUE)
!      IF (POWLOG.EQV..TRUE.) CYCPOW=IVALUE
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'MAXITER', 1, MITERLOG, WEIGHT, IVALUE)
!      IF (MITERLOG.EQV..TRUE.) MAXSTEPS=IVALUE
!      ! End New.
!      !
!    ENDDO
!    CLOSE(1)
!    !
!    !
!IF (mpi__myrank.EQ.0) THEN
!PRINT*, 'POW: ', CYCPOW
!PRINT*, 'MAXSTEPS: ', MAXSTEPS
!ENDIF
!    !
!!    ! TO BE REMOVED SOON:
!!    INV_TEMP=INV_ATMPAR(1)
!!    INV_P0NN=INV_ATMPAR(8)  
!!    INV_PGAS=INV_ATMPAR(2)
!!    INV_RHON=INV_ATMPAR(3)
!!    INV_BXNN=INV_ATMPAR(4)  
!!    INV_BYNN=INV_ATMPAR(5)  
!!    INV_BZNN=INV_ATMPAR(6)  
!!    INV_VLOS=INV_ATMPAR(7)
!!    !
!!    TBMAX=NSLB_MAX(1)
!!    P0BMAX=NSLB_MAX(8)
!!    PGBMAX=NSLB_MAX(2)
!!    RHOBMAX=NSLB_MAX(3)
!!    BXBMAX=NSLB_MAX(4)
!!    BYBMAX=NSLB_MAX(5)
!!    BZBMAX=NSLB_MAX(6)
!!    VZBMAX=NSLB_MAX(7)
!!    ! END TBR.
!    !
!    LINES_READ=LINES_READ+NUM_LINES_OBS
!    !
!  END SUBROUTINE SET_VAR_INVERSION
  !
  !________________________________________________
  !
  SUBROUTINE NSET_STK_INVERSION(START_LINE)
    !
    USE CODE_MODES, ONLY: INPUTFILE
    USE INVERT_PARAM
    !
    INTEGER, INTENT(INOUT)           :: START_LINE
    INTEGER                          :: IERR
    INTEGER                          :: I, J
    LOGICAL                          :: CONDIT
    CHARACTER*800                    :: LINE
    !
    REAL(SP)                         :: MINW
    INTEGER                          :: WEIGHT
    REAL(SP)                         :: DWEIGHT
    !
    CONDIT=.TRUE.
    !DO WHILE (CONDIT.EQV..TRUE.)
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)
      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)
      ! Check if the section has finished:
      DO J=1,LEN(LINE)
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      IF (CONDIT.EQV..FALSE.) EXIT
      !
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKI', 1, INV_STK(1), WSTK(1), WEIGHT, NSTKINV)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKQ', 1, INV_STK(2), WSTK(2), WEIGHT, NSTKINV)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKU', 1, INV_STK(3), WSTK(3), WEIGHT, NSTKINV)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKV', 1, INV_STK(4), WSTK(4), WEIGHT, NSTKINV)
      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'AUTO', 1, AUTOWEIGHT, DWEIGHT, WEIGHT)
      !
    ENDDO
    !
    MINW=SQRT(WSTK(1)**2+WSTK(2)**2+WSTK(3)**2+WSTK(4)**2)/SQRT(REAL(NSTKINV))
    !
    WSTK(1)=WSTK(1)/MINW
    WSTK(2)=WSTK(2)/MINW
    WSTK(3)=WSTK(3)/MINW
    WSTK(4)=WSTK(4)/MINW
    !
    !
  END SUBROUTINE NSET_STK_INVERSION
  !
  !
  !
!  SUBROUTINE SET_STK_INVERSION(LINES_READ)
!    !
!    USE CODE_MODES, ONLY: INPUTFILE
!    USE INVERT_PARAM
!    !
!    INTEGER, INTENT(INOUT)           :: LINES_READ
!    INTEGER                          :: NUM_LINES_OBS
!    INTEGER                          :: IERR
!    INTEGER                          :: I, J
!    LOGICAL                          :: CONDIT
!    CHARACTER*800                    :: LINE
!    !
!    REAL(SP)                         :: MINW
!    INTEGER                          :: WEIGHT
!    REAL(SP)                         :: DWEIGHT
!    !
!    ! Open
!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!    ! Read until last line
!    DO I=1,LINES_READ
!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!    ENDDO
!    !
!    NUM_LINES_OBS=0
!    CONDIT=.TRUE.
!    DO WHILE (CONDIT.EQV..TRUE.)
!      READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!      IF (IERR.LT.0) EXIT
!      ! CHECK IF THE SECTION HAS FINISHED:
!      DO J=1,LEN(LINE)
!        IF (LINE(J:J).EQ.':') THEN
!          CONDIT=.FALSE.
!          EXIT
!        ENDIF
!      ENDDO
!      IF (CONDIT.EQV..FALSE.) EXIT
!      !
!      NUM_LINES_OBS=NUM_LINES_OBS+1
!      !
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKI', 1, INV_STK(1), WSTK(1), WEIGHT, NSTKINV)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKQ', 1, INV_STK(2), WSTK(2), WEIGHT, NSTKINV)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKU', 1, INV_STK(3), WSTK(3), WEIGHT, NSTKINV)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'STKV', 1, INV_STK(4), WSTK(4), WEIGHT, NSTKINV)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'AUTO', 1, AUTOWEIGHT, DWEIGHT, WEIGHT)
!      !
!    ENDDO
!    CLOSE(1)
!    !
!    LINES_READ=LINES_READ+NUM_LINES_OBS
!    !
!    MINW=SQRT(WSTK(1)**2+WSTK(2)**2+WSTK(3)**2+WSTK(4)**2)/SQRT(REAL(NSTKINV))
!    !
!    WSTK(1)=WSTK(1)/MINW
!    WSTK(2)=WSTK(2)/MINW
!    WSTK(3)=WSTK(3)/MINW
!    WSTK(4)=WSTK(4)/MINW
!    !
!!!    ! TO BE REMOVED SOON:
!!    WSTKI=WSTK(1)
!!    WSTKQ=WSTK(2)
!!    WSTKU=WSTK(3)
!!    WSTKV=WSTK(4)
!!!    ! END TBRS
!    !
!  END SUBROUTINE SET_STK_INVERSION
  !
  !________________________________________________
  !
  SUBROUTINE GET_NUMBER(INLINE,IVAL,DVAL)
    !
    CHARACTER(*), INTENT(IN) :: INLINE 
    !
    INTEGER, INTENT(INOUT) :: IVAL
    REAL(DP), INTENT(INOUT) :: DVAL
    !
    INTEGER :: ISEXP
    INTEGER :: ISREAL
    INTEGER :: J
    !
    ISREAL=0
    ISEXP=0
    DO J=1,LEN(INLINE)
      IF(INLINE(J:J).EQ.'.') ISREAL=1
      IF(INLINE(J:J).EQ.'E') ISEXP=J
      IF(INLINE(J:J).EQ.'e') ISEXP=J
    ENDDO
    !
    IF (ISREAL.NE.1) THEN
      !IF (ILOGIC.EQV..TRUE.) THEN
        IF (ISEXP.EQ.0) THEN
          READ(INLINE, '(i14)') IVAL
          DVAL=REAL(IVAL)
        ELSE
          READ(INLINE(1:ISEXP-1), '(i14)') IVAL
          DVAL=REAL(IVAL)
          READ(INLINE(ISEXP+1:LEN(INLINE)), '(i14)') IVAL
          DVAL=DVAL*10.0E0**(IVAL)
        ENDIF
      !ENDIF
    ELSE ! The number is real:
      !IF (ILOGIC.EQV..TRUE.) THEN
        IF (ISEXP.EQ.0) THEN
          READ(INLINE, '(f14.4)') DVAL
        ELSE
          READ(INLINE(1:ISEXP-1), '(f14.4)') DVAL
          READ(INLINE(ISEXP+1:LEN(INLINE)), '(i14)') IVAL
          DVAL=DVAL*10.0D0**(IVAL)
        ENDIF
      !ENDIF
    ENDIF
    !
  END SUBROUTINE GET_NUMBER
  !
  !________________________________________________
  !
!!!!  SUBROUTINE SPLIT_READ_LINE(LINE, MATCH, NMAXCOL, TYPES, DEF, LOGIC, DVALUE, IVALUE, COUNTER)
!!!!    !
!!!!    CHARACTER(*), INTENT(IN)    :: LINE, MATCH
!!!!    INTEGER, INTENT(IN)         :: NMAXCOL
!!!!    INTEGER, INTENT(IN), DIMENSION(NMAXCOL) :: TYPES
!!!!    REAL(DP), INTENT(IN), DIMENSION(NMAXCOL) :: DEF
!!!!    !
!!!!    LOGICAL, INTENT(INOUT)      :: LOGIC
!!!!    !
!!!!    REAL(SP), INTENT(INOUT), DIMENSION(NMAXCOL) :: DVALUE
!!!!    INTEGER, INTENT(INOUT), DIMENSION(NMAXCOL) :: IVALUE
!!!!    INTEGER, INTENT(INOUT), OPTIONAL :: COUNTER
!!!!
!!!!  END SUBROUTINE SPLIT_READ_LINE
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_LINE_WEIGHTS(LINE, MATCH, DEF, LOGIC, DVALUE, IVALUE, COUNTER)
    !
    CHARACTER(*), INTENT(IN)    :: LINE, MATCH
    INTEGER, INTENT(IN)         :: DEF
    !
    LOGICAL, INTENT(INOUT)      :: LOGIC
    REAL(SP), INTENT(INOUT)     :: DVALUE
    INTEGER, INTENT(INOUT)      :: IVALUE
    INTEGER, INTENT(INOUT),OPTIONAL      :: COUNTER
    !
    INTEGER                     :: I, NCOL, START, PRE, ISREAL, ISEXP, DIVALUE, J
    LOGICAL                     :: ILOGIC
    REAL(DP)                     :: VALUED
    !
    NCOL=0
    START=0
    PRE=0
    ILOGIC=.FALSE.
    !DVALUE=0.D0
    !
    DO I=1,LEN(LINE)
! If we have an empty space, then we might need to interpret the thing:
!      IF (mpi__myrank.EQ.1) PRINT*, I, LINE(I:I)
      IF(LINE(I:I).EQ.' ') THEN
! We need to do something if there are previous non empty characters (pre=1)
        IF (PRE.EQ.1) THEN
          ! Assign term to each column:
          SELECT CASE(NCOL)
            CASE(1)
              IF (LINE(START:I).EQ.MATCH) THEN
                LOGIC=.TRUE.
                ILOGIC=.TRUE.
                IF (PRESENT(COUNTER)) COUNTER=COUNTER+1
! Standard value for weights, overwritten if supplied
                DVALUE=REAL(DEF)
                DIVALUE=DEF
              ELSE
                EXIT
              ENDIF
            CASE(2)
! FIRST, WE NEED TO CHECK IF THE NUMBER IS GIVEN AS INTEGER OR OTHERWISE:
              IF (ILOGIC.EQV..TRUE.) THEN
                CALL GET_NUMBER(LINE(START:I),DIVALUE,VALUED)
                DVALUE=REAL(VALUED)
              ENDIF
!!!!!!!!!!!!!!!!!!!!              ISREAL=0
!!!!!!!!!!!!!!!!!!!!              ISEXP=0
!!!!!!!!!!!!!!!!!!!!              DO J=1,LEN(LINE(START:I))
!!!!!!!!!!!!!!!!!!!!                IF(LINE(START+J:START+J).EQ.'.') ISREAL=1
!!!!!!!!!!!!!!!!!!!!                IF(LINE(START+J:START+J).EQ.'E') ISEXP=START+J
!!!!!!!!!!!!!!!!!!!!                IF(LINE(START+J:START+J).EQ.'e') ISEXP=START+J
!!!!!!!!!!!!!!!!!!!!              ENDDO
!!!!!!!!!!!!!!!!!!!!              !PRINT*, 'NUMBER?: ', LINE(START:I), ' ; Is it real?: ', ISREAL
!!!!!!!!!!!!!!!!!!!!              IF (ISREAL.NE.1) THEN
!!!!!!!!!!!!!!!!!!!!                IF (ILOGIC.EQV..TRUE.) THEN
!!!!!!!!!!!!!!!!!!!!                  IF (ISEXP.EQ.0) THEN
!!!!!!!!!!!!!!!!!!!!                    READ(LINE(START:I), '(i14)') DIVALUE
!!!!!!!!!!!!!!!!!!!!                    DVALUE=REAL(DIVALUE)
!!!!!!!!!!!!!!!!!!!!                  ELSE
!!!!!!!!!!!!!!!!!!!!                    READ(LINE(START:ISEXP-1), '(i14)') DIVALUE
!!!!!!!!!!!!!!!!!!!!                    DVALUE=REAL(DIVALUE)
!!!!!!!!!!!!!!!!!!!!                    READ(LINE(ISEXP+1:I), '(i14)') DIVALUE
!!!!!!!!!!!!!!!!!!!!                    DVALUE=DVALUE*10.0E0**(DIVALUE)
!!!!!!!!!!!!!!!!!!!!                  ENDIF
!!!!!!!!!!!!!!!!!!!!                ENDIF
!!!!!!!!!!!!!!!!!!!!              ELSE ! The number is real:
!!!!!!!!!!!!!!!!!!!!                IF (ILOGIC.EQV..TRUE.) THEN
!!!!!!!!!!!!!!!!!!!!                  IF (ISEXP.EQ.0) THEN
!!!!!!!!!!!!!!!!!!!!                    READ(LINE(START:I), '(f14.4)') DVALUE
!!!!!!!!!!!!!!!!!!!!                  ELSE
!!!!!!!!!!!!!!!!!!!!                    READ(LINE(START:ISEXP-1), '(f14.4)') DVALUE
!!!!!!!!!!!!!!!!!!!!                    !DVALUE=DBLE(DIVALUE)
!!!!!!!!!!!!!!!!!!!!                    READ(LINE(ISEXP+1:I), '(i14)') DIVALUE
!!!!!!!!!!!!!!!!!!!!                    DVALUE=DVALUE*10.0D0**(DIVALUE)
!!!!!!!!!!!!!!!!!!!!                  ENDIF
!!!!!!!!!!!!!!!!!!!!                ENDIF
!!!!!!!!!!!!!!!!!!!!              ENDIF
              
          END SELECT
          !END ASSIGNMENT
        ENDIF
        PRE=0
        CYCLE
      ENDIF
! IF IT IS NOT AN EMPTY SPACE AND THE PREVIOUS CHARACTER WAS AN EMPTY...
! ...CHARACTER, THEN ADD A COLUMN:
      IF (PRE.EQ.0) THEN
        NCOL=NCOL+1
        PRE=1
        START=I
      ENDIF

    ENDDO
    IF (ILOGIC.EQV..TRUE.) IVALUE=NINT(DVALUE)
    !
  END SUBROUTINE SPLIT_READ_LINE_WEIGHTS
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_LINE_PSF(LINE, MATCH, LOGIC, SVALUE)
    !
    CHARACTER(*), INTENT(IN)    :: LINE, MATCH
    !
    LOGICAL, INTENT(INOUT)      :: LOGIC
    CHARACTER*800, INTENT(INOUT)      :: SVALUE
    !
    INTEGER                     :: I, NCOL, START, PRE, ISREAL, DIVALUE, J
    LOGICAL                     :: ILOGIC
    CHARACTER*800               :: DUMVALUE
    !
    NCOL=0
    START=0
    PRE=0
    ILOGIC=.FALSE.
    !DVALUE=0.D0
    !
    DO I=1,LEN(LINE)
! IF WE HAVE AN EMPTY SPACE, THEN WE MIGHT NEED TO INTERPRET THE THING:
      IF(LINE(I:I).EQ.' ') THEN
! WE NEED TO DO SOMETHING IF THERE ARE PREVIOUS NON EMPTY CHARACTERS (PRE=1)
        IF (PRE.EQ.1) THEN
          !READ(LINE(START:I), '(A40)') STERM
          ! ASSIGN TERM TO EACH COLUMN:
          SELECT CASE(NCOL)
            CASE(1)
              IF (LINE(START:I).EQ.MATCH) THEN
                LOGIC=.TRUE.
                ILOGIC=.TRUE.
              ELSE
                EXIT
              ENDIF
            CASE(2)
! STORE NAME
              DUMVALUE=TRIM(LINE(START:I))
          END SELECT
          !END ASSIGNMENT
        ENDIF
        PRE=0
        CYCLE
      ENDIF
! IF IT IS NOT AN EMPTY SPACE AND THE PREVIOUS CHARACTER WAS AN EMPTY...
! ...CHARACTER, THEN ADD A COLUMN:
      IF (PRE.EQ.0) THEN
        NCOL=NCOL+1
        PRE=1
        START=I
      ENDIF
    ENDDO
    !
    IF (ILOGIC.EQV..TRUE.) SVALUE=DUMVALUE
    !
  END SUBROUTINE SPLIT_READ_LINE_PSF
  !
  !________________________________________________
  !
  SUBROUTINE NSET_CONTINUUM_NORMALIZATION(START_LINE)
    !
    USE CODE_MODES, ONLY: HSRA_NORMALIZATION, INPUTFILE
    !
    INTEGER, INTENT(INOUT)           :: START_LINE
    INTEGER                          :: NUM_LINES_OBS
    INTEGER                          :: IERR
    INTEGER                          :: I, J
    LOGICAL                          :: CONDIT
    CHARACTER*800                    :: LINE
    CHARACTER(LEN=1)                 :: FIRSTCHAR
    !
    NUM_LINES_OBS=0
    CONDIT=.TRUE.
    !DO WHILE (CONDIT.EQV..TRUE.)
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)
      !
      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)
      !
      IF (IERR.LT.0) EXIT
      ! CHECK IF THE SECTION HAS FINISHED:
      DO J=1,LEN(LINE)
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      !
      IF (CONDIT.EQV..TRUE.) NUM_LINES_OBS=NUM_LINES_OBS+1
      IF (CONDIT.EQV..TRUE.) FIRSTCHAR=TRIM(LINE)
    ENDDO
    ! CHECK IF ONE OF THE KEYWORDS IS PRESENT:
    IF (FIRSTCHAR.EQ.'T') HSRA_NORMALIZATION = .TRUE.
    IF (FIRSTCHAR.EQ.'F') HSRA_NORMALIZATION = .FALSE.
    !
    IF ((FIRSTCHAR.NE.'T').AND.(FIRSTCHAR.NE.'F')) THEN
      WRITE(*,*) '"CONTINUUM NORMALIZATION:" FIELD IS WRONG...'
      WRITE(*,*) '          IT MUST BE EITHER TRUE OR FALSE...'
      WRITE(*,*) '                                      NOT:'//TRIM(FIRSTCHAR)//'.'
      STOP
    ENDIF
    !
  END SUBROUTINE NSET_CONTINUUM_NORMALIZATION
  !
  !________________________________________________
  !
!  SUBROUTINE SET_CONTINUUM_NORMALIZATION(LINES_READ)
!    !
!    USE CODE_MODES, ONLY: HSRA_NORMALIZATION, INPUTFILE
!    !
!    INTEGER, INTENT(INOUT)           :: LINES_READ
!    INTEGER                          :: NUM_LINES_OBS
!    INTEGER                          :: IERR
!    INTEGER                          :: I, J
!    LOGICAL                          :: CONDIT
!    CHARACTER*800                    :: LINE
!    CHARACTER(LEN=1)                 :: FIRSTCHAR
!    ! Open
!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!    ! Read until last line
!    DO I=1,LINES_READ
!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!    ENDDO
!    !
!    NUM_LINES_OBS=0
!    CONDIT=.TRUE.
!    DO WHILE (CONDIT.EQV..TRUE.)
!      READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!      WRITE(*,*) LINE, ' L'
!      IF (IERR.LT.0) EXIT
!      ! CHECK IF THE SECTION HAS FINISHED:
!      DO J=1,LEN(LINE)
!        IF (LINE(J:J).EQ.':') THEN
!          CONDIT=.FALSE.
!          EXIT
!        ENDIF
!      ENDDO
!      !
!      IF (CONDIT.EQV..TRUE.) NUM_LINES_OBS=NUM_LINES_OBS+1
!      IF (CONDIT.EQV..TRUE.) FIRSTCHAR=TRIM(LINE)
!    ENDDO
!    CLOSE(1)
!    ! CHECK IF ONE OF THE KEYWORDS IS PRESENT:
!    IF (FIRSTCHAR.EQ.'T') HSRA_NORMALIZATION = .TRUE.
!    IF (FIRSTCHAR.EQ.'F') HSRA_NORMALIZATION = .FALSE.
!    !
!    IF ((FIRSTCHAR.NE.'T').AND.(FIRSTCHAR.NE.'F')) THEN
!      WRITE(*,*) '"CONTINUUM NORMALIZATION:" FIELD IS WRONG...'
!      WRITE(*,*) '          IT MUST BE EITHER TRUE OR FALSE...'
!      WRITE(*,*) '                                      NOT:'//TRIM(FIRSTCHAR)//'.'
!      STOP
!    ENDIF
!    !
!    LINES_READ=LINES_READ+NUM_LINES_OBS
!    !
!  END SUBROUTINE SET_CONTINUUM_NORMALIZATION
  !
  !________________________________________________
  !
  SUBROUTINE NSET_COUPLED_INVERSION(START_LINE)
    !
    USE COUPLED_PARAM, ONLY: COU_BLCKSZ, COU_NMTHRD, COU_PSFRAD, COU_PSFFNAME
    USE user_mpi, ONLY: mpi__size
    USE LINES_DATABASE, ONLY: NCOLUMNS
    USE CODE_MODES, ONLY: COUPLED, INPUTFILE
    !
    INTEGER, INTENT(INOUT)           :: START_LINE
    INTEGER                          :: NUM_LINES_OBS
    INTEGER                          :: IERR
    INTEGER                          :: I, J
    INTEGER                          :: NUMC
    LOGICAL                          :: CONDIT
    LOGICAL                          :: ISFILENAME
    CHARACTER*800                    :: LINE
    CHARACTER(LEN=1)                 :: FIRSTCHAR
    LOGICAL                          :: DUMLOG
    REAL(SP)                         :: DUMSP
    !
    ! Define defaults:
    COU_BLCKSZ=10 ! MvN's advice. Very large -> SVD will not converge
    COU_NMTHRD=2*mpi__size
    COU_PSFRAD=3
    !
    !
    NUM_LINES_OBS=0
    CONDIT=.TRUE.
    ISFILENAME=.FALSE.
    COU_PSFFNAME=""
    !
    !DO WHILE (CONDIT.EQV..TRUE.)
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)
      !
      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)
      ! Check if the section has finished:
      DO J=1,LEN(LINE)
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      !
      IF (CONDIT.EQV..FALSE.) EXIT
      !
      NUM_LINES_OBS=NUM_LINES_OBS+1
      !
      NUMC=NCOLUMNS(LINE)
      IF (NUMC.EQ.1) THEN
        ! CHECK IF IS IT COUPLED:
        FIRSTCHAR=TRIM(LINE)
        IF ((FIRSTCHAR.EQ.'y').OR.(FIRSTCHAR.EQ.'Y')) THEN
          COUPLED=.TRUE.
        ENDIF
      ELSE IF (NUMC.EQ.2) THEN
        ! READ THE FILENAME OF THE PSF FILE:
        CALL SPLIT_READ_LINE_PSF(LINE, 'FILENAME', ISFILENAME, COU_PSFFNAME)
        CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BLCKSZ', 0, DUMLOG, DUMSP, COU_BLCKSZ)
        CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'NMTHRD', 0, DUMLOG, DUMSP, COU_NMTHRD)
        CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PSFRAD', 0, DUMLOG, DUMSP, COU_PSFRAD)
      ELSE
IF (mpi__myrank.EQ.0) THEN
        PRINT*, ''
        PRINT*, ''
        PRINT*, '"COUPLED INVERSION:" FIELD IN INPUT FILE MUST! CONTAIN...'
        PRINT*, '                                            TRUE/FALSE...'
        PRINT*, '                 FILENAME "name_of_the_file_with_the_psf"'
        PRINT*, ''
        STOP
ENDIF
      ENDIF
      !
    ENDDO
    !
    IF ((COUPLED.EQV..TRUE.).AND.(ISFILENAME.EQV..FALSE.)) THEN
IF (mpi__myrank.EQ.0) THEN
      PRINT*, ''
      PRINT*, '"COUPLED INVERSION:" FIELD MUST HAVE TWO LINES.'
      PRINT*, '   AN EXAMPLE OF HOW THIS PIECE OF INPUT FILE SHOULD LOOK LIKE:'
      PRINT*, ' "...'
      PRINT*, ' COUPLED INVERSION:'
      PRINT*, ' yes/no'
      PRINT*, ' FILENAME name_of_file_with_spatial_psf'
      PRINT*, ' ..."'
      PRINT*, ' IF no, THEN THE SECOND LINE IS NOT NEEDED ANYMORE BUT IT MAY...'
      PRINT*, '...WELL BE THERE, IT DOES NOT MATTER'
      PRINT*, ''
      STOP
ENDIF
    ENDIF
!/IF (mpi__myrank.EQ.0) THEN
!/    PRINT*, 'PSF FILE: ', TRIM(COU_PSFFNAME)
!/    PRINT*, 'COU_BLCKSZ: ', COU_BLCKSZ
!/    PRINT*, 'COU_NMTHRD: ', COU_NMTHRD, mpi__size
!/    PRINT*, 'COU_PSFRAD: ', COU_PSFRAD
!/ENDIF
    !
  END SUBROUTINE NSET_COUPLED_INVERSION
  !
  !________________________________________________
  !
!  SUBROUTINE SET_COUPLED_INVERSION(LINES_READ)
!    !
!    USE COUPLED_PARAM, ONLY: COU_BLCKSZ, COU_NMTHRD, COU_PSFRAD, COU_PSFFNAME
!    USE user_mpi, ONLY: mpi__size
!    USE LINES_DATABASE, ONLY: NCOLUMNS
!    USE CODE_MODES, ONLY: COUPLED, INPUTFILE
!    !
!    INTEGER, INTENT(INOUT)           :: LINES_READ
!    INTEGER                          :: NUM_LINES_OBS
!    INTEGER                          :: IERR
!    INTEGER                          :: I, J
!    INTEGER                          :: NUMC
!    LOGICAL                          :: CONDIT
!    LOGICAL                          :: ISFILENAME
!    CHARACTER*800                    :: LINE
!    CHARACTER(LEN=1)                 :: FIRSTCHAR
!    LOGICAL                          :: DUMLOG
!    REAL(SP)                         :: DUMSP
!    !
!    ! Define defaults:
!    COU_BLCKSZ=10 ! MvN's advice. Very large -> SVD will not converge
!    COU_NMTHRD=2*mpi__size
!    COU_PSFRAD=3
!    !
!    !
!    ! Open
!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!    ! Read until last line
!    DO I=1,LINES_READ
!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!    ENDDO
!    !
!    NUM_LINES_OBS=0
!    CONDIT=.TRUE.
!    ISFILENAME=.FALSE.
!    COU_PSFFNAME=""
!    DO WHILE (CONDIT.EQV..TRUE.)
!      READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!      !WRITE(*,*) LINE, ' L'
!      IF (IERR.LT.0) EXIT
!      ! CHECK IF THE SECTION HAS FINISHED:
!      DO J=1,LEN(LINE)
!        IF (LINE(J:J).EQ.':') THEN
!          CONDIT=.FALSE.
!          EXIT
!        ENDIF
!      ENDDO
!      !
!      IF (CONDIT.EQV..FALSE.) EXIT
!      !
!      NUM_LINES_OBS=NUM_LINES_OBS+1
!      !
!      NUMC=NCOLUMNS(LINE)
!      IF (NUMC.EQ.1) THEN
!        ! CHECK IF IS IT COUPLED:
!        FIRSTCHAR=TRIM(LINE)
!        IF ((FIRSTCHAR.EQ.'y').OR.(FIRSTCHAR.EQ.'Y')) THEN
!          COUPLED=.TRUE.
!        ENDIF
!      ELSE IF (NUMC.EQ.2) THEN
!        ! READ THE FILENAME OF THE PSF FILE:
!        CALL SPLIT_READ_LINE_PSF(LINE, 'FILENAME', ISFILENAME, COU_PSFFNAME)
!        CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BLCKSZ', 0, DUMLOG, DUMSP, COU_BLCKSZ)
!        CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'NMTHRD', 0, DUMLOG, DUMSP, COU_NMTHRD)
!        CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PSFRAD', 0, DUMLOG, DUMSP, COU_PSFRAD)
!      ELSE
!IF (mpi__myrank.EQ.0) THEN
!        PRINT*, ''
!        PRINT*, ''
!        PRINT*, '"COUPLED INVERSION:" FIELD IN INPUT FILE MUST! CONTAIN...'
!        PRINT*, '                                            TRUE/FALSE...'
!        PRINT*, '                 FILENAME "name_of_the_file_with_the_psf"'
!        PRINT*, ''
!        STOP
!ENDIF
!      ENDIF
!      !
!    ENDDO
!    CLOSE(1)
!    !
!    IF ((COUPLED.EQV..TRUE.).AND.(ISFILENAME.EQV..FALSE.)) THEN
!IF (mpi__myrank.EQ.0) THEN
!      PRINT*, ''
!      PRINT*, '"COUPLED INVERSION:" FIELD MUST HAVE TWO LINES.'
!      PRINT*, '   AN EXAMPLE OF HOW THIS PIECE OF INPUT FILE SHOULD LOOK LIKE:'
!      PRINT*, ' "...'
!      PRINT*, ' COUPLED INVERSION:'
!      PRINT*, ' yes/no'
!      PRINT*, ' FILENAME name_of_file_with_spatial_psf'
!      PRINT*, ' ..."'
!      PRINT*, ' IF no, THEN THE SECOND LINE IS NOT NEEDED ANYMORE BUT IT MAY...'
!      PRINT*, '...WELL BE THERE, IT DOES NOT MATTER'
!      PRINT*, ''
!      STOP
!ENDIF
!    ENDIF
!!><    PRINT*, 'COUPLED: ', COUPLED
!!><    PRINT*, 'ISFILENAME: ', ISFILENAME
!!><    PRINT*, 'COU_PSFFNAME: ', TRIM(COU_PSFFNAME), LEN(COU_PSFFNAME)
!IF (mpi__myrank.EQ.0) THEN
!    PRINT*, 'COU_BLCKSZ: ', COU_BLCKSZ
!    PRINT*, 'COU_NMTHRD: ', COU_NMTHRD, mpi__size
!    PRINT*, 'COU_PSFRAD: ', COU_PSFRAD
!ENDIF
!!><    STOP
!    !
!    LINES_READ=LINES_READ+NUM_LINES_OBS
!    !
!  END SUBROUTINE SET_COUPLED_INVERSION
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_LSF_LINE(INLINE, LSFSIGMA, LSFW0)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    !
    CHARACTER(*), INTENT(IN)       :: INLINE
    REAL(DP), INTENT(INOUT)        :: LSFSIGMA, LSFW0
    !
    INTEGER                         :: I, J, NCOL, ENDENTRY, IISIGMA
    !
    LOGICAL                         :: PRECHAR, SPECIALCHAR, ISREAL
    !
    INTEGER, ALLOCATABLE            :: INIT(:), FINI(:)
    !
    NCOL=NCOLUMNS(INLINE)
    !
    ! WHERE DO COLUMNS BEGIN?:
    CALL ALLOCATE_1D_IP(INIT,NCOL,'INIT ARRAY IN READ_LINES.')
    CALL ALLOCATE_1D_IP(FINI,NCOL,'FINI ARRAY IN READ_LINES.')
    !
    NCOL=0
    !
    PRECHAR=.FALSE.
    SPECIALCHAR=.FALSE.
    !
    DO I=1,LEN(INLINE)
      IF (INLINE(I:I).NE.' ') THEN
        ! AM I SPECIAL?
        IF (INLINE(I:I).EQ.'#') THEN
          SPECIALCHAR=.TRUE.
          CYCLE
        ENDIF
        !WAS THE PREVIOUS THING A NON EMPTY SPACE?
        IF ((PRECHAR.EQV..FALSE.).AND.(SPECIALCHAR.EQV..FALSE.)) THEN
            NCOL=NCOL+1
            INIT(NCOL)=I
            IF (NCOL.GT.1) FINI(NCOL-1)=I-1
        ELSE
          SPECIALCHAR=.FALSE.
        ENDIF
        ! UPDATE PRECHAR:
        PRECHAR=.TRUE.
      ELSE
        PRECHAR=.FALSE.
      ENDIF
      IF (SPECIALCHAR.EQV..TRUE.) EXIT
    ENDDO
    IF (SPECIALCHAR.EQV..TRUE.) THEN
      FINI(NCOL)=I-2
    ELSE
      FINI(NCOL)=I-1!LEN(INLINE)
    ENDIF
    !
    IF (NCOL.NE.2) THEN
IF (mpi__myrank.EQ.0) THEN
      PRINT*, '   ***   '
      PRINT*, ' Wrong format! '
      PRINT*, INLINE
      PRINT*, '   ***   '
      STOP
ENDIF
    ENDIF
    !
    ! NOW, WE ASSIGN THE VARIOUS THINGS:
    DO I=1,NCOL
!      PRINT*, INLINE(INIT(I):FINI(I)), I
      SELECT CASE(I)
         CASE(1)
         ! CHECK IF INPUT IS REAL OR INTEGER:
           ISREAL=.FALSE.
           DO J=1,LEN(INLINE(INIT(I):FINI(I)))-1
             IF(INLINE(INIT(I)+J:INIT(I)+J).EQ.'.') ISREAL=.TRUE.
           ENDDO
           IF (ISREAL.EQV..FALSE.) THEN
             READ(INLINE(INIT(I):FINI(I)), '(i14)') IISIGMA
             LSFSIGMA=DBLE(IISIGMA)
           ELSE
             READ(INLINE(INIT(I):FINI(I)), '(d12.5)') LSFSIGMA
           ENDIF
         CASE(2)
         ! CHECK IF INPUT IS REAL OR INTEGER:
           ISREAL=.FALSE.
           DO J=1,LEN(INLINE(INIT(I):FINI(I)))-1
             IF(INLINE(INIT(I)+J:INIT(I)+J).EQ.'.') ISREAL=.TRUE.
           ENDDO
           IF (ISREAL.EQV..FALSE.) THEN
             READ(INLINE(INIT(I):FINI(I)), '(i14)') IISIGMA
             LSFW0=DBLE(IISIGMA)
           ELSE
!PRINT*, INLINE(INIT(I):FINI(I))
             READ(INLINE(INIT(I):FINI(I)), '(d12.5)') LSFW0
           ENDIF
      END SELECT

    ENDDO
!PRINT*, 'INIT: ', INIT
!PRINT*, 'FINI: ', FINI
!PRINT*, INW, ISTART, ISTEP
    !
    DEALLOCATE(INIT)
    DEALLOCATE(FINI)

    !
  END SUBROUTINE SPLIT_LSF_LINE
  !
  !________________________________________________
  !
  SUBROUTINE NSET_LSF(START_LINE)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_DP
    USE FORWARD_PARAM, ONLY: LSF_SIGMA, LSF_W0, BLENDSID, BLENDSMAX, BLENDSDIFF &
        , NUML, LSF_FILE
    USE CODE_MODES, ONLY: MLSF, INPUTFILE, MVERBOSE, DATAPATH
    USE user_mpi, ONLY: mpi__myrank, MPI_COMM_WORLD
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    !
    INTEGER, INTENT(INOUT)           :: START_LINE
    INTEGER                          :: NUM_LINES_OBS
    INTEGER                          :: IERR
    INTEGER                          :: I, J, NCOL
    LOGICAL                          :: CONDIT
    CHARACTER*800                    :: LINE
    CHARACTER*5                      :: empty='EMPTY'
    CHARACTER*5                      :: check='     '
    LOGICAL :: empty_field
    LOGICAL                          :: EXISTS
    !
    REAL(DP)                         :: IS, IW
    !
    CALL ALLOCATE_1D_DP(LSF_SIGMA, NUML, 'LSF_SIGMA')
    CALL ALLOCATE_1D_DP(LSF_W0, NUML, 'LSF_W0')
    LSF_SIGMA(:)=-1.0D0
    LSF_W0(:)=-1.0D0
    ALLOCATE(LSF_FILE(NUML))
    DO I=1,NUML
      LSF_FILE(I)=''
    ENDDO
    !
    NUM_LINES_OBS=0
    CONDIT=.TRUE.
    I=0
    !DO WHILE (CONDIT.EQV..TRUE.)
    DO WHILE (START_LINE.LT.INPUT_READ_NUMBER)
      START_LINE=START_LINE+1
      LINE=INPUT_READ_LINES(START_LINE)
      ! Check if the section has finished:
      DO J=1,800
        IF (LINE(J:J).EQ.':') THEN
          CONDIT=.FALSE.
          START_LINE=START_LINE-1
          EXIT
        ENDIF
      ENDDO
      IF (CONDIT.EQV..FALSE.) EXIT
      !
      NCOL=NCOLUMNS(LINE)
      IF (NCOL.EQ.2) THEN
        CALL SPLIT_LSF_LINE(LINE, IS, IW)
      ELSE IF (NCOL.EQ.1) THEN
        ! It is either empty or a file name:
        PRINT*, ''
      ELSE
        PRINT*, ''
        PRINT*, ' ** Wrong input format in: "LINE SPREAD FUNCTION"'
        PRINT*, ' -- Line: "'//TRIM(ADJUSTL(LINE))//'"!...'
        PRINT*, ' ...does not have either 1 or 2 fields.'
        PRINT*, ''
        PRINT*, ' Stopping'
        PRINT*, ''
        STOP
      ENDIF
      !
      NUM_LINES_OBS=NUM_LINES_OBS+1
      I=I+1
      !
      IF (NUM_LINES_OBS.LE.NUML) THEN
        IF (NCOL.EQ.2) THEN
          IF (IS.EQ.0) THEN
            PRINT*, ''
            PRINT*, ' Error!'
            PRINT*, ' LSF standard deviation for spectral region: ', I
            PRINT*, ' ...is 0. If gaussian LSF is provided, this value...'
            PRINT*, ' ...cannot be zero. Check line: "'//TRIM(LINE)//'"'
            PRINT*, ''
            PRINT*, ''
            STOP
          ENDIF
          LSF_SIGMA(I)=DBLE(IS)
          LSF_W0(I)=DBLE(IW)
        ELSE IF (NCOL.EQ.1) THEN
          ! It is either empty or a file name:
          LINE=ADJUSTL(LINE)
          empty_field = .FALSE.
          IF (LEN(TRIM(LINE)).EQ.5) THEN
            empty_field = .TRUE.
            DO J=1,LEN(TRIM(LINE))
              check(J:J)=LINE(J:J)
              IF (IACHAR(LINE(J:J)).GT.IACHAR('Z') ) check(J:J)=ACHAR(IACHAR(LINE(J:J))+(IACHAR('A')-IACHAR('a')))
              IF (check(J:J).NE.empty(J:J)) empty_field=.FALSE.
            ENDDO
          ENDIF
          IF (.NOT.empty_field) THEN
            DO J=1,LEN(TRIM(LINE))
              LSF_FILE(I)(J:J)=LINE(J:J)
            ENDDO
            !
            ! Check if file exists:
            INQUIRE(FILE=TRIM(DATAPATH)//"/./"//TRIM(LSF_FILE(I)), EXIST=EXISTS)
            IF (EXISTS.EQV..FALSE.) THEN
              PRINT*, ""
              PRINT*, "Error!"
              PRINT*, "I cannot find LSF file: "//TRIM(LSF_FILE(I))//" in directory: "//TRIM(DATAPATH)
              PRINT*, ""
              CALL MPI_ABORT(MPI_COMM_WORLD)
            ENDIF ! Writing error

          ENDIF

        ENDIF
      ENDIF

      !
    ENDDO
    !
    IF (NUM_LINES_OBS.NE.NUML) THEN
IF (mpi__myrank.EQ.0) THEN
      PRINT*, '   ***   '
      PRINT*, ' The number of LSF (', NUM_LINES_OBS, ') and LINE (' &
          , NUML, ') do not match!'
      PRINT*, '   ***   '
      STOP
ENDIF
    ENDIF
    !
    MLSF=.TRUE.
    !LINES_READ=LINES_READ+NUM_LINES_OBS
    !
    IF ((mpi__myrank.EQ.0).AND.(MVERBOSE.GT.0)) THEN
      PRINT*, 'NUML: ', NUML
      PRINT*, 'BLENDSID: '
      DO I=1,NUML
        PRINT*, '(', I, '): ', BLENDSID(I,:)
      ENDDO
      PRINT*, 'BLENDSMAX: ', BLENDSMAX
      PRINT*, 'BLENDSDIFF: '
      DO I=1,NUML
        PRINT*, '(', I, '): ', BLENDSDIFF(I,:)
      ENDDO
      !
      PRINT*, '   IS   ', '   IW   ', '   FILE   '
      DO I=1,NUML
        PRINT*, LSF_SIGMA(I), LSF_W0(I), TRIM(LSF_FILE(I))
      ENDDO
    ENDIF

    !
  END SUBROUTINE NSET_LSF
  !________________________________________________
  !
!  SUBROUTINE SET_LSF(LINES_READ)
!    !
!    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_DP
!    USE FORWARD_PARAM, ONLY: LSF_SIGMA, LSF_W0, BLENDSID, BLENDSMAX, BLENDSDIFF &
!        , NUML
!    USE CODE_MODES, ONLY: MLSF, INPUTFILE, MVERBOSE
!    USE user_mpi, ONLY: mpi__myrank
!    !
!    !
!    INTEGER, INTENT(INOUT)           :: LINES_READ
!    INTEGER                          :: NUM_LINES_OBS
!    INTEGER                          :: IERR
!    INTEGER                          :: I, J
!    LOGICAL                          :: CONDIT
!    CHARACTER*800                    :: LINE
!    !
!    REAL(DP)                         :: IS, IW
!    !
!    CALL ALLOCATE_1D_DP(LSF_SIGMA, NUML, 'LSF_SIGMA')
!    CALL ALLOCATE_1D_DP(LSF_W0, NUML, 'LSF_W0')
!    !
!    ! Open
!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!    ! Read until last line
!    DO I=1,LINES_READ
!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!    ENDDO
!    !
!    NUM_LINES_OBS=0
!    CONDIT=.TRUE.
!    I=0
!    DO WHILE (CONDIT.EQV..TRUE.)
!      READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!      IF (IERR.LT.0) EXIT
!      ! CHECK IF THE SECTION HAS FINISHED:
!      DO J=1,800
!        IF (LINE(J:J).EQ.':') THEN
!          CONDIT=.FALSE.
!          EXIT
!        ENDIF
!      ENDDO
!      IF (CONDIT.EQV..FALSE.) EXIT
!      !
!      CALL SPLIT_LSF_LINE(LINE, IS, IW)
!      NUM_LINES_OBS=NUM_LINES_OBS+1
!      I=I+1
!      !
!      IF (NUM_LINES_OBS.LE.NUML) THEN
!        LSF_SIGMA(I)=DBLE(IS)
!        LSF_W0(I)=DBLE(IW)
!      ENDIF
!      !
!    ENDDO
!    CLOSE(1)
!    !
!    IF (NUM_LINES_OBS.NE.NUML) THEN
!IF (mpi__myrank.EQ.0) THEN
!      PRINT*, '   ***   '
!      PRINT*, ' The number of LSF (', NUM_LINES_OBS, ') and LINE (' &
!          , NUML, ') do not match!'
!      PRINT*, '   ***   '
!      STOP
!ENDIF
!    ENDIF
!    !
!    MLSF=.TRUE.
!    LINES_READ=LINES_READ+NUM_LINES_OBS
!    !
!    IF ((mpi__myrank.EQ.0).AND.(MVERBOSE.GT.0)) THEN
!      PRINT*, 'NUML: ', NUML
!      PRINT*, 'BLENDSID: '
!      DO I=1,NUML
!        PRINT*, '(', I, '): ', BLENDSID(I,:)
!      ENDDO
!      PRINT*, 'BLENDSMAX: ', BLENDSMAX
!      PRINT*, 'BLENDSDIFF: '
!      DO I=1,NUML
!        PRINT*, '(', I, '): ', BLENDSDIFF(I,:)
!      ENDDO
!      !
!      PRINT*, 'IS: ', LSF_SIGMA
!      PRINT*, 'IW: ', LSF_W0
!    ENDIF
!    !
!  END SUBROUTINE SET_LSF
  !
  !________________________________________________
  !
  SUBROUTINE SPLIT_READ_LINE_STRING(LINE, MATCH, SVALUE)
    !
    CHARACTER(*), INTENT(IN)    :: LINE, MATCH
    !
    CHARACTER(*), INTENT(INOUT) :: SVALUE
    !
    INTEGER                     :: I, NCOL, START, PRE, J
    LOGICAL                     :: ILOGIC
    !
    NCOL=0
    START=0
    PRE=0
    ILOGIC=.FALSE.
    !
    DO I=1,LEN(LINE)
! If we have an empty space, then we might need to interpret the thing:
      IF(LINE(I:I).EQ.' ') THEN
! We need to do something if there are previous non empty characters (pre=1)
        IF (PRE.EQ.1) THEN
          !
          SELECT CASE(NCOL)
            CASE(1)
              IF (LINE(START:I).EQ.MATCH) THEN
                ILOGIC=.TRUE.
              ELSE
                EXIT
              ENDIF
            CASE(2)
! Actual path:
              IF (ILOGIC.EQV..TRUE.) THEN
                IF (LEN(TRIM(LINE(START:I))).GT.0) READ(LINE(START:I), '(A)') SVALUE
              ENDIF
          END SELECT
          !End assignment
        ENDIF
        PRE=0
        CYCLE
      ENDIF
! If it is not an empty space and the previous character was an empty...
! ...character, then add a column:
      IF (PRE.EQ.0) THEN
        NCOL=NCOL+1
        PRE=1
        START=I
      ENDIF
      !
    ENDDO
    !
  END SUBROUTINE SPLIT_READ_LINE_STRING
  !
  !________________________________________________
  !
!  SUBROUTINE SET_MISC(LINES_READ)
!    !
!    USE user_mpi
!    USE INVERT_PARAM, ONLY: ASSIST_T, ASSIST_P, ASSIST_B, ASSIST_V &
!        , RSVDTOL, RSIGMAP, ISIGMA
!    USE FORWARD_PARAM, ONLY: FULL_STOKES
!    USE CODE_MODES, ONLY: INPUTFILE, DATAPATH, OUTPPATH, LINEPATH &
!        , MODLPATH, MVERBOSE
!    !
!    INTEGER, INTENT(INOUT)           :: LINES_READ
!    INTEGER                          :: NUM_LINES_OBS
!    INTEGER                          :: IERR
!    INTEGER                          :: I, J
!    LOGICAL                          :: CONDIT
!    CHARACTER*800                    :: LINE
!    REAL(SP)                         :: DUMSP
!    INTEGER                          :: DUMIP
!    LOGICAL                          :: DUMSVDL, DUMSIGL
!    !
!    !
!    ! Open
!    OPEN(UNIT=1, FILE=TRIM(INPUTFILE), STATUS='OLD', IOSTAT=IERR, FORM='FORMATTED')
!    ! Read until last line
!    DO I=1,LINES_READ
!       READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!    ENDDO
!    !
!    NUM_LINES_OBS=0
!    CONDIT=.TRUE.
!    I=0
!    !
!    DO WHILE (CONDIT.EQV..TRUE.)
!      !
!      DUMSVDL=.FALSE.
!      DUMSIGL=.FALSE.
!      !
!      READ(UNIT=1, FMT='(A)',IOSTAT=IERR) LINE
!      IF (IERR.LT.0) EXIT
!      ! CHECK IF THE SECTION HAS FINISHED:
!      DO J=1,800
!        IF (LINE(J:J).EQ.':') THEN
!          CONDIT=.FALSE.
!          EXIT
!        ENDIF
!      ENDDO
!      IF (CONDIT.EQV..FALSE.) EXIT
!      !
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'TEMP', 0, ASSIST_T, DUMSP, DUMIP)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'PGAS', 0, ASSIST_P, DUMSP, DUMIP)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'BVEC', 0, ASSIST_B, DUMSP, DUMIP)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'VLOS', 0, ASSIST_V, DUMSP, DUMIP)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'FULL_STOKES', 0, FULL_STOKES, DUMSP, DUMIP)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'SVDTOL', 0, DUMSVDL, DUMSP, DUMIP)
!      IF (DUMSVDL.EQV..TRUE.) RSVDTOL=DBLE(DUMSP)
!      CALL SPLIT_READ_LINE_WEIGHTS(LINE, 'NOISE', 0, DUMSIGL, DUMSP, DUMIP)
!      IF (DUMSIGL.EQV..TRUE.) ISIGMA=DBLE(DUMSP)
!      CALL SPLIT_READ_LINE_STRING(LINE, 'DATAPATH', DATAPATH)
!      CALL SPLIT_READ_LINE_STRING(LINE, 'MODPATH', MODLPATH)
!      CALL SPLIT_READ_LINE_STRING(LINE, 'OUTPATH', OUTPPATH)
!      CALL SPLIT_READ_LINE_STRING(LINE, 'SLDPATH', LINEPATH)
!      NUM_LINES_OBS=NUM_LINES_OBS+1
!      !
!    ENDDO
!    CLOSE(1)
!    !
!    LINES_READ=LINES_READ+NUM_LINES_OBS
!    !
!    !IF (ABS(RSIGMAP).GT.1.0D-5) ISIGMA=DBLE(RSIGMAP)
!    IF (mpi__myrank.EQ.0) THEN
!IF (MVERBOSE.GT.0) THEN
!      PRINT*, 'Assist T: ', ASSIST_T
!      PRINT*, 'Assist P: ', ASSIST_P
!      PRINT*, 'Assist B: ', ASSIST_B
!      PRINT*, 'Assist V: ', ASSIST_V
!      PRINT*, 'DATAPATH: ', TRIM(DATAPATH)
!      PRINT*, 'MODPATH: ', TRIM(MODLPATH)
!      PRINT*, 'OUTPATH: ', TRIM(OUTPPATH)
!      PRINT*, 'SLDPATH: ', TRIM(LINEPATH)
!ENDIF
!    ENDIF
!    !
!    ! Check SVDTOL and NOISE ARE PROPERLY SUPPLIED:
!!><    IF (DUMSVDL.EQV..TRUE.) THEN
!!><      IF (mpi__myrank.EQ.0) THEN
!!><        PRINT*, ' *** '
!!><        PRINT*, ' SVDTOL in "MISC SETUP:" field is not properly supplied'
!!><        PRINT*, ' --- '
!!><      ENDIF
!!><      CALL MPI_BARRIER(MPI_COMM_WORLD, mpi__ierror)
!!><      STOP
!!><    ENDIF ! SVDTOL
!!><    IF (DUMSIGL.EQV..TRUE.) THEN
!!><PRINT*, DUMSIGL, DUMSP, DUMIP, RSIGMAP
!!><      IF (mpi__myrank.EQ.0) THEN
!!><        PRINT*, ' *** '
!!><        PRINT*, ' NOISE in "MISC SETUP:" field is not properly supplied'
!!><        PRINT*, ' --- '
!!><      ENDIF
!!><      CALL MPI_BARRIER(MPI_COMM_WORLD, mpi__ierror)
!!><      STOP
!!><    ENDIF ! NOISE
!!><CALL MPI_BARRIER(MPI_COMM_WORLD,mpi__ierror)
!!><STOP
!    !
!  END SUBROUTINE SET_MISC
  !
  !________________________________________________
  !
  SUBROUTINE CLOSE_READ_INPUT_FILE()
    !
    USE FORWARD_PARAM, ONLY: IND_LINE, WAVE_INI, DELTA_LAMBDA, NUMWAVE
    !
    DEALLOCATE(IND_LINE)
    DEALLOCATE(WAVE_INI)
    DEALLOCATE(DELTA_LAMBDA)
    DEALLOCATE(NUMWAVE)
    !
  END SUBROUTINE CLOSE_READ_INPUT_FILE
  !
  !________________________________________________
  !
  SUBROUTINE UPDATE_FREE_VAR(LOGIC, VMAX, VABS)
    !
    USE INVERT_PARAM, ONLY: NFREEV, NFREEP, TFREEP, INFREEPMAX
    !
    LOGICAL, INTENT(IN)    :: LOGIC
    INTEGER, INTENT(IN)    :: VMAX, VABS
    !
    IF (LOGIC.EQV..TRUE.) THEN
      NFREEV=NFREEV+1
      NFREEP=NFREEP+VMAX
      TFREEP=TFREEP+VABS
      IF (VMAX.GT.INFREEPMAX) INFREEPMAX=VMAX
    ENDIF
    !
  END SUBROUTINE UPDATE_FREE_VAR
  !
  !________________________________________________
  !
  SUBROUTINE FREE_VAR_SPACE(NZ)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP
    USE INVERT_PARAM, ONLY: INV_ATMPAR, NSLB_MAX, NFREEV, NFREEP &
        , TFREEP, INFREEPMAX, NSLAB_PER_FREEV, REG_TYPE_FREEV
    !
    INTEGER, INTENT(IN) :: NZ
    !
    INTEGER             :: I
    !
    NFREEV=0
    NFREEP=0
    TFREEP=0
    INFREEPMAX=0
    DO I=1,SIZE(INV_ATMPAR)
      IF (I.EQ.8) THEN
        CALL UPDATE_FREE_VAR(INV_ATMPAR(I), NSLB_MAX(I), 1)
      ELSE
        CALL UPDATE_FREE_VAR(INV_ATMPAR(I), NSLB_MAX(I), NZ)
      ENDIF
    ENDDO
    !
    CALL ALLOCATE_1D_IP(NSLAB_PER_FREEV, NFREEV, 'initial_duties: NSLAB_PER_FREEV')
    CALL ALLOCATE_1D_IP(REG_TYPE_FREEV, NFREEV, 'initial_duties: REG_TYPE_FREEV')
    REG_TYPE_FREEV(:)=1
    !
  END SUBROUTINE FREE_VAR_SPACE
   SUBROUTINE NSET_NLTE(START_LINE)
    USE FORWARD_PARAM, ONLY: NUML, BLENDSID, BLENDSMAX, POP_FILE, IND_LINE &
         , LINE_NUM, POPL, POPU, IND_FOUND
    USE CODE_MODES, ONLY: INPUTFILE, MVERBOSE, DATAPATH
    USE user_mpi, ONLY: mpi__myrank, MPI_COMM_WORLD
    USE LINES_DATABASE, ONLY: NCOLUMNS
    USE GRID_PARAM, ONLY: NZ
    !
    INTEGER, INTENT(INOUT)                    :: START_LINE
    INTEGER                                   :: I, J, K, M, NUM_LINES_OBS, NCOL, LINDEX, TAR
    LOGICAL                                   :: CONDIT, EXISTS
    CHARACTER*800                             :: LINE
    CHARACTER*100                             :: FILE
    INTEGER, ALLOCATABLE                      :: IND_FOUND_DUMMY(:)
    CHARACTER*800, ALLOCATABLE                :: POP_FILE_DUMMY(:)
    ALLOCATE(POP_FILE(NUML),IND_FOUND(NUML))
    ALLOCATE(POP_FILE_DUMMY(NUML),IND_FOUND_DUMMY(NUML))
    IND_FOUND(:)=-1
    DO J=1,NUML
       POP_FILE(J)='LTE'
    ENDDO
    !
    NUM_LINES_OBS=1
    I=0
    CONDIT=.TRUE.
    !
    DO WHILE (START_LINE < INPUT_READ_NUMBER)
       START_LINE=START_LINE+1
       LINE=INPUT_READ_LINES(START_LINE)
       ! Check if the section has finished:
       DO J=1,800
          IF (LINE(J:J) == ':') THEN
             CONDIT=.FALSE.
             START_LINE=START_LINE-1
             EXIT
          ENDIF
       ENDDO
       IF (.NOT.CONDIT) EXIT
       NCOL=NCOLUMNS(LINE)
       !PRINT*,NCOL,':',LINE
       SELECT CASE(NCOL)
          CASE(2)
             CALL SPLIT_NLTE_LINE_2C(LINE,LINDEX,FILE)
             POP_FILE(NUM_LINES_OBS)=FILE
             IND_FOUND(NUM_LINES_OBS)=LINDEX
             !PRINT*,FILE,LINDEX,NUM_LINES_OBS
         CASE DEFAULT 
             PRINT*, ''
             PRINT*, ' ** Wrong input format in: "NLTE DEPARTURE COEFFICIENTS"'
             PRINT*, ' -- Line: "'//TRIM(ADJUSTL(LINE))//'"!...'
             PRINT*, ' ...does not have 2 fields.'
             PRINT*, ''
             PRINT*, ' Stopping'
             PRINT*, ''
             STOP
       END SELECT
       ! Checking if file exists
       IF (TRIM(FILE) /= 'LTE') THEN
          INQUIRE(FILE=TRIM(DATAPATH)//"/./"//TRIM(FILE), EXIST=EXISTS)
          IF (.NOT.EXISTS) THEN
             PRINT*, ""
             PRINT*, "Error!"
             PRINT*, "I cannot find NLTE file: "//TRIM(FILE)//" in directory: "//TRIM(DATAPATH)
             PRINT*, ""
             STOP
             !CALL MPI_ABORT(MPI_COMM_WORLD)
          ENDIF ! Writing error
       ENDIF
       NUM_LINES_OBS=NUM_LINES_OBS+1
    ENDDO
    ! Now we need to check that NLTE lines are not blended
    DO I=1,NUML
       IF (IND_FOUND(I) /= -1) THEN
          DO J=1,NUML
             DO K=1,BLENDSMAX+1
                IF (BLENDSID(J,K) /= -1) THEN
                   IF (LINE_NUM(BLENDSID(J,K)) == IND_FOUND(I)) THEN
                      !PRINT*,IND_FOUND(I),LINE_NUM(BLENDSID(J,K))
                      !PRINT*,BLENDSID(J,:)
                      DO M=1,BLENDSMAX+1
                         IF (LINE_NUM(BLENDSID(J,M)) /= IND_FOUND(I)) THEN
                            IF (BLENDSID(J,M) /= -1) THEN
                               IF (MPI__MYRANK ==0) THEN
                                  PRINT*,'NLTE lines cannot be blended. Line index:',IND_FOUND(I)
                                  PRINT*,'They must appear as isolated entried inside LINES:'
                                  PRINT*,'Modify '//TRIM(INPUTFILE)
                               ENDIF
                               !PRINT*,BLENDSID(J,:)
                               !CALL MPI_ABORT(MPI_COMM_WORLD)
                               STOP
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! Now we check that all lines under NLTE DEPARTURE COEFFICIENTS:
    ! are also included in LINES:
    DO I=1,NUML
       IF (IND_FOUND(I) /= -1) THEN
          TAR=0
          DO J=1,NUML
             IF (IND_FOUND(I) == LINE_NUM(BLENDSID(J,1))) THEN
                TAR=TAR+1
                EXIT
             ENDIF
          ENDDO
          IF (TAR == 0) THEN
             IF (MPI__MYRANK == 0) THEN
                WRITE(*,FMT='(A,1X,I3,1X,A)') 'Why do you provide departure coefficients for line #' &
                     ,IND_FOUND(I),' but this line is not under LINES?'
                STOP
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    ! Now we reorder POP_FILE and IND_FOUND so that it has the same
    ! order as in LINES:
    IND_FOUND_DUMMY(:)=-1
    POP_FILE_DUMMY='LTE'
    DO J=1,NUML
       IF (BLENDSID(J,1) /= -1) THEN
          DO K=1,NUML
             IF (IND_FOUND(K) /= -1) THEN
                IF (IND_FOUND(K) == LINE_NUM(BLENDSID(J,1))) THEN
                   IND_FOUND_DUMMY(J)=IND_FOUND(K)
                   POP_FILE_DUMMY(J)=POP_FILE(K)
                   EXIT
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    ! After this IND_FOUND_DUMMY and POP_FILE_DUMMY will
    ! have the exact same order as BLENDSID
    POP_FILE=POP_FILE_DUMMY
    IND_FOUND=IND_FOUND_DUMMY
    DEALLOCATE(POP_FILE_DUMMY,IND_FOUND_DUMMY)
    ! Now we read departure coefficients
    DO I=1,NUML
       IF (IND_FOUND(I) /= -1) THEN
          IF (MPI__MYRANK == 0) THEN
             WRITE(*,FMT='(A,1X,I3,1X,A)') 'Spectral line #',IND_FOUND(I) &
                  ,'will be treated under NLTE using departure coefficients from ' &
                  ,TRIM(POP_FILE(I))
          ENDIF
       ELSE
          DO J=1,BLENDSMAX+1
             IF (MPI__MYRANK == 0) THEN
                IF (BLENDSID(I,J) /= -1) WRITE(*,FMT='(A,1X,I3,1X,A)') 'Spectral line ' &
                     ,LINE_NUM(BLENDSID(I,J)),' will be treated under LTE'
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    !
    DO I=1,NUML
       IF (IND_FOUND(I) /= -1) CALL READ_DEPARTURE_COEFFICIENTS(I)
    ENDDO
    !
  END SUBROUTINE NSET_NLTE
  !
  !_____________________________________________
  !
  SUBROUTINE READ_DEPARTURE_COEFFICIENTS(LINDEX)
    !
    USE GRID_PARAM, ONLY : NZ
    USE FORWARD_PARAM, ONLY: NUML, POP_FILE, POPL, POPU, IND_LINE, LINE_NUM, IND_FOUND
    USE CODE_MODES, ONLY: DATAPATH
    !
    INTEGER, INTENT(IN)        :: LINDEX
    REAL(DP)                   :: BETAL, BETAU
    INTEGER                    :: I
    !
    OPEN(UNIT=99,FILE=TRIM(DATAPATH)//TRIM(POP_FILE(LINDEX)),STATUS='OLD',FORM='UNFORMATTED')
    DO I=1,NZ
       READ(UNIT=99) BETAL, BETAU
       POPL(LINDEX,I)=BETAL
       POPU(LINDEX,I)=BETAU
    ENDDO
    CLOSE(UNIT=99)
    !IF (mpi__myrank.EQ.0) WRITE(*,FMT='(A,I3,A)') 'LINE #',IND_FOUND(LINDEX) &
    !     ,' will be treated under NLTE:',TRIM(POP_FILE(LINDEX))
    !
  END SUBROUTINE READ_DEPARTURE_COEFFICIENTS
  !________________________________________________


  !________________________________________________
  !
  SUBROUTINE SPLIT_NLTE_LINE_2C(INLINE, LINDEX, FILE)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP
    USE LINES_DATABASE, ONLY: NCOLUMNS
    !
    CHARACTER(*), INTENT(IN)       :: INLINE
    INTEGER                        :: LINDEX
    CHARACTER*100                  :: FILE
    REAL(DP)                       :: LSFSIGMA, LSFW0
    !
    INTEGER                         :: I, J, NCOL, ENDENTRY, IISIGMA
    !
    LOGICAL                         :: PRECHAR, SPECIALCHAR, ISREAL
    !
    INTEGER, ALLOCATABLE            :: INIT(:), FINI(:)
    !
    NCOL=NCOLUMNS(INLINE)
    !
    ! WHERE DO COLUMNS BEGIN?:
    CALL ALLOCATE_1D_IP(INIT,NCOL,'INIT ARRAY IN READ_LINES.')
    CALL ALLOCATE_1D_IP(FINI,NCOL,'FINI ARRAY IN READ_LINES.')
    !
    NCOL=0
    !
    PRECHAR=.FALSE.
    SPECIALCHAR=.FALSE.
    !
    DO I=1,LEN(INLINE)
       IF (INLINE(I:I).NE.' ') THEN
          ! AM I SPECIAL?
          IF (INLINE(I:I).EQ.'#') THEN
             SPECIALCHAR=.TRUE.
             CYCLE
          ENDIF
          !WAS THE PREVIOUS THING A NON EMPTY SPACE?
          IF ((PRECHAR.EQV..FALSE.).AND.(SPECIALCHAR.EQV..FALSE.)) THEN
             NCOL=NCOL+1
             INIT(NCOL)=I
             IF (NCOL.GT.1) FINI(NCOL-1)=I-1
          ELSE
             SPECIALCHAR=.FALSE.
          ENDIF
          ! UPDATE PRECHAR:
          PRECHAR=.TRUE.
       ELSE
          PRECHAR=.FALSE.
       ENDIF
       IF (SPECIALCHAR.EQV..TRUE.) EXIT
    ENDDO
    IF (SPECIALCHAR.EQV..TRUE.) THEN
       FINI(NCOL)=I-2
    ELSE
       FINI(NCOL)=I-1!LEN(INLINE)
    ENDIF
    !
    IF (NCOL.NE.2) THEN
       IF (mpi__myrank.EQ.0) THEN
          PRINT*, '   ***   '
          PRINT*, ' Wrong format! '
          PRINT*, INLINE
          PRINT*, '   ***   '
          STOP
       ENDIF
    ENDIF
    !
    ! NOW, WE ASSIGN THE VARIOUS THINGS:
    DO I=1,NCOL
       !PRINT*,INIT(I),FINI(I)
       !PRINT*,INLINE(INIT(I):FINI(I))
       SELECT CASE(I)
       CASE(1)
          READ(INLINE(INIT(I):FINI(I)), '(i14)') LINDEX
          !print*,lindex,i
       CASE(2)
          READ(INLINE(INIT(I):FINI(I)), '(A)') FILE
          !print*,TRIM(file),i
       END SELECT
       !READ(*,*)
    ENDDO

!PRINT*, 'INIT: ', INIT
!PRINT*, 'FINI: ', FINI
!PRINT*, INW, ISTART, ISTEP
    !
    DEALLOCATE(INIT)
    DEALLOCATE(FINI)

    !
  END SUBROUTINE SPLIT_NLTE_LINE_2C
  !________________________________________________
  !
  !================================================
  !
END MODULE INPUT_FILE
!
