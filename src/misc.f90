!
MODULE MISC
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP
  USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_SP,ALLOCATE_3D_SP &
      ,ALLOCATE_L1D_SP,ALLOCATE_5D_SP,ALLOCATE_4D_SP
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, PARAMETER :: LI=KIND(10000000)
  !INTEGER(LI), parameter :: maxrec=327155710
  INTEGER(LI), parameter :: bmaxrec=1073741784
  INTEGER(LI), parameter :: nmaxrec=bmaxrec/4
  !
  PUBLIC :: READ_SP1D
  PUBLIC :: READ_MODEL_3D
  PUBLIC :: READ_MODEL
  PRIVATE :: WRITE_MODEL3D_V1
  PRIVATE :: WRITE_MODEL3D_V2
  PRIVATE :: WRITE_MODEL3D_V3
  PUBLIC :: WRITE_MODEL
  PRIVATE :: WRITE_PROFILES3D_V1
  PRIVATE :: WRITE_PROFILES3D_V2
  PRIVATE :: WRITE_PROFILES3D_V3
  PUBLIC :: WRITE_PROFILES
  PUBLIC :: READ_PROFILES
  PUBLIC :: READ_PROFILE_3D
  PRIVATE :: WRITE_RFS3D_V1
  PRIVATE :: WRITE_RFS3D_V2
  PRIVATE :: WRITE_RFS3D_V3
  PUBLIC :: WRITE_LOGTAUS
  PRIVATE :: WRITE_LOGTAUS_V3
  PUBLIC :: WRITE_RFS
  PUBLIC :: WRITE_RFS2
  PUBLIC :: READ_PSF_2D
  PUBLIC :: WRITE_BIN
  PUBLIC :: READ_LSF_FILE
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! read_sp1d
  ! read_model_3d
  ! read_model
  ! read_profiles
  ! read_profile_3d
  ! read_psf_2d
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_SP1D (FNAME, SZ, TOREAD)
    !
    CHARACTER(*), INTENT(IN)              :: FNAME
    INTEGER,DIMENSION(1),INTENT(IN)       :: SZ
    REAL(SP),DIMENSION(SZ(1)),INTENT(INOUT)  :: TOREAD
    INTEGER :: IERR
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="UNFORMATTED",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      READ(1) TOREAD
    CLOSE(UNIT=1)
    !
  END SUBROUTINE READ_SP1D
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_MODEL_3D (FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,VV)
    !
    CHARACTER(*), INTENT(IN)                       :: FNAME
    INTEGER, INTENT(IN)                            :: NX, NY, NZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: TEM,PG,RHO,BX,BY,BZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: VZ,PEL,MW,TAU
    REAL(DP), INTENT(INOUT), DIMENSION(NX)   :: X
    REAL(DP), INTENT(INOUT), DIMENSION(NY)   :: Y
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)   :: Z
    INTEGER, INTENT(IN)                      :: VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOREAD
    !
    INTEGER                                  :: NTOT
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_READ
    !
    INTEGER                                  :: RFID
    INTEGER                                  :: NSDIMS
    INTEGER                                  :: NPAR, NIX, NIY, NIZ, RNREC
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER                                  :: I
    INTEGER :: IERR
    !
    IF (VV.EQ.1) THEN
      ! READ THE FIRST 5 ELEMENTS:
      !
      NTOT=5
      !
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      ! VERSION: ALREADY CHECKED!
      ! ID
      RFID=INT(TOREAD(3))
      IF ( ABS(RFID-130904).LT.0.01) THEN
        WRITE(*,*) 'READING MODEL ATMOSPHERE FILE: '//TRIM(FNAME)
      ELSE
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      ! NUMBER OF RECORDS:
      RNREC=INT(TOREAD(4))
      ! NDIMS:
      NSDIMS=INT(TOREAD(5))
      IF (NSDIMS.NE.4) THEN
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      DEALLOCATE(TOREAD)
      !
      NTOT=5+NSDIMS
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      NPAR=INT(TOREAD(6))
      NIX=INT(TOREAD(7))
      NIY=INT(TOREAD(8))
      NIZ=INT(TOREAD(9))
      !
      DEALLOCATE(TOREAD)
      !
      IF ((NPAR.NE.13).OR.(NIX.NE.NX).OR.(NIY.NE.NY).OR.(NIZ.NE.NZ)) THEN
        WRITE(*,*) 'EITHER THERE ARE NOT ENOUGH ATMOSPHERE PARAMETERS'
        WRITE(*,*) 'OR THE NUMBER OF Z GRIDS IN THE FILE '//TRIM(FNAME)
        WRITE(*,*) 'AND THE ONE FROM THE INPUTFILE DO NOT MATCH OR BOTH'
        STOP
      ENDIF
      !
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA:
      CALL ALLOCATE_3D_SP(TMP,NIX,NIY,NIZ,'READ_ALLOCATION TMP')
      !
      OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
        DO I=1,RNREC
          IF (I.EQ.1) THEN
            NTOT=5+NSDIMS
          ELSE
            NTOT=INT(NIX*NIY*NIZ)
          ENDIF
          ALLOCATE(TOREAD(NTOT))
          READ(1) TOREAD
          IF (I.NE.1) THEN
            SELECT CASE (I)
              CASE(2)
                TEM(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(3)
                PG(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(4)
                RHO(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(5)
                BX(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(6)
                BY(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(7)
                BZ(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(8)
                VZ(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(9)
                PEL(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(10)
                MW(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(11)
                TAU(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
              CASE(12)
                TMP(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
                X(:)=DBLE(TMP(:,1,1))
              CASE(13)
                TMP(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
                Y(:)=DBLE(TMP(1,:,1))
              CASE(14)
                TMP(:,:,:)=RESHAPE(TOREAD,(/NIX,NIY,NIZ/))
                Z(:)=DBLE(TMP(1,1,:))
            ENDSELECT
          ENDIF
          DEALLOCATE(TOREAD)
        ENDDO
      CLOSE(UNIT=1)
      DEALLOCATE(TMP)
      ! END VERSION 1
!--------------------------------------------------------------------
    ELSE IF (VV.EQ.2) THEN
      ! READ THE FIRST 5 ELEMENTS:
      !
      NTOT=5
      !
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      ! VERSION: ALREADY CHECKED!
      ! ID
      RFID=INT(TOREAD(3))
      IF ( ABS(RFID-130904).LT.0.01) THEN
        WRITE(*,*) 'READING MODEL ATMOSPHERE FILE: '//TRIM(FNAME)
      ELSE
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      ! NUMBER OF RECORDS:
      RNREC=INT(TOREAD(4))
      ! NDIMS:
      NSDIMS=INT(TOREAD(5))
      IF (NSDIMS.NE.4) THEN
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      DEALLOCATE(TOREAD)
      !
      NTOT=5+NSDIMS
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      NPAR=INT(TOREAD(6))
      NIX=INT(TOREAD(7))
      NIY=INT(TOREAD(8))
      NIZ=INT(TOREAD(9))
      !
      DEALLOCATE(TOREAD)
      !
      IF ((NPAR.NE.13).OR.(NIX.NE.NX).OR.(NIY.NE.NY).OR.(NIZ.NE.NZ)) THEN
        WRITE(*,*) 'EITHER THERE ARE NOT ENOUGH ATMOSPHERE PARAMETERS'
        WRITE(*,*) 'OR THE NUMBER OF X,Y,AND/OR Z GRIDS IN THE FILE '
        WRITE(*,*) TRIM(FNAME)//'AND THE ONE FROM THE INPUTFILE DO NOT'
        WRITE(*,*) ' MATCH OR BOTH'
        WRITE(*,*) 'NX', NIX, NX
        WRITE(*,*) 'NY', NIY, NY
        WRITE(*,*) 'NZ', NIZ, NZ
        WRITE(*,*) 'NPAR', NPAR, 13
        STOP
      ENDIF
      !
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA AND STORE IT:
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIZ*NPAR
      CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
      !
      TO_BE_READ=LNTOT
      LOFFSET=1
      !
      OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
        DO I=1,RNREC
          IF (I.EQ.1) THEN
            LNTOT=5
            LNTOT=LNTOT+NSDIMS
          ELSE
            LNTOT=1
            IF (TO_BE_READ.GT.nmaxrec) THEN 
              LNTOT=nmaxrec
            ELSE
              LNTOT=TO_BE_READ
            ENDIF
          ENDIF
          ALLOCATE(TOREAD(LNTOT))
          READ(1) TOREAD
          IF (I.NE.1) THEN
!PRINT*, 'PRE:', SIZE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SIZE(TOREAD)
            TMP1D(LOFFSET:LOFFSET+LNTOT-1)=TOREAD
            LOFFSET=LOFFSET+LNTOT
            TO_BE_READ=TO_BE_READ-LNTOT
          ENDIF
          DEALLOCATE(TOREAD)
        ENDDO
      CLOSE(UNIT=1)
      ! SAVE IT IN THE ATMOSPHERE:
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA:
      CALL ALLOCATE_3D_SP(TMP,NIX,NIY,NIZ,'READ_ALLOCATION TMP')
      !
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIZ
      LOFFSET=1
      !
      DO I=1,NPAR
        TMP(:,:,:)=0.E0
        TMP(:,:,:)=RESHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1),(/NIX,NIY,NIZ/))
        !PRINT*, LOFFSET,LOFFSET+LNTOT,SHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SHAPE(TMP)&
        !    , SUM(TMP/REAL(LNTOT)), I
        SELECT CASE (I)
          CASE(1)
            TEM(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(2)
            PG(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(3)
            RHO(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(4)
            BX(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(5)
            BY(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(6)
            BZ(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(7)
            VZ(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(8)
            PEL(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(9)
            MW(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(10)
            TAU(:,:,:)=RESHAPE(TMP, (/ NZ,NY,NX /),ORDER=(/ 3,2,1 /))
          CASE(11)
            X(:)=DBLE(TMP(:,1,1))
          CASE(12)
            Y(:)=DBLE(TMP(1,:,1))
          CASE(13)
            Z(:)=DBLE(TMP(1,1,:))
        ENDSELECT
        LOFFSET=LOFFSET+LNTOT
      ENDDO
      !
      DEALLOCATE(TMP)
      DEALLOCATE(TMP1D)
!PRINT*, 'INTERNAL READ Z: ', z(1), Z(NIZ)
!PRINT*, 'INTERNAL READ Y: ', Y(1), Y(NIY)
!PRINT*, 'INTERNAL READ X: ', X(1), X(NIX)
      !END VERSION 2
!--------------------------------------------------------------------
    ELSE IF (VV.EQ.3) THEN
      ! READ THE FIRST 5 ELEMENTS:
      !
      NTOT=5
      !
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      ! VERSION: ALREADY CHECKED!
      ! ID
      RFID=INT(TOREAD(3))
      IF ( ABS(RFID-130904).LT.0.01) THEN
        WRITE(*,*) 'READING MODEL ATMOSPHERE FILE: '//TRIM(FNAME)
      ELSE
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      ! NUMBER OF RECORDS:
      RNREC=INT(TOREAD(4))
      ! NDIMS:
      NSDIMS=INT(TOREAD(5))
      IF (NSDIMS.NE.4) THEN
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      DEALLOCATE(TOREAD)
      !
      NTOT=5+NSDIMS
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      NPAR=INT(TOREAD(6))
      NIX=INT(TOREAD(7))
      NIY=INT(TOREAD(8))
      NIZ=INT(TOREAD(9))
      !
      DEALLOCATE(TOREAD)
      !
      IF ((NPAR.NE.13).OR.(NIX.NE.NX).OR.(NIY.NE.NY).OR.(NIZ.NE.NZ)) THEN
        WRITE(*,*) 'EITHER THERE ARE NOT ENOUGH ATMOSPHERE PARAMETERS'
        WRITE(*,*) 'OR THE NUMBER OF X,Y,AND/OR Z GRIDS IN THE FILE '
        WRITE(*,*) TRIM(FNAME)//'AND THE ONE FROM THE INPUTFILE DO NOT'
        WRITE(*,*) ' MATCH OR BOTH'
        WRITE(*,*) 'NX', NIX, NX
        WRITE(*,*) 'NY', NIY, NY
        WRITE(*,*) 'NZ', NIZ, NZ
        WRITE(*,*) 'NPAR', NPAR, 13
        STOP
      ENDIF
      !
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA AND STORE IT:
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIZ*NPAR
      CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
      !
      TO_BE_READ=LNTOT
      LOFFSET=1
!PRINT*, ' -> Reading model: ', nmaxrec
      !
      OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
        DO I=1,RNREC
          IF (I.EQ.1) THEN
            LNTOT=5
            LNTOT=LNTOT+NSDIMS
          ELSE
            LNTOT=1
            IF (TO_BE_READ.GT.nmaxrec) THEN 
              LNTOT=nmaxrec
            ELSE
              LNTOT=TO_BE_READ
            ENDIF
          ENDIF
          ALLOCATE(TOREAD(LNTOT))
          READ(1) TOREAD
          IF (I.NE.1) THEN
!PRINT*, 'PRE:', SIZE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SIZE(TOREAD)
            TMP1D(LOFFSET:LOFFSET+LNTOT-1)=TOREAD
            LOFFSET=LOFFSET+LNTOT
            TO_BE_READ=TO_BE_READ-LNTOT
          ENDIF
          DEALLOCATE(TOREAD)
        ENDDO
      CLOSE(UNIT=1)
      ! SAVE IT IN THE ATMOSPHERE:
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA:
      CALL ALLOCATE_3D_SP(TMP,NIZ,NIY,NIX,'READ_ALLOCATION TMP')
      !
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIZ
      LOFFSET=1
      !
      DO I=1,NPAR
        TMP(:,:,:)=0.E0
        TMP(:,:,:)=RESHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1),(/NIZ,NIY,NIX/))
        !PRINT*, I, SUM(TMP)/REAL(LNTOT)
        SELECT CASE (I)
          CASE(1)
            TEM(:,:,:)=TMP
          CASE(2)
            PG(:,:,:)=TMP
          CASE(3)
            RHO(:,:,:)=TMP
          CASE(4)
            BX(:,:,:)=TMP
          CASE(5)
            BY(:,:,:)=TMP
          CASE(6)
            BZ(:,:,:)=TMP
          CASE(7)
            VZ(:,:,:)=TMP
          CASE(8)
            PEL(:,:,:)=TMP
          CASE(9)
            MW(:,:,:)=TMP
          CASE(10)
            TAU(:,:,:)=TMP
          CASE(11)
            X(:)=DBLE(TMP(1,1,:))
          CASE(12)
            Y(:)=DBLE(TMP(1,:,1))
          CASE(13)
            Z(:)=DBLE(TMP(:,1,1))
        ENDSELECT
        LOFFSET=LOFFSET+LNTOT
      ENDDO
      !
      DEALLOCATE(TMP)
      DEALLOCATE(TMP1D)
!PRINT*, 'INTERNAL READ Z: ', z(1), Z(NIZ)
!PRINT*, 'INTERNAL READ Y: ', Y(1), Y(NIY)
!PRINT*, 'INTERNAL READ X: ', X(1), X(NIX)
      !END VERSION 3
    ENDIF
    !
  END SUBROUTINE READ_MODEL_3D
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_MODEL (FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z)
    !
    CHARACTER(*), INTENT(IN)                       :: FNAME
    INTEGER, INTENT(IN)                            :: NX, NY, NZ
    REAL(SP), INTENT(INOUT), DIMENSION(NX,NY,NZ)   :: TEM,PG,RHO,BX,BY,BZ
    REAL(SP), INTENT(INOUT), DIMENSION(NX,NY,NZ)   :: VZ,PEL,MW,TAU
    REAL(DP), INTENT(INOUT), DIMENSION(NX)   :: X
    REAL(DP), INTENT(INOUT), DIMENSION(NY)   :: Y
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)   :: Z
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOREAD
    !
    INTEGER                                  :: NTOT
    !
    INTEGER                                  :: POSV, NEGV, MEDV
    !
    ! READ THE FIRST 4 ELEMENTS:
    !
    NTOT=2
    !
    CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
    !
    CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
    !
    ! VERSION:
    POSV=INT(TOREAD(1))
    NEGV=INT(TOREAD(2))
    !
    MEDV=(POSV+NEGV)/2
    !
    POSV=POSV-MEDV
    NEGV=NEGV-MEDV
    !
    IF ( (POSV+NEGV).GT.0.01) THEN
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'NON RECOGNIZED FORMAT. STOPPING'
      STOP
    ELSE
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'VERSION ', (POSV - NEGV) / 2
    ENDIF
    !
    IF (MEDV.EQ.3000) THEN
      CALL READ_MODEL_3D(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,(POSV - NEGV) / 2)
    ENDIF
    !
  END SUBROUTINE READ_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_MODEL3D_V1(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                       :: FNAME
    INTEGER, INTENT(IN)                            :: NX, NY, NZ
    REAL(SP), INTENT(INOUT), DIMENSION(NX,NY,NZ)   :: TEM,PG,RHO,BX,BY,BZ
    REAL(SP), INTENT(INOUT), DIMENSION(NX,NY,NZ)   :: VZ,PEL,MW,TAU
    REAL(DP), INTENT(INOUT), DIMENSION(NX)   :: X
    REAL(DP), INTENT(INOUT), DIMENSION(NY)   :: Y
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)   :: Z
    INTEGER, INTENT(IN)                      :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOWRITE
    !
    INTEGER                                  :: NTOT
    INTEGER                                  :: I, J, K
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    INTEGER :: IERR
    !
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_MODEL3D_V1 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=130904.E0
    TOWRITE(4)=14.E0
    TOWRITE(5)=4.E0
    TOWRITE(6)=13.E0
    TOWRITE(7)=REAL(NX)
    TOWRITE(8)=REAL(NY)
    TOWRITE(9)=REAL(NZ)
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error writing: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      NTOT=INT(NX*NY*NZ)
      CALL ALLOCATE_3D_SP(TMP,NX,NY,NZ,'WRITE_MODEL3D_V1 TMP')
      DO I=1,13
        ALLOCATE(TOWRITE(NTOT))
        SELECT CASE (I)
          CASE(1)
            TOWRITE(:)=RESHAPE(TEM(:,:,:),(/NTOT/))
          CASE(2)
            TOWRITE(:)=RESHAPE(PG(:,:,:),(/NTOT/))
          CASE(3)
            TOWRITE(:)=RESHAPE(RHO(:,:,:),(/NTOT/))
          CASE(4)
            TOWRITE(:)=RESHAPE(BX(:,:,:),(/NTOT/))
          CASE(5)
            TOWRITE(:)=RESHAPE(BY(:,:,:),(/NTOT/))
          CASE(6)
            TOWRITE(:)=RESHAPE(BZ(:,:,:),(/NTOT/))
          CASE(7)
            TOWRITE(:)=RESHAPE(VZ(:,:,:),(/NTOT/))
          CASE(8)
            TOWRITE(:)=RESHAPE(PEL(:,:,:),(/NTOT/))
          CASE(9)
            TOWRITE(:)=RESHAPE(MW(:,:,:),(/NTOT/))
          CASE(10)
            TOWRITE(:)=RESHAPE(TAU(:,:,:),(/NTOT/))
          CASE(11)
            DO J=1,NY
              DO K=1,NZ
                TMP(:,J,K)=REAL(X(:))
              ENDDO
            ENDDO
            TOWRITE=RESHAPE(TMP(:,:,:),(/NTOT/))
          CASE(12)
            DO J=1,NX
              DO K=1,NZ
                TMP(J,:,K)=REAL(Y(:))
              ENDDO
            ENDDO
            TOWRITE=RESHAPE(TMP(:,:,:),(/NTOT/))
          CASE(13)
            DO J=1,NX
              DO K=1,NY
                TMP(J,K,:)=REAL(Z(:))
              ENDDO
            ENDDO
            TOWRITE=RESHAPE(TMP(:,:,:),(/NTOT/))
        ENDSELECT
        WRITE(1) TOWRITE
        DEALLOCATE(TOWRITE)
      ENDDO
    CLOSE(UNIT=1)
    !
  END SUBROUTINE WRITE_MODEL3D_V1
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_MODEL3D_V2(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                       :: FNAME
    INTEGER, INTENT(IN)                            :: NX, NY, NZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: TEM,PG,RHO,BX,BY,BZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: VZ,PEL,MW,TAU
    REAL(DP), INTENT(INOUT), DIMENSION(NX)   :: X
    REAL(DP), INTENT(INOUT), DIMENSION(NY)   :: Y
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)   :: Z
    INTEGER, INTENT(IN)                      :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOWRITE
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    INTEGER                                  :: NTOT, NPAR, NRECS
    INTEGER                                  :: I, J, K
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER :: IERR
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_MODEL3D_V2 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=130904.E0
    !
    NPAR=13
    LNTOT=1
    LNTOT=LNTOT*NX*NY*NZ*NPAR
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    !
    TOWRITE(4)=NRECS+1.E0   ! NREC
    TOWRITE(5)=4.E0
    TOWRITE(6)=REAL(NPAR)
    TOWRITE(7)=REAL(NX)
    TOWRITE(8)=REAL(NY)
    TOWRITE(9)=REAL(NZ)
    !
    ! STORE ATMOSPHERE IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
    !
    NTOT=INT(NX*NY*NZ)
    CALL ALLOCATE_3D_SP(TMP,NX,NY,NZ,'WRITE_MODEL3D_V2 TMP')
    !
    LOFFSET=1
    DO I=1,NPAR
      SELECT CASE (I)
        CASE(1)
          TMP=RESHAPE(TEM, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(2)
          TMP=RESHAPE(PG, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(3)
          TMP=RESHAPE(RHO, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(4)
          TMP=RESHAPE(BX, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(5)
          TMP=RESHAPE(BY, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(6)
          TMP=RESHAPE(BZ, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(7)
          TMP=RESHAPE(VZ, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(8)
          TMP=RESHAPE(PEL, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(9)
          TMP=RESHAPE(MW, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(10)
          TMP=RESHAPE(TAU, (/NX,NY,NZ/), ORDER=(/ 3,2,1 /))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(11)
          DO J=1,NY
            DO K=1,NZ
              TMP(:,J,K)=REAL(X(:))
            ENDDO
          ENDDO
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(12)
          DO J=1,NX
            DO K=1,NZ
              TMP(J,:,K)=REAL(Y(:))
            ENDDO
          ENDDO
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(13)
          DO J=1,NX
            DO K=1,NY
              TMP(J,K,:)=REAL(Z(:))
            ENDDO
          ENDDO
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
      ENDSELECT
      LOFFSET=LOFFSET+NTOT
    ENDDO
    DEALLOCATE(TMP)
    !
    !ONCE TMP1D IS FILLED, WE WRITE THE MODEL:
    !
    OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error writing: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN 
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    !
    DEALLOCATE(TMP1D)
    !
  END SUBROUTINE WRITE_MODEL3D_V2
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_MODEL3D_V3(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                       :: FNAME
    INTEGER, INTENT(IN)                            :: NX, NY, NZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: TEM,PG,RHO,BX,BY,BZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: VZ,PEL,MW,TAU
    REAL(DP), INTENT(INOUT), DIMENSION(NX)   :: X
    REAL(DP), INTENT(INOUT), DIMENSION(NY)   :: Y
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)   :: Z
    INTEGER, INTENT(IN)                      :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOWRITE
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    INTEGER                                  :: NTOT, NPAR, NRECS
    INTEGER                                  :: I, J, K
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER :: IERR
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_MODEL3D_V2 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=130904.E0
    !
    NPAR=13
    LNTOT=1
    LNTOT=LNTOT*NX*NY*NZ*NPAR
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    !
    TOWRITE(4)=NRECS+1.E0   ! NREC
    TOWRITE(5)=4.E0
    TOWRITE(6)=REAL(NPAR)
    TOWRITE(7)=REAL(NX)
    TOWRITE(8)=REAL(NY)
    TOWRITE(9)=REAL(NZ)
    !
    ! STORE ATMOSPHERE IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
    !
    NTOT=INT(NX*NY*NZ)
    CALL ALLOCATE_3D_SP(TMP,NZ,NY,NX,'WRITE_MODEL3D_V2 TMP')
    !
    LOFFSET=1
    DO I=1,NPAR
      SELECT CASE (I)
        CASE(1)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TEM(:,:,:),(/NTOT/))
        CASE(2)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(PG(:,:,:),(/NTOT/))
        CASE(3)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(RHO(:,:,:),(/NTOT/))
        CASE(4)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(BX(:,:,:),(/NTOT/))
        CASE(5)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(BY(:,:,:),(/NTOT/))
        CASE(6)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(BZ(:,:,:),(/NTOT/))
        CASE(7)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(VZ(:,:,:),(/NTOT/))
        CASE(8)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(PEL(:,:,:),(/NTOT/))
        CASE(9)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(MW(:,:,:),(/NTOT/))
        CASE(10)
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TAU(:,:,:),(/NTOT/))
        CASE(11)
          DO J=1,NY
            DO K=1,NZ
              TMP(K,J,:)=REAL(X(:))
            ENDDO
          ENDDO
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(12)
          DO J=1,NX
            DO K=1,NZ
              TMP(K,:,J)=REAL(Y(:))
            ENDDO
          ENDDO
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
        CASE(13)
          DO J=1,NX
            DO K=1,NY
              TMP(:,K,J)=REAL(Z(:))
            ENDDO
          ENDDO
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
      ENDSELECT
      LOFFSET=LOFFSET+NTOT
    ENDDO
    DEALLOCATE(TMP)
    !
    !ONCE TMP1D IS FILLED, WE WRITE THE MODEL:
    !
    OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error writing: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN 
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    !
    DEALLOCATE(TMP1D)
    !
  END SUBROUTINE WRITE_MODEL3D_V3
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_MODEL (FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                       :: FNAME
    INTEGER, INTENT(IN)                            :: NX, NY, NZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: TEM,PG,RHO,BX,BY,BZ
    REAL(SP), INTENT(INOUT), DIMENSION(NZ,NY,NX)   :: VZ,PEL,MW,TAU
    REAL(DP), INTENT(INOUT), DIMENSION(NX)   :: X
    REAL(DP), INTENT(INOUT), DIMENSION(NY)   :: Y
    REAL(DP), INTENT(INOUT), DIMENSION(NZ)   :: Z
    INTEGER, INTENT(IN)                      :: SS, VV
    !
    WRITE(*,*) 'WRITING ATMOSPHERE IN FILE: ' // TRIM(FNAME)
    !
    IF (SS.EQ.3000) THEN
      IF (VV.EQ.1) THEN
        CALL WRITE_MODEL3D_V1(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
      ELSE IF (VV.EQ.2) THEN
        CALL WRITE_MODEL3D_V2(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
      ELSE IF (VV.EQ.3) THEN
        CALL WRITE_MODEL3D_V3(FNAME,NX,NY,NZ,TEM,PG,RHO,BX,BY,BZ,VZ,PEL,MW,TAU,X,Y,Z,SS,VV)
      ENDIF
    ENDIF
    !
  END SUBROUTINE WRITE_MODEL
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_PROFILES3D_V1(FNAME,SZ,IND,WAVE,STOKES,SS,VV)

    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    REAL(DP), INTENT(IN), DIMENSION(SZ(3))                      :: IND
    !INTEGER, INTENT(IN), DIMENSION(SZ(3))                      :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(3))                     :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: STOKES
    !
    INTEGER, INTENT(IN)                                        :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                        :: TOWRITE
    !
    INTEGER                                                    :: NTOT, I
    INTEGER :: IERR
    !
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_PROFILES3D_V1 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=160904.E0
    TOWRITE(4)=7.E0
    TOWRITE(5)=4 ! NDIMS
    TOWRITE(6)=REAL(SZ(1))
    TOWRITE(7)=REAL(SZ(2))
    TOWRITE(8)=REAL(SZ(3))
    TOWRITE(9)=REAL(SZ(4))
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error writing: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      DO I=1,6
        SELECT CASE (I)
          CASE(1)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(IND)
          CASE(2)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(WAVE)
          CASE(3)
            NTOT=INT(SZ(1)*SZ(2)*SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=RESHAPE(STOKES(:,:,:,1),(/NTOT/))
          CASE(4)
            NTOT=INT(SZ(1)*SZ(2)*SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=RESHAPE(STOKES(:,:,:,2),(/NTOT/))
          CASE(5)
            NTOT=INT(SZ(1)*SZ(2)*SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=RESHAPE(STOKES(:,:,:,3),(/NTOT/))
          CASE(6)
            NTOT=INT(SZ(1)*SZ(2)*SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=RESHAPE(STOKES(:,:,:,4),(/NTOT/))
        ENDSELECT
        WRITE(1) TOWRITE
        DEALLOCATE(TOWRITE)
      ENDDO
    CLOSE(UNIT=1)
    !
  END SUBROUTINE WRITE_PROFILES3D_V1
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_PROFILES3D_V2(FNAME,SZ,IND,WAVE,STOKES,SS,VV)

    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    REAL(DP), INTENT(IN), DIMENSION(SZ(2))                      :: IND
    !INTEGER, INTENT(IN), DIMENSION(SZ(2))                      :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(2))                     :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: STOKES
    !
    INTEGER, INTENT(IN)                                        :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                        :: TOWRITE
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    INTEGER                                  :: NTOT, NPAR, NRECS
    INTEGER                                  :: I
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER :: IERR
    !
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_PROFILES3D_V1 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=160904.E0
    NPAR=4
    LNTOT=1
    LNTOT=LNTOT*SZ(1)*SZ(2)*SZ(3)*SZ(4)
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    !
    TOWRITE(4)=NRECS+3.E0   ! NREC
    TOWRITE(5)=4           ! NDIMS
    TOWRITE(6)=REAL(SZ(4)) !x
    TOWRITE(7)=REAL(SZ(3)) !y
    TOWRITE(8)=REAL(SZ(2)) !w
    TOWRITE(9)=REAL(SZ(1)) !Stokes
    !
    !
    ! STORE PROFILES IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP1D')
    CALL ALLOCATE_3D_SP(TMP,SZ(4),SZ(3),SZ(2),'READ_ALLOCATION TMP')
    !
    NTOT=INT(SZ(4)*SZ(3)*SZ(2))
    !
    LOFFSET=1
    DO I=1,NPAR
      TMP(:,:,:)=RESHAPE(STOKES(I,:,:,:), (/SZ(4), SZ(3), SZ(2)/), ORDER=(/3,2,1/))
      TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP(:,:,:),(/NTOT/))
      LOFFSET=LOFFSET+NTOT
    ENDDO
    !
    ! ONCE STORED, WE WRITE IT
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error writing: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! INDEX
      NTOT=INT(SZ(2))
      ALLOCATE(TOWRITE(NTOT))
      TOWRITE(:)=REAL(IND)
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! WAVELENGTH
      NTOT=INT(SZ(2))
      ALLOCATE(TOWRITE(NTOT))
      TOWRITE(:)=REAL(WAVE)
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    !
    DEALLOCATE(TMP)
    DEALLOCATE(TMP1D)
    !
  END SUBROUTINE WRITE_PROFILES3D_V2
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_PROFILES3D_V3(FNAME,SZ,IND,WAVE,STOKES,SS,VV)

    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    REAL(DP), INTENT(IN), DIMENSION(SZ(2))                      :: IND
    !INTEGER, INTENT(IN), DIMENSION(SZ(2))                      :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(2))                     :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: STOKES
    !
    INTEGER, INTENT(IN)                                        :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                        :: TOWRITE
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    INTEGER                                  :: NTOT, NPAR, NRECS
    INTEGER                                  :: I
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER :: IERR
    !
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_PROFILES3D_V1 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=160904.E0
    NPAR=4
    LNTOT=1
    LNTOT=LNTOT*SZ(1)*SZ(2)*SZ(3)*SZ(4)
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    !
    TOWRITE(4)=NRECS+3.E0   ! NREC
    TOWRITE(5)=4           ! NDIMS
    TOWRITE(6)=REAL(SZ(4)) !x
    TOWRITE(7)=REAL(SZ(3)) !y
    TOWRITE(8)=REAL(SZ(2)) !w
    TOWRITE(9)=REAL(SZ(1)) !Stokes
    !
    !
    ! STORE PROFILES IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP1D')
    !CALL ALLOCATE_3D_SP(TMP,SZ(4),SZ(3),SZ(2),'READ_ALLOCATION TMP')
    !
    NTOT=INT(SZ(1)*SZ(2)*SZ(3))
    !
    LOFFSET=1
    DO I=1,SZ(4) !x
      TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(STOKES(:,:,:,I),(/NTOT/))
      LOFFSET=LOFFSET+NTOT
    ENDDO
    !
    ! ONCE STORED, WE WRITE IT
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error writing: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! INDEX
      NTOT=INT(SZ(2))
      ALLOCATE(TOWRITE(NTOT))
      TOWRITE(:)=REAL(IND)
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! WAVELENGTH
      NTOT=INT(SZ(2))
      ALLOCATE(TOWRITE(NTOT))
      TOWRITE(:)=REAL(WAVE)
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    !
    !DEALLOCATE(TMP)
    DEALLOCATE(TMP1D)
    !
  END SUBROUTINE WRITE_PROFILES3D_V3
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_PROFILES(FNAME,SZ,IND,WAVE,STOKES,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    REAL(DP), INTENT(IN), DIMENSION(SZ(2))                      :: IND
    !INTEGER, INTENT(IN), DIMENSION(SZ(2))                      :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(2))                     :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: STOKES
    !
    INTEGER, INTENT(IN)                                        :: SS, VV
    !
    IF (SS.EQ.3000) THEN
      IF (VV.EQ.1) THEN
        CALL WRITE_PROFILES3D_V1(FNAME,SZ,IND,WAVE,STOKES,SS,VV)
      ELSE IF (VV.EQ.2) THEN
        CALL WRITE_PROFILES3D_V2(FNAME,SZ,IND,WAVE,STOKES,SS,VV)
      ELSE IF (VV.EQ.3) THEN
        CALL WRITE_PROFILES3D_V3(FNAME,SZ,IND,WAVE,STOKES,SS,VV)
      ENDIF
    ENDIF
    !
    WRITE(*,*) 'WRITING PROFILES IN FILE: ' // TRIM(FNAME)
    !
    !WRITE(*,*) 'I AM write_profiles, BUT I AM NOT IMPLEMENTED YET O_O'
    !
  END SUBROUTINE WRITE_PROFILES
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_PROFILES (FNAME,SZ,IND,WAVE,STOKES)
    !
    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    !INTEGER, INTENT(INOUT), DIMENSION(SZ(2))                      :: IND
    REAL(DP), INTENT(INOUT), DIMENSION(SZ(2))                      :: IND
    REAL(DP), INTENT(INOUT), DIMENSION(SZ(2))                     :: WAVE
    REAL(SP), INTENT(INOUT), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: STOKES
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOREAD
    !
    INTEGER                                  :: NTOT
    !
    INTEGER                                  :: POSV, NEGV, MEDV
    REAL(DP)               , DIMENSION(SZ(2))                      :: DIND
    !INTEGER               , DIMENSION(SZ(2))                      :: DIND
    REAL(DP)               , DIMENSION(SZ(2))                     :: DWAVE
    !
    ! READ THE FIRST 4 ELEMENTS:
    !
    NTOT=2
    !
    CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
    !
    CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
    !
    ! VERSION:
    POSV=INT(TOREAD(1))
    NEGV=INT(TOREAD(2))
    !
    MEDV=(POSV+NEGV)/2
    !
    POSV=POSV-MEDV
    NEGV=NEGV-MEDV
    !
    IF ( (POSV+NEGV).GT.0.01) THEN
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'NON RECOGNIZED FORMAT. STOPPING'
      STOP
    ELSE
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'VERSION ', (POSV - NEGV) / 2
    ENDIF
    !
    !WRITE(*,*) 'LRWNGFA'
    IF (MEDV.EQ.3000) THEN
      WRITE(*,*) 'I AM HERE!'
      CALL READ_PROFILE_3D(FNAME,SHAPE(STOKES),DIND,DWAVE,STOKES,(POSV - NEGV) / 2)
    ENDIF
    !
  END SUBROUTINE READ_PROFILES
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_PROFILE_3D (FNAME,SZ,IND,WAVE,STOKES,VV)
    !
    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    REAL(DP), INTENT(INOUT), DIMENSION(SZ(2))                      :: IND
    !INTEGER, INTENT(INOUT), DIMENSION(SZ(2))                      :: IND
    REAL(DP), INTENT(INOUT), DIMENSION(SZ(2))                     :: WAVE
    REAL(SP), INTENT(INOUT), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: STOKES
    !
    INTEGER, INTENT(IN)                      :: VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOREAD
    !
    INTEGER                                  :: NTOT
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_READ
    !
    INTEGER                                  :: RFID
    INTEGER                                  :: NSDIMS
    INTEGER                                  :: NIX, NIY, NIW, NIS, RNREC
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER                                  :: I
    INTEGER :: IERR
    !
    IF (VV.EQ.1) THEN
      WRITE(*,*) 'NOT IMPLEMENTED YET!'
      STOP
      ! END VERSION 1
!--------------------------------------------------------------------
    ELSE IF (VV.EQ.2) THEN
      ! READ THE FIRST 5 ELEMENTS:
      !
      NTOT=5
      !
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      ! VERSION: ALREADY CHECKED!
      ! ID
      RFID=INT(TOREAD(3))
      IF ( ABS(RFID-160904).LT.0.01) THEN
        PRINT*, 'READING STOKES PROFILES FILE: '//TRIM(FNAME)
      ELSE
        PRINT*, 'WRONG FORMAT FOR A STOKES PROFILES FILE.'
        PRINT*, 'IS '//TRIM(FNAME)//' A FILE CONTAINING STOKES PROFILES?'
        STOP
      ENDIF
      !
      ! NUMBER OF RECORDS:
      RNREC=INT(TOREAD(4))
      ! NDIMS:
      NSDIMS=INT(TOREAD(5))
      IF (NSDIMS.NE.4) THEN
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      DEALLOCATE(TOREAD)
      !
      NTOT=5+NSDIMS
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      NIX=INT(TOREAD(6))
      NIY=INT(TOREAD(7))
      NIW=INT(TOREAD(8))
      NIS=INT(TOREAD(9))
      !
      DEALLOCATE(TOREAD)
      !
      ! NEW ORDER: IF ((NIS.NE.4).OR.(NIX.NE.SZ(1)).OR.(NIY.NE.SZ(2)).OR.(NIW.NE.SZ(3))) THEN
      IF ((NIS.NE.4).OR.(NIX.NE.SZ(4)).OR.(NIY.NE.SZ(3)).OR.(NIW.NE.SZ(2))) THEN
        WRITE(*,*) 'EITHER THERE ARE NOT ENOUGH ATMOSPHERE PARAMETERS'
        WRITE(*,*) 'OR THE NUMBER OF X,Y,AND/OR Z GRIDS IN THE FILE '
        WRITE(*,*) TRIM(FNAME)//'AND THE ONE FROM THE INPUTFILE DO NOT'
        WRITE(*,*) ' MATCH OR BOTH'
        WRITE(*,*) 'NX', NIX, SZ(4)
        WRITE(*,*) 'NY', NIY, SZ(3)
        WRITE(*,*) 'NW', NIW, SZ(2)
        WRITE(*,*) 'NS', 4, SZ(1)
        STOP
      ENDIF
      !
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA AND STORE IT:
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIW*NIS
      CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
      !
      TO_BE_READ=LNTOT
      LOFFSET=1
      !
      OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
        DO I=1,RNREC
          IF (I.EQ.1) THEN
            LNTOT=5
            LNTOT=LNTOT+NSDIMS
          ELSE IF (I.EQ.2) THEN
              LNTOT=NIW
          ELSE IF (I.EQ.3) THEN
              LNTOT=NIW
          ELSE
            LNTOT=1
            IF (TO_BE_READ.GT.nmaxrec) THEN 
              LNTOT=nmaxrec
            ELSE
              LNTOT=TO_BE_READ
            ENDIF
          ENDIF
          ALLOCATE(TOREAD(LNTOT))
          READ(1) TOREAD
          ! STORE:
          IF (I.EQ.2) THEN
!PRINT*, 'I AM GIVING ERROR HERE: ', SHAPE(IND), SHAPE(TMP1D), SHAPE(TOREAD)
            IND=DBLE(TOREAD)
          ELSE IF (I.EQ.3) THEN
            WAVE=DBLE(TOREAD)
          ELSE IF (I.GT.3) THEN
!PRINT*, 'PRE:', SIZE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SIZE(TOREAD)
            TMP1D(LOFFSET:LOFFSET+LNTOT-1)=TOREAD
            LOFFSET=LOFFSET+LNTOT
            TO_BE_READ=TO_BE_READ-LNTOT
          ENDIF
          DEALLOCATE(TOREAD)
        ENDDO
      CLOSE(UNIT=1)
      ! SAVE IT IN THE ATMOSPHERE:
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA:
      CALL ALLOCATE_3D_SP(TMP,NIX,NIY,NIW,'READ_ALLOCATION TMP')
      !
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIW
      LOFFSET=1
      !
      DO I=1,NIS
        TMP(:,:,:)=0.E0
        TMP(:,:,:)=RESHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1),(/NIX,NIY,NIW/))
        !PRINT*, LOFFSET,LOFFSET+LNTOT,SHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SHAPE(TMP)&
        !    , SUM(TMP/REAL(LNTOT))
        STOKES(I,:,:,:)=RESHAPE(TMP, (/NIW,NIY,NIX/), ORDER=(/3,2,1/))
        LOFFSET=LOFFSET+LNTOT
      ENDDO
      !
      DEALLOCATE(TMP)
      DEALLOCATE(TMP1D)
      !END VERSION 2
!--------------------------------------------------------------------
    ELSE IF (VV.EQ.3) THEN
      ! READ THE FIRST 5 ELEMENTS:
      !
      NTOT=5
      !
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      ! VERSION: ALREADY CHECKED!
      ! ID
      RFID=INT(TOREAD(3))
      IF ( ABS(RFID-160904).LT.0.01) THEN
        PRINT*, 'READING STOKES PROFILES FILE: '//TRIM(FNAME)
      ELSE
        PRINT*, 'WRONG FORMAT FOR A STOKES PROFILES FILE.'
        PRINT*, 'IS '//TRIM(FNAME)//' A FILE CONTAINING STOKES PROFILES?'
        STOP
      ENDIF
      !
      ! NUMBER OF RECORDS:
      RNREC=INT(TOREAD(4))
      ! NDIMS:
      NSDIMS=INT(TOREAD(5))
      IF (NSDIMS.NE.4) THEN
        WRITE(*,*) 'WRONG FORMAT FOR A MODEL ATMOSPHERE FILE.'
        WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING MODEL ATMOSPHERE?'
        STOP
      ENDIF
      !
      DEALLOCATE(TOREAD)
      !
      NTOT=5+NSDIMS
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      NIX=INT(TOREAD(6))
      NIY=INT(TOREAD(7))
      NIW=INT(TOREAD(8))
      NIS=INT(TOREAD(9))
      !
      DEALLOCATE(TOREAD)
      !
      ! NEW ORDER: IF ((NIS.NE.4).OR.(NIX.NE.SZ(1)).OR.(NIY.NE.SZ(2)).OR.(NIW.NE.SZ(3))) THEN
      IF ((NIS.NE.4).OR.(NIX.NE.SZ(4)).OR.(NIY.NE.SZ(3)).OR.(NIW.NE.SZ(2))) THEN
        WRITE(*,*) 'EITHER THERE ARE NOT ENOUGH ATMOSPHERE PARAMETERS'
        WRITE(*,*) 'OR THE NUMBER OF X,Y,AND/OR Z GRIDS IN THE FILE '
        WRITE(*,*) TRIM(FNAME)//'AND THE ONE FROM THE INPUTFILE DO NOT'
        WRITE(*,*) ' MATCH OR BOTH'
        WRITE(*,*) 'NX', NIX, SZ(4)
        WRITE(*,*) 'NY', NIY, SZ(3)
        WRITE(*,*) 'NW', NIW, SZ(2)
        WRITE(*,*) 'NS', 4, SZ(1)
        STOP
      ENDIF
      !
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA AND STORE IT:
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY*NIW*NIS
      CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
      !
      TO_BE_READ=LNTOT
      LOFFSET=1
!PRINT*, " *Reading* ", NIX, NIY, NIW, NIS, RNREC
      !
      OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
        DO I=1,RNREC
          IF (I.EQ.1) THEN
            LNTOT=5
            LNTOT=LNTOT+NSDIMS
          ELSE IF (I.EQ.2) THEN
              LNTOT=NIW
          ELSE IF (I.EQ.3) THEN
              LNTOT=NIW
          ELSE
            LNTOT=1
            IF (TO_BE_READ.GT.nmaxrec) THEN 
              LNTOT=nmaxrec
            ELSE
              LNTOT=TO_BE_READ
            ENDIF
          ENDIF
!PRINT*, "     -> ", I, LNTOT
          ALLOCATE(TOREAD(LNTOT))
!PRINT*, "     -> 0)", SIZE(TOREAD)
          READ(1,IOSTAT=IERR) TOREAD
          IF (IERR.NE.0) THEN
            PRINT*, '   ***   '
            PRINT*, ' Error reading: '//TRIM(FNAME)
            PRINT*, ' Error code: ', IERR
            PRINT*, SUM(ABS(TOREAD))
            PRINT*, '   ___   '
          ENDIF 

!PRINT*, "     -> 0.1)", SIZE(TOREAD)
          ! STORE:
          IF (I.EQ.2) THEN
!PRINT*, 'I AM GIVING ERROR HERE: ', SHAPE(IND), SHAPE(TMP1D), SHAPE(TOREAD)
            IND=DBLE(TOREAD)
          ELSE IF (I.EQ.3) THEN
            WAVE=DBLE(TOREAD)
          ELSE IF (I.GT.3) THEN
!PRINT*, 'PRE:', SIZE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SIZE(TOREAD)
!PRINT*, "     -> a)"
            TMP1D(LOFFSET:LOFFSET+LNTOT-1)=TOREAD
!PRINT*, "     -> b)"
            LOFFSET=LOFFSET+LNTOT
!PRINT*, "     -> c)"
            TO_BE_READ=TO_BE_READ-LNTOT
          ENDIF
          DEALLOCATE(TOREAD)
        ENDDO
      CLOSE(UNIT=1)
      ! SAVE IT IN THE ATMOSPHERE:
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA:
      !CALL ALLOCATE_3D_SP(TMP,NIX,NIY,NIW,'READ_ALLOCATION TMP')
      !
      LNTOT=1
      !LNTOT=LNTOT*NIX*NIY*NIW
      LNTOT=LNTOT*NIS*NIW*NIY
      LOFFSET=1
      !
      DO I=1,NIX
        !#TMP(:,:,:)=0.E0
        !#TMP(:,:,:)=RESHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1),(/NIW,NIY,NIW/))
        !#PRINT*, LOFFSET,LOFFSET+LNTOT,SHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1)), SHAPE(TMP)&
        !#    , SUM(TMP/REAL(LNTOT))
        !#STOKES(I,:,:,:)=RESHAPE(TMP, (/NIW,NIY,NIX/), ORDER=(/3,2,1/))
        STOKES(:,:,:,I)=RESHAPE(TMP1D(LOFFSET:LOFFSET+LNTOT-1) &
            , (/NIS,NIW,NIY/))
        LOFFSET=LOFFSET+LNTOT
      ENDDO
      !PRINT*, ': ', SUM(SUM(SUM(STOKES,DIM=2),DIM=2),DIM=2)
      !
      !DEALLOCATE(TMP)
      DEALLOCATE(TMP1D)
      !END VERSION 3
    ENDIF
    !
  END SUBROUTINE READ_PROFILE_3D
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_RFS3D_V1(FNAME,SZ,Z,IND,WAVE,DSTOKES,SS,VV)

    CHARACTER(*), INTENT(IN)                                               :: FNAME
    INTEGER,DIMENSION(6),INTENT(IN)                                        :: SZ
    !
    REAL(DP),DIMENSION(SZ(3)),INTENT(IN)                                   :: Z
    REAL(DP), INTENT(IN), DIMENSION(SZ(4))                                  :: IND
    !INTEGER, INTENT(IN), DIMENSION(SZ(4))                                  :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(4))                                 :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4),SZ(5),SZ(6))   :: DSTOKES
    !
    INTEGER, INTENT(IN)                                                    :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                                    :: TOWRITE
    !
    INTEGER                                                                :: NTOT, I, NLOOP, NDER
!
    INTEGER :: IERR
    !
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=11
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_RFS3D_V1 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=18060904.E0
    TOWRITE(4)=REAL(1+1+1+1+SZ(1))         ! NRECS
    ! SUCH HUGE NUMBER OF RECORDS IS NEEDED BECAUSE EACH RECORD IS LIMITED TO 2GB OF SIZE
    TOWRITE(5)=6.E0          ! NDIMS
    TOWRITE(6)=REAL(SZ(1))   ! NX
    TOWRITE(7)=REAL(SZ(2))   ! NY
    TOWRITE(8)=REAL(SZ(3))   ! NZ
    TOWRITE(9)=REAL(SZ(4))   ! NW
    TOWRITE(10)=REAL(SZ(5))   ! N STOKES (4)
    TOWRITE(11)=REAL(8) ! ONLY: dS/dT_p, dS/dP_t, dS/dR_t, dS/dBx, dS/dBy, dS/dBz, dS/dVz, dS/dP0
    !
    NLOOP=INT(TOWRITE(4))-1
    NDER=INT(TOWRITE(11))
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      DO I=1,NLOOP
        SELECT CASE (I)
          CASE(1)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(Z)
          CASE(2)
            NTOT=INT(SZ(4))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(IND)
          CASE(3)
            NTOT=INT(SZ(4))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(WAVE)
          CASE DEFAULT
            NTOT=INT(SZ(2)*SZ(3)*SZ(4)*SZ(5)*NDER)
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=RESHAPE(DSTOKES(I-3,:,:,:,:,1:NDER),(/NTOT/))
        ENDSELECT
        WRITE(1) TOWRITE
        DEALLOCATE(TOWRITE)
      ENDDO
    CLOSE(UNIT=1)
    !
  END SUBROUTINE WRITE_RFS3D_V1
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_RFS3D_V2(FNAME,SZ,Z,IND,WAVE,DSTOKES,SS,VV)

    CHARACTER(*), INTENT(IN)                                               :: FNAME
    INTEGER,DIMENSION(6),INTENT(IN)                                        :: SZ
    !
    REAL(DP),DIMENSION(SZ(4)),INTENT(IN)                                   :: Z
    !INTEGER, INTENT(IN), DIMENSION(SZ(3))                                  :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(3))                                  :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(3))                                 :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4),SZ(5),SZ(6))   :: DSTOKES
    !
    INTEGER, INTENT(IN)                                                    :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                                    :: TOWRITE
    !
    INTEGER                                                                :: NTOT, I, J, NLOOP, NDER, NRECS
!
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    REAL(SP), DIMENSION(:,:,:,:,:), ALLOCATABLE          :: TMP5D
    REAL(SP), DIMENSION(:,:,:,:), ALLOCATABLE          :: TMP4D
    !
    INTEGER :: IERR
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=11
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_RFS3D_V2 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=18060904.E0
    ! NUMBER OF RECORDS TO BE STORED:
    LNTOT=1
    LNTOT=LNTOT*SZ(1)*SZ(2)*SZ(3)*SZ(4)*SZ(5)*SZ(6)
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    TOWRITE(4)=REAL(1+1+1+1+NRECS)         ! NRECS
    ! SUCH HUGE NUMBER OF RECORDS IS NEEDED BECAUSE EACH RECORD IS LIMITED TO 2GB OF SIZE
    TOWRITE(5)=6.E0          ! NDIMS
    TOWRITE(6)=REAL(SZ(6))   ! NX
    TOWRITE(7)=REAL(SZ(5))   ! NY
    TOWRITE(8)=REAL(SZ(4))   ! NZ
    TOWRITE(9)=REAL(SZ(3))   ! NW
    TOWRITE(10)=REAL(SZ(1))  ! N STOKES (4)
    TOWRITE(11)=REAL(SZ(2))  ! ONLY: dS/dT_r, dS/dP_t, dS/dR_t, dS/dBx, dS/dBy, dS/dBz, dS/dVz, dS/dP0
    !
    NLOOP=3!TOWRITE(4)-1
    NDER=INT(TOWRITE(11))
    ! STORE PROFILES IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
    !
    LOFFSET=1
    IF (REAL(LNTOT/SZ(2)).LT.REAL(nmaxrec)) THEN
      !PRINT*, '?'
      NTOT=SZ(6)*SZ(5)*SZ(4)*SZ(3)*SZ(1)   ! nx*ny*nz*nw*ns
      !PRINT*, '?', NTOT, SZ
      CALL ALLOCATE_5D_SP(TMP5D,SZ(6),SZ(5),SZ(4),SZ(3),SZ(1),'TMP5D')
      DO I=1,SZ(2)   ! np
        TMP5D=RESHAPE(DSTOKES(:,I,:,:,:,:), (/SZ(6),SZ(5),SZ(4),SZ(3),SZ(1)/)&
            , ORDER=(/5,4,3,2,1/))
        TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP5D(:,:,:,:,:),(/NTOT/))
        LOFFSET=LOFFSET+NTOT
      ENDDO
      DEALLOCATE(TMP5D)
    ELSE IF (REAL(LNTOT/SZ(2)/SZ(1)).LT.REAL(nmaxrec)) THEN
      !PRINT*, '!'
      NTOT=SZ(6)*SZ(5)*SZ(4)*SZ(3) ! nx*ny*nz*nw
      CALL ALLOCATE_4D_SP(TMP4D,SZ(6),SZ(5),SZ(4),SZ(3),'TMP4D')
      DO I=1,SZ(2)   ! np
        DO J=1,SZ(1) ! ns
          TMP4D=RESHAPE(DSTOKES(J,I,:,:,:,:), (/SZ(6),SZ(5),SZ(4),SZ(3)/)&
              , ORDER=(/4,3,2,1/))
          TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(TMP4D,(/NTOT/))
          !TMP1D(LOFFSET:LOFFSET+NTOT-1)=RESHAPE(DSTOKES(:,:,:,:,J,I),(/NTOT/))
          LOFFSET=LOFFSET+NTOT
        ENDDO
      ENDDO
      DEALLOCATE(TMP4D)
    ELSE
      PRINT*, 'NOT IMPLEMENTED FEATURE IN: write_rfs3d_v2'
      STOP
    ENDIF
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      DO I=1,NLOOP
        SELECT CASE (I)
          CASE(1)
            NTOT=INT(SZ(4))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(Z)
          CASE(2)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(IND)
          CASE(3)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(WAVE)
        ENDSELECT
        WRITE(1) TOWRITE
        DEALLOCATE(TOWRITE)
      ENDDO
      ! NOW, STORE THE RESPONSE FUNCTION ITSELF:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    DEALLOCATE(TMP1D)
!PRINT*, 'I AM DONE WITH WRITE_RFS3D_V2'
    !
  END SUBROUTINE WRITE_RFS3D_V2
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_RFS3D_V3(FNAME,SZ,Z,LOGV,IND,WAVE,DSTOKES,SS,VV)

    CHARACTER(*), INTENT(IN)                                               :: FNAME
    INTEGER,DIMENSION(6),INTENT(IN)                                        :: SZ
    !
    REAL(DP),DIMENSION(SZ(4)),INTENT(IN)                                   :: Z
    LOGICAL,DIMENSION(SZ(1)),INTENT(IN)                                    :: LOGV
    REAL(DP), INTENT(IN), DIMENSION(SZ(3))                                 :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(3))                                 :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4),SZ(5),SZ(6))   :: DSTOKES
    !
    INTEGER, INTENT(IN)                                                    :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                                    :: TOWRITE
    !
    INTEGER                                  :: NTOT, I, J, NLOOP, NDER, NRECS
    INTEGER                                  :: K,W,P,S,NPARTOSAV
    REAL(SP),DIMENSION(SZ(1))                 :: LOGVSP
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    REAL(SP), DIMENSION(:,:,:,:,:), ALLOCATABLE          :: TMP5D
    REAL(SP), DIMENSION(:,:,:,:), ALLOCATABLE          :: TMP4D
    !
    INTEGER :: IERR
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=11
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_RFS3D_V2 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=18060904.E0
    ! NUMBER OF RECORDS TO BE STORED:
    LNTOT=1
    !
    ! Take into account that not all the parameters are necessarily calculated:
    NPARTOSAV=0
    LOGVSP(:)=0.0E0
    DO P=1,SZ(1)
      IF (LOGV(P).EQV..TRUE.) THEN
        NPARTOSAV=NPARTOSAV+1
        LOGVSP(P)=1.0E0
      ENDIF
    ENDDO
    !
    LNTOT=LNTOT*NPARTOSAV*SZ(2)*SZ(3)*SZ(4)*SZ(5)*SZ(6)
    !
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    TOWRITE(4)=REAL(1+1+1+1+1+NRECS)         ! NRECS: HEADER, Z, IND, WAVE, LOGV + DATA
    ! SUCH HUGE NUMBER OF RECORDS IS NEEDED BECAUSE EACH RECORD IS LIMITED TO 2GB OF SIZE
    TOWRITE(5)=6.E0          ! NDIMS
    TOWRITE(6)=REAL(SZ(6))   ! NX
    TOWRITE(7)=REAL(SZ(5))   ! NY
    TOWRITE(8)=REAL(SZ(4))   ! NZ
    TOWRITE(9)=REAL(SZ(3))   ! NW
    TOWRITE(10)=REAL(SZ(2))  ! N STOKES (4)
    TOWRITE(11)=REAL(SZ(1))  ! ONLY: dS/dT_r, dS/dP_t, dS/dR_t, dS/dBx, dS/dBy, dS/dBz, dS/dVz, dS/dP0
    !
    NLOOP=4!TOWRITE(4)-1
    NDER=INT(TOWRITE(11))
    !
    ! STORE PROFILES IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
    !
    LOFFSET=1
    !
    DO I=1,SZ(6)
      DO J=1,SZ(5)
        DO K=1,SZ(4)
          DO W=1,SZ(3)
            DO S=1,SZ(2)
              DO P=1,SZ(1)
                IF (LOGV(P).EQV..TRUE.) THEN
                  TMP1D(LOFFSET)=DSTOKES(P,S,W,K,J,I)
                  LOFFSET=LOFFSET+1
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      DO I=1,NLOOP
        SELECT CASE (I)
          CASE(1)
            NTOT=INT(SZ(4))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(Z)
          CASE(2)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(IND)
          CASE(3)
            NTOT=INT(SZ(3))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(WAVE)
          CASE(4)
            NTOT=INT(SZ(1))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(LOGVSP)
        ENDSELECT
        WRITE(1) TOWRITE
        DEALLOCATE(TOWRITE)
      ENDDO
      ! NOW, STORE THE RESPONSE FUNCTION ITSELF:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    DEALLOCATE(TMP1D)
!PRINT*, 'I AM DONE WITH WRITE_RFS3D_V2'
    !
  END SUBROUTINE WRITE_RFS3D_V3
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_RFS2(FNAME,SZ,Z,LOGV,IND,WAVE,DSTOKES,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                                            :: FNAME
    INTEGER,DIMENSION(6),INTENT(IN)                                     :: SZ
    !
    REAL(DP),DIMENSION(SZ(4)),INTENT(IN)                                :: Z
    REAL(DP),DIMENSION(SZ(3)),INTENT(IN)                                :: IND
    LOGICAL,DIMENSION(SZ(1)),INTENT(IN)                                 :: LOGV
    REAL(DP),DIMENSION(SZ(3)),INTENT(IN)                                :: WAVE
    REAL(SP),DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4),SZ(5),SZ(6)),INTENT(IN)  :: DSTOKES
    !
    INTEGER, INTENT(IN)                                                 :: SS, VV
    !
    IF (SS.EQ.3000) THEN
      IF ((VV.EQ.1).OR.(VV.EQ.2)) THEN
        PRINT*, 'I MUST NOT BE HERE! STOP'
        STOP
      ELSE IF (VV.EQ.3) THEN
        CALL WRITE_RFS3D_V3(FNAME,SZ,Z,LOGV,IND,WAVE,DSTOKES,SS,VV)
      ENDIF
    ENDIF
    !
    WRITE(*,*) 'WRITING RESPONSE FUNCTIONS IN FILE: ' // TRIM(FNAME)
    !
  END SUBROUTINE WRITE_RFS2
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_RFS(FNAME,SZ,Z,IND,WAVE,DSTOKES,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                                            :: FNAME
    INTEGER,DIMENSION(6),INTENT(IN)                                     :: SZ
    !
    REAL(DP),DIMENSION(SZ(4)),INTENT(IN)                                :: Z
    REAL(DP),DIMENSION(SZ(3)),INTENT(IN)                                 :: IND
    !INTEGER,DIMENSION(SZ(3)),INTENT(IN)                                 :: IND
    REAL(DP),DIMENSION(SZ(3)),INTENT(IN)                                :: WAVE
    REAL(SP),DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4),SZ(5),SZ(6)),INTENT(IN)  :: DSTOKES
    !
    INTEGER, INTENT(IN)                                                 :: SS, VV
    !
    IF (SS.EQ.3000) THEN
      IF (VV.EQ.1) THEN
        PRINT*, 'I MUST NOT BE HERE! STOP'
        STOP
        CALL WRITE_RFS3D_V1(FNAME,SZ,Z,IND,WAVE,DSTOKES,SS,VV)
      ELSE IF ((VV.EQ.2).OR.(VV.EQ.3)) THEN
        CALL WRITE_RFS3D_V2(FNAME,SZ,Z,IND,WAVE,DSTOKES,SS,VV)
      ENDIF
    ENDIF
    !
    WRITE(*,*) 'WRITING RESPONSE FUNCTIONS IN FILE: ' // TRIM(FNAME)
    !
    !WRITE(*,*) 'I AM write_rfs, BUT I AM NOT IMPLEMENTED YET O_O'
    !
  END SUBROUTINE WRITE_RFS
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_LOGTAUS(FNAME,SZ,Z,IND,WAVE,TAUS,SS,VV)
    !
    CHARACTER(*), INTENT(IN)                                :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                         :: SZ
    !
    REAL(DP),DIMENSION(SZ(2)),INTENT(IN)                    :: Z
    REAL(DP),DIMENSION(SZ(1)),INTENT(IN)                    :: IND
    REAL(DP),DIMENSION(SZ(1)),INTENT(IN)                    :: WAVE
    REAL(SP),DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4)),INTENT(IN)  :: TAUS
    !
    INTEGER, INTENT(IN)                                     :: SS, VV
    !
    WRITE(*,*) 'WRITING SPECTRAL LINE OPTICAL DEPTHS IN FILE: ' // TRIM(FNAME)
    IF (SS.EQ.3000) THEN
      IF ((VV.EQ.1).OR.(VV.EQ.1)) THEN
        PRINT*, 'I MUST NOT BE HERE! STOP'
        STOP
      ELSE IF (VV.EQ.3) THEN
        CALL WRITE_LOGTAUS_V3(FNAME,SZ,Z,IND,WAVE,TAUS,SS,VV)
      ENDIF
    ENDIF
    !
  END SUBROUTINE WRITE_LOGTAUS
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_LOGTAUS_V3(FNAME,SZ,Z,IND,WAVE,TAUS,SS,VV)

    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER,DIMENSION(4),INTENT(IN)                            :: SZ
    !
    REAL(DP),DIMENSION(SZ(2)),INTENT(IN)                       :: Z
    REAL(DP), INTENT(IN), DIMENSION(SZ(1))                     :: IND
    REAL(DP), INTENT(IN), DIMENSION(SZ(1))                     :: WAVE
    REAL(SP), INTENT(IN), DIMENSION(SZ(1),SZ(2),SZ(3),SZ(4))   :: TAUS
    !
    INTEGER, INTENT(IN)                                        :: SS, VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE                        :: TOWRITE
    !
    INTEGER                                  :: NTOT, I, J, NLOOP, NDER, NRECS
    INTEGER                                  :: K,W,P,S,NPARTOSAV
    REAL(SP),DIMENSION(SZ(1))                 :: LOGVSP
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_WRITTEN
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    REAL(SP), DIMENSION(:,:,:,:,:), ALLOCATABLE          :: TMP5D
    REAL(SP), DIMENSION(:,:,:,:), ALLOCATABLE          :: TMP4D
    !
    INTEGER :: IERR
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=9
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_RFS3D_V2 ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=12200904.E0
    ! NUMBER OF RECORDS TO BE STORED:
    LNTOT=1
    LNTOT=LNTOT*SZ(1)*SZ(2)*SZ(3)*SZ(4)
    !
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    TOWRITE(4)=REAL(1+1+1+1+NRECS)         ! NRECS: HEADER, Z, IND, WAVE + DATA
    ! SUCH HUGE NUMBER OF RECORDS IS NEEDED BECAUSE EACH RECORD IS LIMITED TO 2GB OF SIZE
    TOWRITE(5)=4.E0          ! NDIMS
    TOWRITE(6)=REAL(SZ(4))   ! NX
    TOWRITE(7)=REAL(SZ(3))   ! NY
    TOWRITE(8)=REAL(SZ(2))   ! NZ
    TOWRITE(9)=REAL(SZ(1))   ! NW
    !
    NLOOP=3
    !
    ! STORE PROFILES IN A 1D LARRAY
    CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
    !
    LOFFSET=1
    !
    DO I=1,SZ(4)
      DO J=1,SZ(3)
        DO K=1,SZ(2)
          !DO W=1,SZ(1)
          TMP1D(LOFFSET:LOFFSET+SZ(1)-1)=TAUS(:,K,J,I)
          LOFFSET=LOFFSET+SZ(1)
          !ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! SEQUENTIAL ACCESS TO WRITE EACH PARAMETER:
      DO I=1,NLOOP
        SELECT CASE (I)
          CASE(1)
            NTOT=INT(SZ(2))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(Z)
          CASE(2)
            NTOT=INT(SZ(1))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(IND)
          CASE(3)
            NTOT=INT(SZ(1))
            ALLOCATE(TOWRITE(NTOT))
            TOWRITE(:)=REAL(WAVE)
        ENDSELECT
        WRITE(1) TOWRITE
        DEALLOCATE(TOWRITE)
      ENDDO
      ! NOW, STORE THE RESPONSE FUNCTION ITSELF:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) TMP1D(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
    DEALLOCATE(TMP1D)
    !
  END SUBROUTINE WRITE_LOGTAUS_V3
  !
  !------------------------------------------------
  !
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_PSF_2D(FNAME,PSF)
    !
    CHARACTER(*), INTENT(IN)                                   :: FNAME
    !
    REAL(SP), INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE            :: PSF
    !
    INTEGER                                  :: VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOREAD
    !
    INTEGER                                  :: NTOT
    INTEGER                                  :: POSV
    INTEGER                                  :: NEGV
    INTEGER                                  :: MEDV
    INTEGER                                  :: RFID
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_READ
    !
    INTEGER                                  :: NSDIMS
    INTEGER                                  :: NIX, NIY, NIW, NIS, RNREC
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER                                  :: I
    INTEGER :: IERR
    !
    !
    NTOT=5
    !
    CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
    !
    CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
    !
    POSV=INT(TOREAD(1))
    NEGV=INT(TOREAD(2))
    !
    MEDV=(POSV+NEGV)/2
    !
    POSV=POSV-MEDV
    NEGV=NEGV-MEDV
    !
    VV=(POSV - NEGV) / 2
    !
    IF ( (POSV+NEGV).GT.0.01) THEN
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'NON RECOGNIZED FORMAT. STOPPING'
      STOP
    ELSE
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'VERSION ', (POSV - NEGV) / 2
    ENDIF
    !
    ! ID
    RFID=NINT(TOREAD(3))
    IF ( ABS(RFID-161906).LT.0.01) THEN
      !PRINT*, 'READING PSF FILE: '//TRIM(FNAME)
    ELSE
      PRINT*, 'WRONG FORMAT FOR AN SPATIAL PSF FILE.'
      PRINT*, 'IS '//TRIM(FNAME)//' A FILE CONTAINING AN SPATIAL PSF?'
      STOP
    ENDIF
    !
    ! NUMBER OF RECORDS:
    RNREC=INT(TOREAD(4))
    ! NDIMS:
    NSDIMS=INT(TOREAD(5))
    IF (NSDIMS.NE.2) THEN
      WRITE(*,*) 'WRONG FORMAT FOR AN SPATIAL PSF FILE.'
      WRITE(*,*) 'IS '//TRIM(FNAME)//' A FILE CONTAINING AN SPATIAL PSF FILE?'
      STOP
    ENDIF
    !
    DEALLOCATE(TOREAD)
    IF (VV.EQ.1) THEN
      WRITE(*,*) 'NOT IMPLEMENTED YET!'
      STOP
      ! END VERSION 1
!--------------------------------------------------------------------
    ELSE IF (VV.EQ.2) THEN
      ! READ THE FIRST 5 ELEMENTS PLUS DIMENSIONS:
      !
      NTOT=5+NSDIMS
      CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
      !
      CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
      !
      NIX=INT(TOREAD(6))
      NIY=INT(TOREAD(7))
!PRINT*, NIX, NIY
      !
      DEALLOCATE(TOREAD)
      !
      ! NOW WE HAVE TO ACCESS ITERATIVELY TO THE DATA AND STORE IT:
      LNTOT=1
      LNTOT=LNTOT*NIX*NIY
      CALL ALLOCATE_L1D_SP(TMP1D,LNTOT,'READ_ALLOCATION TMP')
      !
      TO_BE_READ=LNTOT
      LOFFSET=1
      !
      OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
        DO I=1,RNREC
          IF (I.EQ.1) THEN
            LNTOT=5
            LNTOT=LNTOT+NSDIMS
          ELSE
            LNTOT=1
            IF (TO_BE_READ.GT.nmaxrec) THEN
              LNTOT=nmaxrec
            ELSE
              LNTOT=TO_BE_READ
            ENDIF
          ENDIF
          ALLOCATE(TOREAD(LNTOT))
          READ(1) TOREAD
          ! STORE:
          IF (I.GT.1) THEN
            TMP1D(LOFFSET:LOFFSET+LNTOT-1)=TOREAD
            LOFFSET=LOFFSET+LNTOT
            TO_BE_READ=TO_BE_READ-LNTOT
          ENDIF
          DEALLOCATE(TOREAD)
        ENDDO
      CLOSE(UNIT=1)
      !
      ! SAVE IT IN THE PSF:
      !
      IF (ALLOCATED(PSF)) DEALLOCATE(PSF)
      ALLOCATE(PSF(NIY,NIX))
      PSF(:,:)=RESHAPE(RESHAPE(TMP1D(:),(/NIX,NIY/)), (/NIY,NIX/), ORDER=(/2,1/))
      !
      DEALLOCATE(TMP1D)
      !END VERSION 2
    ENDIF
    !
  END SUBROUTINE READ_PSF_2D
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_BIN(FNAME,RK,SH,SZ,ARRAY,SS,VV)

    CHARACTER(*), INTENT(IN)                   :: FNAME
    INTEGER,INTENT(IN)                         :: RK
    INTEGER,DIMENSION(RK),INTENT(IN)           :: SH
    INTEGER(kind=8),INTENT(IN)                 :: SZ
    !
    REAL(SP), INTENT(IN), DIMENSION(SZ)        :: ARRAY
    !
    INTEGER, INTENT(IN)                        :: SS, VV
    !
!
    !
    INTEGER                                    :: NTOT, I, NRECS
    INTEGER(kind=8)                            :: LNTOT, LOFFSET, TO_BE_WRITTEN
    REAL(SP), ALLOCATABLE, DIMENSION(:)        :: TOWRITE
    INTEGER :: IERR
    !
    ! FIRST RECORD:
    !             POSV, NEGV, ID, NRECS, NDIMS, SIZEDIM1, SIZEDIM2, SIZEDIM3, SIZEDIM4
    NTOT=5+RK
    CALL ALLOCATE_1D_SP(TOWRITE,NTOT,'WRITE_BIN ALLOCATION')
    !
    TOWRITE(1)=REAL(SS+VV)
    TOWRITE(2)=REAL(SS-VV)
    TOWRITE(3)=160206.E0
    ! NUMBER OF RECORDS TO BE STORED:
    LNTOT=SZ
    NRECS=INT(CEILING(REAL(LNTOT)/REAL(nmaxrec)))
    TOWRITE(4)=REAL(1+NRECS)         ! NRECS
    ! SUCH HUGE NUMBER OF RECORDS IS NEEDED BECAUSE EACH RECORD IS LIMITED TO 2GB OF SIZE
    TOWRITE(5)=REAL(RK)          ! NDIMS
    DO I=1,RK
      TOWRITE(5+I)=REAL(SH(I))
    ENDDO
    !
    LOFFSET=1
    !
    OPEN(UNIT=1,FILE=FNAME,FORM="unformatted",IOSTAT=IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '   ***   '
        PRINT*, ' Error reading: '//TRIM(FNAME)
        PRINT*, ' Error code: ', IERR
        PRINT*, '   ___   '
      ENDIF 
      !FIRST ACCESS:
      WRITE(1) TOWRITE
      DEALLOCATE(TOWRITE)
      !
      ! SEQUENTIAL ACCESS TO WRITE EACH RECORD:
      TO_BE_WRITTEN=LNTOT
      LOFFSET=1
      DO I=1,NRECS
        LNTOT=1
        IF (TO_BE_WRITTEN.GT.nmaxrec) THEN
          LNTOT=nmaxrec
        ELSE
          LNTOT=TO_BE_WRITTEN
        ENDIF
        WRITE(1) ARRAY(LOFFSET:LOFFSET+LNTOT-1)
        LOFFSET=LOFFSET+LNTOT
        TO_BE_WRITTEN=TO_BE_WRITTEN-LNTOT
      ENDDO
    CLOSE(UNIT=1)
!PRINT*, 'I AM DONE WITH WRITE_BIN'
    !
  END SUBROUTINE WRITE_BIN
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_LSF_FILE(FNAME, NN, ARRAY, II)
    !
    CHARACTER(*), INTENT(IN)                                   :: FNAME
    INTEGER, INTENT(IN) :: NN, II
    REAL(DP), INTENT(INOUT), DIMENSION(NN) :: ARRAY
    !
    INTEGER                                  :: VV
    !
    REAL(SP), DIMENSION(:), ALLOCATABLE      :: TOREAD
    !
    INTEGER                                  :: NTOT
    INTEGER                                  :: POSV
    INTEGER                                  :: NEGV
    INTEGER                                  :: MEDV
    INTEGER                                  :: RFID
    !
    INTEGER(kind=8)                          :: LNTOT, LOFFSET, TO_BE_READ
    !
    INTEGER                                  :: NSDIMS
    INTEGER                                  :: NIW, RNREC
    !
    REAL(SP), DIMENSION(:,:,:), ALLOCATABLE      :: TMP
    REAL(SP), DIMENSION(:), ALLOCATABLE          :: TMP1D
    !
    INTEGER                                  :: I
    INTEGER :: IERR
    !
    !PRINT*, ''
    !PRINT*, ' + READ_LSF_FILE. Reading: '//TRIM(FNAME)//' +'
    !
    NTOT=5
    !
    CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 1')
    !
    CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
    !
    POSV=INT(TOREAD(1))
    NEGV=INT(TOREAD(2))
    !
    MEDV=(POSV+NEGV)/2
    !
    POSV=POSV-MEDV
    NEGV=NEGV-MEDV
    !
    VV=(POSV - NEGV) / 2
    !
    IF ( (POSV+NEGV).GT.0.01) THEN
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'NON RECOGNIZED FORMAT. STOPPING'
      STOP
    ELSE
      WRITE(*,*) 'READING FILE: '//TRIM(FNAME)
      WRITE(*,*) 'VERSION ', (POSV - NEGV) / 2
    ENDIF
    !
    ! ID
    RFID=NINT(TOREAD(3))
    IF ( ABS(RFID-3330904).LT.0.01) THEN
      !PRINT*, 'READING PSF FILE: '//TRIM(FNAME)
    ELSE
      PRINT*, 'WRONG FORMAT FOR A LINE SPREAD FUNCTION FILE.'
      PRINT*, 'IS '//TRIM(FNAME)//' A FILE CONTAINING A LINE SPREAD FUNCTION?'
      STOP
    ENDIF
    !
    ! NUMBER OF RECORDS:
    RNREC=INT(TOREAD(4))
    ! NDIMS:
    NSDIMS=INT(TOREAD(5))
    IF (NSDIMS.NE.1) THEN
      PRINT*, 'WRONG FORMAT FOR A LINE SPREAD FUNCTION FILE.'
      PRINT*, 'IS '//TRIM(FNAME)//' A FILE CONTAINING A LINE SPREAD FUNCTION?'
      STOP
    ENDIF
    !
    DEALLOCATE(TOREAD)
    !
    ! First, read the whole header:
    NTOT=5+NSDIMS
    CALL ALLOCATE_1D_SP(TOREAD,NTOT,'READ_ALLOCATION 2')
    !
    CALL READ_SP1D(FNAME, SHAPE(TOREAD), TOREAD)
    !
    NIW=INT(TOREAD(6))
    !
    DEALLOCATE(TOREAD)
    !
    ! Check wavelength size supplied as compared to the one from input file:
    IF (NIW.NE.NN) THEN
      PRINT*, ''
      PRINT*, ' Error! Wavelength size for spectral region #=', II
      PRINT*, ' Number of wavelengths for this spectral region from input file: ', NN
      PRINT*, ' Number of wavelengths found in '//TRIM(FNAME)//': ', NN
      PRINT*, ''
      PRINT*, ''
      STOP
    ENDIF

    !
    TO_BE_READ=NIW
    LOFFSET=1
    !
    OPEN(UNIT=1,FILE=TRIM(FNAME),FORM="unformatted",IOSTAT=IERR)
    IF (IERR.NE.0) THEN
      PRINT*, '   ***   '
      PRINT*, ' Error reading: '//TRIM(FNAME)
      PRINT*, ' Error code: ', IERR
      PRINT*, '   ___   '
    ENDIF 
      DO I=1,RNREC
        IF (I.EQ.1) THEN
          LNTOT=5
          LNTOT=LNTOT+NSDIMS
        ELSE
          LNTOT=1
          IF (TO_BE_READ.GT.nmaxrec) THEN
            LNTOT=nmaxrec
          ELSE
            LNTOT=TO_BE_READ
          ENDIF
        ENDIF
        ALLOCATE(TOREAD(LNTOT))
        READ(1) TOREAD
        ! STORE:
        IF (I.GT.1) THEN
          ARRAY(LOFFSET:LOFFSET+LNTOT-1)=DBLE(TOREAD)
          LOFFSET=LOFFSET+LNTOT
          TO_BE_READ=TO_BE_READ-LNTOT
        ENDIF
        DEALLOCATE(TOREAD)
      ENDDO
    CLOSE(UNIT=1)
    !
  END SUBROUTINE READ_LSF_FILE
  !
  !================================================
  !
END MODULE MISC
!
