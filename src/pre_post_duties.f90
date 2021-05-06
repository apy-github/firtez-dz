!
MODULE PRE_POST_DUTIES
  !
  !================================================
  !
  USE CONS_PARAM, ONLY: SP, DP, IMAT
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  PUBLIC :: START_STEP
  PUBLIC :: END_STEP
  !
  PRIVATE :: DEALLOCATE_AND_WRITE
  PRIVATE :: WRITE_STOKES_PROFILES
  PRIVATE :: WRITE_MODEL_ATMOSPHERE
  PRIVATE :: WRITE_RESPONSE_FUNCTIONS
  PRIVATE :: ALLOCATE_AND_READ
  PRIVATE :: READ_PSF
  PRIVATE :: SPATIAL_CONV_INIT
  PRIVATE :: WAVE_INIT
  PRIVATE :: SET_INVERSION_CONSTANTS
  PRIVATE :: LSF_INIT
  PRIVATE :: ALLOCATE_PHYS_ARRAYS
  PRIVATE :: ALLOCATE_STOKES
  PRIVATE :: ALLOCATE_INVERSION
  PRIVATE :: CLOSE_READ_INPUT_FILE
  PRIVATE :: CLOSE_WAVE_INIT
  PRIVATE :: LSF_END
  PRIVATE :: SPATIAL_CONV_END
  PRIVATE :: DEALLOCATE_PHYS_ARRAYS
  PRIVATE :: DEALLOCATE_STOKES
  PRIVATE :: DEALLOCATE_INVERSION
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! start_step
  ! allocate_and_read
  ! read_psf
  ! spatial_conv_init
  ! wave_init
  ! set_inversion_constants
  ! lsf_init
  ! allocate_phys_arrays
  ! allocate_stokes
  ! allocate_inversion
  ! close_read_input_file
  ! close_wave_init
  ! deallocate_phys_arrays
  ! deallocate_stokes
  ! deallocate_inversion
  !
  !________________________________________________
  !
  SUBROUTINE START_STEP()
    !
    USE ATOM_DATABASE, ONLY: ATOM_INIT
    USE LOG, ONLY: LOG_INIT
    USE LINES_DATABASE, ONLY: NREAD_LINES_DATABASE
    USE INPUT_FILE, ONLY: READ_INPUT_FILE, READ_MINPUT_FILE
    USE user_mpi, only: mpi__myrank
    !
    ! FIRST STEP IS DIFFERENT FOR SLAVES AND MASTER...
    ! ...FIRST, THERE IS A COMMON STEP IN WHICH WE READ...
    ! ...THE INPUT FILE AND ALL THE PREVIOUS DUTIES...
    ! ...NEEDED FOR THE READING
    !
    ! Initialize log file
    IF (mpi__myrank.eq.0) CALL LOG_INIT()
    ! Initialize atomic database
    CALL ATOM_INIT()
    !
    ! Prepare to read input file
    CALL READ_MINPUT_FILE()
    !
    CALL NREAD_LINES_DATABASE()
    ! Prepare to read input file
    CALL READ_INPUT_FILE()
    !
    CALL WAVE_INIT()
    !
    ! NOW, WE ALLOCATE THE VARIOUS ARRAYS NEEDED...
    ! ...DEPENDING ON WHETHER I AM THE MASTER...
    ! ...OR A SLAVE:
    CALL ALLOCATE_AND_READ()
    !
    !WRITE(*,*) 'I AM FIRST_STEP', mpi__myrank, NX, NY, NZ
    !
  END SUBROUTINE START_STEP
  !
  !------------------------------------------------
  !
  SUBROUTINE END_STEP()
    !
    CALL DEALLOCATE_AND_WRITE()
    !
  END SUBROUTINE END_STEP
  !
  !------------------------------------------------
  !
  !
  ! Now, private subroutines:
  !
  SUBROUTINE DEALLOCATE_AND_WRITE()
    !
    USE INPUT_FILE, ONLY: CLOSE_READ_INPUT_FILE
    USE user_mpi, ONLY: mpi__myrank
    USE TIME_PARAM
    USE LINES_DATABASE, ONLY: FINI_LINEDATABASE_VARS
    USE FORWARD_PARAM, ONLY: N_FWD_MODELLING
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      CALL SYSTEM_CLOCK(COUNT=TSTART,COUNT_RATE=TRATE,COUNT_MAX=TMAX)
      TIME2B=DBLE(TSTART)
    ENDIF ! Master
    !
    ! Write Stokes profiles if necessary:
    CALL WRITE_STOKES_PROFILES()
    !
    ! Write model atmosphere if necessary (always, actually):
    CALL WRITE_MODEL_ATMOSPHERE()
    !
    ! Write response functions if necessary:
    CALL WRITE_RESPONSE_FUNCTIONS()
    !
    CALL WRITE_LINTAU()
    !
    CALL DEALLOCATE_INVERSION()
    !
    CALL DEALLOCATE_PHYS_ARRAYS()
    !
    CALL DEALLOCATE_STOKES()
    !
    ! We deallocate the various allocated arrays:
    !
    CALL FINI_LINEDATABASE_VARS()
    !
    ! Deallocate input file reading variables:
    CALL CLOSE_READ_INPUT_FILE()
    !
    ! Deallocate wavelength grid
    CALL CLOSE_WAVE_INIT()
    IF (mpi__myrank.EQ.0) THEN
      !
      CALL SYSTEM_CLOCK(COUNT=TSTART,COUNT_RATE=TRATE,COUNT_MAX=TMAX)
      TIME2=DBLE(TSTART)
PRINT*, '*********************************************************************'
PRINT*, ''
PRINT*, 'Number of forward modelling: ', N_FWD_MODELLING
PRINT*, (TIME2-TIME1)/1000., (TIME2B-TIME1B)/1000.
PRINT*, '*********************************************************************'
    ENDIF
    !
  END SUBROUTINE DEALLOCATE_AND_WRITE
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_LINTAU()
    !
    USE MISC, ONLY: WRITE_LOGTAUS
    USE user_mpi, ONLY: mpi__myrank
    USE PHYS_PARAM, ONLY: TAU3DLIN
    USE CODE_MODES, ONLY: MTAULIN, NAMEMODEL, OUTPPATH
    USE GRID_PARAM, ONLY: ZZ
    USE FORWARD_PARAM, ONLY: INDEX, WAVE
    !
    IF (mpi__myrank.EQ.0) THEN
      IF (MTAULIN.EQV..TRUE.) THEN
        CALL WRITE_LOGTAUS(TRIM(OUTPPATH)//'tl_'//TRIM(NAMEMODEL) &
            , SHAPE(TAU3DLIN), ZZ,INDEX,WAVE, TAU3DLIN, 3000, 3)
      ENDIF
    ENDIF
    !
  END SUBROUTINE WRITE_LINTAU
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_STOKES_PROFILES()
    !
    USE MISC, ONLY: WRITE_PROFILES, WRITE_BIN
    USE user_mpi, ONLY: mpi__myrank
    USE USER_FFTW3, ONLY: FFTK2D, FFTY2D, GPLAN2D &
        , FFTW_EXECUTE_DFT_R2C, FFTW_EXECUTE_DFT_C2R
    USE PHYS_PARAM
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_4D_SP
    USE CODE_MODES, ONLY: MSYNTHESIS, MINVERSION, NAMEPROFILE &
        , COUPLED, OUTPPATH
    USE GRID_PARAM, ONLY: NX, NY
    USE INVERT_PARAM, ONLY: INV_STK, NSTKINV
    USE FORWARD_PARAM, ONLY: FULL_STOKES, NUMW, WAVE, INDEX
    !
    INTEGER  :: OFF, I, K
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      IF ( (MSYNTHESIS.EQV..TRUE.) .OR. (MINVERSION.EQV..TRUE.) ) THEN
        !
        PRINT*, 'MASTER WRITING PROFILES: '
        IF (MINVERSION.EQV..TRUE.) THEN
          !
          IF (ALLOCATED(SYN)) DEALLOCATE(SYN)
          CALL ALLOCATE_4D_SP(SYN,4,NUMW,NY,NX,'SYN')
          OFF=0
          DO I=1,4
            IF (INV_STK(I).EQV..TRUE.) THEN
              SYN(I,:,:,:)=SYN3D(OFF+1:OFF+NUMW,:,:)
              OFF=OFF+NUMW
            ENDIF
          ENDDO
          CALL WRITE_PROFILES(TRIM(OUTPPATH)//'out_'//TRIM(NAMEPROFILE), SHAPE(SYN)&
              , INDEX, WAVE, SYN, 3000, 3)
          DEALLOCATE(SYN)

          IF (COUPLED.EQV..TRUE.) THEN

            DO K=1,NSTKINV*NUMW
              !
              ! DFT synthetic data:
              CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),SYN3D(K,:,:),FFTY2D)
              !
              ! Multiply transformed arrays:
              FFTY2D=FFTY2D*FFTK2D
              !
              ! DFT back the product:
              CALL FFTW_EXECUTE_DFT_C2R(GPLAN2D(2),FFTY2D,SYN3D(K,:,:))
              !
            ENDDO
            !
            IF (ALLOCATED(SYN)) DEALLOCATE(SYN)
            CALL ALLOCATE_4D_SP(SYN,4,NUMW,NY,NX,'SYN')
            OFF=0
            DO I=1,4
              IF (INV_STK(I).EQV..TRUE.) THEN
                SYN(I,:,:,:)=SYN3D(OFF+1:OFF+NUMW,:,:)
                OFF=OFF+NUMW
              ENDIF
            ENDDO

            CALL WRITE_PROFILES(TRIM(OUTPPATH)//'cout_'//TRIM(NAMEPROFILE) &
                , SHAPE(SYN), INDEX, WAVE, SYN, 3000, 3)
            DEALLOCATE(SYN)

          ENDIF ! End coupled inversion
          !
        ELSE ! I am not inversion, then:
          !
          ! Transform, from 3D to 4D:
          IF (ALLOCATED(SYN)) DEALLOCATE(SYN)
          CALL ALLOCATE_4D_SP(SYN,4,NUMW,NY,NX,'SYN')
          OFF=0
          DO I=1,4
            IF (INV_STK(I).EQV..TRUE.) THEN
              SYN(I,:,:,:)=SYN3D(OFF+1:OFF+NUMW,:,:)
              OFF=OFF+NUMW
            ENDIF
          ENDDO
          !
          CALL WRITE_PROFILES(TRIM(OUTPPATH)//TRIM(NAMEPROFILE), SHAPE(SYN)&
              , INDEX, WAVE, SYN, 3000, 3)
          !
          ! Check whether we are asked to save FULL_STOKES, and act accordingly:
          IF (FULL_STOKES.EQV..TRUE.) THEN
            CALL WRITE_BIN(TRIM(OUTPPATH)//'fs_'//TRIM(NAMEPROFILE) &
                , SIZE(SHAPE(SYN5D)) &
                , SHAPE(SYN5D), SIZE(SYN5D,KIND=8), SYN5D, 3000, 3)
          ENDIF ! FULL_STOKES
          !
        ENDIF ! End synthesis
        !
      ENDIF ! End either synthesis or inversion
      !
    ENDIF ! Master
    !
  END SUBROUTINE WRITE_STOKES_PROFILES
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_MODEL_ATMOSPHERE()
    !
    USE user_mpi, ONLY: mpi__myrank
    USE MISC, ONLY: WRITE_MODEL
    USE PHYS_PARAM
    USE CODE_MODES, ONLY: MINVERSION, NAMEMODEL, OUTPPATH
    USE GRID_PARAM, ONLY: NX, NY, NZ, XX, YY, ZZ
    !
    INTEGER  :: OFF, I, K
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      IF (MINVERSION.EQV..TRUE.) THEN
        !
        PRINT*, 'WRITING MODEL ERRORS (skip):'
        !CALL WRITE_MODEL(TRIM(OUTPPATH)//'eout_'//TRIM(NAMEMODEL), NX, NY, NZ&
        !        , ETEM3D, EPG3D, ERHO3D, EBX3D, EBY3D, EBZ3D, EVZ3D&
        !        , PEL3D, MW3D, TAU3D5, XX, YY, ZZ, 3000, 3)
        !
      ENDIF ! Inversion
      !
      PRINT*, 'WRITING MODEL:'
      CALL WRITE_MODEL(TRIM(OUTPPATH)//'out_'//TRIM(NAMEMODEL), NX, NY, NZ&
              , TEM3D, PG3D, RHO3D, BX3D, BY3D, BZ3D, VZ3D&
              , PEL3D, MW3D, TAU3D5, XX, YY, ZZ, 3000, 3)
      !
    ENDIF ! Master
    !
  END SUBROUTINE WRITE_MODEL_ATMOSPHERE
  !
  !------------------------------------------------
  !
  SUBROUTINE WRITE_RESPONSE_FUNCTIONS()
    !
    USE user_mpi, ONLY: mpi__myrank
    USE MISC, ONLY: WRITE_RFS, WRITE_RFS2
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_6D_SP
    USE PHYS_PARAM, ONLY: WDSYN, DSYN
    USE CODE_MODES, ONLY: SAVE_RFS, NAMEMODEL, OUTPPATH
    USE GRID_PARAM, ONLY: NX, NY, NZ, ZZ
    USE INVERT_PARAM, ONLY: INV_STK, INV_ATMPAR
    USE FORWARD_PARAM, ONLY: INDEX, WAVE, NUMW
    !
    INTEGER  :: OFF, I, J, POFF, K
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      ! Response functions are only written if explicitly stated:
      IF (SAVE_RFS.EQV..TRUE.) THEN
        PRINT*, 'WRITING response functions:!'
        PRINT*, 'Not just now!'
        CALL ALLOCATE_6D_SP(WDSYN,8,4,NUMW,NZ,NY,NX,'SYN')
        OFF=0
        DO I=1,4
          POFF=0
          IF (INV_STK(I).EQV..TRUE.) THEN
            DO J=1,SIZE(INV_ATMPAR)
!PRINT*, INV_ATMPAR, ' <<<<<<<<<<<<<<<<'
              IF (INV_ATMPAR(J).EQV..TRUE.) THEN
                WDSYN(J,I,:,:,:,:)=DSYN(OFF+1:OFF+NUMW,POFF+1:POFF+NZ,:,:)
                POFF=POFF+NZ
              ENDIF
            ENDDO
            OFF=OFF+NUMW
          ENDIF
        ENDDO
        CALL WRITE_RFS2(TRIM(OUTPPATH)//'rfout_'//TRIM(NAMEMODEL)&
            ,SHAPE(WDSYN),ZZ,INV_ATMPAR,INDEX,WAVE,WDSYN,3000,3)
        !CALL WRITE_RFS(TRIM(OUTPPATH)//'rfout_'//TRIM(NAMEMODEL)&
        !    ,SHAPE(WDSYN),ZZ,INDEX,WAVE,WDSYN,3000,3)
        DEALLOCATE(WDSYN)
      ENDIF ! Save_rfs
      !
    ENDIF ! Master
    !
  END SUBROUTINE WRITE_RESPONSE_FUNCTIONS
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_AND_READ()
    !
    USE MISC, ONLY: READ_MODEL, READ_PROFILES
    USE NORMALIZATION, ONLY: GET_NORM
    !
    USE COUPLED_INVERSION, ONLY: INIT_COUPLED_INVERSION_STATIC_VARS
    USE FORWARD_PARAM, ONLY: NUMW, GAIN1D, N_FWD_MODELLING &
         , WAVE, INDEX
    USE INVERT_PARAM, ONLY: NFREQ, NSTKINV
    USE TIME_PARAM
    USE PHYS_PARAM
    USE user_mpi, ONLY: mpi__myrank
    USE GRID_PARAM, ONLY: NX, NY, NZ, XX, YY, ZZ
    USE INVERT_PARAM, ONLY: INV_STK
    USE CODE_MODES, ONLY: COUPLED, HSRA_NORMALIZATION, MINVERSION, NAMEMODEL &
         , NAMEPROFILE, DATAPATH, MODLPATH!, MSMOOTHING
    !USE EXTEND_2D, ONLY: SMOOTHING
    !
    INTEGER  :: OFF, I,J,K
    !
    ! NUMBER OF FREQUENCIES OBSERVED:
    NFREQ=NUMW * NSTKINV
    !
    !
    ! Clearer alternative:
    CALL ALLOCATE_PHYS_ARRAYS()
    !     and inside, distinguish between master and slave
    ! CALL ALLOCATE_STOKES
    !     idem
    ! CALL ALLOCATE_INVERSION
    !     idem
    !     Actually, this might be inside ALLOCATE_PHYS_ARRAYS?
    !     ?Add here specific TAGNAMES for the various MPI messages?
    !
    ! Normalization and coupled inversion, the latter might be moved to...
    ! ...ALLOCATE_INVERSION?
    ! 
    IF (mpi__myrank.EQ.0) THEN
      !
      ! MASTER
      !
      N_FWD_MODELLING = 0
      !
! Timing 1:
      CALL SYSTEM_CLOCK(COUNT=TSTART,COUNT_RATE=TRATE,COUNT_MAX=TMAX)
      TIME1=DBLE(TSTART)
! End timing 1.
      !
      !CALL ALLOCATE_PHYS3D_ARRAYS(MINVERSION)
      CALL ALLOCATE_STOKES(MINVERSION)
      !
      CALL READ_MODEL(TRIM(MODLPATH)//TRIM(NAMEMODEL), NX, NY, NZ&
              , TEM3D, PG3D, RHO3D, BX3D, BY3D, BZ3D, VZ3D&
              , PEL3D, MW3D, TAU3D5, XX, YY, ZZ)
 !>< ! Avoid T below 2000
 !><     IF (ANY(TEM3D.LT.2000)) THEN
 !><       DO I=1,NX
 !><         DO J=1,NY
 !><           DO K=1,NZ
 !><             IF (TEM3D(K,J,I).LT.2000.0E0) TEM3D(K,J,I)=2000.0E0
 !><           ENDDO
 !><         ENDDO
 !><       ENDDO
 !><     ENDIF
      !
      ! In inversion mode, profiles must be read:
      IF (MINVERSION.EQV..TRUE.) THEN
        !
        ! Load observed profiles:
        CALL READ_PROFILES(TRIM(DATAPATH)//TRIM(NAMEPROFILE), SHAPE(OBS)&
            , INDEX, WAVE, OBS)
        !
        ! Reshape them to a 3D array:
        OFF=0
        DO I=1,4
          IF (INV_STK(I).EQV..TRUE.) THEN
            OBS3D(OFF+1:OFF+NUMW,:,:)=OBS(I,:,:,:)
            OFF=OFF+NUMW
          ENDIF
        ENDDO
        !
        DEALLOCATE(OBS)
        !
      ENDIF
      !
    ELSE
      !
      ! SLAVES:
      !
      !CALL ALLOCATE_PHYS1D_ARRAYS()
      !
      ! Since typically this is gonna be much faster than...
      ! ... master's initial duties, i.e. allocating huge...
      ! ... arrays into memory, each slave is going to...
      ! ... synthesis hsra reference:
      !
      CALL GET_NORM(HSRA_NORMALIZATION,GAIN1D)
      !
    ENDIF
    !
    CALL ALLOCATE_INVERSION()
    !
    ! If we do set spatial coupling:
    !
    IF (COUPLED.EQV..TRUE.) THEN
      CALL READ_PSF()
      CALL INIT_COUPLED_INVERSION_STATIC_VARS()
 !   ELSE
 !     CALL INIT_SMOOTH_2D()
    ENDIF
    !
    !MSMOOTHING=.TRUE.
    !IF (MSMOOTHING.EQV..TRUE.) THEN
    !  CALL SMOOTHING()
    !ENDIF
    !
    IF (mpi__myrank.EQ.0) THEN
! Timing 2:
      CALL SYSTEM_CLOCK(COUNT=TSTART,COUNT_RATE=TRATE,COUNT_MAX=TMAX)
      TIME1B=DBLE(TSTART)
! End timing 2.
    ENDIF
  END SUBROUTINE ALLOCATE_AND_READ
  !
  !------------------------------------------------
  !
  SUBROUTINE READ_PSF()
    !
    USE COUPLED_PARAM, ONLY: COU_PSF,COU_PSFFNAME,COU_NPX,COU_NPY&
         ,COU_SPSF,COU_NSPX,COU_NSPY,COU_PSFRAD
    USE MISC, ONLY: READ_PSF_2D
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_2D_DP
    USE user_mpi, ONLY: mpi__myrank
    USE GRID_PARAM, ONLY: NX, NY
    USE CODE_MODES, ONLY: DATAPATH
    !
    REAL(SP), DIMENSION(:,:), ALLOCATABLE :: PSF
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: DPSF
    INTEGER, DIMENSION(:), ALLOCATABLE    :: DIMS
    INTEGER                               :: X0, Y0, XS, YS
    INTEGER                               :: XI0, YI0
    !
    ! Read PSF file: both slaves and master do this step
    CALL READ_PSF_2D(TRIM(DATAPATH)//COU_PSFFNAME, PSF)
    CALL ALLOCATE_1D_IP(DIMS, SIZE(SHAPE(PSF)), '')
    DIMS(:)=SHAPE(PSF)
    !
    ! Master, has to store the PSF in a different variable to ...
    ! ... give it to MvN's routines:
    IF (mpi__myrank.EQ.0) THEN
      COU_NPY=DIMS(1)
      COU_NPX=DIMS(2)
      CALL ALLOCATE_2D_DP(COU_PSF, COU_NPY, COU_NPX, '')
      COU_PSF(:,:)=PSF(:,:)
      !DEALLOCATE(DIMS)
      !DEALLOCATE(PSF)
      ! ALLOCATE SMALL PSF:
      COU_NSPX=COU_PSFRAD*2+1
      COU_NSPY=COU_PSFRAD*2+1
      CALL ALLOCATE_2D_DP(COU_SPSF, COU_NSPY, COU_NSPX, '')
      X0=COU_NPX/2+1
      Y0=COU_NPY/2+1
      COU_SPSF(:,:)=COU_PSF(Y0-COU_PSFRAD:Y0+COU_PSFRAD &
          ,X0-COU_PSFRAD:X0+COU_PSFRAD)
    ENDIF
!    ELSE ! Slaves:
! Temporarily, all:
      ! Slaves can calculate the plan for the direct and inverse ...
      ! ... transformation and also, they can calculate the discrete ...
      ! ...fourier transform of the psf:
      ! Internal PSF
      CALL ALLOCATE_2D_DP(DPSF, NY, NX, '')
      !X0=DIMS(2)/2+1
      !Y0=DIMS(1)/2+1
      X0=NX/2+1
      Y0=NY/2+1
      !DPSF(:,:)=PSF(Y0-NY/2:Y0+NY/2,X0-NX/2:X0+NX/2)
      XS=MINVAL((/DIMS(2)/2,NX/2/))
      YS=MINVAL((/DIMS(1)/2,NY/2/))
      !
      ! Supplied PSF
      XI0=DIMS(2)/2+1
      YI0=DIMS(1)/2+1
!PRINT*, X0, XS, Y0, YS
!PRINT*, Y0-YS,Y0+YS,X0-XS,X0+XS
!PRINT*, SHAPE(DPSF), SHAPE(PSF(YI0-YS:YI0+YS,XI0-XS:XI0+XS))
!PRINT*, SUM(PSF(YI0-YS:YI0+YS,XI0-XS:XI0+XS))
      DPSF(Y0-YS:Y0+YS,X0-XS:X0+XS)=PSF(YI0-YS:YI0+YS,XI0-XS:XI0+XS)
      !
      CALL SPATIAL_CONV_INIT(DPSF)
      !
!    ENDIF
    !
  END SUBROUTINE READ_PSF
  !
  !------------------------------------------------
  !
  SUBROUTINE SPATIAL_CONV_INIT(DPSF)
    !
    USE USER_FFTW3, ONLY: FFTK2D, FFTY2D, GPLAN2D &
        , MAKE_2D_PLAN, MAKE_2D_IPLAN, FFTW_EXECUTE_DFT_R2C
    USE GRID_PARAM, ONLY: NX, NY
    !
    REAL(DP), DIMENSION(NY,NX),INTENT(INOUT)  :: DPSF
    !
    INTEGER       :: NY2, NX2
    !
    NY2=NY/2+1
    NX2=NX/2+1
    !
    DPSF=CSHIFT(DPSF,SHIFT=NY2-1,DIM=1)
    DPSF=CSHIFT(DPSF,SHIFT=NX2-1,DIM=2)
    ALLOCATE(GPLAN2D(2))
    ALLOCATE(FFTK2D(NY2,NX))
    ALLOCATE(FFTY2D(NY2,NX))
    !
    DPSF=DPSF/DBLE(NX*NY) ! Normalization factor not to be applied afterwards
    !
    ! Define a forward plan:
    CALL MAKE_2D_PLAN(GPLAN2D(1), SHAPE(DPSF), DPSF, SHAPE(FFTK2D), FFTK2D)
    ! Define an inverse plan:
    CALL MAKE_2D_IPLAN(GPLAN2D(2), SHAPE(FFTK2D), FFTK2D, SHAPE(DPSF), DPSF)
    !
    ! Store kernel transformation:
    CALL FFTW_EXECUTE_DFT_R2C(GPLAN2D(1),DPSF,FFTK2D)
    !
  END SUBROUTINE SPATIAL_CONV_INIT
  !
  !------------------------------------------------
  !
  SUBROUTINE WAVE_INIT()
    !
    USE FORWARD_PARAM, ONLY: PIXEL_INI, PIXEL_END, WAVE, INDEX &
        , IND_LINE, NUMW, NUMWAVE, NUML, LINE_L0, DELTA_LAMBDA, WAVE_INI
    !
    INTEGER                                    :: I, J
    !
    ! Total number of wavelenghs
    NUMW=SUM(NUMWAVE)
    ! Allocate localizers
    ALLOCATE( PIXEL_INI(NUML), PIXEL_END(NUML) )
    PIXEL_INI(1)=1
    PIXEL_END(1)=NUMWAVE(1)
    DO I=2,NUML
       PIXEL_INI(I)=PIXEL_END(I-1)+1
       PIXEL_END(I)=PIXEL_END(I-1)+NUMWAVE(I)
    ENDDO
    ! Allocate and create wavelength array
    ALLOCATE(WAVE(NUMW))
    ALLOCATE(INDEX(NUMW))
    DO I=1,NUML
       DO J=PIXEL_INI(I),PIXEL_END(I)
          WAVE(J)=WAVE_INI(I)+DELTA_LAMBDA(I)*(J-PIXEL_INI(I))
          INDEX(J)=LINE_L0(IND_LINE(I))
       ENDDO
    ENDDO
    ! Message in log file
    !CALL LOGW_CHAR(25,'Wavelength array created.')
    !------------------------------------------------
  END SUBROUTINE WAVE_INIT
  !
  !------------------------------------------------
  !
  SUBROUTINE SET_INVERSION_CONSTANTS(SLAVE)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_SP, ALLOCATE_1D_DP
    USE GRID_PARAM, ONLY: NZ
    USE INVERT_PARAM, ONLY: NFREEP, NFREEV, MAXITER, CYCPOW, INFREEPMAX &
        , YTOFIT, YGUESS, ISIGMAP, MPERTURBATION, CURIC, NFREQ, TOFFSET
    !
    INTEGER, OPTIONAL   :: SLAVE
    !
    CURIC=1
    !
 !   CALL FREE_VAR_SPACE(NZ)
    !
    IF (NFREEP.EQ.NFREEV) THEN
      MAXITER=1
    ELSE
      !MAXITER=CEILING(DLOG(DBLE(INFREEPMAX))/DLOG(DBLE(CYCPOW)))+1
      MAXITER=FLOOR(DLOG(DBLE(INFREEPMAX))/DLOG(DBLE(CYCPOW)))
      ! To allow some cycles in which T preceeds the othres, we...
      ! ...add a few additional cycles:
TOFFSET=1
      !MAXITER=MAXITER+TOFFSET
    ENDIF
    !
    CALL ALLOCATE_1D_SP(MPERTURBATION, NFREEV*NZ, 'MPERTURBATION')
    !
    IF (PRESENT(SLAVE).EQV..TRUE.) THEN
      CALL ALLOCATE_1D_SP(YTOFIT,NFREQ,'YTOFIT')
      CALL ALLOCATE_1D_SP(YGUESS,NFREQ,'YGUESS')
      CALL ALLOCATE_1D_DP(ISIGMAP,NFREQ,'ISIGMAP')
    ENDIF
    !
  END SUBROUTINE SET_INVERSION_CONSTANTS
  !
  !------------------------------------------------
  !
  SUBROUTINE REFERENCE_LSF_INIT()
    !
    USE USER_FFTW3, ONLY: GPLAN1D, FFTINI, FFTEND, FFTK1D, FFTY1D &
        , FFTW_EXECUTE_DFT_R2C, MAKE_1D_IPLAN, MAKE_1D_PLAN
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_1D_DP
    USE CODE_MODES, ONLY: MLSF
    USE FORWARD_PARAM, ONLY: PIXEL_INI, PIXEL_END, WAVE, LSF_KERNEL &
        , NUML, NUMW, LSF_SIGMA, LSF_W0, LSF_FILE
    !
    INTEGER    :: I, J, SKI
    INTEGER    :: WI, WF
    INTEGER    :: TOTSIZ, FFTSIZ
    !
    REAL(DP),ALLOCATABLE,DIMENSION(:)   :: XK, YK
    REAL(DP)    :: WOFF
    !
    IF (MLSF.EQV..TRUE.) THEN
      !
      !PRINT*, MLSF, NUML, WAVE
      !PRINT*, PIXEL_INI
      !PRINT*, PIXEL_END
      !
      ! First: estimate transformed dimensions and ...
      ! ... allocate various fftf and plans:
      !
      ALLOCATE(GPLAN1D(2,NUML))
      CALL ALLOCATE_1D_IP(FFTINI, NUML, 'FFTSTR')
      CALL ALLOCATE_1D_IP(FFTEND, NUML, 'FFTEND')
      !
      FFTINI(1)=1
      !
      TOTSIZ=0
      DO I=1,NUML
        !
        WI=PIXEL_INI(I)
        WF=PIXEL_END(I)
        SKI=WF-WI+1
        !
        FFTSIZ=SKI/2+1
        TOTSIZ=TOTSIZ+FFTSIZ
        !
        FFTEND(I)=FFTINI(I)+FFTSIZ-1
        IF (I.LT.NUML) FFTINI(I+1)=FFTEND(I)+1
        !
      ENDDO
      !
      ALLOCATE(FFTK1D(TOTSIZ))
      ALLOCATE(FFTY1D(TOTSIZ))
      FFTK1D(:)=0.0D0
      FFTY1D(:)=0.0D0
      !
      ! Second, Build real kernels:
      CALL ALLOCATE_1D_DP(LSF_KERNEL, NUMW, 'LSF_KERNEL')
      !
      DO I=1,NUML
        !
        ! BUILD AN X ARRAY
        IF (ALLOCATED(XK)) DEALLOCATE(XK)
        IF (ALLOCATED(YK)) DEALLOCATE(YK)
        !
        WI=PIXEL_INI(I)
        WF=PIXEL_END(I)
        SKI=WF-WI+1
        !
        CALL ALLOCATE_1D_DP(XK, SKI, 'XK')
        CALL ALLOCATE_1D_DP(YK, SKI, 'YK')
        !
        DO J=1,SKI
          XK(J)=DBLE(J)
        ENDDO
        ! Center it
        XK=XK-SUM(XK)/DBLE(SKI)
        ! ... in pixels
        XK=XK*(WAVE(WI+1)-WAVE(WI))
        ! ... in mA (assuming constant step)
        ! Deal with even-ity
        WOFF=0.D0
        IF(MOD(SKI,2).EQ.0) WOFF=(XK(2)-XK(1))/2.0D0
        !
        YK=DEXP(-((XK-LSF_W0(I)-WOFF)**2/(2.0D0*LSF_SIGMA(I)**2)))
        YK=YK/SUM(YK)
        YK=CSHIFT(YK,SHIFT=SKI/2,DIM=1)
        !
        FFTSIZ=FFTEND(I)-FFTINI(I)+1
        LSF_KERNEL(WI:WF)=YK/DBLE(SIZE(YK))
        !
        ! Direct plan
        CALL MAKE_1D_PLAN(GPLAN1D(1,I), SKI, LSF_KERNEL(WI:WF) &
            , FFTSIZ, FFTK1D(FFTINI(I):FFTEND(I)))
        ! Inverse plan
        CALL MAKE_1D_IPLAN(GPLAN1D(2,I), FFTSIZ, FFTK1D(FFTINI(I):FFTEND(I)) &
            , SKI, LSF_KERNEL(WI:WF))
        !
        ! Since kernel is always the same, store its dft already:
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN1D(1,I),LSF_KERNEL(WI:WF) &
            , FFTK1D(FFTINI(I):FFTEND(I)))
        !
        IF (ALLOCATED(XK)) DEALLOCATE(XK)
        IF (ALLOCATED(YK)) DEALLOCATE(YK)
        !
      ENDDO
      !
      IF (ALLOCATED(XK)) DEALLOCATE(XK)
      IF (ALLOCATED(YK)) DEALLOCATE(YK)
    ENDIF
    !
  END SUBROUTINE REFERENCE_LSF_INIT  !
  !------------------------------------------------
  !
  SUBROUTINE LSF_INIT()
    !
    USE USER_FFTW3, ONLY: GPLAN1D, FFTINI, FFTEND, FFTK1D, FFTY1D &
        , FFTW_EXECUTE_DFT_R2C, MAKE_1D_IPLAN, MAKE_1D_PLAN
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_IP, ALLOCATE_1D_DP
    USE CODE_MODES, ONLY: MLSF, DATAPATH
    USE FORWARD_PARAM, ONLY: PIXEL_INI, PIXEL_END, WAVE, LSF_KERNEL &
        , NUML, NUMW, LSF_SIGMA, LSF_W0, LSF_FILE, LSF_VALID
    USE MISC, ONLY: READ_LSF_FILE
    USE USER_MPI, ONLY: mpi__myrank
    !
    INTEGER    :: I, J, SKI
    INTEGER    :: WI, WF
    INTEGER    :: TOTSIZ, FFTSIZ
    !
    REAL(DP),ALLOCATABLE,DIMENSION(:)   :: XK, YK
    REAL(DP)    :: WOFF
    !
    IF (MLSF.EQV..TRUE.) THEN
      ALLOCATE(LSF_VALID(NUML))
      LSF_VALID(:)=.FALSE.
      !
      !PRINT*, MLSF, NUML, WAVE
      !PRINT*, PIXEL_INI
      !PRINT*, PIXEL_END
      !
      ! First: estimate transformed dimensions and ...
      ! ... allocate various fftf and plans:
      !
      ALLOCATE(GPLAN1D(2,NUML))
      CALL ALLOCATE_1D_IP(FFTINI, NUML, 'FFTSTR')
      CALL ALLOCATE_1D_IP(FFTEND, NUML, 'FFTEND')
      !
      FFTINI(1)=1
      !
      TOTSIZ=0
      DO I=1,NUML
        !
        WI=PIXEL_INI(I)
        WF=PIXEL_END(I)
        SKI=WF-WI+1
        !
        FFTSIZ=SKI/2+1
        TOTSIZ=TOTSIZ+FFTSIZ
        !
        FFTEND(I)=FFTINI(I)+FFTSIZ-1
        IF (I.LT.NUML) FFTINI(I+1)=FFTEND(I)+1
        !
      ENDDO
      !
      ALLOCATE(FFTK1D(TOTSIZ))
      ALLOCATE(FFTY1D(TOTSIZ))
      FFTK1D(:)=0.0D0
      FFTY1D(:)=0.0D0
      !
      ! Second, Build real kernels:
      CALL ALLOCATE_1D_DP(LSF_KERNEL, NUMW, 'LSF_KERNEL')
      !
      DO I=1,NUML
        !
        ! BUILD AN X ARRAY
        IF (ALLOCATED(XK)) DEALLOCATE(XK)
        IF (ALLOCATED(YK)) DEALLOCATE(YK)
        !
        WI=PIXEL_INI(I)
        WF=PIXEL_END(I)
        SKI=WF-WI+1
        !
        CALL ALLOCATE_1D_DP(YK, SKI, 'YK')

IF ( (LSF_SIGMA(I).NE.-1) .AND. (LSF_W0(I).NE.-1) ) THEN
        CALL ALLOCATE_1D_DP(XK, SKI, 'XK')
        !
        DO J=1,SKI
          XK(J)=DBLE(J)
        ENDDO
        ! Center it
        XK=XK-SUM(XK)/DBLE(SKI)
        ! ... in pixels
        XK=XK*(WAVE(WI+1)-WAVE(WI))
        ! ... in mA (assuming constant step)
        ! Deal with even-ity
        WOFF=0.D0
        IF(MOD(SKI,2).EQ.0) WOFF=(XK(2)-XK(1))/2.0D0
        !
        YK=DEXP(-((XK-LSF_W0(I)-WOFF)**2/(2.0D0*LSF_SIGMA(I)**2)))
        YK=YK/SUM(YK)
        YK=CSHIFT(YK,SHIFT=SKI/2,DIM=1)
        !
ELSE

  IF (LEN(TRIM(LSF_FILE(I))).NE.0) THEN
        CALL READ_LSF_FILE(TRIM(DATAPATH)//"/./"//TRIM(LSF_FILE(I)), SKI, YK, I)
  ENDIF

ENDIF

IF (ABS(SUM(YK)-1.0D0).LT.0.01D0) THEN
if (mpi__myrank.eq.1) PRINT*, ' Valid LSF for spectral region: ', I
LSF_VALID(I)=.TRUE.
        YK=YK/SUM(YK)
        FFTSIZ=FFTEND(I)-FFTINI(I)+1
        !
        LSF_KERNEL(WI:WF)=YK/DBLE(SIZE(YK))
        !
        ! Direct plan
        CALL MAKE_1D_PLAN(GPLAN1D(1,I), SKI, LSF_KERNEL(WI:WF) &
            , FFTSIZ, FFTK1D(FFTINI(I):FFTEND(I)))
        ! Inverse plan
        CALL MAKE_1D_IPLAN(GPLAN1D(2,I), FFTSIZ, FFTK1D(FFTINI(I):FFTEND(I)) &
            , SKI, LSF_KERNEL(WI:WF))
        !
        ! Since kernel is always the same, store its dft already:
        CALL FFTW_EXECUTE_DFT_R2C(GPLAN1D(1,I),LSF_KERNEL(WI:WF) &
            , FFTK1D(FFTINI(I):FFTEND(I)))
ELSE
if (mpi__myrank.eq.1) PRINT*, ' Spectral region: ', I, ' has not LSF'
ENDIF
        !
        IF (ALLOCATED(XK)) DEALLOCATE(XK)
        IF (ALLOCATED(YK)) DEALLOCATE(YK)
        !
      ENDDO
      !
      IF (ALLOCATED(XK)) DEALLOCATE(XK)
      IF (ALLOCATED(YK)) DEALLOCATE(YK)
    ENDIF
    !
  END SUBROUTINE LSF_INIT
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_PHYS_ARRAYS()
    !
    USE ABSORPTION_MATRIX, ONLY: ABS_MAT_INIT
    USE CODE_MODES, ONLY: MINVERSION
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_SP, ALLOCATE_2D_SP &
        , ALLOCATE_3D_SP, ALLOCATE_2D_DP, ALLOCATE_1D_DP &
        , ALLOCATE_1D_IP, ALLOCATE_2D_IP, ALLOCATE_4D_DP &
        , ALLOCATE_4D_SP
    USE PHYS_PARAM
    USE user_mpi, ONLY: mpi__myrank
    USE CODE_MODES, ONLY: SAVE_RFS, MTAULIN
    USE GRID_PARAM, ONLY: XX, YY, ZZ, NX, NY, NZ, DX, DY, DZ
    USE INVERT_PARAM, ONLY: AM_I_DONE, IMASK, MAXITER, NFREQ, INV_ATMPAR, NFREEV
    USE FORWARD_PARAM, ONLY: ATM_ARGS, BLENDSMAX, NUML, DER_ARGS &
        , GAIN1D, KC, KC5, KLIN, NUMW, KLINTAU, TAULIN, FULL_STOKES
    !
    INTEGER           :: I
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      MAXITER=0
      !
      CALL ALLOCATE_1D_DP(XX,NX,'XX')
      CALL ALLOCATE_1D_DP(YY,NY,'YY')
      CALL ALLOCATE_1D_DP(ZZ,NZ,'ZZ')
      !
      CALL ALLOCATE_2D_SP(MODEL1D_SND,NZ,ATM_ARGS,'MODEL1D_SND')

      CALL ALLOCATE_3D_SP(MODEL2D_SND,NZ,ATM_ARGS,2,'MODEL2D_SND')
      !
      CALL ALLOCATE_3D_SP(TEM3D,NZ,NY,NX,'TEM3D')
      CALL ALLOCATE_3D_SP(BX3D,NZ,NY,NX,'BX3D')
      CALL ALLOCATE_3D_SP(BY3D,NZ,NY,NX,'BY3D')
      CALL ALLOCATE_3D_SP(BZ3D,NZ,NY,NX,'BZ3D')
      CALL ALLOCATE_3D_SP(VX3D,NZ,NY,NX,'VX3D')
      CALL ALLOCATE_3D_SP(VY3D,NZ,NY,NX,'VY3D')
      CALL ALLOCATE_3D_SP(VZ3D,NZ,NY,NX,'VZ3D')
      CALL ALLOCATE_3D_SP(PG3D,NZ,NY,NX,'PG3D')
      CALL ALLOCATE_3D_SP(RHO3D,NZ,NY,NX,'RHO3D')
      CALL ALLOCATE_3D_SP(PEL3D,NZ,NY,NX,'PEL3D')
      CALL ALLOCATE_3D_SP(MW3D,NZ,NY,NX,'MW3D')
      CALL ALLOCATE_3D_SP(TAU3D5,NZ,NY,NX,'TAU3D5')
      !
      CALL ALLOCATE_3D_SP(ETEM3D,NZ,NY,NX,'ETEM3D')
      CALL ALLOCATE_3D_SP(EBX3D,NZ,NY,NX,'EBX3D')
      CALL ALLOCATE_3D_SP(EBY3D,NZ,NY,NX,'EBY3D')
      CALL ALLOCATE_3D_SP(EBZ3D,NZ,NY,NX,'EBZ3D')
      CALL ALLOCATE_3D_SP(EVX3D,NZ,NY,NX,'EVX3D')
      CALL ALLOCATE_3D_SP(EVY3D,NZ,NY,NX,'EVY3D')
      CALL ALLOCATE_3D_SP(EVZ3D,NZ,NY,NX,'EVZ3D')
      CALL ALLOCATE_3D_SP(EPG3D,NZ,NY,NX,'EPG3D')
      CALL ALLOCATE_3D_SP(ERHO3D,NZ,NY,NX,'ERHO3D')
      !CALL ALLOCATE_3D_SP(EPEL3D,NZ,NY,NX,'PEL3D')
      !CALL ALLOCATE_3D_SP(EMW3D,NZ,NY,NX,'MW3D')
      !CALL ALLOCATE_3D_SP(ETAU3D5,NZ,NY,NX,'TAU3D5')
      !
      !
!>< NOT USEFUL      ! Geometrical arrays
!>< NOT USEFUL      DO I=1,NX
!>< NOT USEFUL         XX(I)=DBLE(DBLE(I-1)*DX)
!>< NOT USEFUL      ENDDO
!>< NOT USEFUL      DO I=1,NYDBLE
!>< NOT USEFUL         YY(I)=DBLE(DBLE(I-1)*DY)
!>< NOT USEFUL      ENDDO
!>< NOT USEFUL      DO I=1,NZDBLE
!>< NOT USEFUL         ZZ(I)=DBLE(DBLE(I-1)*DZ)
!>< NOT USEFUL      ENDDO
      !
      CALL ALLOCATE_2D_SP(IMASK,NY,NX,'IMASK')
      !
      IF (MINVERSION.EQV..TRUE.) THEN
        CALL ALLOCATE_3D_SP(BEST_TEM3D,NZ,NY,NX,'BEST_TEM3D')
        CALL ALLOCATE_3D_SP(BEST_BX3D,NZ,NY,NX,'BEST_BX3D')
        CALL ALLOCATE_3D_SP(BEST_BY3D,NZ,NY,NX,'BEST_BY3D')
        CALL ALLOCATE_3D_SP(BEST_BZ3D,NZ,NY,NX,'BEST_BZ3D')
        CALL ALLOCATE_3D_SP(BEST_VX3D,NZ,NY,NX,'BEST_VX3D')
        CALL ALLOCATE_3D_SP(BEST_VY3D,NZ,NY,NX,'BEST_VY3D')
        CALL ALLOCATE_3D_SP(BEST_VZ3D,NZ,NY,NX,'BEST_VZ3D')
        CALL ALLOCATE_3D_SP(BEST_PG3D,NZ,NY,NX,'BEST_PG3D')
        CALL ALLOCATE_3D_SP(BEST_RHO3D,NZ,NY,NX,'BEST_RHO3D')
        CALL ALLOCATE_3D_SP(BEST_PEL3D,NZ,NY,NX,'BEST_PEL3D')
        CALL ALLOCATE_3D_SP(BEST_MW3D,NZ,NY,NX,'BEST_MW3D')
        CALL ALLOCATE_3D_SP(BEST_TAU3D5,NZ,NY,NX,'BEST_TAU3D5')
      ENDIF
      !
      IF (MTAULIN.EQV..TRUE.) THEN
        CALL ALLOCATE_4D_SP(TAU3DLIN,NUMW,NZ,NY,NX,'TAU3DLIN')
      ENDIF
    ELSE ! Master. Slaves:
      !
      CALL ALLOCATE_1D_DP(ZZ,NZ,'ZZ')
      CALL ALLOCATE_2D_SP(MODEL1D_RCV,NZ,ATM_ARGS,'MODEL1D_RCV')
      CALL ALLOCATE_3D_SP(MODEL2D_RCV,NZ,ATM_ARGS,2,'MODEL2D_RCV')
      CALL ALLOCATE_2D_DP(SYN1D,4,NUMW,'SYN1D')
      IF (FULL_STOKES.EQV..TRUE.) THEN
        CALL ALLOCATE_3D_SP(FS_SYN1D,4,NUMW,NZ,'SYN1D')
      ENDIF
      CALL ALLOCATE_1D_DP(RSYN1D,NFREQ,'SYN1D')
      CALL ALLOCATE_1D_DP(BEST_RSYN1D,NFREQ,'SYN1D')
      CALL ALLOCATE_1D_DP(ROBS1D,NFREQ,'SYN1D')
      ! This array has always the same size, that is why it is allocated...
      ! ... here:
      CALL ALLOCATE_4D_DP(DSYN1D,4,ATM_ARGS,NUMW,NZ,'DSYN1D')
      CALL ALLOCATE_4D_DP(EVOLG,4,4,NUMW,NZ,'EVOLG')
      ! Identity matrix
      DO I=1,4
         IMAT(I,I)=1.0D0
      ENDDO
  !    ! Geometrical array
  !!    DO I=1,NZ
  !!       ZZ(I)=REAL((I-1)*DZ)
  !!    ENDDO
      ! Allocate elements of the absorption matrix
      CALL ABS_MAT_INIT()
      ! Allocate continuum- and line- absorption coefficients
      CALL ALLOCATE_1D_DP(KC5,NZ,'KC5')
      CALL ALLOCATE_1D_DP(KC,NUML,'KC')
      CALL ALLOCATE_2D_DP(KLIN,NUML,BLENDSMAX+1,'KLIN')
      !
      CALL ALLOCATE_2D_SP(GAIN1D,4,NUMW,'GAIN1D')
      !
      IF (MTAULIN.EQV..TRUE.) THEN
        CALL ALLOCATE_2D_DP(KLINTAU,NUMW,NZ,'KLINTAU')
        CALL ALLOCATE_2D_SP(TAULIN,NUMW,NZ,'KLINTAU')
      ENDIF
!><      CALL ALLOCATE_1D_IP(CALC_RFSP, ATM_ARGS, 'CALC_RFSP')
!><      IF (INV_ATMPAR(1).EQV..TRUE.) CALC_RFSP(1)=1
!><      IF (INV_ATMPAR(2).EQV..TRUE.) CALC_RFSP(2)=1
!><      IF (INV_ATMPAR(3).EQV..TRUE.) CALC_RFSP(3)=1
!><      IF (INV_ATMPAR(4).EQV..TRUE.) CALC_RFSP(4)=1
!><      IF (INV_ATMPAR(5).EQV..TRUE.) CALC_RFSP(5)=1
!><      IF (INV_ATMPAR(6).EQV..TRUE.) CALC_RFSP(6)=1
!><      IF (INV_ATMPAR(7).EQV..TRUE.) CALC_RFSP(7)=1
!><      IF (INV_ATMPAR(8).EQV..TRUE.) CALC_RFSP(8)=1
!><      IF (SAVE_RFS.EQV..TRUE.) CALC_RFSP(:)=CALC_RFSP(:)*0+1
!><      IF (SAVE_RFS.EQV..TRUE.) INV_ATMPAR(:)=.TRUE.
      !
      CALL LSF_INIT()
    ENDIF ! Slave
    CALL ALLOCATE_2D_IP(AM_I_DONE,NY,NX,'AM_I_DONE')
    !
  END SUBROUTINE ALLOCATE_PHYS_ARRAYS
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_STOKES(MINV)
    !
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_3D_DP, ALLOCATE_4D_SP, ALLOCATE_5D_SP
    USE PHYS_PARAM
    USE GRID_PARAM, ONLY: NX, NY, NZ
    USE INVERT_PARAM, ONLY: NFREQ
    USE FORWARD_PARAM, ONLY: FULL_STOKES, NUMW
    !
    LOGICAL, INTENT(IN)   :: MINV
    !
    PRINT*, '---'
    PRINT*, 'ALLOCATING STOKES VARS'!SYNTHETIC, OBSERVED AND RF ARRAYS'
    PRINT*, '(This might take a while)'
    !
    !
    !----------------------------------------
    ! If Non-LTE
    ! CALL ALLOCATE_5D_SP(SYN,NX,NY,NZ,NUMW,4,'SYN')
    !
    ! apy: Allow whole Stokes: (only for forward solver, right now:)
    IF (FULL_STOKES.EQV..TRUE.) THEN
      CALL ALLOCATE_5D_SP(SYN5D,4,NUMW,NZ,NY,NX,'SYN5D')
    ENDIF
    ! If LTE
    CALL ALLOCATE_3D_DP(SYN3D,NFREQ,NY,NX,'SYN')
    !----------------------------------------
!    IF (MRESPFUNCT.EQV..TRUE.) THEN
!      CALL ALLOCATE_6D_SP(DSYN,4,DER_ARGS,NUMW,NZ,NY,NX,'DSYN')
!    ENDIF
    !
    IF (MINV.EQV..TRUE.) THEN
      ! For reading only:
      CALL ALLOCATE_4D_SP(OBS,4,NUMW,NY,NX,'OBS')
      ! To use as storage
      CALL ALLOCATE_3D_DP(OBS3D,NFREQ,NY,NX,'OBS')
      !
      CALL ALLOCATE_3D_DP(BEST_SYN,NFREQ,NY,NX,'BETS_SYN')
      CALL ALLOCATE_3D_DP(ISIGMAP3D,NFREQ,NY,NX,'SYN')
      !
    ENDIF
    !
    PRINT*, ''
    !
  END SUBROUTINE ALLOCATE_STOKES
  !
  !------------------------------------------------
  !
  SUBROUTINE ALLOCATE_INVERSION()
    !
    USE user_mpi, ONLY: mpi__myrank
    USE ALLOCATE_UTILS, ONLY: ALLOCATE_1D_DP, ALLOCATE_3D_DP
    USE CODE_MODES, ONLY: MINVERSION
    USE GRID_PARAM, ONLY: NX, NY
    USE INVERT_PARAM, ONLY: INV_MAS, INV_SLA, MSGSIZE, NFREQ
    !
    !
    IF (MINVERSION.EQV..TRUE.) THEN
      !
      MSGSIZE=9+NFREQ
      IF (mpi__myrank.EQ.0) THEN
        ! Master:
        ! In addition, in inversion mode, we have to set some additional variables:
        CALL SET_INVERSION_CONSTANTS()
        !
        ! 1- STEPS_WITHOUT_IMPROVEMENT
        ! 2- CHICUR
        ! 3- CHIPRE
        ! 4- PRECYC
        ! 5- WSTKI ! MPERT
        ! 6- WSTKQ ! SVDTOL
        ! 7- WSTKU ! IFREEP
        ! 8- WSTKV ! NJEVALS
        ! 9- LAMBDA
        ! 10-() TOFIT
        ! 11-() YGUESS
        CALL ALLOCATE_3D_DP(INV_MAS,MSGSIZE,NY,NX,'INV_MAS')
        !
        ! ALL PIXELS START WITH DEFAULT VALUES:
        ! 1- STEPS_WITHOUT_IMPROVEMENT
        INV_MAS(1,:,:)=0
        ! 2- CHICUR
        INV_MAS(2,:,:)=1.E29
        ! 3- CHIPRE
        INV_MAS(3,:,:)=1.E29
        ! 4- PRECYC
        INV_MAS(4,:,:)=0
        ! 5->8: WSTK
        !><! 5- MPERT
        !><INV_MAS(5,:,:)=1.D-1
        !><! 6- SVDTOL
        !><INV_MAS(6,:,:)=4
        !><! 7- IFREEP (IT WILL BE GIVEN)
        !><INV_MAS(7,:,:)=0
        !><! 8- NJEVALS
        !><INV_MAS(8,:,:)=1
        ! 9- LAMBDA
        INV_MAS(9,:,:)=1.D-3
        !
        !><OBS><CALL ALLOCATE_3D_SP(DMODEL,NFREEV*NZ,NY,NX,'DMODEL')
      ELSE ! Master. Slaves:
        ! In addition, in inversion mode, we have to set some additional variables:
        CALL SET_INVERSION_CONSTANTS(1)
        !
        ! 1- STEPS_WITHOUT_IMPROVEMENT
        ! 2- CHICUR
        ! 3- CHIPRE
        ! 4- PRECYC
        ! 5- MPERT
        ! 6- SVDTOL
        ! 7- IFREEP
        ! 8- NJEVALS
        ! 9- LAMBDA
        ! 10-() TOFIT
        ! 11-() YGUESS
        CALL ALLOCATE_1D_DP(INV_SLA,MSGSIZE,'INV_SLA')
        !
      ENDIF ! Slaves.
    ENDIF ! Inversion.
    !
  END SUBROUTINE ALLOCATE_INVERSION
  !
  !------------------------------------------------
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
  !------------------------------------------------
  !
  SUBROUTINE CLOSE_WAVE_INIT()
    !
    USE FORWARD_PARAM, ONLY: WAVE, INDEX, PIXEL_INI, PIXEL_END
    !
    DEALLOCATE(PIXEL_INI)
    DEALLOCATE(PIXEL_END)
    DEALLOCATE(WAVE)
    DEALLOCATE(INDEX)
    !
  END SUBROUTINE CLOSE_WAVE_INIT
  !
  !------------------------------------------------
  !
  SUBROUTINE LSF_END()
    !
    USE USER_FFTW3, ONLY: GPLAN1D, FFTK1D, FFTY1D
    !
    !PRINT*, 'LSF_END'
    DEALLOCATE(GPLAN1D)
    DEALLOCATE(FFTK1D)
    DEALLOCATE(FFTY1D)
    !
  END SUBROUTINE LSF_END
  !
  !------------------------------------------------
  !
  SUBROUTINE SPATIAL_CONV_END()
    !
    USE USER_FFTW3, ONLY: GPLAN2D, FFTK2D, FFTY2D
    !
    DEALLOCATE(GPLAN2D)
    DEALLOCATE(FFTK2D)
    DEALLOCATE(FFTY2D)
    !
  END SUBROUTINE SPATIAL_CONV_END
  !
  !------------------------------------------------
  !
  SUBROUTINE DEALLOCATE_PHYS_ARRAYS()
    !
    USE user_mpi, ONLY: mpi__myrank
    USE ABSORPTION_MATRIX, ONLY: ABS_MAT_END
    USE PHYS_PARAM
    USE CODE_MODES, ONLY: MLSF, MTAULIN
    USE GRID_PARAM, ONLY: XX, YY, ZZ
    USE FORWARD_PARAM, ONLY: KC5, KC, KLIN, KLINTAU, TAULIN
    !
    IF (mpi__myrank.EQ.0) THEN
      ! Master:
      !
      DEALLOCATE(XX)
      DEALLOCATE(YY)
      DEALLOCATE(ZZ)
      !
      DEALLOCATE(MODEL1D_SND)
      DEALLOCATE(MODEL2D_SND)
      !
      DEALLOCATE(TEM3D)
      DEALLOCATE(BX3D)
      DEALLOCATE(BY3D)
      DEALLOCATE(BZ3D)
      DEALLOCATE(VX3D)
      DEALLOCATE(VY3D)
      DEALLOCATE(VZ3D)
      DEALLOCATE(PG3D)
      DEALLOCATE(RHO3D)
      DEALLOCATE(PEL3D)
      DEALLOCATE(MW3D)
      DEALLOCATE(TAU3D5)
      !
      DEALLOCATE(ETEM3D)
      DEALLOCATE(EBX3D)
      DEALLOCATE(EBY3D)
      DEALLOCATE(EBZ3D)
      DEALLOCATE(EVX3D)
      DEALLOCATE(EVY3D)
      DEALLOCATE(EVZ3D)
      DEALLOCATE(EPG3D)
      DEALLOCATE(ERHO3D)
      !DEALLOCATE(EPEL3D)
      !DEALLOCATE(EMW3D)
      !DEALLOCATE(ETAU3D5)
      !
      IF (ALLOCATED(BEST_TEM3D)) DEALLOCATE(BEST_TEM3D)
      IF (ALLOCATED(BEST_BX3D)) DEALLOCATE(BEST_BX3D)
      IF (ALLOCATED(BEST_BY3D)) DEALLOCATE(BEST_BY3D)
      IF (ALLOCATED(BEST_BZ3D)) DEALLOCATE(BEST_BZ3D)
      IF (ALLOCATED(BEST_VX3D)) DEALLOCATE(BEST_VX3D)
      IF (ALLOCATED(BEST_VY3D)) DEALLOCATE(BEST_VY3D)
      IF (ALLOCATED(BEST_VZ3D)) DEALLOCATE(BEST_VZ3D)
      IF (ALLOCATED(BEST_PG3D)) DEALLOCATE(BEST_PG3D)
      IF (ALLOCATED(BEST_RHO3D)) DEALLOCATE(BEST_RHO3D)
      IF (ALLOCATED(BEST_PEL3D)) DEALLOCATE(BEST_PEL3D)
      IF (ALLOCATED(BEST_MW3D)) DEALLOCATE(BEST_MW3D)
      IF (ALLOCATED(BEST_TAU3D5)) DEALLOCATE(BEST_TAU3D5)
      !
      IF (MTAULIN.EQV..TRUE.) THEN
        DEALLOCATE(TAU3DLIN)
      ENDIF
      !
    ELSE ! Master ends. Slaves start:
      ! Slaves:
      !
      DEALLOCATE(ZZ)
      DEALLOCATE(MODEL1D_RCV)
      DEALLOCATE(SYN1D)
      DEALLOCATE(RSYN1D)
      DEALLOCATE(DSYN1D)
      DEALLOCATE(EVOLG)
      ! Deallocate elements of the absorption matrix
      CALL ABS_MAT_END()
      ! Deallocate continuum- and line- absorption coefficients
      DEALLOCATE(KC5)
      DEALLOCATE(KC)
      DEALLOCATE(KLIN)
      !
      IF (MTAULIN.EQV..TRUE.) THEN
        DEALLOCATE(KLINTAU)
        DEALLOCATE(TAULIN)
      ENDIF
      !
      IF (MLSF.EQV..TRUE.) CALL LSF_END()
      !
    ENDIF ! Slaves end.
    !
  END SUBROUTINE DEALLOCATE_PHYS_ARRAYS
  !
  !------------------------------------------------
  !
  SUBROUTINE DEALLOCATE_STOKES()
    !
    USE PHYS_PARAM
    USE user_mpi, ONLY: mpi__myrank
    !
    IF (mpi__myrank.EQ.0) THEN
      !
      !----------------------------------------
      IF (ALLOCATED(SYN3D)) DEALLOCATE(SYN3D)
      !----------------------------------------
      IF (ALLOCATED(OBS)) DEALLOCATE(OBS)
      IF (ALLOCATED(OBS3D)) DEALLOCATE(OBS3D)
      IF (ALLOCATED(DSYN)) DEALLOCATE(DSYN)
      IF (ALLOCATED(BEST_SYN)) DEALLOCATE(BEST_SYN)
      IF (ALLOCATED(BEST_DSYN)) DEALLOCATE(BEST_DSYN)
      IF (ALLOCATED(ISIGMAP3D)) DEALLOCATE(ISIGMAP3D)
      !
    ENDIF
    !
  END SUBROUTINE DEALLOCATE_STOKES
  !
  !------------------------------------------------
  !
  SUBROUTINE DEALLOCATE_INVERSION()
    !
    USE user_mpi, ONLY: mpi__myrank
    USE COUPLED_INVERSION, ONLY: END_COUPLED_INVERSION_STATIC_VARS
    USE CODE_MODES, ONLY: COUPLED, MINVERSION
    USE INVERT_PARAM, ONLY: INV_MAS, INV_SLA, IMASK, YTOFIT &
        , INV_SLA, AM_I_DONE, YGUESS
    !
    IF (MINVERSION.EQV..TRUE.) THEN
      !
      IF (mpi__myrank.EQ.0) THEN
        DEALLOCATE(INV_MAS)
        DEALLOCATE(IMASK)
      ELSE
        DEALLOCATE(INV_SLA)
        DEALLOCATE(YTOFIT)
        DEALLOCATE(YGUESS)
      ENDIF
      !
      DEALLOCATE(AM_I_DONE)
      !
    ENDIF
    IF (COUPLED.EQV..TRUE.) THEN
      CALL END_COUPLED_INVERSION_STATIC_VARS()
    ENDIF
    !
  END SUBROUTINE DEALLOCATE_INVERSION
  !
  !================================================
  !
END MODULE PRE_POST_DUTIES
!
