!
MODULE CODE_MODES
  !
  !================================================
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  CHARACTER*800, PUBLIC  :: NAMEPROFILE, NAMEMODEL
  CHARACTER*800, PUBLIC  :: INPUTFILE
  !
  CHARACTER*800, PUBLIC  :: DATAPATH
  CHARACTER*800, PUBLIC  :: MODLPATH
  CHARACTER*800, PUBLIC  :: OUTPPATH
  CHARACTER*800, PUBLIC  :: LINEPATH
  !
  LOGICAL, PUBLIC        :: BOX, LINES, SIR, SPINOR, MURAM
  !
  LOGICAL, PUBLIC        :: MSYNTHESIS, MINVERSION, MGETTAU, MGETHEQ
  LOGICAL, PUBLIC        :: MRESPFUNCT, SAVE_RFS
  LOGICAL, PUBLIC        :: HYDROSTATIC
  !
  LOGICAL, PUBLIC        :: HSRA_NORMALIZATION
  LOGICAL, PUBLIC        :: COUPLED
  !
  LOGICAL, PUBLIC        :: TWODSPATIALSPD
  LOGICAL, PUBLIC        :: MSMOOTHING
  !
  LOGICAL, PUBLIC        :: MLSF
  INTEGER, PUBLIC        :: MVERBOSE
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
END MODULE CODE_MODES
!
