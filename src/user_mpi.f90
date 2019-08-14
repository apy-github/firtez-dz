!
MODULE USER_MPI
  !
  !================================================
  !
  USE mpi
  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, PUBLIC :: mpi__ierror, mpi__size, mpi__myrank 
  !
  !MPI3  TYPE(MPI_REQUEST)   :: mpi__request2, mpi__request3
  !MPI3  TYPE(MPI_STATUS)    :: mpi__status
  INTEGER, PUBLIC,DIMENSION(MPI_STATUS_SIZE)    :: mpi__status
  !
  PUBLIC :: MPI_COMM_WORLD
  PUBLIC :: MPI_DOUBLE_PRECISION
  PUBLIC :: MPI_REAL
  PUBLIC :: MPI_INTEGER
  PUBLIC :: MPI_LOGICAL
  PUBLIC :: START_MPI, END_MPI, NICE_WAITING
  PUBLIC :: MPI_BCAST
  PRIVATE
  !
  !************************************************
  !
  CONTAINS
  !
  ! start_mpi
  ! nice_waiting
  !
  !------------------------------------------------
  !
  SUBROUTINE START_MPI()
  !
    CALL MPI_INIT(mpi__ierror)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi__size, mpi__ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi__myrank, mpi__ierror)
  !
  END SUBROUTINE START_MPI
  !
  !------------------------------------------------
  !
  SUBROUTINE NICE_WAITING(SOU, TAGN)
    !
    INTEGER, INTENT(IN)  :: SOU, TAGN
    !
    LOGICAL        :: READY  
    !
    READY=.FALSE.
    DO WHILE (READY.EQV..FALSE.)
      !CALL SLEEP(1)
#ifdef __INTEL_COMPILER
      CALL SYSTEM('sleep 0.1s')
#else
      CALL EXECUTE_COMMAND_LINE('sleep 0.1s')
#endif
      CALL MPI_Iprobe(SOU, TAGN,MPI_COMM_WORLD,READY,MPI__STATUS,mpi__ierror)
    ENDDO
    !
  END SUBROUTINE NICE_WAITING
  !
  !------------------------------------------------
  !
  SUBROUTINE END_MPI()
  !
    CALL MPI_FINALIZE(mpi__ierror)
  END SUBROUTINE END_MPI
  !
  !================================================
  !
END MODULE USER_MPI
!
