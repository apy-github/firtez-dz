!
MODULE USER_MPI
  !
  !================================================
  !
  !USE mpi
  use mpi_f08

  !
  IMPLICIT NONE
  !
  !::::::::::::::::::::::::::::::::::::::::::::::::
  !
  INTEGER, PUBLIC :: mpi__ierror, mpi__size, mpi__myrank 
  !
  !MPI3  TYPE(MPI_REQUEST)   :: mpi__request2, mpi__request3
  !MPI3  TYPE(MPI_STATUS)    :: mpi__status
!!  INTEGER, PUBLIC,DIMENSION(MPI_STATUS_SIZE)    :: mpi__status
    type(mpi_request)   :: mpi__request
    type(mpi_request), dimension(:), allocatable   :: mpi__vrequest
    type(mpi_status)    :: mpi__status
  !
  PUBLIC :: MPI_COMM_WORLD
  PUBLIC :: MPI_DOUBLE_PRECISION
  PUBLIC :: MPI_REAL
  PUBLIC :: MPI_CHARACTER
  PUBLIC :: MPI_INTEGER
  PUBLIC :: MPI_LOGICAL
  PUBLIC :: START_MPI, END_MPI, NICE_WAITING
  PUBLIC :: keep_workers_waiting, mpi__status
  PUBLIC :: mpi__request
  PUBLIC :: mpi__vrequest
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

  subroutine slow_waiting(sou, tagn)
    !
    integer, intent(in)  :: sou, tagn
    !
    logical        :: ready
    !
    ready=.false.
    do while (ready.eqv..false.)
      !call sleep(1)
#ifdef __intel_compiler
      call system('sleep 0.01s')
#else
      call execute_command_line('sleep 0.01s')
#endif
      call mpi_iprobe(sou, tagn,mpi_comm_world,ready,mpi__status,mpi__ierror)
    enddo
    !
  end subroutine slow_waiting

  subroutine keep_workers_waiting()

    integer :: msg, i

    if (mpi__myrank.eq.0) then

      do i=1,mpi__size-1
        call mpi_isend(msg, 1, mpi_real, i, i, mpi_comm_world &
            , mpi__request,mpi__ierror)
      enddo
      call mpi_wait(mpi__request, mpi__status, mpi__ierror)

    else ! Master. Workers:

      call slow_waiting(0, mpi__myrank)
      call mpi_recv(msg,1,mpi_real,0,mpi__myrank &
          ,mpi_comm_world,mpi__status,mpi__ierror)
    endif

  end subroutine keep_workers_waiting


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
