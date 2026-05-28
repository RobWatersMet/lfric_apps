!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> This is a code that uses the LFRic infrastructure to build a model that
!> includes the GungHo dynamical core and physics parametrisation schemes
!> that are currently provided through the use of unified model code.

!> @brief Main program used to illustrate an atmospheric model built using
!>        LFRic infrastructure

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the atmospheric model.

program lfric_atm

  use cli_mod,                only: parse_command_line
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_counter_mod,     only: init_counters, final_counters
  use driver_log_mod,         only: init_logger, final_logger
  use driver_time_mod,        only: init_time, final_time
  use gungho_mod,             only: gungho_required_namelists
  use driver_modeldb_mod,     only: modeldb_type
  use gungho_driver_mod,      only: initialise, step, finalise
  use lfric_mpi_mod,          only: global_mpi
  use namelist_mod,           only: namelist_type
  use timing_mod,             only: init_timing, final_timing, &
                                    start_timing, stop_timing, &
                                    tik, LPROF
  use io_config_mod,          only: timer_output_path

  implicit none

  ! Model run working data set
  type(modeldb_type) :: modeldb

  character(*), parameter      :: application_name = "lfric_atm"
  character(:), allocatable    :: filename
  integer(tik)                 :: id_setup
  type(namelist_type), pointer :: io_nml

  logical :: lsubroutine_timers

  ! ── debug variables for final_comms issue diagnostic ───────────────────────
  integer            :: dbg_unit, dbg_rank
  character(len=64)  :: dbg_path

  call parse_command_line( filename )

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( application_name, &
                                         table_len=10 )
  call modeldb%values%initialise( 'values', 5 )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository", &
                                                 table_len = 100)
  call modeldb%fields%add_empty_field_collection("prognostic_fields", &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("diagnostic_fields", &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("lbc_fields",        &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("radiation_fields",  &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("ls_fields",         &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("fd_fields",         &
                                                  table_len = 100)

  call modeldb%io_contexts%initialise(application_name, 100)

  call init_comm( application_name, modeldb )

  call init_config( filename, gungho_required_namelists, &
                    modeldb%configuration )
  call init_logger( modeldb%mpi%get_comm(), application_name )

  io_nml => modeldb%configuration%get_namelist('io')
  call io_nml%get_value('subroutine_timers', lsubroutine_timers)
  call init_timing( modeldb%mpi%get_comm(), lsubroutine_timers, application_name, timer_output_path )
  nullify( io_nml )
  if ( LPROF ) call start_timing( id_setup, '__setup__' )


  call init_collections()
  call init_time( modeldb )
  call init_counters( application_name )
  deallocate( filename )

  call initialise( application_name, modeldb )
  if ( LPROF ) call stop_timing( id_setup, '__setup__' )
  do while (modeldb%clock%tick())
    call step( modeldb )
  end do
  call finalise( application_name, modeldb )

  !call final_counters( application_name )
  !call final_time( modeldb )
  !call final_collections()
  !call final_timing( application_name )
  !call final_logger( application_name )
  !call final_config()
  !call final_comm( modeldb )

  ! get rank number from mpi
  dbg_rank = modeldb%mpi%get_comm_rank()
  ! create dbg path name
  write(dbg_path,'("shutdown.",I0.5,".trace")') dbg_rank
  open(newunit=dbg_unit, file=trim(dbg_path), status='replace', &
       action='write', form='formatted')

  ! befor each shutdown stage call mark
  call mark(dbg_unit, 'A before final_counters')
  call final_counters( application_name )
  call mark(dbg_unit, 'B after  final_counters')
  call final_time( modeldb )
  call mark(dbg_unit, 'C after  final_time')
  call final_collections()
  call mark(dbg_unit, 'D after  final_collections')
  call final_timing( application_name )
  call mark(dbg_unit, 'E after  final_timing')
  call final_logger( application_name )
  call mark(dbg_unit, 'F after  final_logger')
  call final_config()
  call mark(dbg_unit, 'G after  final_config')
  ! final_comm is the major suspect - consider not calling it
  call final_comm( modeldb )
  call mark(dbg_unit, 'H after  final_comm (about to end program)')

contains

  subroutine mark(u, msg)
    ! write out the msg and the current memory usage to the debug log which has
    ! a unit of u
    integer,      intent(in) :: u
    character(*), intent(in) :: msg
    integer :: rss_kb
    rss_kb = vmrss_kb()
    write(u,'(A,"  | VmRSS=",I0," kB")') msg, rss_kb
    flush(u)
  end subroutine mark

  integer function vmrss_kb()
    ! get current Resident Set Size - physical RAM used
    integer :: u, ios
    character(len=256) :: line
    vmrss_kb = -1
    open(newunit=u, file='/proc/self/status', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:6) == 'VmRSS:') then
        read(line(7:),*,iostat=ios) vmrss_kb
        exit
      end if
    end do
    close(u)
  end function vmrss_kb

end program lfric_atm
