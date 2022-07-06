!> Contains routines necessary to initialize the SmartRedis client
module MOM_smartredis

use MOM_cpu_clock,        only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_error_handler,    only : MOM_error, FATAL, WARNING, NOTE, MOM_mesg, is_root_pe
use MOM_file_parser,      only : read_param, get_param, log_version, param_file_type
use smartredis_client,    only : client_type

implicit none; private

!> Control structure to store SmartRedis client related parameters and objects
type, public :: smartredis_CS_type
  type(client_type) :: client !< The SmartRedis client itself
  logical           :: use_smartredis !< If True, use SmartRedis within MOM6
  logical           :: colocated !< If True, the orchestrator was setup in 'co-located' mode
  logical           :: cluster   !< If True, the orchestrator has three shards or more
  integer           :: colocated_stride !< Sets which ranks will load the model from the file
                                        !! e.g. mod(rank,colocated_stride) == 0
end type

public :: client_type
public :: smartredis_init

contains

subroutine smartredis_init(param_file, CS, client_in)
  type(param_file_type),       intent(in   ) :: param_file !< Parameter file structure
  type(smartredis_CS_type),    intent(inout) :: CS         !< Control structure for SmartRedis
  type(client_type), optional, intent(in   ) :: client_in !< If present, use a previously initialized
                                                          !! SmartRedis client

  character(len=40) :: mdl = "MOM_SMARTREDIS"
  integer :: id_client_init
  integer :: return_code
  call get_param(param_file, mdl, "USE_SMARTREDIS",  CS%use_smartredis, &
                 "If true, use the data client to connect"//&
                 "with the SmartRedis database", default=.false.)

  if (present(client_in)) then ! The driver (e.g. the NUOPC cap) has already initialized the client

    CS%client = client_in

    if (.not. CS%client%isinitialized() .and. CS%use_smartredis) then
      call MOM_error(FATAL, &
      "If using a SmartRedis client not initialized within MOM, client%initialize must have already been invoked."//&
      " Check that the client has been initialized in the driver before the call to initialize_MOM")
    endif

  elseif (CS%use_smartredis) then ! The client will be initialized within MOM

    call get_param(param_file, mdl, "SMARTREDIS_COLOCATED",  CS%colocated, &
                   "If true, the SmartRedis database is colocated on the simulation nodes.",&
                   default=.false.)
    if (CS%colocated) then
      CS%cluster = .false.
      call get_param(param_file, mdl, "SMARTREDIS_COLOCATED_STRIDE",  CS%colocated_stride, &
                     "If true, the SmartRedis database is colocated on the simulation nodes.",&
                     default=0)
    else
      call get_param(param_file, mdl, "SMARTREDIS_CLUSTER",  CS%cluster, &
                     "If true, the SmartRedis database is distributed over multiple nodes.",&
                     default=.true.)
    endif
    id_client_init = cpu_clock_id('(SMARTREDIS client init)', grain=CLOCK_ROUTINE)
    call MOM_error(NOTE,"SmartRedis Client Initializing")
    call cpu_clock_begin(id_client_init)
    return_code = CS%client%initialize(CS%cluster)
    if (CS%client%SR_error_parser(return_code)) then
      call MOM_error(FATAL, "SmartRedis client failed to initialize")
    endif
    call MOM_error(NOTE,"SmartRedis Client Initialized")
    call cpu_clock_end(id_client_init)

  endif
end subroutine smartredis_init

end module MOM_smartredis

