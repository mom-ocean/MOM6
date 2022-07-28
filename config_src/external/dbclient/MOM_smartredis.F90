!> Contains routines necessary to initialize the Database client
module MOM_dbclient

! This file is part of MOM6. See LICENSE.md for the license.
use MOM_file_parser,      only : param_file_type
use MOM_error_handler,    only : MOM_error, WARNING
use smartredis_client,    only : dbclient_type => client_type

implicit none; private

!> Control structure to store Database client related parameters and objects
type, public :: dbclient_CS_type
  type(dbclient_type) :: client !< The Database client itself
  logical           :: use_dbclient !< If True, use Database within MOM6
  logical           :: colocated !< If True, the orchestrator was setup in 'co-located' mode
  logical           :: cluster   !< If True, the orchestrator has three shards or more
  integer           :: colocated_stride !< Sets which ranks will load the model from the file
                                        !! e.g. mod(rank,colocated_stride) == 0
end type dbclient_CS_type

public :: dbclient_type
public :: dbclient_init

contains

subroutine dbclient_init(param_file, CS, client_in)
  type(param_file_type),       intent(in   ) :: param_file !< Parameter file structure
  type(dbclient_CS_type),    intent(inout) :: CS         !< Control structure for Database
  type(dbclient_type), optional, intent(in   ) :: client_in !< If present, use a previously initialized
                                                          !! Database client

  call MOM_error(WARNING,"dbclient_init was compiled using the dummy module. If this was\n"//&
                       "a mistake, please follow the instructions in:\n"//&
                       "MOM6/config_src/external/dbclient/README.md")
end subroutine dbclient_init

end module MOM_dbclient

