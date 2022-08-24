module database_client_interface

! This file is part of MOM6. See LICENSE.md for the license.
  use iso_fortran_env, only : int8, int16, int32, int64, real32, real64

  implicit none; private

  !> Dummy type for dataset
  type, public :: dataset_type
    private
  end type dataset_type

  !> Stores all data and methods associated with the communication client that is used to communicate with the database
  type, public :: dbclient_type
    private

    contains

    ! Public procedures
    !> Puts a tensor into the database for a variety of datatypes
    generic :: put_tensor => put_tensor_float_1d, put_tensor_float_2d, put_tensor_float_3d, put_tensor_float_4d, &
                             put_tensor_double_1d, put_tensor_double_2d, put_tensor_double_3d, put_tensor_double_4d, &
                             put_tensor_int32_1d, put_tensor_int32_2d, put_tensor_int32_3d, put_tensor_int32_4d
    !> Retrieve the tensor in the database into already allocated memory for a variety of datatypesm
    generic :: unpack_tensor => unpack_tensor_float_1d, unpack_tensor_float_2d, &
                                unpack_tensor_float_3d, unpack_tensor_float_4d, &
                                unpack_tensor_double_1d, unpack_tensor_double_2d, &
                                unpack_tensor_double_3d, unpack_tensor_double_4d, &
                                unpack_tensor_int32_1d, unpack_tensor_int32_2d, &
                                unpack_tensor_int32_3d, unpack_tensor_int32_4d

    !> Decode a response code from an API function
    procedure :: SR_error_parser
    !> Initializes a new instance of the communication client
    procedure :: initialize => initialize_client
    !> Check if a communication client has been initialized
    procedure :: isinitialized
    !> Destructs a new instance of the communication client
    procedure :: destructor
    !> Rename a tensor within the database
    procedure :: rename_tensor
    !> Delete a tensor from the database
    procedure :: delete_tensor
    !> Copy a tensor within the database to a new name
    procedure :: copy_tensor
    !> Set a model from a file
    procedure :: set_model_from_file
    !> Set a model from a file on a system with multiple GPUs
    procedure :: set_model_from_file_multigpu
    !> Set a model from a byte string that has been loaded within the application
    procedure :: set_model
    !> Set a model from a byte string that has been loaded within the application on a system with multiple GPUs
    procedure :: set_model_multigpu
    !> Retrieve the model as a byte string
    procedure :: get_model
    !> Set a script from a specified file
    procedure :: set_script_from_file
    !> Set a script from a specified file on a system with multiple GPUS
    procedure :: set_script_from_file_multigpu
    !> Set a script as a byte or text string
    procedure :: set_script
    !> Set a script as a byte or text string on a system with multiple GPUs
    procedure :: set_script_multigpu
    !> Retrieve the script from the database
    procedure :: get_script
    !> Run a script that has already been stored in the database
    procedure :: run_script
    !> Run a script that has already been stored in the database with multiple GPUs
    procedure :: run_script_multigpu
    !> Run a model that has already been stored in the database
    procedure :: run_model
    !> Run a model that has already been stored in the database with multiple GPUs
    procedure :: run_model_multigpu
    !> Remove a script from the database
    procedure :: delete_script
    !> Remove a script from the database with multiple GPUs
    procedure :: delete_script_multigpu
    !> Remove a model from the database
    procedure :: delete_model
    !> Remove a model from the database with multiple GPUs
    procedure :: delete_model_multigpu
    !> Put a communication dataset into the database
    procedure :: put_dataset
    !> Retrieve a communication dataset from the database
    procedure :: get_dataset
    !> Rename the dataset within the database
    procedure :: rename_dataset
    !> Copy a dataset stored in the database into another name
    procedure :: copy_dataset
    !> Delete the dataset from the database
    procedure :: delete_dataset

    ! Private procedures
    !> Put a 1d, 32-bit real tensor into database
    procedure, private :: put_tensor_float_1d
    !> Put a 2d, 32-bit real tensor into database
    procedure, private :: put_tensor_float_2d
    !> Put a 3d, 32-bit real tensor into database
    procedure, private :: put_tensor_float_3d
    !> Put a 4d, 32-bit real tensor into database
    procedure, private :: put_tensor_float_4d
    !> Put a 1d, 64-bit real tensor into database
    procedure, private :: put_tensor_double_1d
    !> Put a 2d, 64-bit real tensor into database
    procedure, private :: put_tensor_double_2d
    !> Put a 3d, 64-bit real tensor into database
    procedure, private :: put_tensor_double_3d
    !> Put a 4d, 64-bit real tensor into database
    procedure, private :: put_tensor_double_4d
    !> Put a 1d, 32-bit integer tensor into database
    procedure, private :: put_tensor_int32_1d
    !> Put a 2d, 32-bit integer tensor into database
    procedure, private :: put_tensor_int32_2d
    !> Put a 3d, 32-bit integer tensor into database
    procedure, private :: put_tensor_int32_3d
    !> Put a 4d, 32-bit integer tensor into database
    procedure, private :: put_tensor_int32_4d
    !> Unpack a 1d, 32-bit real tensor from the database
    procedure, private :: unpack_tensor_float_1d
    !> Unpack a 2d, 32-bit real tensor from the database
    procedure, private :: unpack_tensor_float_2d
    !> Unpack a 3d, 32-bit real tensor from the database
    procedure, private :: unpack_tensor_float_3d
    !> Unpack a 4d, 32-bit real tensor from the database
    procedure, private :: unpack_tensor_float_4d
    !> Unpack a 1d, 64-bit real tensor from the database
    procedure, private :: unpack_tensor_double_1d
    !> Unpack a 2d, 64-bit real tensor from the database
    procedure, private :: unpack_tensor_double_2d
    !> Unpack a 3d, 64-bit real tensor from the database
    procedure, private :: unpack_tensor_double_3d
    !> Unpack a 4d, 64-bit real tensor from the database
    procedure, private :: unpack_tensor_double_4d
    !> Unpack a 1d, 32-bit integer tensor from the database
    procedure, private :: unpack_tensor_int32_1d
    !> Unpack a 2d, 32-bit integer tensor from the database
    procedure, private :: unpack_tensor_int32_2d
    !> Unpack a 3d, 32-bit integer tensor from the database
    procedure, private :: unpack_tensor_int32_3d
    !> Unpack a 4d, 32-bit integer tensor from the database
    procedure, private :: unpack_tensor_int32_4d

  end type dbclient_type

  contains

  !> Decode a response code from an API function
  function SR_error_parser(self, response_code) result(is_error)
    class(dbclient_type),       intent(in) :: self    !< Receives the initialized client
    integer, intent(in) :: response_code !< The response code to decode
    logical                              :: is_error      !< Indicates whether this is an error response

    is_error = .true.
  end function SR_error_parser

  !> Initializes a new instance of a communication client
  function initialize_client(self, cluster)
    integer           :: initialize_client
    class(dbclient_type), intent(inout) :: self    !< Receives the initialized client
    logical, optional,  intent(in   ) :: cluster !< If true, client uses a database cluster (Default: .false.)

    initialize_client = -1
  end function initialize_client

  !> Check whether the client has been initialized
  logical function isinitialized(this)
    class(dbclient_type) :: this
    isinitialized = .false.
  end function isinitialized

  !> A destructor for the communication client
  function destructor(self)
    integer           :: destructor
    class(dbclient_type), intent(inout) :: self

    destructor = -1
  end function destructor

  !> Put a 32-bit real 1d tensor into the database
  function put_tensor_float_1d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:), intent(in) :: data !< Data to be sent
    class(dbclient_type),            intent(in) :: self !< Fortran communication client
    character(len=*),                intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),           intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_float_1d

  !> Put a 32-bit real 2d tensor into the database
  function put_tensor_float_2d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),              intent(in) :: self !< Fortran communication client
    character(len=*),                  intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),             intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_float_2d

  !> Put a 32-bit real 3d tensor into the database
  function put_tensor_float_3d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:,:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_float_3d

  !> Put a 32-bit real 4d tensor into the database
  function put_tensor_float_4d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:,:,:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_float_4d

  !> Put a 64-bit real 1d tensor into the database
  function put_tensor_double_1d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:), intent(in) :: data !< Data to be sent
    class(dbclient_type),            intent(in) :: self !< Fortran communication client
    character(len=*),                intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),           intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_double_1d

  !> Put a 64-bit real 2d tensor into the database
  function put_tensor_double_2d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),              intent(in) :: self !< Fortran communication client
    character(len=*),                  intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),             intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_double_2d

  !> Put a 64-bit real 3d tensor into the database
  function put_tensor_double_3d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:,:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_double_3d

  !> Put a 64-bit real 4d tensor into the database
  function put_tensor_double_4d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:,:,:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_double_4d

  !> Put a 32-bit integer 1d tensor into the database
  function put_tensor_int32_1d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:), intent(in) :: data !< Data to be sent
    class(dbclient_type),            intent(in) :: self !< Fortran communication client
    character(len=*),                intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),           intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_int32_1d

  !> Put a 32-bit integer 2d tensor into the database
  function put_tensor_int32_2d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),              intent(in) :: self !< Fortran communication client
    character(len=*),                  intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),             intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_int32_2d

  !> Put a 32-bit integer 3d tensor into the database
  function put_tensor_int32_3d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:,:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_int32_3d

  !> Put a 32-bit integer 4d tensor into the database
  function put_tensor_int32_4d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:,:,:,:), intent(in) :: data !< Data to be sent
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_int32_4d

  !> Unpack a 32-bit real 1d tensor from the database
  function unpack_tensor_float_1d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:), intent(  out) :: data !< Data to be received
    class(dbclient_type),            intent(in) :: self !< Fortran communication client
    character(len=*),                intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),           intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_float_1d

  !> Unpack a 32-bit real 2d tensor from the database
  function unpack_tensor_float_2d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),              intent(in) :: self !< Fortran communication client
    character(len=*),                  intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),             intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_float_2d

  !> Unpack a 32-bit real 3d tensor from the database
  function unpack_tensor_float_3d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:,:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_float_3d

  !> Unpack a 32-bit real 4d tensor from the database
  function unpack_tensor_float_4d(self, name, data, dims) result(code)
    real(kind=real32), dimension(:,:,:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_float_4d

  !> Unpack a 64-bit real 1d tensor from the database
  function unpack_tensor_double_1d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:), intent(  out) :: data !< Data to be received
    class(dbclient_type),            intent(in) :: self !< Fortran communication client
    character(len=*),                intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),           intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_double_1d

  !> Unpack a 64-bit real 2d tensor from the database
  function unpack_tensor_double_2d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),              intent(in) :: self !< Fortran communication client
    character(len=*),                  intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),             intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_double_2d

  !> Unpack a 64-bit real 3d tensor from the database
  function unpack_tensor_double_3d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:,:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_double_3d

  !> Unpack a 64-bit real 4d tensor from the database
  function unpack_tensor_double_4d(self, name, data, dims) result(code)
    real(kind=real64), dimension(:,:,:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_double_4d

  !> Unpack a 32-bit integer 1d tensor from the database
  function unpack_tensor_int32_1d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:), intent(  out) :: data !< Data to be received
    class(dbclient_type),            intent(in) :: self !< Fortran communication client
    character(len=*),                intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),           intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_int32_1d

  !> Unpack a 32-bit integer 2d tensor from the database
  function unpack_tensor_int32_2d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),              intent(in) :: self !< Fortran communication client
    character(len=*),                  intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),             intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_int32_2d

  !> Unpack a 32-bit integer 3d tensor from the database
  function unpack_tensor_int32_3d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:,:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_int32_3d

  !> Unpack a 32-bit integer 4d tensor from the database
  function unpack_tensor_int32_4d(self, name, data, dims) result(code)
    integer(kind=int32), dimension(:,:,:,:), intent(  out) :: data !< Data to be received
    class(dbclient_type),                intent(in) :: self !< Fortran communication client
    character(len=*),                    intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),               intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function unpack_tensor_int32_4d

  !> Move a tensor to a new name
  function rename_tensor(self, old_name, new_name) result(code)
    class(dbclient_type), intent(in) :: self     !< The initialized Fortran communication client
    character(len=*),   intent(in) :: old_name !< The current name for the tensor
                                               !! excluding null terminating character
    character(len=*),   intent(in) :: new_name !< The new tensor name
    integer        :: code

    code = -1
  end function rename_tensor

  !> Delete a tensor
  function delete_tensor(self, name) result(code)
    class(dbclient_type), intent(in) :: self !< The initialized Fortran communication client
    character(len=*),   intent(in) :: name !< The name associated with the tensor
    integer        :: code

    code = -1
  end function delete_tensor

  !> Copy a tensor to the destination name
  function copy_tensor(self, src_name, dest_name) result(code)
    class(dbclient_type), intent(in) :: self      !< The initialized Fortran communication client
    character(len=*),   intent(in) :: src_name  !< The name associated with the tensor
                                                !! excluding null terminating character
    character(len=*),   intent(in) :: dest_name !< The new tensor name
    integer        :: code

    code = -1
  end function copy_tensor

  !> Retrieve the model from the database
  function get_model(self, name, model) result(code)
    class(dbclient_type),               intent(in  ) :: self  !< An initialized communication client
    character(len=*),                 intent(in  ) :: name  !< The name associated with the model
    character(len=*),                 intent( out) :: model !< The model as a continuous buffer
    integer                        :: code

    code = -1
  end function get_model

  !> Load the machine learning model from a file and set the configuration
  function set_model_from_file(self, name, model_file, backend, device, batch_size, min_batch_size, tag, &
      inputs, outputs) result(code)
    class(dbclient_type),                       intent(in) :: self           !< An initialized communication client
    character(len=*),                         intent(in) :: name           !< The name to use to place the model
    character(len=*),                         intent(in) :: model_file     !< The file storing the model
    character(len=*),                         intent(in) :: backend        !< The name of the backend
                                                                           !! (TF, TFLITE, TORCH, ONNX)
    character(len=*),                         intent(in) :: device         !< The name of the device
                                                                           !! (CPU, GPU, GPU:0, GPU:1...)
    integer,                        optional, intent(in) :: batch_size     !< The batch size for model execution
    integer,                        optional, intent(in) :: min_batch_size !< The minimum batch size for model execution
    character(len=*),               optional, intent(in) :: tag            !< A tag to attach to the model for
                                                                           !! information purposes
    character(len=*), dimension(:), optional, intent(in) :: inputs         !< One or more names of model
                                                                           !! input nodes (TF models)
    character(len=*), dimension(:), optional, intent(in) :: outputs        !< One or more names of model
                                                                           !! output nodes (TF models)
    integer                              :: code

    code = -1
  end function set_model_from_file

  !> Load the machine learning model from a file and set the configuration for use in multi-GPU systems
  function set_model_from_file_multigpu(self, name, model_file, backend, first_gpu, num_gpus, batch_size, &
                                        min_batch_size, tag, inputs, outputs) result(code)
    class(dbclient_type),                       intent(in) :: self           !< An initialized communication client
    character(len=*),                         intent(in) :: name           !< The name to use to place the model
    character(len=*),                         intent(in) :: model_file     !< The file storing the model
    character(len=*),                         intent(in) :: backend        !< The name of the backend
                                                                           !! (TF, TFLITE, TORCH, ONNX)
    integer,                                  intent(in) :: first_gpu      !< The first GPU (zero-based)
                                                                           !! to use with the model
    integer,                                  intent(in) :: num_gpus       !< The number of GPUs to use with the model
    integer,                        optional, intent(in) :: batch_size     !< The batch size for model execution
    integer,                        optional, intent(in) :: min_batch_size !< The minimum batch size for model execution
    character(len=*),               optional, intent(in) :: tag            !< A tag to attach to the model for
                                                                           !! information purposes
    character(len=*), dimension(:), optional, intent(in) :: inputs         !< One or more names of model
                                                                           !! input nodes (TF models)
    character(len=*), dimension(:), optional, intent(in) :: outputs        !< One or more names of model
                                                                           !! output nodes (TF models)
    integer                              :: code

    code = -1
  end function set_model_from_file_multigpu

  !> Establish a model to run
  function set_model(self, name, model, backend, device, batch_size, min_batch_size, tag, &
      inputs, outputs) result(code)
    class(dbclient_type),             intent(in) :: self           !< An initialized communication client
    character(len=*),               intent(in) :: name           !< The name to use to place the model
    character(len=*),               intent(in) :: model          !< The binary representation of the model
    character(len=*),               intent(in) :: backend        !< The name of the backend (TF, TFLITE, TORCH, ONNX)
    character(len=*),               intent(in) :: device         !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
    integer,                        intent(in) :: batch_size     !< The batch size for model execution
    integer,                        intent(in) :: min_batch_size !< The minimum batch size for model execution
    character(len=*),               intent(in) :: tag            !< A tag to attach to the model for
                                                                 !! information purposes
    character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of model input nodes (TF models)
    character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of model output nodes (TF models)
    integer                    :: code

    code = -1
  end function set_model

  !> Set a model from a byte string to run on a system with multiple GPUs
  function set_model_multigpu(self, name, model, backend, first_gpu, num_gpus, batch_size, min_batch_size, tag, &
      inputs, outputs) result(code)
    class(dbclient_type),             intent(in) :: self           !< An initialized communication client
    character(len=*),               intent(in) :: name           !< The name to use to place the model
    character(len=*),               intent(in) :: model          !< The binary representation of the model
    character(len=*),               intent(in) :: backend        !< The name of the backend (TF, TFLITE, TORCH, ONNX)
    integer,                        intent(in) :: first_gpu      !< The first GPU (zero-based) to use with the model
    integer,                        intent(in) :: num_gpus       !< The number of GPUs to use with the model
    integer,                        intent(in) :: batch_size     !< The batch size for model execution
    integer,                        intent(in) :: min_batch_size !< The minimum batch size for model execution
    character(len=*),               intent(in) :: tag            !< A tag to attach to the model for
                                                                 !! information purposes
    character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of model input nodes (TF models)
    character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of model output nodes (TF models)
    integer                    :: code

    code = -1
  end function set_model_multigpu

  !> Run a model in the database using the specified input and output tensors
  function run_model(self, name, inputs, outputs) result(code)
    class(dbclient_type),             intent(in) :: self    !< An initialized communication client
    character(len=*),               intent(in) :: name    !< The name to use to place the model
    character(len=*), dimension(:), intent(in) :: inputs  !< One or more names of model input nodes (TF models)
    character(len=*), dimension(:), intent(in) :: outputs !< One or more names of model output nodes (TF models)
    integer                    :: code

    code = -1
  end function run_model

  !> Run a model in the database using the specified input and output tensors in a multi-GPU system
  function run_model_multigpu(self, name, inputs, outputs, offset, first_gpu, num_gpus) result(code)
    class(dbclient_type),             intent(in) :: self    !< An initialized communication client
    character(len=*),               intent(in) :: name    !< The name to use to place the model
    character(len=*), dimension(:), intent(in) :: inputs  !< One or more names of model input nodes (TF models)
    character(len=*), dimension(:), intent(in) :: outputs !< One or more names of model output nodes (TF models)
    integer,                        intent(in) :: offset  !< Index of the current image, such as a processor ID
                                                          !! or MPI rank
    integer,                        intent(in) :: first_gpu !< The first GPU (zero-based) to use with the model
    integer,                        intent(in) :: num_gpus  !< The number of GPUs to use with the model
    integer                    :: code

    code = -1
  end function run_model_multigpu

  !> Remove a model from the database
  function delete_model(self, name) result(code)
    class(dbclient_type),             intent(in) :: self    !< An initialized communication client
    character(len=*),               intent(in) :: name    !< The name to use to remove the model
    integer                    :: code

    code = -1
  end function delete_model

  !> Remove a model from the database
  function delete_model_multigpu(self, name, first_gpu, num_gpus) result(code)
    class(dbclient_type),             intent(in) :: self    !< An initialized communication client
    character(len=*),               intent(in) :: name    !< The name to use to remove the model
    integer,                        intent(in) :: first_gpu !< The first GPU (zero-based) to use with the model
    integer,                        intent(in) :: num_gpus !< The number of GPUs to use with the model
    integer                    :: code

    code = -1
  end function delete_model_multigpu

  !> Retrieve the script from the database
  function get_script(self, name, script) result(code)
    class(dbclient_type), intent(in  ) :: self   !< An initialized communication client
    character(len=*),   intent(in  ) :: name   !< The name to use to place the script
    character(len=*),   intent( out) :: script !< The script as a continuous buffer
    integer          :: code

    code = -1
  end function get_script

  !> Set a script (from file) in the database for future execution
  function set_script_from_file(self, name, device, script_file) result(code)
    class(dbclient_type), intent(in) :: self        !< An initialized communication client
    character(len=*),   intent(in) :: name        !< The name to use to place the script
    character(len=*),   intent(in) :: device      !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
    character(len=*),   intent(in) :: script_file !< The file storing the script
    integer        :: code

    code = -1
  end function set_script_from_file

  !> Set a script (from file) in the database for future execution in a multi-GPU system
  function set_script_from_file_multigpu(self, name, script_file, first_gpu, num_gpus) result(code)
    class(dbclient_type), intent(in) :: self        !< An initialized communication client
    character(len=*),   intent(in) :: name        !< The name to use to place the script
    character(len=*),   intent(in) :: script_file !< The file storing the script
    integer,            intent(in) :: first_gpu   !< The first GPU (zero-based) to use with the model
    integer,            intent(in) :: num_gpus    !< The number of GPUs to use with the model
    integer        :: code

    code = -1
  end function set_script_from_file_multigpu

  !> Set a script (from buffer) in the database for future execution
  function set_script(self, name, device, script) result(code)
    class(dbclient_type), intent(in) :: self   !< An initialized communication client
    character(len=*),   intent(in) :: name   !< The name to use to place the script
    character(len=*),   intent(in) :: device !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
    character(len=*),   intent(in) :: script !< The file storing the script
    integer        :: code

    code = -1
  end function set_script

  !> Set a script (from buffer) in the database for future execution in a multi-GPU system
  function set_script_multigpu(self, name, script, first_gpu, num_gpus) result(code)
    class(dbclient_type), intent(in) :: self   !< An initialized communication client
    character(len=*),   intent(in) :: name   !< The name to use to place the script
    character(len=*),   intent(in) :: script !< The file storing the script
    integer,            intent(in) :: first_gpu !< The first GPU (zero-based) to use with the model
    integer,            intent(in) :: num_gpus  !< The number of GPUs to use with the model
    integer        :: code

    code = -1
  end function set_script_multigpu

  function run_script(self, name, func, inputs, outputs) result(code)
    class(dbclient_type),             intent(in) :: self           !< An initialized communication client
    character(len=*),               intent(in) :: name           !< The name to use to place the script
    character(len=*),               intent(in) :: func           !< The name of the function in the script to call
    character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of script
                                                                 !! input nodes (TF scripts)
    character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of script
                                                                 !! output nodes (TF scripts)
    integer                    :: code

    code = -1
  end function run_script

  function run_script_multigpu(self, name, func, inputs, outputs, offset, first_gpu, num_gpus) result(code)
    class(dbclient_type),             intent(in) :: self           !< An initialized communication client
    character(len=*),               intent(in) :: name           !< The name to use to place the script
    character(len=*),               intent(in) :: func           !< The name of the function in the script to call
    character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of script
                                                                 !! input nodes (TF scripts)
    character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of script
                                                                 !! output nodes (TF scripts)
    integer,                        intent(in) :: offset  !< Index of the current image, such as a processor ID
                                                          !! or MPI rank
    integer,                        intent(in) :: first_gpu !< The first GPU (zero-based) to use with the model
    integer,                        intent(in) :: num_gpus  !< The number of GPUs to use with the model
    integer                    :: code

    code = -1
  end function run_script_multigpu

  !> Remove a script from the database
  function delete_script(self, name) result(code)
    class(dbclient_type),             intent(in) :: self    !< An initialized communication client
    character(len=*),               intent(in) :: name    !< The name to use to delete the script
    integer                    :: code

    code = -1
  end function delete_script

  !> Remove a script_multigpu from the database
  function delete_script_multigpu(self, name, first_gpu, num_gpus) result(code)
    class(dbclient_type),             intent(in) :: self    !< An initialized communication client
    character(len=*),               intent(in) :: name    !< The name to use to delete the script_multigpu
    integer,                        intent(in) :: first_gpu !< The first GPU (zero-based) to use with the model
    integer,                        intent(in) :: num_gpus !< The number of GPUs to use with the model
    integer                    :: code

    code = -1
  end function delete_script_multigpu

  !> Store a dataset in the database
  function put_dataset(self, dataset) result(code)
    class(dbclient_type), intent(in) :: self    !< An initialized communication client
    type(dataset_type), intent(in) :: dataset !< Dataset to store in the dataset
    integer        :: code

    code = -1
  end function put_dataset

  !> Retrieve a dataset from the database
  function get_dataset(self, name, dataset) result(code)
    class(dbclient_type), intent(in )  :: self    !< An initialized communication client
    character(len=*),   intent(in )  :: name    !< Name of the dataset to get
    type(dataset_type), intent( out) :: dataset !< receives the dataset
    integer          :: code

    code = -1
  end function get_dataset

  !> Rename a dataset stored in the database
  function rename_dataset(self, name, new_name) result(code)
    class(dbclient_type), intent(in) :: self     !< An initialized communication client
    character(len=*),   intent(in) :: name     !< Original name of the dataset
    character(len=*),   intent(in) :: new_name !< New name of the dataset
    integer        :: code

    code = -1
  end function rename_dataset

  !> Copy a dataset within the database to a new name
  function copy_dataset(self, name, new_name) result(code)
    class(dbclient_type), intent(in) :: self     !< An initialized communication client
    character(len=*),   intent(in) :: name     !< Source name of the dataset
    character(len=*),   intent(in) :: new_name !< Name of the new dataset
    integer        :: code

    code = -1
  end function copy_dataset

  !> Delete a dataset stored within a database
  function delete_dataset(self, name) result(code)
    class(dbclient_type), intent(in) :: self !< An initialized communication client
    character(len=*),   intent(in) :: name !< Name of the dataset to delete
    integer        :: code

    code = -1
  end function delete_dataset

  !> Appends a dataset to the aggregation list When appending a dataset to an aggregation list, the list will
  !! automatically be created if it does not exist (i.e. this is the first entry in the list). Aggregation
  !! lists work by referencing the dataset by storing its key, so appending a dataset to an aggregation list
  !! does not create a copy of the dataset.  Also, for this reason, the dataset must have been previously
  !! placed into the database with a separate call to put_dataset().
  function append_to_list(self, list_name, dataset) result(code)
    class(dbclient_type), intent(in) :: self       !< An initialized communication client
    character(len=*),   intent(in) :: list_name  !< Name of the dataset to get
    type(dataset_type), intent(in) :: dataset    !< Dataset to append to the list
    integer        :: code

    code = -1
  end function append_to_list

  !> Delete an aggregation list
  function delete_list(self, list_name) result(code)
    class(dbclient_type),   intent(in) :: self       !< An initialized communication client
    character(len=*),     intent(in) :: list_name  !< Name of the aggregated dataset list to delete
    integer          :: code

    code = -1
  end function delete_list

  !> Copy an aggregation list
  function copy_list(self, src_name, dest_name) result(code)
    class(dbclient_type),   intent(in) :: self      !< An initialized communication client
    character(len=*),     intent(in) :: src_name  !< Name of the dataset to copy
    character(len=*),     intent(in) :: dest_name !< The new list name
    integer          :: code

    code = -1
  end function copy_list

  !> Rename an aggregation list
  function rename_list(self, src_name, dest_name) result(code)
    class(dbclient_type),   intent(in) :: self      !< An initialized communication client
    character(len=*),     intent(in) :: src_name  !< Name of the dataset to rename
    character(len=*),     intent(in) :: dest_name !< The new list name
    integer          :: code

    code = -1
  end function rename_list

  end module database_client_interface

