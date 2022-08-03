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
    !> Puts a tensor into the database (overloaded)
    generic :: put_tensor => put_tensor_i8, put_tensor_i16, put_tensor_i32, put_tensor_i64, &
                             put_tensor_float, put_tensor_double
    !> Retrieve the tensor in the database into already allocated memory (overloaded)
    generic :: unpack_tensor => unpack_tensor_i8, unpack_tensor_i16, unpack_tensor_i32, unpack_tensor_i64, &
                                unpack_tensor_float, unpack_tensor_double

    !> Decode a response code from an API function
    procedure :: SR_error_parser
    !> Initializes a new instance of the communication client
    procedure :: initialize => initialize_client
    !> Check if a communication client has been initialized
    procedure :: isinitialized
    !> Destructs a new instance of the communication client
    procedure :: destructor
    !> Check the database for the existence of a specific model
    procedure :: model_exists
    !> Check the database for the existence of a specific tensor
    procedure :: tensor_exists
    !> Check the database for the existence of a specific key
    procedure :: key_exists
    !> Check the database for the existence of a specific dataset
    procedure :: dataset_exists
    !> Poll the database and return if the model exists
    procedure :: poll_model
    !> Poll the database and return if the tensor exists
    procedure :: poll_tensor
    !> Poll the database and return if the datasaet exists
    procedure :: poll_dataset
    !> Poll the database and return if the key exists
    procedure :: poll_key
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

    !> If true, preprend the ensemble id for tensor-related keys
    procedure :: use_tensor_ensemble_prefix
    !> If true, preprend the ensemble id for model-related keys
    procedure :: use_model_ensemble_prefix
    !> If true, preprend the ensemble id for dataset list-related keys
    procedure :: use_list_ensemble_prefix
    !> Specify a specific source of data (e.g. another ensemble member)
    procedure :: set_data_source

    !> Append a dataset to a list for aggregation
    procedure :: append_to_list
    !> Delete an aggregation list
    procedure :: delete_list
    !> Copy an aggregation list
    procedure :: copy_list
    !> Rename an existing aggregation list
    procedure :: rename_list
    !> Retrieve the number of datasets in the list
    procedure :: get_list_length
    !> Repeatedly check the length of the list until it is a given size
    procedure :: poll_list_length
    !> Repeatedly check the length of the list until it greater than or equal to the given size
    procedure :: poll_list_length_gte
    !> Repeatedly check the length of the list until it less than or equal to the given size
    procedure :: poll_list_length_lte
    !> Retrieve vector of datasetes from the list
    procedure :: get_datasets_from_list

    ! Private procedures
    procedure, private :: put_tensor_i8     !< Put 8-bit integer tensor into database
    procedure, private :: put_tensor_i16    !< Put 16-bit integer tensor into database
    procedure, private :: put_tensor_i32    !< Put 32-bit integer tensor into database
    procedure, private :: put_tensor_i64    !< Put 64-bit tensor into database
    procedure, private :: put_tensor_float  !< Put 32-bit real tensor into database
    procedure, private :: put_tensor_double !< Put 64-bit real tensor into database
    procedure, private :: unpack_tensor_i8     !< Unpack a 8-bit integer tensor into memory
    procedure, private :: unpack_tensor_i16    !< Unpack a 16-bit integer tensor into memory
    procedure, private :: unpack_tensor_i32    !< Unpack a 32-bit integer tensor into memory
    procedure, private :: unpack_tensor_i64    !< Unpack a 64-bit integer tensor into memory
    procedure, private :: unpack_tensor_float  !< Unpack a 32-bit real tensor into memory
    procedure, private :: unpack_tensor_double !< Unpack a 64-bit real tensor into memory

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

  !> Check if the specified key exists in the database
  function key_exists(self, key, exists)
    class(dbclient_type),   intent(in)  :: self   !< The client
    character(len=*),     intent(in)  :: key    !< The key to check
    logical, intent(out) :: exists !< Receives whether the key exists
    integer           :: key_exists

    key_exists = -1
  end function key_exists

  !> Check if the specified model exists in the database
  function model_exists(self, model_name, exists) result(code)
    class(dbclient_type),   intent(in)  :: self       !< The client
    character(len=*),     intent(in)  :: model_name !< The model to check
    logical, intent(out) :: exists     !< Receives whether the model exists
    integer           :: code

    code = -1
  end function model_exists

  !> Check if the specified tensor exists in the database
  function tensor_exists(self, tensor_name, exists) result(code)
    class(dbclient_type),   intent(in)  :: self        !< The client
    character(len=*),     intent(in)  :: tensor_name !< The tensor to check
    logical, intent(out) :: exists      !< Receives whether the model exists
    integer           :: code

    code = -1
  end function tensor_exists

  !> Check if the specified dataset exists in the database
  function dataset_exists(this, dataset_name, exists) result(code)
    class(dbclient_type),   intent(in)  :: this          !< The client
    character(len=*),     intent(in)  :: dataset_name  !< The dataset to check
    logical, intent(out) :: exists        !< Receives whether the model exists
    integer           :: code

    code = -1
  end function dataset_exists

  !> Repeatedly poll the database until the tensor exists or the number of tries is exceeded
  function poll_tensor(self, tensor_name, poll_frequency_ms, num_tries, exists) result(code)
    class(dbclient_type),   intent(in)  :: self              !< The client
    character(len=*),     intent(in)  :: tensor_name       !< name in the database to poll
    integer,              intent(in)  :: poll_frequency_ms !< Frequency at which to poll the database (ms)
    integer,              intent(in)  :: num_tries         !< Number of times to poll the database before failing
    logical, intent(out) :: exists            !< Receives whether the tensor exists
    integer           :: code

    code = -1
  end function poll_tensor

  !> Repeatedly poll the database until the dataset exists or the number of tries is exceeded
  function poll_dataset(self, dataset_name, poll_frequency_ms, num_tries, exists)
    integer           :: poll_dataset
    class(dbclient_type),   intent(in)  :: self              !< The client
    character(len=*),     intent(in)  :: dataset_name      !< Name in the database to poll
    integer,              intent(in)  :: poll_frequency_ms !< Frequency at which to poll the database (ms)
    integer,              intent(in)  :: num_tries         !< Number of times to poll the database before failing
    logical, intent(out) :: exists            !< Receives whether the tensor exists

    poll_dataset = -1
  end function poll_dataset

  !> Repeatedly poll the database until the model exists or the number of tries is exceeded
  function poll_model(self, model_name, poll_frequency_ms, num_tries, exists) result(code)
    class(dbclient_type),   intent(in)  :: self              !< The client
    character(len=*),     intent(in)  :: model_name        !< Name in the database to poll
    integer,              intent(in)  :: poll_frequency_ms !< Frequency at which to poll the database (ms)
    integer,              intent(in)  :: num_tries         !< Number of times to poll the database before failing
    logical, intent(out) :: exists            !< Receives whether the model exists
    integer           :: code

    code = -1
  end function poll_model

  !> Repeatedly poll the database until the key exists or the number of tries is exceeded
  function poll_key(self, key, poll_frequency_ms, num_tries, exists) result(code)
    class(dbclient_type),   intent(in)  :: self               !< The client
    character(len=*),     intent(in)  :: key                !< Key in the database to poll
    integer,              intent(in)  :: poll_frequency_ms  !< Frequency at which to poll the database (ms)
    integer,              intent(in)  :: num_tries          !< Number of times to poll the database before failing
    logical, intent(out) :: exists             !< Receives whether the key exists
    integer           :: code

    code = -1
  end function poll_key

  !> Put a tensor whose Fortran type is the equivalent 'int8' C-type
  function put_tensor_i8(self, name, data, dims) result(code)
    integer(kind=int8), dimension(..), target, intent(in) :: data !< Data to be sent
    class(dbclient_type),                    intent(in) :: self !< Fortran communication client
    character(len=*),                      intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_i8

  !> Put a tensor whose Fortran type is the equivalent 'int16' C-type
  function put_tensor_i16(self, name, data, dims) result(code)
    integer(kind=int16), dimension(..), target, intent(in) :: data !< Data to be sent
    class(dbclient_type),                    intent(in) :: self !< Fortran communication client
    character(len=*),                      intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_i16

  !> Put a tensor whose Fortran type is the equivalent 'int32' C-type
  function put_tensor_i32(self, name, data, dims) result(code)
    integer(kind=int32), dimension(..), target, intent(in) :: data !< Data to be sent
    class(dbclient_type),                    intent(in) :: self !< Fortran communication client
    character(len=*),                      intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_i32

  !> Put a tensor whose Fortran type is the equivalent 'int64' C-type
  function put_tensor_i64(self, name, data, dims) result(code)
    integer(kind=int64), dimension(..), target, intent(in) :: data !< Data to be sent
    class(dbclient_type),                    intent(in) :: self !< Fortran communication client
    character(len=*),                      intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_i64

  !> Put a tensor whose Fortran type is the equivalent 'float' C-type
  function put_tensor_float(self, name, data, dims) result(code)
    real(kind=real32), dimension(..), target, intent(in) :: data !< Data to be sent
    class(dbclient_type),                    intent(in) :: self !< Fortran communication client
    character(len=*),                      intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_float

  !> Put a tensor whose Fortran type is the equivalent 'double' C-type
  function put_tensor_double(self, name, data, dims) result(code)
    real(kind=real64), dimension(..), target, intent(in) :: data !< Data to be sent
    class(dbclient_type),                    intent(in) :: self !< Fortran communication client
    character(len=*),                      intent(in) :: name !< The unique name used to store in the database
    integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
    integer                           :: code

    code = -1
  end function put_tensor_double

  !> Put a tensor whose Fortran type is the equivalent 'int8' C-type
  function unpack_tensor_i8(self, name, result, dims) result(code)
    integer(kind=int8), dimension(..), target, intent(out) :: result !< Data to be sent
    class(dbclient_type),                   intent(in) :: self  !< Pointer to the initialized client
    character(len=*),                     intent(in) :: name  !< The name to use to place the tensor
    integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor
    integer                          :: code

    code = -1
  end function unpack_tensor_i8

  !> Put a tensor whose Fortran type is the equivalent 'int16' C-type
  function unpack_tensor_i16(self, name, result, dims) result(code)
    integer(kind=int16), dimension(..), target, intent(out) :: result !< Data to be sent
    class(dbclient_type),                   intent(in) :: self  !< Pointer to the initialized client
    character(len=*),                     intent(in) :: name  !< The name to use to place the tensor
    integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor
    integer                          :: code

    code = -1
  end function unpack_tensor_i16

  !> Put a tensor whose Fortran type is the equivalent 'int32' C-type
  function unpack_tensor_i32(self, name, result, dims) result(code)
    integer(kind=int32), dimension(..), target, intent(out) :: result !< Data to be sent
    class(dbclient_type),                   intent(in) :: self  !< Pointer to the initialized client
    character(len=*),                     intent(in) :: name  !< The name to use to place the tensor
    integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor
    integer                          :: code

    code = -1
  end function unpack_tensor_i32

  !> Put a tensor whose Fortran type is the equivalent 'int64' C-type
  function unpack_tensor_i64(self, name, result, dims) result(code)
    integer(kind=int64), dimension(..), target, intent(out) :: result !< Data to be sent
    class(dbclient_type),                   intent(in) :: self  !< Pointer to the initialized client
    character(len=*),                     intent(in) :: name  !< The name to use to place the tensor
    integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor
    integer                          :: code

    code = -1
  end function unpack_tensor_i64

  !> Put a tensor whose Fortran type is the equivalent 'float' C-type
  function unpack_tensor_float(self, name, result, dims) result(code)
    real(kind=real32), dimension(..), target, intent(out) :: result !< Data to be sent
    class(dbclient_type),                   intent(in) :: self  !< Pointer to the initialized client
    character(len=*),                     intent(in) :: name  !< The name to use to place the tensor
    integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor
    integer                          :: code

    code = -1
  end function unpack_tensor_float

  !> Put a tensor whose Fortran type is the equivalent 'double' C-type
  function unpack_tensor_double(self, name, result, dims) result(code)
    real(kind=real64), dimension(..), target, intent(out) :: result !< Data to be sent
    class(dbclient_type),                   intent(in) :: self  !< Pointer to the initialized client
    character(len=*),                     intent(in) :: name  !< The name to use to place the tensor
    integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor
    integer                          :: code

    code = -1
  end function unpack_tensor_double

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

  !> Set the data source (i.e. name prefix for get functions)
  function set_data_source(self, source_id) result(code)
    class(dbclient_type), intent(in) :: self      !< An initialized communication client
    character(len=*),   intent(in) :: source_id !< The name prefix
    integer        :: code

    code = -1
  end function set_data_source

  !> Set whether names of model and script entities should be prefixed (e.g. in an ensemble) to form database names.
  !! Prefixes will only be used if they were previously set through the environment variables SSKEYOUT and SSKEYIN.
  !! Keys of entities created before client function is called will not be affected. By default, the client does not
  !! prefix model and script names.
  function use_model_ensemble_prefix(self, use_prefix) result(code)
    class(dbclient_type),   intent(in) :: self       !< An initialized communication client
    logical,              intent(in) :: use_prefix !< The prefix setting
    integer          :: code

    code = -1
  end function use_model_ensemble_prefix


  !> Set whether names of tensor and dataset entities should be prefixed (e.g. in an ensemble) to form database keys.
  !! Prefixes will only be used if they were previously set through the environment variables SSKEYOUT and SSKEYIN.
  !! Keys of entities created before client function is called will not be affected. By default, the client prefixes
  !! tensor and dataset keys with the first prefix specified with the SSKEYIN and SSKEYOUT environment variables.
  function use_tensor_ensemble_prefix(self, use_prefix) result(code)
    class(dbclient_type),   intent(in) :: self       !< An initialized communication client
    logical,              intent(in) :: use_prefix !< The prefix setting
    integer          :: code

    code = -1
  end function use_tensor_ensemble_prefix

  !> Control whether aggregation lists are prefixed
  function use_list_ensemble_prefix(self, use_prefix) result(code)
    class(dbclient_type),   intent(in) :: self       !< An initialized communication client
    logical,              intent(in) :: use_prefix !< The prefix setting
    integer          :: code

    code = -1
  end function use_list_ensemble_prefix

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

  !> Get the length of the aggregation list
  function get_list_length(self, list_name, result_length) result(code)
    class(dbclient_type),   intent(in   ) :: self           !< An initialized communication client
    character(len=*),     intent(in   ) :: list_name      !< Name of the dataset to get
    integer,              intent(  out) :: result_length  !< The length of the list
    integer             :: code

    code = -1
  end function get_list_length

  !> Get the length of the aggregation list
  function poll_list_length(self, list_name, list_length, poll_frequency_ms, num_tries, poll_result) result(code)
    class(dbclient_type),   intent(in   ) :: self               !< An initialized communication client
    character(len=*),     intent(in   ) :: list_name          !< Name of the dataset to get
    integer,              intent(in   ) :: list_length        !< The desired length of the list
    integer,              intent(in   ) :: poll_frequency_ms  !< Frequency at which to poll the database (ms)
    integer,              intent(in   ) :: num_tries          !< Number of times to poll the database before failing
    logical, intent(  out) :: poll_result        !< True if the list is the requested length,
                                                              !! False if not after num_tries.
    integer             :: code

    code = -1
  end function poll_list_length

  !> Get the length of the aggregation list
  function poll_list_length_gte(self, list_name, list_length, poll_frequency_ms, num_tries, poll_result) result(code)
    class(dbclient_type),   intent(in   ) :: self               !< An initialized communication client
    character(len=*),     intent(in   ) :: list_name          !< Name of the dataset to get
    integer,              intent(in   ) :: list_length        !< The desired length of the list
    integer,              intent(in   ) :: poll_frequency_ms  !< Frequency at which to poll the database (ms)
    integer,              intent(in   ) :: num_tries          !< Number of times to poll the database before failing
    logical, intent(  out) :: poll_result        !< True if the list is the requested length,
                                                              !! False if not after num_tries.
    integer          :: code

    code = -1
  end function poll_list_length_gte

  !> Get the length of the aggregation list
  function poll_list_length_lte(self, list_name, list_length, poll_frequency_ms, num_tries, poll_result) result(code)
    class(dbclient_type),   intent(in) :: self                !< An initialized communication client
    character(len=*),     intent(in) :: list_name           !< Name of the dataset to get
    integer,              intent(in)  :: list_length        !< The desired length of the list
    integer,              intent(in)  :: poll_frequency_ms  !< Frequency at which to poll the database (ms)
    integer,              intent(in)  :: num_tries          !< Number of times to poll the database before failing
    logical, intent(  out) :: poll_result        !< True if the list is the requested length,
                                                              !! False if not after num_tries.

    integer          :: code

    code = -1
  end function poll_list_length_lte

  !> Get datasets from an aggregation list. Note that this will deallocate an existing list.
  !! NOTE: This potentially be less performant than get_datasets_from_list_range due to an
  !! extra query to the database to get the list length. This is for now necessary because
  !! difficulties in allocating memory for Fortran alloctables from within C.
  function get_datasets_from_list(self, list_name, datasets, num_datasets) result(code)
    class(dbclient_type),   intent(in) :: self       !< An initialized communication client
    character(len=*),     intent(in) :: list_name  !< Name of the dataset to get
    type(dataset_type), dimension(:), allocatable, intent(  out) :: datasets !< The array of datasets included
    integer          :: code
                                                                             !! in the list
    integer,              intent(out) :: num_datasets !< The numbr of datasets returned

    code = -1
  end function get_datasets_from_list

  !> Get datasets from an aggregation list over a given range by index. Note that this will deallocate an existing list
  function get_datasets_from_list_range(self, list_name, start_index, end_index, datasets) result(code)
    class(dbclient_type),   intent(in) :: self        !< An initialized communication client
    character(len=*),     intent(in) :: list_name   !< Name of the dataset to get
    integer,              intent(in) :: start_index !< The starting index of the range (inclusive,
                                                    !! starting at zero).  Negative values are
                                                    !! supported.  A negative value indicates offsets
                                                    !! starting at the end of the list. For example, -1 is
                                                    !! the last element of the list.
    integer,              intent(in) :: end_index   !< The ending index of the range (inclusive,
                                                    !! starting at zero).  Negative values are
                                                    !! supported.  A negative value indicates offsets
                                                    !! starting at the end of the list. For example, -1 is
                                                    !! the last element of the list.

    type(dataset_type), dimension(:), allocatable, intent(  out) :: datasets !< The array of datasets included
    integer          :: code
                                                                             !! in the list

    code = -1
  end function get_datasets_from_list_range

  end module database_client_interface

