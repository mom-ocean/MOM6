! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements the general purpose Artificial Neural Network (ANN).
module MOM_ANN

! This file is part of MOM6. See LICENSE.md for the license

use MOM_io, only : MOM_read_data, field_exists
use MOM_error_handler, only : MOM_error, FATAL, MOM_mesg
use numerical_testing_type, only : testing

implicit none ; private

!#include <MOM_memory.h>

public ANN_init, ANN_allocate, ANN_apply, ANN_end, ANN_unit_tests
public ANN_apply_vector_orig, ANN_apply_vector_oi, ANN_apply_array_sio
public set_layer, set_input_normalization, set_output_normalization
public ANN_random, randomize_layer

!> Applies ANN to x, returning results in y
interface ANN_apply
  module procedure ANN_apply_vector_oi
  module procedure ANN_apply_array_sio
end interface ANN_apply

!> Type for a single Linear layer of ANN,
!! i.e. stores the matrix A and bias b
!! for matrix-vector multiplication
!! y = A*x + b.
type, private :: layer_type ; private
  integer :: output_width        !< Number of rows in matrix A
  integer :: input_width         !< Number of columns in matrix A
  logical :: activation = .True. !< If true, apply the default activation function

  real, allocatable :: A(:,:) !< Matrix in column-major order
                              !! of size A(output_width, input_width) [nondim]
  real, allocatable :: b(:)   !< bias vector of size output_width [nondim]
end type layer_type

!> Control structure/type for ANN
type, public :: ANN_CS ; private
  ! Parameters
  integer :: num_layers          !< Number of layers in the ANN, including the input and output.
                                 !! For example, for ANN with one hidden layer, num_layers = 3.
  integer, allocatable &
          :: layer_sizes(:)      !< Array of length num_layers, storing the number of neurons in
                                 !! each layer.

  type(layer_type), allocatable &
          :: layers(:)           !< Array of length num_layers-1, where each element is the Linear
                                 !! transformation between layers defined by Matrix A and vias b.

  real, allocatable :: &
    input_means(:), &  !< Array of length layer_sizes(1) containing the mean of each input feature
                       !! prior to normalization by input_norms [arbitrary].
    input_norms(:), &  !< Array of length layer_sizes(1) containing the *inverse* of the standard
                       !! deviation for each input feature used to normalize (multiply) before
                       !! feeding into the ANN [arbitrary]
    output_means(:), & !< Array of length layer_sizes(num_layers) containing the mean of each
                       !! output prior to normalization by output_norms [arbitrary].
    output_norms(:)    !< Array of length layer_sizes(num_layers) containing the standard deviation
                       !! each output of the ANN will be multiplied [arbitrary]

  integer, public :: parameters = 0 !< Count of number of parameters
end type ANN_CS

contains

!> Initialization of ANN. Allocates memory and reads ANN parameters from NetCDF file.
!! The NetCDF file must contain:
!! Integer num_layers.
!! Integer arrays: layer_sizes, input_norms, output_norms
!! Matrices and biases for Linear layers can be Real(4) or Real(8) and
!! are named as: A0, b0 for the first layer; A1, b1 for the second layer and so on.
subroutine ANN_init(CS, NNfile)
  type(ANN_CS), intent(inout)  :: CS     !< ANN control structure.
  character(*), intent(in)     :: NNfile !< The name of NetCDF file having neural network parameters
  ! Local variables
  integer :: i
  integer :: num_layers ! Number of layers, including input and output layers
  integer, allocatable :: layer_sizes(:) ! Number of neurons in each layer
  character(len=1) :: layer_num_str
  character(len=3) :: fieldname

  call MOM_mesg('ANN: init from ' // trim(NNfile), 2)

  ! Read the number of layers
  call MOM_read_data(NNfile, "num_layers", num_layers)

  ! Read size of layers
  allocate( layer_sizes(num_layers) )
  call MOM_read_data(NNfile, "layer_sizes", layer_sizes)

  ! Allocates the memory for storing normalization, weights and biases
  call ANN_allocate(CS, num_layers, layer_sizes)
  deallocate( layer_sizes )

  ! Read normalization factors
  if (field_exists(NNfile, 'input_means')) &
    call MOM_read_data(NNfile, 'input_means', CS%input_means)
  if (field_exists(NNfile, 'input_norms')) then
    call MOM_read_data(NNfile, 'input_norms', CS%input_norms)
    ! We calculate the reciprocal here to avoid repeated divisions later
    CS%input_norms(:) = 1.  / CS%input_norms(:)
  endif
  if (field_exists(NNfile, 'output_means')) &
    call MOM_read_data(NNfile, 'output_means', CS%output_means)
  if (field_exists(NNfile, 'output_norms')) &
    call MOM_read_data(NNfile, 'output_norms', CS%output_norms)

  ! Allocate and read matrix A and bias b for each layer
  do i = 1,CS%num_layers-1
    CS%layers(i)%input_width = CS%layer_sizes(i)
    CS%layers(i)%output_width = CS%layer_sizes(i+1)

    ! Reading matrix A
    write(layer_num_str, '(I0)') i-1
    fieldname = trim('A') // trim(layer_num_str)
    call MOM_read_data(NNfile, fieldname, CS%layers(i)%A, &
                        (/1,1,1,1/),(/CS%layers(i)%output_width,CS%layers(i)%input_width,1,1/))

    ! Reading bias b
    fieldname = trim('b') // trim(layer_num_str)
    call MOM_read_data(NNfile, fieldname, CS%layers(i)%b)
  enddo

  ! No activation function for the last layer
  CS%layers(CS%num_layers-1)%activation = .False.

  if (field_exists(NNfile, 'x_test') .and. field_exists(NNfile, 'y_test') ) &
  call ANN_test(CS, NNfile)

  call MOM_mesg('ANN: have been read from ' // trim(NNfile), 2)

end subroutine ANN_init

!> Allocate an ANN
!!
!! This creates the memory for storing weights and intermediate work arrays, but does not set
!! the values of weights or biases (not even initializing with zeros).
subroutine ANN_allocate(CS, num_layers, layer_sizes)
  type(ANN_CS), intent(inout) :: CS !< ANN control structure
  integer,      intent(in)    :: num_layers !< The number of layers, including the input and output layer
  integer,      intent(in)    :: layer_sizes(num_layers) !< The number of neurons in each layer
  ! Local variables
  integer :: l ! Layer number

  ! Assert that there is always an input and output layer
  if (num_layers < 2) call MOM_error(FATAL, "The number of layers in an ANN must be >=2")

  CS%num_layers = num_layers

  ! Layers
  allocate( CS%layer_sizes(CS%num_layers) )
  CS%layer_sizes(:) = layer_sizes(:)

  ! Input and output normalization values
  allocate( CS%input_means(CS%layer_sizes(1)), source=0. ) ! Assume zero mean by default
  allocate( CS%input_norms(CS%layer_sizes(1)), source=1. ) ! Assume unit variance by default
  allocate( CS%output_means(CS%layer_sizes(CS%num_layers)), source=0. ) ! Assume zero mean by default
  allocate( CS%output_norms(CS%layer_sizes(CS%num_layers)), source=1. ) ! Assume unit variance by default

  ! Allocate the Linear transformations between layers
  allocate(CS%layers(CS%num_layers-1))
  CS%parameters = 2 * CS%layer_sizes(1) ! For input normalization

  ! Allocate matrix A and bias b for each layer
  do l = 1, CS%num_layers-1
    CS%layers(l)%input_width = CS%layer_sizes(l)
    CS%layers(l)%output_width = CS%layer_sizes(l+1)

    allocate( CS%layers(l)%A(CS%layers(l)%output_width, CS%layers(l)%input_width) )
    allocate( CS%layers(l)%b(CS%layers(l)%output_width) )

    CS%parameters = CS%parameters &
       + CS%layer_sizes(l) * CS%layer_sizes(l+1) & ! For weights
       + CS%layer_sizes(l+1) ! For bias
  enddo
  CS%parameters = CS%parameters &
     + 2 * CS%layer_sizes(CS%num_layers) ! For output normalization

end subroutine ANN_allocate

!> Test ANN by comparing the prediction with the test data.
subroutine ANN_test(CS, NNfile)
  type(ANN_CS), intent(inout) :: CS     !< ANN control structure.
  character(*), intent(in)    :: NNfile !< The name of NetCDF file having neural network parameters
  ! Local variables
  real, dimension(:), allocatable :: x_test, y_test, y_pred ! [arbitrary]
  real :: relative_error ! [arbitrary]
  character(len=200) :: relative_error_str

  ! Allocate data
  allocate(x_test(CS%layer_sizes(1)))
  allocate(y_test(CS%layer_sizes(CS%num_layers)))
  allocate(y_pred(CS%layer_sizes(CS%num_layers)))

  ! Read test vectors
  call MOM_read_data(NNfile, 'x_test', x_test)
  call MOM_read_data(NNfile, 'y_test', y_test)

  ! Compute prediction
  call ANN_apply_vector_oi(x_test, y_pred, CS)

  relative_error = maxval(abs(y_pred(:) - y_test(:))) / maxval(abs(y_test(:)))

  if (relative_error > 1e-5) then
    write(relative_error_str, '(ES12.4)') relative_error
    call MOM_error(FATAL, 'Relative error in ANN prediction is too large: ' // trim(relative_error_str))
  endif

  deallocate(x_test)
  deallocate(y_test)
  deallocate(y_pred)
end subroutine ANN_test

!> Deallocates memory of ANN
subroutine ANN_end(CS)
  type(ANN_CS), intent(inout) :: CS !< ANN control structure.
  ! Local variables
  integer :: i

  deallocate(CS%layer_sizes)
  deallocate(CS%input_means)
  deallocate(CS%input_norms)
  deallocate(CS%output_means)
  deallocate(CS%output_norms)

  do i = 1, CS%num_layers-1
    deallocate(CS%layers(i)%A)
    deallocate(CS%layers(i)%b)
  enddo
  deallocate(CS%layers)

end subroutine ANN_end

!> The default activation function
pure elemental function activation_fn(x) result (y)
  real, intent(in) :: x !< Scalar input value [nondim]
  real             :: y !< Scalar output value [nondim]

  y = max(x, 0.0) ! ReLU activation

end function activation_fn

!> Single application of ANN inference using vector input and output
!!
!! This implementation is the simplest using allocation and de-allocation
!! of temporary arrays
subroutine ANN_apply_vector_orig(x, y, CS)
  type(ANN_CS), intent(in)    :: CS                               !< ANN instance
  real,         intent(in)    :: x(CS%layer_sizes(1))             !< Inputs [arbitrary]
  real,         intent(inout) :: y(CS%layer_sizes(CS%num_layers)) !< Outputs [arbitrary]
  ! Local variables
  real, allocatable :: x_1(:), x_2(:) ! intermediate states [nondim]
  integer :: i, o ! Input, output indices

  ! Normalize input
  allocate(x_1(CS%layer_sizes(1)))
  do i = 1,CS%layer_sizes(1)
    x_1(i) = ( x(i) - CS%input_means(i) ) * CS%input_norms(i)
  enddo

  ! Apply Linear layers
  do i = 1, CS%num_layers-1
    allocate(x_2(CS%layer_sizes(i+1)))
    call layer_apply_orig(x_1, x_2, CS%layers(i))
    deallocate(x_1)
    allocate(x_1(CS%layer_sizes(i+1)))
    x_1(:) = x_2(:)
    deallocate(x_2)
  enddo

  ! Un-normalize output
  do o = 1, CS%layer_sizes(CS%num_layers)
    y(o) = ( x_1(o) * CS%output_norms(o) ) + CS%output_means(o)
  enddo

  deallocate(x_1)

  contains

  !> Applies linear layer to input data x and stores the result in y with
  !! y = A*x + b with optional application of the activation function so the
  !! overall operations is ReLU(A*x + b)
  subroutine layer_apply_orig(x, y, layer)
    type(layer_type), intent(in)    :: layer                 !< Linear layer
    real,             intent(in)    :: x(layer%input_width)  !< Input vector [nondim]
    real,             intent(inout) :: y(layer%output_width) !< Output vector [nondim]
    ! Local variables
    integer :: i, o ! Input, output indices

    ! Add bias
    y(:) = layer%b(:)
    ! Multiply by kernel
    do i=1,layer%input_width
      do o=1,layer%output_width
        y(o) = y(o) + x(i) * layer%A(o, i)
      enddo
    enddo
    ! Apply activation function
    if (layer%activation) y(:) = activation_fn(y(:))

  end subroutine layer_apply_orig
end subroutine ANN_apply_vector_orig

!> Single application of ANN inference using vector input and output
!!
!! This implementation avoids repeated reallocation of work arrays and uses the
!! output index for the fastest (inner-most) loop in the layer matrix multiply.
subroutine ANN_apply_vector_oi(x, y, CS)
  type(ANN_CS), intent(in)    :: CS                               !< ANN instance
  real,         intent(in)    :: x(CS%layer_sizes(1))             !< Inputs [arbitrary]
  real,         intent(inout) :: y(CS%layer_sizes(CS%num_layers)) !< Outputs [arbitrary]
  ! Local variables
  real, allocatable :: x_1(:), x_2(:) ! intermediate states [nondim]
  integer :: i, o ! Input, output indices

  allocate( x_1( maxval( CS%layer_sizes(:) ) ) )
  allocate( x_2( maxval( CS%layer_sizes(:) ) ) )

  ! Normalize input
  do i = 1,CS%layer_sizes(1)
    x_1(i) = ( x(i) - CS%input_means(i) ) * CS%input_norms(i)
  enddo

  ! Apply Linear layers
  do i = 1, CS%num_layers-2, 2
    call layer_apply_oi(x_1, x_2, CS%layers(i))
    call layer_apply_oi(x_2, x_1, CS%layers(i+1))
  enddo
  if (mod(CS%num_layers,2)==0) then
    call layer_apply_oi(x_1, x_2, CS%layers(CS%num_layers-1))
    ! Un-normalize output
    do o = 1, CS%layer_sizes(CS%num_layers)
      y(o) = x_2(o) * CS%output_norms(o) + CS%output_means(o)
    enddo
  else
    ! Un-normalize output
    do o = 1, CS%layer_sizes(CS%num_layers)
      y(o) = x_1(o) * CS%output_norms(o) + CS%output_means(o)
    enddo
  endif

  deallocate(x_1, x_2)

  contains

  !> Applies linear layer to input data x and stores the result in y with
  !! y = A*x + b with optional application of the activation function so the
  !! overall operations is ReLU(A*x + b)
  subroutine layer_apply_oi(x, y, layer)
    type(layer_type), intent(in)    :: layer                 !< Linear layer
    real,             intent(in)    :: x(layer%input_width)  !< Input vector [nondim]
    real,             intent(inout) :: y(layer%output_width) !< Output vector [nondim]
    ! Local variables
    integer :: i, o ! Input, output indices

    ! Add bias
    y(:) = layer%b(:)
    ! Multiply by kernel
    do i=1,layer%input_width
      do o=1,layer%output_width
        y(o) = y(o) + x(i) * layer%A(o, i)
      enddo
    enddo
    ! Apply activation function
    if (layer%activation) y(:) = activation_fn(y(:))

  end subroutine layer_apply_oi
end subroutine ANN_apply_vector_oi

!> Single application of ANN inference using array input and output
!! with (space,feature) indexing
!!
!! This implementation uses the space index for the fastest (inner-most) loop
!! in the layer matrix multiply, with the input index as the next fastest loop,
!! and uses the weights matrix A(output,index). It also applies the activation
!! function within the outer loop of the matrix multiply.
subroutine ANN_apply_array_sio(nij, x, y, CS)
  type(ANN_CS), intent(in)    :: CS !< ANN control structure
  integer,      intent(in)    :: nij !< Size of spatial dimension
  real,         intent(in)    :: x(nij, CS%layer_sizes(1)) !< input [arbitrary]
  real,         intent(inout) :: y(nij, CS%layer_sizes(CS%num_layers)) !< output [arbitrary]
  ! Local variables
  real, allocatable :: x_1(:,:), x_2(:,:) ! intermediate states [nondim]
  integer :: l, i, o ! Layer, input, output index

  allocate( x_1( nij, maxval( CS%layer_sizes(:) ) ) )
  allocate( x_2( nij, maxval( CS%layer_sizes(:) ) ) )

  ! Normalize input
  do i = 1, CS%layer_sizes(1)
    x_1(:,i) = ( x(:,i) - CS%input_means(i) ) * CS%input_norms(i)
  enddo

  ! Apply Linear layers
  do l = 1, CS%num_layers-2, 2
    call layer_apply_sio(nij, x_1, x_2, CS%layers(l))
    call layer_apply_sio(nij, x_2, x_1, CS%layers(l+1))
  enddo
  if (mod(CS%num_layers,2)==0) then
    call layer_apply_sio(nij, x_1, x_2, CS%layers(CS%num_layers-1))
    ! Un-normalize output
    do o = 1, CS%layer_sizes(CS%num_layers)
      y(:,o) = x_2(:,o) * CS%output_norms(o) + CS%output_means(o)
    enddo
  else
    ! Un-normalize output
    do o = 1, CS%layer_sizes(CS%num_layers)
      y(:,o) = x_1(:,o) * CS%output_norms(o) + CS%output_means(o)
    enddo
  endif

  deallocate(x_1, x_2)

  contains

  !> Applies linear layer to input data x and stores the result in y with
  !! y = A*x + b with optional application of the activation function so the
  !! overall operations is ReLU(A*x + b)
  subroutine layer_apply_sio(nij, x, y, layer)
    type(layer_type), intent(in)    :: layer !< Linear layer
    integer,          intent(in)    :: nij   !< Size of spatial dimension
    real,             intent(in)    :: x(nij, layer%input_width) !< Input vector [nondim]
    real,             intent(inout) :: y(nij, layer%output_width) !< Output vector [nondim]
    ! Local variables
    integer :: i, o ! Input, output indices

    do o = 1, layer%output_width
      ! Add bias
      y(:,o) = layer%b(o)
      ! Multiply by kernel
      do i = 1, layer%input_width
        y(:,o) = y(:,o) + x(:,i) * layer%A(o, i)
      enddo
      ! Apply activation function
      if (layer%activation) y(:,o) = activation_fn(y(:,o))
    enddo

  end subroutine layer_apply_sio
end subroutine ANN_apply_array_sio

!> Sets weights and bias for a single layer
subroutine set_layer(ANN, layer, weights, biases, activation)
  type(ANN_CS), intent(inout) :: ANN !< ANN control structure
  integer,      intent(in)    :: layer !< The number of the layer being adjusted
  real,         intent(in)    :: weights(:,:) !< The weights to assign
  real,         intent(in)    :: biases(:) !< The biases to assign
  logical,      intent(in)    :: activation !< Turn on the activation function

  if ( layer >= ANN%num_layers ) &
      call MOM_error(FATAL, "MOM_ANN, set_layer: layer is out of range")
  if ( layer < 1 ) &
      call MOM_error(FATAL, "MOM_ANN, set_layer: layer should be >= 1")

  if ( size(biases) /= size(ANN%layers(layer)%b) ) &
      call MOM_error(FATAL, "MOM_ANN, set_layer: mismatch in size of biases")
  ANN%layers(layer)%b(:) = biases(:)

  if ( size(weights,1) /= size(ANN%layers(layer)%A,1) ) &
      call MOM_error(FATAL, "MOM_ANN, set_layer: mismatch in size of weights (first dim)")
  if ( size(weights,2) /= size(ANN%layers(layer)%A,2) ) &
      call MOM_error(FATAL, "MOM_ANN, set_layer: mismatch in size of weights (second dim)")
  ANN%layers(layer)%A(:,:) = weights(:,:)

  ANN%layers(layer)%activation = activation
end subroutine set_layer

!> Sets input normalization
subroutine set_input_normalization(ANN, means, norms)
  type(ANN_CS),   intent(inout) :: ANN !< ANN control structure
  real, optional, intent(in)    :: means(:) !< The mean of each input
  real, optional, intent(in)    :: norms(:) !< The standard deviation of each input

  if (present(means)) then
    if ( size(means) /= size(ANN%input_means) ) &
        call MOM_error(FATAL, "MOM_ANN, set_input_normalization: mismatch in size of means")
    ANN%input_means(:) = means(:)
  endif

  if (present(norms)) then
    if ( size(norms) /= size(ANN%input_norms) ) &
        call MOM_error(FATAL, "MOM_ANN, set_input_normalization: mismatch in size of norms")
    ANN%input_norms(:) = norms(:)
  endif

end subroutine set_input_normalization

!> Sets output normalization
subroutine set_output_normalization(ANN, means, norms)
  type(ANN_CS),   intent(inout) :: ANN !< ANN control structure
  real, optional, intent(in)    :: means(:) !< The mean of each output
  real, optional, intent(in)    :: norms(:) !< The standard deviation of each output

  if (present(means)) then
    if ( size(means) /= size(ANN%output_means) ) &
        call MOM_error(FATAL, "MOM_ANN, set_output_normalization: mismatch in size of means")
    ANN%output_means(:) = means(:)
  endif

  if (present(norms)) then
    if ( size(norms) /= size(ANN%output_norms) ) &
        call MOM_error(FATAL, "MOM_ANN, set_output_normalization: mismatch in size of norms")
    ANN%output_norms(:) = norms(:)
  endif

end subroutine set_output_normalization

!> Create a random ANN
subroutine ANN_random(ANN, nlayers, widths)
  type(ANN_CS), intent(inout) :: ANN !< ANN control structure
  integer,      intent(in)    :: nlayers !< Number of layers
  integer,      intent(in)    :: widths(nlayers) !< Width of each layer
  ! Local variables
  integer :: l

  call ANN_allocate(ANN, nlayers, widths)

  do l = 1, nlayers-1
    call randomize_layer(ANN, nlayers, l, widths)
  enddo

end subroutine ANN_random

!> Fill a layer with random numbers
subroutine randomize_layer(ANN, nlayers, layer, widths)
  type(ANN_CS), intent(inout) :: ANN !< ANN control structure
  integer,      intent(in)    :: nlayers !< Number of layers
  integer,      intent(in)    :: layer !< Layer number to randomize
  integer,      intent(in)    :: widths(nlayers) !< Width of each layer
  ! Local variables
  real :: weights(widths(layer+1),widths(layer)) ! Weights
  real :: biases(widths(layer+1)) ! Biases

  call random_number(weights)
  weights(:,:) = 2. * weights(:,:) - 1.

  call random_number(biases)
  biases(:) = 2. * biases(:) - 1.

  call set_layer(ANN, layer, weights, biases, layer<nlayers-1)

end subroutine randomize_layer

!> Runs unit tests on ANN functions.
!!
!! Should only be called from a single/root thread.
!! Returns True if a test fails, otherwise False.
logical function ANN_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  type(ANN_CS) :: ANN ! An ANN
  type(testing) :: test ! Manage tests
  real, allocatable :: x(:), y(:), y_good(:), x2(:,:), y2(:,:) ! Inputs, outputs [arbitrary]
  integer, parameter :: max_rand_nlay = 10 ! Deepest random ANN to generate
  integer :: widths(max_rand_nlay) ! Number of layers for random ANN
  integer :: nlay ! Number of layers for random ANN
  integer :: i, iter ! Loop counters
  logical :: rand_res ! Status of random tests

  ANN_unit_tests = .false. ! Start by assuming all is well
  call test%set(verbose=verbose) ! Pass verbose mode to test

  ! Identity ANN for one input
  allocate( y(1) )
  call ANN_allocate(ANN, 2, [1,1])
  call set_layer(ANN, 1, reshape([1.],[1,1]), [0.], .false.)
  call ANN_apply([1.], y, ANN)
  call test%real_scalar(y(1), 1., 'Scalar identity')
  deallocate( y )
  call ANN_end(ANN)

  ! Summation ANN
  allocate( y(1) )
  call ANN_allocate(ANN, 2, [4,1])
  call set_layer(ANN, 1, reshape([1.,1.,1.,1.], [1,4]), [0.], .false.)
  call ANN_apply([-1.,0.,1.,2.], y, ANN)
  call test%real_scalar(y(1), 2., 'Summation')
  deallocate( y )
  call ANN_end(ANN)

  ! Identity ANN for vector input/output
  call ANN_allocate(ANN, 2, [3,3])
  allocate( y(3) )
  call set_layer(ANN, 1, reshape([1.,0.,0., &
                                  0.,1.,0., &
                                  0.,0.,1.], [3,3]), [0.,0.,0.], .false.)
  call ANN_apply([-1.,0.,1.], y, ANN)
  call test%real_arr(3, y, [-1.,0.,1.], 'Vector identity')
  deallocate( y )
  call ANN_end(ANN)

  ! Rectifying ANN for vector input/output
  allocate( y(3) )
  call ANN_allocate(ANN, 2, [3,3])
  call set_layer(ANN, 1, reshape([1.,0.,0., &
                                  0.,1.,0., &
                                  0.,0.,1.], [3,3]), [0.,0.,0.], .true.)
  call ANN_apply([-1.,0.,1.], y, ANN)
  call test%real_arr(3, y, [0.,0.,1.], 'Rectifier')
  deallocate( y )
  call ANN_end(ANN)

  ! The next 3 tests re-use the same network with 4 inputs, a 4-wide hidden layer, and one output
  allocate( y(1) )
  call ANN_allocate(ANN, 3, [4,4,1])

  ! 1 hidden layer: rectifier followed by summation
  ! Inputs: [-1,0,1,2]
  ! Rectified: [0,0,1,2]
  ! Sum: 3
  ! Outputs: 3
  call set_layer(ANN, 1, reshape([1.,0.,0.,0., &
                                  0.,1.,0.,0., &
                                  0.,0.,1.,0., &
                                  0.,0.,0.,1.], [4,4]), [0.,0.,0.,0.], .true.)
  call set_layer(ANN, 2, reshape([1.,1.,1.,1.], [1,4]), [0.], .false.)
  call ANN_apply_vector_orig([-1.,0.,1.,2.], y, ANN)
  call test%real_scalar(y(1), 3., 'Rectifier+summation')

  ! as above but with biases
  ! Inputs: [-2,-1,0,1]
  ! After bias: [-1,0,1,2] with b=1
  ! Rectified: [0,0,1,2]
  ! Sum: 3
  ! After bias: 6 with b=3
  ! Outputs: 6
  call set_layer(ANN, 1, reshape([1.,0.,0.,0., &
                                  0.,1.,0.,0., &
                                  0.,0.,1.,0., &
                                  0.,0.,0.,1.], [4,4]), [1.,1.,1.,1.], .true.)
  call set_layer(ANN, 2, reshape([1.,1.,1.,1.], [1,4]), [3.], .false.)
  call ANN_apply_vector_orig([-2.,-1.,0.,1.], y, ANN)
  call test%real_scalar(y(1), 6., 'Rectifier+summation+bias')

  ! as above but with normalization of inputs and outputs
  ! Inputs: [0,2,4,6]
  ! Normalized inputs: [-2,-1,0,1] (using mean=-4, norm=2)
  ! Normalized outputs: 6
  ! De-normalized output: 2 (using mean=-10, norm=2)
  call set_input_normalization(ANN, means=[4.,4.,4.,4.], norms=[0.5,0.5,0.5,0.5])
  call set_output_normalization(ANN, norms=[2.], means=[-10.])
  call ANN_apply_vector_orig([0.,2.,4.,6.], y, ANN)
  call test%real_scalar(y(1), 2., 'Rectifier+summation+bias+norms')

  deallocate( y )
  call ANN_end(ANN)

  ! as above with a 1x1 4th identity layer (to check loop combinations)
  allocate( y(1) )
  call ANN_allocate(ANN, 4, [4,4,1,1])
  call set_layer(ANN, 1, reshape([1.,0.,0.,0., &
                                  0.,1.,0.,0., &
                                  0.,0.,1.,0., &
                                  0.,0.,0.,1.], [4,4]), [1.,1.,1.,1.], .true.)
  call set_layer(ANN, 2, reshape([1.,1.,1.,1.], [1,4]), [3.], .false.)
  call set_layer(ANN, 3, reshape([1.],[1,1]), [0.], .false.)
  call set_input_normalization(ANN, means=[4.,4.,4.,4.], norms=[0.5,0.5,0.5,0.5])
  call set_output_normalization(ANN, norms=[2.], means=[-10.])
  call ANN_apply_vector_orig([0.,2.,4.,6.], y, ANN)
  call test%real_scalar(y(1), 2., 'Rectifier+summation+bias+norms 4-layer')

  ! as above with v2 of ANN_apply
  call ANN_apply_vector_oi([0.,2.,4.,6.], y, ANN)
  call test%real_scalar(y(1), 2., 'Rectifier+summation+bias+norms 4-layer v2')
  deallocate( y )

  allocate( y2(1,2) )
  ! as above with v5 of ANN_apply applied to 2d inputs, x(space,feature)
  call ANN_apply_array_sio(2, reshape([0.,1.,2.,3.,4.,5.,6.,7.],[2,4]), y2, ANN)
  call test%real_arr(2, y2, [2.,5.], 'Rectifier+summation+bias+norms 4-layer array v2')
  deallocate( y2 )

  call ANN_end(ANN)

  ! The following block checks that for random ANN (weights and layers widths)
  ! each of the various implementations of inference give identical results.
  ! This helped catch loop and allocation errors.
  rand_res = .false.
  do iter = 1, 1000
    allocate( y(max_rand_nlay+1) )
    call random_number(y) ! Vector of random numbers 0..1
    nlay = 2 + floor( y(max_rand_nlay+1) * ( max_rand_nlay - 1 ) ) ! 2 < nlay < max_rand_nlay
    widths(:) = 1 + floor( y(1:nlay) * 8 ) ! 1 < layer width < 8
    deallocate( y )
    call ANN_random(ANN, nlay, widths)
    allocate( x(widths(1)), y(widths(nlay)), y_good(widths(nlay)) )
    call ANN_apply_vector_orig(x, y_good, ANN)
    call ANN_apply_vector_oi(x, y, ANN)
    rand_res = rand_res .or. maxval( abs( y(:) - y_good(:) ) ) > 0. ! Check results from v2 = v1
    allocate( x2(20,widths(1)), y2(20,widths(nlay)) ) ! 2D input, output
    do i = 1, 20
      x2(i,:) = x(:)
    enddo
    call ANN_apply_array_sio(20, x2, y2, ANN)
    rand_res = rand_res .or. maxval( abs( maxval(y2(:,:),1) - y_good(:) ) ) > 0. ! Check results from array v2 = v1
    rand_res = rand_res .or. maxval( abs( minval(y2(:,:),1) - y_good(:) ) ) > 0. ! Check results from array v2 = v1
    deallocate( x, y, y_good, x2, y2 )
    call ANN_end(ANN)
  enddo
  call test%test(rand_res, 'Equivalence between inference variants with random results')

  ANN_unit_tests = test%summarize('ANN_unit_tests')

end function ANN_unit_tests

!> \namespace mom_ann
!!
!! The mom_ann module is a pure fortran implementation of fully-connected feed-forward
!! networks to facilitate easy evaluation of data-driven functions in MOM6. For performant
!! implementations or for novel architectires, using machine-learning libraries (e.g. via
!! mom_database_comms) are necessary, or at least likely to be more efficient.
!!
!! The artificial neural network (ANN) understood by this MOM6 module has \f$ N \f$ layers,
!! including the input-layer and output-layer, thus requireing \f$ N \geq 2\f$.
!!
!! The output values (neurons or nodes) of any layer other than the input layer (i.e. \f$ l>1 \f$) are
!! \f[
!! y_{l,j} = f_l( b_{l,j} + A_{l,j,i} x_{l-1,i} )
!! \f]
!! where \f$ f(x) = max(0, x) \f$ is the ReLU activation function, \f$b_{l,j}\f$ is a bias for each neuron,
!! \f$A_{l,j,i}\f$ are a rectangular matrix of weights for each layer, and \f$x_{l-1,i}\f$ are the outputs
!! of the previous layer, \f$l-1\f$. The subscript on \f$ f_l() \f$ indicates the activation function is
!! optional for each layer.
!!
!! Currently, the performance of various implementations is dependent on the shape/size of the network and
!! the size of input data. For this reason we provide several versions that all yield the same result but
!! for differently shaped inputs.
!!
!! \image html https://upload.wikimedia.org/wikipedia/commons/4/46/Colored_neural_network.svg
!! Fig: A three layer network with 3 inputs, 2 outputs, and 1 hidden layer. There are two rectanglar
!! matrices of weights (black arrows). The bias for each neuron is implied."

end module MOM_ANN
