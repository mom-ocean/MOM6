!> Provides buffers that can dynamically grow as needed. These are primarily intended for the
!! diagnostics which need to store intermediate or partial states of state variables
module MOM_diag_buffers

use MOM_io, only : stdout, stderr

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public :: diag_buffer_unit_tests_2d, diag_buffer_unit_tests_3d

type, abstract :: buffer_base
end type buffer_base

!> Holds a 2d field
type, extends(buffer_base) :: buffer_2d
  real, dimension(:,:), allocatable :: field !< The actual 2d field to be stored [arbitrary]
end type buffer_2d

!> Holds a 3d field
type, extends(buffer_base) :: buffer_3d
  real, dimension(:,:,:), allocatable :: field !< The actual 3d field to be stored [arbitrary]
end type buffer_3d

!> The base class for the diagnostic buffers in this module
type, abstract :: diag_buffer_base ; private
  integer :: is !< The start slot of the array i-direction
  integer :: js !< The start slot of the array j-direction
  integer :: ie !< The end slot of the array i-direction
  integer :: je !< The end slot of the array j-direction
  real :: fill_value = 0. !< Set the fill value to use when growing the buffer [arbitrary]

  integer, allocatable, dimension(:) :: ids  !< List of diagnostic ids whose slot corresponds to the row in the buffer
  integer :: length = 0 !< The number of slots in the buffer

  contains

  procedure(a_grow), deferred :: grow !< Increase the size of the buffer
  procedure, public :: set_fill_value !< Set the fill value to use when growing the buffer
  procedure, public :: check_capacity_by_id !< Check the size size of the buffer and increase if necessary
  procedure, public :: set_horizontal_extents !< Define the horizontal extents of the arrays
  procedure, public :: mark_available !< Mark that a slot in the buffer can be reused
  procedure, public :: grow_ids !< Increase the size of the vector storing diagnostic ids
  procedure, public :: find_buffer_slot !< Find the slot corresponding to a specific diagnostic id
end type diag_buffer_base

!> Dynamically growing buffer for 2D arrays.
type, extends(diag_buffer_base), public :: diag_buffer_2d ; private
  type(buffer_2d), public, dimension(:), allocatable :: buffer !< The actual 2D buffer which will dynamically grow

  contains

  procedure, public :: grow => grow_2d !< Increase the size of the buffer
  procedure, public :: store => store_2d !< Store a field in the buffer, increasing as necessary
  procedure, public :: set_extents_from_array => set_extents_from_array_2d !< Set extents from array bounds
end type diag_buffer_2d

!> Dynamically growing buffer for 3D arrays.
type, extends(diag_buffer_base), public :: diag_buffer_3d ; private
  type(buffer_3d), public, dimension(:), allocatable :: buffer !< The actual 2D buffer which will dynamically grow
  integer :: ks !< The start slot in the k-dimension
  integer :: ke !< The last slot in the k-dimension

  contains

  procedure, public :: set_vertical_extent !< Set the vertical extents of the buffer
  procedure, public :: grow => grow_3d !< Increase the size of the buffer
  procedure, public :: store => store_3d !< Store a field in the buffer, increasing as necessary
  procedure, public :: set_extents_from_array => set_extents_from_array_3d !< Set extents from array bounds
end type diag_buffer_3d

contains

!> Signature for the grow methods on n-dimension diagnostic buffer types
subroutine a_grow(this)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
end subroutine

!> Set the fill value to use when growing the buffer
subroutine set_fill_value(this, fill_value)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
  real,                    intent(in)    :: fill_value !< The fill value to use when growing the buffer [arbitrary]

  this%fill_value = fill_value
end subroutine set_fill_value

!> Mark a slot in the buffer as unused based on a diagnostic id. For example,
!! the data in that slot has already been consumed and can thus be overwritten
subroutine mark_available(this, id)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
  integer,                 intent(in)    :: id   !< The diagnostic id
  integer :: slot

  slot = this%find_buffer_slot(id)
  this%ids(slot) = 0
end subroutine mark_available

!> Return the slot of the buffer corresponding to the diagnostic id
pure function find_buffer_slot(this, id) result(slot)
  class(diag_buffer_base), intent(in) :: this !< The diagnostic buffer
  integer, intent(in) :: id !< The diagnostic id

  integer, dimension(1) :: temp
  integer :: slot !< The slot in the buffer corresponding to the diagnostic id

  if (allocated(this%ids)) then
    !NOTE: Alternatively could do slot = SUM(findloc(...))
    temp = findloc(this%ids(:), id)
    slot = temp(1)
  else
    slot = 0
  endif

end function find_buffer_slot

!> Grow the ids array by one
subroutine grow_ids(this)
  class(diag_buffer_base), intent(inout) :: this !< This buffer

  integer, allocatable, dimension(:) :: temp
  integer :: n

  n = this%length

  allocate(temp(n+1))
  if (n>0) temp(1:n) = this%ids(:)
  call move_alloc(temp, this%ids)
end subroutine grow_ids

!> Check whether the id already has a slot reserved. If not, find a new empty slot and if
!! need be, grow the buffer.
impure function check_capacity_by_id(this, id) result(slot)
  class(diag_buffer_base), intent(inout) :: this !< This 2d buffer
  integer,                 intent(in)    :: id   !< The diagnostic id
  integer :: slot

  slot = this%find_buffer_slot(id)
  if (slot==0) then
    ! Check to see if there is an open slot
    if (allocated(this%ids)) slot = this%find_buffer_slot(0)
    ! If slot is still 0, then the buffer must grow
    if (slot==0) then
      call this%grow()
      slot = this%length
    endif
    this%ids(slot) = id
  endif
end function check_capacity_by_id

!> Set the horizontal extents of the buffer
subroutine set_horizontal_extents(this, is, ie, js, je)
  class(diag_buffer_base), intent(inout) :: this !< The diagnostic buffer
  integer,               intent(in)    :: is !< The start slot of the array i-direction
  integer,               intent(in)    :: ie !< The end slot of the array i-direction
  integer,               intent(in)    :: js !< The start slot of the array j-direction
  integer,               intent(in)    :: je !< The end slot of the array j-direction

  this%is = is ; this%ie = ie ; this%js = js ; this%je = je
end subroutine set_horizontal_extents

!> Set the vertical extent of the buffer
subroutine set_vertical_extent(this, ks, ke)
  class(diag_buffer_3d), intent(inout) :: this !< The diagnostic buffer
  integer,               intent(in)    :: ks !< The start slot of the array i-direction
  integer,               intent(in)    :: ke !< The end slot of the array i-direction

  this%ks = ks ; this%ke = ke
end subroutine set_vertical_extent

!> Set the extents of a 2D buffer from the bounds of a 2D array
subroutine set_extents_from_array_2d(this, array, fill_value_in)
  class(diag_buffer_2d), intent(inout) :: this !< The diagnostic buffer
  real, dimension(:,:), intent(in)     :: array !< The array whose bounds define the buffer extents
  real, optional,       intent(in)     :: fill_value_in !< Optional fill value

  call this%set_horizontal_extents(lbound(array,1), ubound(array,1), &
                                   lbound(array,2), ubound(array,2))
  if (present(fill_value_in)) call this%set_fill_value(fill_value_in)
end subroutine set_extents_from_array_2d

!> Set the extents of a 3D buffer from the bounds of a 3D array
subroutine set_extents_from_array_3d(this, array, fill_value_in)
  class(diag_buffer_3d), intent(inout) :: this !< The diagnostic buffer
  real, dimension(:,:,:), intent(in)   :: array !< The array whose bounds define the buffer extents
  real, optional,         intent(in)   :: fill_value_in !< Optional fill value

  call this%set_horizontal_extents(lbound(array,1), ubound(array,1), &
                                   lbound(array,2), ubound(array,2))
  call this%set_vertical_extent(lbound(array,3), ubound(array,3))
  if (present(fill_value_in)) call this%set_fill_value(fill_value_in)
end subroutine set_extents_from_array_3d

!> Grow a 2d diagnostic buffer
subroutine grow_2d(this)
  class(diag_buffer_2d), intent(inout) :: this

  integer :: i, n
  integer :: is, ie, js, je
  type(buffer_2d), dimension(:), allocatable :: new_buffer

  ! Grow the ID array
  call this%grow_ids()

  is = this%is ; ie = this%ie ; js = this%js ; je = this%je
  n = this%length

  allocate(new_buffer(n+1))
  do i=1,n
    allocate(new_buffer(i)%field(is:ie,js:je))
    new_buffer(i)%field(:,:) = this%buffer(i)%field(:,:)
  enddo
  allocate(new_buffer(n+1)%field(is:ie,js:je), source=this%fill_value)
  call move_alloc(new_buffer, this%buffer)
  this%length = n+1

end subroutine grow_2d

!> Store a 2D array into this buffer
subroutine store_2d(this, data, id)
  class(diag_buffer_2d), intent(inout) :: this !< This 2d buffer
  real, dimension(:,:),  intent(in)    :: data !< The data to be stored in the buffer [arbitrary]
  integer,               intent(in)    :: id !< The diagnostic id

  integer :: slot

  slot = this%check_capacity_by_id(id)
  this%buffer(slot)%field(:,:) = data(:,:)
end subroutine store_2d

!> Grow a 2d diagnostic buffer
subroutine grow_3d(this)
  class(diag_buffer_3d), intent(inout) :: this

  integer :: i, n
  integer :: is, ie, js, je, ks, ke
  type(buffer_3d), dimension(:), allocatable :: new_buffer

  ! Grow the ID array
  call this%grow_ids()

  is = this%is ; ie = this%ie ; js = this%js ; je = this%je ; ks = this%ks ; ke = this%ke
  n = this%length

  allocate(new_buffer(n+1))
  do i=1,n
    allocate(new_buffer(i)%field(is:ie,js:je,ks:ke))
    new_buffer(i)%field(:,:,:) = this%buffer(i)%field(:,:,:)
  enddo
  allocate(new_buffer(n+1)%field(is:ie,js:je,ks:ke), source=this%fill_value)
  call move_alloc(new_buffer, this%buffer)
  this%length = n+1

end subroutine grow_3d

!> Store a 3d array into this buffer
subroutine store_3d(this, data, id)
  class(diag_buffer_3d),  intent(inout) :: this !< This 3d buffer
  real, dimension(:,:,:), intent(in)    :: data !< The data to be stored in the buffer [arbitrary]
  integer,                intent(in)    :: id !< The diagnostic id

  integer :: slot

  ! Find the first slot in the ids array that is 0, i.e. this is a portion of the buffer that can be reused
  slot = this%check_capacity_by_id(id)
  this%buffer(slot)%field(:,:,:) = data(:,:,:)
end subroutine store_3d

!> Unit tests for the 2d version of the diag buffer
function diag_buffer_unit_tests_2d(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

  fail = .false.
  write(stdout,*) '==== MOM_diag_buffers: diag_buffers_unit_tests_2d ==='
  fail = fail .or. new_buffer_2d()
  fail = fail .or. grow_buffer_2d()
  fail = fail .or. fill_value_2d()
  fail = fail .or. store_buffer_2d()
  fail = fail .or. reuse_buffer_2d()

  contains

  !> Ensure properties of a newly initialized buffer
  function new_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    local_fail = .false.
    local_fail = local_fail .or. allocated(buffer%buffer)
    if (verbose) write(stdout,*) "new_buffer_2d: ", local_fail
    local_fail = local_fail .or. allocated(buffer%ids)
    if (verbose) write(stdout,*) "new_buffer_2d: ", local_fail
    local_fail = local_fail .or. buffer%length /= 0
    if (verbose) write(stdout,*) "new_buffer_2d: ", local_fail
  end function new_buffer_2d

  !> Test the growing of a buffer
  function grow_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6
    integer :: i

    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    ! Grow the buffer 3 times
    do i=1,3
      call buffer%grow()
      local_fail = local_fail .or. (buffer%length /= i)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 1) /= is)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 1) /= ie)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 2) /= js)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 2) /= je)
    enddo
    if (verbose) write(stdout,*) "grow_buffer_2d: ", local_fail
  end function grow_buffer_2d

  !> Test that growing new buffer fills the array with a set fill value
  function fill_value_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6
    real, parameter :: fill_value = -123.456
    integer :: i


    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    call buffer%set_fill_value(fill_value)
    ! Grow the buffer 3 times
    call buffer%grow()
    if (any(buffer%buffer(1)%field(:,:) /= fill_value)) local_fail = .true.
    if (verbose) write(stdout,*) "fill_value_2d: ", local_fail
  end function fill_value_2d

  !> Test storing a buffer based on a unique id
  function store_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, nlen=3
    integer :: i, slot
    real, allocatable, dimension(:,:,:) :: test_2d

    local_fail = .false.

    allocate(test_2d(nlen, is:ie, js:je))
    call random_number(test_2d)
    buffer%is = is
    buffer%ie = ie
    buffer%js = js
    buffer%je = je

    do i=1,nlen
      call buffer%store(test_2d(i,:,:), i*3)
      slot = buffer%find_buffer_slot(i*3)
      local_fail = local_fail .or. ANY(buffer%buffer(slot)%field(:,:) /= test_2d(i,:,:))
    enddo

    if (verbose) write(stdout,*) "store_buffer_2d: ", local_fail
  end function store_buffer_2d

  !> Test the reuse of a buffer. Fill it first like store_buffer_2d. Then,
  !! loop through again, but use the slots of the buffer in the following
  !! order: 2, 1, 3
  function reuse_buffer_2d() result(local_fail)
    type(diag_buffer_2d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, nlen=3
    integer :: i, new_i, id, new_id
    real, dimension(nlen, is:ie, js:je) :: test_2d_first, test_2d_second
    integer, dimension(nlen) :: reorder = [2,1,3]

    local_fail = .false.
    call random_number(test_2d_first)
    call random_number(test_2d_second)

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)

    do i=1,nlen
      call buffer%store(test_2d_first(i,:,:), id=i*3)
    enddo

    do i=1,nlen
      new_i = reorder(i)
      ! id and new_id are multiplied by primes to make sure they are unique
      id = reorder(i)*3
      new_id = i*7
      call buffer%mark_available(id=reorder(i)*3)
      call buffer%store(test_2d_second(i,:,:), id=new_id)
      local_fail = local_fail .or. buffer%find_buffer_slot(new_id) /= new_i
      test_2d_first(new_i,:,:) = test_2d_second(i,:,:)
    enddo
    local_fail = local_fail .or. any(buffer%ids /= [14, 7, 21])
    do i=1,nlen
      local_fail = local_fail .or. any(buffer%buffer(i)%field(:,:) /= test_2d_first(i,:,:))
    enddo
    if (verbose) write(stdout,*) "reuse_buffer_2d: ", local_fail
  end function reuse_buffer_2d

end function diag_buffer_unit_tests_2d

!> Test the 3d version of the buffer
function diag_buffer_unit_tests_3d(verbose) result(fail)
  logical, intent(in) :: verbose !< If true, write results to stdout
  logical :: fail !< True if any of the unit tests fail

  fail = .false.
  write(stdout,*) '==== MOM_diag_buffers: diag_buffers_unit_tests_3d ==='
  fail = fail .or. new_buffer_3d()
  fail = fail .or. grow_buffer_3d()
  fail = fail .or. fill_value_3d()
  fail = fail .or. store_buffer_3d()
  fail = fail .or. reuse_buffer_3d()

  contains

  !> Ensure properties of a newly initialized buffer
  function new_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    local_fail = .false.
    local_fail = local_fail .or. allocated(buffer%buffer)
    local_fail = local_fail .or. allocated(buffer%ids)
    local_fail = local_fail .or. buffer%length /= 0
    if (verbose) write(stdout,*) "new_buffer_3d: ", local_fail
  end function new_buffer_3d

  !> Test the growing of a buffer
  function grow_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10
    integer :: i

    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    call buffer%set_vertical_extent(ks=ks, ke=ke)
    ! Grow the buffer 3 times
    do i=1,3
      call buffer%grow()
      local_fail = local_fail .or. (buffer%length /= i)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 1) /= is)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 1) /= ie)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 2) /= js)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 2) /= je)
      local_fail = local_fail .or. (lbound(buffer%buffer(i)%field, 3) /= ks)
      local_fail = local_fail .or. (ubound(buffer%buffer(i)%field, 3) /= ke)
    if (verbose) write(stdout,*) "grow_buffer_3d: ", local_fail
    enddo
    if (verbose) write(stdout,*) "grow_buffer_3d: ", local_fail
  end function grow_buffer_3d

  !> Test that growing new buffer fills the array with a set fill value
  function fill_value_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail
    integer, parameter :: is=1, ie=2, js=3, je=6
    real, parameter :: fill_value = -123.456
    integer :: i


    local_fail = .false.

    call buffer%set_horizontal_extents(is=is, ie=ie, js=js, je=je)
    call buffer%set_fill_value(fill_value)
    ! Grow the buffer 3 times
    call buffer%grow()
    if (any(buffer%buffer(1)%field(:,:,:) /= fill_value)) local_fail = .true.
    if (verbose) write(stdout,*) "fill_value_3d: ", local_fail
  end function fill_value_3d

  !> Test storing a buffer based on a unique id
  function store_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10, nlen=3
    integer :: i, slot
    real, dimension(nlen,is:ie,js:je,ks:ke) :: test_3d

    local_fail = .false.
    call random_number(test_3d)
    buffer%is = is
    buffer%ie = ie
    buffer%js = js
    buffer%je = je
    buffer%ks = ks
    buffer%ke = ke

    do i=1,nlen
      call buffer%store(test_3d(i,:,:,:), i*3)
      slot = buffer%find_buffer_slot(i*3)
      local_fail = local_fail .or. ANY(buffer%buffer(slot)%field(:,:,:) /= test_3d(i,:,:,:))
    enddo

    if (verbose) write(stdout,*) "store_buffer_3d: ", local_fail
  end function store_buffer_3d

  !> Test the reuse of a buffer. Fill it first like store_buffer_3d. Then,
  !! loop through again, but use the slots of the buffer in the following
  !! order: 2, 1, 3
  function reuse_buffer_3d() result(local_fail)
    type(diag_buffer_3d) :: buffer
    logical :: local_fail !< True if any of the unit tests fail

    integer, parameter :: is=1, ie=2, js=3, je=6, ks=1, ke=10, nlen=3
    integer :: i, new_i, id, new_id
    real, dimension(nlen, is:ie, js:je, ks:ke) :: test_3d_first, test_3d_second
    integer, dimension(nlen) :: reorder = [2,1,3]

    local_fail = .false.
    call random_number(test_3d_first)
    call random_number(test_3d_second)

    buffer%is = is
    buffer%ie = ie
    buffer%js = js
    buffer%je = je
    buffer%ks = ks
    buffer%ke = ke

    do i=1,nlen
      call buffer%store(test_3d_first(i,:,:,:), id=i*3)
    enddo

    do i=1,nlen
      new_i = reorder(i)
      ! id and new_id are multiplied by primes to make sure they are unique
      id = reorder(i)*3
      new_id = i*7
      call buffer%mark_available(id=reorder(i)*3)
      call buffer%store(test_3d_second(i,:,:,:), id=new_id)
      local_fail = local_fail .or. buffer%find_buffer_slot(new_id) /= new_i
      test_3d_first(new_i,:,:,:) = test_3d_second(i,:,:,:)
    enddo
    local_fail = local_fail .or. any(buffer%ids /= [14, 7, 21])
    do i=1,nlen
      local_fail = local_fail .or. any(buffer%buffer(i)%field(:,:,:) /= test_3d_first(i,:,:,:))
    enddo
    if (verbose) write(stdout,*) "reuse_buffer_3d: ", local_fail
  end function reuse_buffer_3d

end function diag_buffer_unit_tests_3d

end module MOM_diag_buffers

