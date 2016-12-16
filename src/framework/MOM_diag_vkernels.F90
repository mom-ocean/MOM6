!> Provides kernels for single-column interpolation, re-integration (re-mapping of integrated quantities)
!! and intensive-variable remapping in the vertical
module MOM_diag_vkernels

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public diag_vkernels_unit_tests
public interpolate_column
public reintegrate_column

! Module data used only for debugging
logical :: verbose_diag_vkernels_unit_tests = .true. ! If true, always write out unit tests results

contains

!> Linearly interpolate interface data, u_src, from grid h_src to a grid h_dest
subroutine interpolate_column(nsrc, h_src, u_src, ndest, h_dest, missing_value, u_dest)
  integer,                  intent(in)    :: nsrc !< Number of source cells
  real, dimension(nsrc),    intent(in)    :: h_src !< Thickness of source cells
  real, dimension(nsrc+1),  intent(in)    :: u_src !< Values at source cell interfaces
  integer,                  intent(in)    :: ndest !< Number of destination cells
  real, dimension(ndest),   intent(in)    :: h_dest !< Thickness of destination cells
  real,                     intent(in)    :: missing_value !< Value to assign in vanished cells
  real, dimension(ndest+1), intent(inout) :: u_dest !< Interpolated value at destination cell interfaces
  ! Local variables
  real :: x_dest ! Relative position of target interface
  real :: dh ! Source cell thickness
  real :: u1, u2 ! Values to interpolate between
  real :: weight_a, weight_b ! Weights for interpolation
  integer :: k_src, k_dest ! Index of cell in src and dest columns
  logical :: still_vanished ! Used for figuring out what to mask as missing

  ! Initial values for the loop
  still_vanished = .true.

  ! The following forces the "do while" loop to do one cycle that will set u1, u2, dh.
  k_src = 0
  dh = 0.
  x_dest = 0.

  do k_dest=1, ndest+1
    do while (dh<=x_dest .and. k_src<nsrc)
      ! Move positions pointers forward until the interval 0 .. dh spans x_dest.
      x_dest = x_dest - dh
      k_src = k_src + 1
      dh = h_src(k_src) ! Thickness of layer k_src

      ! Values that will be used for the interpolation.
      u1 = u_src(k_src) ! Value on left of source cell
      u2 = u_src(k_src+1) ! Value on right of source cell
    enddo

    if (dh>0.) then
      weight_a = max(0., ( dh - x_dest ) / dh) ! Weight of u1
      weight_b = min(1., x_dest / dh) ! Weight of u2
      u_dest(k_dest) = weight_a * u1 + weight_b * u2 ! Linear interpolation between u1 and u2
    else
      u_dest(k_dest) = 0.5 * ( u1 + u2 ) ! For a vanished layer we need to do something reasonable...
    endif

    ! Mask vanished layers at the surface which would be under an ice-shelf.
    ! TODO: Need to figure out what to do for an isopycnal coordinate diagnostic that could
    !       also have vanished layers at the surface.
    if (k_dest<=ndest) then
      x_dest = x_dest + h_dest(k_dest) ! Position of interface k_dest+1
      if (still_vanished .and. h_dest(k_dest)==0.) then
        ! When the layer k_dest is vanished and all layers above are also vanished, the k_dest
        ! interface value should be missing.
        u_dest(k_dest) = missing_value
      else
        still_vanished = .false.
      endif
    endif

  enddo

  ! Mask vanished layers on topography
  still_vanished = .true.
  do k_dest=ndest, 1, -1
    if (still_vanished .and. h_dest(k_dest)==0.) then
      ! When the layer k_dest is vanished and all layers below are also vanished, the k_dest+1
      ! interface value should be missing.
      u_dest(k_dest+1) = missing_value
    else
      exit
    endif
  enddo

end subroutine interpolate_column

!> Conservatively calculate integrated data, uh_dest, on grid h_dest, from layer-integrated data, uh_src, on grid h_src
subroutine reintegrate_column(nsrc, h_src, uh_src, ndest, h_dest, missing_value, uh_dest)
  integer,                intent(in)    :: nsrc !< Number of source cells
  real, dimension(nsrc),  intent(in)    :: h_src !< Thickness of source cells
  real, dimension(nsrc),  intent(in)    :: uh_src !< Values at source cell interfaces
  integer,                intent(in)    :: ndest !< Number of destination cells
  real, dimension(ndest), intent(in)    :: h_dest !< Thickness of destination cells
  real,                   intent(in)    :: missing_value !< Value to assign in vanished cells
  real, dimension(ndest), intent(inout) :: uh_dest !< Interpolated value at destination cell interfaces
  ! Local variables
  real :: x_dest ! Relative position of target interface
  real :: h_src_rem, h_dest_rem, dh ! Incremental thicknesses
  real :: uh_src_rem, uh_dest_rem, duh ! Incremental amounts of stuff
  integer :: k_src, k_dest ! Index of cell in src and dest columns
  integer :: iter
  logical :: src_ran_out, src_exists

  uh_dest(:) = missing_value

  k_src = 0
  k_dest = 0
  h_dest_rem = 0.
  h_src_rem = 0.
  src_ran_out = .false.
  src_exists = .false.

  do while(.true.)
    if (h_src_rem==0. .and. k_src<nsrc) then
      ! Supply is empty so move to the next source cell
      k_src = k_src + 1
      h_src_rem = h_src(k_src)
      uh_src_rem = uh_src(k_src)
      if (h_src_rem==0.) cycle
      src_exists = .true. ! This stops us masking out the entire column
    endif
    if (h_dest_rem==0. .and. k_dest<ndest) then
      ! Sink has no capacity so move to the next destination cell
      k_dest = k_dest + 1
      h_dest_rem = h_dest(k_dest)
      uh_dest(k_dest) = 0.
      if (h_dest_rem==0.) cycle
    endif
    if (k_src==nsrc .and. h_src_rem==0.) then
      if (src_ran_out) exit ! This is the second time implying there is no more src
      src_ran_out = .true.
      cycle
    endif
    duh = 0.
    if (h_src_rem<h_dest_rem) then
      ! The source cell is fully within the destination cell
      dh = h_src_rem
      if (dh>0.) duh = uh_src_rem
      h_src_rem = 0.
      uh_src_rem = 0.
      h_dest_rem = max(0., h_dest_rem - dh)
    elseif (h_src_rem>h_dest_rem) then
      ! Only part of the source cell can be used up
      dh = h_dest_rem
      duh = (dh / h_src_rem) * uh_src_rem
      h_src_rem = max(0., h_src_rem - dh)
      uh_src_rem = uh_src_rem - duh
      h_dest_rem = 0.
    else ! h_src_rem==h_dest_rem
      ! The source cell exactly fits the destination cell
      duh = uh_src_rem
      h_src_rem = 0.
      uh_src_rem = 0.
      h_dest_rem = 0.
    endif
    uh_dest(k_dest) = uh_dest(k_dest) + duh
    if (k_dest==ndest .and. (k_src==nsrc .or. h_dest_rem==0.)) exit
  enddo

  if (.not. src_exists) uh_dest(1:ndest) = missing_value

end subroutine reintegrate_column

!> Returns true if any unit tests for module MOM_diag_vkernels fail
logical function diag_vkernels_unit_tests()
  ! Local variables
  real, parameter :: missing_value=-9.999999999E9 ! Value to use for vanished layers
  logical :: fail

  write(0,*) '============ MOM_diag_kernels: diag_vkernels_unit_tests =================='
  write(0,*) '- - - - - - - - - - interpolation tests  - - - - - - - - - - - - - - - - -'

  fail = test_interp('Identity: 3 layer', &
                     3, (/1.,2.,3./), (/1.,2.,3.,4./), &
                     3, (/1.,2.,3./), (/1.,2.,3.,4./) )
  diag_vkernels_unit_tests = fail

  fail = test_interp('A: 3 layer to 2', &
                     3, (/1.,1.,1./), (/1.,2.,3.,4./), &
                     2, (/1.5,1.5/), (/1.,2.5,4./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('B: 2 layer to 3', &
                     2, (/1.5,1.5/), (/1.,4.,7./), &
                     3, (/1.,1.,1./), (/1.,3.,5.,7./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('C: 3 layer (vanished middle) to 2', &
                     3, (/1.,0.,2./), (/1.,2.,2.,3./), &
                     2, (/1.,2./), (/1.,2.,3./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('D: 3 layer (deep) to 3', &
                     3, (/1.,2.,3./), (/1.,2.,4.,7./), &
                     2, (/2.,2./), (/1.,3.,5./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('E: 3 layer to 3 (deep)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     3, (/2.,3.,4./), (/1.,3.,6.,8./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('F: 3 layer to 4 with vanished top/botton', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,5.,0./), (/missing_value,1.,3.,8.,missing_value/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('Fs: 3 layer to 4 with vanished top/botton (shallow)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,4.,0./), (/missing_value,1.,3.,7.,missing_value/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_interp('Fd: 3 layer to 4 with vanished top/botton (deep)', &
                     3, (/1.,2.,4./), (/1.,2.,4.,8./), &
                     4, (/0.,2.,6.,0./), (/missing_value,1.,3.,8.,missing_value/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  write(0,*) '- - - - - - - - - - reintegration tests  - - - - - - - - - - - - - - - - -'

  fail = test_reintegrate('Identity: 3 layer', &
                     3, (/1.,2.,3./), (/-5.,2.,1./), &
                     3, (/1.,2.,3./), (/-5.,2.,1./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('A: 3 layer to 2', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,3./), (/-4.,2./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('A: 3 layer to 2 (deep)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,4./), (/-4.,2./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('A: 3 layer to 2 (shallow)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     2, (/3.,2./), (/-4.,1.5/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('B: 3 layer to 4 with vanished top/bottom', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     4, (/0.,3.,3.,0./), (/0.,-4.,2.,0./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('C: 3 layer to 4 with vanished top//middle/bottom', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     5, (/0.,3.,0.,3.,0./), (/0.,-4.,0.,2.,0./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('D: 3 layer to 3 (vanished)', &
                     3, (/2.,2.,2./), (/-5.,2.,1./), &
                     3, (/0.,0.,0./), (/0.,0.,0./) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('D: 3 layer (vanished) to 3', &
                     3, (/0.,0.,0./), (/-5.,2.,1./), &
                     3, (/2.,2.,2./), (/missing_value, missing_value, missing_value/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('D: 3 layer (vanished) to 3 (vanished)', &
                     3, (/0.,0.,0./), (/-5.,2.,1./), &
                     3, (/0.,0.,0./), (/missing_value, missing_value, missing_value/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  fail = test_reintegrate('D: 3 layer (vanished) to 3 (vanished)', &
                     3, (/0.,0.,0./), (/0.,0.,0./), &
                     3, (/0.,0.,0./), (/missing_value, missing_value, missing_value/) )
  diag_vkernels_unit_tests = diag_vkernels_unit_tests .or. fail

  write(0,*) '=========================================================================='

  contains

  !> Returns true if a test of interpolate_column() produces the wrong answer
  logical function test_interp(msg, nsrc, h_src, u_src, ndest, h_dest, u_true)
    character(len=*),         intent(in) :: msg !< Message to label test
    integer,                  intent(in) :: nsrc !< Number of source cells
    real, dimension(nsrc),    intent(in) :: h_src !< Thickness of source cells
    real, dimension(nsrc+1),  intent(in) :: u_src !< Values at source cell interfaces
    integer,                  intent(in) :: ndest !< Number of destination cells
    real, dimension(ndest),   intent(in) :: h_dest !< Thickness of destination cells
    real, dimension(ndest+1), intent(in) :: u_true !< Correct value at destination cell interfaces
    ! Local variables
    real, dimension(ndest+1) :: u_dest ! Interpolated value at destination cell interfaces
    integer :: k
    real :: error
    logical :: print_results

    ! Interpolate from src to dest
    call interpolate_column(nsrc, h_src, u_src, ndest, h_dest, missing_value, u_dest)

    test_interp = .false.
    do k=1,ndest+1
      if (u_dest(k)/=u_true(k)) test_interp = .true.
    enddo
    if (verbose_diag_vkernels_unit_tests .or. test_interp) then
      write(0,'(2a)') ' Test: ',msg
      write(0,'(a3,3(a24))') 'k','u_result','u_true','error'
      do k=1,ndest+1
        error = u_dest(k)-u_true(k)
        if (error==0.) then
          write(0,'(i3,3(1pe24.16))') k,u_dest(k),u_true(k),u_dest(k)-u_true(k)
        else
          write(0,'(i3,3(1pe24.16),x,a)') k,u_dest(k),u_true(k),u_dest(k)-u_true(k),'<--- WRONG!'
        endif
      enddo
    endif
  end function test_interp

  !> Returns true if a test of reintegrate_column() produces the wrong answer
  logical function test_reintegrate(msg, nsrc, h_src, uh_src, ndest, h_dest, uh_true)
    character(len=*),       intent(in) :: msg !< Message to label test
    integer,                intent(in) :: nsrc !< Number of source cells
    real, dimension(nsrc),  intent(in) :: h_src !< Thickness of source cells
    real, dimension(nsrc),  intent(in) :: uh_src !< Values of source cell stuff
    integer,                intent(in) :: ndest !< Number of destination cells
    real, dimension(ndest), intent(in) :: h_dest !< Thickness of destination cells
    real, dimension(ndest), intent(in) :: uh_true !< Correct value of destination cell stuff
    ! Local variables
    real, dimension(ndest) :: uh_dest ! Reintegrated value on destination cells
    integer :: k
    real :: error
    logical :: print_results

    ! Interpolate from src to dest
    call reintegrate_column(nsrc, h_src, uh_src, ndest, h_dest, missing_value, uh_dest)

    test_reintegrate = .false.
    do k=1,ndest
      if (uh_dest(k)/=uh_true(k)) test_reintegrate = .true.
    enddo
    if (verbose_diag_vkernels_unit_tests .or. test_reintegrate) then
      write(0,'(2a)') ' Test: ',msg
      write(0,'(a3,3(a24))') 'k','uh_result','uh_true','error'
      do k=1,ndest
        error = uh_dest(k)-uh_true(k)
        if (error==0.) then
          write(0,'(i3,3(1pe24.16))') k,uh_dest(k),uh_true(k),uh_dest(k)-uh_true(k)
        else
          write(0,'(i3,3(1pe24.16),x,a)') k,uh_dest(k),uh_true(k),uh_dest(k)-uh_true(k),'<--- WRONG!'
        endif
      enddo
    endif
  end function test_reintegrate

end function diag_vkernels_unit_tests

end module MOM_diag_vkernels
