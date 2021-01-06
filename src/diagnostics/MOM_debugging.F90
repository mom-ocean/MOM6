!> Provides checksumming functions for debugging
!!
!! This module contains subroutines that perform various error checking and
!! debugging functions for MOM6.  This routine is similar to it counterpart in
!! the SIS2 code, except for the use of the ocean_grid_type and by keeping them
!! separate we retain the ability to set up MOM6 and SIS2 debugging separately.
module MOM_debugging

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_checksums, only : hchksum, Bchksum, qchksum, uvchksum, hchksum_pair
use MOM_checksums, only : is_NaN, chksum, MOM_checksums_init
use MOM_coms, only : PE_here, root_PE, num_PEs
use MOM_coms, only : min_across_PEs, max_across_PEs, reproducing_sum
use MOM_domains, only : pass_vector, pass_var, pe_here
use MOM_domains, only : BGRID_NE, AGRID, To_All, Scalar_Pair
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : log_version, param_file_type, get_param
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type

implicit none ; private

public :: check_redundant_C, check_redundant_B, check_redundant_T, check_redundant
public :: vec_chksum, vec_chksum_C, vec_chksum_B, vec_chksum_A
public :: MOM_debugging_init, totalStuff, totalTandS
public :: check_column_integral, check_column_integrals

! These interfaces come from MOM_checksums.
public :: hchksum, Bchksum, qchksum, is_NaN, chksum, uvchksum, hchksum_pair

!> Check for consistency between the duplicated points of a C-grid vector
interface check_redundant
  module procedure check_redundant_vC3d, check_redundant_vC2d
end interface check_redundant
!> Check for consistency between the duplicated points of a C-grid vector
interface check_redundant_C
  module procedure check_redundant_vC3d, check_redundant_vC2d
end interface check_redundant_C
!> Check for consistency between the duplicated points of a B-grid vector or scalar
interface check_redundant_B
  module procedure check_redundant_vB3d, check_redundant_vB2d
  module procedure check_redundant_sB3d, check_redundant_sB2d
end interface check_redundant_B
!> Check for consistency between the duplicated points of an A-grid vector or scalar
interface check_redundant_T
  module procedure check_redundant_sT3d, check_redundant_sT2d
  module procedure check_redundant_vT3d, check_redundant_vT2d
end interface check_redundant_T

!> Do checksums on the components of a C-grid vector
interface vec_chksum
  module procedure chksum_vec_C3d, chksum_vec_C2d
end interface vec_chksum
!> Do checksums on the components of a C-grid vector
interface vec_chksum_C
  module procedure chksum_vec_C3d, chksum_vec_C2d
end interface vec_chksum_C
!> Do checksums on the components of a B-grid vector
interface vec_chksum_B
  module procedure chksum_vec_B3d, chksum_vec_B2d
end interface vec_chksum_B
!> Do checksums on the components of an A-grid vector
interface vec_chksum_A
  module procedure chksum_vec_A3d, chksum_vec_A2d
end interface vec_chksum_A

! Note: these parameters are module data but ONLY used when debugging and
!       so can violate the thread-safe requirement of no module/global data.
integer :: max_redundant_prints = 100 !< Maximum number of times to write redundant messages
integer :: redundant_prints(3) = 0 !< Counters for controlling redundant printing
logical :: debug = .false. !< Write out verbose debugging data
logical :: debug_chksums = .true. !< Perform checksums on arrays
logical :: debug_redundant = .true. !< Check redundant values on PE boundaries

contains

!> MOM_debugging_init initializes the MOM_debugging module, and sets
!! the parameterts that control which checks are active for MOM6.
subroutine MOM_debugging_init(param_file)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_debugging" ! This module's name.

  call log_version(param_file, mdl, version, debugging=.true.)
  call get_param(param_file, mdl, "DEBUG", debug, &
                 "If true, write out verbose debugging data.", &
                 default=.false., debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_CHKSUMS", debug_chksums, &
                 "If true, checksums are performed on arrays in the "//&
                 "various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)
  call get_param(param_file, mdl, "DEBUG_REDUNDANT", debug_redundant, &
                 "If true, debug redundant data points during calls to "//&
                 "the various vec_chksum routines.", default=debug, &
                 debuggingParam=.true.)

  call MOM_checksums_init(param_file)

end subroutine MOM_debugging_init

!> Check for consistency between the duplicated points of a 3-D C-grid vector
subroutine check_redundant_vC3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                    intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),               intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:),   intent(in)    :: u_comp !< The u-component of the vector
                                                               !! to be checked for consistency
  real, dimension(G%isd:,G%JsdB:,:),   intent(in)    :: v_comp !< The u-component of the vector
                                                               !! to be checked for consistency
  integer,                   optional, intent(in)    :: is     !< The starting i-index to check
  integer,                   optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                   optional, intent(in)    :: js     !< The starting j-index to check
  integer,                   optional, intent(in)    :: je     !< The ending j-index to check
  integer,                   optional, intent(in)    :: direction !< the direction flag to be
                                                               !! passed to pass_vector
  ! Local variables
  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(u_comp,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_vC2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction)
  enddo
end subroutine  check_redundant_vC3d

!> Check for consistency between the duplicated points of a 2-D C-grid vector
subroutine check_redundant_vC2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector
                                                           !! to be checked for consistency
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The u-component of the vector
                                                           !! to be checked for consistency
  integer,               optional, intent(in)    :: is     !< The starting i-index to check
  integer,               optional, intent(in)    :: ie     !< The ending i-index to check
  integer,               optional, intent(in)    :: js     !< The starting j-index to check
  integer,               optional, intent(in)    :: je     !< The ending j-index to check
  integer,               optional, intent(in)    :: direction !< the direction flag to be
                                                           !! passed to pass_vector
  ! Local variables
  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: u_resym(G%IsdB:G%IedB,G%jsd:G%jed)
  real :: v_resym(G%isd:G%ied,G%JsdB:G%JedB)
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, direction)

  do I=IsdB,IedB ; do j=jsd,jed ; u_resym(I,j) = u_comp(I,j) ; enddo ; enddo
  do i=isd,ied ; do J=JsdB,JedB ; v_resym(i,J) = v_comp(i,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call pass_vector(u_resym, v_resym, G%Domain, direction)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_resym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(3) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(3) = redundant_prints(3) + 1
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(3) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(3) = redundant_prints(3) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_vC2d

!> Check for consistency between the duplicated points of a 3-D scalar at corner points
subroutine check_redundant_sB3d(mesg, array, G, is, ie, js, je)
  character(len=*),                     intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),                intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:),   intent(in)    :: array !< The array to be checked for consistency
  integer,                    optional, intent(in)    :: is    !< The starting i-index to check
  integer,                    optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                    optional, intent(in)    :: js    !< The starting j-index to check
  integer,                    optional, intent(in)    :: je    !< The ending j-index to check

  ! Local variables
  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(array,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_sB2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                              G, is, ie, js, je)
  enddo
end subroutine  check_redundant_sB3d

!> Check for consistency between the duplicated points of a 2-D scalar at corner points
subroutine check_redundant_sB2d(mesg, array, G, is, ie, js, je)
  character(len=*),                 intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: array !< The array to be checked for consistency
  integer,                optional, intent(in)    :: is    !< The starting i-index to check
  integer,                optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                optional, intent(in)    :: js    !< The starting j-index to check
  integer,                optional, intent(in)    :: je    !< The ending j-index to check
  ! Local variables
  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: a_resym(G%IsdB:G%IedB,G%JsdB:G%JedB)
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call pass_vector(a_nonsym, a_nonsym, G%Domain_aux, &
                   direction=To_All+Scalar_Pair, stagger=BGRID_NE)

  do I=IsdB,IedB ; do J=JsdB,JedB ; a_resym(I,J) = array(I,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    a_resym(i,j) = a_nonsym(i,j)
  enddo ; enddo
  call pass_vector(a_resym, a_resym, G%Domain, direction=To_All+Scalar_Pair, &
                   stagger=BGRID_NE)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_resym(i,j) /= array(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_resym(i,j),array(i,j)-a_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_sB2d

!> Check for consistency between the duplicated points of a 3-D B-grid vector
subroutine check_redundant_vB3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                    intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),               intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in)    :: u_comp !< The u-component of the vector
                                                           !! to be checked for consistency
  real, dimension(G%IsdB:,G%JsdB:,:),  intent(in)    :: v_comp !< The v-component of the vector
                                                           !! to be checked for consistency
  integer,                   optional, intent(in)    :: is     !< The starting i-index to check
  integer,                   optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                   optional, intent(in)    :: js     !< The starting j-index to check
  integer,                   optional, intent(in)    :: je     !< The ending j-index to check
  integer,                   optional, intent(in)    :: direction !< the direction flag to be
                                                               !! passed to pass_vector
  ! Local variables
  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(u_comp,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_vB2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction)
  enddo
end subroutine  check_redundant_vB3d

!> Check for consistency between the duplicated points of a 2-D B-grid vector
subroutine check_redundant_vB2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                                direction)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: u_comp !< The u-component of the vector
                                                            !! to be checked for consistency
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector
                                                            !! to be checked for consistency
  integer,                optional, intent(in)    :: is     !< The starting i-index to check
  integer,                optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                optional, intent(in)    :: js     !< The starting j-index to check
  integer,                optional, intent(in)    :: je     !< The ending j-index to check
  integer,                optional, intent(in)    :: direction !< the direction flag to be
                                                            !! passed to pass_vector
  ! Local variables
  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: u_resym(G%IsdB:G%IedB,G%JsdB:G%JedB)
  real :: v_resym(G%IsdB:G%IedB,G%JsdB:G%JedB)
  character(len=128) :: mesg2
  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.(present(is) .or. present(ie) .or. present(js) .or. present(je))) then
    ! This only works with symmetric memory, so otherwise return.
    if ((isd == IsdB) .and. (jsd == JsdB)) return
  endif

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call pass_vector(u_nonsym, v_nonsym, G%Domain_aux, direction, stagger=BGRID_NE)

  do I=IsdB,IedB ; do J=JsdB,JedB
    u_resym(I,J) = u_comp(I,J) ; v_resym(I,J) = v_comp(I,J)
  enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call pass_vector(u_resym, v_resym, G%Domain, direction, stagger=BGRID_NE)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (u_resym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo
  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(2) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(2) = redundant_prints(2) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_vB2d

!> Check for consistency between the duplicated points of a 3-D scalar at tracer points
subroutine check_redundant_sT3d(mesg, array, G, is, ie, js, je)
  character(len=*),                     intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),                intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:),     intent(in)    :: array !< The array to be checked for consistency
  integer,                    optional, intent(in)    :: is    !< The starting i-index to check
  integer,                    optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                    optional, intent(in)    :: js    !< The starting j-index to check
  integer,                    optional, intent(in)    :: je    !< The ending j-index to check
  ! Local variables
  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(array,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_sT2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                              G, is, ie, js, je)
  enddo
end subroutine  check_redundant_sT3d


!> Check for consistency between the duplicated points of a 2-D scalar at tracer points
subroutine check_redundant_sT2d(mesg, array, G, is, ie, js, je)
  character(len=*),                 intent(in)    :: mesg  !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G     !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:),   intent(in)    :: array !< The array to be checked for consistency
  integer,                optional, intent(in)    :: is    !< The starting i-index to check
  integer,                optional, intent(in)    :: ie    !< The ending i-index to check
  integer,                optional, intent(in)    :: js    !< The starting j-index to check
  integer,                optional, intent(in)    :: je    !< The ending j-index to check
  ! Local variables
  real :: a_nonsym(G%isd:G%ied,G%jsd:G%jed)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  is_ch = G%isc ; ie_ch = G%iec ; js_ch = G%jsc ; je_ch = G%jec
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  ! This only works on points outside of the standard computational domain.
  if ((is_ch == G%isc) .and. (ie_ch == G%iec) .and. &
      (js_ch == G%jsc) .and. (je_ch == G%jec)) return

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  call pass_var(a_nonsym, G%Domain)

  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_nonsym(i,j) /= array(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_nonsym(i,j),array(i,j)-a_nonsym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_sT2d

!> Check for consistency between the duplicated points of a 3-D A-grid vector
subroutine check_redundant_vT3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                    intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),               intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:),    intent(in)    :: u_comp !< The u-component of the vector
                                                           !! to be checked for consistency
  real, dimension(G%isd:,G%jsd:,:),    intent(in)    :: v_comp !< The v-component of the vector
                                                           !! to be checked for consistency
  integer,                   optional, intent(in)    :: is     !< The starting i-index to check
  integer,                   optional, intent(in)    :: ie     !< The ending i-index to check
  integer,                   optional, intent(in)    :: js     !< The starting j-index to check
  integer,                   optional, intent(in)    :: je     !< The ending j-index to check
  integer,                   optional, intent(in)    :: direction !< the direction flag to be
                                                           !! passed to pass_vector
  ! Local variables
  character(len=24) :: mesg_k
  integer :: k

  do k=1,size(u_comp,3)
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_vT2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction)
  enddo
end subroutine  check_redundant_vT3d

!> Check for consistency between the duplicated points of a 2-D A-grid vector
subroutine check_redundant_vT2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:),  intent(in)    :: u_comp !< The u-component of the vector
                                                           !! to be checked for consistency
  real, dimension(G%isd:,G%jsd:),  intent(in)    :: v_comp !< The v-component of the vector
                                                           !! to be checked for consistency
  integer,               optional, intent(in)    :: is     !< The starting i-index to check
  integer,               optional, intent(in)    :: ie     !< The ending i-index to check
  integer,               optional, intent(in)    :: js     !< The starting j-index to check
  integer,               optional, intent(in)    :: je     !< The ending j-index to check
  integer,               optional, intent(in)    :: direction !< the direction flag to be
                                                           !! passed to pass_vector
  ! Local variables
  real :: u_nonsym(G%isd:G%ied,G%jsd:G%jed)
  real :: v_nonsym(G%isd:G%ied,G%jsd:G%jed)
  character(len=128) :: mesg2

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  is_ch = G%isc ; ie_ch = G%iec ; js_ch = G%jsc ; je_ch = G%jec
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je

  ! This only works on points outside of the standard computational domain.
  if ((is_ch == G%isc) .and. (ie_ch == G%iec) .and. &
      (js_ch == G%jsc) .and. (je_ch == G%jec)) return

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  call pass_vector(u_nonsym, v_nonsym, G%Domain, direction, stagger=AGRID)

  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_nonsym(i,j) /= u_comp(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_nonsym(i,j),u_comp(i,j)-u_nonsym(i,j),i,j,pe_here()
      write(0,'(A130)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_nonsym(i,j) /= v_comp(i,j) .and. &
        redundant_prints(1) < max_redundant_prints) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_nonsym(i,j),v_comp(i,j)-v_nonsym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(0,'(A155)') trim(mesg)//trim(mesg2)
      redundant_prints(1) = redundant_prints(1) + 1
    endif
  enddo ; enddo

end subroutine  check_redundant_vT2d

!> Do a checksum and redundant point check on a 3d C-grid vector.
subroutine chksum_vec_C3d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                  intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),             intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                 optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                 optional, intent(in)    :: scalars !< If true this is a pair of
                                                             !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call uvchksum(mesg, u_comp, v_comp, G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_C3d

!> Do a checksum and redundant point check on a 2d C-grid vector.
subroutine chksum_vec_C2d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),           intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector
  integer,               optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,               optional, intent(in)    :: scalars !< If true this is a pair of
                                                           !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call uvchksum(mesg, u_comp, v_comp, G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_C(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_C(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_C2d

!> Do a checksum and redundant point check on a 3d B-grid vector.
subroutine chksum_vec_B3d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                   intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),              intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%IsdB:,G%JsdB:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                  optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                  optional, intent(in)    :: scalars !< If true this is a pair of
                                                              !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call Bchksum(u_comp, mesg//"(u)", G%HI, halos)
    call Bchksum(v_comp, mesg//"(v)", G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_B3d

! Do a checksum and redundant point check on a 2d B-grid vector.
subroutine chksum_vec_B2d(mesg, u_comp, v_comp, G, halos, scalars, symmetric)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%IsdB:,G%JsdB:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  logical,                optional, intent(in)    :: symmetric !< If true, do the checksums on the
                                                            !! full symmetric computational domain.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call Bchksum(u_comp, mesg//"(u)", G%HI, halos, symmetric=symmetric)
    call Bchksum(v_comp, mesg//"(v)", G%HI, halos, symmetric=symmetric)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_B(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_B(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_B2d

!> Do a checksum and redundant point check on a 3d C-grid vector.
subroutine chksum_vec_A3d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),                 intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),            intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%jsd:,:), intent(in)    :: v_comp !< The v-component of the vector
  integer,                optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,                optional, intent(in)    :: scalars !< If true this is a pair of
                                                            !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call hchksum(u_comp, mesg//"(u)", G%HI, halos)
    call hchksum(v_comp, mesg//"(v)", G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_A3d

!> Do a checksum and redundant point check on a 2d C-grid vector.
subroutine chksum_vec_A2d(mesg, u_comp, v_comp, G, halos, scalars)
  character(len=*),               intent(in)    :: mesg   !< An identifying message
  type(ocean_grid_type),          intent(inout) :: G      !< The ocean's grid structure
  real, dimension(G%isd:,G%jsd:), intent(in)    :: u_comp !< The u-component of the vector
  real, dimension(G%isd:,G%jsd:), intent(in)    :: v_comp !< The v-component of the vector
  integer,              optional, intent(in)    :: halos  !< The width of halos to check (default 0)
  logical,              optional, intent(in)    :: scalars !< If true this is a pair of
                                                          !! scalars that are being checked.
  ! Local variables
  logical :: are_scalars
  are_scalars = .false. ; if (present(scalars)) are_scalars = scalars

  if (debug_chksums) then
    call hchksum(u_comp, mesg//"(u)", G%HI, halos)
    call hchksum(v_comp, mesg//"(v)", G%HI, halos)
  endif
  if (debug_redundant) then
    if (are_scalars) then
      call check_redundant_T(mesg, u_comp, v_comp, G, direction=To_All+Scalar_Pair)
    else
      call check_redundant_T(mesg, u_comp, v_comp, G)
    endif
  endif

end subroutine chksum_vec_A2d

!> This function returns the sum over computational domain of all
!! processors of hThick*stuff, where stuff is a 3-d array at tracer points.
function totalStuff(HI, hThick, areaT, stuff)
  type(hor_index_type),               intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights
  real, dimension(HI%isd:,HI%jsd:),   intent(in) :: areaT  !< The array of cell areas [m2]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: stuff  !< The array of stuff to be summed
  real                                         :: totalStuff !< the globally integrated amoutn of stuff
  ! Local variables
  real, dimension(HI%isc:HI%iec, HI%jsc:HI%jec) :: tmp_for_sum
  integer :: i, j, k, nz

  nz = size(hThick,3)
  tmp_for_sum(:,:) = 0.0
  do k=1,nz ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    tmp_for_sum(i,j) = tmp_for_sum(i,j) + hThick(i,j,k) * stuff(i,j,k) * areaT(i,j)
  enddo ; enddo ; enddo
  totalStuff = reproducing_sum(tmp_for_sum)

end function totalStuff

!> This subroutine display the total thickness, temperature and salinity
!! as well as the change since the last call.
subroutine totalTandS(HI, hThick, areaT, temperature, salinity, mesg)
  type(hor_index_type),               intent(in) :: HI     !< A horizontal index type
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: hThick !< The array of thicknesses to use as weights
  real, dimension(HI%isd:,HI%jsd:),   intent(in) :: areaT  !< The array of cell areas [m2]
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: temperature !< The temperature field to sum
  real, dimension(HI%isd:,HI%jsd:,:), intent(in) :: salinity    !< The salinity field to sum
  character(len=*),                   intent(in) :: mesg        !< An identifying message
  ! NOTE: This subroutine uses "save" data which is not thread safe and is purely for
  ! extreme debugging without a proper debugger.
  real, save :: totalH = 0., totalT = 0., totalS = 0.
  ! Local variables
  logical, save :: firstCall = .true.
  real, dimension(HI%isc:HI%iec, HI%jsc:HI%jec) :: tmp_for_sum
  real :: thisH, thisT, thisS, delH, delT, delS
  integer :: i, j, k, nz

  nz = size(hThick,3)
  tmp_for_sum(:,:) = 0.0
  do k=1,nz ; do j=HI%jsc,HI%jec ; do i=HI%isc,HI%iec
    tmp_for_sum(i,j) = tmp_for_sum(i,j) + hThick(i,j,k) * areaT(i,j)
  enddo ; enddo ; enddo
  thisH = reproducing_sum(tmp_for_sum)
  thisT = totalStuff(HI, hThick, areaT, temperature)
  thisS = totalStuff(HI, hThick, areaT, salinity)

  if (is_root_pe()) then
    if (firstCall) then
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Totals H,T,S:',thisH,thisT,thisS,' ',mesg
      firstCall = .false.
    else
      delH = thisH - totalH
      delT = thisT - totalT
      delS = thisS - totalS
      totalH = thisH ; totalT = thisT ; totalS = thisS
      write(0,*) 'Tot/del H,T,S:',thisH,thisT,thisS,delH,delT,delS,' ',mesg
    endif
  endif

end subroutine totalTandS

!> Returns false if the column integral of a given quantity is within roundoff
logical function check_column_integral(nk, field, known_answer)
  integer,             intent(in) :: nk           !< Number of levels in column
  real, dimension(nk), intent(in) :: field        !< Field to be summed
  real, optional,      intent(in) :: known_answer !< If present is the expected sum,
                                                  !! If missing, assumed zero
  ! Local variables
  real    :: u_sum, error, expected
  integer :: k

  u_sum = field(1)
  error = 0.

  ! Reintegrate and sum roundoff errors
  do k=2,nk
    u_sum = u_sum + field(k)
    error = error + EPSILON(u_sum)*MAX(ABS(u_sum),ABS(field(k)))
  enddo

  ! Assign expected answer to either the optional input or 0
  if (present(known_answer)) then
    expected = known_answer
  else
    expected = 0.
  endif

  ! Compare the column integrals against calculated roundoff error
  if (abs(u_sum-expected) > error) then
    check_column_integral = .true.
  else
    check_column_integral = .false.
  endif

end function check_column_integral

!> Returns false if the column integrals of two given quantities are within roundoff of each other
logical function check_column_integrals(nk_1, field_1, nk_2, field_2, missing_value)
  integer,               intent(in) :: nk_1           !< Number of levels in field 1
  integer,               intent(in) :: nk_2           !< Number of levels in field 2
  real, dimension(nk_1), intent(in) :: field_1        !< First field to be summed
  real, dimension(nk_2), intent(in) :: field_2        !< Second field to be summed
  real, optional,        intent(in) :: missing_value  !< If column contains missing values,
                                                      !! mask them from the sum
  ! Local variables
  real    :: u1_sum, error1, u2_sum, error2, misval
  integer :: k

  ! Assign missing value
  if (present(missing_value)) then
    misval = missing_value
  else
    misval = 0.
  endif

  u1_sum = field_1(1)
  error1 = 0.

  ! Reintegrate and sum roundoff errors
  do k=2,nk_1
    if (field_1(k)/=misval) then
      u1_sum = u1_sum + field_1(k)
      error1 = error1 + EPSILON(u1_sum)*MAX(ABS(u1_sum),ABS(field_1(k)))
    endif
  enddo

  u2_sum = field_2(1)
  error2 = 0.

  ! Reintegrate and sum roundoff errors
  do k=2,nk_2
    if (field_2(k)/=misval) then
      u2_sum = u2_sum + field_2(k)
      error2 = error2 + EPSILON(u2_sum)*MAX(ABS(u2_sum),ABS(field_2(k)))
    endif
  enddo

  ! Compare the column integrals against calculated roundoff error
  if (abs(u1_sum-u2_sum) > (error1+error2)) then
    check_column_integrals = .true.
  else
    check_column_integrals = .false.
  endif

end function check_column_integrals

end module MOM_debugging
