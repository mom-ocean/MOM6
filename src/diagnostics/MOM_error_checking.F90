module MOM_error_checking
!***********************************************************************
!*                   GNU General Public License                        *
!* This file is a part of MOM.                                         *
!*                                                                     *
!* MOM is free software; you can redistribute it and/or modify it and  *
!* are expected to follow the terms of the GNU General Public License  *
!* as published by the Free Software Foundation; either version 2 of   *
!* the License, or (at your option) any later version.                 *
!*                                                                     *
!* MOM is distributed in the hope that it will be useful, but WITHOUT  *
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  *
!* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    *
!* License for more details.                                           *
!*                                                                     *
!* For the full text of the GNU General Public License,                *
!* write to: Free Software Foundation, Inc.,                           *
!*           675 Mass Ave, Cambridge, MA 02139, USA.                   *
!* or see:   http://www.gnu.org/licenses/gpl.html                      *
!***********************************************************************

!********+*********+*********+*********+*********+*********+*********+**
!*                                                                     *
!*  By Robert Hallberg, August 2006                                    *
!*                                                                     *
!*    This file contains various subroutines that are useful for       *
!*  catching certain types of errors.                                  *
!*                                                                     *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_domains, only : pe_here, BGRID_NE, To_All, Scalar_Pair
use MOM_domains, only : create_group_pass, do_group_pass, group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_grid, only : ocean_grid_type

implicit none ; private

#include <MOM_memory.h>

public :: check_redundant

interface check_redundant
  module procedure check_redundant_v3d, check_redundant_v2d
  module procedure check_redundant_s3d, check_redundant_s2d
end interface check_redundant

contains

subroutine check_redundant_v3d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction, stagger)
  character(len=*),                    intent(in)    :: mesg
  real, dimension(NIMEMB_,NJMEM_,NKMEM_), intent(in) :: u_comp
  real, dimension(NIMEM_,NJMEMB_,NKMEM_), intent(in) :: v_comp
  type(ocean_grid_type),               intent(inout) :: G
  integer,                   optional, intent(in)    :: is, ie, js, je
  integer,                   optional, intent(in)    :: direction
  integer,                   optional, intent(in)    :: stagger
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.
!  (in/opt)  stagger - the stagger flag to be passed to pass_vector.

  character(len=24) :: mesg_k
  integer :: k
  
  do k=1,G%ke
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_v2d(trim(mesg)//trim(mesg_k), u_comp(:,:,k), &
             v_comp(:,:,k), G, is, ie, js, je, direction, stagger)
  enddo
end subroutine  check_redundant_v3d

subroutine check_redundant_v2d(mesg, u_comp, v_comp, G, is, ie, js, je, &
                               direction, stagger)
  character(len=*),                intent(in)    :: mesg
  real, dimension(NIMEMB_,NJMEM_), intent(in)    :: u_comp
  real, dimension(NIMEM_,NJMEMB_), intent(in)    :: v_comp
  type(ocean_grid_type),           intent(inout) :: G
  integer,               optional, intent(in)    :: is, ie, js, je
  integer,               optional, intent(in)    :: direction
  integer,               optional, intent(in)    :: stagger
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  direction - the direction flag to be passed to pass_vector.
!  (in/opt)  stagger - the stagger flag to be passed to pass_vector.

  real :: u_nonsym(SZI_(G),SZJ_(G))
  real :: v_nonsym(SZI_(G),SZJ_(G))
  real :: u_resym(SZIB_(G),SZJ_(G))
  real :: v_resym(SZI_(G),SZJB_(G))
  character(len=128) :: mesg2
  type(group_pass_type), save :: pass_nonsym_uv, pass_resym_uv !For group halo pass

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if ((isd == IsdB) .and. (jsd == JsdB)) return

  do i=isd,ied ; do j=jsd,jed
    u_nonsym(i,j) = u_comp(i,j) ; v_nonsym(i,j) = v_comp(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")

  call create_group_pass(pass_nonsym_uv, u_nonsym, v_nonsym, G%Domain_aux, &
                         direction, stagger)
  call create_group_pass(pass_resym_uv, u_resym, v_resym, G%Domain,  &
                         direction, stagger)
  call do_group_pass(pass_nonsym_uv, G%Domain_aux)

  do I=IsdB,IedB ; do j=jsd,jed ; u_resym(I,j) = u_comp(I,j) ; enddo ; enddo
  do i=isd,ied ; do J=JsdB,JedB ; v_resym(i,J) = v_comp(i,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    u_resym(i,j) = u_nonsym(i,j) ; v_resym(i,j) = v_nonsym(i,j)
  enddo ; enddo
  call do_group_pass(pass_resym_uv, G%Domain)

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je
  
  do i=is_ch,ie_ch ; do j=js_ch+1,je_ch
    if (u_resym(i,j) /= u_comp(i,j)) then
      write(mesg2,'(" redundant u-components",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           u_comp(i,j), u_resym(i,j),u_comp(i,j)-u_resym(i,j),i,j,pe_here()
      write(*,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo
  do i=is_ch+1,ie_ch ; do j=js_ch,je_ch
    if (v_resym(i,j) /= v_comp(i,j)) then
      write(mesg2,'(" redundant v-comps",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," x,y = ",2(1pe12.4)" on pe ",i4)') &
           v_comp(i,j), v_resym(i,j),v_comp(i,j)-v_resym(i,j),i,j, &
           G%geoLonBu(i,j), G%geoLatBu(i,j), pe_here()
      write(*,'(A155)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_v2d

subroutine check_redundant_s3d(mesg, array, G, is, ie, js, je, stagger)
  character(len=*),                     intent(in)    :: mesg
  real, dimension(NIMEMB_,NJMEMB_,NKMEM_), intent(in) :: array
  type(ocean_grid_type),                intent(inout) :: G
  integer,                    optional, intent(in)    :: is, ie, js, je
  integer,                    optional, intent(in)    :: stagger
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  stagger - the stagger flag to be passed to pass_vector.

  character(len=24) :: mesg_k
  integer :: k
  
  do k=1,G%ke
    if (k < 10) then ; write(mesg_k,'(" Layer",i2," ")') k
    elseif (k < 100) then ; write(mesg_k,'(" Layer",i3," ")') k
    elseif (k < 1000) then ; write(mesg_k,'(" Layer",i4," ")') k
    else ; write(mesg_k,'(" Layer",i9," ")') k ; endif

    call check_redundant_s2d(trim(mesg)//trim(mesg_k), array(:,:,k), &
                             G, is, ie, js, je, stagger)
  enddo
end subroutine  check_redundant_s3d


subroutine check_redundant_s2d(mesg, array, G, is, ie, js, je, stagger)
  character(len=*),                intent(in)    :: mesg
  real, dimension(NIMEMB_,NJMEMB_), intent(in)   :: array
  type(ocean_grid_type),           intent(inout) :: G
  integer,               optional, intent(in)    :: is, ie, js, je
  integer,               optional, intent(in)    :: stagger
! Arguments: u_comp - The u-component of the vector being checked.
!  (in)      v_comp - The v-component of the vector being checked.
!  (in)      mesg - A message indicating what is being checked.
!  (in)      G - The ocean's grid structure.
!  (in/opt)  is, ie, js, je - the i- and j- range of indices to check.
!  (in/opt)  stagger - the stagger flag to be passed to pass_vector.

  real :: a_nonsym(SZI_(G),SZJ_(G))
  real :: a_resym(SZIB_(G),SZJB_(G))
  character(len=128) :: mesg2
  type(group_pass_type), save :: pass_a_nonsym, pass_a_resym !For group halo pass

  integer :: i, j, is_ch, ie_ch, js_ch, je_ch
  integer :: Isq, Ieq, Jsq, Jeq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if ((isd == IsdB) .and. (jsd == JsdB)) return

  do i=isd,ied ; do j=jsd,jed
    a_nonsym(i,j) = array(i,j)
  enddo ; enddo

  if (.not.associated(G%Domain_aux)) call MOM_error(FATAL," check_redundant"//&
    " called with a non-associated auxiliary domain the grid type.")
  call create_group_pass(pass_a_nonsym, a_nonsym, a_nonsym, G%Domain_aux, &
                         direction=To_All+Scalar_Pair, stagger=BGRID_NE)
  call create_group_pass(pass_a_resym, a_resym, a_resym, G%Domain,     &
                         direction=To_All+Scalar_Pair, stagger=BGRID_NE)
  call do_group_pass(pass_a_nonsym, G%Domain_aux)

  do I=IsdB,IedB ; do J=JsdB,JedB ; a_resym(I,J) = array(I,J) ; enddo ; enddo
  do i=isd,ied ; do j=jsd,jed
    a_resym(i,j) = a_nonsym(i,j)
  enddo ; enddo
  call do_group_pass(pass_a_resym, G%Domain) 

  is_ch = Isq ; ie_ch = Ieq ; js_ch = Jsq ; je_ch = Jeq
  if (present(is)) is_ch = is ; if (present(ie)) ie_ch = ie
  if (present(js)) js_ch = js ; if (present(js)) je_ch = je
  
  do i=is_ch,ie_ch ; do j=js_ch,je_ch
    if (a_resym(i,j) /= array(i,j)) then
      write(mesg2,'(" Redundant points",2(1pe12.4)," differ by ", &
                    & 1pe12.4," at i,j = ",2i4," on pe ",i4)') &
           array(i,j), a_resym(i,j),array(i,j)-a_resym(i,j),i,j,pe_here()
      write(*,'(A130)') trim(mesg)//trim(mesg2)
    endif
  enddo ; enddo

end subroutine  check_redundant_s2d

end module MOM_error_checking
