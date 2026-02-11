! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A null version of K-d tree from geoKdTree
module kdtree
  implicit none
  private

  public :: kd_root

  !> A K-d tree tpe
  type kd_root
    integer :: dummy !< To stop a compiler from doing nothing
  end type kd_root
end module kdtree
