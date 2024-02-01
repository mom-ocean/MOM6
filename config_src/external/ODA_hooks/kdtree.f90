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
