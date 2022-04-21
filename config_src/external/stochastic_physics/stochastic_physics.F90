! The are stubs for ocean stochastic physics
! the fully functional code is available at
! http://github.com/noaa-psd/stochastic_physics
module stochastic_physics

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, WARNING

implicit none ; private

public :: init_stochastic_physics_ocn
public :: run_stochastic_physics_ocn

contains

!> Initializes the stochastic physics perturbations.
subroutine init_stochastic_physics_ocn(delt, geoLonT, geoLatT, nx, ny, nz, pert_epbl_in, do_sppt_in, &
                                       mpiroot, mpicomm, iret)
  real,    intent(in)    :: delt !< timestep in seconds between calls to run_stochastic_physics_ocn [s]
  integer, intent(in)    :: nx   !< number of gridpoints in the x-direction of the compute grid
  integer, intent(in)    :: ny   !< number of gridpoints in the y-direction of the compute grid
  integer, intent(in)    :: nz   !< number of gridpoints in the z-direction of the compute grid
  real,    intent(in)    :: geoLonT(nx,ny) !< Longitude in degrees
  real,    intent(in)    :: geoLatT(nx,ny) !< Latitude in degrees
  logical, intent(in)    :: pert_epbl_in !< logical flag, if true generate random pattern for ePBL perturbations
  logical, intent(in)    :: do_sppt_in   !< logical flag, if true generate random pattern for SPPT perturbations
  integer, intent(in)    :: mpiroot !< root processor
  integer, intent(in)    :: mpicomm !< mpi communicator
  integer, intent(out)   :: iret    !< return code

  iret=0
  if (pert_epbl_in) then
    call MOM_error(WARNING, 'init_stochastic_physics_ocn: pert_epbl needs to be false if using the stub')
    iret=-1
  endif
  if (do_sppt_in) then
    call MOM_error(WARNING, 'init_stochastic_physics_ocn: do_sppt needs to be false if using the stub')
    iret=-1
  endif

  ! This stub function does not actually do anything.
  return
end subroutine init_stochastic_physics_ocn

!> Determines the stochastic physics perturbations.
subroutine run_stochastic_physics_ocn(sppt_wts, t_rp1, t_rp2)
  real, intent(inout) :: sppt_wts(:,:) !< array containing random weights for SPPT range [0,2]
  real, intent(inout) :: t_rp1(:,:)    !< array containing random weights for ePBL
                                       !! perturbations (KE generation) range [0,2]
  real, intent(inout) :: t_rp2(:,:)    !< array containing random weights for ePBL
                                       !! perturbations (KE dissipation) range [0,2]

  ! This stub function does not actually do anything.
  return
end subroutine run_stochastic_physics_ocn

end module stochastic_physics
