! The are stubs for ocean stochastic physics
! the fully functional code is available at
! http://github.com/noaa-psd/stochastic_physics
module stochastic_physics

implicit none

private

public :: init_stochastic_physics_ocn
public :: run_stochastic_physics_ocn

contains

!!!!!!!!!!!!!!!!!!!!
subroutine init_stochastic_physics_ocn(delt,geoLonT,geoLatT,nx,ny,nz,pert_epbl_in,do_sppt_in, &
                                       mpiroot, mpicomm, iret)
implicit none
real,intent(in)       :: delt !< timestep in seconds between calls to run_stochastic_physics_ocn
integer,intent(in)    :: nx   !< number of gridpoints in the x-direction of the compute grid
integer,intent(in)    :: ny   !< number of gridpoints in the y-direction of the compute grid
integer,intent(in)    :: nz   !< number of gridpoints in the z-direction of the compute grid
real,intent(in)       :: geoLonT(nx,ny) !< Longitude in degrees
real,intent(in)       :: geoLatT(nx,ny) !< Latitude in degrees
logical,intent(in)    :: pert_epbl_in !< logical flag, if true generate random pattern for ePBL perturbations
logical,intent(in)    :: do_sppt_in   !< logical flag, if true generate random pattern for SPPT perturbations
integer,intent(in)    :: mpiroot !< root processor
integer,intent(in)    :: mpicomm !< mpi communicator
integer, intent(out)  :: iret    !< return code

iret=0
if (pert_epbl_in .EQV. .true. ) then
   print*,'pert_epbl needs to be false if using the stub'
   iret=-1
endif
if (do_sppt_in.EQV. .true. ) then
   print*,'do_sppt needs to be false if using the stub'
   iret=-1
endif
return
end subroutine init_stochastic_physics_ocn

subroutine run_stochastic_physics_ocn(sppt_wts,t_rp1,t_rp2)
implicit none
real, intent(inout) :: sppt_wts(:,:) !< array containing random weights for SPPT range [0,2]
real, intent(inout) :: t_rp1(:,:)    !< array containing random weights for ePBL 
                                     !! perturbations (KE generation) range [0,2]
real, intent(inout) :: t_rp2(:,:)    !< array containing random weights for ePBL
                                     !! perturbations (KE dissipation) range [0,2]
return
end subroutine run_stochastic_physics_ocn

end module stochastic_physics

module get_stochy_pattern_mod

private

public  :: write_stoch_restart_ocn

contains
subroutine write_stoch_restart_ocn(sfile)

character(len=*) :: sfile   !< name of restart file
return
end subroutine write_stoch_restart_ocn

end module get_stochy_pattern_mod
