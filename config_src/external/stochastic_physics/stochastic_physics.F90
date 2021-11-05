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
real,intent(in)  :: delt
integer,intent(in) :: nx,ny,nz
real,intent(in) :: geoLonT(nx,ny),geoLatT(nx,ny)
logical,intent(in) :: pert_epbl_in,do_sppt_in
integer,intent(in)    :: mpiroot, mpicomm
integer, intent(out) :: iret

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
real, intent(inout) :: sppt_wts(:,:),t_rp1(:,:),t_rp2(:,:)
return
end subroutine run_stochastic_physics_ocn

end module stochastic_physics

module get_stochy_pattern_mod

private

public  :: write_stoch_restart_ocn

contains
subroutine write_stoch_restart_ocn(sfile)

character(len=*) :: sfile
return
end subroutine write_stoch_restart_ocn

end module get_stochy_pattern_mod
