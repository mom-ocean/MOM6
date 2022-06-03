! The are stubs for ocean stochastic physics
! the fully functional code is available at
! http://github.com/noaa-psd/stochastic_physics
module get_stochy_pattern_mod

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public  :: write_stoch_restart_ocn

contains

!> Write the restart file for the stochastic physics perturbations.
subroutine write_stoch_restart_ocn(sfile)
  character(len=*) :: sfile   !< name of restart file

  ! This stub function does not actually do anything.
  return
end subroutine write_stoch_restart_ocn

end module get_stochy_pattern_mod
