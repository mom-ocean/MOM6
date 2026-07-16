!>@brief The module 'stochy_data_mod' contains the initilization routine that read the stochastic phyiscs
!! namelist and determins the number of random patterns.
module stochy_data_mod

implicit none
public :: stoch_restfile

character(len=128) :: stoch_restfile = './INPUT/ocn_stoch.res.nc' !< default restart file name

end module stochy_data_mod
