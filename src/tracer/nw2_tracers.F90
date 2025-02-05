!> Ideal tracers designed to help diagnose a tracer diffusivity tensor in NeverWorld2
module nw2_tracers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_interface_heights, only : thickness_to_dz
use MOM_io, only : file_exists, MOM_read_data, slasher, vardesc, var_desc
use MOM_restart, only : query_initialized, set_initialized, MOM_restart_CS
use MOM_time_manager, only : time_type, time_type_to_real
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_nw2_tracers
public initialize_nw2_tracers
public nw2_tracer_column_physics
public nw2_tracers_end

!> The control structure for the nw2_tracers package
type, public :: nw2_tracers_CS ; private
  integer :: ntr = 0  !< The number of tracers that are actually used.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL()   !< The array of tracers used in this package, in [conc] (g m-3)?
  real, allocatable , dimension(:) :: restore_rate !< The rate at which the tracer is damped toward
                                             !! its target profile [T-1 ~> s-1]
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart controls structure
end type nw2_tracers_CS

contains

!> Register the NW2 tracer fields to be used with MOM.
logical function register_nw2_tracers(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(nw2_tracers_CS),       pointer    :: CS !< The control structure returned by a previous
                                               !! call to register_nw2_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "nw2_tracers" ! This module's name.
  character(len=8)  :: var_name ! The variable's name.
  real, pointer :: tr_ptr(:,:,:) => NULL() ! The tracer concentration [conc]
  integer :: isd, ied, jsd, jed, nz, m, ig
  integer :: n_groups ! Number of groups of three tracers (i.e. # tracers/3)
  real, allocatable, dimension(:) :: timescale_in_days ! Damping timescale [days]
  type(vardesc) :: tr_desc ! Descriptions and metadata for the tracers
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_nw2_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NW2_TRACER_GROUPS", n_groups, &
                 "The number of tracer groups where a group is of three tracers "//&
                 "initialized and restored to sin(x), y and z, respectively. Each "//&
                 "group is restored with an independent restoration rate.", &
                 default=3)
  allocate(timescale_in_days(n_groups))
  timescale_in_days = (/365., 730., 1460./)
  call get_param(param_file, mdl, "NW2_TRACER_RESTORE_TIMESCALE", timescale_in_days, &
                 "A list of timescales, one for each tracer group.", &
                 units="days")

  CS%ntr = 3 * n_groups
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)
  allocate(CS%restore_rate(CS%ntr))

  do m=1,CS%ntr
    write(var_name(1:8),'(a6,i2.2)') 'tracer',m
    tr_desc = var_desc(var_name, "1", "Ideal Tracer", caller=mdl)
    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, tr_desc=tr_desc, &
                         registry_diags=.true., restart_CS=restart_CS, mandatory=.false.)
    ig = int( (m+2)/3 ) ! maps (1,2,3)->1, (4,5,6)->2, ...
    CS%restore_rate(m) = 1.0 / ( timescale_in_days(ig) * 86400.0*US%s_to_T )
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_nw2_tracers = .true.
end function register_nw2_tracers

!> Sets the NW2 traces to their initial values and sets up the tracer output
subroutine initialize_nw2_tracers(restart, day, G, GV, US, h, tv, diag, CS)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various
                                                         !! thermodynamic variables
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(nw2_tracers_CS),               pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_nw2_tracer.
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta ! Interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz  ! Vertical extent of layers [Z ~> m]
  real :: rscl ! z* scaling factor [nondim]
  character(len=8)  :: var_name ! The variable's name.
  integer :: i, j, k, m

  if (.not.associated(CS)) return

  CS%Time => day
  CS%diag => diag

  ! Calculate z* interface positions
  call thickness_to_dz(h, tv, dz, G, GV, US)

  if (GV%Boussinesq) then
    ! First calculate interface positions in z-space (m)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,GV%ke+1) = - G%mask2dT(i,j) * G%bathyT(i,j)
    enddo ; enddo
    do k=GV%ke,1,-1 ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k)
    enddo ; enddo ; enddo
    ! Re-calculate for interface positions in z*-space (m)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%bathyT(i,j)>0.) then
        rscl = G%bathyT(i,j) / ( eta(i,j,1) + G%bathyT(i,j) )
        do K=GV%ke, 1, -1
          eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k) * rscl
        enddo
      endif
    enddo ; enddo
  else
    call MOM_error(FATAL, "NW2 tracers assume Boussinesq mode")
  endif

  do m=1,CS%ntr
    ! Initialize only if this is not a restart or we are using a restart
    ! in which the tracers were not present
    write(var_name(1:8),'(a6,i2.2)') 'tracer',m
    if ((.not.restart) .or. &
        (.not. query_initialized(CS%tr(:,:,:,m), var_name, CS%restart_CSp))) then
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
          CS%tr(i,j,k,m) = nw2_tracer_dist(m, G, GV, eta, i, j, k)
      enddo ; enddo ; enddo
      call set_initialized(CS%tr(:,:,:,m), var_name, CS%restart_CSp)
    endif ! restart
  enddo ! Tracer loop

end subroutine initialize_nw2_tracers

!> Applies diapycnal diffusion, aging and regeneration at the surface to the NW2 tracers
subroutine nw2_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, tv, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(nw2_tracers_CS),    pointer    :: CS   !< The control structure returned by a previous
                                              !! call to register_nw2_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]
!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.
! This is a simple example of a set of advected passive tracers.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)
  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: eta ! Interface heights [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: dz  ! Vertical extent of layers [Z ~> m]
  integer :: i, j, k, m
  real :: dt_x_rate ! dt * restoring rate [nondim]
  real :: rscl ! z* scaling factor [nondim]
  real :: target_value ! tracer target value for damping [conc]

! if (.not.associated(CS)) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

  ! Calculate z* interface positions
  call thickness_to_dz(h_new, tv, dz, G, GV, US)

  if (GV%Boussinesq) then
    ! First calculate interface positions in z-space [Z ~> m]
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,GV%ke+1) = - G%mask2dT(i,j) * G%bathyT(i,j)
    enddo ; enddo
    do k=GV%ke,1,-1 ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k)
    enddo ; enddo ; enddo
    ! Re-calculate for interface positions in z*-space [Z ~> m]
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (G%bathyT(i,j)>0.) then
        rscl = G%bathyT(i,j) / ( eta(i,j,1) + G%bathyT(i,j) )
        do K=GV%ke, 1, -1
          eta(i,j,K) = eta(i,j,K+1) + G%mask2dT(i,j) * dz(i,j,k) * rscl
        enddo
      endif
    enddo ; enddo
  else
    call MOM_error(FATAL, "NW2 tracers assume Boussinesq mode")
  endif

  do m=1,CS%ntr
    dt_x_rate = dt * CS%restore_rate(m)
    !$OMP parallel do default(shared) private(target_value)
    do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      target_value = nw2_tracer_dist(m, G, GV, eta, i, j, k)
      CS%tr(i,j,k,m) = CS%tr(i,j,k,m) + G%mask2dT(i,j) * dt_x_rate * ( target_value - CS%tr(i,j,k,m) )
    enddo ; enddo ; enddo
  enddo

end subroutine nw2_tracer_column_physics

!> The target value of a NeverWorld2 tracer label m [conc] at non-dimensional
!! position x=lon/Lx, y=lat/Ly, z=eta/H
real function nw2_tracer_dist(m, G, GV, eta, i, j, k)
  integer, intent(in) :: m !< Indicates the NW2 tracer
  type(ocean_grid_type),   intent(in) :: G   !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV  !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                           intent(in) :: eta !< Interface position [Z ~> m]
  integer, intent(in) :: i !< Cell index i
  integer, intent(in) :: j !< Cell index j
  integer, intent(in) :: k !< Layer index k
  ! Local variables
  real :: pi ! 3.1415... [nondim]
  real :: x, y, z ! non-dimensional relative positions [nondim]
  pi = 2.*acos(0.)
  x = ( G%geolonT(i,j) - G%west_lon ) / G%len_lon ! 0 ... 1
  y = -G%geolatT(i,j) / G%south_lat ! -1 ... 1
  z = - 0.5 * ( eta(i,j,K) + eta(i,j,K+1) ) / GV%max_depth ! 0 ... 1
  select case ( mod(m-1,3) )
  case (0) ! sin(2 pi x/L)
    nw2_tracer_dist = sin( 2.0 * pi * x )
  case (1) ! y/L
    nw2_tracer_dist = y
  case (2) ! -z/L
    nw2_tracer_dist = -z
  case default
    stop 'This should not happen. Died in nw2_tracer_dist()!'
  end select
  nw2_tracer_dist = nw2_tracer_dist * G%mask2dT(i,j)
end function nw2_tracer_dist

!> Deallocate any memory associated with this tracer package
subroutine nw2_tracers_end(CS)
  type(nw2_tracers_CS), pointer :: CS !< The control structure returned by a previous
                                      !! call to register_nw2_tracers.

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine nw2_tracers_end

!> \namespace nw2_tracers
!!
!!  TBD

end module nw2_tracers
