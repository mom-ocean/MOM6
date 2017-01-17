module MOM_regularize_layers
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
!*  By Robert Hallberg and Alistair Adcroft, 2011.                     *
!*                                                                     *
!*    This file contains the code to do vertical remapping of mass,    *
!*  temperature and salinity in MOM. Other tracers and the horizontal  *
!*  velocity components will be remapped outside of this subroutine    *
!*  using the values that are stored in ea and eb.                     *
!*    The code that is here now only applies in very limited cases     *
!*  where the mixed- and buffer-layer structures are problematic, but  *
!*  future additions will include the ability to emulate arbitrary     *
!*  vertical coordinates.                                              *
!*                                                                     *
!*  Macros written all in capital letters are defined in MOM_memory.h. *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, T, S, ea, eb, etc.                    *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_domains,       only : create_group_pass, do_group_pass
use MOM_domains,       only : group_pass_type
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>
#undef  DEBUG_CODE

public regularize_layers, regularize_layers_init

type, public :: regularize_layers_CS ; private
  logical :: regularize_surface_layers ! If true, vertically restructure the
                             ! near-surface layers when they have too much
                             ! lateral variations to allow for sensible lateral
                             ! barotropic transports.
  logical :: reg_sfc_detrain
  real    :: h_def_tol1      ! The value of the relative thickness deficit at
                             ! which to start modifying the structure, 0.5 by
                             ! default (or a thickness ratio of 5.83).
  real    :: h_def_tol2      ! The value of the relative thickness deficit at
                             ! which to the structure modification is in full
                             ! force, now 20% of the way from h_def_tol1 to 1.
  real    :: h_def_tol3      ! The values of the relative thickness defitic at
  real    :: h_def_tol4      ! which to start detrainment from the buffer layers
                             ! to the interior, and at which to do this at full
                             ! intensity.  Now 30% and 50% of the way from
                             ! h_def_tol1 to 1.
  real    :: Hmix_min        ! The minimum mixed layer thickness in m.
  type(time_type), pointer :: Time ! A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  logical :: debug           ! If true, do more thorough checks for debugging purposes.

  type(group_pass_type) :: pass_h ! For group pass

  integer :: id_def_rat = -1
  logical :: allow_clocks_in_omp_loops  ! If true, clocks can be called
                                        ! from inside loops that can be threaded.
                                        ! To run with multiple threads, set to False.
#ifdef DEBUG_CODE
  integer :: id_def_rat_2 = -1, id_def_rat_3 = -1
  integer :: id_def_rat_u = -1, id_def_rat_v = -1
  integer :: id_e1 = -1, id_e2 = -1, id_e3 = -1
  integer :: id_def_rat_u_1b = -1, id_def_rat_v_1b = -1
  integer :: id_def_rat_u_2 = -1, id_def_rat_u_2b = -1
  integer :: id_def_rat_v_2 = -1, id_def_rat_v_2b = -1
  integer :: id_def_rat_u_3 = -1, id_def_rat_u_3b = -1
  integer :: id_def_rat_v_3 = -1, id_def_rat_v_3b = -1
#endif
end type regularize_layers_CS

integer :: id_clock_pass, id_clock_EOS

contains

subroutine regularize_layers(h, tv, dt, ea, eb, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h
  type(thermo_var_ptrs),                    intent(inout) :: tv
  real,                                     intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: ea, eb
  type(regularize_layers_CS),               pointer       :: CS

!    This subroutine partially steps the bulk mixed layer model.
!  The following processes are executed, in the order listed.

! Arguments: h - Layer thickness, in m or kg m-2. (Intent in/out)
!                The units of h are referred to as H below.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - Time increment, in s.
!  (in/out)  ea - The amount of fluid moved downward into a layer; this should
!                 be increased due to mixed layer detrainment, in the same units
!                 as h - usually m or kg m-2 (i.e., H).
!  (in/out)  eb - The amount of fluid moved upward into a layer; this should
!                 be increased due to mixed layer entrainment, in the same units
!                 as h - usually m or kg m-2 (i.e., H).
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 regularize_layers_init.

  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_regularize_layers: "//&
         "Module must be initialized before it is used.")

  if (CS%regularize_surface_layers) then
    call cpu_clock_begin(id_clock_pass)
    call create_group_pass(CS%pass_h,h,G%Domain)
    call do_group_pass(CS%pass_h,G%Domain)
    call cpu_clock_end(id_clock_pass)
  endif

  if (CS%regularize_surface_layers) then
    call regularize_surface(h, tv, dt, ea, eb, G, GV, CS)
  endif

end subroutine regularize_layers

subroutine regularize_surface(h, tv, dt, ea, eb, G, GV, CS)
  type(ocean_grid_type),                    intent(inout) :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: h
  type(thermo_var_ptrs),                    intent(inout) :: tv
  real,                                     intent(in)    :: dt
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: ea, eb
  type(regularize_layers_CS),               pointer       :: CS

!    This subroutine ensures that there is a degree of horizontal smoothness
!  in the depths of the near-surface interfaces.

! Arguments: h - Layer thickness, in m or kg m-2. (Intent in/out)
!                The units of h are referred to as H below.
!  (in/out)  tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - Time increment, in s.
!  (in/out)  ea - The amount of fluid moved downward into a layer; this should
!                 be increased due to mixed layer detrainment, in the same units
!                 as h - usually m or kg m-2 (i.e., H).
!  (in/out)  eb - The amount of fluid moved upward into a layer; this should
!                 be increased due to mixed layer entrainment, in the same units
!                 as h - usually m or kg m-2 (i.e., H).
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 regularize_layers_init.

  real, dimension(SZIB_(G),SZJ_(G)) :: &
    def_rat_u   ! The ratio of the thickness deficit to the minimum depth, ND.
  real, dimension(SZI_(G),SZJB_(G)) :: &
    def_rat_v   ! The ratio of the thickness deficit to the minimum depth, ND.
  real, dimension(SZI_(G),SZJ_(G)) :: &
    def_rat_h   ! The ratio of the thickness deficit to the minimum depth, ND.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    e           ! The interface depths, in H, positive upward.

#ifdef DEBUG_CODE
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    def_rat_u_1b, def_rat_u_2, def_rat_u_2b, def_rat_u_3, def_rat_u_3b
  real, dimension(SZI_(G),SZJB_(G)) :: &
    def_rat_v_1b, def_rat_v_2, def_rat_v_2b, def_rat_v_3, def_rat_v_3b
  real, dimension(SZI_(G),SZJB_(G)) :: &
    def_rat_h2, def_rat_h3
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1) :: &
    ef          ! The filtered interface depths, in H, positive upward.
#endif

  real, dimension(SZI_(G),SZK_(G)+1) :: &
    e_filt, e_2d  ! The interface depths, in H, positive upward.
  real, dimension(SZI_(G),SZK_(G)) :: &
    h_2d, &     !   A 2-d version of h, in H.
    T_2d, &     !   A 2-d version of tv%T, in deg C.
    S_2d, &     !   A 2-d version of tv%S, in PSU.
    Rcv, &      !   A 2-d version of the coordinate density, in kg m-3.
    h_2d_init, &  ! The initial value of h_2d, in H.
    T_2d_init, &  ! THe initial value of T_2d, in deg C.
    S_2d_init, &  ! The initial value of S_2d, in PSU.
    d_eb, &     !   The downward increase across a layer in the entrainment from
                ! below, in H.  The sign convention is that positive values of
                ! d_eb correspond to a gain in mass by a layer by upward motion.
    d_ea        !   The upward increase across a layer in the entrainment from
                ! above, in H.  The sign convention is that positive values of
                ! d_ea mean a net gain in mass by a layer from downward motion.
  real, dimension(SZI_(G)) :: &
    p_ref_cv, & !   Reference pressure for the potential density which defines
                ! the coordinate variable, set to P_Ref, in Pa.
    Rcv_tol, &  !   A tolerence, relative to the target density differences
                ! between layers, for detraining into the interior, ND.
    h_add_tgt, h_add_tot, &
    h_tot1, Th_tot1, Sh_tot1, &
    h_tot3, Th_tot3, Sh_tot3, &
    h_tot2, Th_tot2, Sh_tot2
  real, dimension(SZK_(G)) :: &
    h_prev_1d     ! The previous thicknesses, in H.
  real :: I_dtol  ! The inverse of the tolerance changes, nondim.
  real :: I_dtol34 ! The inverse of the tolerance changes, nondim.
  real :: h1, h2  ! Temporary thicknesses, in H.
  real :: e_e, e_w, e_n, e_s  ! Temporary interface heights, in H.
  real :: wt    ! The weight of the filted interfaces in setting the targets, ND.
  real :: scale ! A scaling factor, ND.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  real, dimension(SZK_(G)+1) :: &
    int_flux, int_Tflux, int_Sflux, int_Rflux
  real :: h_add
  real :: h_det_tot
  real :: max_def_rat
  real :: Rcv_min_det  ! The lightest (min) and densest (max) coordinate density
  real :: Rcv_max_det  ! that can detrain into a layer, in kg m-3.

  real :: int_top, int_bot
  real :: h_predicted
  real :: h_prev
  real :: h_deficit

  logical :: cols_left, ent_any, more_ent_i(SZI_(G)), ent_i(SZI_(G))
  logical :: det_any, det_i(SZI_(G))
  logical :: do_j(SZJ_(G)), do_i(SZI_(G)), find_i(SZI_(G))
  logical :: debug = .false.
  logical :: fatal_error
  character(len=256) :: mesg    ! Message for error messages.
  integer :: i, j, k, is, ie, js, je, nz, nkmb, nkml, k1, k2, k3, ks, nz_filt

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_regularize_layers: "//&
         "Module must be initialized before it is used.")

  if (GV%nkml<1) return
  nkmb = GV%nk_rho_varies ; nkml = GV%nkml
  if (.not.ASSOCIATED(tv%eqn_of_state)) call MOM_error(FATAL, &
    "MOM_regularize_layers: This module now requires the use of temperature and "//&
    "an equation of state.")

  h_neglect = GV%H_subroundoff
  debug = (debug .or. CS%debug)
#ifdef DEBUG_CODE
  debug = .true.
  if (CS%id_def_rat_2 > 0) then ! Calculate over a slightly larger domain.
    is = G%isc-1 ; ie = G%iec+1 ; js = G%jsc-1 ; je = G%jec+1
  endif
#endif

  I_dtol = 1.0 / max(CS%h_def_tol2 - CS%h_def_tol1, 1e-40)
  I_dtol34 = 1.0 / max(CS%h_def_tol4 - CS%h_def_tol3, 1e-40)

  p_ref_cv(:) = tv%P_Ref

  do j=js-1,je+1 ; do i=is-1,ie+1
    e(i,j,1) = 0.0
  enddo ; enddo
  do K=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    e(i,j,K+1) = e(i,j,K) - h(i,j,k)
  enddo ; enddo ; enddo

#ifdef DEBUG_CODE
  call find_deficit_ratios(e, def_rat_u, def_rat_v, G, GV, CS, def_rat_u_1b, def_rat_v_1b, 1, h)
#else
  call find_deficit_ratios(e, def_rat_u, def_rat_v, G, GV, CS, h=h)
#endif
  ! Determine which columns are problematic
  do j=js,je ; do_j(j) = .false. ; enddo
  do j=js,je ; do i=is,ie
    def_rat_h(i,j) = max(def_rat_u(I-1,j), def_rat_u(I,j), &
                         def_rat_v(i,J-1), def_rat_v(i,J))
    if (def_rat_h(i,j) > CS%h_def_tol1) do_j(j) = .true.
  enddo ; enddo

#ifdef DEBUG_CODE
  if ((CS%id_def_rat_3 > 0) .or. (CS%id_e3 > 0) .or. &
      (CS%id_def_rat_u_3 > 0) .or. (CS%id_def_rat_u_3b > 0) .or. &
      (CS%id_def_rat_v_3 > 0) .or. (CS%id_def_rat_v_3b > 0) ) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      ef(i,j,1) = 0.0
    enddo ; enddo
    do K=2,nz+1 ; do j=js,je ; do i=is,ie
      if (G%mask2dCu(I,j) <= 0.0) then ; e_e = e(i,j,K) ; else
        e_e = max(e(i+1,j,K) + min(e(i,j,K) - e(i+1,j,nz+1), 0.0), e(i,j,nz+1))
      endif
      if (G%mask2dCu(I-1,j) <= 0.0) then ; e_w = e(i,j,K) ; else
        e_w = max(e(i-1,j,K) + min(e(i,j,K) - e(i-1,j,nz+1), 0.0), e(i,j,nz+1))
      endif
      if (G%mask2dCv(i,J) <= 0.0) then ; e_n = e(i,j,K) ; else
        e_n = max(e(i,j+1,K) + min(e(i,j,K) - e(i,j+1,nz+1), 0.0), e(i,j,nz+1))
      endif
      if (G%mask2dCv(i,J-1) <= 0.0) then ; e_s = e(i,j,K) ; else
        e_s = max(e(i,j-1,K) + min(e(i,j,K) - e(i,j-1,nz+1), 0.0), e(i,j,nz+1))
      endif

      wt = 1.0
      ef(i,j,k) = (1.0 - 0.5*wt) * e(i,j,K) + &
                  wt * 0.125 * ((e_e + e_w) + (e_n + e_s))
    enddo ; enddo ; enddo
    call find_deficit_ratios(ef, def_rat_u_3, def_rat_v_3, G, GV, CS, def_rat_u_3b, def_rat_v_3b)

    ! Determine which columns are problematic
    do j=js,je ; do i=is,ie
      def_rat_h3(i,j) = max(def_rat_u_3(I-1,j), def_rat_u_3(I,j), &
                            def_rat_v_3(i,J-1), def_rat_v_3(i,J))
    enddo ; enddo

    if (CS%id_e3 > 0) call post_data(CS%id_e3, ef, CS%diag)
    if (CS%id_def_rat_3 > 0) call post_data(CS%id_def_rat_3, def_rat_h3, CS%diag)
    if (CS%id_def_rat_u_3 > 0) call post_data(CS%id_def_rat_u_3, def_rat_u_3, CS%diag)
    if (CS%id_def_rat_u_3b > 0) call post_data(CS%id_def_rat_u_3b, def_rat_u_3b, CS%diag)
    if (CS%id_def_rat_v_3 > 0) call post_data(CS%id_def_rat_v_3, def_rat_v_3, CS%diag)
    if (CS%id_def_rat_v_3b > 0) call post_data(CS%id_def_rat_v_3b, def_rat_v_3b, CS%diag)
  endif
#endif


  ! Now restructure the layers.
!$OMP parallel do default(none) shared(is,ie,js,je,nz,do_j,def_rat_h,CS,nkmb,G,GV,&
!$OMP                                  e,I_dtol,h,tv,debug,h_neglect,p_ref_cv,ea, &
!$OMP                                  eb,id_clock_EOS,nkml)                      &
!$OMP                          private(d_ea,d_eb,max_def_rat,do_i,nz_filt,e_e,e_w,&
!$OMP                                  e_n,e_s,wt,e_filt,e_2d,h_2d,T_2d,S_2d,     &
!$OMP                                  h_2d_init,T_2d_init,S_2d_init,ent_any,     &
!$OMP                                  more_ent_i,ent_i,h_add_tgt,h_add_tot,      &
!$OMP                                  cols_left,h_add,h_prev,ks,det_any,det_i,   &
!$OMP                                  Rcv_tol,Rcv,k1,k2,h_det_tot,Rcv_min_det,   &
!$OMP                                  Rcv_max_det,h_deficit,h_tot3,Th_tot3,      &
!$OMP                                  Sh_tot3,scale,int_top,int_flux,int_Rflux,  &
!$OMP                                  int_Tflux,int_Sflux,int_bot,h_prev_1d,     &
!$OMP                                  h_tot1,Th_tot1,Sh_tot1,h_tot2,Th_tot2,     &
!$OMP                                  Sh_tot2,h_predicted,fatal_error,mesg )
  do j=js,je ; if (do_j(j)) then

!  call cpu_clock_begin(id_clock_EOS)
!  call calculate_density_derivs(T(:,1), S(:,1), p_ref_cv, dRcv_dT, dRcv_dS, &
!                                is, ie-is+1, tv%eqn_of_state)
!  call cpu_clock_end(id_clock_EOS)

    do k=1,nz ; do i=is,ie ; d_ea(i,k) = 0.0 ; d_eb(i,k) = 0.0 ; enddo ; enddo

    max_def_rat = 0.0
    do i=is,ie
      do_i(i) = def_rat_h(i,j) > CS%h_def_tol1
      if (def_rat_h(i,j) > max_def_rat) max_def_rat = def_rat_h(i,j)
    enddo
    nz_filt = nkmb+1 ; if (max_def_rat > CS%h_def_tol3) nz_filt = nz+1

    ! Find a 2-D 1-2-1 filtered version of e to target.  Area weights are
    ! deliberately omitted here.  This is slightly more complicated than a
    ! simple filter so that the effects of topography are eliminated.
    do K=1,nz_filt ; do i=is,ie ; if (do_i(i)) then
      if (G%mask2dCu(I,j) <= 0.0) then ; e_e = e(i,j,K) ; else
        e_e = max(e(i+1,j,K) + min(e(i,j,K) - e(i+1,j,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom)

      endif
      if (G%mask2dCu(I-1,j) <= 0.0) then ; e_w = e(i,j,K) ; else
        e_w = max(e(i-1,j,K) + min(e(i,j,K) - e(i-1,j,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom)
      endif
      if (G%mask2dCv(i,J) <= 0.0) then ; e_n = e(i,j,K) ; else
        e_n = max(e(i,j+1,K) + min(e(i,j,K) - e(i,j+1,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom)
      endif
      if (G%mask2dCv(i,J-1) <= 0.0) then ; e_s = e(i,j,K) ; else
        e_s = max(e(i,j-1,K) + min(e(i,j,K) - e(i,j-1,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom)
      endif

      wt = max(0.0, min(1.0, I_dtol*(def_rat_h(i,j)-CS%h_def_tol1)))

      e_filt(i,k) = (1.0 - 0.5*wt) * e(i,j,K) + &
                  wt * 0.125 * ((e_e + e_w) + (e_n + e_s))
      e_2d(i,k) = e(i,j,K)
    endif ; enddo ; enddo
    do k=1,nz ; do i=is,ie
      h_2d(i,k) = h(i,j,k)
      T_2d(i,k) = tv%T(i,j,k) ; S_2d(i,k) = tv%S(i,j,k)
    enddo ; enddo

    if (debug) then
      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        h_2d_init(i,k) = h(i,j,k)
        T_2d_init(i,k) = tv%T(i,j,k) ; S_2d_init(i,k) = tv%S(i,j,k)
      endif ; enddo ; enddo
    endif

    ! First, try to entrain from the interior.
    ent_any = .false.
    do i=is,ie
      more_ent_i(i) = .false. ; ent_i(i) = .false.
      h_add_tgt(i) = 0.0 ; h_add_tot(i) = 0.0
      if (do_i(i) .and. (e_2d(i,nkmb+1) > e_filt(i,nkmb+1))) then
        more_ent_i(i) = .true. ; ent_i(i) = .true. ; ent_any = .true.
        h_add_tgt(i) = e_2d(i,nkmb+1) - e_filt(i,nkmb+1)
      endif
    enddo

    if (ent_any) then
      do k=nkmb+1,nz
        cols_left = .false.
        do i=is,ie ; if (more_ent_i(i)) then
          if (h_2d(i,k) - GV%Angstrom > h_neglect) then
            if (e_2d(i,nkmb+1)-e_filt(i,nkmb+1) > h_2d(i,k) - GV%Angstrom) then
              h_add = h_2d(i,k) - GV%Angstrom
              h_2d(i,k) = GV%Angstrom
            else
              h_add = e_2d(i,nkmb+1)-e_filt(i,nkmb+1)
              h_2d(i,k) = h_2d(i,k) - h_add
            endif
            d_eb(i,k-1) = d_eb(i,k-1) + h_add
            h_add_tot(i) = h_add_tot(i) + h_add
            e_2d(i,nkmb+1) = e_2d(i,nkmb+1) - h_add
            h_prev = h_2d(i,nkmb)
            h_2d(i,nkmb) = h_2d(i,nkmb) + h_add

            T_2d(i,nkmb) = (h_prev*T_2d(i,nkmb) + h_add*T_2d(i,k)) / h_2d(i,nkmb)
            S_2d(i,nkmb) = (h_prev*S_2d(i,nkmb) + h_add*S_2d(i,k)) / h_2d(i,nkmb)

            if ((e_2d(i,nkmb+1) <= e_filt(i,nkmb+1)) .or. &
                (h_add_tot(i) > 0.6*h_add_tgt(i))) then  !### 0.6 is adjustable?.
              more_ent_i(i) = .false.
            else
              cols_left = .true.
            endif
          else
            cols_left = .true.
          endif
        endif ; enddo
        if (.not.cols_left) exit
      enddo

      ks = min(k-1,nz-1)
      do k=ks,nkmb,-1 ; do i=is,ie ; if (ent_i(i)) then
        d_eb(i,k) = d_eb(i,k) + d_eb(i,k+1)
      endif ; enddo ; enddo
    endif ! ent_any

    !   This is where code to detrain to the interior will go.
    ! The buffer layers can only detrain water into layers when the buffer
    ! layer potential density is between (c*Rlay(k-1) + (1-c)*Rlay(k)) and
    ! (c*Rlay(k+1) + (1-c)*Rlay(k)), where 0.5 <= c < 1.0.
    !    Do not detrain if the 2-layer deficit ratio is not significant.
    !    Detrainment must be able to come from all mixed and buffer layers.
    !    All water is moved out of the buffer layers below before moving from
    !  a shallower layer (characteristics do not cross).
    det_any = .false.
    if ((max_def_rat > CS%h_def_tol3) .and. (CS%reg_sfc_detrain)) then
      do i=is,ie
        det_i(i) = .false. ; Rcv_tol(i) = 0.0
        if (do_i(i) .and. (e_2d(i,nkmb+1) < e_filt(i,nkmb+1)) .and. &
            (def_rat_h(i,j) > CS%h_def_tol3)) then
          det_i(i) = .true. ; det_any = .true.
          Rcv_tol(i) = min((def_rat_h(i,j) - CS%h_def_tol3), 1.0)
        endif
      enddo
    endif
    if (det_any) then
      call cpu_clock_begin(id_clock_EOS)
      do k=1,nkmb
        call calculate_density(T_2d(:,k),S_2d(:,k),p_ref_cv,Rcv(:,k), &
                               is,ie-is+1,tv%eqn_of_state)
      enddo
      call cpu_clock_end(id_clock_EOS)

      do i=is,ie ; if (det_i(i)) then
        k1 = nkmb ; k2 = nz
        h_det_tot = 0.0
        do ! This loop is terminated by exits.
          if (k1 <= 1) exit
          if (k2 <= nkmb) exit
          ! ### The 0.6 here should be adjustable?  It gives 20% overlap for now.
          Rcv_min_det = GV%Rlay(k2) + 0.6*Rcv_tol(i)*(GV%Rlay(k2-1)-GV%Rlay(k2))
          if (k2 < nz) then
            Rcv_max_det = GV%Rlay(k2) + 0.6*Rcv_tol(i)*(GV%Rlay(k2+1)-GV%Rlay(k2))
          else
            Rcv_max_det = GV%Rlay(nz) + 0.6*Rcv_tol(i)*(GV%Rlay(nz)-GV%Rlay(nz-1))
          endif
          if (Rcv(i,k1) > Rcv_max_det) &
            exit ! All shallower interior layers are too light for detrainment.

          h_deficit = (e_filt(i,k2)-e_filt(i,k2+1)) - h_2d(i,k2)
          if ((e_filt(i,k2) > e_2d(i,k1+1)) .and. (h_deficit > 0.0) .and. &
              (Rcv(i,k1) < Rcv_max_det) .and. (Rcv(i,k1) > Rcv_min_det)) then
            ! Detrainment will occur.
            h_add = min(e_filt(i,k2) - e_2d(i,k2), h_deficit )
            if (h_add < h_2d(i,k1)) then
              ! Only part of layer k1 detrains.
              if (h_add > 0.0) then
                h_prev = h_2d(i,k2)
                h_2d(i,k2) = h_2d(i,k2) + h_add
                e_2d(i,k2) = e_2d(i,k2+1) + h_2d(i,k2)
                d_ea(i,k2) = d_ea(i,k2) + h_add
                ! ### THIS IS UPWIND.  IT SHOULD BE HIGHER ORDER...
                T_2d(i,k2) = (h_prev*T_2d(i,k2) + h_add*T_2d(i,k1)) / h_2d(i,k2)
                S_2d(i,k2) = (h_prev*S_2d(i,k2) + h_add*S_2d(i,k1)) / h_2d(i,k2)
                h_det_tot = h_det_tot + h_add

                h_2d(i,k1) = h_2d(i,k1) - h_add
                do k3=k1,nkmb ; e_2d(i,k3+1) = e_2d(i,k3) - h_2d(i,k3) ; enddo
                do k3=k1+1,nkmb ; d_ea(i,k3) = d_ea(i,k3) + h_add ; enddo
              else
                if (h_add < 0.0) &
                  call MOM_error(FATAL, "h_add is negative.  Some logic is wrong.")
                h_add = 0.0 ! This usually should not happen...
              endif

              ! Move up to the next target layer.
              k2 = k2-1
              if (k2>nkmb+1) e_2d(i,k2) = e_2d(i,k2) + h_det_tot
            else
              h_add = h_2d(i,k1)
              h_prev = h_2d(i,k2)
              h_2d(i,k2) = h_2d(i,k2) + h_add
              e_2d(i,k2) = e_2d(i,k2+1) + h_2d(i,k2)
              d_ea(i,k2) = d_ea(i,k2) + h_add
              T_2d(i,k2) = (h_prev*T_2d(i,k2) + h_add*T_2d(i,k1)) / h_2d(i,k2)
              S_2d(i,k2) = (h_prev*S_2d(i,k2) + h_add*S_2d(i,k1)) / h_2d(i,k2)
              h_det_tot = h_det_tot + h_add

              h_2d(i,k1) = 0.0
              do k3=k1,nkmb ; e_2d(i,k3+1) = e_2d(i,k3) - h_2d(i,k3) ; enddo
              do k3=k1+1,nkmb ; d_ea(i,k3) = d_ea(i,k3) + h_add ; enddo

              ! Move up to the next source layer.
              k1 = k1-1
            endif

          else
            ! Move up to the next target layer.
            k2 = k2-1
            if (k2>nkmb+1) e_2d(i,k2) = e_2d(i,k2) + h_det_tot
          endif

        enddo ! exit terminated loop.
      endif ; enddo
      ! ### This could be faster if the deepest k with nonzero d_ea were kept.
      do k=nz-1,nkmb+1,-1 ; do i=is,ie ; if (det_i(i)) then
        d_ea(i,k) = d_ea(i,k) + d_ea(i,k+1)
      endif ; enddo ; enddo
    endif  ! Detrainment to the interior.
    if (debug) then
      do i=is,ie ; h_tot3(i) = 0.0 ; Th_tot3(i) = 0.0 ; Sh_tot3(i) = 0.0 ; enddo
      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        h_tot3(i) = h_tot3(i) + h_2d(i,k)
        Th_tot3(i) = Th_tot3(i) + h_2d(i,k) * T_2d(i,k)
        Sh_tot3(i) = Sh_tot3(i) + h_2d(i,k) * S_2d(i,k)
      endif ; enddo ; enddo
    endif

    do i=is,ie ; if (do_i(i)) then
      ! Rescale the interface targets so the depth at the bottom of the deepest
      ! buffer layer matches.
      scale = e_2d(i,nkmb+1) / e_filt(i,nkmb+1)
      do k=2,nkmb+1 ; e_filt(i,k) = e_filt(i,k) * scale ; enddo

      ! Ensure that layer 1 only has water from layers 1 to nkml and rescale
      ! the remaining layer thicknesses if necessary.
      if (e_filt(i,2) < e_2d(i,nkml)) then
        scale = (e_2d(i,nkml) - e_filt(i,nkmb+1)) / &
                ((e_filt(i,2) - e_filt(i,nkmb+1)) + h_neglect)
        do k=3,nkmb
          e_filt(i,k) = e_filt(i,nkmb+1) + scale * (e_filt(i,k) - e_filt(i,nkmb+1))
        enddo
        e_filt(i,2) = e_2d(i,nkml)
      endif

      ! Map the water back into the layers.
      k1 = 1 ; k2 = 1
      int_top = 0.0
      do k=1,nkmb+1
        int_flux(k) = 0.0 ; int_Rflux(k) = 0.0
        int_Tflux(k) = 0.0 ; int_Sflux(k) = 0.0
      enddo
      do k=1,2*nkmb
        int_bot = max(e_2d(i,k1+1),e_filt(i,k2+1))
        h_add = int_top - int_bot

        if (k2 > k1) then
          do k3=k1+1,k2
            d_ea(i,k3) = d_ea(i,k3) + h_add
            int_flux(k3) = int_flux(k3) + h_add
            int_Tflux(k3) = int_Tflux(k3) + h_add*T_2d(i,k1)
            int_Sflux(k3) = int_Sflux(k3) + h_add*S_2d(i,k1)
          enddo
        elseif (k1 > k2) then
          do k3=k2,k1-1
            d_eb(i,k3) = d_eb(i,k3) + h_add
            int_flux(k3+1) = int_flux(k3+1) - h_add
            int_Tflux(k3+1) = int_Tflux(k3+1) - h_add*T_2d(i,k1)
            int_Sflux(k3+1) = int_Sflux(k3+1) - h_add*S_2d(i,k1)
          enddo
        endif

        if (int_bot <= e_filt(i,k2+1)) then
          ! Increment the target layer.
          k2 = k2 + 1
        elseif (int_bot <= e_2d(i,k1+1)) then
          ! Increment the source layer.
          k1 = k1 + 1
        else
          call MOM_error(FATAL, &
            "Regularize_surface: Could not increment target or source.")
        endif
        if ((k1 > nkmb) .or. (k2 > nkmb)) exit
        int_top = int_bot
      enddo
      if (k2 < nkmb) &
        call MOM_error(FATAL, "Regularize_surface: Did not assign fluid to layer nkmb.")

      ! Note that movement of water across the base of the bottommost buffer
      ! layer has already been dealt with separately.
      do k=1,nkmb ; h_prev_1d(k) = h_2d(i,k) ; enddo
      h_2d(i,1) = h_2d(i,1) - int_flux(2)
      do k=2,nkmb-1
        h_2d(i,k) = h_2d(i,k) + (int_flux(k) - int_flux(k+1))
      enddo
      ! Note that movement of water across the base of the bottommost buffer
      ! layer has already been dealt with separately.
      h_2d(i,nkmb) = h_2d(i,nkmb) + int_flux(nkmb)

      T_2d(i,1) = (T_2d(i,1)*h_prev_1d(1) - int_Tflux(2)) / h_2d(i,1)
      S_2d(i,1) = (S_2d(i,1)*h_prev_1d(1) - int_Sflux(2)) / h_2d(i,1)
      do k=2,nkmb-1
        T_2d(i,k) = (T_2d(i,k)*h_prev_1d(k) + (int_Tflux(k) - int_Tflux(k+1))) / h_2d(i,k)
        S_2d(i,k) = (S_2d(i,k)*h_prev_1d(k) + (int_Sflux(k) - int_Sflux(k+1))) / h_2d(i,k)
      enddo
      T_2d(i,nkmb) = (T_2d(i,nkmb)*h_prev_1d(nkmb) + int_Tflux(nkmb) ) / h_2d(i,nkmb)
      S_2d(i,nkmb) = (S_2d(i,nkmb)*h_prev_1d(nkmb) + int_Sflux(nkmb) ) / h_2d(i,nkmb)

    endif ; enddo ! i-loop

    ! Copy the interior thicknesses and other fields back to the 3-d arrays.
    do k=1,nz ; do i=is,ie ; if (do_i(i)) then
      h(i,j,k) = h_2d(i,k)
      tv%T(i,j,k) = T_2d(i,k) ; tv%S(i,j,k) = S_2d(i,k)
      ea(i,j,k) = ea(i,j,k) + d_ea(i,k)
      eb(i,j,k) = eb(i,j,k) + d_eb(i,k)
    endif ; enddo ; enddo

    if (debug) then
      do i=is,ie ; h_tot1(i) = 0.0 ; Th_tot1(i) = 0.0 ; Sh_tot1(i) = 0.0 ; enddo
      do i=is,ie ; h_tot2(i) = 0.0 ; Th_tot2(i) = 0.0 ; Sh_tot2(i) = 0.0 ; enddo

      do k=1,nz ; do i=is,ie ; if (do_i(i)) then
        h_tot1(i) = h_tot1(i) + h_2d_init(i,k)
        h_tot2(i) = h_tot2(i) + h(i,j,k)

        Th_tot1(i) = Th_tot1(i) + h_2d_init(i,k) * T_2d_init(i,k)
        Th_tot2(i) = Th_tot2(i) + h(i,j,k) * tv%T(i,j,k)
        Sh_tot1(i) = Sh_tot1(i) + h_2d_init(i,k) * S_2d_init(i,k)
        Sh_tot2(i) = Sh_tot2(i) + h(i,j,k) * tv%S(i,j,k)
        if (h(i,j,k) < 0.0) &
          call MOM_error(FATAL,"regularize_surface: Negative thicknesses.")
        if (k==1) then ; h_predicted = h_2d_init(i,k) + (d_eb(i,k) - d_ea(i,k+1))
        elseif (k==nz) then ; h_predicted = h_2d_init(i,k) + (d_ea(i,k) - d_eb(i,k-1))
        else
          h_predicted = h_2d_init(i,k) + ((d_ea(i,k) - d_eb(i,k-1)) + &
                                          (d_eb(i,k) - d_ea(i,k+1)))
        endif
        if (abs(h(i,j,k) - h_predicted) > MAX(1e-9*abs(h_predicted),GV%Angstrom)) &
          call MOM_error(FATAL, "regularize_surface: d_ea mismatch.")
      endif ; enddo ; enddo
      do i=is,ie ; if (do_i(i)) then
        fatal_error = .false.
        if (abs(h_tot1(i) - h_tot2(i)) > 1e-12*h_tot1(i)) then
          write(mesg,'(ES11.4," became ",ES11.4," diff ",ES11.4)') &
                h_tot1(i), h_tot2(i), (h_tot1(i) - h_tot2(i))
          call MOM_error(WARNING, "regularize_surface: Mass non-conservation."//&
                          trim(mesg), .true.)
          fatal_error = .true.
        endif
        if (abs(Th_tot1(i) - Th_tot2(i)) > 1e-12*(Th_tot1(i)+10.0*h_tot1(i))) then
          write(mesg,'(ES11.4," became ",ES11.4," diff ",ES11.4," int diff ",ES11.4)') &
                Th_tot1(i), Th_tot2(i), (Th_tot1(i) - Th_tot2(i)), (Th_tot1(i) - Th_tot3(i))
          call MOM_error(WARNING, "regularize_surface: Heat non-conservation."//&
                          trim(mesg), .true.)
          fatal_error = .true.
        endif
        if (abs(Sh_tot1(i) - Sh_tot2(i)) > 1e-12*(Sh_tot1(i)+10.0*h_tot1(i))) then
          write(mesg,'(ES11.4," became ",ES11.4," diff ",ES11.4," int diff ",ES11.4)') &
                Sh_tot1(i), Sh_tot2(i), (Sh_tot1(i) - Sh_tot2(i)), (Sh_tot1(i) - Sh_tot3(i))
          call MOM_error(WARNING, "regularize_surface: Salinity non-conservation."//&
                          trim(mesg), .true.)
          fatal_error = .true.
        endif
        if (fatal_error) then
          write(mesg,'("Error at lat/lon ",2(ES11.4))') G%geoLatT(i,j), G%geoLonT(i,j)
          call MOM_error(FATAL, "regularize_surface: Terminating with fatal error.  "//&
                          trim(mesg))
        endif
      endif ; enddo
    endif

  endif ; enddo ! j-loop.

  if (CS%id_def_rat > 0) call post_data(CS%id_def_rat, def_rat_h, CS%diag)

#ifdef DEBUG_CODE
  if (CS%id_e1 > 0) call post_data(CS%id_e1, e, CS%diag)
  if (CS%id_def_rat_u > 0) call post_data(CS%id_def_rat_u, def_rat_u, CS%diag)
  if (CS%id_def_rat_u_1b > 0) call post_data(CS%id_def_rat_u_1b, def_rat_u_1b, CS%diag)
  if (CS%id_def_rat_v > 0) call post_data(CS%id_def_rat_v, def_rat_v, CS%diag)
  if (CS%id_def_rat_v_1b > 0) call post_data(CS%id_def_rat_v_1b, def_rat_v_1b, CS%diag)

  if ((CS%id_def_rat_2 > 0) .or. &
      (CS%id_def_rat_u_2 > 0) .or. (CS%id_def_rat_u_2b > 0) .or. &
      (CS%id_def_rat_v_2 > 0) .or. (CS%id_def_rat_v_2b > 0) ) then
    do j=js-1,je+1 ; do i=is-1,ie+1
      e(i,j,1) = 0.0
    enddo ; enddo
    do K=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
      e(i,j,K+1) = e(i,j,K) - h(i,j,k)
    enddo ; enddo ; enddo

    call find_deficit_ratios(e, def_rat_u_2, def_rat_v_2, G, GV, CS, def_rat_u_2b, def_rat_v_2b, h=h)

    ! Determine which columns are problematic
    do j=js,je ; do i=is,ie
      def_rat_h2(i,j) = max(def_rat_u_2(I-1,j), def_rat_u_2(I,j), &
                            def_rat_v_2(i,J-1), def_rat_v_2(i,J))
    enddo ; enddo

    if (CS%id_def_rat_2 > 0) call post_data(CS%id_def_rat_2, def_rat_h2, CS%diag)
    if (CS%id_e2 > 0) call post_data(CS%id_e2, e, CS%diag)
    if (CS%id_def_rat_u_2 > 0) call post_data(CS%id_def_rat_u_2, def_rat_u_2, CS%diag)
    if (CS%id_def_rat_u_2b > 0) call post_data(CS%id_def_rat_u_2b, def_rat_u_2b, CS%diag)
    if (CS%id_def_rat_v_2 > 0) call post_data(CS%id_def_rat_v_2, def_rat_v_2, CS%diag)
    if (CS%id_def_rat_v_2b > 0) call post_data(CS%id_def_rat_v_2b, def_rat_v_2b, CS%diag)
  endif
#endif

end subroutine regularize_surface

subroutine find_deficit_ratios(e, def_rat_u, def_rat_v, G, GV, CS, &
                               def_rat_u_2lay, def_rat_v_2lay, halo, h)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(in) :: e
  real, dimension(SZIB_(G),SZJ_(G)),         intent(out) :: def_rat_u
  real, dimension(SZI_(G),SZJB_(G)),         intent(out) :: def_rat_v
  type(regularize_layers_CS),                pointer     :: CS
  real, dimension(SZIB_(G),SZJ_(G)), optional, intent(out) :: def_rat_u_2lay
  real, dimension(SZI_(G),SZJB_(G)), optional, intent(out) :: def_rat_v_2lay
  integer,                         optional, intent(in)  :: halo
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), optional, intent(in)  :: h
!    This subroutine determines the amount by which the harmonic mean
!  thickness at velocity points differ from the arithmetic means, relative to
!  the the arithmetic means, after eliminating thickness variations that are
!  solely due to topography and aggregating all interior layers into one.

! Arguments: e - Interface depths, in m or kg m-2.
!  (out)     def_rat_u - The thickness deficit ratio at u points, nondim.
!  (out)     def_rat_v - The thickness deficit ratio at v points, nondim.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 regularize_layers_init.
!  (out,opt) def_rat_u_2lay - The thickness deficit ratio at u points when the
!                 mixed and buffer layers are aggregated into 1 layer, nondim.
!  (out,opt) def_rat_v_2lay - The thickness deficit ratio at v pointswhen the
!                 mixed and buffer layers are aggregated into 1 layer, nondim.
!  (in,opt)  halo - An extra-wide halo size, 0 by default.
!  (in,opt)  h - The layer thicknesse; if not present take vertical differences of e.
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    h_def_u, &  ! The vertically summed thickness deficits at u-points, in H.
    h_norm_u, & ! The vertically summed arithmetic mean thickness by which
                ! h_def_u is normalized, in H.
    h_def2_u
  real, dimension(SZI_(G),SZJB_(G)) :: &
    h_def_v, &  ! The vertically summed thickness deficits at v-points, in H.
    h_norm_v, & ! The vertically summed arithmetic mean thickness by which
                ! h_def_v is normalized, in H.
    h_def2_v
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  real :: h1, h2  ! Temporary thicknesses, in H.
  integer :: i, j, k, is, ie, js, je, nz, nkmb

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif
  nkmb = GV%nk_rho_varies
  h_neglect = GV%H_subroundoff

  ! Determine which zonal faces are problematic.
  do j=js,je ; do I=is-1,ie
    ! Aggregate all water below the mixed and buffer layers for the purposes of
    ! this diagnostic.
    h1 = e(i,j,nkmb+1)-e(i,j,nz+1) ; h2 = e(i+1,j,nkmb+1)-e(i+1,j,nz+1)
    if (e(i,j,nz+1) < e(i+1,j,nz+1)) then
      if (h1 > h2) h1 = max(e(i,j,nkmb+1)-e(i+1,j,nz+1), h2)
    elseif (e(i+1,j,nz+1) < e(i,j,nz+1)) then
      if (h2 > h1) h2 = max(e(i+1,j,nkmb+1)-e(i,j,nz+1), h1)
    endif
    h_def_u(I,j) = 0.5*(h1-h2)**2 / ((h1 + h2) + h_neglect)
    h_norm_u(I,j) = 0.5*(h1+h2)
  enddo ; enddo
  if (present(def_rat_u_2lay)) then ; do j=js,je ; do I=is-1,ie
    ! This is a particular metric of the aggregation into two layers.
    h1 = e(i,j,1)-e(i,j,nkmb+1) ; h2 = e(i+1,j,1)-e(i+1,j,nkmb+1)
    if (e(i,j,nkmb+1) < e(i+1,j,nz+1)) then
      if (h1 > h2) h1 = max(e(i,j,1)-e(i+1,j,nz+1), h2)
    elseif (e(i+1,j,nkmb+1) < e(i,j,nz+1)) then
      if (h2 > h1) h2 = max(e(i+1,j,1)-e(i,j,nz+1), h1)
    endif
    h_def2_u(I,j) = h_def_u(I,j) + 0.5*(h1-h2)**2 / ((h1 + h2) + h_neglect)
  enddo ; enddo ; endif
  do k=1,nkmb ; do j=js,je ; do I=is-1,ie
    if (present(h)) then
      h1 = h(i,j,k) ; h2 = h(i+1,j,k)
    else
      h1 = e(i,j,K)-e(i,j,K+1) ; h2 = e(i+1,j,K)-e(i+1,j,K+1)
    endif
    ! Thickness deficits can not arise simply because a layer's bottom is bounded
    ! by the bathymetry.
    if (e(i,j,K+1) < e(i+1,j,nz+1)) then
      if (h1 > h2) h1 = max(e(i,j,K)-e(i+1,j,nz+1), h2)
    elseif (e(i+1,j,K+1) < e(i,j,nz+1)) then
      if (h2 > h1) h2 = max(e(i+1,j,K)-e(i,j,nz+1), h1)
    endif
    h_def_u(I,j) = h_def_u(I,j) + 0.5*(h1-h2)**2 / ((h1 + h2) + h_neglect)
    h_norm_u(I,j) = h_norm_u(I,j) + 0.5*(h1+h2)
  enddo ; enddo ; enddo
  if (present(def_rat_u_2lay)) then ; do j=js,je ; do I=is-1,ie
    def_rat_u(I,j) = G%mask2dCu(I,j) * h_def_u(I,j) / &
                     (max(CS%Hmix_min, h_norm_u(I,j)) + h_neglect)
    def_rat_u_2lay(I,j) = G%mask2dCu(I,j) * h_def2_u(I,j) / &
                          (max(CS%Hmix_min, h_norm_u(I,j)) + h_neglect)
  enddo ; enddo ; else ; do j=js,je ; do I=is-1,ie
    def_rat_u(I,j) = G%mask2dCu(I,j) * h_def_u(I,j) / &
                     (max(CS%Hmix_min, h_norm_u(I,j)) + h_neglect)
  enddo ; enddo ; endif

  ! Determine which meridional faces are problematic.
  do J=js-1,je ; do i=is,ie
    ! Aggregate all water below the mixed and buffer layers for the purposes of
    ! this diagnostic.
    h1 = e(i,j,nkmb+1)-e(i,j,nz+1) ; h2 = e(i,j+1,nkmb+1)-e(i,j+1,nz+1)
    if (e(i,j,nz+1) < e(i,j+1,nz+1)) then
      if (h1 > h2) h1 = max(e(i,j,nkmb+1)-e(i,j+1,nz+1), h2)
    elseif (e(i,j+1,nz+1) < e(i,j,nz+1)) then
      if (h2 > h1) h2 = max(e(i,j+1,nkmb+1)-e(i,j,nz+1), h1)
    endif
    h_def_v(i,J) = 0.5*(h1-h2)**2 / ((h1 + h2) + h_neglect)
    h_norm_v(i,J) = 0.5*(h1+h2)
  enddo ; enddo
  if (present(def_rat_v_2lay)) then ; do J=js-1,je ; do i=is,ie
    ! This is a particular metric of the aggregation into two layers.
    h1 = e(i,j,1)-e(i,j,nkmb+1) ; h2 = e(i,j+1,1)-e(i,j+1,nkmb+1)
    if (e(i,j,nkmb+1) < e(i,j+1,nz+1)) then
      if (h1 > h2) h1 = max(e(i,j,1)-e(i,j+1,nz+1), h2)
    elseif (e(i,j+1,nkmb+1) < e(i,j,nz+1)) then
      if (h2 > h1) h2 = max(e(i,j+1,1)-e(i,j,nz+1), h1)
    endif
    h_def2_v(i,J) = h_def_v(i,J) + 0.5*(h1-h2)**2 / ((h1 + h2) + h_neglect)
  enddo ; enddo ; endif
  do k=1,nkmb ; do J=js-1,je ; do i=is,ie
    if (present(h)) then
      h1 = h(i,j,k) ; h2 = h(i,j+1,k)
    else
      h1 = e(i,j,K)-e(i,j,K+1) ; h2 = e(i,j+1,K)-e(i,j+1,K+1)
    endif
    ! Thickness deficits can not arise simply because a layer's bottom is bounded
    ! by the bathymetry.
    if (e(i,j,K+1) < e(i,j+1,nz+1)) then
      if (h1 > h2) h1 = max(e(i,j,K)-e(i,j+1,nz+1), h2)
    elseif (e(i,j+1,K+1) < e(i,j,nz+1)) then
      if (h2 > h1) h2 = max(e(i,j+1,K)-e(i,j,nz+1), h1)
    endif
    h_def_v(i,J) = h_def_v(i,J) + 0.5*(h1-h2)**2 / ((h1 + h2) + h_neglect)
    h_norm_v(i,J) = h_norm_v(i,J) + 0.5*(h1+h2)
  enddo ; enddo ; enddo
  if (present(def_rat_v_2lay)) then ; do J=js-1,je ; do i=is,ie
    def_rat_v(i,J) = G%mask2dCv(i,J) * h_def_v(i,J) / &
                      (max(CS%Hmix_min, h_norm_v(i,J)) + h_neglect)
    def_rat_v_2lay(i,J) = G%mask2dCv(i,J) * h_def2_v(i,J) / &
                      (max(CS%Hmix_min, h_norm_v(i,J)) + h_neglect)
  enddo ; enddo ; else ; do J=js-1,je ; do i=is,ie
    def_rat_v(i,J) = G%mask2dCv(i,J) * h_def_v(i,J) / &
                      (max(CS%Hmix_min, h_norm_v(i,J)) + h_neglect)
  enddo ; enddo ; endif

end subroutine find_deficit_ratios

subroutine regularize_layers_init(Time, G, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(regularize_layers_CS), pointer    :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                  for this module
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod = "MOM_regularize_layers"  ! This module's name.
  logical :: use_temperature
  integer :: isd, ied, jsd, jed
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (associated(CS)) then
    call MOM_error(WARNING, "regularize_layers_init called with an associated"// &
                            "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%diag => diag
  CS%Time => Time

! Set default, read and log parameters
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "REGULARIZE_SURFACE_LAYERS", CS%regularize_surface_layers, &
                 "If defined, vertically restructure the near-surface \n"//&
                 "layers when they have too much lateral variations to \n"//&
                 "allow for sensible lateral barotropic transports.", &
                 default=.false.)
  if (CS%regularize_surface_layers) then
    call get_param(param_file, mod, "REGULARIZE_SURFACE_DETRAIN", CS%reg_sfc_detrain, &
                 "If true, allow the buffer layers to detrain into the \n"//&
                 "interior as a part of the restructuring when \n"//&
                 "REGULARIZE_SURFACE_LAYERS is true.", default=.true.)
  endif

  call get_param(param_file, mod, "HMIX_MIN", CS%Hmix_min, &
                 "The minimum mixed layer depth if the mixed layer depth \n"//&
                 "is determined dynamically.", units="m", default=0.0)
  call get_param(param_file, mod, "REG_SFC_DEFICIT_TOLERANCE", CS%h_def_tol1, &
                 "The value of the relative thickness deficit at which \n"//&
                 "to start modifying the layer structure when \n"//&
                 "REGULARIZE_SURFACE_LAYERS is true.", units="nondim", &
                 default=0.5)
  CS%h_def_tol2 = 0.2 + 0.8*CS%h_def_tol1
  CS%h_def_tol3 = 0.3 + 0.7*CS%h_def_tol1
  CS%h_def_tol4 = 0.5 + 0.5*CS%h_def_tol1

  call get_param(param_file, mod, "DEBUG", CS%debug, default=.false.)
!  if (.not. CS%debug) &
!    call get_param(param_file, mod, "DEBUG_CONSERVATION", CS%debug, &
!                 "If true, monitor conservation and extrema.", default=.false.)

  call get_param(param_file, mod, "ALLOW_CLOCKS_IN_OMP_LOOPS", &
                 CS%allow_clocks_in_omp_loops, &
                 "If true, clocks can be called from inside loops that can \n"//&
                 "be threaded. To run with multiple threads, set to False.", &
                 default=.true.)

  CS%id_def_rat = register_diag_field('ocean_model', 'deficit_ratio', diag%axesT1, &
      Time, 'Max face thickness deficit ratio', 'Nondim')

#ifdef DEBUG_CODE
  CS%id_def_rat_2 = register_diag_field('ocean_model', 'deficit_rat2', diag%axesT1, &
      Time, 'Corrected thickness deficit ratio', 'Nondim')
  CS%id_def_rat_3 = register_diag_field('ocean_model', 'deficit_rat3', diag%axesT1, &
      Time, 'Filtered thickness deficit ratio', 'Nondim')
  CS%id_e1 = register_diag_field('ocean_model', 'er_1', diag%axesTi, &
      Time, 'Intial interface depths before remapping', 'm')
  CS%id_e2 = register_diag_field('ocean_model', 'er_2', diag%axesTi, &
      Time, 'Intial interface depths after remapping', 'm')
  CS%id_e3 = register_diag_field('ocean_model', 'er_3', diag%axesTi, &
      Time, 'Intial interface depths filtered', 'm')

  CS%id_def_rat_u = register_diag_field('ocean_model', 'defrat_u', diag%axesCu1, &
      Time, 'U-point thickness deficit ratio', 'Nondim')
  CS%id_def_rat_u_1b = register_diag_field('ocean_model', 'defrat_u_1b', diag%axesCu1, &
      Time, 'U-point 2-layer thickness deficit ratio', 'Nondim')
  CS%id_def_rat_u_2 = register_diag_field('ocean_model', 'defrat_u_2', diag%axesCu1, &
      Time, 'U-point corrected thickness deficit ratio', 'Nondim')
  CS%id_def_rat_u_2b = register_diag_field('ocean_model', 'defrat_u_2b', diag%axesCu1, &
      Time, 'U-point corrected 2-layer thickness deficit ratio', 'Nondim')
  CS%id_def_rat_u_3 = register_diag_field('ocean_model', 'defrat_u_3', diag%axesCu1, &
      Time, 'U-point filtered thickness deficit ratio', 'Nondim')
  CS%id_def_rat_u_3b = register_diag_field('ocean_model', 'defrat_u_3b', diag%axesCu1, &
      Time, 'U-point filtered 2-layer thickness deficit ratio', 'Nondim')

  CS%id_def_rat_v = register_diag_field('ocean_model', 'defrat_v', diag%axesCv1, &
      Time, 'V-point thickness deficit ratio', 'Nondim')
  CS%id_def_rat_v_1b = register_diag_field('ocean_model', 'defrat_v_1b', diag%axesCv1, &
      Time, 'V-point 2-layer thickness deficit ratio', 'Nondim')
  CS%id_def_rat_v_2 = register_diag_field('ocean_model', 'defrat_v_2', diag%axesCv1, &
      Time, 'V-point corrected thickness deficit ratio', 'Nondim')
  CS%id_def_rat_v_2b = register_diag_field('ocean_model', 'defrat_v_2b', diag%axesCv1, &
      Time, 'V-point corrected 2-layer thickness deficit ratio', 'Nondim')
  CS%id_def_rat_v_3 = register_diag_field('ocean_model', 'defrat_v_3', diag%axesCv1, &
      Time, 'V-point filtered thickness deficit ratio', 'Nondim')
  CS%id_def_rat_v_3b = register_diag_field('ocean_model', 'defrat_v_3b', diag%axesCv1, &
      Time, 'V-point filtered 2-layer thickness deficit ratio', 'Nondim')
#endif

  if(CS%allow_clocks_in_omp_loops) then
    id_clock_EOS = cpu_clock_id('(Ocean regularize_layers EOS)', grain=CLOCK_ROUTINE)
  endif
  id_clock_pass = cpu_clock_id('(Ocean regularize_layers halo updates)', grain=CLOCK_ROUTINE)

end subroutine regularize_layers_init

end module MOM_regularize_layers
