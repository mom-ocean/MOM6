!> Provides regularization of layers in isopycnal mode
module MOM_regularize_layers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs, EOS_domain

implicit none ; private

#include <MOM_memory.h>

public regularize_layers, regularize_layers_init

!> This control structure holds parameters used by the MOM_regularize_layers module
type, public :: regularize_layers_CS ; private
  logical :: regularize_surface_layers !< If true, vertically restructure the
                             !! near-surface layers when they have too much
                             !! lateral variations to allow for sensible lateral
                             !! barotropic transports.
  logical :: reg_sfc_detrain !< If true, allow the buffer layers to detrain into the
                             !! interior as a part of the restructuring when
                             !! regularize_surface_layers is true
  real    :: density_match_tol !< A relative tolerance for how well the densities must match
                             !! with the target densities during detrainment when regularizing
                             !! the near-surface layers [nondim]
  real    :: h_def_tol1      !< The value of the relative thickness deficit at
                             !! which to start modifying the structure, 0.5 by
                             !! default (or a thickness ratio of 5.83) [nondim].
  real    :: h_def_tol2      !< The value of the relative thickness deficit at
                             !! which to the structure modification is in full
                             !! force, now 20% of the way from h_def_tol1 to 1 [nondim].
  real    :: h_def_tol3      !< The value of the relative thickness deficit at which to start
                             !! detrainment from the buffer layers to the interior, now 30% of
                             !! the way from h_def_tol1 to 1 [nondim].
  real    :: h_def_tol4      !< The value of the relative thickness deficit at which to do
                             !! detrainment from the buffer layers to the interior at full
                             !! force, now 50% of the way from h_def_tol1 to 1 [nondim].
  real    :: Hmix_min        !< The minimum mixed layer thickness [H ~> m or kg m-2].
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                             !! regulate the timing of diagnostic output.
  logical :: answers_2018    !< If true, use the order of arithmetic and expressions that recover the
                             !! answers from the end of 2018.  Otherwise, use updated and more robust
                             !! forms of the same expressions.
  logical :: debug           !< If true, do more thorough checks for debugging purposes.

  integer :: id_def_rat = -1 !< A diagnostic ID
  logical :: allow_clocks_in_omp_loops  !< If true, clocks can be called from inside loops that
                             !! can be threaded. To run with multiple threads, set to False.
end type regularize_layers_CS

!>@{ Clock IDs
!! \todo Should these be global?
integer :: id_clock_pass, id_clock_EOS
!>@}

contains

!> This subroutine partially steps the bulk mixed layer model.
!! The following processes are executed, in the order listed.
subroutine regularize_layers(h, tv, dt, ea, eb, G, GV, US, CS)
  type(ocean_grid_type),      intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),      intent(inout) :: tv !< A structure containing pointers to any
                                                  !! available thermodynamic fields. Absent fields
                                                  !! have NULL ptrs.
  real,                       intent(in)    :: dt !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ea !< The amount of fluid moved downward into a
                                                  !! layer; this should be increased due to mixed
                                                  !! layer detrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eb !< The amount of fluid moved upward into a layer
                                                  !! this should be increased due to mixed layer
                                                  !! entrainment [H ~> m or kg m-2].
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(regularize_layers_CS), pointer       :: CS !< The control structure returned by a previous
                                                  !! call to regularize_layers_init.
  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_regularize_layers: "//&
         "Module must be initialized before it is used.")

  if (CS%regularize_surface_layers) then
    call pass_var(h, G%Domain, clock=id_clock_pass)
    call regularize_surface(h, tv, dt, ea, eb, G, GV, US, CS)
  endif

end subroutine regularize_layers

!> This subroutine ensures that there is a degree of horizontal smoothness
!! in the depths of the near-surface interfaces.
subroutine regularize_surface(h, tv, dt, ea, eb, G, GV, US, CS)
  type(ocean_grid_type),      intent(inout) :: G  !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2].
  type(thermo_var_ptrs),      intent(inout) :: tv !< A structure containing pointers to any
                                                  !! available thermodynamic fields. Absent fields
                                                  !! have NULL ptrs.
  real,                       intent(in)    :: dt !< Time increment [T ~> s].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: ea !< The amount of fluid moved downward into a
                                                  !! layer; this should be increased due to mixed
                                                  !! layer detrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(inout) :: eb !< The amount of fluid moved upward into a layer
                                                  !! this should be increased due to mixed layer
                                                  !! entrainment [H ~> m or kg m-2].
  type(unit_scale_type),      intent(in)    :: US !< A dimensional unit scaling type
  type(regularize_layers_CS), pointer       :: CS !< The control structure returned by a previous
                                                  !! call to regularize_layers_init.
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    def_rat_u   ! The ratio of the thickness deficit to the minimum depth [nondim].
  real, dimension(SZI_(G),SZJB_(G)) :: &
    def_rat_v   ! The ratio of the thickness deficit to the minimum depth [nondim].
  real, dimension(SZI_(G),SZJ_(G)) :: &
    def_rat_h   ! The ratio of the thickness deficit to the minimum depth [nondim].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: &
    e           ! The interface depths [H ~> m or kg m-2], positive upward.

  real, dimension(SZI_(G),SZK_(GV)+1) :: &
    e_filt, e_2d  ! The interface depths [H ~> m or kg m-2], positive upward.
  real, dimension(SZI_(G),SZK_(GV)) :: &
    h_2d, &     !   A 2-d version of h [H ~> m or kg m-2].
    T_2d, &     !   A 2-d version of tv%T [degC].
    S_2d, &     !   A 2-d version of tv%S [ppt].
    Rcv, &      !   A 2-d version of the coordinate density [R ~> kg m-3].
    h_2d_init, &  ! The initial value of h_2d [H ~> m or kg m-2].
    T_2d_init, &  ! THe initial value of T_2d [degC].
    S_2d_init, &  ! The initial value of S_2d [ppt].
    d_eb, &     !   The downward increase across a layer in the entrainment from
                ! below [H ~> m or kg m-2].  The sign convention is that positive values of
                ! d_eb correspond to a gain in mass by a layer by upward motion.
    d_ea        !   The upward increase across a layer in the entrainment from
                ! above [H ~> m or kg m-2].  The sign convention is that positive values of
                ! d_ea mean a net gain in mass by a layer from downward motion.
  real, dimension(SZI_(G)) :: &
    p_ref_cv, & !   Reference pressure for the potential density which defines
                ! the coordinate variable, set to P_Ref [R L2 T-2 ~> Pa].
    Rcv_tol, &  !   A tolerence, relative to the target density differences
                ! between layers, for detraining into the interior [nondim].
    h_add_tgt, h_add_tot, &
    h_tot1, Th_tot1, Sh_tot1, &
    h_tot3, Th_tot3, Sh_tot3, &
    h_tot2, Th_tot2, Sh_tot2
  real, dimension(SZK_(GV)) :: &
    h_prev_1d     ! The previous thicknesses [H ~> m or kg m-2].
  real :: I_dtol  ! The inverse of the tolerance changes [nondim].
  real :: I_dtol34 ! The inverse of the tolerance changes [nondim].
  real :: h1, h2  ! Temporary thicknesses [H ~> m or kg m-2].
  real :: e_e, e_w, e_n, e_s  ! Temporary interface heights [H ~> m or kg m-2].
  real :: wt    ! The weight of the filted interfaces in setting the targets [nondim].
  real :: scale ! A scaling factor [nondim].
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real, dimension(SZK_(GV)+1) :: &
    int_flux, int_Tflux, int_Sflux, int_Rflux
  real :: h_add
  real :: h_det_tot
  real :: max_def_rat
  real :: Rcv_min_det  ! The lightest (min) and densest (max) coordinate density
  real :: Rcv_max_det  ! that can detrain into a layer [R ~> kg m-3].

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
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz, nkmb, nkml, k1, k2, k3, ks, nz_filt, kmax_d_ea

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_regularize_layers: "//&
         "Module must be initialized before it is used.")

  if (GV%nkml<1) return
  nkmb = GV%nk_rho_varies ; nkml = GV%nkml
  if (.not.associated(tv%eqn_of_state)) call MOM_error(FATAL, &
    "MOM_regularize_layers: This module now requires the use of temperature and "//&
    "an equation of state.")

  h_neglect = GV%H_subroundoff
  debug = (debug .or. CS%debug)

  I_dtol = 1.0 / max(CS%h_def_tol2 - CS%h_def_tol1, 1e-40)
  I_dtol34 = 1.0 / max(CS%h_def_tol4 - CS%h_def_tol3, 1e-40)

  p_ref_cv(:) = tv%P_Ref
  EOSdom(:) = EOS_domain(G%HI)

  do j=js-1,je+1 ; do i=is-1,ie+1
    e(i,j,1) = 0.0
  enddo ; enddo
  do K=1,nz ; do j=js-1,je+1 ; do i=is-1,ie+1
    e(i,j,K+1) = e(i,j,K) - h(i,j,k)
  enddo ; enddo ; enddo

  call find_deficit_ratios(e, def_rat_u, def_rat_v, G, GV, CS, h=h)

  ! Determine which columns are problematic
  do j=js,je ; do_j(j) = .false. ; enddo
  do j=js,je ; do i=is,ie
    def_rat_h(i,j) = max(def_rat_u(I-1,j), def_rat_u(I,j), &
                         def_rat_v(i,J-1), def_rat_v(i,J))
    if (def_rat_h(i,j) > CS%h_def_tol1) do_j(j) = .true.
  enddo ; enddo

  ! Now restructure the layers.
  !$OMP parallel do default(private) shared(is,ie,js,je,nz,do_j,def_rat_h,CS,nkmb,G,GV,US, &
  !$OMP                                     e,I_dtol,h,tv,debug,h_neglect,p_ref_cv,ea, &
  !$OMP                                     eb,id_clock_EOS,nkml,EOSdom)
  do j=js,je ; if (do_j(j)) then

!  call cpu_clock_begin(id_clock_EOS)
!  call calculate_density_derivs(T(:,1), S(:,1), p_ref_cv, dRcv_dT, dRcv_dS, tv%eqn_of_state, EOSdom)
!  call cpu_clock_end(id_clock_EOS)

    do k=1,nz ; do i=is,ie ; d_ea(i,k) = 0.0 ; d_eb(i,k) = 0.0 ; enddo ; enddo
    kmax_d_ea = 0

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
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom_H)

      endif
      if (G%mask2dCu(I-1,j) <= 0.0) then ; e_w = e(i,j,K) ; else
        e_w = max(e(i-1,j,K) + min(e(i,j,K) - e(i-1,j,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom_H)
      endif
      if (G%mask2dCv(i,J) <= 0.0) then ; e_n = e(i,j,K) ; else
        e_n = max(e(i,j+1,K) + min(e(i,j,K) - e(i,j+1,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom_H)
      endif
      if (G%mask2dCv(i,J-1) <= 0.0) then ; e_s = e(i,j,K) ; else
        e_s = max(e(i,j-1,K) + min(e(i,j,K) - e(i,j-1,nz+1), 0.0), &
                  e(i,j,nz+1) + (nz+1-k)*GV%Angstrom_H)
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
          if (h_2d(i,k) - GV%Angstrom_H > h_neglect) then
            if (e_2d(i,nkmb+1)-e_filt(i,nkmb+1) > h_2d(i,k) - GV%Angstrom_H) then
              h_add = h_2d(i,k) - GV%Angstrom_H
              h_2d(i,k) = GV%Angstrom_H
              e_2d(i,nkmb+1) = e_2d(i,nkmb+1) - h_add
            else
              h_add = e_2d(i,nkmb+1) - e_filt(i,nkmb+1)
              h_2d(i,k) = h_2d(i,k) - h_add
              if (CS%answers_2018) then
                e_2d(i,nkmb+1) = e_2d(i,nkmb+1) - h_add
              else
                e_2d(i,nkmb+1) = e_filt(i,nkmb+1)
              endif
            endif
            d_eb(i,k-1) = d_eb(i,k-1) + h_add
            h_add_tot(i) = h_add_tot(i) + h_add
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
          ! The CS%density_match_tol default value of 0.6 gives 20% overlap in acceptable densities.
          Rcv_tol(i) = CS%density_match_tol * min((def_rat_h(i,j) - CS%h_def_tol3), 1.0)
        endif
      enddo
    endif
    if (det_any) then
      call cpu_clock_begin(id_clock_EOS)
      do k=1,nkmb
        call calculate_density(T_2d(:,k), S_2d(:,k), p_ref_cv, Rcv(:,k), tv%eqn_of_state, EOSdom)
      enddo
      call cpu_clock_end(id_clock_EOS)

      do i=is,ie ; if (det_i(i)) then
        k1 = nkmb ; k2 = nz
        h_det_tot = 0.0
        do ! This loop is terminated by exits.
          if (k1 <= 1) exit
          if (k2 <= nkmb) exit
          Rcv_min_det = (GV%Rlay(k2) + Rcv_tol(i)*(GV%Rlay(k2-1)-GV%Rlay(k2)))
          if (k2 < nz) then
            Rcv_max_det = (GV%Rlay(k2) + Rcv_tol(i)*(GV%Rlay(k2+1)-GV%Rlay(k2)))
          else
            Rcv_max_det = (GV%Rlay(nz) + Rcv_tol(i)*(GV%Rlay(nz)-GV%Rlay(nz-1)))
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
                kmax_d_ea = max(kmax_d_ea, k2)
                ! This is upwind.  It should perhaps be higher order...
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
              kmax_d_ea = max(kmax_d_ea, k2)
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
      do k=kmax_d_ea-1,nkmb+1,-1 ; do i=is,ie ; if (det_i(i)) then
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

      ! Map the water back into the layers.  There are not mixed or buffer layers that are exceedingly
      ! small compared to the others, so the code here is less prone to roundoff than elsewhere in MOM6.
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
        if (abs(h(i,j,k) - h_predicted) > MAX(1e-9*abs(h_predicted),GV%Angstrom_H)) &
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

end subroutine regularize_surface

!>  This subroutine determines the amount by which the harmonic mean
!! thickness at velocity points differ from the arithmetic means, relative to
!! the the arithmetic means, after eliminating thickness variations that are
!! solely due to topography and aggregating all interior layers into one.
subroutine find_deficit_ratios(e, def_rat_u, def_rat_v, G, GV, CS, &
                               def_rat_u_2lay, def_rat_v_2lay, halo, h)
  type(ocean_grid_type),      intent(in)  :: G         !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in)  :: GV        !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1), &
                              intent(in)  :: e         !< Interface depths [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),          &
                              intent(out) :: def_rat_u !< The thickness deficit ratio at u points,
                                                       !! [nondim].
  real, dimension(SZI_(G),SZJB_(G)),          &
                              intent(out) :: def_rat_v !< The thickness deficit ratio at v points,
                                                       !! [nondim].
  type(regularize_layers_CS), pointer     :: CS        !< The control structure returned by a
                                                       !! previous call to regularize_layers_init.
  real, dimension(SZIB_(G),SZJ_(G)),          &
                    optional, intent(out) :: def_rat_u_2lay !< The thickness deficit ratio at u
                                                       !! points when the mixed and buffer layers
                                                       !! are aggregated into 1 layer [nondim].
  real, dimension(SZI_(G),SZJB_(G)),          &
                    optional, intent(out) :: def_rat_v_2lay !< The thickness deficit ratio at v
                                                       !! pointswhen the mixed and buffer layers
                                                       !! are aggregated into 1 layer [nondim].
  integer,          optional, intent(in)  :: halo      !< An extra-wide halo size, 0 by default.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  &
                    optional, intent(in)  :: h         !< Layer thicknesses [H ~> m or kg m-2].
                                                       !! If h is not present, vertical differences
                                                       !! in interface heights are used instead.
  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G)) :: &
    h_def_u, &  ! The vertically summed thickness deficits at u-points [H ~> m or kg m-2].
    h_norm_u, & ! The vertically summed arithmetic mean thickness by which
                ! h_def_u is normalized [H ~> m or kg m-2].
    h_def2_u
  real, dimension(SZI_(G),SZJB_(G)) :: &
    h_def_v, &  ! The vertically summed thickness deficits at v-points [H ~> m or kg m-2].
    h_norm_v, & ! The vertically summed arithmetic mean thickness by which
                ! h_def_v is normalized [H ~> m or kg m-2].
    h_def2_v
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: Hmix_min  ! A local copy of CS%Hmix_min [H ~> m or kg m-2].
  real :: h1, h2  ! Temporary thicknesses [H ~> m or kg m-2].
  integer :: i, j, k, is, ie, js, je, nz, nkmb

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif
  nkmb = GV%nk_rho_varies
  h_neglect = GV%H_subroundoff
  Hmix_min = CS%Hmix_min

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
                     (max(Hmix_min, h_norm_u(I,j)) + h_neglect)
    def_rat_u_2lay(I,j) = G%mask2dCu(I,j) * h_def2_u(I,j) / &
                          (max(Hmix_min, h_norm_u(I,j)) + h_neglect)
  enddo ; enddo ; else ; do j=js,je ; do I=is-1,ie
    def_rat_u(I,j) = G%mask2dCu(I,j) * h_def_u(I,j) / &
                     (max(Hmix_min, h_norm_u(I,j)) + h_neglect)
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
                      (max(Hmix_min, h_norm_v(i,J)) + h_neglect)
    def_rat_v_2lay(i,J) = G%mask2dCv(i,J) * h_def2_v(i,J) / &
                      (max(Hmix_min, h_norm_v(i,J)) + h_neglect)
  enddo ; enddo ; else ; do J=js-1,je ; do i=is,ie
    def_rat_v(i,J) = G%mask2dCv(i,J) * h_def_v(i,J) / &
                      (max(Hmix_min, h_norm_v(i,J)) + h_neglect)
  enddo ; enddo ; endif

end subroutine find_deficit_ratios

!> Initializes the regularize_layers control structure
subroutine regularize_layers_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< The current model time.
  type(ocean_grid_type),   intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for
                                                 !! run-time parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< A structure that is used to regulate
                                                 !! diagnostic output.
  type(regularize_layers_CS), pointer    :: CS   !< A pointer that is set to point to the
                                                 !! control structure for this module.
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_regularize_layers"  ! This module's name.
  logical :: use_temperature
  logical :: default_2018_answers
  logical :: just_read
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
  call get_param(param_file, mdl, "REGULARIZE_SURFACE_LAYERS", CS%regularize_surface_layers, &
                 default=.false., do_not_log=.true.)
  call log_version(param_file, mdl, version, "", all_default=.not.CS%regularize_surface_layers)
  call get_param(param_file, mdl, "REGULARIZE_SURFACE_LAYERS", CS%regularize_surface_layers, &
                 "If defined, vertically restructure the near-surface "//&
                 "layers when they have too much lateral variations to "//&
                 "allow for sensible lateral barotropic transports.", &
                 default=.false.)
  just_read = .not.CS%regularize_surface_layers
  if (CS%regularize_surface_layers) then
    call get_param(param_file, mdl, "REGULARIZE_SURFACE_DETRAIN", CS%reg_sfc_detrain, &
                 "If true, allow the buffer layers to detrain into the "//&
                 "interior as a part of the restructuring when "//&
                 "REGULARIZE_SURFACE_LAYERS is true.", default=.true., do_not_log=just_read)
    call get_param(param_file, mdl, "REG_SFC_DENSE_MATCH_TOLERANCE", CS%density_match_tol, &
                 "A relative tolerance for how well the densities must match with the target "//&
                 "densities during detrainment when regularizing the near-surface layers.  The "//&
                 "default of 0.6 gives 20% overlaps in density", &
                 units="nondim", default=0.6, do_not_log=just_read)
    call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.false., do_not_log=just_read)
    call get_param(param_file, mdl, "REGULARIZE_LAYERS_2018_ANSWERS", CS%answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the answers "//&
                 "from the end of 2018.  Otherwise, use updated and more robust forms of the "//&
                 "same expressions.", default=default_2018_answers, do_not_log=just_read)
  endif

  call get_param(param_file, mdl, "HMIX_MIN", CS%Hmix_min, &
                 "The minimum mixed layer depth if the mixed layer depth is determined "//&
                 "dynamically.", units="m", default=0.0, scale=GV%m_to_H, do_not_log=just_read)
  call get_param(param_file, mdl, "REG_SFC_DEFICIT_TOLERANCE", CS%h_def_tol1, &
                 "The value of the relative thickness deficit at which "//&
                 "to start modifying the layer structure when "//&
                 "REGULARIZE_SURFACE_LAYERS is true.", units="nondim", &
                 default=0.5, do_not_log=just_read)
  CS%h_def_tol2 = 0.2 + 0.8*CS%h_def_tol1
  CS%h_def_tol3 = 0.3 + 0.7*CS%h_def_tol1
  CS%h_def_tol4 = 0.5 + 0.5*CS%h_def_tol1

  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
!  if (.not. CS%debug) &
!    call get_param(param_file, mdl, "DEBUG_CONSERVATION", CS%debug, &
!                 "If true, monitor conservation and extrema.", default=.false., do_not_log=just_read)

  call get_param(param_file, mdl, "ALLOW_CLOCKS_IN_OMP_LOOPS", CS%allow_clocks_in_omp_loops, &
                 "If true, clocks can be called from inside loops that can "//&
                 "be threaded. To run with multiple threads, set to False.", &
                 default=.true., do_not_log=just_read)

  if (.not.CS%regularize_surface_layers) return

  CS%id_def_rat = register_diag_field('ocean_model', 'deficit_ratio', diag%axesT1, &
      Time, 'Max face thickness deficit ratio', 'nondim')

  if (CS%allow_clocks_in_omp_loops) then
    id_clock_EOS = cpu_clock_id('(Ocean regularize_layers EOS)', grain=CLOCK_ROUTINE)
  endif
  id_clock_pass = cpu_clock_id('(Ocean regularize_layers halo updates)', grain=CLOCK_ROUTINE)

end subroutine regularize_layers_init

end module MOM_regularize_layers
