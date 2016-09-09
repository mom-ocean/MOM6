module MOM_entrain_diffusive
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
!*  By Robert Hallberg, September 1997 - July 2000                     *
!*                                                                     *
!*    This file contains the subroutines that implement diapycnal      *
!*  mixing and advection in isopycnal layers.  The main subroutine,    *
!*  calculate_entrainment, returns the entrainment by each layer       *
!*  across the interfaces above and below it.  These are calculated    *
!*  subject to the constraints that no layers can be driven to neg-    *
!*  ative thickness and that the each layer maintains its target       *
!*  density, using the scheme described in Hallberg (MWR 2000). There  *
!*  may or may not be a bulk mixed layer above the isopycnal layers.   *
!*  The solution is iterated until the change in the entrainment       *
!*  between successive iterations is less than some small tolerance.   *
!*                                                                     *
!*    The dual-stream entrainment scheme of MacDougall and Dewar       *
!*  (JPO 1997) is used for combined diapycnal advection and diffusion, *
!*  modified as described in Hallberg (MWR 2000) to be solved          *
!*  implicitly in time.  Any profile of diffusivities may be used.     *
!*  Diapycnal advection is fundamentally the residual of diapycnal     *
!*  diffusion, so the fully implicit upwind differencing scheme that   *
!*  is used is entirely appropriate.  The downward buoyancy flux in    *
!*  each layer is determined from an implicit calculation based on     *
!*  the previously calculated flux of the layer above and an estim-    *
!*  ated flux in the layer below.  This flux is subject to the foll-   *
!*  owing conditions:  (1) the flux in the top and bottom layers are   *
!*  set by the boundary conditions, and (2) no layer may be driven     *
!*  below an Angstrom thickness.  If there is a bulk mixed layer, the  *
!*  mixed and buffer layers are treated as Eulerian layers, whose      *
!*  thicknesses only change due to entrainment by the interior layers. *
!*                                                                     *
!*    In addition, the model may adjust the fluxes to drive the layer  *
!*  densities (sigma 2?) back toward their targer values.              *
!*                                                                     *
!*     A small fragment of the grid is shown below:                    *
!*                                                                     *
!*    j+1  x ^ x ^ x   At x:  q                                        *
!*    j+1  > o > o >   At ^:  v                                        *
!*    j    x ^ x ^ x   At >:  u                                        *
!*    j    > o > o >   At o:  h, buoy, T, S, ea, eb, etc.              *
!*    j-1  x ^ x ^ x                                                   *
!*        i-1  i  i+1  At x & ^:                                       *
!*           i  i+1    At > & o:                                       *
!*                                                                     *
!*  The boundaries always run through q grid points (x).               *
!*                                                                     *
!********+*********+*********+*********+*********+*********+*********+**

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl, time_type
use MOM_error_handler, only : MOM_error, is_root_pe, FATAL, WARNING, NOTE
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_EOS, only : calculate_density, calculate_density_derivs

implicit none ; private

#include <MOM_memory.h>

public entrainment_diffusive, entrain_diffusive_init, entrain_diffusive_end

type, public :: entrain_diffusive_CS ; private
  logical :: bulkmixedlayer  ! If true, a refined bulk mixed layer is used with
                             ! GV%nk_rho_varies variable density mixed & buffer
                             ! layers.
  logical :: correct_density ! If true, the layer densities are restored toward
                             ! their target variables by the diapycnal mixing.
  integer :: max_ent_it      ! The maximum number of iterations that may be
                             ! used to calculate the diapycnal entrainment.
  real    :: Tolerance_Ent   ! The tolerance with which to solve for entrainment
                             ! values, in m.
  type(diag_ctrl), pointer :: diag ! A structure that is used to regulate the
                             ! timing of diagnostic output.
  integer :: id_Kd = -1, id_diff_work = -1
end type entrain_diffusive_CS

contains

subroutine entrainment_diffusive(u, v, h, tv, fluxes, dt, G, GV, CS, ea, eb, &
                                 kb_out, Kd_Lay, Kd_int)
  type(ocean_grid_type),                     intent(in)  :: G
  type(verticalGrid_type),                   intent(in)  :: GV
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(in)  :: u
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(in)  :: v
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)  :: h
  type(thermo_var_ptrs),                     intent(in)  :: tv
  type(forcing),                             intent(in)  :: fluxes
  real,                                      intent(in)  :: dt
  type(entrain_diffusive_CS),                pointer     :: CS
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(out) :: ea, eb
  integer, dimension(SZI_(G),SZJ_(G)),        optional, intent(inout) :: kb_out
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   optional, intent(in) :: Kd_Lay
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), optional, intent(in) :: Kd_int

!   This subroutine calculates ea and eb, the rates at which a layer
! entrains from the layers above and below.  The entrainment rates
! are proportional to the buoyancy flux in a layer and inversely
! proportional to the density differences between layers.  The
! scheme that is used here is described in detail in Hallberg, Mon.
! Wea. Rev. 2000.

! Arguments: u - Zonal velocity, in m s-1.
!  (in)      v - Meridional velocity, in m s-1.
!  (in)      h - Layer thickness, in m or kg m-2.
!  (in)      fluxes - A structure of surface fluxes that may be used.
!  (in)      kb_out - The index of the lightest layer denser than the
!                     buffer layers.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      dt - The time increment in s.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - The control structure returned by a previous call to
!                 entrain_diffusive_init.
!  (out)     ea - The amount of fluid entrained from the layer above within
!                 this time step, in the same units as h, m or kg m-2.
!  (out)     eb - The amount of fluid entrained from the layer below within
!                 this time step, in the same units as h, m or kg m-2.
!  (out,opt) kb - The index of the lightest layer denser than the buffer layer.
! At least one of the two arguments must be present.
!  (in,opt)  Kd_Lay - The diapycnal diffusivity of layers, in m2 s-1.
!  (in,opt)  Kd_int - The diapycnal diffusivity of interfaces, in m2 s-1.

! In the comments below, H is used as shorthand for the units of h, m or kg m-2.
  real, dimension(SZI_(G),SZK_(G)) :: &
    dtKd    ! The layer diapycnal diffusivity times the time step, translated
            ! into the same unints as h, m2 or kg2 m-4 (i.e. H2).
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    dtKd_int ! The diapycnal diffusivity at the interfaces times the time step,
            ! translated into the same unints as h, m2 or kg2 m-4 (i.e. H2).
  real, dimension(SZI_(G),SZK_(G)) :: &
    F, &    ! The density flux through a layer within a time step divided by the
            ! density difference across the interface below the layer, in H.
    maxF, & ! maxF is the maximum value of F that will not deplete all of the
            ! layers above or below a layer within a timestep, in H.
    minF, & ! minF is the minimum flux that should be expected in the absence of
            ! interactions between layers, in H.
    Fprev, &! The previous estimate of F, in H.
    dFdfm, &! The partial derivative of F with respect to changes in F of the
            ! neighboring layers.  Nondimensional.
    h_guess ! An estimate of the layer thicknesses after entrainment, but
            ! before the entrainments are adjusted to drive the layer
            ! densities toward their target values, in H.
  real, dimension(SZI_(G),SZK_(G)+1) :: &
    Ent_bl  ! The average entrainment upward and downward across
            ! each interface around the buffer layers, in H.
  real, allocatable, dimension(:,:,:) :: &
    Kd_eff, &     ! The effective diffusivity that actually applies to each
                  ! layer after the effects of boundary conditions are
                  ! considered, in m2 s-1.
    diff_work     ! The work actually done by diffusion across each
                  ! interface, in W m-2.  Sum vertically for the total work.

  real :: hm, fm, fr, fk  ! Work variables with units of H, H, H, and H2.

  real :: b1(SZI_(G))         ! b1 and c1 are variables used by the
  real :: c1(SZI_(G),SZK_(G)) ! tridiagonal solver.

  real, dimension(SZI_(G)) :: &
    htot, &       ! The total thickness above or below a layer in H.
    Rcv, &        ! Value of the coordinate variable (potential density)
                  ! based on the simulated T and S and P_Ref, kg m-3.
    pres, &       ! Reference pressure (P_Ref) in Pa.
    eakb, &       ! The entrainment from above by the layer below the buffer
                  ! layer (i.e. layer kb), in H.
    ea_kbp1, &    ! The entrainment from above by layer kb+1, in H.
    eb_kmb, &     ! The entrainment from below by the deepest buffer layer, in H.
    dS_kb, &      ! The reference potential density difference across the
                  ! interface between the buffer layers and layer kb, in kg m-3.
    dS_anom_lim, &! The amount by which dS_kb is reduced when limits are
                  ! applied, in kg m-3.
    I_dSkbp1, &   ! The inverse of the potential density difference across the
                  ! interface below layer kb, in m3 kg-1.
    dtKd_kb, &    ! The diapycnal diffusivity in layer kb times the time step,
                  ! in units of H2.
    maxF_correct, & ! An amount by which to correct maxF due to excessive
                  ! surface heat loss, in H.
    zeros, &      ! An array of all zeros. (Usually used with units of H.)
    max_eakb, &   ! The maximum value of eakb that might be realized, in H.
    min_eakb, &   ! The minimum value of eakb that might be realized, in H.
    err_max_eakb0, & ! The value of error returned by determine_Ea_kb
    err_min_eakb0, & ! when eakb = min_eakb and max_eakb and ea_kbp1 = 0.
    err_eakb0, &  ! A value of error returned by determine_Ea_kb.
    F_kb, &       ! The value of F in layer kb, or equivalently the entrainment
                  ! from below by layer kb, in H.
    dFdfm_kb, &   ! The partial derivative of F with fm, nondim. See dFdfm.
    maxF_kb, &    ! The maximum value of F_kb that might be realized, in H.
    eakb_maxF, &  ! The value of eakb that gives F_kb=maxF_kb, in H.
    F_kb_maxEnt   ! The value of F_kb when eakb = max_eakb, in H.
  real, dimension(SZI_(G),SZK_(G)) :: &
    Sref, &  ! The reference potential density of the mixed and buffer layers,
             ! and of the two lightest interior layers (kb and kb+1) copied
             ! into layers kmb+1 and kmb+2, in kg m-3.
    h_bl     ! The thicknesses of the mixed and buffer layers, and of the two
             ! lightest interior layers (kb and kb+1) copied into layers kmb+1
             ! and kmb+2, in H.

  real, dimension(SZI_(G),SZK_(G)) :: &
    ds_dsp1, &      ! The coordinate variable (sigma-2) difference across an
                    ! interface divided by the difference across the interface
                    ! below it. Nondimensional.
    dsp1_ds, &      ! The inverse coordinate variable (sigma-2) difference
                    ! across an interface times the difference across the
                    ! interface above it. Nondimensional.
    I2p2dsp1_ds, &  ! 1 / (2 + 2 * ds_k+1 / ds_k). Nondimensional.
    grats           ! 2*(2 + ds_k+1 / ds_k + ds_k / ds_k+1) =
                    !       4*ds_Lay*(1/ds_k + 1/ds_k+1). Nondim.

  real :: dRHo      ! The change in locally referenced potential density between
                    ! the layers above and below an interface, in kg m-3.
  real :: g_2dt     ! 0.5 * G_Earth / dt, in m s-3.
  real, dimension(SZI_(G)) :: &
    pressure, &      ! The pressure at an interface, in Pa.
    T_eos, S_eos, &  ! The potential temperature and salinity at which to
                     ! evaluate dRho_dT and dRho_dS, in degC and PSU.
    dRho_dT, dRho_dS ! The partial derivatives of potential density with
                     ! temperature and salinity, in kg m-3 K-1 and kg m-3 psu-1.

  real :: tolerance  ! The tolerance within which E must be converged, in H.
  real :: Angstrom   ! The minimum layer thickness, in H.
  real :: h_neglect  ! A thickness that is so small it is usually lost
                     ! in roundoff and can be neglected, in H.
  real :: F_cor      ! A correction to the amount of F that is used to
                     ! entrain from the layer above, in H.
  real :: Kd_here    ! The effective diapycnal diffusivity, in H2 s-1.
  real :: h_avail    ! The thickness that is available for entrainment, in H.
  real :: dS_kb_eff  ! The value of dS_kb after limiting is taken into account.
  real :: Rho_cor    ! The depth-integrated potential density anomaly that
                     ! needs to be corrected for, in kg m-2.
  real :: ea_cor     ! The corrective adjustment to eakb, in H.
  real :: h1         ! The layer thickness after entrainment through the
                     ! interface below is taken into account, in H.
  real :: Idt        ! The inverse of the time step, in s-1.
  real :: H_to_m, m_to_H  ! Local copies of unit conversion factors.

  logical :: do_any
  logical :: do_i(SZI_(G)), did_i(SZI_(G)), reiterate, correct_density
  integer :: it, i, j, k, is, ie, js, je, nz, K2, kmb
  integer :: kb(SZI_(G))  ! The value of kb in row j.
  integer :: kb_min       ! The minimum value of kb in the current j-row.
  integer :: kb_min_act   ! The minimum active value of kb in the current j-row.
  integer :: is1, ie1     ! The minimum and maximum active values of i in the current j-row.
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  Angstrom = GV%Angstrom
  h_neglect = GV%H_subroundoff

  if (.not. associated(CS)) call MOM_error(FATAL, &
         "MOM_entrain_diffusive: Module must be initialized before it is used.")

  if (.not.(present(Kd_Lay) .or. present(Kd_int))) call MOM_error(FATAL, &
      "MOM_entrain_diffusive: Either Kd_Lay or Kd_int must be present in call.")

  if ((.not.CS%bulkmixedlayer .and. .not.ASSOCIATED(fluxes%buoy)) .and. &
      (ASSOCIATED(fluxes%lprec) .or. ASSOCIATED(fluxes%evap) .or. &
       ASSOCIATED(fluxes%sens) .or. ASSOCIATED(fluxes%sw))) then
    if (is_root_pe()) call MOM_error(NOTE, "Calculate_Entrainment: &
          &The code to handle evaporation and precipitation without &
          &a bulk mixed layer has not been implemented.")
    if (is_root_pe()) call MOM_error(FATAL, &
         "Either define BULKMIXEDLAYER in MOM_input or use fluxes%buoy &
         &and a linear equation of state to drive the model.")
  endif

  H_to_m = GV%H_to_m ; m_to_H = GV%m_to_H
  tolerance = m_to_H * CS%Tolerance_Ent
  g_2dt = 0.5 * GV%g_Earth / dt
  kmb = GV%nk_rho_varies
  K2 = max(kmb+1,2) ; kb_min = K2
  if (.not. CS%bulkmixedlayer) then
    kb(:) = 1
    ! These lines fill in values that are arbitrary, but needed because
    ! they are used to normalize the buoyancy flux in layer nz.
    do i=is,ie ; ds_dsp1(i,nz) = 2.0 ; dsp1_ds(i,nz) = 0.5 ; enddo
  else
    kb(:) = 0
    do i=is,ie ; ds_dsp1(i,nz) = 0.0 ; dsp1_ds(i,nz) = 0.0 ; enddo
  endif

  if (CS%id_diff_work > 0) allocate(diff_work(G%isd:G%ied,G%jsd:G%jed,nz+1))
  if (CS%id_Kd > 0)        allocate(Kd_eff(G%isd:G%ied,G%jsd:G%jed,nz))

  correct_density = (CS%correct_density .and. associated(tv%eqn_of_state))
  if (correct_density) then
    pres(:) = tv%P_Ref
  else
    pres(:) = 0.0
  endif

!$OMP parallel do default(none) shared(is,ie,js,je,nz,Kd_Lay,G,GV,dt,Kd_int,CS,h,tv, &
!$OMP                                  kmb,Angstrom,fluxes,K2,h_neglect,tolerance, &
!$OMP                                  ea,eb,correct_density,Kd_eff,diff_work,     &
!$OMP                                  g_2dt, kb_out, m_to_H, H_to_m)              &
!$OMP                     firstprivate(kb,ds_dsp1,dsp1_ds,pres,kb_min)             &
!$OMP                          private(dtKd,dtKd_int,do_i,Ent_bl,dtKd_kb,h_bl,     &
!$OMP                                  I2p2dsp1_ds,grats,htot,max_eakb,I_dSkbp1,   &
!$OMP                                  zeros,maxF_kb,maxF,ea_kbp1,eakb,Sref,       &
!$OMP                                  maxF_correct,do_any,                        &
!$OMP                                  err_min_eakb0,err_max_eakb0,eakb_maxF,      &
!$OMP                                  min_eakb,err_eakb0,F,minF,hm,fk,F_kb_maxent,&
!$OMP                                  F_kb,is1,ie1,kb_min_act,dFdfm_kb,b1,dFdfm,  &
!$OMP                                  Fprev,fm,fr,c1,reiterate,eb_kmb,did_i,      &
!$OMP                                  h_avail,h_guess,dS_kb,Rcv,F_cor,dS_kb_eff,  &
!$OMP                                  Rho_cor,ea_cor,h1,Idt,Kd_here,pressure,     &
!$OMP                                  T_eos,S_eos,dRho_dT,dRho_dS,dRho,dS_anom_lim)
  do j=js,je
    do i=is,ie ; kb(i) = 1 ; enddo

    if (present(Kd_Lay)) then
      do k=1,nz ; do i=is,ie
        dtKd(i,k) = m_to_H**2 * (dt*Kd_Lay(i,j,k))
      enddo ; enddo
      if (present(Kd_int)) then
        do K=1,nz+1 ; do i=is,ie
          dtKd_int(i,K) = m_to_H**2 * (dt*Kd_int(i,j,K))
        enddo ; enddo
      else
        do K=2,nz ; do i=is,ie
          dtKd_int(i,K) = m_to_H**2 * (0.5*dt*(Kd_Lay(i,j,k-1) + Kd_Lay(i,j,k)))
        enddo ; enddo
      endif
    else ! Kd_int must be present, or there already would have been an error.
      do k=1,nz ; do i=is,ie
        dtKd(i,k) = m_to_H**2 * (0.5*dt*(Kd_int(i,j,K)+Kd_int(i,j,K+1)))
      enddo ; enddo
      dO K=1,nz+1 ; do i=is,ie
        dtKd_int(i,K) = m_to_H**2 * (dt*Kd_int(i,j,K))
      enddo ; enddo
    endif

    do i=is,ie ; do_i(i) = (G%mask2dT(i,j) > 0.5) ; enddo
    do i=is,ie ; ds_dsp1(i,nz) = 0.0 ; enddo
    do i=is,ie ; dsp1_ds(i,nz) = 0.0 ; enddo

    do k=2,nz-1 ; do i=is,ie
      ds_dsp1(i,k) = GV%g_prime(k) / GV%g_prime(k+1)
    enddo ; enddo

    if (CS%bulkmixedlayer) then
      !   This subroutine determines the averaged entrainment across each
      ! interface and causes thin and relatively light interior layers to be
      ! entrained by the deepest buffer layer.  This also determines kb.
      call set_Ent_bl(h, dtKd_int, tv, kb, kmb, do_i, G, GV, CS, j, Ent_bl, Sref, h_bl)

      do i=is,ie
        dtKd_kb(i) = 0.0 ; if (kb(i) < nz) dtKd_kb(i) = dtKd(i,kb(i))
      enddo
    else
      do i=is,ie ; Ent_bl(i,Kmb+1) = 0.0 ; enddo
    endif

    do k=2,nz-1 ; do i=is,ie
      dsp1_ds(i,k) = 1.0 / ds_dsp1(i,k)
      I2p2dsp1_ds(i,k) = 0.5/(1.0+dsp1_ds(i,k))
      grats(i,k) = 2.0*(2.0+(dsp1_ds(i,k)+ds_dsp1(i,k)))
    enddo ; enddo

!   Determine the maximum flux, maxF, for each of the isopycnal layers.
!   Also determine when the fluxes start entraining
! from various buffer or mixed layers, where appropriate.
    if (CS%bulkmixedlayer) then
      kb_min = nz
      do i=is,ie
        htot(i) = h(i,j,1) - Angstrom
      enddo
      do k=2,kmb ; do i=is,ie
        htot(i) = htot(i) + (h(i,j,k) - Angstrom)
      enddo ; enddo
      do i=is,ie
        max_eakb(i) = MAX(Ent_bl(i,Kmb+1) + 0.5*htot(i), htot(i))
        I_dSkbp1(i) = 1.0 / (Sref(i,kmb+2) - Sref(i,kmb+1))
        zeros(i) = 0.0
      enddo

      !   Find the maximum amount of entrainment from below that the top
      ! interior layer could exhibit (maxF(i,kb)), approximating that
      ! entrainment by (eakb*max(dS_kb/dSkbp1,0)).  eakb is in the range
      ! from 0 to max_eakb.
      call find_maxF_kb(h_bl, Sref, Ent_bl, I_dSkbp1, zeros, max_eakb, kmb, &
                        is, ie, G, GV, CS, maxF_kb, eakb_maxF, do_i, F_kb_maxent)
      do i=is,ie ; if (kb(i) <= nz) then
        maxF(i,kb(i)) = MAX(0.0, maxF_kb(i), F_kb_maxent(i))
        if ((maxF_kb(i) > F_kb_maxent(i)) .and. (eakb_maxF(i) >= htot(i))) &
          max_eakb(i) = eakb_maxF(i)
      endif ; enddo

      do i=is,ie ; ea_kbp1(i) = 0.0 ; eakb(i) = max_eakb(i) ; enddo
      call determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, ea_kbp1, &
                           max_eakb, max_eakb, kmb, is, ie, do_i, G, GV, CS, eakb, &
                           error=err_max_eakb0, F_kb=F_kb)

      !   The maximum value of F(kb) between htot and max_eakb determines
      ! what maxF(kb+1) should be.
      do i=is,ie ; min_eakb(i) = MIN(htot(i), max_eakb(i)) ; enddo
      call find_maxF_kb(h_bl, Sref, Ent_bl, I_dSkbp1, min_eakb, max_eakb, &
                        kmb, is, ie, G, GV, CS, F_kb_maxEnt, do_i_in = do_i)

      do i=is,ie
        if ((.not.do_i(i)) .or. (err_max_eakb0(i) >= 0.0)) then
          eakb(i) = 0.0 ; min_eakb(i) = 0.0
        else ! If error_max_eakb0 < 0 the buffer layers are always all entrained.
          eakb(i) = max_eakb(i) ; min_eakb(i) = max_eakb(i)
        endif
      enddo

      !   Find the amount of entrainment of the buffer layers that would occur
      ! if there were no entrainment by the deeper interior layers.  Also find
      ! how much entrainment of the deeper layers would occur.
      call determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, ea_kbp1, &
                           zeros, max_eakb, kmb, is, ie, do_i, G, GV, CS, min_eakb, &
                           error=err_min_eakb0, F_kb=F_kb, err_max_eakb0=err_max_eakb0)
      ! Error_min_eakb0 should be ~0 unless error_max_eakb0 < 0.
      do i=is,ie ; if ((kb(i)<nz) .and. (kb_min>kb(i))) kb_min = kb(i) ; enddo
    else
      ! Without a bulk mixed layer, surface fluxes are applied in this
      ! subroutine.  (Otherwise, they are handled in mixedlayer.)
      !   Initially the maximum flux in layer zero is given by the surface
      ! buoyancy flux.  It will later be limited if the surface flux is
      ! too large.  Here buoy is the surface buoyancy flux.
      do i=is,ie
        maxF(i,1) = 0.0
        htot(i) = h(i,j,1) - Angstrom
      enddo
      if (ASSOCIATED(fluxes%buoy)) then ; do i=is,ie
        maxF(i,1) = (dt*fluxes%buoy(i,j)) / GV%g_prime(2)
      enddo ; endif
    endif

! The following code calculates the maximum flux, maxF, for the interior
! layers.
    do k=kb_min,nz-1 ; do i=is,ie
      if ((k == kb(i)+1) .and. CS%bulkmixedlayer) then
        maxF(i,k) = ds_dsp1(i,k)*(F_kb_maxEnt(i) + htot(i))
        htot(i) = htot(i) + (h(i,j,k) - Angstrom)
      elseif (k > kb(i)) then
        maxF(i,k) = ds_dsp1(i,k)*(maxF(i,k-1) + htot(i))
        htot(i) = htot(i) + (h(i,j,k) - Angstrom)
      endif
    enddo ; enddo
    do i=is,ie
      maxF(i,nz) = 0.0
      if (.not.CS%bulkmixedlayer) then
        maxF_correct(i) = MAX(0.0, -(maxF(i,nz-1) + htot(i)))
      endif
      htot(i) = h(i,j,nz) - Angstrom
    enddo
    if (.not.CS%bulkmixedlayer) then
      do_any = .false. ; do i=is,ie ; if (maxF_correct(i) > 0.0) do_any = .true. ; enddo
      if (do_any) then
        do k=nz-1,1,-1 ; do i=is,ie
          maxF(i,k) = maxF(i,k) + maxF_correct(i)
          maxF_correct(i) = maxF_correct(i) * dsp1_ds(i,k)
        enddo ; enddo
      endif
    endif
    do k=nz-1,kb_min,-1 ; do i=is,ie ; if (do_i(i)) then
      if (k>=kb(i)) then
        maxF(i,k) = MIN(maxF(i,k),dsp1_ds(i,k+1)*maxF(i,k+1) + htot(i))
        htot(i) = htot(i) + (h(i,j,k) - Angstrom)
        if ( (k == kb(i)) .and. ((maxF(i,k) < F_kb(i)) .or. &
            (maxF(i,k) < maxF_kb(i)) .and. (eakb_maxF(i) <= max_eakb(i))) ) then
          !   In this case, too much was being entrained by the topmost interior
          ! layer, even with the minimum initial estimate.  The buffer layer
          ! will always entrain the maximum amount.
          F_kb(i) = maxF(i,k)
          if ((F_kb(i) <= maxF_kb(i)) .and. (eakb_maxF(i) <= max_eakb(i))) then
            eakb(i) = eakb_maxF(i)
          else
            eakb(i) = max_eakb(i)
          endif
          call F_kb_to_ea_kb(h_bl, Sref, Ent_bl, I_dSkbp1, F_kb, kmb, i, &
                             G, GV, CS, eakb, Angstrom)
          if ((eakb(i) < max_eakb(i)) .or. (eakb(i) < min_eakb(i))) then
            call determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, zeros, &
                                 eakb, eakb, kmb, i, i, do_i, G, GV, CS, eakb, &
                                 error=err_eakb0)
            if (eakb(i) < max_eakb(i)) then
              max_eakb(i) = eakb(i) ; err_max_eakb0(i) = err_eakb0(i)
            endif
            if (eakb(i) < min_eakb(i)) then
              min_eakb(i) = eakb(i) ; err_min_eakb0(i) = err_eakb0(i)
            endif
          endif
        endif
      endif
    endif ; enddo ; enddo
    if (.not.CS%bulkmixedlayer) then
      do i=is,ie
        maxF(i,1) = MIN(maxF(i,1),dsp1_ds(i,2)*maxF(i,2) + htot(i))
      enddo
    endif

!   The following code provides an initial estimate of the flux in
! each layer, F.  The initial guess for the layer diffusive flux is
! the smaller of a forward discretization or the maximum diffusive
! flux starting from zero thickness in one time step without
! considering adjacent layers.
    do i=is,ie
      F(i,1) = maxF(i,1)
      F(i,nz) = maxF(i,nz) ; minF(i,nz) = 0.0
    enddo
    do k=nz-1,K2,-1
      do i=is,ie
        if ((k==kb(i)) .and. (do_i(i))) then
          eakb(i) = min_eakb(i)
          minF(i,k) = 0.0
        elseif ((k>kb(i)) .and. (do_i(i))) then
!   Here the layer flux is estimated, assuming no entrainment from
! the surrounding layers.  The estimate is a forward (steady) flux,
! limited by the maximum flux for a layer starting with zero
! thickness.  This is often a good guess and leads to few iterations.
          hm = h(i,j,k) + h_neglect
  !   Note: Tried sqrt((0.5*ds_dsp1(i,k))*dtKd(i,k)) for the second limit,
  ! but it usually doesn't work as well.
          F(i,k) = MIN(maxF(i,k), sqrt(ds_dsp1(i,k)*dtKd(i,k)), &
                       0.5*(ds_dsp1(i,k)+1.0) * (dtKd(i,k) / hm))

!   Calculate the minimum flux that can be expected if there is no entrainment
! from the neighboring layers.  The 0.9 is used to give used to give a 10%
! known error tolerance.
          fk = dtKd(i,k) * grats(i,k)
          minF(i,k) = MIN(maxF(i,k), &
                          0.9*(I2p2dsp1_ds(i,k) * fk / (hm + sqrt(hm*hm + fk))))
          if (k==kb(i)) minF(i,k) = 0.0 ! BACKWARD COMPATIBILITY - DELETE LATER?
        else
          F(i,k) = 0.0
          minF(i,k) = 0.0
        endif
      enddo ! end of i loop
    enddo ! end of k loop

    ! This is where the fluxes are actually calculated.

    is1 = ie+1 ; ie1 = is-1
    do i=is,ie ; if (do_i(i)) then ; is1 = i ; exit ; endif ; enddo
    do i=ie,is,-1 ; if (do_i(i)) then ; ie1 = i ; exit ; endif ; enddo

    if (CS%bulkmixedlayer) then
      kb_min_act = nz
      do i=is,ie
        if (do_i(i) .and. (kb(i) < kb_min_act)) kb_min_act = kb(i)
      enddo
      !   Solve for the entrainment rate from above in the topmost interior
      ! layer, eakb, such that
      !   eakb*dS_implicit = dt*Kd*dS_layer_implicit / h_implicit.
      do i=is1,ie1
        ea_kbp1(i) = 0.0
        if (do_i(i) .and. (kb(i) < nz)) &
          ea_kbp1(i) = dsp1_ds(i,kb(i)+1)*F(i,kb(i)+1)
      enddo
      call determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, ea_kbp1, min_eakb, &
                           max_eakb, kmb, is1, ie1, do_i, G, GV, CS, eakb, F_kb=F_kb, &
                           err_max_eakb0=err_max_eakb0, err_min_eakb0=err_min_eakb0, &
                           dFdfm_kb=dFdfm_kb)
    else
      kb_min_act = kb_min
    endif

    do it=0,CS%max_ent_it-1
      do i=is1,ie1 ; if (do_i(i)) then
        if (.not.CS%bulkmixedlayer) F(i,1) = MIN(F(i,1),maxF(i,1))
        b1(i) = 1.0
      endif ; enddo ! end of i loop

     ! F_kb has already been found for this iteration, either above or at
     ! the end of the code for the previous iteration.
      do k=kb_min_act,nz-1 ; do i=is1,ie1 ; if (do_i(i) .and. (k>=kb(i))) then
        ! Calculate the flux in layer k.
        if (CS%bulkmixedlayer .and. (k==kb(i))) then
          F(i,k) = F_kb(i)
          dFdfm(i,k) = dFdfm_kb(i)
        else ! k > kb(i)
          Fprev(i,k) = F(i,k)
          fm = (F(i,k-1) - h(i,j,k)) + dsp1_ds(i,k+1)*F(i,k+1)
          fk = grats(i,k)*dtKd(i,k)
          fr = sqrt(fm*fm + fk)

          if (fm>=0) then
            F(i,k) = MIN(maxF(i,k), I2p2dsp1_ds(i,k) * (fm+fr))
          else
            F(i,k) = MIN(maxF(i,k), I2p2dsp1_ds(i,k) * (fk / (-1.0*fm+fr)))
          endif

          if ((F(i,k) >= maxF(i,k)) .or. (fr == 0.0)) then
            dFdfm(i,k) = 0.0
          else
            dFdfm(i,k) = I2p2dsp1_ds(i,k) * ((fr + fm) / fr)
          endif

          if (k > K2) then
            ! This is part of a tridiagonal solver for the actual flux.
            c1(i,k) = dFdfm(i,k-1)*(dsp1_ds(i,k)*b1(i))
            b1(i) = 1.0 / (1.0 - c1(i,k)*dFdfm(i,k))
            F(i,k) = MIN(b1(i)*(F(i,k)-Fprev(i,k)) + Fprev(i,k), maxF(i,k))
            if (F(i,k) >= maxF(i,k)) dFdfm(i,k) = 0.0
          endif
        endif
      endif ; enddo ; enddo

      do k=nz-2,kb_min_act,-1 ; do i=is1,ie1
        if (do_i(i) .and. (k > kb(i))) &
          F(i,k) = MIN((F(i,k)+c1(i,k+1)*(F(i,k+1)-Fprev(i,k+1))),maxF(i,k))
      enddo ; enddo

      if (CS%bulkmixedlayer) then
        do i=is1,ie1
          if (do_i(i) .and. (kb(i) < nz)) then
            ! F will be increased to minF later.
            ea_kbp1(i) = dsp1_ds(i,kb(i)+1)*max(F(i,kb(i)+1), minF(i,kb(i)+1))
          else
            ea_kbp1(i) = 0.0
          endif
        enddo
        call determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, ea_kbp1, min_eakb, &
                             max_eakb, kmb, is1, ie1, do_i, G, GV, CS, eakb, F_kb=F_kb, &
                             err_max_eakb0=err_max_eakb0, err_min_eakb0=err_min_eakb0, &
                             dFdfm_kb=dFdfm_kb)
        do i=is1,ie1
          if (do_i(i) .and. (kb(i) < nz)) F(i,kb(i)) = F_kb(i)
        enddo
      endif

! Determine whether to do another iteration.
      if (it < CS%max_ent_it-1) then

        reiterate = .false.
        if (CS%bulkmixedlayer) then ; do i=is1,ie1 ; if (do_i(i)) then
          eb_kmb(i) = max(2.0*Ent_bl(i,Kmb+1) - eakb(i), 0.0)
        endif ; enddo ; endif
        do i=is1,ie1
          did_i(i) = do_i(i) ; do_i(i) = .false.
        enddo
        do k=kb_min_act,nz-1 ; do i=is1,ie1
          if (did_i(i) .and. (k >= kb(i))) then
            if (F(i,k) < minF(i,k)) then
              F(i,k) = minF(i,k)
              do_i(i) = .true. ; reiterate = .true.
            elseif (k > kb(i)) then
              if ((abs(F(i,k) - Fprev(i,k)) > tolerance) .or. &
                  ((h(i,j,k) + ((1.0+dsp1_ds(i,k))*F(i,k) - &
                   (F(i,k-1) + dsp1_ds(i,k+1)*F(i,k+1)))) < 0.9*Angstrom)) then
                do_i(i) = .true. ; reiterate = .true.
              endif
            else ! (k == kb(i))
! A more complicated test is required for the layer beneath the buffer layer,
! since its flux may be partially used to entrain directly from the mixed layer.
! Negative fluxes should not occur with the bulk mixed layer.
              if (h(i,j,k) + ((F(i,k) + eakb(i)) - &
                  (eb_kmb(i) + dsp1_ds(i,k+1)*F(i,k+1))) < 0.9*Angstrom) then
                do_i(i) = .true. ; reiterate = .true.
              endif
            endif
          endif
        enddo ; enddo
        if (.not.reiterate) exit
      endif ! end of if (it < CS%max_ent_it-1)
    enddo ! end of it loop
! This is the end of the section that might be iterated.


    if (it == (CS%max_ent_it)) then
      !   Limit the flux so that the layer below is not depleted.
      ! This should only be applied to the last iteration.
      do i=is1,ie1 ; if (do_i(i)) then
        F(i,nz-1) = MAX(F(i,nz-1), MIN(minF(i,nz-1), 0.0))
        if (kb(i) >= nz-1) then ; ea_kbp1(i) = 0.0 ; endif
      endif ; enddo
      do k=nz-2,kb_min_act,-1 ; do i=is1,ie1 ; if (do_i(i)) then
        if (k>kb(i)) then
          F(i,k) = MIN(MAX(minF(i,k),F(i,k)), (dsp1_ds(i,k+1)*F(i,k+1) + &
               MAX(((F(i,k+1)-dsp1_ds(i,k+2)*F(i,k+2)) + &
                    (h(i,j,k+1) - Angstrom)), 0.5*(h(i,j,k+1)-Angstrom))))
        elseif (k==kb(i)) then
          ea_kbp1(i) = dsp1_ds(i,k+1)*F(i,k+1)
          h_avail = dsp1_ds(i,k+1)*F(i,k+1) + MAX(0.5*(h(i,j,k+1)-Angstrom), &
               ((F(i,k+1)-dsp1_ds(i,k+2)*F(i,k+2)) + (h(i,j,k+1) - Angstrom)))
          if ((F(i,k) > 0.0) .and. (F(i,k) > h_avail)) then
            F_kb(i) = MAX(0.0, h_avail)
            F(i,k) = F_kb(i)
            if ((F_kb(i) < maxF_kb(i)) .and. (eakb_maxF(i) <= eakb(i))) &
              eakb(i) = eakb_maxF(i)
            call F_kb_to_ea_kb(h_bl, Sref, Ent_bl, I_dSkbp1, F_kb, kmb, i, &
                               G, GV, CS, eakb)
          endif
        endif
      endif ; enddo ; enddo


      if (CS%bulkmixedlayer) then ; do i=is1,ie1
        if (do_i(i) .and. (kb(i) < nz)) then
          h_avail = eakb(i) + MAX(0.5*(h_bl(i,kmb+1)-Angstrom), &
                (F_kb(i)-ea_kbp1(i)) + (h_bl(i,kmb+1)-Angstrom))
          ! Ensure that 0 < eb_kmb < h_avail.
          Ent_bl(i,Kmb+1) = MIN(Ent_bl(i,Kmb+1),0.5*(eakb(i) + h_avail))

          eb_kmb(i) = max(2.0*Ent_bl(i,Kmb+1) - eakb(i), 0.0)
        endif
      enddo ; endif

     !  Limit the flux so that the layer above is not depleted.
      do k=kb_min_act+1,nz-1 ; do i=is1,ie1 ; if (do_i(i)) then
        if ((.not.CS%bulkmixedlayer) .or. (k > kb(i)+1)) then
          F(i,k) = MIN(F(i,k), ds_dsp1(i,k)*( ((F(i,k-1) + &
              dsp1_ds(i,k-1)*F(i,k-1)) - F(i,k-2)) + (h(i,j,k-1) - Angstrom)))
          F(i,k) = MAX(F(i,k),MIN(minF(i,k),0.0))
        else if (k == kb(i)+1) then
          F(i,k) = MIN(F(i,k), ds_dsp1(i,k)*( ((F(i,k-1) + eakb(i)) - &
              eb_kmb(i)) + (h(i,j,k-1) - Angstrom)))
          F(i,k) = MAX(F(i,k),MIN(minF(i,k),0.0))
        endif
      endif ; enddo ; enddo
    endif ! (it == (CS%max_ent_it))

    call F_to_ent(F, h, kb, kmb, j, G, GV, CS, dsp1_ds, eakb, Ent_bl, ea, eb)

    !   Calculate the layer thicknesses after the entrainment to constrain the
    ! corrective fluxes.
    if (correct_density) then
      do i=is,ie
        h_guess(i,1) = (h(i,j,1) - Angstrom) + (eb(i,j,1) - ea(i,j,2))
        h_guess(i,nz) = (h(i,j,nz) - Angstrom) + (ea(i,j,nz) - eb(i,j,nz-1))
        if (h_guess(i,1) < 0.0) h_guess(i,1) = 0.0
        if (h_guess(i,nz) < 0.0) h_guess(i,nz) = 0.0
      enddo
      do k=2,nz-1 ; do i=is,ie
        h_guess(i,k) = (h(i,j,k) - Angstrom) + ((ea(i,j,k) - eb(i,j,k-1)) + &
                   (eb(i,j,k) - ea(i,j,k+1)))
        if (h_guess(i,k) < 0.0) h_guess(i,k) = 0.0
      enddo ; enddo
      if (CS%bulkmixedlayer) then
        call determine_dSkb(h_bl, Sref, Ent_bl, eakb, is, ie, kmb, G, GV, &
                            .true., dS_kb, dS_anom_lim=dS_anom_lim)
        do k=nz-1,kb_min,-1
          call calculate_density(tv%T(is:ie,j,k), tv%S(is:ie,j,k), pres(is:ie), &
                                 Rcv(is:ie), 1, ie-is+1, tv%eqn_of_state)
          do i=is,ie
            if ((k>kb(i)) .and. (F(i,k) > 0.0)) then
              ! Within a time step, a layer may entrain no more than its
              ! thickness for correction.  This limitation should apply
              ! extremely rarely, but precludes undesirable behavior.
              !  Note: Corrected a sign/logic error & factor of 2 error, and
              !    the layers tracked the target density better, mostly due to
              !    the factor of 2 error.
              F_cor = h(i,j,k) * MIN(1.0 , MAX(-ds_dsp1(i,k), &
                          (GV%Rlay(k) - Rcv(i)) / (GV%Rlay(k+1)-GV%Rlay(k))) )

              ! Ensure that (1) Entrainments are positive, (2) Corrections in
              ! a layer cannot deplete the layer itself (very generously), and
              ! (3) a layer can take no more than a quarter the mass of its
              ! neighbor.
              if (F_cor > 0.0) then
                F_cor = MIN(F_cor, 0.9*F(i,k), ds_dsp1(i,k)*0.5*h_guess(i,k), &
                            0.25*h_guess(i,k+1))
              else
                F_cor = -MIN(-F_cor, 0.9*F(i,k), 0.5*h_guess(i,k), &
                             0.25*ds_dsp1(i,k)*h_guess(i,k-1) )
              endif

              ea(i,j,k) = ea(i,j,k) - dsp1_ds(i,k)*F_cor
              eb(i,j,k) = eb(i,j,k) + F_cor
            else if ((k==kb(i)) .and. (F(i,k) > 0.0)) then
              !   Rho_cor is the density anomaly that needs to be corrected,
              ! taking into account that the true potential density of the
              ! deepest buffer layer is not exactly what is returned as dS_kb.
              dS_kb_eff = 2.0*dS_kb(i) - dS_anom_lim(i) ! Could be negative!!!
              Rho_cor = h(i,j,k) * (GV%Rlay(k)-Rcv(i)) + eakb(i)*dS_anom_lim(i)

              ! Ensure that  -.9*eakb < ea_cor < .9*eakb
              if (abs(Rho_cor) < abs(0.9*eakb(i)*dS_kb_eff)) then
                ea_cor = -Rho_cor / dS_kb_eff
              else
                ea_cor = sign(0.9*eakb(i),-Rho_cor*dS_kb_eff)
              endif

              if (ea_cor > 0.0) then
                ! Ensure that -F_cor < 0.5*h_guess
                ea_cor = MIN(ea_cor, 0.5*(max_eakb(i) - eakb(i)), &
                             0.5*h_guess(i,k) / (dS_kb(i) * I_dSkbp1(i)))
              else
                ! Ensure that -ea_cor < 0.5*h_guess & F_cor < 0.25*h_guess(k+1)
                ea_cor = -MIN(-ea_cor, 0.5*h_guess(i,k), &
                              0.25*h_guess(i,k+1) / (dS_kb(i) * I_dSkbp1(i)))
              endif

              ea(i,j,k) = ea(i,j,k) + ea_cor
              eb(i,j,k) = eb(i,j,k) - (dS_kb(i) * I_dSkbp1(i)) * ea_cor
            else if (k < kb(i)) then
              ! Repetative, unless ea(kb) has been corrected.
              ea(i,j,k) = ea(i,j,k+1)
            endif
          enddo
        enddo
        do k=kb_min-1,K2,-1 ; do i=is,ie
          ea(i,j,k) = ea(i,j,k+1)
        enddo ; enddo

        ! Repetative, unless ea(kb) has been corrected.
        k=kmb
        do i=is,ie
          ! Do not adjust eb through the base of the buffer layers, but it
          ! may be necessary to change entrainment from above.
          h1 = (h(i,j,k) - Angstrom) + (eb(i,j,k) - ea(i,j,k+1))
          ea(i,j,k) = MAX(Ent_bl(i,K), Ent_bl(i,K)-0.5*h1, -h1)
        enddo
        do k=kmb-1,2,-1 ; do i=is,ie
          ! Determine the entrainment from below for each buffer layer.
          eb(i,j,k) = max(2.0*Ent_bl(i,K+1) - ea(i,j,k+1), 0.0)

          ! Determine the entrainment from above for each buffer layer.
          h1 = (h(i,j,k) - Angstrom) + (eb(i,j,k) - ea(i,j,k+1))
          ea(i,j,k) = MAX(Ent_bl(i,K), Ent_bl(i,K)-0.5*h1, -h1)
        enddo ; enddo
        do i=is,ie
          eb(i,j,1) = max(2.0*Ent_bl(i,2) - ea(i,j,2), 0.0)
        enddo

      else ! not bulkmixedlayer
        do k=K2,nz-1
          call calculate_density(tv%T(is:ie,j,k), tv%S(is:ie,j,k), pres(is:ie), &
                                 Rcv(is:ie), 1, ie-is+1, tv%eqn_of_state)
          do i=is,ie ; if (F(i,k) > 0.0) then
            ! Within a time step, a layer may entrain no more than
            ! its thickness for correction.  This limitation should
            ! apply extremely rarely, but precludes undesirable
            ! behavior.
            F_cor = h(i,j,k) * MIN(dsp1_ds(i,k) , MAX(-1.0, &
                       (GV%Rlay(k) - Rcv(i)) / (GV%Rlay(k+1)-GV%Rlay(k))) )

            ! Ensure that (1) Entrainments are positive, (2) Corrections in
            ! a layer cannot deplete the layer itself (very generously), and
            ! (3) a layer can take no more than a quarter the mass of its
            ! neighbor.
            if (F_cor >= 0.0) then
              F_cor = MIN(F_cor, 0.9*F(i,k), 0.5*dsp1_ds(i,k)*h_guess(i,k), &
                          0.25*h_guess(i,k+1))
            else
              F_cor = -MIN(-F_cor, 0.9*F(i,k), 0.5*h_guess(i,k), &
                           0.25*ds_dsp1(i,k)*h_guess(i,k-1) )
            endif

            ea(i,j,k) = ea(i,j,k) - dsp1_ds(i,k)*F_cor
            eb(i,j,k) = eb(i,j,k) + F_cor
          endif ; enddo
        enddo
      endif

    endif   ! correct_density

    if (CS%id_Kd > 0) then
      Idt = 1.0 / dt
      do k=2,nz-1 ; do i=is,ie
        if (k<kb(i)) then ; Kd_here = 0.0 ; else
          Kd_here = F(i,k) * ( h(i,j,k) + ((ea(i,j,k) - eb(i,j,k-1)) + &
              (eb(i,j,k) - ea(i,j,k+1))) ) / (I2p2dsp1_ds(i,k) * grats(i,k))
        endif

        Kd_eff(i,j,k) = H_to_m**2 * (MAX(dtKd(i,k),Kd_here)*Idt)
      enddo ; enddo
      do i=is,ie
        Kd_eff(i,j,1) = H_to_m**2 * (dtKd(i,1)*Idt)
        Kd_eff(i,j,nz) = H_to_m**2 * (dtKd(i,nz)*Idt)
      enddo
    endif

    if (CS%id_diff_work > 0) then
      do i=is,ie ; diff_work(i,j,1) = 0.0 ; diff_work(i,j,nz+1) = 0.0 ; enddo
      if (associated(tv%eqn_of_state)) then
        if (associated(fluxes%p_surf)) then
          do i=is,ie ; pressure(i) = fluxes%p_surf(i,j) ; enddo
        else
          do i=is,ie ; pressure(i) = 0.0 ; enddo
        endif
        do K=2,nz
          do i=is,ie ; pressure(i) = pressure(i) + GV%H_to_Pa*h(i,j,k-1) ; enddo
          do i=is,ie
            if (k==kb(i)) then
              T_eos(i) = 0.5*(tv%T(i,j,kmb) + tv%T(i,j,k))
              S_eos(i) = 0.5*(tv%S(i,j,kmb) + tv%S(i,j,k))
            else
              T_eos(i) = 0.5*(tv%T(i,j,k-1) + tv%T(i,j,k))
              S_eos(i) = 0.5*(tv%S(i,j,k-1) + tv%S(i,j,k))
            endif
          enddo
          call calculate_density_derivs(T_eos, S_eos, pressure, &
                  dRho_dT, dRho_dS, is, ie-is+1, tv%eqn_of_state)
          do i=is,ie
            if ((k>kmb) .and. (k<kb(i))) then ; diff_work(i,j,K) = 0.0
            else
              if (k==kb(i)) then
                dRho = dRho_dT(i) * (tv%T(i,j,k)-tv%T(i,j,kmb)) + &
                       dRho_dS(i) * (tv%S(i,j,k)-tv%S(i,j,kmb))
              else
                dRho = dRho_dT(i) * (tv%T(i,j,k)-tv%T(i,j,k-1)) + &
                       dRho_dS(i) * (tv%S(i,j,k)-tv%S(i,j,k-1))
              endif
              diff_work(i,j,K) = g_2dt * dRho * &
                   (ea(i,j,k) * (h(i,j,k) + ea(i,j,k)) + &
                    eb(i,j,k-1)*(h(i,j,k-1) + eb(i,j,k-1)))
            endif
          enddo
        enddo
      else
        do K=2,nz ; do i=is,ie
          diff_work(i,j,K) = g_2dt * (GV%Rlay(k)-GV%Rlay(k-1)) * &
               (ea(i,j,k) * (h(i,j,k) + ea(i,j,k)) + &
                eb(i,j,k-1)*(h(i,j,k-1) + eb(i,j,k-1)))
        enddo ; enddo
      endif
    endif

    if (present(kb_out)) then
      do i=is,ie ; kb_out(i,j) = kb(i) ; enddo
    endif

  enddo ! end of j loop

! Offer diagnostic fields for averaging.
  if (CS%id_Kd > 0) call post_data(CS%id_Kd, Kd_eff, CS%diag)
  if (CS%id_Kd > 0) deallocate(Kd_eff)
  if (CS%id_diff_work > 0) call post_data(CS%id_diff_work, diff_work, CS%diag)
  if (CS%id_diff_work > 0) deallocate(diff_work)

end subroutine entrainment_diffusive


subroutine F_to_ent(F, h, kb, kmb, j, G, GV, CS, dsp1_ds, eakb, Ent_bl, ea, eb, do_i_in)
  type(ocean_grid_type),                    intent(in)    :: G
  type(verticalGrid_type),                  intent(in)    :: GV
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: F
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h
  integer, dimension(SZI_(G)),              intent(in)    :: kb
  integer,                                  intent(in)    :: kmb, j
  type(entrain_diffusive_CS),               intent(in)    :: CS
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: dsp1_ds
  real, dimension(SZI_(G)),                 intent(in)    :: eakb
  real, dimension(SZI_(G),SZK_(G)),         intent(in)    :: Ent_bl
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(inout) :: ea, eb
  logical, dimension(SZI_(G)),    optional, intent(in)    :: do_i_in
!   This subroutine calculates the actual entrainments (ea and eb) and the
! amount of surface forcing that is applied to each layer if there is no bulk
! mixed layer.

  real :: h1        ! The thickness in excess of the minimum that will remain
                    ! after exchange with the layer below, in m or kg m-2.
  logical :: do_i(SZI_(G))
  integer :: i, k, is, ie, nz

  is = G%isc ; ie = G%iec ; nz = G%ke

  if (present(do_i_in)) then
    do i=is,ie ; do_i(i) = do_i_in(i) ; enddo
    do i=G%isc,G%iec ; if (do_i(i)) then
      is = i ; exit
    endif ; enddo
    do i=G%iec,G%isc,-1 ; if (do_i(i)) then
      ie = i ; exit
    endif ; enddo
  else
    do i=is,ie ; do_i(i) = .true. ; enddo
  endif

  do i=is,ie
    ea(i,j,nz) = 0.0 ; eb(i,j,nz) = 0.0
  enddo
  if (CS%bulkmixedlayer) then
    do i=is,ie
      eb(i,j,kmb) = max(2.0*Ent_bl(i,Kmb+1) - eakb(i), 0.0)
    enddo
    do k=nz-1,kmb+1,-1 ; do i=is,ie ; if (do_i(i)) then
      if (k > kb(i)) then
        ! With a bulk mixed layer, surface buoyancy fluxes are applied
        ! elsewhere, so F should always be nonnegative.
        ea(i,j,k) = dsp1_ds(i,k)*F(i,k)
        eb(i,j,k) = F(i,k)
      else if (k == kb(i)) then
        ea(i,j,k) = eakb(i)
        eb(i,j,k) = F(i,k)
      elseif (k == kb(i)-1) then
        ea(i,j,k) = ea(i,j,k+1)
        eb(i,j,k) = eb(i,j,kmb)
      else
        ea(i,j,k) = ea(i,j,k+1)
        !   Add the entrainment of the thin interior layers to eb going
        ! up into the buffer layer.
        eb(i,j,k) = eb(i,j,k+1) + max(0.0, h(i,j,k+1) - GV%Angstrom)
      endif
    endif ; enddo ; enddo
    k = kmb
    do i=is,ie ; if (do_i(i)) then
      ! Adjust the previously calculated entrainment from below by the deepest
      ! buffer layer to account for entrainment of thin interior layers .
      if (kb(i) > kmb+1) &
        eb(i,j,k) = eb(i,j,k+1) + max(0.0, h(i,j,k+1) - GV%Angstrom)

      ! Determine the entrainment from above for each buffer layer.
      h1 = (h(i,j,k) - GV%Angstrom) + (eb(i,j,k) - ea(i,j,k+1))
      ea(i,j,k) = MAX(Ent_bl(i,K), Ent_bl(i,K)-0.5*h1, -h1)
    endif ; enddo
    do k=kmb-1,2,-1 ; do i=is,ie ; if (do_i(i)) then
      ! Determine the entrainment from below for each buffer layer.
      eb(i,j,k) = max(2.0*Ent_bl(i,K+1) - ea(i,j,k+1), 0.0)

      ! Determine the entrainment from above for each buffer layer.
      h1 = (h(i,j,k) - GV%Angstrom) + (eb(i,j,k) - ea(i,j,k+1))
      ea(i,j,k) = MAX(Ent_bl(i,K), Ent_bl(i,K)-0.5*h1, -h1)
!       if (h1 >= 0.0) then ;                     ea(i,j,k) = Ent_bl(i,K)
!       elseif (Ent_bl(i,K)+0.5*h1 >= 0.0) then ; ea(i,j,k) = Ent_bl(i,K)-0.5*h1
!       else ;                                    ea(i,j,k) = -h1 ; endif
    endif ; enddo ; enddo
    do i=is,ie ; if (do_i(i)) then
      eb(i,j,1) = max(2.0*Ent_bl(i,2) - ea(i,j,2), 0.0)
      ea(i,j,1) = 0.0
    endif ; enddo
  else                                          ! not BULKMIXEDLAYER
    ! Calculate the entrainment by each layer from above and below.
    ! Entrainment is always positive, but F may be negative due to
    ! surface buoyancy fluxes.
    do i=is,ie
      ea(i,j,1) = 0.0 ; eb(i,j,1) = MAX(F(i,1),0.0)
      ea(i,j,2) = dsp1_ds(i,2)*F(i,2) - MIN(F(i,1),0.0)
    enddo

    do k=2,nz-1 ; do i=is,ie
      eb(i,j,k) = MAX(F(i,k),0.0)
      ea(i,j,k+1) = dsp1_ds(i,k+1)*F(i,k+1) - (F(i,k)-eb(i,j,k))
      if (ea(i,j,k+1) < 0.0) then
        eb(i,j,k) = eb(i,j,k) - ea(i,j,k+1)
        ea(i,j,k+1) = 0.0
      endif
    enddo ; enddo
  endif                                         ! end BULKMIXEDLAYER
end subroutine F_to_ent

subroutine set_Ent_bl(h, dtKd_int, tv, kb, kmb, do_i, G, GV, CS, j, Ent_bl, Sref, h_bl)
  type(ocean_grid_type),                    intent(in)  :: G
  type(verticalGrid_type),                  intent(in)  :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)  :: h
  real, dimension(SZI_(G),SZK_(G)+1),       intent(in)  :: dtKd_int
  type(thermo_var_ptrs),                    intent(in)  :: tv
  integer, dimension(SZI_(G)),              intent(inout) :: kb
  integer,                                  intent(in)  :: kmb
  logical, dimension(SZI_(G)),              intent(in)  :: do_i
  type(entrain_diffusive_CS),               pointer     :: CS
  integer,                                  intent(in)  :: j
  real, dimension(SZI_(G),SZK_(G)+1),       intent(out) :: Ent_bl
  real, dimension(SZI_(G),SZK_(G)),         intent(out) :: Sref, h_bl
! Arguments: h - Layer thickness, in m or kg m-2 (abbreviated as H below).
!  (in)      dtKd_int - The diapycnal diffusivity across each interface times
!                       the time step, in H2.
!  (in)      tv - A structure containing pointers to any available
!                 thermodynamic fields. Absent fields have NULL ptrs.
!  (in)      kb - The index of the lightest layer denser than the
!                 buffer layer or 1 if there is no buffer layer.
!  (in)      do_i - A logical variable indicating which i-points to work on.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - This module's control structure.
!  (in)      j - The meridional index upon which to work.
!  (out)     Ent_bl - The average entrainment upward and downward across
!                     each interface around the buffer layers, in H.
!  (out)     Sref - The coordinate potential density - 1000 for each layer,
!                   in kg m-3.
!  (out)     h_bl - The thickness of each layer, in H.

!   This subroutine sets the average entrainment across each of the interfaces
! between buffer layers within a timestep. It also causes thin and relatively
! light interior layers to be entrained by the deepest buffer layer.
!   Also find the initial coordinate potential densities (Sref) of each layer.

! Does there need to be limiting when the layers below are all thin?
  real, dimension(SZI_(G)) :: &
    b1, d1, &   ! Variables used by the tridiagonal solver, in H-1 and ND.
    Rcv, &      ! Value of the coordinate variable (potential density)
                ! based on the simulated T and S and P_Ref, kg m-3.
    pres, &     ! Reference pressure (P_Ref) in Pa.
    frac_rem, & ! The fraction of the diffusion remaining, ND.
    h_interior  ! The interior thickness available for entrainment, in H.
  real, dimension(SZI_(G), SZK_(G)) :: &
    S_est ! An estimate of the coordinate potential density - 1000 after
          ! entrainment for each layer, in kg m-3.
  real :: max_ent  ! The maximum possible entrainment, in H.
  real :: dh       ! An available thickness, in H.
  real :: Kd_x_dt  ! The diffusion that remains after thin layers are
                   ! entrained, in H2.
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  integer :: i, k, is, ie, nz
  is = G%isc ; ie = G%iec ; nz = G%ke

!  max_ent = 1.0e14*GV%Angstrom ! This is set to avoid roundoff problems.
  max_ent = 1.0e4*GV%m_to_H
  h_neglect = GV%H_subroundoff

  do i=is,ie ; pres(i) = tv%P_Ref ; enddo
  do k=1,kmb
    call calculate_density(tv%T(is:ie,j,k), tv%S(is:ie,j,k), pres(is:ie), &
                           Rcv(is:ie), 1, ie-is+1, tv%eqn_of_state)
    do i=is,ie
      h_bl(i,k) = h(i,j,k) + h_neglect
      Sref(i,k) = Rcv(i) - 1000.0
    enddo
  enddo

  do i=is,ie
    h_interior(i) = 0.0 ; Ent_bl(i,1) = 0.0
!     if (kb(i) > nz) Ent_bl(i,Kmb+1) = 0.0
  enddo

  do k=2,kmb ; do i=is,ie
    if (do_i(i)) then
      Ent_bl(i,K) = min(2.0 * dtKd_int(i,K) / (h(i,j,k-1) + h(i,j,k) + h_neglect), &
                        max_ent)
    else ; Ent_bl(i,K) = 0.0 ; endif
  enddo ; enddo

  !   Determine the coordinate density of the bottommost buffer layer if there
  ! is no entrainment from the layers below.  This is a partial solver, based
  ! on the first pass of a tridiagonal solver, as the values in the upper buffer
  ! layers are not needed.

  do i=is,ie
    b1(i) = 1.0 / (h_bl(i,1) + Ent_bl(i,2))
    d1(i) = h_bl(i,1) * b1(i)  ! = 1.0 - Ent_bl(i,2)*b1(i)
    S_est(i,1) = (h_bl(i,1)*Sref(i,1)) * b1(i)
  enddo
  do k=2,kmb-1 ; do i=is,ie
    b1(i) = 1.0 / ((h_bl(i,k) + Ent_bl(i,K+1)) + d1(i)*Ent_bl(i,K))
    d1(i) = (h_bl(i,k) + d1(i)*Ent_bl(i,K)) * b1(i)  ! = 1.0 - Ent_bl(i,K+1)*b1(i)
    S_est(i,k) = (h_bl(i,k)*Sref(i,k) + Ent_bl(i,K)*S_est(i,k-1)) * b1(i)
  enddo ; enddo
  do i=is,ie
    S_est(i,kmb) = (h_bl(i,kmb)*Sref(i,kmb) + Ent_bl(i,Kmb)*S_est(i,kmb-1)) / &
                   (h_bl(i,kmb) + d1(i)*Ent_bl(i,Kmb))
    frac_rem(i) = 1.0
  enddo

  !   Entrain any thin interior layers that are lighter (in the coordinate
  ! potential density) than the deepest buffer layer will be, and adjust kb.
  do i=is,ie ; kb(i) = nz+1 ; if (do_i(i)) kb(i) = kmb+1 ; enddo

  do k=kmb+1,nz ; do i=is,ie ; if (do_i(i)) then
    if ((k == kb(i)) .and. (S_est(i,kmb) > (GV%Rlay(k) - 1000.0))) then
      if (4.0*dtKd_int(i,Kmb+1)*frac_rem(i) > &
          (h_bl(i,kmb) + h(i,j,k)) * (h(i,j,k) - GV%Angstrom)) then
        ! Entrain this layer into the buffer layer and move kb down.
        dh = max((h(i,j,k) - GV%Angstrom), 0.0)
        if (dh > 0.0) then
          frac_rem(i) = frac_rem(i) - ((h_bl(i,kmb) + h(i,j,k)) * dh) / &
                                       (4.0*dtKd_int(i,Kmb+1))
          Sref(i,kmb) = (h_bl(i,kmb)*Sref(i,kmb) + dh*(GV%Rlay(k)-1000.0)) / &
                        (h_bl(i,kmb) + dh)
          h_bl(i,kmb) = h_bl(i,kmb) + dh
          S_est(i,kmb) = (h_bl(i,kmb)*Sref(i,kmb) + Ent_bl(i,Kmb)*S_est(i,kmb-1)) / &
                         (h_bl(i,kmb) + d1(i)*Ent_bl(i,Kmb))
        endif
        kb(i) = kb(i) + 1
      endif
    endif
  endif ; enddo ; enddo

 !    This is where variables are be set up with a different vertical grid
 !  in which the (newly?) massless layers are taken out.
  do k=nz,kmb+1,-1 ; do i=is,ie
    if (k >= kb(i)) h_interior(i) = h_interior(i) + (h(i,j,k)-GV%Angstrom)
    if (k==kb(i)) then
      h_bl(i,kmb+1) = h(i,j,k) ; Sref(i,kmb+1) = GV%Rlay(k) - 1000.0
    elseif (k==kb(i)+1) then
      h_bl(i,kmb+2) = h(i,j,k) ; Sref(i,kmb+2) = GV%Rlay(k) - 1000.0
    endif
  enddo ; enddo
  do i=is,ie ; if (kb(i) >= nz) then
    h_bl(i,kmb+1) = h(i,j,nz)
    Sref(i,kmb+1) = GV%Rlay(nz) - 1000.0
    h_bl(i,kmb+2) = GV%Angstrom
    Sref(i,kmb+2) = Sref(i,kmb+1) + (GV%Rlay(nz) - GV%Rlay(nz-1))
  endif ; enddo

  !   Perhaps we should revisit the way that the average entrainment between the
  ! buffer layer and the interior is calculated so that it is not unduly
  ! limited when the layers are less than sqrt(Kd * dt) thick?
  do i=is,ie ; if (do_i(i)) then
    Kd_x_dt = frac_rem(i) * dtKd_int(i,Kmb+1)
    if ((Kd_x_dt <= 0.0) .or. (h_interior(i) <= 0.0)) then
      Ent_bl(i,Kmb+1) = 0.0
    else
      !   If the combined layers are exceptionally thin, use sqrt(Kd*dt) as the
      ! estimate of the thickness in the denominator of the thickness diffusion.
      Ent_bl(i,Kmb+1) = MIN(0.5*h_interior(i), sqrt(Kd_x_dt), &
                            Kd_x_dt / (0.5*(h_bl(i,kmb) + h_bl(i,kmb+1))))
    endif
  else
    Ent_bl(i,Kmb+1) = 0.0
  endif ; enddo

end subroutine set_Ent_bl

subroutine determine_dSkb(h_bl, Sref, Ent_bl, E_kb, is, ie, kmb, G, GV, limit, &
                          dSkb, ddSkb_dE, dSlay, ddSlay_dE, dS_anom_lim, do_i_in)
  type(ocean_grid_type),            intent(in)    :: G
  type(verticalGrid_type),          intent(in)    :: GV
  real, dimension(SZI_(G),SZK_(G)), intent(in)    :: h_bl, Sref, Ent_bl
  real, dimension(SZI_(G)),         intent(in)    :: E_kb
  integer,                          intent(in)    :: is, ie, kmb
  logical,                          intent(in)    :: limit
  real, dimension(SZI_(G)),         intent(inout) :: dSkb
  real, dimension(SZI_(G)), optional, intent(inout) :: ddSkb_dE, dSlay, ddSlay_dE
  real, dimension(SZI_(G)), optional, intent(inout) :: dS_anom_lim
  logical, dimension(SZI_(G)), optional, intent(in) :: do_i_in
! Arguments: h_bl - Layer thickness, in m or kg m-2 (abbreviated as H below).
!  (in)      Sref - Reference potential vorticity (in kg m-3?)
!  (in)      Ent_bl - The average entrainment upward and downward across
!                     each interface around the buffer layers, in H.
!  (in)      E_kb - The entrainment by the top interior layer, in H.
!  (in)      is, ie - The range of i-indices to work on.
!  (in)      kmb - The number of mixed and buffer layers.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      limit - If true, limit dSkb and dSlay to avoid negative values.
!  (out)     dSkb - The limited potential density difference across the
!                   interface between the bottommost buffer layer and the
!                   topmost interior layer.  dSkb > 0.
!  (out,opt) dSlay - The limited potential density difference across the
!                        topmost interior layer. 0 < dSkb
!  (out,opt) ddSkb_dE - The partial derivative of dSkb with E, in kg m-3 H-1.
!  (out,opt) ddSlay_dE - The partial derivative of dSlay with E, in kg m-3 H-1.
!  (in,opt)  do_i_in - If present, determines which columns are worked on.
! Note that dSkb, ddSkb_dE, dSlay, ddSlay_dE, and dS_anom_lim are declared
! intent(inout) because they should not change where do_i_in is false.

!   This subroutine determines the reference density difference between the
! bottommost buffer layer and the first interior after the mixing between mixed
! and buffer layers and mixing with the layer below. Within the mixed and buffer
! layers, entrainment from the layer above is increased when it is necessary to
! keep the layers from developing a negative thickness; otherwise it equals
! Ent_bl.  At each interface, the upward and downward fluxes average out to
! Ent_bl, unless entrainment by the layer below is larger than twice Ent_bl.
!   The density difference across the first interior layer may also be returned.
! It could also be limited to avoid negative values or values that greatly
! exceed the density differences across an interface.
!   Additionally, the partial derivatives of dSkb and dSlay with E_kb could
! also be returned.
  real, dimension(SZI_(G),SZK_(G)) :: &
    b1, c1, &       ! b1 and c1 are variables used by the tridiagonal solver.
    S, dS_dE, &     ! The coordinate density and its derivative with R.
    ea, dea_dE, &   ! The entrainment from above and its derivative with R.
    eb, deb_dE      ! The entrainment from below and its derivative with R.
  real :: deriv_dSkb(SZI_(G))
  real :: d1(SZI_(G))  ! d1 = 1.0-c1 is also used by the tridiagonal solver.
  real :: src       ! A source term for dS_dR.
  real :: h1        ! The thickness in excess of the minimum that will remain
                    ! after exchange with the layer below, in m or kg m-2.
  logical, dimension(SZI_(G)) :: do_i
  real :: h_neglect ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected, in H.
  real :: h_tr      ! h_tr is h at tracer points with a tiny thickness
                    ! added to ensure positive definiteness, in m or kg m-2.
  real :: b_denom_1 ! The first term in the denominator of b1 in m or kg m-2.
  real :: rat
  real :: dS_kbp1, IdS_kbp1
  real :: deriv_dSLay
  real :: Inv_term     ! Nondimensional.
  real :: f1, df1_drat ! Nondimensional temporary variables.
  real :: z, dz_drat, f2, df2_dz, expz ! Nondimensional temporary variables.
  real :: eps_dSLay, eps_dSkb ! Small nondimensional constants.
  integer :: i, k

  if (present(ddSlay_dE) .and. .not.present(dSlay)) call MOM_error(FATAL, &
      "In deterimine_dSkb, ddSLay_dE may only be present if dSlay is.")

  h_neglect = GV%H_subroundoff

  do i=is,ie
    ea(i,kmb+1) = E_kb(i) ; dea_dE(i,kmb+1) = 1.0
    S(i,kmb+1) = Sref(i,kmb+1) ; dS_dE(i,kmb+1) = 0.0
    b1(i,kmb+1) = 0.0
    d1(i) = 1.0
    do_i(i) = .true.
  enddo
  if (present(do_i_in)) then
    do i=is,ie ; do_i(i) = do_i_in(i) ; enddo
  endif
  do k=kmb,1,-1 ; do i=is,ie
    if (do_i(i)) then
      ! The do_i test here is only for efficiency.
    ! Determine the entrainment from below for each buffer layer.
      if (2.0*Ent_bl(i,K+1) > ea(i,k+1)) then
        eb(i,k) = 2.0*Ent_bl(i,K+1) - ea(i,k+1) ; deb_dE(i,k) = -dea_dE(i,k+1)
      else
        eb(i,k) = 0.0 ; deb_dE(i,k) = 0.0
      endif

      ! Determine the entrainment from above for each buffer layer.
      h1 = (h_bl(i,k) - GV%Angstrom) + (eb(i,k) - ea(i,k+1))
      if (h1 >= 0.0) then
        ea(i,k) = Ent_bl(i,K) ; dea_dE(i,k) = 0.0
      elseif (Ent_bl(i,K) + 0.5*h1 >= 0.0) then
        ea(i,k) = Ent_bl(i,K) - 0.5*h1
        dea_dE(i,k) = 0.5*(dea_dE(i,k+1) - deb_dE(i,k))
      else
        ea(i,k) = -h1
        dea_dE(i,k) = dea_dE(i,k+1) - deb_dE(i,k)
      endif
    else
      ea(i,k) = 0.0 ; dea_dE(i,k) = 0.0 ; eb(i,k) = 0.0 ; deb_dE(i,k) = 0.0
    endif

    ! This is the first-pass of a tridiagonal solver for S.
    h_tr = h_bl(i,k) + h_neglect
    c1(i,k) = ea(i,k+1) * b1(i,k+1)
    b_denom_1 = (h_tr + d1(i)*eb(i,k))
    b1(i,k) = 1.0 / (b_denom_1 + ea(i,k))
    d1(i) = b_denom_1 * b1(i,k)

    S(i,k) = (h_tr*Sref(i,k) + eb(i,k)*S(i,k+1)) * b1(i,k)
  enddo ; enddo
  do k=2,kmb ; do i=is,ie
    S(i,k) = S(i,k) + c1(i,k-1)*S(i,k-1)
  enddo ; enddo

  if (present(ddSkb_dE) .or. present(ddSlay_dE)) then
    ! These two tridiagonal solvers cannot be combined because the solutions for
    ! S are required as a source for dS_dE.
    do k=kmb,2,-1 ; do i=is,ie
      if (do_i(i) .and. (dea_dE(i,k) - deb_dE(i,k) > 0.0)) then
        src = (((S(i,k+1) - Sref(i,k)) * (h_bl(i,k) + h_neglect) + &
                (S(i,k+1) - S(i,k-1)) * ea(i,k)) * deb_dE(i,k) - &
               ((Sref(i,k) - S(i,k-1)) * h_bl(i,k) + &
                (S(i,k+1) - S(i,k-1)) * eb(i,k)) * dea_dE(i,k)) / &
              ((h_bl(i,k) + h_neglect + ea(i,k)) + eb(i,k))
      else ; src = 0.0 ; endif
      dS_dE(i,k) = (src + eb(i,k)*dS_dE(i,k+1)) * b1(i,k)
    enddo ; enddo
    do i=is,ie
      if (do_i(i) .and. (deb_dE(i,1) < 0.0)) then
        src = (((S(i,2) - Sref(i,1)) * (h_bl(i,1) + h_neglect)) * deb_dE(i,1)) / &
              (h_bl(i,1) + h_neglect + eb(i,1))
      else ; src = 0.0 ; endif
      dS_dE(i,1) = (src + eb(i,1)*dS_dE(i,2)) * b1(i,1)
    enddo
    do k=2,kmb ; do i=is,ie
      dS_dE(i,k) = dS_dE(i,k) + c1(i,k-1)*dS_dE(i,k-1)
    enddo ; enddo
  endif

  ! Now, apply any limiting and return the requested variables.

  eps_dSkb = 1.0e-6   ! Should be a small, nondimensional, positive number.
  if (.not.limit) then
    do i=is,ie ; if (do_i(i)) then
      dSkb(i) = Sref(i,kmb+1) - S(i,kmb)
    endif ; enddo
    if (present(ddSkb_dE)) then ; do i=is,ie ; if (do_i(i)) then
      ddSkb_dE(i) = -1.0*dS_dE(i,kmb)
    endif ; enddo ; endif

    if (present(dSlay)) then ; do i=is,ie ; if (do_i(i)) then
      dSlay(i) = 0.5 * (Sref(i,kmb+2) - S(i,kmb))
    endif ; enddo ; endif
    if (present(ddSlay_dE)) then ; do i=is,ie ; if (do_i(i)) then
      ddSlay_dE(i) = -0.5*dS_dE(i,kmb)
    endif ; enddo ; endif
  else
    do i=is,ie ; if (do_i(i)) then
      ! Need to ensure that 0 < dSkb <= S_kb - Sbl
      if (Sref(i,kmb+1) - S(i,kmb) < eps_dSkb*(Sref(i,kmb+2) - Sref(i,kmb+1))) then
        dSkb(i) = eps_dSkb * (Sref(i,kmb+2) - Sref(i,kmb+1)) ; deriv_dSkb(i) = 0.0
      else
        dSkb(i) = Sref(i,kmb+1) - S(i,kmb) ; deriv_dSkb(i) = -1.0
      endif
      if (present(ddSkb_dE)) ddSkb_dE(i) = deriv_dSkb(i)*dS_dE(i,kmb)
    endif ; enddo

    if (present(dSLay)) then
      dz_drat = 1000.0    ! The limit of large dz_drat the same as choosing a
                          ! Heaviside function.
      eps_dSLay = 1.0e-10 ! Should be ~= GV%Angstrom / sqrt(Kd*dt)
      do i=is,ie ; if (do_i(i)) then
        dS_kbp1 = Sref(i,kmb+2) - Sref(i,kmb+1)
        IdS_kbp1 = 1.0 / (Sref(i,kmb+2) - Sref(i,kmb+1))
        rat = (Sref(i,kmb+1) - S(i,kmb)) * IdS_kbp1
        ! Need to ensure that 0 < dSLay <= 2*dSkb
        if (rat < 0.5) then
          ! The coefficients here are chosen so that at rat = 0.5, the value (1.5)
          ! and first derivative (-0.5) match with the "typical" case (next).
          ! The functional form here is arbitrary.
          !   f1 provides a reasonable profile that matches the value and derivative
          ! of the "typical" case at rat = 0.5, and has a maximum of less than 2.
          Inv_term = 1.0 / (1.0-rat)
          f1 = 2.0 - 0.125*(Inv_term**2)
          df1_drat = - 0.25*(Inv_term**3)

          !   f2 ensures that dSLay goes to 0 rapidly if rat is significantly
          ! negative.
          z = dz_drat * rat + 4.0 ! The 4 here gives f2(0) = 0.982.
          if (z >= 18.0) then ; f2 = 1.0 ; df2_dz = 0.0
          elseif (z <= -58.0) then ; f2 = eps_dSLay ; df2_dz = 0.0
          else
            expz = exp(z) ; Inv_term = 1.0 / (1.0 + expz)
            f2 = (eps_dSLay + expz) * Inv_term
            df2_dz = (1.0 - eps_dSLay) * expz * Inv_term**2
          endif

          dSLay(i) = dSkb(i) * f1 * f2
          deriv_dSLay = deriv_dSkb(i) * (f1 * f2) - (dSkb(i)*IdS_kbp1) * &
                            (df1_drat*f2 + f1 * dz_drat * df2_dz)
        elseif (dSkb(i) <= 3.0*dS_kbp1) then
          ! This is the "typical" case.
          dSLay(i) = 0.5 * (dSkb(i) + dS_kbp1)
          deriv_dSLay = 0.5 * deriv_dSkb(i) ! = -0.5
        else
          dSLay(i) = 2.0*dS_kbp1
          deriv_dSLay = 0.0
        endif
        if (present(ddSlay_dE)) ddSlay_dE(i) = deriv_dSLay*dS_dE(i,kmb)
      endif ; enddo
    endif ! present(dSlay)
  endif ! Not limited.

  if (present(dS_anom_lim)) then ; do i=is,ie ; if (do_i(i)) then
    dS_anom_lim(i) = max(0.0, eps_dSkb * (Sref(i,kmb+2) - Sref(i,kmb+1)) - &
                              (Sref(i,kmb+1) - S(i,kmb)) )
  endif ; enddo ; endif

end subroutine determine_dSkb


subroutine F_kb_to_ea_kb(h_bl, Sref, Ent_bl, I_dSkbp1, F_kb, kmb, i, &
                         G, GV, CS, ea_kb, tol_in)
  type(ocean_grid_type),         intent(in)    :: G
  type(verticalGrid_type),       intent(in)    :: GV
  real, dimension(SZI_(G),SZK_(G)), intent(in) :: h_bl, Sref, Ent_bl
  real, dimension(SZI_(G)),      intent(in)    :: I_dSkbp1, F_kb
  integer,                       intent(in)    :: kmb, i
  type(entrain_diffusive_CS),    pointer       :: CS
  real, dimension(SZI_(G)),      intent(inout) :: ea_kb
  real,                optional, intent(in)    :: tol_in

  !   Given an entrainment from below for layer kb, determine a consistent
  ! entrainment from above, such that dSkb * ea_kb = dSkbp1 * F_kb.  The input
  ! value of ea_kb is both the maximum value that can be obtained and the first
  ! guess of the iterations.  Also, make sure that ea_kb is an under-estimate
  real :: max_ea, min_ea
  real :: err, err_min, err_max
  real :: derr_dea
  real :: val, tolerance, tol1
  real :: ea_prev
  real :: dS_kbp1
  logical :: bisect_next, Newton
  real, dimension(SZI_(G)) :: dS_kb
  real, dimension(SZI_(G)) :: maxF, ent_maxF, zeros
  real, dimension(SZI_(G)) :: ddSkb_dE
  integer :: it
  integer, parameter :: MAXIT = 30

  dS_kbp1 = Sref(i,kmb+2) - Sref(i,kmb+1)
  max_ea = ea_kb(i) ; min_ea = 0.0
  val = dS_kbp1 * F_kb(i)
  err_min = -val

  tolerance = GV%m_to_H * CS%Tolerance_Ent
  if (present(tol_in)) tolerance = tol_in
  bisect_next = .true.

  call determine_dSkb(h_bl, Sref, Ent_bl, ea_kb, i, i, kmb, G, GV, .true., &
                      dS_kb, ddSkb_dE)

  err = dS_kb(i) * ea_kb(i) - val
  derr_dea = dS_kb(i) + ddSkb_dE(i) * ea_kb(i)
  ! Return if Newton's method on the first guess would give a tolerably small
  ! change in the value of ea_kb.
  if ((err <= 0.0) .and. (abs(err) <= tolerance*abs(derr_dea))) return

  if (err == 0.0) then ; return ! The exact solution on the first guess...
  elseif (err > 0.0) then ! The root is properly bracketed.
    max_ea = ea_kb(i) ; err_max = err
    !   Use Newton's method (if it stays bounded) or the false position method
    ! to find the next value.
    if ((derr_dea > 0.0) .and. (derr_dea*(ea_kb(i) - min_ea) > err) .and. &
        (derr_dea*(max_ea - ea_kb(i)) > -1.0*err)) then
      ea_kb(i) = ea_kb(i) - err / derr_dea
    else ! Use the bisection for the next guess.
      ea_kb(i) = 0.5*(max_ea+min_ea)
    endif
  else
    !   Try to bracket the root first.  If unable to bracket the root, return
    ! the maximum.
    zeros(i) = 0.0
    call find_maxF_kb(h_bl, Sref, Ent_bl, I_dSkbp1, zeros, ea_kb, &
                      kmb, i, i, G, GV, CS, maxF, ent_maxF, F_thresh = F_kb)
    err_max = dS_kbp1 * maxF(i) - val
    ! If err_max is negative, there is no good solution, so use the maximum
    ! value of F in the valid range.
    if (err_max <= 0.0) then
      ea_kb(i) = ent_maxF(i) ; return
    else
      max_ea = ent_maxF(i)
      ea_kb(i) = 0.5*(max_ea+min_ea) ! Use bisection for the next guess.
    endif
  endif

  ! Exit if the range between max_ea and min_ea already acceptable.
  ! if (abs(max_ea - min_ea) < 0.1*tolerance) return

  do it = 1, MAXIT
    call determine_dSkb(h_bl, Sref, Ent_bl, ea_kb, i, i, kmb, G, GV, .true., &
                        dS_kb, ddSkb_dE)

    err = dS_kb(i) * ea_kb(i) - val
    derr_dea = dS_kb(i) + ddSkb_dE(i) * ea_kb(i)

    ea_prev = ea_kb(i)
    ! Use Newton's method or the false position method to find the next value.
    Newton = .false.
    if (err > 0.0) then
      max_ea = ea_kb(i) ; err_max = err
      if ((derr_dea > 0.0) .and. (derr_dea*(ea_kb(i)-min_ea) > err)) Newton = .true.
    else
      min_ea = ea_kb(i) ; err_min = err
      if ((derr_dea > 0.0) .and. (derr_dea*(ea_kb(i)-max_ea) < err)) Newton = .true.
    endif

    if (Newton) then
      ea_kb(i) = ea_kb(i) - err / derr_dea
    elseif (bisect_next) then ! Use bisection to reduce the range.
      ea_kb(i) = 0.5*(max_ea+min_ea)
      bisect_next = .false.
    else  ! Use the false-position method for the next guess.
      ea_kb(i) = min_ea + (max_ea-min_ea) * (err_min/(err_min - err_max))
      bisect_next = .true.
    endif

    tol1 = tolerance ; if (err > 0.0) tol1 = 0.099*tolerance
    if (dS_kb(i) <= dS_kbp1) then
      if (abs(ea_kb(i) - ea_prev) <= tol1) return
    else
      if (dS_kbp1*abs(ea_kb(i) - ea_prev) <= dS_kb(i)*tol1) return
    endif
  enddo

end subroutine F_kb_to_ea_kb


subroutine determine_Ea_kb(h_bl, dtKd_kb, Sref, I_dSkbp1, Ent_bl, ea_kbp1, &
                           min_eakb, max_eakb, kmb, is, ie, do_i, G, GV, CS, Ent, &
                           error, err_min_eakb0, err_max_eakb0, F_kb, dFdfm_kb)
  type(ocean_grid_type),            intent(in)  :: G
  type(verticalGrid_type),          intent(in)  :: GV
  real, dimension(SZI_(G),SZK_(G)), intent(in)  :: h_bl, Sref, Ent_bl
  real, dimension(SZI_(G)),         intent(in)  :: I_dSkbp1, dtKd_kb, ea_kbp1
  real, dimension(SZI_(G)),         intent(in)  :: min_eakb, max_eakb
  integer,                          intent(in)  :: kmb, is, ie
  logical, dimension(SZI_(G)),      intent(in)  :: do_i
  type(entrain_diffusive_CS),       pointer     :: CS
  real, dimension(SZI_(G)),         intent(inout) :: Ent
  real, dimension(SZI_(G)),     intent(out), optional :: error
  real, dimension(SZI_(G)),     intent(in),  optional :: err_min_eakb0, err_max_eakb0
  real, dimension(SZI_(G)),     intent(out), optional :: F_kb, dFdfm_kb
! Arguments: h_bl - Layer thickness, with the top interior layer at k-index
!                   kmb+1, in units of m or kg m-2 (abbreviated as H below).
!  (in)      dtKd_kb - The diapycnal diffusivity in the top interior layer times
!                      the time step, in H2.
!  (in)      Sref - The coordinate reference potential density, with the
!                   value of the topmost interior layer at layer kmb+1,
!                   in units of kg m-3.
!  (in)      I_dSkbp1 - The inverse of the difference in reference potential
!                       density across the base of the uppermost interior layer,
!                       in units of m3 kg-1.
!  (in)      Ent_bl - The average entrainment upward and downward across
!                     each interface around the buffer layers, in H.
!  (in)      ea_kbp1 - The entrainment from above by layer kb+1, in H.
!  (in)      min_eakb - The minimum permissible rate of entrainment, in H.
!  (in)      max_eakb - The maximum permissible rate of entrainment, in H.
!  (in)      is, ie - The range of i-indices to work on.
!  (in)      do_i - A logical variable indicating which i-points to work on.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - This module's control structure.
!  (in/out)  Ent - The entrainment rate of the uppermost interior layer, in H.
!                  The input value is the first guess.
!  (out,opt) error - The error (locally defined in this routine) associated with
!                    the returned solution.
!  (in,opt)  error_min_eakb0, error_max_eakb0 - The errors (locally defined)
!                    associated with min_eakb and max_eakb when ea_kbp1 = 0,
!                    returned from a previous call to this routine.
!  (out,opt) F_kb - The entrainment from below by the uppermost interior layer
!                   corresponding to the returned value of Ent, in H.
!  (out,out) dFdfm_kb - The partial derivative of F_kb with ea_kbp1, nondim.

!  This subroutine determines the entrainment from above by the top interior
! layer (labeled kb elsewhere) given an entrainment by the layer below it,
! constrained to be within the provided bounds.
  real, dimension(SZI_(G)) :: &
    dS_kb, &                !   The coordinate-density difference between the
                            ! layer kb and deepest buffer layer, limited to
                            ! ensure that it is positive, in kg m-3.
    dS_Lay, &               !   The coordinate-density difference across layer
                            ! kb, limited to ensure that it is positive and not
                            ! too much bigger than dS_kb or dS_kbp1, in kg m-3.
    ddSkb_dE, ddSlay_dE, &  ! The derivatives of dS_kb and dS_Lay with E,
                            ! in units of kg m-3 H-1.
    derror_dE, &            ! The derivative of err with E, in H.
    err, &                  ! The "error" whose zero is being sought, in H2.
    E_min, E_max, &         ! The minimum and maximum values of E, in H.
    error_minE, error_maxE  ! err when E = E_min or E = E_max, in H2.
  real :: err_est           ! An estimate of what err will be, in H2.
  real :: eL                ! 1 or 0, depending on whether increases in E lead
                            ! to decreases in the entrainment from below by the
                            ! deepest buffer layer.
  real :: fa, fk, fm, fr    ! Temporary variables used to calculate err, in ND, H2, H, H.
  real :: tolerance         ! The tolerance within which E must be converged, in H.
  real :: E_prev            ! The previous value of E, in H.
  logical, dimension(SZI_(G)) :: false_position ! If true, the false position
                            ! method might be used for the next iteration.
  logical, dimension(SZI_(G)) :: redo_i ! If true, more work is needed on this column.
  logical :: do_any
  real, parameter :: LARGE_VAL = 1.0e30
  integer :: i, it
  integer, parameter :: MAXIT = 30

  if (.not.CS%bulkmixedlayer) then
    call MOM_error(FATAL, "determine_Ea_kb should not be called "//&
                           "unless BULKMIXEDLAYER is defined.")
  endif
  tolerance = GV%m_to_H * CS%Tolerance_Ent

  do i=is,ie ; redo_i(i) = do_i(i) ; enddo

  do i=is,ie ; if (do_i(i)) then
    ! The first guess of Ent was the value from the previous iteration.

    !   These were previously calculated and provide good limits and estimates
    ! of the errors there. By construction the errors increase with R*ea_kbp1.
    E_min(i) = min_eakb(i) ; E_max(i) = max_eakb(i)
    error_minE(i) = -LARGE_VAL ; error_maxE(i) = LARGE_VAL
    false_position(i) = .true. ! Used to alternate between false_position and
                               ! bisection when Newton's method isn't working.
    if (present(err_min_eakb0)) error_minE(i) = err_min_eakb0(i) - E_min(i) * ea_kbp1(i)
    if (present(err_max_eakb0)) error_maxE(i) = err_max_eakb0(i) - E_max(i) * ea_kbp1(i)

    if ((error_maxE(i) <= 0.0) .or. (error_minE(i) >= 0.0)) then
      ! The root is not bracketed and one of the limiting values should be used.
      if (error_maxE(i) <= 0.0) then
        ! The errors decrease with E*ea_kbp1, so E_max is the best solution.
        Ent(i) = E_max(i) ; err(i) = error_maxE(i)
      else  ! error_minE >= 0 is equivalent to ea_kbp1 = 0.0.
        Ent(i) = E_min(i) ; err(i) = error_minE(i)
      endif
      derror_dE(i) = 0.0
      redo_i(i) = .false.
    endif
  endif ; enddo   ! End of i-loop

  do it = 1,MAXIT
    do_any = .false. ; do i=is,ie ; if (redo_i(i)) do_any = .true. ; enddo
    if (.not.do_any) exit
    call determine_dSkb(h_bl, Sref, Ent_bl, Ent, is, ie, kmb, G, GV, .true., dS_kb, &
                        ddSkb_dE, dS_lay, ddSlay_dE, do_i_in = redo_i)
    do i=is,ie ; if (redo_i(i)) then
      !  The correct root is bracketed between E_min and E_max.
      ! Note the following limits:  Ent >= 0 ; fa > 1 ; fk > 0
      eL = 0.0 ; if (2.0*Ent_bl(i,Kmb+1) >= Ent(i)) eL = 1.0
      fa = (1.0 + eL) + dS_kb(i)*I_dSkbp1(i)
      fk = dtKd_kb(i) * (dS_Lay(i)/dS_kb(i))
      fm = (ea_kbp1(i) - h_bl(i,kmb+1)) + eL*2.0*Ent_bl(i,Kmb+1)
      if (fm > -GV%Angstrom) fm = fm + GV%Angstrom  ! This could be smooth if need be.
      err(i) = (fa * Ent(i)**2 - fm * Ent(i)) - fk
      derror_dE(i) = ((2.0*fa + (ddSkb_dE(i)*I_dSkbp1(i))*Ent(i))*Ent(i) - fm) - &
          dtKd_kb(i) * (ddSlay_dE(i)*dS_kb(i) - ddSkb_dE(i)*dS_Lay(i))/(dS_kb(i)**2)

      if (err(i) == 0.0) then
        redo_i(i) = .false. ; cycle
      elseif (err(i) > 0.0) then
        E_max(i) = Ent(i) ; error_maxE(i) = err(i)
      else
        E_min(i) = Ent(i) ; error_minE(i) = err(i)
      endif

      E_prev = Ent(i)
      if ((it == 1) .or. (derror_dE(i) <= 0.0)) then
        !   Assuming that the coefficients of the quadratic equation are correct
        ! will usually give a very good first guess.  Also, if derror_dE < 0.0,
        ! R is on the wrong side of the approximate parabola.  In either case,
        ! try assuming that the error is approximately a parabola and solve.
        fr = sqrt(fm**2 + 4.0*fa*fk)
        if (fm >= 0.0) then
          Ent(i) = (fm + fr) / (2.0 * fa)
        else
          Ent(i) = (2.0 * fk) / (fr - fm)
        endif
        ! But make sure that the root stays bracketed, bisecting if needed.
        if ((Ent(i) > E_max(i)) .or. (Ent(i) < E_min(i))) &
          Ent(i) = 0.5*(E_max(i) + E_min(i))
      elseif (((E_max(i)-Ent(i))*derror_dE(i) > -err(i)) .and. &
              ((Ent(i)-E_min(i))*derror_dE(i) > err(i)) ) then
        ! Use Newton's method for the next estimate, provided it will
        ! remain bracketed between Rmin and Rmax.
        Ent(i) = Ent(i) - err(i) / derror_dE(i)
      elseif (false_position(i) .and. &
              (error_maxE(i) - error_minE(i) < 0.9*LARGE_VAL)) then
        ! Use the false postion method if there are decent error estimates.
        Ent(i) = E_min(i) + (E_max(i)-E_min(i)) * &
                (-error_minE(i)/(error_maxE(i) - error_minE(i)))
        false_position(i) = .false.
      else ! Bisect as a last resort or if the false position method was used last.
        Ent(i) = 0.5*(E_max(i) + E_min(i))
        false_position(i) = .true.
      endif

      if (abs(E_prev - Ent(i)) < tolerance) then
        err_est = err(i) + (Ent(i) - E_prev) * derror_dE(i)
        if ((it > 1) .or. (err_est*err(i) <= 0.0) .or. &
            (abs(err_est) < abs(tolerance*derror_dE(i)))) redo_i(i) = .false.
      endif

    endif ; enddo   ! End of i-loop
  enddo ! End of iterations to determine Ent(i).

  ! Update the value of dS_kb for consistency with Ent.
  if (present(F_kb) .or. present(dFdfm_kb)) &
    call determine_dSkb(h_bl, Sref, Ent_bl, Ent, is, ie, kmb, G, GV, .true., &
                        dS_kb, do_i_in = do_i)

  if (present(F_kb)) then ; do i=is,ie ; if (do_i(i)) then
    F_kb(i) = Ent(i) * (dS_kb(i) * I_dSkbp1(i))
  endif ; enddo ; endif
  if (present(error)) then ; do i=is,ie ; if (do_i(i)) then
    error(i) = err(i)
  endif ; enddo ; endif
  if (present(dFdfm_kb)) then ; do i=is,ie ; if (do_i(i)) then
    !   derror_dE and ddSkb_dE are _not_ recalculated here, since dFdfm_kb is
    ! only used in Newton's method, and slightly increasing the accuracy of the
    ! estimate is unlikely to speed convergence.
    if (derror_dE(i) > 0.0) then
      dFdfm_kb(i) = ((dS_kb(i) + Ent(i) * ddSkb_dE(i)) * I_dSkbp1(i)) * &
                    (Ent(i) / derror_dE(i))
    else ! Use Adcroft's division by 0 convention.
      dFdfm_kb(i) = 0.0
    endif
  endif ; enddo ; endif

end subroutine determine_Ea_kb

subroutine find_maxF_kb(h_bl, Sref, Ent_bl, I_dSkbp1, min_ent_in, max_ent_in, &
                        kmb, is, ie, G, GV, CS, maxF, ent_maxF, do_i_in, &
                        F_lim_maxent, F_thresh)
  type(ocean_grid_type),        intent(in)  :: G
  type(verticalGrid_type),      intent(in)  :: GV
  real, dimension(SZI_(G),SZK_(G)), intent(in) :: h_bl, Sref, Ent_bl
  real, dimension(SZI_(G)),     intent(in)  :: I_dSkbp1, min_ent_in, max_ent_in
  integer,                      intent(in)  :: kmb, is, ie
  type(entrain_diffusive_CS),   pointer     :: CS
  real, dimension(SZI_(G)),     intent(out) :: maxF
  real, dimension(SZI_(G)),     intent(out), optional :: ent_maxF
  logical, dimension(SZI_(G)),  intent(in),  optional :: do_i_in
  real, dimension(SZI_(G)),     intent(out), optional :: F_lim_maxent
  real, dimension(SZI_(G)),     intent(in),  optional :: F_thresh
! Arguments: h_bl - Layer thickness, in m or kg m-2 (abbreviated as H below).
!  (in)      Sref - Reference potential density (in kg m-3?)
!  (in)      Ent_bl - The average entrainment upward and downward across
!                     each interface around the buffer layers, in H.
!  (in)      I_dSkbp1 - The inverse of the difference in reference potential
!                       density across the base of the uppermost interior layer,
!                       in units of m3 kg-1.
!  (in)      min_ent_in - The minimum value of ent to search, in H.
!  (in)      max_ent_in - The maximum value of ent to search, in H.
!  (in)      is, ie - The range of i-indices to work on.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      CS - This module's control structure.
!  (out)     maxF - The maximum value of F = ent*ds_kb*I_dSkbp1 found in the
!                   range min_ent < ent < max_ent, in H.
!  (out,opt) ent_maxF - The value of ent at that maximum, in H.
!  (in, opt) do_i_in - A logical array indicating which columns to work on.
!  (out,opt) F_lim_maxent - If present, do not apply the limit in finding the
!                           maximum value, but return the limited value at
!                           ent=max_ent_in in this array, in H.
!  (in, opt) F_thresh - If F_thresh is present, return the first value found
!                       that has F > F_thresh, or the maximum.

! Maximize F = ent*ds_kb*I_dSkbp1 in the range min_ent < ent < max_ent.
! ds_kb may itself be limited to positive values in determine_dSkb, which gives
! the prospect of two local maxima in the range - one at max_ent_in with that
! minimum value of ds_kb, and the other due to the unlimited (potentially
! negative) value.  It is faster to find the true maximum by first finding the
! unlimited maximum and comparing it to the limited value at max_ent_in.
  real, dimension(SZI_(G)) :: &
    ent, &
    minent, maxent, ent_best, &
    F_max_ent_in, &
    F_maxent, F_minent, F, F_best, &
    dF_dent, dF_dE_max, dF_dE_min, dF_dE_best, &
    dS_kb, dS_kb_lim, ddSkb_dE, dS_anom_lim, &
    chg_prev, chg_pre_prev
  real :: dF_dE_mean, maxslope, minslope
  real :: tolerance
  real :: ratio_select_end
  real :: rat, max_chg, min_chg, chg1, chg2, chg
  logical, dimension(SZI_(G)) :: do_i, last_it, need_bracket, may_use_best
  logical :: doany, OK1, OK2, bisect, new_min_bound
  integer :: i, it, is1, ie1
  integer, parameter :: MAXIT = 20

  tolerance = GV%m_to_H * CS%Tolerance_Ent

  if (present(do_i_in)) then
    do i=is,ie ; do_i(i) = do_i_in(i) ; enddo
  else
    do i=is,ie ; do_i(i) = .true. ; enddo
  endif

  ! The most likely value is at max_ent.
  call determine_dSkb(h_bl, Sref, Ent_bl, max_ent_in, is, ie, kmb, G, GV, .false., &
                      dS_kb, ddSkb_dE , dS_anom_lim=dS_anom_lim)
  ie1 = is-1 ; doany = .false.
  do i=is,ie
    dS_kb_lim(i) = dS_kb(i) + dS_anom_lim(i)
    F_max_ent_in(i) = max_ent_in(i)*dS_kb_lim(i)*I_dSkbp1(i)
    maxent(i) = max_ent_in(i) ; minent(i) = min_ent_in(i)
    if ((abs(maxent(i) - minent(i)) < tolerance) .or. (.not.do_i(i))) then
      F_best(i) = max_ent_in(i)*dS_kb(i)*I_dSkbp1(i)
      ent_best(i) = max_ent_in(i) ; ent(i) = max_ent_in(i)
      do_i(i) = .false.
    else
      F_maxent(i) = maxent(i) * dS_kb(i) * I_dSkbp1(i)
      dF_dE_max(i) = (dS_kb(i) + maxent(i)*ddSkb_dE(i)) * I_dSkbp1(i)
      doany = .true. ; last_it(i) = .false. ; need_bracket(i) = .true.
    endif
  enddo

  if (doany) then
    ie1 = is-1 ; do i=is,ie ; if (do_i(i)) ie1 = i ; enddo
    do i=ie1,is,-1 ; if (do_i(i)) is1 = i ; enddo
    ! Find the value of F and its derivative at min_ent.
    call determine_dSkb(h_bl, Sref, Ent_bl, minent, is1, ie1, kmb, G, GV, .false., &
                        dS_kb, ddSkb_dE, do_i_in = do_i)
    do i=is1,ie1 ; if (do_i(i)) then
      F_minent(i) = minent(i) * dS_kb(i) * I_dSkbp1(i)
      dF_dE_min(i) = (dS_kb(i) + minent(i)*ddSkb_dE(i)) * I_dSkbp1(i)
    endif ; enddo

    ratio_select_end = 0.9
    do it=1,MAXIT
      ratio_select_end = 0.5*ratio_select_end
      do i=is1,ie1 ; if (do_i(i)) then
        if (need_bracket(i)) then
          dF_dE_mean = (F_maxent(i) - F_minent(i)) / (maxent(i) - minent(i))
          maxslope = MAX(dF_dE_mean, dF_dE_min(i), dF_dE_max(i))
          minslope = MIN(dF_dE_mean, dF_dE_min(i), dF_dE_max(i))
          if (F_minent(i) >= F_maxent(i)) then
            if (dF_dE_min(i) > 0.0) then ; rat = 0.02 ! A small step should bracket the soln.
            elseif (maxslope < ratio_select_end*minslope) then
              ! The maximum of F is at minent.
              F_best(i) = F_minent(i) ; ent_best(i) = minent(i) ; rat = 0.0
              do_i(i) = .false.
            else ; rat = 0.382 ; endif ! Use the golden ratio
          else
            if (dF_dE_max(i) < 0.0) then ; rat = 0.98 ! A small step should bracket the soln.
            elseif (minslope > ratio_select_end*maxslope) then
              ! The maximum of F is at maxent.
              F_best(i) = F_maxent(i) ; ent_best(i) = maxent(i) ; rat = 1.0
              do_i(i) = .false.
            else ; rat = 0.618 ; endif ! Use the golden ratio
          endif

          if (rat >= 0.0) ent(i) = rat*maxent(i) + (1.0-rat)*minent(i)
          if (((maxent(i) - minent(i)) < tolerance) .or. (it==MAXIT)) &
            last_it(i) = .true.
        else ! The maximum is bracketed by minent, ent_best, and maxent.
          chg1 = 2.0*(maxent(i) - minent(i)) ; chg2 = chg1
          if (dF_dE_best(i) > 0) then
            max_chg = maxent(i) - ent_best(i) ; min_chg = 0.0
          else
            max_chg = 0.0 ; min_chg = minent(i) - ent_best(i) ! < 0
          endif
          if (max_chg - min_chg < 2.0*tolerance) last_it(i) = .true.
          if (dF_dE_max(i) /= dF_dE_best(i)) &
            chg1 = (maxent(i) - ent_best(i))*dF_dE_best(i) / &
                   (dF_dE_best(i) - dF_dE_max(i))
          if (dF_dE_min(i) /= dF_dE_best(i)) &
            chg2 = (minent(i) - ent_best(i))*dF_dE_best(i) / &
                   (dF_dE_best(i) - dF_dE_min(i))
          OK1 = ((chg1 < max_chg) .and. (chg1 > min_chg))
          OK2 = ((chg2 < max_chg) .and. (chg2 > min_chg))
          if (.not.(OK1 .or. OK2)) then ; bisect = .true. ; else
            if (OK1 .and. OK2) then ! Take the acceptable smaller change.
              chg = chg1 ; if (abs(chg2) < abs(chg1)) chg = chg2
            elseif (OK1) then ; chg = chg1
            else ; chg = chg2 ; endif
            if (abs(chg) > 0.5*abs(chg_pre_prev(i))) then ; bisect = .true.
            else ; bisect = .false. ; endif
          endif
          chg_pre_prev(i) = chg_prev(i)
          if (bisect) then
            if (dF_dE_best(i) > 0.0) then
              ent(i) = 0.5*(maxent(i) + ent_best(i))
              chg_prev(i) = 0.5*(maxent(i) - ent_best(i))
            else
              ent(i) = 0.5*(minent(i) + ent_best(i))
              chg_prev(i) = 0.5*(minent(i) - ent_best(i))
            endif
          else
            if (abs(chg) < tolerance) chg = SIGN(tolerance,chg)
            ent(i) = ent_best(i) + chg
            chg_prev(i) = chg
          endif
        endif
      endif ; enddo

      if (mod(it,3) == 0) then  ! Re-determine the loop bounds.
        ie1 = is-1 ; do i=is1,ie ; if (do_i(i)) ie1 = i ; enddo
        do i=ie1,is,-1 ; if (do_i(i)) is1 = i ; enddo
      endif

      call determine_dSkb(h_bl, Sref, Ent_bl, ent, is1, ie1, kmb, G, GV, .false., &
                          dS_kb, ddSkb_dE, do_i_in = do_i)
      do i=is1,ie1 ; if (do_i(i)) then
        F(i) = ent(i)*dS_kb(i)*I_dSkbp1(i)
        dF_dent(i) = (dS_kb(i) + ent(i)*ddSkb_dE(i)) * I_dSkbp1(i)
      endif ; enddo

      if (present(F_thresh)) then ; do i=is1,ie1 ; if (do_i(i)) then
        if (F(i) >= F_thresh(i)) then
          F_best(i) = F(i) ; ent_best(i) = ent(i) ; do_i(i) = .false.
        endif
      endif ; enddo ; endif

      doany = .false.
      do i=is1,ie1 ; if (do_i(i)) then
        if (.not.last_it(i)) doany = .true.
        if (last_it(i)) then
          if (need_bracket(i)) then
            if ((F(i) > F_maxent(i)) .and. (F(i) > F_minent(i))) then
              F_best(i) = F(i) ; ent_best(i) = ent(i)
            elseif (F_maxent(i) > F_minent(i)) then
              F_best(i) = F_maxent(i) ; ent_best(i) = maxent(i)
            else
              F_best(i) = F_minent(i) ; ent_best(i) = minent(i)
            endif
          elseif (F(i) > F_best(i)) then
            F_best(i) = F(i) ; ent_best(i) = ent(i)
          endif
          do_i(i) = .false.
        elseif (need_bracket(i)) then
          if ((F(i) > F_maxent(i)) .and. (F(i) > F_minent(i))) then
            need_bracket(i) = .false. ! The maximum is now bracketed.
            chg_prev(i) = (maxent(i) - minent(i))
            chg_pre_prev(i) = 2.0*chg_prev(i)
            ent_best(i) = ent(i) ; F_best(i) = F(i) ; dF_dE_best(i) = dF_dent(i)
          elseif ((F(i) <= F_maxent(i)) .and. (F(i) > F_minent(i))) then
            new_min_bound = .true.  ! We have a new minimum bound.
          elseif ((F(i) <= F_maxent(i)) .and. (F(i) > F_minent(i))) then
            new_min_bound = .false. ! We have a new maximum bound.
          else ! This case would bracket a minimum.  Wierd.
             ! Unless the derivative indicates that there is a maximum near the
             ! lower bound, try keeping the end with the larger value of F;
             ! in a tie keep the minimum as the answer here will be compared
             ! with the maximum input value later.
             new_min_bound = .true.
             if (dF_dE_min(i) > 0.0 .or. (F_minent(i) >= F_maxent(i))) &
               new_min_bound = .false.
          endif
          if (need_bracket(i)) then ! Still not bracketed.
            if (new_min_bound) then
              minent(i) = ent(i) ; F_minent(i) = F(i) ; dF_dE_min(i) = dF_dent(i)
            else
              maxent(i) = ent(i) ; F_maxent(i) = F(i) ; dF_dE_max(i) = dF_dent(i)
            endif
          endif
        else  ! The root was previously bracketed.
          if (F(i) >= F_best(i)) then ! There is a new maximum.
            if (ent(i) > ent_best(i)) then ! Replace minent with ent_prev.
              minent(i) = ent_best(i) ; F_minent(i) = F_best(i) ; dF_dE_min(i) = dF_dE_best(i)
            else ! Replace maxent with ent_best.
              maxent(i) = ent_best(i) ; F_maxent(i) = F_best(i) ; dF_dE_max(i) = dF_dE_best(i)
            endif
            ent_best(i) = ent(i) ; F_best(i) = F(i) ; dF_dE_best(i) = dF_dent(i)
          else
            if (ent(i) < ent_best(i)) then ! Replace the minent with ent.
              minent(i) = ent(i) ; F_minent(i) = F(i) ; dF_dE_min(i) = dF_dent(i)
            else ! Replace maxent with ent_prev.
              maxent(i) = ent(i) ; F_maxent(i) = F(i) ; dF_dE_max(i) = dF_dent(i)
            endif
          endif
          if ((maxent(i) - minent(i)) <= tolerance) do_i(i) = .false. ! Done.
        endif ! need_bracket.
      endif ; enddo
      if (.not.doany) exit
    enddo
  endif

  if (present(F_lim_maxent)) then
    ! Return the unlimited maximum in maxF, and the limited value of F at maxent.
    do i=is,ie
      maxF(i) = F_best(i)
      F_lim_maxent(i) = F_max_ent_in(i)
      if (present(ent_maxF)) ent_maxF(i) = ent_best(i)
    enddo
  else
    ! Now compare the two? potential maxima using the limited value of dF_kb.
    doany = .false.
    do i=is,ie
      may_use_best(i) = (ent_best(i) /= max_ent_in(i))
      if (may_use_best(i)) doany = .true.
    enddo
    if (doany) then
      ! For efficiency, could save previous value of dS_anom_lim_best?
      call determine_dSkb(h_bl, Sref, Ent_bl, ent_best, is, ie, kmb, G, GV, .true., &
                          dS_kb_lim)
      do i=is,ie
        F_best(i) = ent_best(i)*dS_kb_lim(i)*I_dSkbp1(i)
        ! The second test seems necessary because of roundoff differences that
        ! can arise during compilation.
        if ((F_best(i) > F_max_ent_in(i)) .and. (may_use_best(i))) then
          maxF(i) = F_best(i)
          if (present(ent_maxF)) ent_maxF(i) = ent_best(i)
        else
          maxF(i) = F_max_ent_in(i)
          if (present(ent_maxF)) ent_maxF(i) = max_ent_in(i)
        endif
      enddo
    else
      ! All of the maxima are at the maximum entrainment.
      do i=is,ie ; maxF(i) = F_max_ent_in(i) ; enddo
      if (present(ent_maxF)) then
        do i=is,ie ; ent_maxF(i) = max_ent_in(i) ; enddo
      endif
    endif
  endif

end subroutine find_maxF_kb

subroutine entrain_diffusive_init(Time, G, GV, param_file, diag, CS)
  type(time_type),         intent(in)    :: Time
  type(ocean_grid_type),   intent(in)    :: G
  type(verticalGrid_type),               intent(in)    :: GV
  type(param_file_type),   intent(in)    :: param_file
  type(diag_ctrl), target, intent(inout) :: diag
  type(entrain_diffusive_CS), pointer     :: CS
! Arguments: Time - The current model time.
!  (in)      G - The ocean's grid structure.
!  (in)      GV - The ocean's vertical grid structure.
!  (in)      param_file - A structure indicating the open file to parse for
!                         model parameter values.
!  (in)      diag - A structure that is used to regulate diagnostic output.
!  (in/out)  CS - A pointer that is set to point to the control structure
!                 for this module
  real :: decay_length, dt, Kd
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mod  = "MOM_entrain_diffusive" ! This module's name.

  if (associated(CS)) then
    call MOM_error(WARNING, "entrain_diffusive_init called with an associated "// &
                            "control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  CS%bulkmixedlayer = (GV%nkml > 0)

! Set default, read and log parameters
  call log_version(param_file, mod, version, "")
  call get_param(param_file, mod, "CORRECT_DENSITY", CS%correct_density, &
                 "If true, and USE_EOS is true, the layer densities are \n"//&
                 "restored toward their target values by the diapycnal \n"//&
                 "mixing, as described in Hallberg (MWR, 2000).", &
                 default=.true.)
  call get_param(param_file, mod, "MAX_ENT_IT", CS%max_ent_it, &
                 "The maximum number of iterations that may be used to \n"//&
                 "calculate the interior diapycnal entrainment.", default=5)
! In this module, KD is only used to set the default for TOLERANCE_ENT. (m2 s-1)
  call get_param(param_file, mod, "KD", Kd, fail_if_missing=.true.)
  call get_param(param_file, mod, "DT", dt, &
                 "The (baroclinic) dynamics time step.", units = "s", &
                 fail_if_missing=.true.)
! CS%Tolerance_Ent = MAX(100.0*GV%Angstrom,1.0e-4*sqrt(dt*Kd)) !
  call get_param(param_file, mod, "TOLERANCE_ENT", CS%Tolerance_Ent, &
                 "The tolerance with which to solve for entrainment values.", &
                 units="m", default=MAX(100.0*GV%Angstrom,1.0e-4*sqrt(dt*Kd)))

  CS%id_Kd = register_diag_field('ocean_model', 'Kd_effective', diag%axesTL, Time, &
      'Diapycnal diffusivity as applied', 'meter2 second-1')
  CS%id_diff_work = register_diag_field('ocean_model', 'diff_work', diag%axesTi, Time, &
      'Work actually done by diapycnal diffusion across each interface', 'W m-2')

end subroutine entrain_diffusive_init

subroutine entrain_diffusive_end(CS)
  type(entrain_diffusive_CS), pointer :: CS

  if (associated(CS)) deallocate(CS)

end subroutine entrain_diffusive_end

end module MOM_entrain_diffusive
