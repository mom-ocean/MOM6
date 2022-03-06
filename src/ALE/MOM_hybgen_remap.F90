!> This module contains the hybgen remapping routines from HYCOM, with minor
!! modifications to follow the MOM6 coding conventions
module MOM_hybgen_remap

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

public hybgen_plm_coefs, hybgen_ppm_coefs, hybgen_weno_coefs

contains

!> Set up the coefficients for PLM remapping of a set of scalars
subroutine hybgen_plm_coefs(si, dpi, slope, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: si(nk,ns) !< The cell-averaged input scalar fields [A]
  real,    intent(in)  :: dpi(nk)   !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: slope(nk,ns) !< The PLM slope times cell width [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: piecewise linear across each input cell with
!             monotonized central-difference limiter.
!
!     van Leer, B., 1977, J. Comp. Phys., 23 276-299.
!
!  2) input arguments:
!       si    - initial scalar fields in pi-layer space
!       dpi   - initial layer thicknesses (dpi(k) = pi(k+1)-pi(k))
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       slope - coefficients for hybgen_plm_remap
!                profile(y) = si+slope*(y-1),  -0.5 <= y <= 0.5
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  real :: qcen   ! A layer's thickness divided by the distance between the centers
                 ! of the adjacent cells, usually ~0.5, but always <= 1 [nondim]
  real :: zbot, zcen, ztop ! Tracer slopes times the layer thickness [A]
  integer :: i, k

  do i=1,ns
    slope(1, i) = 0.0
    slope(nk,i) = 0.0
  enddo !i
  do k= 2,nk-1
    if (dpi(k) <= thin) then  !use PCM
      do i=1,ns ; slope(k,i) = 0.0 ; enddo
    else
! ---     use qcen in place of 0.5 to allow for non-uniform grid
      qcen = dpi(k) / (dpi(k)+0.5*(dpi(k-1)+dpi(k+1)))  !dpi(k)>thin
      do i=1,ns
! ---       PLM (non-zero slope, but no new extrema)
! ---       layer value is si-0.5*slope at top    interface,
! ---                  and si+0.5*slope at bottom interface.
!
! ---       monotonized central-difference limiter (van Leer, 1977,
! ---       JCP 23 pp 276-299).  For a discussion of PLM limiters, see
! ---       Finite Volume Methods for Hyperbolic Problems by R.J. Leveque.
        ztop = 2.0*(si(k,  i)-si(k-1,i))
        zbot = 2.0*(si(k+1,i)-si(k,  i))
        zcen = qcen*(si(k+1,i)-si(k-1,i))
        if     (ztop*zbot > 0.0) then !ztop,zbot are the same sign
          slope(k,i) = sign(min(abs(zcen),abs(zbot),abs(ztop)), zbot)
        else
          slope(k,i) = 0.0  !local extrema, so no slope
        endif
      enddo !i
    endif  !PCM:PLM
  enddo !k

  if (present(PCM_lay)) then
    do k=1,nk ; if (PCM_lay(k)) then
      do i=1,ns ; slope(k,i) = 0.0 ; enddo
    endif ; enddo
  endif

end subroutine hybgen_plm_coefs


!> Set up the coefficients for PPM remapping of a set of scalars
subroutine hybgen_ppm_coefs(s, h_src, edges, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The scalar fields to work on
  real,    intent(in)  :: s(nk,ns)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2,ns) !< The PPM interpolation edge values of the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic piecewise parabolic across each input cell
!
!     Colella, P. & P.R. Woodward, 1984, J. Comp. Phys., 54, 174-201.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       edges - cell edge scalar values for the PPM reconstruction
!                edges.1 is value at interface above
!                edges.2 is value at interface below
!
!  4) Tim Campbell, Mississippi State University, October 2002.
!     Alan J. Wallcraft,  Naval Research Laboratory,  Aug. 2007.
!-----------------------------------------------------------------------
!
  real :: dp(nk) ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  logical :: PCM_layer(nk) ! True for layers that should use PCM remapping, either because they are
                           ! very thin, or because this is specified by PCM_lay.
  real :: da        ! Difference between the unlimited scalar edge value estimates [A]
  real :: a6        ! Scalar field differences that are proportional to the curvature [A]
  real :: slk, srk  ! Differences between adjacent cell averages of scalars [A]
  real :: sck       ! Scalar differences across a cell.
  real :: as(nk)    ! Scalar field difference across each cell [A]
  real :: al(nk), ar(nk)   ! Scalar field at the left and right edges of a cell [A]
  real :: h112(nk+1), h122(nk+1)  ! Combinations of thicknesses [H ~> m or kg m-2]
  real :: I_h12(nk+1) ! Inverses of combinations of thickesses [H-1 ~> m-1 or m2 kg-1]
  real :: h2_h123(nk)  ! A ratio of a layer thickness of the sum of 3 adjacent thicknesses [nondim]
  real :: I_h0123(nk)     ! Inverse of the sum of 4 adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: h01_h112(nk+1) ! A ratio of sums of adjacent thicknesses [nondim], 2/3 in the limit of uniform thicknesses.
  real :: h23_h122(nk+1) ! A ratio of sums of adjacent thicknesses [nondim], 2/3 in the limit of uniform thicknesses.
  integer :: k, i

  ! This PPM remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
  do k=1,nk ; dp(k) = max(h_src(k), thin) ; enddo

  ! Specify the layers that will use PCM remapping.
  if (present(PCM_lay)) then
    do k=1,nk ; PCM_layer(k) = (PCM_lay(k) .or. dp(k) <= thin) ; enddo
  else
    do k=1,nk ; PCM_layer(k) = (dp(k) <= thin) ; enddo
  endif

  !compute grid metrics
  do k=2,nk
    h112(K) = 2.*dp(k-1) + dp(k)
    h122(K) = dp(k-1) + 2.*dp(k)
    I_h12(K) = 1.0 / (dp(k-1) + dp(k))
  enddo !k
  do k=2,nk-1
    h2_h123(k) = dp(k) / (dp(k) + (dp(k-1)+dp(k+1)))
  enddo
  do K=3,nk-1
    I_h0123(K) = 1.0 / ((dp(k-2) + dp(k-1)) + (dp(k) + dp(k+1)))

    h01_h112(K) = (dp(k-2) + dp(k-1)) / (2.0*dp(k-1) + dp(k))
    h23_h122(K) = (dp(k) + dp(k+1))   / (dp(k-1) + 2.0*dp(k))
  enddo

  do i=1,ns
    !Compute average slopes: Colella, Eq. (1.8)
    as(1) = 0.
    do k=2,nk-1
      if (PCM_layer(k)) then  !use PCM
        as(k) = 0.0
      else
        slk = s(k,  i)-s(k-1,i)
        srk = s(k+1,i)-s(k,  i)
        if (slk*srk > 0.) then
          sck = h2_h123(k)*( h112(K)*srk*I_h12(K+1) + h122(K+1)*slk*I_h12(K) )
          as(k) = sign(min(abs(2.0*slk), abs(sck), abs(2.0*srk)), sck)
        else
          as(k) = 0.
        endif
      endif  !PCM:PPM
    enddo !k
    as(nk) = 0.
    !Compute "first guess" edge values: Colella, Eq. (1.6)
    al(1) = s(1,i)  ! 1st layer PCM
    ar(1) = s(1,i)  ! 1st layer PCM
    al(2) = s(1,i)  ! 1st layer PCM
    do K=3,nk-1
      ! This is a 4th order explicit edge value estimate.
      al(k) = (dp(k)*s(k-1,i) + dp(k-1)*s(k,i)) * I_h12(K) &
            + I_h0123(K)*( 2.*dp(k)*dp(k-1)*I_h12(K)*(s(k,i)-s(k-1,i)) * &
                           ( h01_h112(K) - h23_h122(K) ) &
                    + (dp(k)*as(k-1)*h23_h122(K) - dp(k-1)*as(k)*h01_h112(K)) )
     ar(k-1) = al(k)
    enddo !k
    ar(nk-1) = s(nk,i) ! last layer PCM
    al(nk)  = s(nk,i)  ! last layer PCM
    ar(nk)  = s(nk,i)  ! last layer PCM
    !Impose monotonicity: Colella, Eq. (1.10)
    do k=2,nk-1
      if ((PCM_layer(k)) .or. ((s(k+1,i)-s(k,i))*(s(k,i)-s(k-1,i)) <= 0.)) then !local extremum
        al(k) = s(k,i)
        ar(k) = s(k,i)
      else
        da = ar(k)-al(k)
        a6 = 6.0*s(k,i) - 3.0*(al(k)+ar(k))
        if (da*a6 > da*da) then !peak in right half of zone
          al(k) = 3.0*s(k,i) - 2.0*ar(k)
        elseif (da*a6 < -da*da) then !peak in left half of zone
          ar(k) = 3.0*s(k,i) - 2.0*al(k)
        endif
      endif
    enddo !k
    !Set coefficients
    do k=1,nk
      edges(k,1,i) = al(k)
      edges(k,2,i) = ar(k)
    enddo !k
  enddo !i

end subroutine hybgen_ppm_coefs


!> Set up the coefficients for PPM remapping of a set of scalars
subroutine hybgen_weno_coefs(s, h_src, edges, nk, ns, thin, PCM_lay)
  integer, intent(in)  :: nk        !< The number of input layers
  integer, intent(in)  :: ns        !< The number of scalar fields to work on
  real,    intent(in)  :: s(nk,ns)  !< The input scalar fields [A]
  real,    intent(in)  :: h_src(nk) !< The input grid layer thicknesses [H ~> m or kg m-2]
  real,    intent(out) :: edges(nk,2,ns) !< The WENO interpolation edge values of the scalar fields [A]
  real,    intent(in)  :: thin      !< A negligible layer thickness that can be ignored [H ~> m or kg m-2]
  logical, optional, intent(in)  :: PCM_lay(nk) !< If true for a layer, use PCM remapping for that layer

!-----------------------------------------------------------------------
!  1) coefficients for remapping from one set of vertical cells to another.
!     method: monotonic WENO-like alternative to PPM across each input cell
!             a second order polynomial approximation of the profiles
!             using a WENO reconciliation of the slopes to compute the
!             interfacial values
!
!     This scheme might have ben developed by Shchepetkin. A.F., personal communication.
!     See also Engwirda, D., and M. Kelley, A WENO-type slope-limiter for a family of piecewise
!       polynomial methods, arXive:1606.08188v1, 27 June 2016.
!
!  2) input arguments:
!       s     - initial scalar fields in pi-layer space
!       h_src - initial layer thicknesses (>=0)
!       nk    - number of layers
!       ns    - number of fields
!       thin  - layer thickness (>0) that can be ignored
!       PCM_lay - use PCM for selected layers (optional)
!
!  3) output arguments:
!       edges - cell edge scalar values for the WENO reconstruction
!                edges.1 is value at interface above
!                edges.2 is value at interface below
!
!  4) Laurent Debreu, Grenoble.
!     Alan J. Wallcraft,  Naval Research Laboratory,  July 2008.
!-----------------------------------------------------------------------
!
!  real, parameter :: dsmll=1.0e-8  ! This has units of [A2], and hence can not be a parameter.
!
  real :: curv_cell   ! An estimate of the tracer curvature centered on a cell times the grid
                      ! spacing [A H-1 ~> A m-1 or A kg m-2]
  real :: seh1, seh2  ! Tracer slopes at the cell edges times the cell grid spacing [A]
  real :: q01, q02    ! Various tracer differences between a cell average and the edge values [A]
  real :: q001, q002  ! Tracer slopes at the cell edges times the cell grid spacing [A]
  real :: ds2a, ds2b  ! Squared tracer differences between a cell average and the edge values [A2]
  logical :: PCM_layer(nk) ! True for layers that should use PCM remapping, either because they are
                      ! very thin, or because this is specified by PCM_lay.
  real :: dp(nk)      ! Input grid layer thicknesses, but with a minimum thickness given by thin [H ~> m or kg m-2]
  real :: qdpkm(nk)   ! Inverse of the sum of two adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: qdpkmkp(nk) ! Inverse of the sum of three adjacent thicknesses [H-1 ~> m-1 or m2 kg-1]
  real :: dpkm2kp(nk) ! Twice the distance between the centers of the layers two apart [H ~> m or kg m-2]
  real :: zw(nk,2)    ! Squared combinations of the differences between the the cell average tracer
                      ! concentrations and the left and right edges [A2]
  real :: min_ratio   ! The minimum ratio of the values of zw used to interpolate the edge values [nondim]
  real :: wt1         ! The weight of the upper layer in the interpolated shared edge value [nondim]
  real :: slope_edge(nk+1)  ! Tracer slopes at the edges [A H-1 ~> A m-1 or A kg m-2]
  real :: val_edge(nk+1)    ! A weighted average edge concentration [A]
  integer :: i, k

  min_ratio = 1.0e-8

  ! The WENO remapper is not currently written to work with massless layers, so set
  ! the thicknesses for very thin layers to some minimum value.
  do k=1,nk ; dp(k) = max(h_src(k), thin) ; enddo

  ! Specify the layers that will use PCM remapping.
  if (present(PCM_lay)) then
    do k=1,nk ; PCM_layer(k) = (PCM_lay(k) .or. dp(k) <= thin) ; enddo
  else
    do k=1,nk ; PCM_layer(k) = (dp(k) <= thin) ; enddo
  endif

  !compute grid metrics
  do k=2,nk-1
    qdpkm(  K) = 1.0 / (dp(k-1) + dp(k))
    qdpkmkp(k) = 1.0 / (dp(k-1) + dp(k) + dp(k+1))
    dpkm2kp(k) = dp(k-1) + 2.0*dp(k) + dp(k+1)
  enddo !k
  qdpkm(nk) = 1.0 / (dp(nk-1) + dp(nk))

  do i=1,ns
    do K=2,nk
      slope_edge(K) = qdpkm(K) * (s(k,i)-s(k-1,i))
    enddo !k
    k = 1  !PCM first layer
    edges(k,1,i) = s(k,i)
    edges(k,2,i) = s(k,i)
    zw(k,1) = 0.0
    zw(k,2) = 0.0
    do k=2,nk-1
      if ((slope_edge(K)*slope_edge(K+1) < 0.0) .or. PCM_layer(k)) then  !use PCM
        edges(k,1,i) = s(k,i)
        edges(k,2,i) = s(k,i)
        zw(k,1) = 0.0
        zw(k,2) = 0.0
      else
        seh1 = dp(k)*slope_edge(K+1)
        seh2 = dp(k)*slope_edge(K)
        q01 = dpkm2kp(k)*slope_edge(K+1)
        q02 = dpkm2kp(k)*slope_edge(K)
        if (abs(seh1) > abs(q02)) then
          seh1 = q02
        endif
        if (abs(seh2) > abs(q01)) then
          seh2 = q01
        endif
        curv_cell = (seh1 - seh2) * qdpkmkp(k)
        q001 = seh1 - curv_cell*dp(k+1)
        q002 = seh2 + curv_cell*dp(k-1)
        ! q001 = (seh1 * (dp(k-1) + dp(k)) + seh2 * dp(k+1)) * qdpkmkp(k)
        ! q002 = (seh2 * (dp(k+1) + dp(k)) + seh1 * dp(k-1)) * qdpkmkp(k)

        edges(k,2,i) = s(k,i) + q001
        edges(k,1,i) = s(k,i) - q002
        zw(k,1) = (2.0*q001 - q002)**2
        zw(k,2) = (2.0*q002 - q001)**2
      endif  !PCM:WENO
    enddo !k
    k = nk  !PCM last layer
    edges(k,1,i) = s(k,i)
    edges(k,2,i) = s(k,i)
    zw(k,  1) = 0.0
    zw(k,  2) = 0.0

    do k=2,nk
      ! This was the original code based on that in Hycom, but because zw has
      ! dimensions of [A2], it can not use a constant (hard coded) value of dsmll.
      !   ds2a = max(zw(k-1,2), dsmll)
      !   ds2b = max(zw(k,  1), dsmll)
      !   val_edge(K) = (ds2b*edges(k-1,2,i)+ds2a*edges(k,1,i)) / (ds2b+ds2a)
      ! Use a weighted average of the two layers' estimated edge values as the actual edge value.
      if (zw(k,1) + zw(k-1,2) <= 0.0) then
        wt1 = 0.5
      elseif (zw(k,1) <= min_ratio * (zw(k,1) + zw(k-1,2))) then
        wt1 = min_ratio
      elseif (zw(k-1,2) <= min_ratio * (zw(k,1) + zw(k-1,2))) then
        wt1 = (1.0 - min_ratio)
      else
        wt1 = zw(k,1) / (zw(k,1) + zw(k-1,2))
      endif
      val_edge(k) = wt1*edges(k-1,2,i) + (1.0-wt1)*edges(k,1,i)
    enddo !k
    val_edge(   1) = 2.0*s( 1,i)-val_edge( 2)  !not used?
    val_edge(nk+1) = 2.0*s(nk,i)-val_edge(nk)  !not used?

    do k=2,nk-1
      if (.not.PCM_layer(k)) then  !don't use PCM
        q01 = val_edge(K+1) - s(k,i)
        q02 = s(k,i) - val_edge(K)
        if (q01*q02 < 0.0) then
          q01 = 0.0
          q02 = 0.0
        elseif (abs(q01) > abs(2.0*q02)) then
          q01 = 2.0*q02
        elseif (abs(q02) > abs(2.0*q01)) then
          q02 = 2.0*q01
        endif
        edges(k,1,i) = s(k,i) - q02
        edges(k,2,i) = s(k,i) + q01
      endif  ! PCM:WENO
    enddo !k
  enddo !i

end subroutine hybgen_weno_coefs

end module MOM_hybgen_remap
