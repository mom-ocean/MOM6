use netcdf

implicit none

integer, parameter :: dp = selected_real_kind(10, 100)
  !! Double precision (8-byte)

integer, parameter :: nx = 14, ny = 10
  !! Grid size
real(kind=dp), parameter :: depth0 = 100._dp
  !! Uniform depth
real(kind=dp), parameter :: ds = 0.01_dp
  !! Grid resolution at the equator in degrees
real(kind=dp), parameter :: Re = 6.378e6_dp
  !! Radius of earth
real(kind=dp), parameter :: rad_per_deg = (4. * atan(1._dp)) / 180._dp
  !! Degress to radians (= pi/180.)

integer :: ncid
integer :: nx_id, ny_id, nxp_id, nyp_id, ntile_id, string_id
integer :: depth_id, h2_id
integer :: x_id, y_id, dx_id, dy_id, area_id, angle_id, tile_id

! Fields on model grid
real(kind=dp) :: depth(nx, ny)

! Grid fields (defined on supergrid)
real(kind=dp) :: xg(0:2*nx), yg(0:2*ny)
real(kind=dp) :: x(0:2*nx, 0:2*ny), y(0:2*nx, 0:2*ny)
real(kind=dp) :: dx(0:2*nx-1, 0:2*ny)
real(kind=dp) :: dy(0:2*nx, 0:2*ny-1)
real(kind=dp) :: area(0:2*nx-1, 0:2*ny-1)
real(kind=dp) :: angle_dx(0:2*nx, 0:2*ny)

integer :: i, j, rc


! Topography
rc = nf90_create('topog.nc', NF90_CLOBBER, ncid)

rc = nf90_def_dim(ncid, 'ny', ny, ny_id)
rc = nf90_def_dim(ncid, 'nx', nx, nx_id)
rc = nf90_def_dim(ncid, 'ntiles', 1, ntile_id)

rc = nf90_def_var(ncid, 'depth', NF90_DOUBLE, [nx_id, ny_id], depth_id)
rc = nf90_def_var(ncid, 'h2', NF90_DOUBLE, [nx_id, ny_id], h2_id)

rc = nf90_enddef(ncid)

depth(:,:) = depth0
rc = nf90_put_var(ncid, depth_id, depth)

rc = nf90_close(ncid)


! Horizontal grid
rc = nf90_create('ocean_hgrid.nc', NF90_CLOBBER, ncid)

rc = nf90_def_dim(ncid, 'ny', 2*ny, ny_id)
rc = nf90_def_dim(ncid, 'nx', 2*nx, nx_id)
rc = nf90_def_dim(ncid, 'nyp', 2*ny+1, nyp_id)
rc = nf90_def_dim(ncid, 'nxp', 2*nx+1, nxp_id)
rc = nf90_def_dim(ncid, 'string', 5, string_id)

rc = nf90_def_var(ncid, 'y', NF90_DOUBLE, [nxp_id, nyp_id], y_id)
rc = nf90_def_var(ncid, 'x', NF90_DOUBLE, [nxp_id, nyp_id], x_id)
rc = nf90_def_var(ncid, 'dy', NF90_DOUBLE, [nxp_id, ny_id], dy_id)
rc = nf90_def_var(ncid, 'dx', NF90_DOUBLE, [nx_id, nyp_id], dx_id)
rc = nf90_def_var(ncid, 'area', NF90_DOUBLE, [nx_id, ny_id], area_id)
rc = nf90_def_var(ncid, 'angle_dx', NF90_DOUBLE, [nxp_id, nyp_id], angle_id)
rc = nf90_def_var(ncid, 'tile', NF90_CHAR, string_id, tile_id)

rc = nf90_put_att(ncid, y_id, 'units', 'degrees')
rc = nf90_put_att(ncid, x_id, 'units', 'degrees')
rc = nf90_put_att(ncid, dy_id, 'units', 'meters')
rc = nf90_put_att(ncid, dx_id, 'units', 'meters')
rc = nf90_put_att(ncid, area_id, 'units', 'm2')
rc = nf90_put_att(ncid, angle_id, 'units', 'degrees')

rc = nf90_enddef(ncid)

xg = ds * [(i, i=0, 2*nx)]
yg = ds * [(j, j=0, 2*ny)]

! NOTE: sin() and cos() are compiler-dependent

x(:,:) = spread(xg(:), 2, 2*ny+1)
y(:,:) = spread(yg(:), 1, 2*nx+1)
dx(:,:) = rad_per_deg * Re * (x(1:,:) - x(:2*nx-1,:)) &
    * cos(0.5 * rad_per_deg * (y(1:,:) + y(:2*nx-1,:)))
dy(:,:) = rad_per_deg * Re * (y(:,1:) - y(:,:2*ny-1))

area(:,:) = rad_per_deg * Re * Re &
  * spread(sin(rad_per_deg * yg(1:)) - sin(rad_per_deg * yg(:2*ny-1)), 1, 2*nx) &
  * spread(xg(1:) - xg(:2*nx-1), 2, 2*ny)

angle_dx(:,:) = 0.

rc = nf90_put_var(ncid, x_id, x)
rc = nf90_put_var(ncid, y_id, y)
rc = nf90_put_var(ncid, dx_id, dx)
rc = nf90_put_var(ncid, dy_id, dy)
rc = nf90_put_var(ncid, area_id, area)
rc = nf90_put_var(ncid, angle_id, angle_dx)
rc = nf90_put_var(ncid, tile_id, 'tile1')

rc = nf90_close(ncid)
end
