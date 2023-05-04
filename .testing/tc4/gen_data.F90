use netcdf
implicit none

integer, parameter :: dp = selected_real_kind(10, 100)
  !! Double precision (8-byte)

integer, parameter :: nz = 3
  !! Number of vertical layers
real(kind=dp), parameter :: salt0 = 35._dp
  !! Background salinity
real(kind=dp), parameter :: dampTime = 20._dp
  !! Sponge damping timescale [days]
real(kind=dp), parameter :: secs_per_day = 86400._dp
  !! Seconds per calendar day

integer :: ncid

integer :: x_id, y_id
integer :: lon_dimid, lat_dimid, depth_dimid, time_dimid
integer :: lon_id, lat_id, depth_id, time_id, temp_id, salt_id, idamp_id
integer :: field_dimids(2)
integer :: nx, ny

integer :: i, rc

real(kind=dp), allocatable :: x(:,:), y(:,:), z(:)
  !! Temperature grid positions
real(kind=dp), allocatable :: zbot(:,:)
  !! Bottom topography
real(kind=dp) :: zbot0
  !! Maximum topographic depth
real(kind=dp), allocatable :: temp(:,:,:), salt(:,:,:)
  !! Initial temperature and salinity fields
real(kind=dp), allocatable :: Idamp(:,:)
  !! Sponge dampening rate

! Read the domain grid
rc = nf90_open('ocean_hgrid.nc', NF90_NOWRITE, ncid)

rc = nf90_inq_varid(ncid, 'x', x_id)
rc = nf90_inq_varid(ncid, 'y', y_id)

rc = nf90_inquire_variable(ncid, x_id, dimids=field_dimids)
rc = nf90_inquire_dimension(ncid, field_dimids(1), len=nx)
rc = nf90_inquire_dimension(ncid, field_dimids(2), len=ny)

! Extract center ("T") points of supergrid
nx = nx / 2
ny = ny / 2
allocate(x(nx, ny), y(nx, ny))
rc = nf90_get_var(ncid, x_id, x, start=[2,2], stride=[2,2])
rc = nf90_get_var(ncid, y_id, y, start=[2,2], stride=[2,2])

rc = nf90_close(ncid)


! Read the topographic domain
rc = nf90_open('topog.nc', NF90_NOWRITE, ncid)

rc = nf90_inq_varid(ncid, 'depth', depth_id)
rc = nf90_inquire_variable(ncid, depth_id, dimids=field_dimids)
rc = nf90_inquire_dimension(ncid, field_dimids(1), len=nx)
rc = nf90_inquire_dimension(ncid, field_dimids(2), len=ny)

allocate(zbot(nx, ny))
rc = nf90_get_var(ncid, depth_id, zbot)
rc = nf90_close(ncid)


! Construct the vertical axis
allocate(z(nz))
z = [(i, i=0,nz-1)] * maxval(zbot) / nz

allocate(temp(nx, ny, nz), salt(nx, ny, nz))
call t_fc(x, y, z, temp)
salt(:,:,:) = salt0


! Write T/S initial state
rc = nf90_create('temp_salt_ic.nc', NF90_CLOBBER, ncid)

rc = nf90_def_dim(ncid, 'lon', nx, lon_dimid)
rc = nf90_def_dim(ncid, 'lat', ny, lat_dimid)
rc = nf90_def_dim(ncid, 'depth', nz, depth_dimid)
rc = nf90_def_dim(ncid, 'Time', NF90_UNLIMITED, time_dimid)

rc = nf90_def_var(ncid, 'depth', NF90_DOUBLE, [depth_dimid], depth_id)
rc = nf90_def_var(ncid, 'lon', NF90_DOUBLE, [lon_dimid], lon_id)
rc = nf90_def_var(ncid, 'lat', NF90_DOUBLE, [lat_dimid], lat_id)
rc = nf90_def_var(ncid, 'Time', NF90_DOUBLE, [time_dimid], time_id)

rc = nf90_put_att(ncid, time_id, 'calendar', 'noleap')
rc = nf90_put_att(ncid, time_id, 'units', 'days since 0001-01-01 00:00:00.0')
! NOTE: nf90_put_att() truncates empty strings, so use nf90_put_att_any()
rc = nf90_put_att_any(ncid, time_id, 'modulo', NF90_CHAR, 1, ' ')

rc = nf90_def_var(ncid, 'ptemp', NF90_DOUBLE, &
    [lon_dimid, lat_dimid, depth_dimid, time_dimid], temp_id)
rc = nf90_def_var_fill(ncid, temp_id, 0, -1e20_dp)

rc = nf90_def_var(ncid, 'salt', NF90_DOUBLE, &
    [lon_dimid, lat_dimid, depth_dimid, time_dimid], salt_id)
rc = nf90_def_var_fill(ncid, salt_id, 0, -1e20_dp)

rc = nf90_enddef(ncid)

rc = nf90_put_var(ncid, lon_id, x(:,1))
rc = nf90_put_var(ncid, lat_id, y(1,:))
rc = nf90_put_var(ncid, depth_id, z)
rc = nf90_put_var(ncid, time_id, 0.)
rc = nf90_put_var(ncid, temp_id, temp)
rc = nf90_put_var(ncid, salt_id, salt)

rc = nf90_close(ncid)


! Sponge file
rc = nf90_create('sponge.nc', NF90_CLOBBER, ncid)

rc = nf90_def_dim(ncid, 'lon', nx, lon_dimid)
rc = nf90_def_dim(ncid, 'lat', ny, lat_dimid)

rc = nf90_def_var(ncid, 'lon', NF90_DOUBLE, lon_id)
rc = nf90_def_var(ncid, 'lat', NF90_DOUBLE, lat_id)
rc = nf90_def_var(ncid, 'Idamp', NF90_DOUBLE, [lon_dimid, lat_dimid], Idamp_id)
rc = nf90_def_var_fill(ncid, Idamp_id, 0, -1e20_dp)

rc = nf90_enddef(ncid)

allocate(Idamp(nx, ny))
Idamp = 0.
if (dampTime > 0.) &
  Idamp(:,:) = 1. / (dampTime * secs_per_day)

rc = nf90_put_var(ncid, Idamp_id, Idamp)
rc = nf90_put_var(ncid, lon_id, x(:,1))
rc = nf90_put_var(ncid, lat_id, y(1,:))

rc = nf90_close(ncid)

contains

subroutine t_fc(x, y, z, tl, radius, tmag)
  real(kind=dp), intent(in) :: x(:,:), y(:,:), z(:)
    !! Grid positions
  real(kind=dp), intent(inout) :: tl(:,:,:)
    !! Temperature field on the model grid
  real(kind=dp), intent(in), optional :: radius
    !! Temperature anomaly radius
  real(kind=dp), intent(in), optional :: tmag
    !! Temperature anomaly maximum

  real(kind=dp) :: t_rad, t_max
    !! Temperature field parameters (radius, max value)
  real(kind=dp) :: x0, y0
    !! Center of anomaly (currently midpoint of domain)
  real(kind=dp), allocatable :: r(:,:), zd(:)
    !! Radial and vertical extent of anomaly
  integer :: k, nz
    !! Vertical level indexing

  t_rad = 5._dp
  if (present(radius)) t_rad = radius

  t_max = 1._dp
  if (present(tmag)) t_max = tmag

  ! Reduce supergrid size to T/S grid
  allocate(zd, source=z)

  x0 = x(1 + size(x, 1)/2, 1 + size(x, 2)/2)
  y0 = y(1 + size(y, 1)/2, 1 + size(y, 2)/2)

  tl(:,:,:) = 0.
  nz = size(z)
  if (nz > 1) then
    zd(:) = z(:) / z(nz)
  else
    zd(:) = 0.
  endif

  allocate(r, source=x)
  r(:,:) = hypot(x(:,:) - x0, y(:,:) - y0)
  do k = 1, nz
    tl(:,:,k) = (1. - min(r(:,:) / t_rad, 1.)) * t_max * (1. - zd(k))
  enddo
end subroutine t_fc

end
