import netCDF4 as nc
from netCDF4 import stringtochar
import numpy as np

nx, ny = 14, 10     # Grid size
depth0 = 100.       # Uniform depth
ds = 0.01           # grid resolution at the equator in degrees
Re = 6.378e6        # Radius of earth

topo_ = np.zeros((ny, nx)) + depth0
f_topo = nc.Dataset('topog.nc', 'w', format='NETCDF3_CLASSIC')
ny, nx = topo_.shape
f_topo.createDimension('ny', ny)
f_topo.createDimension('nx', nx)
f_topo.createDimension('ntiles', 1)
f_topo.createVariable('depth', 'f8', ('ny', 'nx'))
f_topo.createVariable('h2', 'f8', ('ny', 'nx'))
f_topo.variables['depth'][:] = topo_
f_topo.sync()
f_topo.close()

x_ = np.arange(0, 2*nx + 1) * ds  # units are degrees E
y_ = np.arange(0, 2*ny + 1) * ds  # units are degrees N
x, y = np.meshgrid(x_, y_)

dx = np.zeros((2*ny + 1, 2*nx))
dy = np.zeros((2*ny, 2*nx + 1))
rad_deg = np.pi / 180.
dx[:] = (rad_deg * Re * (x[:, 1:] - x[:, 0:-1])
         * np.cos(0.5*rad_deg*(y[:, 0:-1] + y[:, 1:])))
dy[:] = rad_deg * Re * (y[1:, :] - y[0:-1, :])

f_sg = nc.Dataset('ocean_hgrid.nc', 'w', format='NETCDF3_CLASSIC')
f_sg.createDimension('ny', 2*ny)
f_sg.createDimension('nx', 2*nx)
f_sg.createDimension('nyp', 2*ny + 1)
f_sg.createDimension('nxp', 2*nx + 1)
f_sg.createDimension('string', 5)
f_sg.createVariable('y', 'f8', ('nyp', 'nxp'))
f_sg.createVariable('x', 'f8', ('nyp', 'nxp'))
dyv = f_sg.createVariable('dy', 'f8', ('ny', 'nxp'))
dxv = f_sg.createVariable('dx', 'f8', ('nyp', 'nx'))
areav = f_sg.createVariable('area', 'f8', ('ny', 'nx'))
dxv.units = 'm'
dyv.units = 'm'
areav.units = 'm2'
f_sg.createVariable('angle_dx', 'f8', ('nyp', 'nxp'))
f_sg.createVariable('tile', 'S1', ('string'))
f_sg.variables['y'].units = 'degrees'
f_sg.variables['x'].units = 'degrees'
f_sg.variables['dy'].units = 'meters'
f_sg.variables['dx'].units = 'meters'
f_sg.variables['area'].units = 'm2'
f_sg.variables['angle_dx'].units = 'degrees'
f_sg.variables['y'][:] = y
f_sg.variables['x'][:] = x
f_sg.variables['dx'][:] = dx
f_sg.variables['dy'][:] = dy

# Compute the area bounded by lines of constant
# latitude-longitud on a sphere in m2.
dlon = x_[1:] - x_[:-1]
dlon = np.tile(dlon[np.newaxis, :], (2*ny, 1))
y1_ = y_[:-1]
y1_ = y1_[:, np.newaxis]*rad_deg
y2_ = y_[1:]
y2_ = y2_[:, np.newaxis]*rad_deg
y1_ = np.tile(y1_, (1, 2*nx))
y2_ = np.tile(y2_, (1, 2*nx))
area = rad_deg * Re * Re * (np.sin(y2_) - np.sin(y1_)) * dlon
f_sg.variables['area'][:] = area
f_sg.variables['angle_dx'][:] = 0.
str_ = stringtochar(np.array(['tile1'], dtype='S5'))
f_sg.variables['tile'][:] = str_
f_sg.sync()
f_sg.close()
