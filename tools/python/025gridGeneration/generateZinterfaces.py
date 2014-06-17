import netCDF4
import numpy

dz=netCDF4.Dataset('../../../examples/ocean_SIS/MOM6z_SIS_025/INPUT/vgrid_cm4.nc').variables['dz'][:]
zi=numpy.zeros(76)
zi[1:]=numpy.cumsum(-dz)
print zi

Zbot = -netCDF4.Dataset('ocean_topog.nc').variables['depth'][:]

Zi =numpy.zeros((76,1080,1440))
for k in range(0,76):
  Zi[k,:,:] = numpy.maximum( Zbot, zi[k] )
