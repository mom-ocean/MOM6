import xarray as xr
import numpy as np

# tc5
ds = xr.Dataset()
ds['h2'] = xr.DataArray(np.zeros((10,10)), dims=('y', 'x'))
ds.to_netcdf('h2_tc5.nc')

