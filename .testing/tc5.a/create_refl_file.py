import xarray as xr
import matplotlib.pyplot as plt
import numpy as np


def create_refl_walls(ds, nx=10, ny=10, rotate=False):
    """ create reflection for E-W walls """
    
    refl = 0 * np.ones((ny, nx))
    refl_angle = -999.9 * np.ones((ny, nx))
    refl_pref =  0* np.ones((ny, nx))
    refl_dbl = 0 * np.ones((ny, nx))
    trans = 0 * np.ones((ny, nx))

    jsouth=0
    jnorth=-1
    #jsouth=2
    #jnorth=-3

    refl[jsouth,:] = 1.
    refl[jnorth,:] = 1.
    refl_pref[jsouth,:] = 1.
    refl_pref[jnorth,:] = 1.
    refl_angle[jsouth,:] = np.pi ##I believe this is correct
    refl_angle[jnorth,:] = 0.
    refl_dbl[jsouth,:] = 0.
    refl_dbl[jnorth,:] = 0.
    trans[jsouth,:] = 0.
    trans[jnorth,:] = 0.

    if rotate:
        refl = refl.transpose()
        refl_pref = refl_pref.transpose()
        refl_angle = refl_angle.transpose()
        refl_dbl = refl_dbl.transpose()
        trans = trans.transpose()

    ds['refl'] = xr.DataArray(refl, dims=('y', 'x'), attrs = {'_FillValue': 1e+20})
    ds['refl_angle'] = xr.DataArray(refl_angle, dims=('y','x'), attrs = {'_FillValue': 1e+20})
    ds['refl_pref'] = xr.DataArray(refl_pref, dims=('y','x'), attrs = {'_FillValue': 1e+20})
    ds['refl_dbl'] = xr.DataArray(refl_dbl, dims=('y','x'), attrs = {'_FillValue': 1e+20})
    ds['trans'] = xr.DataArray(trans, dims=('y','x'), attrs = {'_FillValue': 1e+20})
    return ds


#--------------------- tc5 --------------------------------
ds = xr.Dataset()
ds = create_refl_walls(ds)
ds.to_netcdf('IWcoefs_tc5.nc')
