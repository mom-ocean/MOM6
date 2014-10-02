#!/usr/bin/env python

from midas.rectgrid import *
import numpy as np
import netCDF4 as nc

sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
O=state('/archive/gold/datasets/obs/WOA05_pottemp_salt.nc',fields=['SALT'],z_indices=np.arange(0,1),default_calendar='noleap')
O.grid.cyclic_x=True
O.var_dict['SALT']['Z'] = None
OM=O.horiz_interp('SALT',target=grid,method='bilinear')
OM.grid.D=nc.Dataset('ocean_topog.nc').variables['depth'][:]
OM.grid.wet=np.zeros(OM.grid.D.shape)
OM.grid.wet[OM.grid.D>0.]=1.
OM.fill_interior('SALT',smooth=True,num_pass=10000)
OM.mask_where('SALT','grid.D<=0.')
OM.rename_field('SALT','salt')
OM.var_dict['salt']['xax_data']=grid.x_T[0,:]
OM.var_dict['salt']['yax_data']=grid.y_T[:,grid.im/4]
OM.write_nc('salt_restore.nc',fields=['salt'])
