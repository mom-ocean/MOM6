#!/usr/bin/env python

from midas.rectgrid import *
import netCDF4 as nc
import numpy as np

sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.lath=grid.y_T[:,grid.im/4]
grid.latq=grid.y_T_bounds[:,grid.im/4]
grid.D=nc.Dataset('ocean_topog.nc').variables['depth'][:]
grid.wet=np.zeros(grid.D.shape)
grid.wet[grid.D>0.]=1
S=state(grid=grid)

# Model vertical grid
#dz=nc.Dataset('../../../examples/ocean_SIS/OM4_025/INPUT/vgrid_75_2m.nc').variables['dz'][:]
#nk = dz.shape[0]
#zi=np.zeros(nk+1)
#zi[1:]=np.cumsum(-dz)
# Analysis vertical  grid
zi=-nc.Dataset('../../../examples/ocean_SIS/OM4_025/INPUT/vgrid_75_2m.nc').variables['zw'][:]
nk=zi.shape[0]-1

zb =np.zeros((nk+1,S.grid.jm,S.grid.im))
for k in range(0,nk+1):
  zb[k,:]=zi[k]
  zb[k,:,:] = np.maximum( -S.grid.D, zb[k] )

for n in np.arange(0,12):
   O=state('/archive/gold/datasets/obs/WOA05_pottemp_salt.nc',fields=['SALT','PTEMP'],time_indices=np.arange(n,n+1),default_calendar='noleap',z_orientation=-1)
   O.grid.cyclic_x=True
   O.rename_field('PTEMP','ptemp')
   O.rename_field('SALT','salt')
   OM=O.horiz_interp('salt',target=S.grid,method='bilinear')
   OM=O.horiz_interp('ptemp',target=S.grid,method='bilinear',PrevState=OM)
   OM.adjust_thickness('ptemp')
   OM.adjust_thickness('salt')
   OM.fill_interior('salt',smooth=True,num_pass=10000)
   OM.fill_interior('ptemp',smooth=True,num_pass=10000)

   OM.remap_ALE(fields=['ptemp','salt'],z_bounds=zb,zbax_data=-zi,method='ppm_h4',bndy_extrapolation=False)
   OM.rename_field('ptemp_remap','ptemp')
   OM.rename_field('salt_remap','salt')
   OM.mask_where('ptemp','grid.wet==0.')
   OM.mask_where('salt','grid.wet==0.')
   OM.ptemp=np.ma.masked_where(OM.var_dict['ptemp']['dz'][np.newaxis,:]<1.e-2, OM.ptemp)
   OM.salt=np.ma.masked_where(OM.var_dict['ptemp']['dz'][np.newaxis,:]<1.e-2, OM.salt)

   if n==0:
      OM.write_nc('WOA05_ptemp_salt_monthly.nc',['ptemp','salt'],append=False,write_interface_positions=True)
   else:
      OM.write_nc('WOA05_ptemp_salt_monthly.nc',['ptemp','salt'],append=True,write_interface_positions=True)
