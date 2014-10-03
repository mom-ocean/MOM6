#!/usr/bin/env python

from midas.rectgrid import *
import netCDF4 as nc
import numpy as np

f=nc.Dataset('DATA/grid_tpxo7_atlas.nc')

lon_u=f.variables['lon_u'][:].T
lat_u=f.variables['lat_u'][:].T
lon_v=f.variables['lon_v'][:].T
lat_v=f.variables['lat_v'][:].T

grid_u=quadmesh(lon=lon_u,lat=lat_u,cyclic=True)
grid_v=quadmesh(lon=lon_v,lat=lat_v,cyclic=True)

Su=state(grid=grid_u)
Sv=state(grid=grid_v)

f=nc.Dataset('DATA/u_tpxo7_atlas.nc')

ua=np.zeros((1,8,721,1440))
va=np.zeros((1,8,721,1440))

for i in np.arange(0,8):
    tmp=f.variables['ua'][i,:]
    ua[0,i,:]=tmp.T
    tmp=f.variables['va'][i,:]
    va[0,i,:]=tmp.T    

ua[ua==0.0]=-1.e20
va[va==0.0]=-1.e20
ua[:,:,-1,:] = ua[:,:,-2,:] # Work around for missing data at N-pole
va[:,:,-1,:] = va[:,:,-2,:] # Work around for missing data at N-pole
ua=np.ma.masked_where(ua==-1.e20,ua)
va=np.ma.masked_where(va==-1.e20,va)


vdict_u={}
vdict_v={}

vdict_u['X']=f.variables['ua'].dimensions[1]
vdict_u['Y']=f.variables['ua'].dimensions[2]
vdict_u['Z']=f.variables['ua'].dimensions[0][0:8]
vdict_u['T']=None
vdict_u['units']='cm s-1'
vdict_u['path']='DATA/u_tpxo7_atlas.nc'
vdict_u['Zdir']=1
vdict_u['Ztype']='Fixed'
vdict_u['Zb']=None
vdict_u['z']=np.arange(0,8)
vdict_u['z_interfaces']=None
vdict_u['zunits']='none'
vdict_u['zax_data']=np.arange(0,8)
vdict_u['xax_data']=grid_u.lonh
vdict_u['yax_data']=grid_u.lath
vdict_u['xunits']='degrees_east'
vdict_u['yunits']='degrees_north'
vdict_u['zunits']='component'
vdict_u['_FillValue']=-1.e20
vdict_u['missing_value']=-1.e20
vdict_u['masked']=True
vdict_v['X']=f.variables['va'].dimensions[1]
vdict_v['Y']=f.variables['va'].dimensions[2]
vdict_v['Z']=f.variables['va'].dimensions[0][0:8]
vdict_v['T']=None
vdict_v['units']='cm s-1'
vdict_v['path']='DATA/u_tpxo7_atlas.nc'
vdict_v['Zdir']=1
vdict_v['Ztype']='Fixed'
vdict_v['Zb']=None
vdict_v['z']=np.arange(0,8)
vdict_v['z_interfaces']=None
vdict_v['zunits']='none'
vdict_v['zax_data']=np.arange(0,8)
vdict_v['xax_data']=grid_v.lonh
vdict_v['yax_data']=grid_v.lath
vdict_v['xunits']='degrees_east'
vdict_v['yunits']='degrees_north'
vdict_v['zunits']='component'
vdict_v['_FillValue']=-1.e20
vdict_v['missing_value']=-1.e20
vdict_v['masked']=True

Su.add_field_from_array(ua,'ua',var_dict=vdict_u)
Sv.add_field_from_array(va,'va',var_dict=vdict_v)

grid=quadmesh(lon=lon_v,lat=lat_u,cyclic=True)
grid.wet=np.ones((grid.jm,grid.im))

S=Su.horiz_interp('ua',target=grid,src_modulo=True,method='bilinear')
S=Sv.horiz_interp('va',target=grid,src_modulo=True,method='bilinear',PrevState=S)

u2mod = (S.ua**2.0 + S.va**2.0)
umod=np.sum(u2mod,axis=1)**0.5
umod=umod[:,np.newaxis,:]
umod = 1.e-2*umod

vdict=S.var_dict['ua'].copy()
vdict['units']='m s-1'
vdict['Z'] = None

S.add_field_from_array(umod,'umod',var_dict=vdict)

sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
output_grid = quadmesh(supergrid=sgrid,cyclic=True)
output_grid.D=nc.Dataset('ocean_topog.nc').variables['depth'][:]
output_grid.wet = np.zeros(output_grid.D.shape)
output_grid.wet[output_grid.D>0.]=1.

S.fill_interior('umod')

T=S.horiz_interp('umod',target=output_grid,src_modulo=True,method='bilinear')

T.fill_interior('umod')

T.var_dict['umod']['yax_data'] = output_grid.y_T[:,output_grid.im/4]

#T.umod = np.ma.filled(T.umod,0.0)6.n

T.rename_field('umod','tideamp')

T.write_nc('tidal_amplitude.nc',['tideamp'])
