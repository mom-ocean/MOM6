#!/usr/bin/env python

from midas import *
import netCDF4 as nc
import numpy as np

f=nc.Dataset('scap_supergrid.nc')
f2=nc.Dataset('antarctic_spherical_supergrid.nc')
g=nc.Dataset('mercator_supergrid.nc')
h=nc.Dataset('ncap_supergrid.nc')

y1=f.variables['y'][:]
y12=f2.variables['y'][:]
y2=g.variables['y'][:]
y3=h.variables['y'][:]
y=np.concatenate((y1,y12[1:,:]),axis=0)
y=np.concatenate((y,y2[1:,:]),axis=0)
y=np.concatenate((y,y3[1:,:]),axis=0)

dy1=f.variables['dy'][:]
dy12=f2.variables['dy'][:]
dy2=g.variables['dy'][:]
dy3=h.variables['dy'][:]
dy=np.concatenate((dy1,dy12),axis=0)
dy=np.concatenate((dy,dy2),axis=0)
dy=np.concatenate((dy,dy3),axis=0)

x1=f.variables['x'][:]
x12=f2.variables['x'][:]
x2=g.variables['x'][:]
x3=h.variables['x'][:]
x=np.concatenate((x1,x12[1:,:]),axis=0)
x=np.concatenate((x,x2[1:,:]),axis=0)
x=np.concatenate((x,x3[1:,:]),axis=0)

dx1=f.variables['dx'][:]
dx12=f2.variables['dx'][:]
dx2=g.variables['dx'][:]
dx3=h.variables['dx'][:]
dx=np.concatenate((dx1,dx12[1:,:]),axis=0)
dx=np.concatenate((dx,dx2[1:,:]),axis=0)
dx=np.concatenate((dx,dx3[1:,:]),axis=0)


area1=f.variables['area'][:]
area12=f2.variables['area'][:]
area2=g.variables['area'][:]
area3=h.variables['area'][:]
area=np.concatenate((area1,area12),axis=0)
area=np.concatenate((area,area2),axis=0)
area=np.concatenate((area,area3),axis=0)

angle_dx1=f.variables['angle_dx'][:-1,:]
angle_dx12=f2.variables['angle_dx'][:-1,:]
angle_dx2=g.variables['angle_dx'][:-1,:]
angle_dx3=h.variables['angle_dx'][:]
angle_dx=np.concatenate((angle_dx1,angle_dx12),axis=0)
angle_dx=np.concatenate((angle_dx,angle_dx2),axis=0)
angle_dx=np.concatenate((angle_dx,angle_dx3),axis=0)

fout=nc.Dataset('supergrid.nc','w',format='NETCDF3_CLASSIC')

ny=area.shape[0]; nx = area.shape[1]
nyp=ny+1; nxp=nx+1

print 'ny,nx= ',ny,nx

nyp=fout.createDimension('nyp',nyp)
nxp=fout.createDimension('nxp',nxp)
ny=fout.createDimension('ny',ny)
nx=fout.createDimension('nx',nx)
string=fout.createDimension('string',255)    

tile=fout.createVariable('tile','S1',('string'))
yv=fout.createVariable('y','f8',('nyp','nxp'))
xv=fout.createVariable('x','f8',('nyp','nxp'))    
yv.units='degrees'
xv.units='degrees'
yv[:]=y
xv[:]=x

tile[0:4]='tile1'
dyv=fout.createVariable('dy','f8',('ny','nxp'))
dyv.units='meters'
dyv[:]=dy
dxv=fout.createVariable('dx','f8',('nyp','nx'))
dxv.units='meters'
dxv[:]=dx
areav=fout.createVariable('area','f8',('ny','nx'))
areav.units='m2'
areav[:]=area
anglev=fout.createVariable('angle_dx','f8',('nyp','nxp'))
anglev.units='degrees'
anglev[:]=angle_dx            

fout.sync()
fout.close()
