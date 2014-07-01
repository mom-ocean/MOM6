#!/usr/bin/env python

import netCDF4
import numpy
import sys
sys.path.append('../../../tools/analysis/')
import m6plot
import matplotlib.pyplot as plt

def subplot( path, axis, title ):
  #x = netCDF4.Dataset(path+'/Initial_state.nc').variables['lonh'][:]
  #e = netCDF4.Dataset(path+'/Initial_state.nc').variables['eta'][:,0,:]
  #S = netCDF4.Dataset(path+'/Initial_state.nc').variables['Salt'][0,:,0,:]
  x = netCDF4.Dataset(path+'/prog.nc').variables['xh'][:]
  e = netCDF4.Dataset(path+'/prog.nc').variables['e'][-1,:,0,:]
  S = netCDF4.Dataset(path+'/prog.nc').variables['salt'][-1,:,0,:]
  m6plot.yzplot(S, x, e, axis=axis, title=title, clim=m6plot.linCI(34,35,.05))
  plt.plot(x, e.T, hold=True, color='k')

ax = plt.subplot(2,2,1); subplot( 'layer', ax, 'Layer' )
ax = plt.subplot(2,2,2); subplot( 'z', ax, 'z*' )
ax = plt.subplot(2,2,3); subplot( 'sigma', ax, r'$\sigma$' )
ax = plt.subplot(2,2,4); subplot( 'rho', ax, r'$\rho$' )

plt.tight_layout()
plt.savefig('flow_downslope.png')
