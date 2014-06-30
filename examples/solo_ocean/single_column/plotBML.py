#!/usr/bin/env python

import netCDF4
import matplotlib.pyplot as plt
import numpy

t = netCDF4.Dataset('prog.nc').variables['Time'][:]
T = netCDF4.Dataset('prog.nc').variables['temp'][:,:,0,0]
e = netCDF4.Dataset('prog.nc').variables['e'][:,:,0,0]
K = netCDF4.Dataset('visc.nc').variables['Kd_interface'][:,:,0,0]
mld = netCDF4.Dataset('visc.nc').variables['MLD_003'][:,0,0]

nk=63
plt.subplot(2,1,1)
plt.contourf(t+0*e[:,:nk].T, e[:,:nk].T, T[:,:nk].T, levels=numpy.arange(18.,28,.5)); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$\theta$ [$\degree$C]')
plt.subplot(2,1,2)
plt.contourf(t+0*e[:,:nk].T, e[:,:nk].T, numpy.log10( numpy.maximum(K[:,:nk].T,1e-8) ), levels=numpy.arange(-5,.25,.25) ); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$\kappa$ [m$^2$s$^{-1}$]')
plt.xlabel('Time [days]')
plt.savefig('BML.png')

