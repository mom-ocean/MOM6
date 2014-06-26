#!/usr/bin/env python

import netCDF4
import matplotlib.pyplot as plt
import numpy

z = netCDF4.Dataset('prog.nc').variables['zl'][:]
zi = netCDF4.Dataset('prog.nc').variables['zi'][:]
t = netCDF4.Dataset('prog.nc').variables['Time'][:]
T = netCDF4.Dataset('prog.nc').variables['temp'][:,:,0,0]
h = netCDF4.Dataset('visc.nc').variables['KPP_OBLdepth'][:,0,0]
K = netCDF4.Dataset('visc.nc').variables['KPP_Kheat'][:,:,0,0]
mld = netCDF4.Dataset('visc.nc').variables['MLD_003'][:,0,0]

nk=27
plt.subplot(2,1,1)
plt.contourf(t, -z[:nk], T[:,:nk].T, levels=numpy.arange(18.,28,.5)); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-h,color='k',hold=True)
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$\theta$ [$\degree$C]')
plt.subplot(2,1,2)
plt.contourf(t, -zi[:nk], numpy.log10( numpy.maximum(K[:,:nk].T,1e-8) ), levels=numpy.arange(-5,.25,.25) ); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-h,color='k',hold=True)
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$\kappa$ [m$^2$s$^{-1}$]')
plt.xlabel('Time [days]')
plt.savefig('KPP.png')

