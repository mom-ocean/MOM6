#!/usr/bin/env python

import netCDF4
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

z = netCDF4.Dataset('prog.nc').variables['zl'][:]
zi = netCDF4.Dataset('prog.nc').variables['zi'][:]
t = netCDF4.Dataset('prog.nc').variables['Time'][:]
T = netCDF4.Dataset('prog.nc').variables['temp'][:,:,0,0]
e = netCDF4.Dataset('prog.nc').variables['e'][:,:,0,0]
h = netCDF4.Dataset('visc.nc').variables['KPP_OBLdepth'][:,0,0]
K = netCDF4.Dataset('visc.nc').variables['KPP_Kheat'][:,:,0,0]
mld = netCDF4.Dataset('visc.nc').variables['MLD_003'][:,0,0]
N2 = netCDF4.Dataset('visc.nc').variables['KPP_N2'][:,:,0,0]

nk=27
cmap = plt.get_cmap()
plt.subplot(2,2,1)
plt.contourf(t, -z[:nk], T[:,:nk].T, levels=numpy.arange(18.,28,.5)); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-h,color='k',hold=True)
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$\theta$ [$\degree$C]')
plt.subplot(2,2,2)
plt.contourf(t, -zi[:nk], numpy.log10( numpy.maximum(K[:,:nk].T,1e-8) ), levels=numpy.arange(-5,.25,.25) ); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-h,color='k',hold=True)
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$log_{10}|\kappa|$ [m$^2$s$^{-1}$]')
plt.subplot(2,2,3)
norm = BoundaryNorm(numpy.arange(-5.,-1.,.25), ncolors=cmap.N)
#plt.pcolormesh(t+0*e[:,:nk+1].T, e[:,:nk+1].T, 0.5*numpy.log10( numpy.maximum(N2[:,:nk].T,0e-12) ), cmap=cmap, norm=norm); plt.xlim(0,365); plt.ylim(-300,0); plt.colorbar()
plt.contourf(t, -zi[:nk], 0.5*numpy.log10( numpy.maximum(N2[:,:nk].T,0e-12) ), levels=numpy.arange(-5.,-1.,.25) ); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-h,color='k',hold=True)
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$log_{10}|N|$ [s$^{-1}$]')
plt.xlabel('Time [days]')
plt.subplot(2,2,4)
plt.plot(t,T[:,0]); plt.ylim(18,30)
plt.ylabel(r'$\theta$ [$\degree$C]'); plt.title('SST')
plt.xlabel('Time [days]')

plt.tight_layout()
plt.savefig('KPP.png')
