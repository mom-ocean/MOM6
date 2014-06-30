#!/usr/bin/env python

import netCDF4
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

t = netCDF4.Dataset('prog.nc').variables['Time'][:]
T = netCDF4.Dataset('prog.nc').variables['temp'][:,:,0,0]
e = netCDF4.Dataset('prog.nc').variables['e'][:,:,0,0]
K = netCDF4.Dataset('visc.nc').variables['Kd_interface'][:,:,0,0]
mld = netCDF4.Dataset('visc.nc').variables['MLD_003'][:,0,0]
N2 = netCDF4.Dataset('visc.nc').variables['KPP_N2'][:,:,0,0]

nk=63
plt.subplot(2,2,1)
cmap = plt.get_cmap();
norm = BoundaryNorm(numpy.arange(16,29.5,.5), ncolors=cmap.N)
plt.pcolormesh(t+0*e[:,:nk+1].T, e[:,:nk+1].T, T[:,:nk].T, cmap=cmap, norm=norm); plt.xlim(0,365); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$\theta$ [$\degree$C]')
plt.subplot(2,2,2)
plt.contourf(t+0*e[:,:nk].T, e[:,:nk].T, numpy.log10( numpy.maximum(K[:,:nk].T,1e-8) ), levels=numpy.arange(-5,.25,.25) ); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$log_{10}|\kappa$| [m$^2$s$^{-1}$]')
plt.subplot(2,2,3)
norm = BoundaryNorm(numpy.arange(-5.,-1.,.25), ncolors=cmap.N)
plt.pcolormesh(t+0*e[:,:nk+1].T, e[:,:nk+1].T, 0.5*numpy.log10( numpy.maximum(N2[:,:nk].T,0e-12) ), cmap=cmap, norm=norm); plt.xlim(0,365); plt.ylim(-300,0); plt.colorbar()
#plt.contourf(t+0*e[:,:nk].T, e[:,:nk].T, 0.5*numpy.log10( numpy.maximum(N2[:,:nk].T,0e-12) ), levels=numpy.arange(-5.,-1.,.25) ); plt.ylim(-300,0); plt.colorbar()
plt.plot(t,-mld,color='m',hold=True)
plt.ylabel('z [m]'); plt.title(r'$log_{10}|N|$ [s$^{-1}$]')
plt.xlabel('Time [days]')
plt.subplot(2,2,4)
plt.plot(t,T[:,0]); plt.ylim(18,30)
plt.ylabel('$\theta$ [$\degree$C]'); plt.title(r'SST [s$\degree$C]')
plt.xlabel('Time [days]')

plt.tight_layout()
plt.savefig('BML.png')

