#!/usr/bin/env python
#
# Written by Matthew Harrison (matthew.harrison@noaa.gov)
# Calculate depth-integrated heat and salt content
# over a limited depth range.  

from midas.rectgrid import *
from midas.rectgrid_gen import *
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import argparse

parser = argparse.ArgumentParser(description='''Script for generating heat and salt content over a specified depth range.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'temp'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-s','--start_depth', type=float, default='0', help='''Starting Depth for analysis.''')
parser.add_argument('-e','--end_depth', type=float, default='-1', help='''Ending Depth for analysis.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-r','--reference_data', type=str, default='.', help='''The reference data containing temperature and salinity on the model grid.''')

cmdLineArgs = parser.parse_args()


path_ann=cmdLineArgs.annual_file
path_obs=cmdLineArgs.reference_data

sgrid=supergrid(file=cmdLineArgs.gridspecdir+'/ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.D=nc.Dataset(cmdLineArgs.gridspecdir+'/ocean_topog.nc').variables['depth'][:]
grid.lath=grid.y_T[:,grid.im/4]
grid.latq=grid.y_T_bounds[:,grid.im/4]
grid.wet=np.zeros(grid.D.shape)
grid.wet[grid.D>0.]=1.

S=state(path_ann,grid=grid,fields=['temp','salt'],interfaces='e',verbose=False)
O=state(path_obs,grid=grid,fields=['ptemp','salt'],interfaces='eta',verbose=False)

O.rename_field('ptemp','temp')

max_depth=numpy.max(-S.var_dict['temp']['z_interfaces'][:,-1,:])
S.grid.D[S.grid.D>max_depth]=max_depth
O.grid.D[O.grid.D>max_depth]=max_depth

if cmdLineArgs.end_depth == -1:
    cmdLineArgs.end_depth = max_depth

if cmdLineArgs.end_depth > max_depth:
    cmdLineArgs.end_depth = max_depth    

if cmdLineArgs.start_depth > cmdLineArgs.end_depth:
    print "Error in command arguments. Start depth must be less then end depth"
    raise
    
if cmdLineArgs.start_depth > 0:
    zbax=numpy.array([0.,-cmdLineArgs.start_depth,-cmdLineArgs.end_depth,-max_depth])
else:
    zbax=numpy.array([0.,-cmdLineArgs.end_depth,-max_depth])
    
zb=zbax[:,numpy.newaxis,numpy.newaxis]
zb=numpy.tile(zb,(1,grid.jm,grid.im))
zi=S.var_dict['temp']['z_interfaces']
zb[-1,:]=zi[0,-1,:]
zb[1,:]=numpy.maximum(zi[0,-1,:],-300.)

S.remap_ALE(['temp','salt'],z_bounds=zb,zbax_data=-zbax,method='pcm')
S.adjust_thickness('temp_remap')
S.adjust_thickness('salt_remap')
S.rename_field('temp_remap','temp')
S.rename_field('salt_remap','salt')

O.remap_ALE(['temp','salt'],z_bounds=zb,zbax_data=-zbax,method='pcm')
O.adjust_thickness('temp_remap')
O.adjust_thickness('salt_remap')
O.rename_field('temp_remap','temp')
O.rename_field('salt_remap','salt')

vdict=S.var_dict['temp'].copy()
S.add_field_from_array(O.temp,'temp_obs',var_dict=vdict)
vdict=S.var_dict['salt'].copy()
S.add_field_from_array(O.salt,'salt_obs',var_dict=vdict)


S.write_nc('tmp.nc',['temp','salt','temp_obs','salt_obs'],write_interface_positions=True)


if cmdLineArgs.start_depth > 0:
    z_indices=numpy.array([1])
else:
    z_indices=numpy.array([0])
    
cp = 3989.0;rho0=1.035e3
S=state('tmp.nc',fields=['temp','salt','temp_obs','salt_obs'],grid=grid,z_indices=z_indices,interfaces='eta',verbose=False)


S.volume_integral('temp','Z',normalize=False)
S.temp_zint = S.temp_zint * rho0 * cp*1.e-10
S.rename_field('temp_zint','hc')
S.var_dict['hc']['units']='10^10 Joules'
S.volume_integral('temp_obs','Z',normalize=False)
S.temp_obs_zint = S.temp_obs_zint * rho0 * cp * 1.e-10
S.rename_field('temp_obs_zint','hc_obs')
S.var_dict['hc_obs']['units']='10^10 Joules'
S.volume_integral('salt','Z',normalize=False)
S.salt_zint = S.salt_zint*rho0/1.e6
S.rename_field('salt_zint','sc')
S.var_dict['sc']['units']='Mg Salt'
S.volume_integral('salt_obs','Z',normalize=False)
S.salt_obs_zint = S.salt_obs_zint*rho0/1.e6
S.rename_field('salt_obs_zint','sc_obs')
S.var_dict['sc_obs']['units']='Mg Salt'
vdict=S.var_dict['hc'].copy()
S.add_field_from_array(S.hc-S.hc_obs,'hc_bias',var_dict=vdict)
vdict=S.var_dict['sc'].copy()
S.add_field_from_array(S.sc-S.sc_obs,'sc_bias',var_dict=vdict)
fnam_out=cmdLineArgs.outdir+'/'+cmdLineArgs.label+'_heat_salt_'+str(int(cmdLineArgs.start_depth))+'_'+str(int(cmdLineArgs.end_depth))+'.nc'
S.write_nc(fnam_out,['hc','sc','hc_bias','sc_bias'],write_interface_positions=True)


fig=plt.figure(1,figsize=(8.5,11))
ax1=fig.add_subplot(211)
ax2=fig.add_subplot(212)
cf=ax1.contourf(S.grid.x_T,S.grid.y_T,numpy.squeeze(S.hc_bias),np.arange(-1,1.02,.02),extend='both',cmap=plt.cm.bwr)
ax1.contour(S.grid.x_T,S.grid.y_T,S.grid.wet,[0.5,0.5],colors='k')
tit=cmdLineArgs.label+' Heat Content (10^10 Joules) Bias z= '+str(cmdLineArgs.start_depth)+' to '+str(cmdLineArgs.end_depth)
ax1.set_title(tit,fontsize=10)
plt.colorbar(cf,ax=ax1)
cf=ax2.contourf(S.grid.x_T,S.grid.y_T,numpy.squeeze(S.sc_bias),np.arange(-0.5,0.55,.05),extend='both',cmap=plt.cm.bwr)
ax2.contour(S.grid.x_T,S.grid.y_T,S.grid.wet,[0.5,0.5],colors='k')
tit=cmdLineArgs.label+' Salt Content (10^6 gSalt) Bias z= '+str(cmdLineArgs.start_depth)+' to '+str(cmdLineArgs.end_depth)
ax2.set_title(tit,fontsize=10)
plt.colorbar(cf,ax=ax2)
fnam_out=cmdLineArgs.outdir+'/'+cmdLineArgs.label+'_heat_salt_'+str(int(cmdLineArgs.start_depth))+'_'+str(int(cmdLineArgs.end_depth))+'.png'
plt.savefig(fnam_out)
