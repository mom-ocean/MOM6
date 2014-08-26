#!/usr/bin/env python

import netCDF4
import numpy
import m6plot
import math

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting annual-average SST bias.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file contained 3D temp.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspec', type=str,
  default='/archive/gold/datasets/MOM6z_SIS_025/siena/mosaic.unpacked',
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
cmdLineArgs = parser.parse_args()

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'ssu' not in rootGroup.variables: raise Exception('Could not find "ssu" in file "%s"'%(cmdLineArgs.annual_file))
if 'ssv' not in rootGroup.variables: raise Exception('Could not find "ssv" in file "%s"'%(cmdLineArgs.annual_file))

x = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
msk = numpy.ma.array(msk, mask=(msk==0))

#[t,z,y,x] corresponds to axis [0,1,2,3] which can be indexed by [-4,-3,-2,-1]

ssu = rootGroup.variables['ssu']
ssu_mean = ssu[:].mean(axis=0)
eke_u = (0.5*(ssu-ssu_mean)**2).mean(axis=0)
eke = (eke_u + numpy.roll(eke_u,1,axis=-1))/2 # U-point to T-point transform

ssv = rootGroup.variables['ssv']
ssv_mean = ssv[:].mean(axis=0)
eke_v = (0.5*(ssv-ssv_mean)**2).mean(axis=0)
eke = eke + (eke_v + numpy.roll(eke_v,1,axis=-2))/2 

ci=m6plot.pmCI(0.0,0.5,0.1)

m6plot.xyplot( 10000*eke, x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Eddy Kinetic Energy annual mean [(cm/s)^2]',
      clim=ci, logscale=True,
      save=cmdLineArgs.outdir+'/EKE_mean.png')

