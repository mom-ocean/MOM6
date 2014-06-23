#!/usr/bin/env python

import netCDF4
import matplotlib.pyplot as plt
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting annual-average SST bias.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file contained 3D salt.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspec', type=str,
  default='/archive/gold/datasets/MOM6z_SIS_025/siena/mosaic.unpacked',
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
cmdLineArgs = parser.parse_args()

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'salt' not in rootGroup.variables: raise Exception('Could not find "salt" in file "%s"'%(cmdLineArgs.annual_file))

x = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
msk = numpy.ma.array(msk, mask=(msk==0))

Tobs = netCDF4.Dataset( '/archive/gold/datasets/MOM6z_SIS_025/siena/obs/WOA05_ptemp_salt_annual_v1.nc' ).variables['salt'][1]

variable = rootGroup.variables['salt']
if variable.shape[0]>1: salt = variable[:,0].mean(axis=0)
else: salt = variable[0,0]

ci=m6plot.pmCI(0.125,2.25,.25)
m6plot.xyplot( salt - Tobs , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='SST bias (w.r.t. WOA\'05) [ppt]',
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/SSS_bias_WOA05.png')

m6plot.xycompare( salt, Tobs , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label,
      title1='SSS [ppt]',
      title2='WOA\'05 SSS [ppt]',
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='both',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/SSS_bias_WOA05.3_panel.png')
