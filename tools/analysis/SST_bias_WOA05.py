#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting annual-average SST bias.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'temp'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-w','--woa', type=str, required=True,
  help='''File containing WOA (or obs) data to compare against.''')
cmdLineArgs = parser.parse_args()

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'temp' not in rootGroup.variables: raise Exception('Could not find "temp" in file "%s"'%(cmdLineArgs.annual_file))

x = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
msk = numpy.ma.array(msk, mask=(msk==0))

Tobs = netCDF4.Dataset( cmdLineArgs.woa ).variables['temp'][0]

variable = rootGroup.variables['temp']
if variable.shape[0]>1: temp = variable[:,0].mean(axis=0)
else: temp = variable[0,0]

ci=m6plot.pmCI(0.25,4.5,.5)
m6plot.xyplot( temp - Tobs , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='SST bias (w.r.t. WOA\'05) [$\degree$C]',
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/SST_bias_WOA05.png')

m6plot.xycompare( temp, Tobs , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label,
      title1='SST [$\degree$C]',
      title2='WOA\'05 SST [$\degree$C]',
      clim=m6plot.linCI(-2,29,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/SST_bias_WOA05.3_panel.png')
