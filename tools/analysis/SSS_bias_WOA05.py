#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting annual-average SSS bias.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'salt'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-w','--woa', type=str, required=True,
  help='''File containing WOA (or obs) data to compare against.''')
cmdLineArgs = parser.parse_args()

x = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
msk = numpy.ma.array(msk, mask=(msk==0))

Sobs = netCDF4.Dataset( cmdLineArgs.woa ).variables['salt']
if len(Sobs.shape)==3: Sobs = Sobs[0]
else: Sobs = Sobs[:,0].mean(axis=0)

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'salt' in rootGroup.variables: varName = 'salt'
elif 'so' in rootGroup.variables: varName = 'so'
else: raise Exception('Could not find "salt" or "so" in file "%s"'%(cmdLineArgs.annual_file))
if rootGroup.variables[varName].shape[0]>1: salt = rootGroup.variables[varName][:,0].mean(axis=0)
else: Smod = rootGroup.variables[varName][0,0]

ci=m6plot.pmCI(0.125,2.25,.25)
m6plot.xyplot( Smod - Sobs , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='SSS bias (w.r.t. WOA\'05) [ppt]',
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/SSS_bias_WOA05.png')

m6plot.xycompare( Smod, Sobs , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label,
      title1='SSS [ppt]',
      title2='WOA\'05 SSS [ppt]',
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='both',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/SSS_bias_WOA05.3_panel.png')
