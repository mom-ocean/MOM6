#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting difference in annual-average SSS between two experiments.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'salt'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-1','--label1', type=str, default='', help='''Lable to give experiment 1.''')
parser.add_argument('-2','--label2', type=str, default='', help='''Lable to give experiment 2.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-r','--ref', type=str, required=True,
  help='''File containing reference experiment to compare against.''')
cmdLineArgs = parser.parse_args()

x = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
msk = numpy.ma.array(msk, mask=(msk==0))

rootGroupRef = netCDF4.Dataset( cmdLineArgs.ref )
if 'salt' in rootGroupRef.variables: varName='salt'
elif 'so' in rootGroupRef.variables: varName='so'
else: raise Exception('Could not find "salt" or "so" in file "%s"'%(cmdLineArgs.ref))
if rootGroupRef.variables[varName].shape[0]>1: Sref = rootGroupRef.variables[varName][:,0].mean(axis=0)
else: Sref = rootGroupRef.variables[varName][0,0]

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'salt' in rootGroup.variables: varName='salt'
elif 'so' in rootGroup.variables: varName='so'
else: raise Exception('Could not find "salt" or "so" in file "%s"'%(cmdLineArgs.annual_file))
if len(rootGroup.variables[varName].shape)==4: Smod = rootGroup.variables[varName][:,0].mean(axis=0)
else: Smod = rootGroup.variables[varName][0]

if len(cmdLineArgs.label1): title1 = cmdLineArgs.label1
else: title1 = rootGroup.title
if len(cmdLineArgs.label2): title2 = cmdLineArgs.label2
else: title2 = rootGroupRef.title

ci=m6plot.pmCI(0.025,0.525,.05)
m6plot.xyplot( Smod - Sref , x, y, area=area,
      suptitle=title1+' - '+title2,
      title='SSS difference [ppt] '+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/SSS_difference.png')

m6plot.xycompare( Smod, Sref , x, y, area=area,
      suptitle='SSS difference [ppt] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/SSS_difference.3_panel.png')
