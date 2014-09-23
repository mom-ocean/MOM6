#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting difference in annual-average MLD between two experiments.''')
parser.add_argument('monthly_file', type=str, help='''Monthly-averaged file containing 2D 'MLD_003'.''')
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
if 'MLD_003' not in rootGroupRef.variables: raise Exception('Could not find "MLD_003" in file "%s"'%(cmdLineArgs.monthly_file))
variable = rootGroupRef.variables['MLD_003']
shape = variable.shape
MLDref = variable[:].reshape(shape[0]/12,12,shape[1],shape[2]).mean(axis=0)

rootGroup = netCDF4.Dataset( cmdLineArgs.monthly_file )
if 'MLD_003' not in rootGroup.variables: raise Exception('Could not find "MLD_003" in file "%s"'%(cmdLineArgs.monthly_file))
variable = rootGroup.variables['MLD_003']
shape = variable.shape
MLDmod = variable[:].reshape(shape[0]/12,12,shape[1],shape[2]).mean(axis=0)

if len(cmdLineArgs.label1): title1 = cmdLineArgs.label1
else: title1 = rootGroup.title
if len(cmdLineArgs.label2): title2 = cmdLineArgs.label2
else: title2 = rootGroupRef.title

ciMin = m6plot.linCI(0,95,5)
ciMax = m6plot.linCI(0,680,20)
diMin = m6plot.pmCI(2.5,32.5,5)
diMax = m6plot.pmCI(0,20,5, 40,100,20)
diMax = numpy.array([-250,-200,-150,-100,-50,-20,-10,-3,3,10,20,50,100,150,200,250])

m6plot.xyplot( MLDmod.min(axis=0) - MLDref.min(axis=0), x, y, area=area,
      suptitle=title1+' - '+title2,
      title='MLD difference [m] '+cmdLineArgs.label,
      clim=diMin, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/MLD_minimum_difference.png')

m6plot.xycompare( MLDmod.min(axis=0), MLDref.min(axis=0), x, y, area=area,
      suptitle='MLD difference [m] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=ciMin, colormap='dunneRainbow', extend='max',
      dlim=diMin, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/MLD_minimum_difference.3_panel.png')

m6plot.xyplot( MLDmod.max(axis=0) - MLDref.max(axis=0), x, y, area=area,
      suptitle=title1+' - '+title2,
      title='MLD difference [m] '+cmdLineArgs.label,
      clim=diMax, colormap='dunnePM', centerlabels=False, extend='both',
      save=cmdLineArgs.outdir+'/MLD_maximum_difference.png')

m6plot.xycompare( MLDmod.max(axis=0), MLDref.max(axis=0), x, y, area=area,
      suptitle='MLD difference [m] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=ciMax, colormap='dunneRainbow', extend='max',
      dlim=diMax, dcolormap='dunnePM', dextend='both', centerdlabels=False,
      save=cmdLineArgs.outdir+'/MLD_maximum_difference.3_panel.png')
