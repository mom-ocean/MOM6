#!/usr/bin/env python

import netCDF4
import numpy
import m6plot
import matplotlib.pyplot as plt

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting annual-min/max mixed layer depth.''')
parser.add_argument('monthly_file', type=str, help='''Monthly-averaged file containing at least 12 months of MLD_003.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-od','--obsdata', type=str,
  default='/archive/gold/datasets/obs/Hosada2010_MLD_climatology.v20140515.nc',
  help='''File containing the observational MLD data (Hosoda et al., 2010).''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
cmdLineArgs = parser.parse_args()

rootGroup = netCDF4.Dataset( cmdLineArgs.monthly_file )
if 'MLD_003' not in rootGroup.variables: raise Exception('Could not find "MLD_003" in file "%s"'%(cmdLineArgs.monthly_file))

x = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)

variable = rootGroup.variables['MLD_003']
shape = variable.shape
MLD = variable[:].reshape(shape[0]/12,12,shape[1],shape[2]).mean(axis=0)

MLD_obs = netCDF4.Dataset(cmdLineArgs.obsdata).variables['MLD'][:]
x_obs = netCDF4.Dataset(cmdLineArgs.obsdata).variables['LONGITUDE'][:]
y_obs = netCDF4.Dataset(cmdLineArgs.obsdata).variables['LATITUDE'][:]

ciMin = m6plot.linCI(0,95,5)
ciMax = m6plot.linCI(0,680,20)

# Plot of shallowest model MLD (summer)
m6plot.xyplot( MLD.min(axis=0), x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Annual-minimum MLD$_{0.03}$ [m]',
      clim=ciMin, extend='max', colormap='dunneRainbow',
      save=cmdLineArgs.outdir+'/MLD_003_minimum.png')

# 2-panel plot of shallowest model MLD + obs (summer)
m6plot.setFigureSize(aspect=[3,3], verticalresolution=976, npanels=0)
ax1 = plt.subplot(2,1,1)
m6plot.xyplot( numpy.roll(MLD_obs.min(axis=0),300,axis=-1), x_obs-300, y_obs,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Hosoda et al., 2010, annual-minimum MLD$_{0.03}$ [m]',
      clim=ciMin, extend='max', colormap='dunneRainbow',
      axis=ax1)
ax2 = plt.subplot(2,1,2)
m6plot.xyplot( MLD.min(axis=0), x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Annual-minimum MLD$_{0.03}$ [m]',
      clim=ciMin, extend='max', colormap='dunneRainbow',
      axis=ax2,
      save=cmdLineArgs.outdir+'/MLD_003_minimum.2_panel.png')

# Plot of deepest model MLD (winter)
m6plot.xyplot( MLD.max(axis=0), x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Annual-maximum MLD$_{0.03}$ [m]',
      clim=ciMax, extend='max', colormap='dunneRainbow',
      save=cmdLineArgs.outdir+'/MLD_003_maximum.png')

# 2-panel plot of deepest model MLD + obs (winter)
m6plot.setFigureSize(aspect=[3,3], verticalresolution=976, npanels=0)
ax1 = plt.subplot(2,1,1)
m6plot.xyplot( numpy.roll(MLD_obs.max(axis=0),300,axis=-1), x_obs-300, y_obs,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Hosoda et al., 2010, annual-maximum MLD$_{0.03}$ [m]',
      clim=ciMax, extend='max', colormap='dunneRainbow',
      axis=ax1)
ax2 = plt.subplot(2,1,2)
m6plot.xyplot( MLD.max(axis=0), x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label, title='Annual-maximum MLD$_{0.03}$ [m]',
      clim=ciMax, extend='max', colormap='dunneRainbow',
      axis=ax2,
      save=cmdLineArgs.outdir+'/MLD_003_maximum.2_panel.png')
