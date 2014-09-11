#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting change in annual-average zonal salinity.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'salt' and 'e'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to title of the main plot.''')
parser.add_argument('-1','--label1', type=str, default='', help='''Label to give experiment 1.''')
parser.add_argument('-2','--label2', type=str, default='', help='''Label to give experiment 2.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-r','--ref', type=str, required=True,
  help='''File containing reference experiment to compare against.''')
cmdLineArgs = parser.parse_args()

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'salt' not in rootGroup.variables: raise Exception('Could not find "salt" in file "%s"'%(cmdLineArgs.annual_file))
if 'e' not in rootGroup.variables: raise Exception('Could not find "e" in file "%s"'%(cmdLineArgs.annual_file))

rootGroupRef = netCDF4.Dataset( cmdLineArgs.ref )
if 'salt' not in rootGroupRef.variables: raise Exception('Could not find "salt" in file "%s"'%(cmdLineArgs.ref))
if 'e' not in rootGroupRef.variables: raise Exception('Could not find "e" in file "%s"'%(cmdLineArgs.ref))

y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][1::2,1::2].max(axis=-1)
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
basin = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/basin_codes.nc').variables['basin'][:]

def zonalAverage(T, eta, area, mask=1.):
  vols = ( mask * area ) * ( eta[:-1] - eta[1:] ) # mask * area * level thicknesses
  return numpy.sum( vols * T, axis=-1 ) / numpy.sum( vols, axis=-1 ), (mask*eta).min(axis=-1)

if rootGroupRef.variables['salt'].shape[0]>1:
  Sref = rootGroupRef.variables['salt'][:,:].mean(axis=0)
  Zref = rootGroupRef.variables['e'][0]
else:
  Sref = rootGroupRef.variables['salt'][0]
  Zref = rootGroupRef.variables['e'][0]

if rootGroup.variables['salt'].shape[0]>1:
  Smod = rootGroup.variables['salt'][:,:].mean(axis=0)
  Zmod = rootGroup.variables['e'][0]
else:
  Smod = rootGroup.variables['salt'][0]
  Zmod = rootGroup.variables['e'][0]

ci=m6plot.pmCI(0.0125,.225,.025)
if len(cmdLineArgs.label1): title1 = cmdLineArgs.label1
else: title1 = rootGroup.title
if len(cmdLineArgs.label2): title2 = cmdLineArgs.label2
else: title2 = rootGroupRef.title

# Global
sPlot, z = zonalAverage(Smod, Zmod, area)
sRefPlot, _ = zonalAverage(Sref, Zref, area)
m6plot.yzplot( sPlot - sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title='''Global zonal-average salinity response [ppt] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/S_global_xave_response.png')

m6plot.yzcompare( sPlot, sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle='Global zonal-average salinity [ppt] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='both',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/S_global_xave_response.3_panel.png')

# Atlantic + Arctic
newMask = 1.*msk; newMask[ (basin!=2) & (basin!=4) ] = 0.
sPlot, z = zonalAverage(Smod, Zmod, area, mask=newMask)
sRefPlot, _ = zonalAverage(Sref, Zref, area, mask=newMask)
m6plot.yzplot( sPlot - sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title='''Atlantic zonal-average salinity response [ppt] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/S_Atlantic_xave_response.png')

m6plot.yzcompare( sPlot, sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle='Atlantic zonal-average salinity [ppt] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='both',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/S_Atlantic_xave_response.3_panel.png')

# Pacific
newMask = 1.*msk; newMask[ (basin!=3) ] = 0.
sPlot, z = zonalAverage(Smod, Zmod, area, mask=newMask)
sRefPlot, _ = zonalAverage(Sref, Zref, area, mask=newMask)
m6plot.yzplot( sPlot - sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title='''Pacific zonal-average salinity response [ppt] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/S_Pacific_xave_response.png')

m6plot.yzcompare( sPlot, sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle='Pacific zonal-average salinity [ppt] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='both',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/S_Pacific_xave_response.3_panel.png')

# Indian
newMask = 1.*msk; newMask[ (basin!=5) ] = 0.
sPlot, z = zonalAverage(Smod, Zmod, area, mask=newMask)
sRefPlot, _ = zonalAverage(Sref, Zref, area, mask=newMask)
m6plot.yzplot( sPlot - sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title='''Indian zonal-average salinity response [ppt] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/S_Indian_xave_response.png')

m6plot.yzcompare( sPlot, sRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle='Indian zonal-average salinity [ppt] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(20,30,10, 31,39,.5), colormap='dunneRainbow', extend='both',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/S_Indian_xave_response.3_panel.png')
