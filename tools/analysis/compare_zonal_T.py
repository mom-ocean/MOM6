#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting change in annual-average zonal temperature.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'temp' and 'e'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the title of the main plot.''')
parser.add_argument('-1','--label1', type=str, default='', help='''Lable to give experiment 1.''')
parser.add_argument('-2','--label2', type=str, default='', help='''Lable to give experiment 2.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-r','--ref', type=str, required=True,
  help='''File containing reference experiment to compare against.''')
cmdLineArgs = parser.parse_args()

y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][1::2,1::2].max(axis=-1)
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
basin = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/basin_codes.nc').variables['basin'][:]

rootGroupRef = netCDF4.Dataset( cmdLineArgs.ref )
if 'temp' in rootGroupRef.variables: varName = 'temp'
elif 'ptemp' in rootGroupRef.variables: varName = 'ptemp'
elif 'thetao' in rootGroupRef.variables: varName = 'thetao'
else: raise Exception('Could not find "temp", "ptemp" or "thetao" in file "%s"'%(cmdLineArgs.ref))
if len(rootGroupRef.variables[varName].shape)==4: Tref = rootGroupRef.variables[varName][:].mean(axis=0)
else: Tref = rootGroupRef.variables[varName][:]
if 'e' in rootGroupRef.variables: Zref = rootGroupRef.variables['e'][0]
else:
  D = -netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_topog.nc').variables['depth'][:]
  zw = -rootGroupRef.variables['zw'][:]
  Zref = numpy.zeros((Tref.shape[0]+1,Tref.shape[1],Tref.shape[2]))
  for k in range(Tref.shape[0]):
    Zref[k+1] = numpy.maximum( Zref[k] - abs(zw[k+1] - zw[k]), D)

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'temp' in rootGroup.variables: varName = 'temp'
elif 'ptemp' in rootGroup.variables: varName = 'ptemp'
elif 'thetao' in rootGroup.variables: varName = 'thetao'
else: raise Exception('Could not find "temp", "ptemp" or "thetao" in file "%s"'%(cmdLineArgs.annual_file))
if len(rootGroup.variables[varName].shape)==4: Tmod = rootGroup.variables[varName][:].mean(axis=0)
else: Tmod = rootGroup.variables[varName][:]
if 'e' in rootGroup.variables: Zmod = rootGroup.variables['e'][0]
else: Zmod = Zref

def zonalAverage(T, eta, area, mask=1.):
  vols = ( mask * area ) * ( eta[:-1] - eta[1:] ) # mask * area * level thicknesses
  return numpy.sum( vols * T, axis=-1 ) / numpy.sum( vols, axis=-1 ), (mask*eta).min(axis=-1)

ci=m6plot.pmCI(0.025,.45,.05)
if len(cmdLineArgs.label1): title1 = cmdLineArgs.label1
else: title1 = rootGroup.title
if len(cmdLineArgs.label2): title2 = cmdLineArgs.label2
else: title2 = rootGroupRef.title

# Global
tPlot, z = zonalAverage(Tmod, Zmod, area)
tRefPlot, _ = zonalAverage(Tref, Zref, area)
m6plot.yzplot( tPlot - tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title=r'''Global zonal-average $\theta$ response [$\degree$C] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/T_global_xave_response.png')

m6plot.yzcompare( tPlot, tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=r'Global zonal-average $\theta$ [$\degree$C] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(-2,29,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/T_global_xave_response.3_panel.png')

# Atlantic + Arctic
newMask = 1.*msk; newMask[ (basin!=2) & (basin!=4) ] = 0.
tPlot, z = zonalAverage(Tmod, Zmod, area, mask=newMask)
tRefPlot, _ = zonalAverage(Tref, Zref, area, mask=newMask)
m6plot.yzplot( tPlot - tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title=r'''Atlantic zonal-average $\theta$ response [$\degree$C] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/T_Atlantic_xave_response.png')

m6plot.yzcompare( tPlot, tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=r'Atlantic zonal-average $\theta$ response [$\degree$C] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(-2,29,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/T_Atlantic_xave_response.3_panel.png')

# Pacific
newMask = 1.*msk; newMask[ (basin!=3) ] = 0.
tPlot, z = zonalAverage(Tmod, Zmod, area, mask=newMask)
tRefPlot, _ = zonalAverage(Tref, Zref, area, mask=newMask)
m6plot.yzplot( tPlot - tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title=r'''Pacific zonal-average $\theta$ response [$\degree$C] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/T_Pacific_xave_response.png')

m6plot.yzcompare( tPlot, tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=r'Pacific zonal-average $\theta$ response [$\degree$C] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(-2,29,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/T_Pacific_xave_response.3_panel.png')

# Indian
newMask = 1.*msk; newMask[ (basin!=5) ] = 0.
tPlot, z = zonalAverage(Tmod, Zmod, area, mask=newMask)
tRefPlot, _ = zonalAverage(Tref, Zref, area, mask=newMask)
m6plot.yzplot( tPlot - tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=title1+' - '+title2, title=r'''Indian zonal-average $\theta$ response [$\degree$C] '''+cmdLineArgs.label,
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/T_Indian_xave_response.png')

m6plot.yzcompare( tPlot, tRefPlot , y, z, splitscale=[0., -1000., -6500.],
      suptitle=r'Indian zonal-average $\theta$ response [$\degree$C] '+cmdLineArgs.label,
      title1=title1, title2=title2,
      clim=m6plot.linCI(-2,29,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/T_Indian_xave_response.3_panel.png')
