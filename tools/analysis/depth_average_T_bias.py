#!/usr/bin/env python

import netCDF4
import numpy
import m6plot

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description='''Script for plotting depth-average temperature bias.''')
parser.add_argument('annual_file', type=str, help='''Annually-averaged file containing 3D 'temp' and 'e'.''')
parser.add_argument('-l','--label', type=str, default='', help='''Label to add to the plot.''')
parser.add_argument('-o','--outdir', type=str, default='.', help='''Directory in which to place plots.''')
parser.add_argument('-g','--gridspecdir', type=str, required=True,
  help='''Directory containing mosaic/grid-spec files (ocean_hgrid.nc and ocean_mask.nc).''')
parser.add_argument('-w','--woa', type=str, required=True,
  help='''File containing WOA (or obs) data to compare against.''')
parser.add_argument('-zt','--ztop', type=float, default=0., help='''Depth (+ve) over which to start the average (default 300).''')
parser.add_argument('-zb','--zbottom', type=float, default=300., help='''Depth (+ve) over which to end the average (default 300).''')
cmdLineArgs = parser.parse_args()

x = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['x'][::2,::2]
y = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['y'][::2,::2]
msk = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_mask.nc').variables['mask'][:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_hgrid.nc').variables['area'][:,:].reshape([msk.shape[0], 2, msk.shape[1], 2]).sum(axis=-3).sum(axis=-1)
depth = netCDF4.Dataset(cmdLineArgs.gridspecdir+'/ocean_topog.nc').variables['depth'][:]
depth = numpy.ma.array(depth, mask=depth<=cmdLineArgs.ztop)
lDepth = cmdLineArgs.zbottom
uDepth = cmdLineArgs.ztop

Tobs = netCDF4.Dataset( cmdLineArgs.woa )
if 'temp' in Tobs.variables: Tobs = Tobs.variables['temp']
else: Tobs = Tobs.variables['ptemp']
if len(Tobs.shape)==3: Tobs = Tobs[:].filled(0.)
else: Tobs = Tobs[:].mean(axis=0).filled(0.)
Zobs = netCDF4.Dataset( cmdLineArgs.woa ).variables['eta'][:]

rootGroup = netCDF4.Dataset( cmdLineArgs.annual_file )
if 'temp' in rootGroup.variables: varName = 'temp'
elif 'ptemp' in rootGroup.variables: varName = 'ptemp'
elif 'thetao' in rootGroup.variables: varName = 'thetao'
else: raise Exception('Could not find "temp", "ptemp" or "thetao" in file "%s"'%(cmdLineArgs.annual_file))
if len(rootGroup.variables[varName].shape)==4: Tmod = rootGroup.variables[varName][:].mean(axis=0).filled(0.)
else: Tmod = rootGroup.variables[varName][:].filled(0.)
if 'e' in rootGroup.variables: Zmod = rootGroup.variables['e'][0]
else: Zmod = Zobs # Using model z-ou:put

def depthAverageT(T, z, d, lowerDepth, upperDepth):
  D = numpy.minimum(d, lowerDepth)
  H = numpy.zeros(T.shape[1:])
  HT = numpy.zeros(T.shape[1:])
  for k in range(T.shape[0]):
    zTop = numpy.minimum( z[k], -upperDepth )
    zBot = numpy.maximum( z[k+1], -D )
    dh = numpy.maximum( zTop -zBot, 0. )
    #dh = numpy.minimum( z[k]-z[k+1], D-H)
    H = H + dh
    HT = HT + dh*T[k]
  return HT/(H+1.e-20)

ci=m6plot.pmCI(0.25,4.5,.5)

tPlot = depthAverageT(Tmod, Zmod, depth, lDepth, uDepth)
tObsPlot = depthAverageT(Tobs, Zobs, depth, lDepth, uDepth)
m6plot.xyplot( tPlot - tObsPlot , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label,
      title=r'''%i-%im depth-average $\theta$ bias (w.r.t. WOA'05) [$\degree$C]'''%(uDepth,lDepth),
      clim=ci, colormap='dunnePM', centerlabels=True, extend='both',
      save=cmdLineArgs.outdir+'/T_%i-%im_zave_bias_WOA05.png'%(uDepth,lDepth))

m6plot.xycompare( tPlot, tObsPlot , x, y, area=area,
      suptitle=rootGroup.title+' '+cmdLineArgs.label,
      title1=r'%i-%im depth-average $\theta$ [$\degree$C]'%(uDepth,lDepth),
      title2=r'''WOA'05 $\theta$ [$\degree$C]''',
      clim=m6plot.linCI(-2,29,.5), colormap='dunneRainbow', extend='max',
      dlim=ci, dcolormap='dunnePM', dextend='both', centerdlabels=True,
      save=cmdLineArgs.outdir+'/T_%i-%im_zave_bias_WOA05.3_panel.png'%(uDepth,lDepth))
