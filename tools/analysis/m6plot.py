"""
A method for producing a standardized pseudo-color plot of 2D data
"""

import os
try: 
  if os.environ['DISPLAY'] != None: pass
except: 
  import matplotlib
  matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm
from matplotlib.ticker import MaxNLocator
import math
import numpy, numpy.matlib
import m6toolbox
import VerticalSplitScale

def xyplot(field, x=None, y=None, area=None,
  xlabel=None, xunits=None, ylabel=None, yunits=None,
  title='', suptitle='',
  clim=None, colormap=None, extend=None, centerlabels=False,
  nbins=None, landcolor=[.5,.5,.5],
  aspect=[16,9], resolution=576, axis=None,
  ignore=None, save=None, debug=False, show=False, interactive=False, logscale=False):
  """
  Renders plot of scalar field, field(x,y).

  Arguments:
  field        Scalar 2D array to be plotted.
  x            x coordinate (1D or 2D array). If x is the same size as field then x is treated as
               the cell center coordinates.
  y            y coordinate (1D or 2D array). If x is the same size as field then y is treated as
               the cell center coordinates.
  area         2D array of cell areas (used for statistics). Default None.
  xlabel       The label for the x axis. Default 'Longitude'.
  xunits       The units for the x axis. Default 'degrees E'.
  ylabel       The label for the y axis. Default 'Latitude'.
  yunits       The units for the y axis. Default 'degrees N'.
  title        The title to place at the top of the panel. Default ''.
  suptitle     The super-title to place at the top of the figure. Default ''.
  clim         A tuple of (min,max) color range OR a list of contour levels. Default None.
  colormap     The name of the colormap to use. Default None.
  extend       Can be one of 'both', 'neither', 'max', 'min'. Default None.
  centerlabels If True, will move the colorbar labels to the middle of the interval. Default False.
  nbins        The number of colors levels (used is clim is missing or only specifies the color range).
  landcolor    An rgb tuple to use for the color of land (no data). Default [.5,.5,.5].
  aspect       The aspect ratio of the figure, given as a tuple (W,H). Default [16,9].
  resolution   The vertical resolution of the figure given in pixels. Default 720.
  axis         The axis handle to plot to. Default None.
  ignore       A value to use as no-data (NaN). Default None.
  save         Name of file to save figure in. Default None.
  debug        If true, report stuff for debugging. Default False.
  show         If true, causes the figure to appear on screen. Used for testing. Default False.
  interactive  If true, adds interactive features such as zoom, close and cursor. Default False.
  logscale     If true, use logaritmic coloring scheme. Default False.
  """

  # Create coordinates if not provided
  xlabel, xunits, ylabel, yunits = createXYlabels(x, y, xlabel, xunits, ylabel, yunits)
  if debug: print 'x,y label/units=',xlabel,xunits,ylabel,yunits
  xCoord, yCoord = createXYcoords(field, x, y)

  # Diagnose statistics
  if ignore!=None: maskedField = numpy.ma.masked_array(field, mask=[field==ignore])
  else: maskedField = field.copy()
  sMin, sMax, sMean, sStd, sRMS = myStats(maskedField, area, debug=debug)
  xLims = boundaryStats(xCoord)
  yLims = boundaryStats(yCoord)

  # Choose colormap
  if nbins==None and (clim==None or len(clim)==2): nbins=35
  if colormap==None: colormap = chooseColorMap(sMin, sMax)
  cmap, norm, extend = chooseColorLevels(sMin, sMax, colormap, clim=clim, nbins=nbins, extend=extend, logscale=logscale)

  if axis==None:
    setFigureSize(aspect, resolution, debug=debug)
    #plt.gcf().subplots_adjust(left=.08, right=.99, wspace=0, bottom=.09, top=.9, hspace=0)
    axis = plt.gca()
  plt.pcolormesh(xCoord, yCoord, maskedField, cmap=cmap, norm=norm)
  if interactive: addStatusBar(xCoord, yCoord, maskedField)
  cb = plt.colorbar(fraction=.08, pad=0.02, extend=extend)
  if centerlabels and len(clim)>2: cb.set_ticks(  0.5*(clim[:-1]+clim[1:]) )
  axis.set_axis_bgcolor(landcolor)
  plt.xlim( xLims )
  plt.ylim( yLims )
  axis.annotate('max=%.5g\nmin=%.5g'%(sMax,sMin), xy=(0.0,1.01), xycoords='axes fraction', verticalalignment='bottom', fontsize=10)
  if area!=None:
    axis.annotate('mean=%.5g\nrms=%.5g'%(sMean,sRMS), xy=(1.0,1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right', fontsize=10)
    axis.annotate(' sd=%.5g\n'%(sStd), xy=(1.0,1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
  if len(xlabel+xunits)>0: plt.xlabel(label(xlabel, xunits))
  if len(ylabel+yunits)>0: plt.ylabel(label(ylabel, yunits))
  if len(title)>0: plt.title(title)
  if len(suptitle)>0: plt.suptitle(suptitle)

  if save!=None: plt.savefig(save)
  if interactive: addInteractiveCallbacks()
  if show: plt.show(block=False)


def xycompare(field1, field2, x=None, y=None, area=None,
  xlabel=None, xunits=None, ylabel=None, yunits=None,
  title1='', title2='', title3='A - B', addplabel=True, suptitle='',
  clim=None, colormap=None, extend=None, centerlabels=False,
  dlim=None, dcolormap=None, dextend=None, centerdlabels=False,
  nbins=None, landcolor=[.5,.5,.5],
  aspect=None, resolution=None, axis=None, npanels=3,
  ignore=None, save=None, debug=False, show=False, interactive=False):
  """
  Renders n-panel plot of two scalar fields, field1(x,y) and field2(x,y).

  Arguments:
  field1        Scalar 2D array to be plotted and compared to field2.
  field2        Scalar 2D array to be plotted and compared to field1.
  x             x coordinate (1D or 2D array). If x is the same size as field then x is treated as
                the cell center coordinates.
  y             y coordinate (1D or 2D array). If x is the same size as field then y is treated as
                the cell center coordinates.
  area          2D array of cell areas (used for statistics). Default None.
  xlabel        The label for the x axis. Default 'Longitude'.
  xunits        The units for the x axis. Default 'degrees E'.
  ylabel        The label for the y axis. Default 'Latitude'.
  yunits        The units for the y axis. Default 'degrees N'.
  title1        The title to place at the top of panel 1. Default ''.
  title2        The title to place at the top of panel 1. Default ''.
  title3        The title to place at the top of panel 1. Default 'A-B'.
  addplabel     Adds a 'A:' or 'B:' to the title1 and title2. Default True.
  suptitle      The super-title to place at the top of the figure. Default ''.
  clim          A tuple of (min,max) color range OR a list of contour levels for the field plots. Default None.
  colormap      The name of the colormap to use for the field plots. Default None.
  extend        Can be one of 'both', 'neither', 'max', 'min'. Default None.
  centerlabels  If True, will move the colorbar labels to the middle of the interval. Default False.
  dlim          A tuple of (min,max) color range OR a list of contour levels for the difference plot. Default None.
  dcolormap     The name of the colormap to use for the difference plot. Default None.
  dextend       For the difference colorbar. Can be one of 'both', 'neither', 'max', 'min'. Default None.
  centerdlabels If True, will move the difference colorbar labels to the middle of the interval. Default False.
  nbins         The number of colors levels (used is clim is missing or only specifies the color range).
  landcolor     An rgb tuple to use for the color of land (no data). Default [.5,.5,.5].
  aspect        The aspect ratio of the figure, given as a tuple (W,H). Default [16,9].
  resolution    The vertical resolution of the figure given in pixels. Default 1280.
  axis          The axis handle to plot to. Default None.
  npanels       Number of panels to display (1, 2 or 3). Default 3.
  ignore        A value to use as no-data (NaN). Default None.
  save          Name of file to save figure in. Default None.
  debug         If true, report stuff for debugging. Default False.
  show          If true, causes the figure to appear on screen. Used for testing. Default False.
  interactive   If true, adds interactive features such as zoom, close and cursor. Default False.
  """

  if (field1.shape)!=(field2.shape): raise Exception('field1 and field2 must be the same shape')

  # Create coordinates if not provided
  xlabel, xunits, ylabel, yunits = createXYlabels(x, y, xlabel, xunits, ylabel, yunits)
  if debug: print 'x,y label/units=',xlabel,xunits,ylabel,yunits
  xCoord, yCoord = createXYcoords(field1, x, y)

  # Diagnose statistics
  if ignore!=None: maskedField1 = numpy.ma.masked_array(field1, mask=[field1==ignore])
  else: maskedField1 = field1.copy()
  s1Min, s1Max, s1Mean, s1Std, s1RMS = myStats(maskedField1, area, debug=debug)
  if ignore!=None: maskedField2 = numpy.ma.masked_array(field2, mask=[field2==ignore])
  else: maskedField2 = field2.copy()
  s2Min, s2Max, s2Mean, s2Std, s2RMS = myStats(maskedField2, area, debug=debug)
  dMin, dMax, dMean, dStd, dRMS = myStats(maskedField1 - maskedField2, area, debug=debug)
  if s1Mean!=None: dRxy = corr(maskedField1 - s1Mean, maskedField2 - s2Mean, area)
  else: dRxy = None
  s12Min = min(s1Min, s2Min); s12Max = max(s1Max, s2Max)
  xLims = boundaryStats(xCoord); yLims = boundaryStats(yCoord)
  if debug:
    print 's1: min, max, mean =', s1Min, s1Max, s1Mean
    print 's2: min, max, mean =', s2Min, s2Max, s2Mean
    print 's12: min, max =', s12Min, s12Max

  # Choose colormap
  if nbins==None and (clim==None or len(clim)==2): cBins=35
  else: cBins=nbins
  if nbins==None and (dlim==None or len(dlim)==2): nbins=35
  if colormap==None: colormap = chooseColorMap(s12Min, s12Max)
  cmap, norm, extend = chooseColorLevels(s12Min, s12Max, colormap, clim=clim, nbins=cBins, extend=extend)

  def annotateStats(axis, sMin, sMax, sMean, sStd, sRMS):
    axis.annotate('max=%.5g\nmin=%.5g'%(sMax,sMin), xy=(0.0,1.025), xycoords='axes fraction', verticalalignment='bottom', fontsize=10)
    if sMean!=None:
      axis.annotate('mean=%.5g\nrms=%.5g'%(sMean,sRMS), xy=(1.0,1.025), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right', fontsize=10)
      axis.annotate(' sd=%.5g\n'%(sStd), xy=(1.0,1.025), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left', fontsize=10)

  if addplabel: preTitleA = 'A: '; preTitleB = 'B: '
  else: preTitleA = ''; preTitleB = ''

  if axis==None:
    setFigureSize(aspect, resolution, npanels=npanels, debug=debug)

  if npanels in [2,3]:
    axis = plt.subplot(npanels,1,1)
    plt.pcolormesh(xCoord, yCoord, maskedField1, cmap=cmap, norm=norm)
    if interactive: addStatusBar(xCoord, yCoord, maskedField1)
    cb1 = plt.colorbar(fraction=.08, pad=0.02, extend=extend)
    if centerlabels and len(clim)>2: cb1.set_ticks(  0.5*(clim[:-1]+clim[1:]) )
    axis.set_axis_bgcolor(landcolor)
    plt.xlim( xLims ); plt.ylim( yLims )
    annotateStats(axis, s1Min, s1Max, s1Mean, s1Std, s1RMS)
    axis.set_xticklabels([''])
    if len(ylabel+yunits)>0: plt.ylabel(label(ylabel, yunits))
    if len(title1)>0: plt.title(preTitleA+title1)

    axis = plt.subplot(npanels,1,2)
    plt.pcolormesh(xCoord, yCoord, maskedField2, cmap=cmap, norm=norm)
    if interactive: addStatusBar(xCoord, yCoord, maskedField2)
    cb2 = plt.colorbar(fraction=.08, pad=0.02, extend=extend)
    if centerlabels and len(clim)>2: cb2.set_ticks(  0.5*(clim[:-1]+clim[1:]) )
    axis.set_axis_bgcolor(landcolor)
    plt.xlim( xLims ); plt.ylim( yLims )
    annotateStats(axis, s2Min, s2Max, s2Mean, s2Std, s2RMS)
    if npanels>2: axis.set_xticklabels([''])
    if len(ylabel+yunits)>0: plt.ylabel(label(ylabel, yunits))
    if len(title2)>0: plt.title(preTitleB+title2)

  if npanels in [1,3]:
    axis = plt.subplot(npanels,1,npanels)
    if dcolormap==None: dcolormap = chooseColorMap(dMin, dMax)
    cmap, norm, extend = chooseColorLevels(dMin, dMax, dcolormap, clim=dlim, nbins=nbins, extend=dextend)
    plt.pcolormesh(xCoord, yCoord, maskedField1 - maskedField2, cmap=cmap, norm=norm)
    if interactive: addStatusBar(xCoord, yCoord, maskedField1 - maskedField2)
    cb3 = plt.colorbar(fraction=.08, pad=0.02, extend=dextend) # was extend!
    if centerdlabels and len(dlim)>2: cb3.set_ticks(  0.5*(dlim[:-1]+dlim[1:]) )
    axis.set_axis_bgcolor(landcolor)
    plt.xlim( xLims ); plt.ylim( yLims )
    annotateStats(axis, dMin, dMax, dMean, dStd, dRMS)
    if len(ylabel+yunits)>0: plt.ylabel(label(ylabel, yunits))
    if len(title3)>0: plt.title(title3)

  if dRxy!=None: axis.annotate(' r(A,B)=%.5g\n'%(dRxy), xy=(1.0,-0.20), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='center', fontsize=10)
  if len(xlabel+xunits)>0: plt.xlabel(label(xlabel, xunits))
  if len(suptitle)>0: plt.suptitle(suptitle)

  if save!=None: plt.savefig(save)
  if interactive: addInteractiveCallbacks()
  if show: plt.show(block=False)


def yzplot(field, y=None, z=None,
  ylabel=None, yunits=None, zlabel=None, zunits=None,
  splitscale=None,
  title='', suptitle='',
  clim=None, colormap=None, extend=None, centerlabels=False,
  nbins=None, landcolor=[.5,.5,.5],
  aspect=[16,9], resolution=576, axis=None,
  ignore=None, save=None, debug=False, show=False, interactive=False):
  """
  Renders section plot of scalar field, field(x,z).

  Arguments:
  field       Scalar 2D array to be plotted.
  y           y (or x) coordinate (1D array). If y is the same size as field then x is treated as
              the cell center coordinates.
  z           z coordinate (1D or 2D array). If z is the same size as field then y is treated as
              the cell center coordinates.
  ylabel      The label for the x axis. Default 'Latitude'.
  yunits      The units for the x axis. Default 'degrees N'.
  zlabel      The label for the z axis. Default 'Elevation'.
  zunits      The units for the z axis. Default 'm'.
  splitscale    A list of depths to define equal regions of projection in the vertical, e.g. [0.,-1000,-6500]
  title       The title to place at the top of the panel. Default ''.
  suptitle    The super-title to place at the top of the figure. Default ''.
  clim        A tuple of (min,max) color range OR a list of contour levels. Default None.
  colormap    The name of the colormap to use. Default None.
  extend      Can be one of 'both', 'neither', 'max', 'min'. Default None.
  centerlabels If True, will move the colorbar labels to the middle of the interval. Default False.
  nbins       The number of colors levels (used is clim is missing or only specifies the color range).
  landcolor   An rgb tuple to use for the color of land (no data). Default [.5,.5,.5].
  aspect      The aspect ratio of the figure, given as a tuple (W,H). Default [16,9].
  resolution  The vertical resolution of the figure given in pixels. Default 720.
  axis         The axis handle to plot to. Default None.
  ignore      A value to use as no-data (NaN). Default None.
  save        Name of file to save figure in. Default None.
  debug       If true, report stuff for debugging. Default False.
  show        If true, causes the figure to appear on screen. Used for testing. Default False.
  interactive If true, adds interactive features such as zoom, close and cursor. Default False.
  """

  # Create coordinates if not provided
  ylabel, yunits, zlabel, zunits = createYZlabels(y, z, ylabel, yunits, zlabel, zunits)
  if debug: print 'y,z label/units=',ylabel,yunits,zlabel,zunits
  if len(y)==z.shape[-1]: y = expand(y)
  elif len(y)==z.shape[-1]+1: y = y
  else: raise Exception('Length of y coordinate should be equal or 1 longer than horizontal length of z')
  if ignore!=None: maskedField = numpy.ma.masked_array(field, mask=[field==ignore])
  else: maskedField = field.copy()
  yCoord, zCoord, field2 = m6toolbox.section2quadmesh(y, z, maskedField)

  # Diagnose statistics
  sMin, sMax, sMean, sStd, sRMS = myStats(maskedField, yzWeight(y, z), debug=debug)
  yLims = numpy.amin(yCoord), numpy.amax(yCoord)
  zLims = boundaryStats(zCoord)

  # Choose colormap
  if nbins==None and (clim==None or len(clim)==2): nbins=35
  if colormap==None: colormap = chooseColorMap(sMin, sMax)
  cmap, norm, extend = chooseColorLevels(sMin, sMax, colormap, clim=clim, nbins=nbins, extend=extend)

  if axis==None:
    setFigureSize(aspect, resolution, debug=debug)
    #plt.gcf().subplots_adjust(left=.10, right=.99, wspace=0, bottom=.09, top=.9, hspace=0)
    axis = plt.gca()
  plt.pcolormesh(yCoord, zCoord, field2, cmap=cmap, norm=norm)
  if interactive: addStatusBar(yCoord, zCoord, field2)
  cb = plt.colorbar(fraction=.08, pad=0.02, extend=extend)
  if centerlabels and len(clim)>2: cb.set_ticks(  0.5*(clim[:-1]+clim[1:]) )
  axis.set_axis_bgcolor(landcolor)
  if splitscale!=None:
    for zzz in splitscale[1:-1]: plt.axhline(zzz,color='k',linestyle='--')
    axis.set_yscale('splitscale', zval=splitscale)
  plt.xlim( yLims ); plt.ylim( zLims )
  axis.annotate('max=%.5g\nmin=%.5g'%(sMax,sMin), xy=(0.0,1.01), xycoords='axes fraction', verticalalignment='bottom', fontsize=10)
  if sMean!=None:
    axis.annotate('mean=%.5g\nrms=%.5g'%(sMean,sRMS), xy=(1.0,1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right', fontsize=10)
    axis.annotate(' sd=%.5g\n'%(sStd), xy=(1.0,1.01), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left', fontsize=10)
  if len(ylabel+yunits)>0: plt.xlabel(label(ylabel, yunits))
  if len(zlabel+zunits)>0: plt.ylabel(label(zlabel, zunits))
  if len(title)>0: plt.title(title)
  if len(suptitle)>0: plt.suptitle(suptitle)

  if save!=None: plt.savefig(save)
  if interactive: addInteractiveCallbacks()
  if show: plt.show(block=False)


def yzcompare(field1, field2, y=None, z=None,
  ylabel=None, yunits=None, zlabel=None, zunits=None,
  splitscale=None,
  title1='', title2='', title3='A - B', addplabel=True, suptitle='',
  clim=None, colormap=None, extend=None, centerlabels=False,
  dlim=None, dcolormap=None, dextend=None, centerdlabels=False,
  nbins=None, landcolor=[.5,.5,.5],
  aspect=None, resolution=None, axis=None, npanels=3,
  ignore=None, save=None, debug=False, show=False, interactive=False):
  """
  Renders n-panel plot of two scalar fields, field1(x,y) and field2(x,y).

  Arguments:
  field1        Scalar 2D array to be plotted and compared to field2.
  field2        Scalar 2D array to be plotted and compared to field1.
  y             y coordinate (1D array). If y is the same size as field then y is treated as
                the cell center coordinates.
  z             z coordinate (1D or 2D array). If z is the same size as field then z is treated as
                the cell center coordinates.
  ylabel        The label for the y axis. Default 'Latitude'.
  yunits        The units for the y axis. Default 'degrees N'.
  zlabel        The label for the z axis. Default 'Elevation'.
  zunits        The units for the z axis. Default 'm'.
  splitscale    A list of depths to define equal regions of projection in the vertical, e.g. [0.,-1000,-6500]
  title1        The title to place at the top of panel 1. Default ''.
  title2        The title to place at the top of panel 1. Default ''.
  title3        The title to place at the top of panel 1. Default 'A-B'.
  addplabel     Adds a 'A:' or 'B:' to the title1 and title2. Default True.
  suptitle      The super-title to place at the top of the figure. Default ''.
  clim          A tuple of (min,max) color range OR a list of contour levels for the field plots. Default None.
  colormap      The name of the colormap to use for the field plots. Default None.
  extend        Can be one of 'both', 'neither', 'max', 'min'. Default None.
  centerlabels  If True, will move the colorbar labels to the middle of the interval. Default False.
  dlim          A tuple of (min,max) color range OR a list of contour levels for the difference plot. Default None.
  dcolormap     The name of the colormap to use for the differece plot. Default None.
  dextend       For the difference colorbar. Can be one of 'both', 'neither', 'max', 'min'. Default None.
  centerdlabels If True, will move the difference colorbar labels to the middle of the interval. Default False.
  nbins         The number of colors levels (used is clim is missing or only specifies the color range).
  landcolor     An rgb tuple to use for the color of land (no data). Default [.5,.5,.5].
  aspect        The aspect ratio of the figure, given as a tuple (W,H). Default [16,9].
  resolution    The vertical resolution of the figure given in pixels. Default 1280.
  axis          The axis handle to plot to. Default None.
  npanels       Number of panels to display (1, 2 or 3). Default 3.
  ignore        A value to use as no-data (NaN). Default None.
  save          Name of file to save figure in. Default None.
  debug         If true, report stuff for debugging. Default False.
  show          If true, causes the figure to appear on screen. Used for testing. Default False.
  interactive   If true, adds interactive features such as zoom, close and cursor. Default False.
  """

  if (field1.shape)!=(field2.shape): raise Exception('field1 and field2 must be the same shape')

  # Create coordinates if not provided
  ylabel, yunits, zlabel, zunits = createYZlabels(y, z, ylabel, yunits, zlabel, zunits)
  if debug: print 'y,z label/units=',ylabel,yunits,zlabel,zunits
  if len(y)==z.shape[-1]: y= expand(y)
  elif len(y)==z.shape[-1]+1: y= y
  else: raise Exception('Length of y coordinate should be equal or 1 longer than horizontal length of z')
  if ignore!=None: maskedField1 = numpy.ma.masked_array(field1, mask=[field1==ignore])
  else: maskedField1 = field1.copy()
  yCoord, zCoord, field1 = m6toolbox.section2quadmesh(y, z, maskedField1)

  # Diagnose statistics
  yzWeighting = yzWeight(y, z)
  s1Min, s1Max, s1Mean, s1Std, s1RMS = myStats(maskedField1, yzWeighting, debug=debug)
  if ignore!=None: maskedField2 = numpy.ma.masked_array(field2, mask=[field2==ignore])
  else: maskedField2 = field2.copy()
  yCoord, zCoord, field2 = m6toolbox.section2quadmesh(y, z, maskedField2)
  s2Min, s2Max, s2Mean, s2Std, s2RMS = myStats(maskedField2, yzWeighting, debug=debug)
  dMin, dMax, dMean, dStd, dRMS = myStats(maskedField1 - maskedField2, yzWeighting, debug=debug)
  dRxy = corr(maskedField1 - s1Mean, maskedField2 - s2Mean, yzWeighting)
  s12Min = min(s1Min, s2Min); s12Max = max(s1Max, s2Max)
  xLims = numpy.amin(yCoord), numpy.amax(yCoord); yLims = boundaryStats(zCoord)
  if debug:
    print 's1: min, max, mean =', s1Min, s1Max, s1Mean
    print 's2: min, max, mean =', s2Min, s2Max, s2Mean
    print 's12: min, max =', s12Min, s12Max

  # Choose colormap
  if nbins==None and (clim==None or len(clim)==2): cBins=35
  else: cBins=nbins
  if nbins==None and (dlim==None or len(dlim)==2): nbins=35
  if colormap==None: colormap = chooseColorMap(s12Min, s12Max)
  cmap, norm, extend = chooseColorLevels(s12Min, s12Max, colormap, clim=clim, nbins=cBins, extend=extend)

  def annotateStats(axis, sMin, sMax, sMean, sStd, sRMS):
    axis.annotate('max=%.5g\nmin=%.5g'%(sMax,sMin), xy=(0.0,1.025), xycoords='axes fraction', verticalalignment='bottom', fontsize=10)
    if sMean!=None:
      axis.annotate('mean=%.5g\nrms=%.5g'%(sMean,sRMS), xy=(1.0,1.025), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='right', fontsize=10)
      axis.annotate(' sd=%.5g\n'%(sStd), xy=(1.0,1.025), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='left', fontsize=10)

  if addplabel: preTitleA = 'A: '; preTitleB = 'B: '
  else: preTitleA = ''; preTitleB = ''

  if axis==None:
    setFigureSize(aspect, resolution, npanels=npanels, debug=debug)
    #plt.gcf().subplots_adjust(left=.13, right=.94, wspace=0, bottom=.05, top=.94, hspace=0.15)

  if npanels in [2, 3]:
    axis = plt.subplot(npanels,1,1)
    plt.pcolormesh(yCoord, zCoord, field1, cmap=cmap, norm=norm)
    if interactive: addStatusBar(yCoord, zCoord, field1)
    cb1 = plt.colorbar(fraction=.08, pad=0.02, extend=extend)
    if centerlabels and len(clim)>2: cb1.set_ticks(  0.5*(clim[:-1]+clim[1:]) )
    axis.set_axis_bgcolor(landcolor)
    plt.xlim( xLims ); plt.ylim( yLims )
    if splitscale!=None:
      for zzz in splitscale[1:-1]: plt.axhline(zzz,color='k',linestyle='--')
      axis.set_yscale('splitscale', zval=splitscale)
    annotateStats(axis, s1Min, s1Max, s1Mean, s1Std, s1RMS)
    axis.set_xticklabels([''])
    if len(zlabel+zunits)>0: plt.ylabel(label(zlabel, zunits))
    if len(title1)>0: plt.title(preTitleA+title1)

    axis = plt.subplot(npanels,1,2)
    plt.pcolormesh(yCoord, zCoord, field2, cmap=cmap, norm=norm)
    if interactive: addStatusBar(yCoord, zCoord, field2)
    cb2 = plt.colorbar(fraction=.08, pad=0.02, extend=extend)
    if centerlabels and len(clim)>2: cb2.set_ticks(  0.5*(clim[:-1]+clim[1:]) )
    axis.set_axis_bgcolor(landcolor)
    plt.xlim( xLims ); plt.ylim( yLims )
    if splitscale!=None:
      for zzz in splitscale[1:-1]: plt.axhline(zzz,color='k',linestyle='--')
      axis.set_yscale('splitscale', zval=splitscale)
    annotateStats(axis, s2Min, s2Max, s2Mean, s2Std, s2RMS)
    if npanels>2: axis.set_xticklabels([''])
    if len(zlabel+zunits)>0: plt.ylabel(label(zlabel, zunits))
    if len(title2)>0: plt.title(preTitleB+title2)

  if npanels in [1, 3]:
    axis = plt.subplot(npanels,1,npanels)
    if dcolormap==None: dcolormap = chooseColorMap(dMin, dMax)
    cmap, norm, extend = chooseColorLevels(dMin, dMax, dcolormap, clim=dlim, nbins=nbins, extend=dextend)
    plt.pcolormesh(yCoord, zCoord, field1 - field2, cmap=cmap, norm=norm)
    if interactive: addStatusBar(yCoord, zCoord, field1 - field2)
    cb3 = plt.colorbar(fraction=.08, pad=0.02, extend=dextend)
    if centerdlabels and len(dlim)>2: cb3.set_ticks(  0.5*(dlim[:-1]+dlim[1:]) )
    axis.set_axis_bgcolor(landcolor)
    plt.xlim( xLims ); plt.ylim( yLims )
    if splitscale!=None:
      for zzz in splitscale[1:-1]: plt.axhline(zzz,color='k',linestyle='--')
      axis.set_yscale('splitscale', zval=splitscale)
    annotateStats(axis, dMin, dMax, dMean, dStd, dRMS)
    if len(zlabel+zunits)>0: plt.ylabel(label(zlabel, zunits))

  axis.annotate(' r(A,B)=%.5g\n'%(dRxy), xy=(1.0,-0.20), xycoords='axes fraction', verticalalignment='bottom', horizontalalignment='center', fontsize=10)
  if len(ylabel+yunits)>0: plt.xlabel(label(ylabel, yunits))
  if len(title3)>0: plt.title(title3)
  if len(suptitle)>0: plt.suptitle(suptitle)

  if save!=None: plt.savefig(save)
  if interactive: addInteractiveCallbacks()
  if show: plt.show(block=False)


def chooseColorMap(sMin, sMax):
  """
  Based on the min/max extremes of the data, choose a colormap that fits the data.
  """
  if sMin<0 and sMax>0: return 'seismic'
  elif sMax>0 and sMin<0.1*sMax: return 'hot'
  elif sMin<0 and sMax>0.1*sMin: return 'hot_r'
  else: return 'spectral'


def chooseColorLevels(sMin, sMax, colorMapName, clim=None, nbins=None, steps=[1,2,2.5,5,10], extend=None, logscale=False):
  """
  If nbins is a positive integer, choose sensible color levels with nbins colors.
  If clim is a 2-element tuple, create color levels within the clim range
  or if clim is a vector, use clim as contour levels.
  If clim provides more than 2 color interfaces, nbins must be absent.
  If clim is absent, the sMin,sMax are used as the color range bounds.
  
  Returns cmap, norm and extend.
  """
  if nbins==None and clim==None: raise Exception('At least one of clim or nbins is required.')
  if clim!=None:
    if len(clim)<2: raise Exception('clim must be at least 2 values long.')
    if nbins==None and len(clim)==2: raise Exception('nbins must be provided when clims specifies a color range.')
    if nbins!=None and len(clim)>2: raise Exception('nbins cannot be provided when clims specifies color levels.')
  if clim==None: levels = MaxNLocator(nbins=nbins, steps=steps).tick_values(sMin, sMax)
  elif len(clim)==2: levels = MaxNLocator(nbins=nbins, steps=steps).tick_values(clim[0], clim[1])
  else: levels = clim

  nColors = len(levels)-1
  if extend==None:
    if sMin<levels[0] and sMax>levels[-1]: extend = 'both'#; eColors=[1,1]
    elif sMin<levels[0] and sMax<=levels[-1]: extend = 'min'#; eColors=[1,0]
    elif sMin>=levels[0] and sMax>levels[-1]: extend = 'max'#; eColors=[0,1]
    else: extend = 'neither'#; eColors=[0,0]
  eColors = [0,0]
  if extend in ['both', 'min']: eColors[0] = 1
  if extend in ['both', 'max']: eColors[1] = 1

  cmap = plt.get_cmap(colorMapName)#,lut=nColors+eColors[0]+eColors[1])
  #cmap0 = cmap(0.)
  #cmap1 = cmap(1.)
  #cmap = ListedColormap(cmap(range(eColors[0],nColors+1-eColors[1]+eColors[0])))#, N=nColors)
  #if eColors[0]>0: cmap.set_under(cmap0)
  #if eColors[1]>0: cmap.set_over(cmap1)
  if logscale: norm = LogNorm(vmin=levels[0], vmax=levels[-1])
  else: norm = BoundaryNorm(levels, ncolors=cmap.N)
  return cmap, norm, extend


def linCI(min, max, ci, *args):
  """
  Returns list of linearly spaced contour intervals from min to max with spacing ci.
  Unline numpy.arange this max is included IF max = min + ci*N for an integer N.
  """
  if len(args): return numpy.concatenate( (numpy.arange(min, max+ci, ci), linCI(*args)) )
  return numpy.arange(min, max+ci, ci)


def pmCI(min, max, ci, *args):
  """
  Returns list of linearly spaced contour intervals from -max to -min then min to max with spacing ci.
  Unline numpy.arange this max is included IF max = min + ci*N for an integer N.
  """
  ci = linCI(min, max, ci, *args)
  if ci[0]>0: return numpy.concatenate( (-ci[::-1],ci) )
  else: return numpy.concatenate( (-ci[::-1],ci[1:]) )


def myStats(s, area, s2=None, debug=False):
  """
  Calculates mean, standard deviation and root-mean-square of s.
  """
  sMin = numpy.ma.min(s); sMax = numpy.ma.max(s)
  if area==None: return sMin, sMax, None, None, None
  weight = area.copy()
  if debug: print 'myStats: sum(area) =',numpy.ma.sum(weight)
  if not numpy.ma.getmask(s).any()==numpy.ma.nomask: weight[s.mask] = 0.
  sumArea = numpy.ma.sum(weight)
  if debug: print 'myStats: sum(area) =',sumArea,'after masking'
  if debug: print 'myStats: sum(s) =',numpy.ma.sum(s)
  if debug: print 'myStats: sum(area*s) =',numpy.ma.sum(weight*s)
  mean = numpy.ma.sum(weight*s)/sumArea
  std = math.sqrt( numpy.ma.sum( weight*((s-mean)**2) )/sumArea )
  rms = math.sqrt( numpy.ma.sum( weight*(s**2) )/sumArea )
  if debug: print 'myStats: mean(s) =',mean
  if debug: print 'myStats: std(s) =',std
  if debug: print 'myStats: rms(s) =',rms
  return sMin, sMax, mean, std, rms


def corr(s1, s2, area):
  """
  Calculates the correlation coefficient between s1 and s2, assuming s1 and s2 have
  not mean. That is s1 = S - mean(S), etc.
  """
  weight = area.copy()
  if not numpy.ma.getmask(s1).any()==numpy.ma.nomask: weight[s1.mask] = 0.
  sumArea = numpy.ma.sum(weight)
  v1 = numpy.ma.sum( weight*(s1**2) )/sumArea
  v2 = numpy.ma.sum( weight*(s2**2) )/sumArea
  if v1==0 or v2==0: return numpy.NaN
  rxy = numpy.ma.sum( weight*(s1*s2) )/sumArea / math.sqrt( v1*v2 )
  return rxy


def createXYcoords(s, x, y):
  """
  Checks that x and y are appropriate 2D corner coordinates
  and tries to make some if they are not.
  """
  nj, ni = s.shape
  if x==None: xCoord = numpy.arange(0., ni+1)
  else: xCoord = numpy.ma.filled(x, 0.)
  if y==None: yCoord = numpy.arange(0., nj+1)
  else: yCoord = numpy.ma.filled(y, 0.)

  # Turn coordinates into 2D arrays if 1D arrays were provided
  if len(xCoord.shape)==1:
    nxy = yCoord.shape
    xCoord = numpy.matlib.repmat(xCoord, nxy[0], 1)
  nxy = xCoord.shape
  if len(yCoord.shape)==1: yCoord = numpy.matlib.repmat(yCoord.T, nxy[-1], 1).T
  if xCoord.shape!=yCoord.shape: raise Exception('The shape of coordinates are mismatched!')

  # Create corner coordinates from center coordinates is center coordinates were provided
  if xCoord.shape!=yCoord.shape: raise Exception('The shape of coordinates are mismatched!')
  if s.shape==xCoord.shape:
    xCoord = expandJ( expandI( xCoord ) )
    yCoord = expandJ( expandI( yCoord ) )
  return xCoord, yCoord


def expandI(a):
  """
  Expands an array by one column, averaging the data to the middle columns and
  extrapolating for the first and last columns. Needed for shifting coordinates
  from centers to corners.
  """
  nj, ni = a.shape
  b = numpy.zeros((nj, ni+1))
  b[:,1:-1] = 0.5*( a[:,:-1] + a[:,1:] )
  b[:,0] = a[:,0] + 0.5*( a[:,0] - a[:,1] )
  b[:,-1] = a[:,-1] + 0.5*( a[:,-1] - a[:,-2] )
  return b


def expandJ(a):
  """
  Expands an array by one row, averaging the data to the middle columns and
  extrapolating for the first and last rows. Needed for shifting coordinates
  from centers to corners.
  """
  nj, ni = a.shape
  b = numpy.zeros((nj+1, ni))
  b[1:-1,:] = 0.5*( a[:-1,:] + a[1:,:] )
  b[0,:] = a[0,:] + 0.5*( a[0,:] - a[1,:] )
  b[-1,:] = a[-1,:] + 0.5*( a[-1,:] - a[-2,:] )
  return b


def expand(a):
  """
  Expands a vector by one element, averaging the data to the middle columns and
  extrapolating for the first and last rows. Needed for shifting coordinates
  from centers to corners.
  """
  b = numpy.zeros((len(a)+1))
  b[1:-1] = 0.5*( a[:-1] + a[1:] )
  b[0] = a[0] + 0.5*( a[0] - a[1] )
  b[-1] = a[-1] + 0.5*( a[-1] - a[-2] )
  return b


def boundaryStats(a):
  """
  Returns the minimum and maximum values of a only on the boundaries of the array.
  """
  amin = numpy.amin(a[0,:])
  amin = min(amin, numpy.amin(a[1:,-1]))
  amin = min(amin, numpy.amin(a[-1,:-1]))
  amin = min(amin, numpy.amin(a[1:-1,0]))
  amax = numpy.amax(a[0,:])
  amax = max(amax, numpy.amax(a[1:,-1]))
  amax = max(amax, numpy.amax(a[-1,:-1]))
  amax = max(amax, numpy.amax(a[1:-1,0]))
  return amin, amax


def setFigureSize(aspect=None, verticalresolution=None, horiztonalresolution=None,
  npanels=1, debug=False):
  """
  Set the figure size based on vertical resolution and aspect ratio (tuple of W,H).
  """
  if (not horiztonalresolution==None) and (not verticalresolution==None):
    if aspect==None: aspect=[horiztonalresolution, verticalresolution]
    else: raise Exception('Aspect-ratio and both h-/v- resolutions can not be specified together')
  if aspect==None: aspect = {1:[16,9], 2:[1,1], 3:[7,10]}[npanels]
  if (not horiztonalresolution==None) and (verticalresolution==None):
    verticalresolution = int(1.*aspect[1]/aspect[0] * horiztonalresolution)
  if verticalresolution==None: verticalresolution = {1:576, 2:720, 3:1200}[npanels]
  width = int(1.*aspect[0]/aspect[1] * verticalresolution) # First guess
  if debug: print 'setFigureSize: first guess width =',width
  width = width + ( width % 2 ) # Make even
  if debug: print 'setFigureSize: corrected width =',width
  if debug: print 'setFigureSize: height =',verticalresolution
  plt.figure(figsize=(width/100., verticalresolution/100.)) # 100 dpi always?
  if npanels==1: plt.gcf().subplots_adjust(left=.08, right=.99, wspace=0, bottom=.09, top=.9, hspace=0)
  elif npanels==2: plt.gcf().subplots_adjust(left=.11, right=.94, wspace=0, bottom=.09, top=.9, hspace=0.15)
  elif npanels==3: plt.gcf().subplots_adjust(left=.11, right=.94, wspace=0, bottom=.05, top=.93, hspace=0.15)
  elif npanels==0: pass
  else: raise Exception('npanels out of range')


def label(label, units):
  """
  Combines a label string and units string together in the form 'label [units]'
  unless one of the other is empty.
  """
  string = unicode(label)
  if len(units)>0: string = string + ' [' + unicode(units) + ']'
  return string


def createXYlabels(x, y, xlabel, xunits, ylabel, yunits):
  """
  Checks that x and y labels are appropriate and tries to make some if they are not.
  """
  if x==None:
    if xlabel==None: xlabel='i'
    if xunits==None: xunits=''
  else:
    if xlabel==None: xlabel='Longitude'
    #if xunits==None: xunits=u'\u00B0E'
    if xunits==None: xunits=r'$\degree$E'
  if y==None:
    if ylabel==None: ylabel='j'
    if yunits==None: yunits=''
  else:
    if ylabel==None: ylabel='Latitude'
    #if yunits==None: yunits=u'\u00B0N'
    if yunits==None: yunits=r'$\degree$N'
  return xlabel, xunits, ylabel, yunits


def addInteractiveCallbacks():
  """
  Adds interactive features to a plot on screen.
  Key 'q' to close window.
  Zoom button to center.
  Zoom wheel to zoom in and out.
  """
  def keyPress(event):
    if event.key=='Q': exit(0) # Exit python
    elif event.key=='q': plt.close() # Close just the active figure
  class hiddenStore:
    def __init__(self,axis):
      self.axis = axis
      self.xMin, self.xMax = axis.get_xlim()
      self.yMin, self.yMax = axis.get_ylim()
  save = hiddenStore(plt.gca())
  def zoom(event): # Scroll wheel up/down
    if event.button == 'up': scaleFactor = 1/1.5 # deal with zoom in
    elif event.button == 'down': scaleFactor = 1.5 # deal with zoom out
    elif event.button == 2: scaleFactor = 1.0
    else: return
    axis = event.inaxes
    axmin,axmax=axis.get_xlim(); aymin,aymax=axis.get_ylim();
    (axmin,axmax),(aymin,aymax) = newLims(
        (axmin,axmax), (aymin,aymax), (event.xdata, event.ydata),
        (save.xMin,save.xMax), (save.yMin,save.yMax), scaleFactor)
    if axmin==None: return
    for axis in plt.gcf().get_axes():
      if axis.get_navigate():
        axis.set_xlim(axmin, axmax); axis.set_ylim(aymin, aymax)
    plt.draw() # force re-draw
  def zoom2(event): zoom(event)
  plt.gcf().canvas.mpl_connect('key_press_event', keyPress)
  plt.gcf().canvas.mpl_connect('scroll_event', zoom)
  plt.gcf().canvas.mpl_connect('button_press_event', zoom2)


def addStatusBar(xCoord, yCoord, zData):
  """
  Reformats status bar message
  """
  class hiddenStore:
    def __init__(self,axis):
      self.axis = axis
      self.xMin, self.xMax = axis.get_xlim()
      self.yMin, self.yMax = axis.get_ylim()
  save = hiddenStore(plt.gca())
  def statusMessage(x,y):
    # THIS NEEDS TESTING FOR ACCURACY, ESPECIALLY IN YZ PLOTS -AJA
    if len(xCoord.shape)==1 and len(yCoord.shape)==1:
      # -2 needed because of coords are for vertices and need to be averaged to centers
      i = min(range(len(xCoord)-2), key=lambda l: abs((xCoord[l]+xCoord[l+1])/2.-x))
      j = min(range(len(yCoord)-2), key=lambda l: abs((yCoord[l]+yCoord[l+1])/2.-y))
    elif len(xCoord.shape)==1 and len(yCoord.shape)==2:
      i = min(range(len(xCoord)-2), key=lambda l: abs((xCoord[l]+xCoord[l+1])/2.-x))
      j = min(range(len(yCoord[:,i])-1), key=lambda l: abs((yCoord[l,i]+yCoord[l+1,i])/2.-y))
    elif len(xCoord.shape)==2 and len(yCoord.shape)==2:
      idx = numpy.abs( numpy.fabs( xCoord[0:-1,0:-1]+xCoord[1:,1:]+xCoord[0:-1,1:]+xCoord[1:,0:-1]-4*x)
          +numpy.fabs( yCoord[0:-1,0:-1]+yCoord[1:,1:]+yCoord[0:-1,1:]+yCoord[1:,0:-1]-4*y) ).argmin()
      j,i = numpy.unravel_index(idx,zData.shape)
    else: raise Exception('Combindation of coordinates shapes is VERY UNUSUAL!')
    if not i==None:
      val = zData[j,i]
      if val is numpy.ma.masked: return 'x,y=%.3f,%.3f  f(%i,%i)=NaN'%(x,y,i+1,j+1)
      else: return 'x,y=%.3f,%.3f  f(%i,%i)=%g'%(x,y,i+1,j+1,val)
    else: return 'x,y=%.3f,%.3f'%(x,y)
  plt.gca().format_coord = statusMessage


def newLims(cur_xlim, cur_ylim, cursor, xlim, ylim, scale_factor):
  cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
  cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
  xdata = cursor[0]; ydata = cursor[1]
  new_xrange = cur_xrange*scale_factor; new_yrange = cur_yrange*scale_factor
  xdata = min( max( xdata, xlim[0]+new_xrange ), xlim[1]-new_xrange )
  xL = max( xlim[0], xdata - new_xrange ); xR = min( xlim[1], xdata + new_xrange )
  if ylim[1]>ylim[0]:
    ydata = min( max( ydata, ylim[0]+new_yrange ), ylim[1]-new_yrange )
    yL = max( ylim[0], ydata - new_yrange ); yR = min( ylim[1], ydata + new_yrange )
  else:
    ydata = min( max( ydata, ylim[1]-new_yrange ), ylim[0]+new_yrange )
    yR = max( ylim[1], ydata + new_yrange ); yL = min( ylim[0], ydata - new_yrange )
  if xL==cur_xlim[0] and xR==cur_xlim[1] and \
     yL==cur_ylim[0] and yR==cur_ylim[1]: return (None, None), (None, None)
  return (xL, xR), (yL, yR)


def createYZlabels(y, z, ylabel, yunits, zlabel, zunits):
  """
  Checks that y and z labels are appropriate and tries to make some if they are not.
  """
  if y==None:
    if ylabel==None: ylabel='j'
    if yunits==None: yunits=''
  else:
    if ylabel==None: ylabel='Latitude'
    #if yunits==None: yunits=u'\u00B0N'
    if yunits==None: yunits=r'$\degree$N'
  if z==None:
    if zlabel==None: zlabel='k'
    if zunits==None: zunits=''
  else:
    if zlabel==None: zlabel='Elevation'
    if zunits==None: zunits='m'
  return ylabel, yunits, zlabel, zunits


def yzWeight(y, z):
  """
  Calculates the weights to use when calculating the statistics of a y-z section.

  y(nj+1) is a 1D vector of column edge positions and z(nk+1,nj) is the interface
  elevations of each column. Returns weight(nk,nj).
  """
  dz = z[:-1,:] - z[1:,:]
  return numpy.matlib.repmat(y[1:] - y[:-1], dz.shape[0], 1) * dz


def dunne_rainbow(N=256):
  """
  Spectral/rainbow colormap from John Dunne.
  """
  cdict = {'red': [(0.00, 0.95, 0.95),
                   (0.09, 0.85, 0.85),
                   (0.18, 0.60, 0.60),
                   (0.32, 0.30, 0.30),
                   (0.45, 0.00, 0.00),
                   (0.60, 1.00, 1.00),
                   (0.85, 1.00, 1.00),
                   (1.00, 0.40, 0.00)],

         'green': [(0.00, 0.75, 0.75),
                   (0.09, 0.85, 0.85),
                   (0.18, 0.60, 0.60),
                   (0.32, 0.20, 0.20),
                   (0.45, 0.60, 0.60),
                   (0.60, 1.00, 1.00),
                   (0.73, 0.70, 0.70),
                   (0.85, 0.00, 0.00),
                   (1.00, 0.00, 0.00)],

         'blue':  [(0.00, 1.00, 1.00),
                   (0.32, 1.00, 1.00),
                   (0.45, 0.30, 0.30),
                   (0.60, 0.00, 0.00),
                   (1.00, 0.00, 0.00)]}
  import matplotlib
  cmap = matplotlib.colors.LinearSegmentedColormap('dunneRainbow', cdict, N=N)
  #cmap.set_under([1,.65,.85]); cmap.set_over([.25,0.,0.])
  cmap.set_under([.95*.9,.75*.9,.9]); cmap.set_over([.3,0.,0.])
  #cmap.set_bad('w')
  matplotlib.cm.register_cmap(cmap=cmap)
  return cmap


def dunne_pm(N=256):
  """
  Plus-minus  colormap from John Dunne.
  """
  cdict = {'red':   [(0.00, 0.3, 0.3),
                     (0.05, 0.5, 0.5),
                     (0.20, 0.0, 0.0),
                     (0.30, 0.4, 0.4),
                     (0.40, 0.8, 0.8),
                     (0.50, 1.0, 1.0),
                     (0.95, 0.6, 0.6),
                     (1.00, 0.4, 0.4)],
  
           'green': [(0.00, 0.0, 0.0),
                     (0.30, 0.5, 0.5),
                     (0.40, 1.0, 1.0),
                     (0.70, 1.0, 1.0),
                     (1.00, 0.0, 0.0)],
  
           'blue':  [(0.00, 0.3, 0.3),
                     (0.05, 0.5, 0.5),
                     (0.20, 1.0, 1.0),
                     (0.50, 1.0, 1.0),
                     (0.60, 0.7, 0.7),
                     (0.70, 0.0, 0.0),
                     (1.00, 0.0, 0.0)]}
  import matplotlib
  cmap = matplotlib.colors.LinearSegmentedColormap('dunnePM', cdict, N=N)
  cmap.set_under([.1,.0,.1]); cmap.set_over([.2,0.,0.])
  #cmap.set_bad('w')
  matplotlib.cm.register_cmap(cmap=cmap)
  return cmap

# Load new named colormaps
c = dunne_rainbow()
c = dunne_pm()

# Test
if __name__ == '__main__':
  import nccf
  file = 'baseline/19000101.ocean_static.nc'
  D,(y,x),_ = nccf.readVar(file,'depth_ocean')
  y,_,_ = nccf.readVar(file,'geolat')
  x,_,_ = nccf.readVar(file,'geolon')
  area,_,_ = nccf.readVar(file,'area_t')
  xyplot(D, x, y, title='Depth', ignore=0, suptitle='Testing', area=area, clim=[0, 5500], nbins=12, debug=True, interactive=True, show=True)#, save='fig_test.png')
  xycompare(D, .9*D, x, y, title1='Depth', ignore=0, suptitle='Testing', area=area, nbins=12)#, save='fig_test2.png')
  annual = 'baseline/19000101.ocean_annual.nc'
  monthly = 'baseline/19000101.ocean_month.nc'
  e,(t,z,y,x),_ = nccf.readVar(annual,'e',0,None,None,1100)
  temp,(t,z,y,x),_ = nccf.readVar(monthly,'temp',0,None,None,1100)
  temp2,(t,z,y,x),_ = nccf.readVar(monthly,'temp',11,None,None,1100)
  yzplot(temp, y, e)
  yzcompare(temp, temp2, y, e, interactive=True)
  yzcompare(temp, temp2, y, e, npanels=2)
  yzcompare(temp, temp2, y, e, npanels=1)
  plt.show()
