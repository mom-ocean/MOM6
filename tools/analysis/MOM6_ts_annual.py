# Script to plot sub-surface ocean temperature drift.

# Analysis: using newer python 2.7.3
"""
module purge
module use -a /home/fms/local/modulefiles
module load gcc
module load netcdf/4.2
module load python/2.7.3 
"""

import os
import math
import numpy as np
from numpy import ma
from netCDF4 import Dataset, MFDataset, num2date, date2num
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Function to convert from page coordinates to non-dimensional coordinates
def page_to_ndc( panel, page ):

  if len(panel) == 4:
    ndc = [ 0.0, 0.0, 0.0, 0.0 ]
    ndc[0] = (panel[0]-page[0])/(page[2]-page[0])
    ndc[1] = (panel[1]-page[1])/(page[3]-page[1])
    ndc[2] = (panel[2]-panel[0])/(page[2]-page[0])
    ndc[3] = (panel[3]-panel[1])/(page[3]-page[1])
    return ndc

  elif len(panel) == 2:
    ndc = [ 0.0, 0.0 ]
    ndc[0] = (panel[0]-page[0])/(page[2]-page[0])
    ndc[1] = (panel[1]-page[1])/(page[3]-page[1])
    return ndc

# -----------------------------------------------------------------------------
# Function to discretize colormap with option to white out certain regions
def cmap_discretize(cmap, N, white=None):

    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)

    # White levels?
    if white != None:
      for i in range(N):
        if white[i] > 0.0:
          colors_rgba[i,:] = 1.0

    # Construct colormap distionary
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]

    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


# -----------------------------------------------------------------------------

# Radius of the earth (shared/constants/constants.F90)
radius = 6371.0e3

# Ocean heat capacity (ocean_core/ocean_parameters.F90)
cp_ocean = 3992.10322329649

# Read 'descriptor' and 'years' from external file
f = open("files.txt")
for line in f.readlines():
  exec(line.lstrip())
f.close()
model_label = "%s (%s)" % (descriptor,years)

# TMPDIR where input files are located
tmpdir = "./"

# Open input files
#fstatic = Dataset(tmpdir+'19000101.ocean_geometry.nc', 'r')
fstatic = Dataset(tmpdir+'ocean_annual.static.nc', 'r')
ftemp = MFDataset(tmpdir+'ocean_annual.*.temp.nc')
fsalt = MFDataset(tmpdir+'ocean_annual.*.salt.nc')

# Time info
time = ftemp.variables["time"]
ntimes = len(time[:])
date = num2date(time,time.units,time.calendar.lower())
year = [d.year for d in date]
time_days = date2num(date,'days since 01-01-0001',time.calendar.lower())

# Grid info
#area = fstatic.variables["Ah"][:]
area = fstatic.variables["area_t"][:]

z = ftemp.variables["zl"][:]
nz = len(z) 

# Input variables
temp = ftemp.variables["temp"]
salt = fsalt.variables["salt"]
# Create arrays to hold derived variables
ztemp = ma.array( np.zeros((ntimes,nz), 'float64'), mask=True )
zsalt = ma.array( np.zeros((ntimes,nz), 'float64'), mask=True )

# Loop over time
#for itime in range(ntimes):
for itime in range(1):

  # Compute vertical profile of zemperature
  tmp = temp[itime,:,:,:]
  contmp = salt[itime,:,:,:]
  for iz in range(nz):
    ztemp[itime,iz] = ma.average(tmp[iz,:,:], weights=area)
    zsalt[itime,iz] = ma.average(contmp[iz,:,:], weights=area)


# Transpose for compatibility with contour plots
ztemp = ztemp.transpose()
zsalt = zsalt.transpose()

# Close files
fstatic.close()
ftemp.close()
fsalt.close()

# -----------------------------------------------------------------------------
# Create plot

# Specify plots position in points: [left bottom right top]
page   = [   0.0,   0.0, 612.0, 792.0 ]    # corresponding to papertype='letter'
plot1a = [  89.0, 497.0, 480.0, 670.0 ]
plot1b = [  89.0, 324.0, 480.0, 497.0 ]
cbar   = [ 506.0, 324.0, 531.0, 670.0 ]
plot2  = [  89.0,  99.0, 480.0, 272.0 ]
plot   = [  89.0,  99.0, 480.0, 670.0 ]

#plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.dpi'] = 72.0
plt.rcParams['figure.figsize'] = [ (page[2]-page[0])/72.0, (page[3]-page[1])/72.0 ]

fig = plt.figure()
ax1a = plt.axes(page_to_ndc(plot,page))              

ax1a.set_ylim(5300,0)
ax1a.set_ylabel('Depth (m)')
ax1a.set_xlabel('Ocean Temp (C)',color='r')
ax1a.plot(ztemp,z,ls='-',color='r')
ax1b = ax1a.twiny()
ax1b.set_xlabel('Ocean Salinity (PSU)',color='b')
ax1b.plot(zsalt,z,ls='-',color='b')


# Figure title
xy = page_to_ndc([280.0,730.0],page)
fig.text(xy[0],xy[1],model_label,ha="center",size="x-large")

# Save figure
fig.savefig("ocean_temp.ps")

