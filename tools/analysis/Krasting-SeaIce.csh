#!/bin/csh 
#------------------------------------
#PBS -N Krasting-SeaIce
#PBS -l size=1
#PBS -l walltime=60:00:00
#PBS -r y
#PBS -j oe
#PBS -o
#PBS -q batch
#----------------------------------
# Script:  Krasting-SeaIce.csh
# Authors: John Krasting
#          (Based on Figures and Data from the ESMValTool Project)
#
# Source data: pp/ice/ts/monthly
#
# Output:  Creates figures in $out_dir/ice_${yr1}_${yr2}
#
# Sample frepp usage (http://www.gfdl.noaa.gov/fms/fre/#analysis):
# <component type="ice">
#    <timeSeries ... >
#       <analysis script="script_name [options]"/>
#    </timeSeries>
# </component>
#
# Overview:
# This analysis script generates polar stereographic maps of both
# Arctic (March) and Antarctic (September) sea ice concentration
# relative to the NSIDC 1981-2000 average.
#
# The script also generates line plots of the annual cycle of both
# sea ice area and extent.  Sea ice extent (also denoted as a red
# line on the maps) is based on a minimum sea ice concentration of
# 15% which is the same threshold used by NSIDC.
#
# The plots in this script are inspired by the ESMValTool that is 
# being developed as part of the EMBRACE project.

#clean up the temp directory
wipetmp

#fields set by frepp
set argu
set descriptor
set in_data_dir
set out_dir
set WORKDIR    
set mode       
set yr1
set yr2
set databegyr
set dataendyr
set datachunk
set staticfile
set script_path
set fremodule
set gridspecfile

set oceanstaticfile = `echo ${staticfile} | sed -e 's/ice/ocean_monthly/g'`

set model = "GFDL"

# make sure valid platform and required modules are loaded
if (`gfdl_platform` == "hpcs-csc") then
   source $MODULESHOME/init/csh
   module use -a /home/fms/local/modulefiles
   module purge
   module load $fremodule
   module load fre-analysis/test
   module unload python

   # This script relies on the ESMF regridder, available in 
   # uv-cdat version 1.5 and higher.  This package is currently
   # not available as a module on GFDL PAN. As a substitute,
   # source the paths of a local uv-cdat build, which relies on Qt.
   
   module load qt
   source /nbhome/jpk/uvcdat-1.5/install/bin/setup_runtime.csh
   which python

   setenv UVCDAT_ANONYMOUS_LOG no

else
   echo "ERROR: invalid platform"
   exit 1
endif

# check again?
if (! $?FRE_ANALYSIS_HOME) then
   echo "ERROR: environment variable FRE_ANALYSIS_HOME not set."
   exit 1
endif

# Create a working directory
set workDir = `mktemp -d -p $TMPDIR`
cd ${workDir}
mkdir -p ${workDir}/data

# Generate a list of files to dmget, then dmget them
cd ${in_data_dir}
set yrStart = `echo $yr1 | awk '{x=$0+0;print x}'`
set yrEnd = `echo $yr2 | awk '{x=$0+0;print x}'`
set fileList = ""
foreach var (CN)
 foreach f (`dir -1 *.${var}.nc`)
  set yr1_1 = `echo $f | cut -f 2 -d '.' | cut -c 1-4 | awk '{x=$0+0;print x}'`
  set yr2_2 = `echo $f | cut -f 2 -d '.' | cut -c 8-11 | awk '{x=$0+0;print x}'`
  if (($yr1_1 >= $yrStart) && ($yr2_2 <= $yrEnd)) then
    set fileList = "$fileList $f"
  else if (($yr1_1 >= $yrStart) && ($yr2_2 <= ($yrEnd + $datachunk))) then
    set fileList = "$fileList $f"
  endif
 end
end
echo "dmgetting $fileList"
dmget $fileList
echo "copying files to work directory"
foreach f ($fileList)
 gcp -v $f ${workDir}/data/
end

# Bring in observed data
dmget /archive/fms/fre-analysis/test/John.Krasting/ESMValTool-NSIDC-SeaIce-Obs-20140521.tar

cd ${workDir}
echo "Copying input data sets to local disk"
gcp -v /archive/fms/fre-analysis/test/John.Krasting/ESMValTool-NSIDC-SeaIce-Obs-20140521.tar ${workDir}
cd ${workDir}/data
tar -xvf ${workDir}/ESMValTool-NSIDC-SeaIce-Obs-20140521.tar
cd ${workDir}

# Get and unpack gridspecfile
cd ${workDir}
mkdir -p gridspec
if ( ${gridspecfile} =~ "*.nc" ) then  
  set gridspecdir = `echo ${gridspecfile} | rev | cut -f 2-100 -d '/' | rev`
  gcp -v ${gridspecdir}/ocean_hgrid.nc ${workDir}/gridspec
else
  dmget ${gridspecfile}
  gcp -v ${gridspecfile} ${workDir}
  cd gridspec
  tar -xvf ${workDir}/`basename ${gridspecfile}`
endif
cd ${workDir}

# Make an input XML of the Sea Ice Data
cd ${workDir}
cdscan -x input.xml data/ice*.nc

#generate a namelist that is required as input for the tool
cat > ${workDir}/Krasting-SeaIce.py <<EOF
import matplotlib as mpl
mpl.use('Agg')
import cdms2, numpy, scipy, MV2, netCDF4, cdutil
import cdms2.coord, cdms2.hgrid
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

expName = '${descriptor}'

#cdms2.setAutoBounds('on')

#-- Read coordinate bounds
lon = netCDF4.Dataset('gridspec/ocean_hgrid.nc').variables['x'][::2,::2]
lat = netCDF4.Dataset('gridspec/ocean_hgrid.nc').variables['y'][::2,::2]
geolon = netCDF4.Dataset('gridspec/ocean_hgrid.nc').variables['x'][1::2,1::2]
geolat = netCDF4.Dataset('gridspec/ocean_hgrid.nc').variables['y'][1::2,1::2]
area_t = netCDF4.Dataset('gridspec/ocean_hgrid.nc').variables['area'][:].reshape((geolon.shape[0],2,geolon.shape[1],2)).sum(axis=1).sum(axis=-1)

#-- Read geolon/geolat
fs = cdms2.open('${oceanstaticfile}')
yh = fs('area_t').getAxis(0)
xh = fs('area_t').getAxis(1)

#-- Read data and compute climatology
f = cdms2.open('input.xml')
CN = f('CN',time=('${yr1}-01-01','${yr2}-12-31'))
year1Label = str(CN.getTime().asComponentTime()[0]).split('-')[0]
year2Label = str(CN.getTime().asComponentTime()[-1]).split('-')[0]
tax = CN.getAxis(0)
sic = MV2.array(CN.sum(axis=1)*100.)
sic.setAxisList((tax,yh,xh))
print "doing average ..."
cdutil.setTimeBoundsMonthly(sic)
sic = cdutil.ANNUALCYCLE.climatology(sic)
print "done"
timeMod = sic.getAxis(0)
sic_area = sic*area_t
sic_ext  = MV2.where(MV2.greater(sic,15.0),100.0,0.0)*area_t

#-- Read NSIDC Climatology (from K. Gottschaldt)
fo = cdms2.open('data/OBS_NSIDC_sat_NH_T2Ms_sic.nc')
lonobs = fo('lon')
latobs = fo('lat')
areacello = fo('areacello')
sic_nsidc = fo('sic',time=('1981-1-1','2000-12-31'))*100.0
cdutil.setTimeBoundsMonthly(sic_nsidc)
sic_nsidc = cdutil.ANNUALCYCLE.climatology(sic_nsidc)
timeObs = sic_nsidc.getAxis(0)
sic_area_nsidc = sic_nsidc*areacello
sic_ext_nsidc  = MV2.where(MV2.greater(sic_nsidc,15.0),100.0,0.0)*areacello

#-- Specify transient grid for obs data
iaxisObs = cdms2.coord.TransientVirtualAxis("i", len(sic_nsidc.getAxis(2)))
jaxisObs = cdms2.coord.TransientVirtualAxis("j", len(sic_nsidc.getAxis(1)))
lonAxisObs = cdms2.coord.TransientAxis2D(lonobs,axes=(jaxisObs, iaxisObs),attributes={'units': 'degrees_east'},id='lon')
latAxisObs = cdms2.coord.TransientAxis2D(latobs,axes=(jaxisObs, iaxisObs),attributes={'units': 'degrees_north'},id='lat')
obsGrid = cdms2.hgrid.TransientCurveGrid(latAxisObs, lonAxisObs, id='lats_lons')
sic_nsidc_grd = MV2.array(sic_nsidc, axes=[timeObs, jaxisObs,iaxisObs], grid=obsGrid, missing=0.0)

#-- Specify transient grid for model data
iaxisMod = cdms2.coord.TransientVirtualAxis("i", len(sic.getAxis(2)))
jaxisMod = cdms2.coord.TransientVirtualAxis("j", len(sic.getAxis(1)))
lonAxisMod = cdms2.coord.TransientAxis2D(geolon,axes=(jaxisMod, iaxisMod),attributes={'units': 'degrees_east'},id='lon')
latAxisMod = cdms2.coord.TransientAxis2D(geolat,axes=(jaxisMod, iaxisMod),attributes={'units': 'degrees_north'},id='lat')
ModGrid = cdms2.hgrid.TransientCurveGrid(latAxisMod, lonAxisMod, id='lats_lons')
sic_grd = MV2.array(sic, axes=[timeMod, jaxisMod,iaxisMod], grid=ModGrid)

#-- Regrid obs to model grid
sic_nsidc_regrid = sic_nsidc_grd.regrid(sic_grd.getGrid(),regridTool='esmf',regridMethod='conserve',fixSrcBounds=True,fixDstBounds=True)
sic_nsidc_regrid.data[:,:,:] = numpy.where(numpy.greater(sic_nsidc_regrid.data,1.e19),0.0,sic_nsidc_regrid.data)
sic_nsidc_regrid.mask = sic.mask

#-- Make a plot
fig = plt.figure(figsize=(11,8.5))
plt.subplot(2,3,1)
ice = sic[2,:,:]
m = Basemap(projection='npstere',boundinglat=55,lon_0=270,resolution='l')
#m = Basemap(projection='kav7',lon_0=-120,resolution=None)
m.drawmapboundary(fill_color='0.3')
#im0 = m.bluemarble(scale=0.5)
im1 = m.pcolormesh(numpy.minimum(lon,60.),lat,ice,shading='flat',cmap='Blues_r',latlon=True)
plt.clim(0.,100.)
im2 = m.contour(geolon,geolat,ice,[15],latlon=True,colors='red',linewidth=5)
plt.clim(0.0,100.)
m.drawparallels(numpy.arange(-90.,99.,30.))
m.drawmeridians(numpy.arange(-180.,180.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
plt.title('Model - Years ' + year1Label + ' to ' + year2Label)

plt.subplot(2,3,2)
ice = sic_nsidc[2,:,:]
m = Basemap(projection='npstere',boundinglat=55,lon_0=270,resolution='l')
#m = Basemap(projection='kav7',lon_0=-120,resolution=None)
m.drawmapboundary(fill_color='0.3')
#im0 = m.bluemarble(scale=0.5)
im1 = m.pcolormesh(lonobs,latobs,ice,shading='flat',cmap='Blues_r',latlon=True)
plt.clim(0.,100.)
im2 = m.contour(lonobs,latobs,ice,[15],latlon=True,colors='red',linewidth=5)
plt.clim(0.0,100.)
m.drawparallels(numpy.arange(-90.,99.,30.))
m.drawmeridians(numpy.arange(-180.,180.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
plt.title('NSIDC - Years 1981 to 2000')

plt.subplot(2,3,3)
ice = sic[2,:,:] - sic_nsidc_regrid[2,:,:]
m = Basemap(projection='npstere',boundinglat=55,lon_0=270,resolution='l')
#m = Basemap(projection='kav7',lon_0=-120,resolution=None)
m.drawmapboundary(fill_color='0.3')
#im0 = m.bluemarble(scale=0.5)
im1 = m.pcolormesh(numpy.minimum(lon,60.),lat,ice,shading='flat',cmap='RdBu',latlon=True)
plt.clim(-100.,100.)
m.drawparallels(numpy.arange(-90.,99.,30.))
m.drawmeridians(numpy.arange(-180.,180.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%", ticks=[-100.,-60., -30., 0., 30., 60., 100.])
plt.title('Difference')

ax4 = plt.subplot(2,3,4)
mons = numpy.linspace(1,13,13)
model = numpy.array(sic_area(latitude=(0,90)).sum(axis=1).sum(axis=1))
model = numpy.concatenate((model[7:12], model[0:8]))*1.e-14
obs = numpy.array(sic_area_nsidc(latitude=(0,90)).sum(axis=1).sum(axis=1))
obs = numpy.concatenate((obs[7:12], obs[0:8]))*1.e-14
max_area_model = str(numpy.around(model.max(),decimals=5))
min_area_model = str(numpy.around(model.min(),decimals=5))
max_area_obs   = str(numpy.around(obs.max(),decimals=5))
min_area_obs   = str(numpy.around(obs.min(),decimals=5))
plt.plot(mons,model,'r',label='model')
plt.plot(mons,obs,'k',label='obs')
ax4.set_xticks(mons)
ax4.set_xlim(1,13)
ax4.set_xticklabels(['A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'])
ax4.set_ylabel('1e6 km^2')
ax4.set_title('Sea Ice Area')
ax4.legend(loc=2,fontsize=10)

ax5= plt.subplot(2,3,5)
mons = numpy.linspace(1,13,13)
model = numpy.array(sic_ext(latitude=(0,90)).sum(axis=1).sum(axis=1))
model = numpy.concatenate((model[7:12], model[0:8]))*1.e-14
obs = numpy.array(sic_ext_nsidc(latitude=(0,90)).sum(axis=1).sum(axis=1))
obs = numpy.concatenate((obs[7:12], obs[0:8]))*1.e-14
max_ext_model = str(numpy.around(model.max(),decimals=5))
min_ext_model = str(numpy.around(model.min(),decimals=5))
max_ext_obs   = str(numpy.around(obs.max(),decimals=5))
min_ext_obs   = str(numpy.around(obs.min(),decimals=5))
plt.plot(mons,model,'r',label='model')
plt.plot(mons,obs,'k',label='obs')
ax5.set_xticks(mons)
ax5.set_xlim(1,13)
ax5.set_xticklabels(['A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A'])
ax5.set_ylabel('1e6 km^2')
ax5.set_title('Sea Ice Extent')
ax5.legend(loc=2,fontsize=10)

plt.subplots_adjust(top=0.8)
fig.text(.5,.92,'Northern Hemisphere Sea Ice',ha='center',fontsize=22,weight='bold')
fig.text(0.5,0.89,'Climatological March Sea Ice Concentration  /  Annual Cycle of Area and Extent',ha='center',fontsize=14)
fig.text(0.5,0.86,expName,ha='center',fontsize=14)
fig.text(0.67,0.39,'Annual Sea Ice Area',fontsize=10)
fig.text(0.67,0.375,'Model Max: ' + max_area_model ,fontsize=10)
fig.text(0.67,0.36, 'Obs Max: ' + max_area_obs ,fontsize=10)
fig.text(0.67,0.345,'Model Min: ' + min_area_model ,fontsize=10)
fig.text(0.67,0.33, 'Obs Min: ' + min_area_obs ,fontsize=10)
fig.text(0.67,0.285,'Annual Sea Ice Extent: ' ,fontsize=10)
fig.text(0.67,0.27, 'Model Max: ' + max_ext_model ,fontsize=10)
fig.text(0.67,0.255,'Obs Max: ' + max_ext_obs ,fontsize=10)
fig.text(0.67,0.24, 'Model Min: ' + min_ext_model ,fontsize=10)
fig.text(0.67,0.225,'Obs Min: ' + min_ext_obs ,fontsize=10)

plt.savefig('seaice.NH.png')

fo.close()

#-- Read NSIDC Climatology (from K. Gottschaldt)
fo = cdms2.open('data/OBS_NSIDC_sat_SH_T2Ms_sic.nc')
lonobs = fo('lon')
latobs = fo('lat')
areacello = fo('areacello')
sic_nsidc = fo('sic',time=('1981-1-1','2000-12-31'))*100.0
cdutil.setTimeBoundsMonthly(sic_nsidc)
sic_nsidc = cdutil.ANNUALCYCLE.climatology(sic_nsidc)
timeObs = sic_nsidc.getAxis(0)
sic_area_nsidc = sic_nsidc*areacello
sic_ext_nsidc  = MV2.where(MV2.greater(sic_nsidc,15.0),100.0,0.0)*areacello

#-- Specify transient grid for obs data
iaxisObs = cdms2.coord.TransientVirtualAxis("i", len(sic_nsidc.getAxis(2)))
jaxisObs = cdms2.coord.TransientVirtualAxis("j", len(sic_nsidc.getAxis(1)))
lonAxisObs = cdms2.coord.TransientAxis2D(lonobs,axes=(jaxisObs, iaxisObs),attributes={'units': 'degrees_east'},id='lon')
latAxisObs = cdms2.coord.TransientAxis2D(latobs,axes=(jaxisObs, iaxisObs),attributes={'units': 'degrees_north'},id='lat')
obsGrid = cdms2.hgrid.TransientCurveGrid(latAxisObs, lonAxisObs, id='lats_lons')
sic_nsidc_grd = MV2.array(sic_nsidc, axes=[timeObs, jaxisObs,iaxisObs], grid=obsGrid, missing=0.0)

#-- Regrid obs to model grid
sic_nsidc_regrid = sic_nsidc_grd.regrid(sic_grd.getGrid(),regridTool='esmf',regridMethod='conserve',fixSrcBounds=True,fixDstBounds=True)
sic_nsidc_regrid.data[:,:,:] = numpy.where(numpy.greater(sic_nsidc_regrid.data,1.e19),0.0,sic_nsidc_regrid.data)
sic_nsidc_regrid.mask = sic.mask

#-- Make a plot
fig = plt.figure(figsize=(11,8.5))
plt.subplot(2,3,1)
ice = sic[8,:,:]
m = Basemap(projection='spstere',boundinglat=-55,lon_0=270,resolution='l')
#m = Basemap(projection='kav7',lon_0=-120,resolution=None)
m.drawmapboundary(fill_color='0.3')
#im0 = m.bluemarble(scale=0.5)
im1 = m.pcolormesh(numpy.minimum(lon,60.),lat,ice,shading='flat',cmap='Blues_r',latlon=True)
plt.clim(0.,100.)
im2 = m.contour(geolon,geolat,ice,[15],latlon=True,colors='red',linewidth=5)
plt.clim(0.0,100.)
m.drawparallels(numpy.arange(-90.,99.,30.))
m.drawmeridians(numpy.arange(-180.,180.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
plt.title('Model - Years ' + year1Label + ' to ' + year2Label)

plt.subplot(2,3,2)
ice = sic_nsidc[8,:,:]
m = Basemap(projection='spstere',boundinglat=-55,lon_0=270,resolution='l')
#m = Basemap(projection='kav7',lon_0=-120,resolution=None)
m.drawmapboundary(fill_color='0.3')
#im0 = m.bluemarble(scale=0.5)
im1 = m.pcolormesh(lonobs,latobs,ice,shading='flat',cmap='Blues_r',latlon=True)
plt.clim(0.,100.)
im2 = m.contour(lonobs,latobs,ice,[15],latlon=True,colors='red',linewidth=5)
plt.clim(0.0,100.)
m.drawparallels(numpy.arange(-90.,99.,30.))
m.drawmeridians(numpy.arange(-180.,180.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
plt.title('NSIDC - Years 1981 to 2000')

plt.subplot(2,3,3)
ice = sic[8,:,:] - sic_nsidc_regrid[8,:,:]
m = Basemap(projection='spstere',boundinglat=-55,lon_0=270,resolution='l')
#m = Basemap(projection='kav7',lon_0=-120,resolution=None)
m.drawmapboundary(fill_color='0.3')
#im0 = m.bluemarble(scale=0.5)
im1 = m.pcolormesh(numpy.minimum(lon,60.),lat,ice,shading='flat',cmap='RdBu',latlon=True)
plt.clim(-100.,100.)
m.drawparallels(numpy.arange(-90.,99.,30.))
m.drawmeridians(numpy.arange(-180.,180.,60.))
cb = m.colorbar(im1,"bottom", size="5%", pad="2%", ticks=[-100.,-60., -30., 0., 30., 60., 100.])
plt.title('Difference')

ax4 = plt.subplot(2,3,4)
mons = numpy.linspace(1,13,13)
model = numpy.array(sic_area(latitude=(-90,0)).sum(axis=1).sum(axis=1))
model = numpy.append(model, model[0])*1.e-14
obs = numpy.array(sic_area_nsidc(latitude=(-90,0)).sum(axis=1).sum(axis=1))
obs = numpy.append(obs, obs[0])*1.e-14
max_area_model = str(numpy.around(model.max(),decimals=5))
min_area_model = str(numpy.around(model.min(),decimals=5))
max_area_obs   = str(numpy.around(obs.max(),decimals=5))
min_area_obs   = str(numpy.around(obs.min(),decimals=5))
plt.plot(mons,model,'r',label='model')
plt.plot(mons,obs,'k',label='obs')
ax4.set_xticks(mons)
ax4.set_xlim(1,13)
ax4.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'])
ax4.set_ylabel('1e6 km^2')
ax4.set_title('Sea Ice Area')
ax4.legend(loc=2,fontsize=10)

ax5= plt.subplot(2,3,5)
mons = numpy.linspace(1,13,13)
model = numpy.array(sic_ext(latitude=(-90,0)).sum(axis=1).sum(axis=1))
model = numpy.append(model, model[0])*1.e-14
obs = numpy.array(sic_ext_nsidc(latitude=(-90,0)).sum(axis=1).sum(axis=1))
obs = numpy.append(obs, obs[0])*1.e-14
max_ext_model = str(numpy.around(model.max(),decimals=5))
min_ext_model = str(numpy.around(model.min(),decimals=5))
max_ext_obs   = str(numpy.around(obs.max(),decimals=5))
min_ext_obs   = str(numpy.around(obs.min(),decimals=5))
plt.plot(mons,model,'r',label='model')
plt.plot(mons,obs,'k',label='obs')
ax5.set_xticks(mons)
ax5.set_xlim(1,13)
ax5.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J'])
ax5.set_ylabel('1e6 km^2')
ax5.set_title('Sea Ice Extent')
ax5.legend(loc=2,fontsize=10)

plt.subplots_adjust(top=0.8)
fig.text(.5,.92,'Southern Hemisphere Sea Ice',ha='center',fontsize=22,weight='bold')
fig.text(0.5,0.89,'Climatological September Sea Ice Concentration  /  Annual Cycle of Area and Extent',ha='center',fontsize=14)
fig.text(0.5,0.86,expName,ha='center',fontsize=14)
fig.text(0.67,0.39,'Annual Sea Ice Area',fontsize=10)
fig.text(0.67,0.375,'Model Max: ' + max_area_model ,fontsize=10)
fig.text(0.67,0.36, 'Obs Max: ' + max_area_obs ,fontsize=10)
fig.text(0.67,0.345,'Model Min: ' + min_area_model ,fontsize=10)
fig.text(0.67,0.33, 'Obs Min: ' + min_area_obs ,fontsize=10)
fig.text(0.67,0.285,'Annual Sea Ice Extent: ' ,fontsize=10)
fig.text(0.67,0.27, 'Model Max: ' + max_ext_model ,fontsize=10)
fig.text(0.67,0.255,'Obs Max: ' + max_ext_obs ,fontsize=10)
fig.text(0.67,0.24, 'Model Min: ' + min_ext_model ,fontsize=10)
fig.text(0.67,0.225,'Obs Min: ' + min_ext_obs ,fontsize=10)

plt.savefig('seaice.SH.png')

EOF

#run the python script
cd ${workDir}
echo "Running Python"
python Krasting-SeaIce.py

#copy output figures to analysis directory
if (! -d ${out_dir}/ice_${yr1}_${yr2}/Krasting.SeaIce) then
 mkdir -p ${out_dir}/ice_${yr1}_${yr2}/Krasting.SeaIce
endif
cd ${workDir}
foreach plot (*.png)
  gcp -v -cd ${plot} ${out_dir}/ice_${yr1}_${yr2}/Krasting.SeaIce/
end

exit
