#------------------------------------------------------------------------------
#  refineDiag_data_stager_globalAve.csh
#
#  2014/05/07 JPK
#
#  DESCRIPTION:
#    This script serves two primary functions:
#
#    1.  It unpacks the history tar file to the /ptmp file system.  It allows
#        for more efficient post-processing when individual components are 
#        called by frepp.  (i.e. when the frepp "atmos_month" post-processing
#        script runs, frepp will copy only the unpacked "*atmos_month*" .nc 
#        files from /ptmp to the $work directory rather than the entire history
#        tar file.
#
#    2.  It performs a global annual average of all 3D variables (time, lat, lon)
#        and stores the values in a sqlite database that resides in a parallel
#        directory to the frepp scripts and stdout
#
#------------------------------------------------------------------------------

echo ""
echo ""
echo ""
echo "  ---------- begin refineDiag_data_stager.csh ----------  "

cd $work/$hsmdate
pwd

#-- Unload any previous versions of Python and load the system default
module unload python
module unload cdat
module load python

#-- Unpack gridSpec file.  Right now this hardcoded and this is bad practice.  
#   It would be much better to have the refineDiag script know about the gridSpec location
#   through an already populated FRE variable.  Will talk to Amy L. about alternatives.
#set gridSpecFile = "/archive/cjg/mdt/awg/input/grid/c96_GIS_025_grid_v20140327.tar"
#set gsArchRoot = `echo ${gridSpecFile} | rev | cut -f 2-100 -d '/' | rev`
#set gsBaseName = `basename ${gridSpecFile} | cut -f 1 -d '.'`
#hsmget -v -a ${gsArchRoot} -p /ptmp/$USER/${gsArchRoot} -w `pwd` ${gsBaseName}/\*

#-- Create a directory to house the sqlite database (if it does not already exist)
set localRoot = `echo $scriptName | rev | cut -f 4-100 -d '/' | rev`
if (! -d ${localRoot}/db) then 
  mkdir -p ${localRoot}/db
endif

#-- If db exists, copy it for safe keeping and prevent file locks in the 
foreach reg (global nh sh tropics)
  cp -f ${localRoot}/db/${reg}AveAtmos.db ${localRoot}/db/.${reg}AveAtmos.db
  cp -f ${localRoot}/db/${reg}AveOcean.db ${localRoot}/db/.${reg}AveOcean.db
end

#-- Cat a Python script that performs the averages and writes to a copy of the DB
#   in case it is locked by another user
cat > global_atmos_ave.py <<EOF
import sqlite3, cdms2, cdutil, MV2, numpy, cdtime
import sys

# Set current year
fYear = "${oname}"

fgs1 = cdms2.open(fYear + '.grid_spec.tile1.nc')
fgs2 = cdms2.open(fYear + '.grid_spec.tile2.nc')
fgs3 = cdms2.open(fYear + '.grid_spec.tile3.nc')
fgs4 = cdms2.open(fYear + '.grid_spec.tile4.nc')
fgs5 = cdms2.open(fYear + '.grid_spec.tile5.nc')
fgs6 = cdms2.open(fYear + '.grid_spec.tile6.nc')

geoLat   = MV2.concatenate((MV2.array(fgs1('grid_latt')), MV2.array(fgs2('grid_latt')), MV2.array(fgs3('grid_latt')), MV2.array(fgs4('grid_latt')), MV2.array(fgs5('grid_latt')), MV2.array(fgs6('grid_latt'))),axis=0)
geoLon   = MV2.concatenate((MV2.array(fgs1('grid_lont')), MV2.array(fgs2('grid_lont')), MV2.array(fgs3('grid_lont')), MV2.array(fgs4('grid_lont')), MV2.array(fgs5('grid_lont')), MV2.array(fgs6('grid_lont'))),axis=0)
cellArea = MV2.concatenate((MV2.array(fgs1('area')), MV2.array(fgs2('area')), MV2.array(fgs3('area')), MV2.array(fgs4('area')), MV2.array(fgs5('area')), MV2.array(fgs6('area'))),axis=0)

#Read in 6 nc files
fdata1 = cdms2.open(fYear + '.atmos_month.tile1.nc')
fdata2 = cdms2.open(fYear + '.atmos_month.tile2.nc')
fdata3 = cdms2.open(fYear + '.atmos_month.tile3.nc')
fdata4 = cdms2.open(fYear + '.atmos_month.tile4.nc')
fdata5 = cdms2.open(fYear + '.atmos_month.tile5.nc')
fdata6 = cdms2.open(fYear + '.atmos_month.tile6.nc')

def areaMean(varName,cellArea,geoLat,geoLon,region='global'):
  var = MV2.concatenate((MV2.array(fdata1(varName)), MV2.array(fdata2(varName)), MV2.array(fdata3(varName)), MV2.array(fdata4(varName)), MV2.array(fdata5(varName)), MV2.array(fdata6(varName))),axis=1)
  var = cdutil.YEAR(var).squeeze()
  if (region == 'tropics'):
    var = MV2.masked_where(MV2.logical_or(geoLat < -30., geoLat > 30.),var)
    cellArea = MV2.masked_where(MV2.logical_or(geoLat < -30., geoLat > 30.),cellArea)
  elif (region == 'nh'):
    var  = MV2.masked_where(MV2.less_equal(geoLat,30.),var)
    cellArea  = MV2.masked_where(MV2.less_equal(geoLat,30.),cellArea)
  elif (region == 'sh'):
    var  = MV2.masked_where(MV2.greater_equal(geoLat,-30.),var)
    cellArea  = MV2.masked_where(MV2.greater_equal(geoLat,-30.),cellArea)
  elif (region == 'global'):
    var  = var
    cellArea = cellArea
  res = MV2.array(var*cellArea).sum()/cellArea.sum()
  return res

varDict = fdata1.variables
globalMeanDic={}
tropicsMeanDic={}
nhMeanDic={}
shMeanDic={}
for varName in varDict:
  if (len(varDict[varName].shape) == 3):
    
    conn = sqlite3.connect("${localRoot}/db/.globalAveAtmos.db")
    c = conn.cursor()
    globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='global')
    sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
    sqlres = c.execute(sql)
    sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
    sqlres = c.execute(sql)
    conn.commit()
    c.close()
    conn.close()
    
    conn = sqlite3.connect("${localRoot}/db/.tropicsAveAtmos.db")
    c = conn.cursor()
    globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='tropics')
    sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
    sqlres = c.execute(sql)
    sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
    sqlres = c.execute(sql)
    conn.commit()
    c.close()
    conn.close()
    
    conn = sqlite3.connect("${localRoot}/db/.nhAveAtmos.db")
    c = conn.cursor()
    globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='nh')
    sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
    sqlres = c.execute(sql)
    sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
    sqlres = c.execute(sql)
    conn.commit()
    c.close()
    conn.close()
    
    conn = sqlite3.connect("${localRoot}/db/.shAveAtmos.db")
    c = conn.cursor()
    globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='sh')
    sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
    sqlres = c.execute(sql)
    sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
    sqlres = c.execute(sql)
    conn.commit()
    c.close()
    conn.close()

EOF

cat > global_ocean_ave.py <<EOF
import sqlite3, cdms2, cdutil, MV2, numpy, cdtime
import sys

# Set current year
fYear = "${oname}"

# Test to see if sqlite databse exits, if not, then create it
dbFile = "${localRoot}/db/.globalAveOcean.db"

#Read in gridSpec files
fgs = cdms2.open(fYear + '.ocean_static.nc')

cellArea = fgs('area_t') 
geoLat   = fgs('geolat') 
geoLon   = fgs('geolon') 

#Read in data nc files
fdata = cdms2.open(fYear + '.ocean_month.nc')

def areaMean(varName,cellArea,geoLat,geoLon,region='global'):
  var = fdata(varName)
  var = cdutil.YEAR(var).squeeze()
  if (region == 'tropics'):
    var = MV2.masked_where(MV2.logical_or(geoLat < -30., geoLat > 30.),var)
    cellArea = MV2.masked_where(MV2.logical_or(geoLat < -30., geoLat > 30.),cellArea)
  elif (region == 'nh'):
    var  = MV2.masked_where(MV2.less_equal(geoLat,30.),var)
    cellArea  = MV2.masked_where(MV2.less_equal(geoLat,30.),cellArea)
  elif (region == 'sh'):
    var  = MV2.masked_where(MV2.greater_equal(geoLat,-30.),var)
    cellArea  = MV2.masked_where(MV2.greater_equal(geoLat,-30.),cellArea)
  elif (region == 'global'):
    var  = var
    cellArea = cellArea
  res = MV2.array(var*cellArea).sum()/cellArea.sum()
  return res

varDict = fdata.variables
globalMeanDic={}
tropicsMeanDic={}
nhMeanDic={}
shMeanDic={}
for varName in varDict:
  if (len(varDict[varName].shape) == 3):
    if (fdata(varName).getAxis(1).id == 'yh' and fdata(varName).getAxis(2).id == 'xh'):
      
      conn = sqlite3.connect("${localRoot}/db/.globalAveOcean.db")
      c = conn.cursor()
      globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='global')
      sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
      sqlres = c.execute(sql)
      sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
      sqlres = c.execute(sql)
      conn.commit()
      c.close()
      conn.close()
      
      conn = sqlite3.connect("${localRoot}/db/.tropicsAveOcean.db")
      c = conn.cursor()
      globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='tropics')
      sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
      sqlres = c.execute(sql)
      sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
      sqlres = c.execute(sql)
      conn.commit()
      c.close()
      conn.close()
      
      conn = sqlite3.connect("${localRoot}/db/.nhAveOcean.db")
      c = conn.cursor()
      globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='nh')
      sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
      sqlres = c.execute(sql)
      sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
      sqlres = c.execute(sql)
      conn.commit()
      c.close()
      conn.close()
      
      conn = sqlite3.connect("${localRoot}/db/.shAveOcean.db")
      c = conn.cursor()
      globalMeanDic[varName] = areaMean(varName,cellArea,geoLat,geoLon,region='sh')
      sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
      sqlres = c.execute(sql)
      sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
      sqlres = c.execute(sql)
      conn.commit()
      c.close()
      conn.close()

EOF

#-- Run the averager script
python global_atmos_ave.py
python global_ocean_ave.py

#-- Copy the database back to its original location
foreach reg (global nh sh tropics)
  cp -f ${localRoot}/db/.${reg}AveAtmos.db ${localRoot}/db/${reg}AveAtmos.db
  cp -f ${localRoot}/db/.${reg}AveOcean.db ${localRoot}/db/${reg}AveOcean.db
end

echo "  ---------- end refineDiag_data_stager.csh ----------  "
echo ""
echo ""
echo ""

exit
