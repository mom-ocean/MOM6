#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

def blend12(d1,d2, f1, f2):
  nj,ni=d1.shape
  x=np.arange(0,nj,dtype=float)/(nj-1)
  x=(x-f1)/(f2-f1)
  x=np.maximum(0.,x)
  x=np.minimum(1.,x)
  weight1d=1. - x
  weight1d=weight1d[:,np.newaxis]
  weight=np.tile( weight1d, (1,ni) )
  return d1*weight+d2*(1-weight)

def blend1234(d1,d2a,d2b,d3,d4a,d4b):
  #d=np.concatenate((d1, blend12(d2a,d2b,0.8,1.0) ),axis=0) # For GIS?
  d=np.concatenate((d1, blend12(d2a,d2b,0.,0.1) ),axis=0) # For CM4
  d=np.concatenate((d,d3),axis=0)
  d=np.concatenate((d, blend12(d4a,d4b,0.3,0.4) ),axis=0)
  return d


#scap_bedmap=nc.Dataset('scap_topog_bedmap2.nc') # For GIS?
scap_bedmap=nc.Dataset('scap_topog_gebco.nc') # For CM4
scap_gebco=nc.Dataset('scap_topog_gebco.nc')
#so_bedmap=nc.Dataset('so_topog_bedmap2.nc') # For GIS?
so_bedmap=nc.Dataset('so_topog_gebco.nc') # For CM4
so_gebco=nc.Dataset('so_topog_gebco.nc')
equator=nc.Dataset('mercator_topog_gebco.nc')
ncap_ibcao=nc.Dataset('ncap_topog.nc')
ncap_gebco=nc.Dataset('ncap_topog_gebco.nc')

mean=blend1234( 
  scap_bedmap.variables['mean'][:],
  so_bedmap.variables['mean'][:],
  so_gebco.variables['mean'][:],
  equator.variables['mean'][:],
  ncap_gebco.variables['mean'][:],
  ncap_ibcao.variables['mean'][:]
  )

std=blend1234( 
  scap_bedmap.variables['std'][:],
  so_bedmap.variables['std'][:],
  so_gebco.variables['std'][:],
  equator.variables['std'][:],
  ncap_gebco.variables['std'][:],
  ncap_ibcao.variables['std'][:]
  )

fout=nc.Dataset('interpolated_topog.nc','w',format='NETCDF3_CLASSIC')

ny=mean.shape[0]; nx = mean.shape[1]

print 'ny,nx= ',ny,nx

yax=fout.createDimension('ny',ny)
xax=fout.createDimension('nx',nx)
fout.createDimension('ntiles',1)

meanv=fout.createVariable('depth','f8',('ny','nx'))
meanv.units='meters'
meanv.standard_name='topographic depth at T-cell centers'
meanv[:]=mean
#maxv=fout.createVariable('max','f8',('ny','nx'))
#maxv.units='meters'
#maxv[:]=max
#minv=fout.createVariable('min','f8',('ny','nx'))
#minv.units='meters'
#minv[:]=min
stdv=fout.createVariable('std','f8',('ny','nx'))
stdv.units='meters'
stdv[:]=std
#countv=fout.createVariable('count','f8',('ny','nx'))
#countv.units='none'
#countv[:]=count

fout.sync()
fout.close()




        
        

        
        
        
