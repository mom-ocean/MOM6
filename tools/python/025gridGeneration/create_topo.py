#!python

#============================================================
# Generate tiles for the northern/southern caps
# and central mercator grid.
#
# python create_topo.py
# Output: mercator_supergrid.nc, ncap_supergrid.nc, scap_supergrid.nc
# These are supergrids (2x grid tracer refinement) containing positions
# cell lengths, areas and angles
#
# Generate topography for grid tiles using BEDMAP for the Antarctic cap
# GEBCO 2 minute data for the Mercator grid and either
# IBCAO or GEBCO for the Northern cap (these files need to be linked to the
# current directory prior to running this command)

# python create_topo.py --tile ncap|scap|mercator 
#
#============================================================



from midas.rectgrid import *
from midas.rectgrid_gen import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.basemap import interp 
import argparse

def csv(value):
   return map(int, value.split(","))

#### Begin User Input

parser = argparse.ArgumentParser()
parser.add_argument('tile',type=str,choices=['ncap','scap','mercator'])
parser.add_argument('--use_ice_sheet', action='store_true', default=False,
        help='Use ice sheet mask for generating coupler mosaic')
parser.add_argument('--use_gebco', action='store_true', default=False,
        help='Use GEBCO for polar regions')


args=parser.parse_args()

tile = args.tile


do_ncap=False
do_mercator=False
do_scap=False
use_ice_sheet_mask=args.use_ice_sheet
use_gebco=args.use_gebco


if tile == 'ncap':
   do_ncap = True
if tile == 'mercator':
   do_mercator = True
if tile == 'scap':   
   do_scap = True


#### End User Input

#### Begin Mercator Grid

mercator=supergrid(file='mercator_supergrid.nc')
mercator_grid=quadmesh(supergrid=mercator,is_latlon=True,cyclic=True)

#### Begin Tripolar Cap

tripolar_n=supergrid(file='ncap_supergrid.nc')
tripolar_n_grid=quadmesh(supergrid=tripolar_n,is_latlon=True,cyclic=True)

#### Begin Antarctic Cap

antarctic_sph=supergrid(file='antarctic_spherical_supergrid.nc')
antarctic_sph_grid=quadmesh(supergrid=antarctic_sph,is_latlon=True,cyclic=True)

antarctic_cap=supergrid(file='scap_supergrid.nc')
antarctic_cap_grid=quadmesh(supergrid=antarctic_cap,is_latlon=True,cyclic=True)


if do_ncap:

######## Interpolate bathymetry from IBCAO to northern cap 
######## on a np stereo projection
      
   if use_gebco:
      ingrid=quadmesh('GEBCO_08_v1.nc',var='depth',simple_grid=True,cyclic=True)
#      np_reg=ingrid.geo_region(y=(mercator.y.max()-1.0,90.0))
      np_reg=ingrid.indexed_region(j=(18360,21599),i=(0,43199))      
      print np_reg['y'][0],np_reg['x'][0]
      print np_reg['y'][-1],np_reg['x'][-1]
      
      TOPO=state('GEBCO_08_v1.nc',grid=ingrid,geo_region=np_reg,fields=['depth'])


      
      TOPO.rename_field('depth','topo')
      TOPO.var_dict['topo']['Ztype']='Fixed'

      import hashlib
      hash=hashlib.md5(TOPO.topo)
      TOPO.write_nc('tmp_new.nc',['topo'])
      print hash.hexdigest()
      
      fnam = 'ncap_topog_gebco.nc'
   
      R=TOPO.subtile('topo',target=tripolar_n_grid)
      R.write_nc(fnam,['mean','max','min','std','count'])

   else:
      xlen=2904000.0*2.0
      x=np.linspace(0.0,xlen,11617)
      X,Y=np.meshgrid(x,x)
      grid_ibcao = quadmesh(lon=X,lat=Y,is_latlon=False,is_cartesian=True)
   
      m = Basemap(projection='stere',width=xlen,height=xlen,lon_0=0.0,lat_0=90.0,resolution='l')

      IBCAO=state('IBCAO_V3_500m_RR.grd',grid=grid_ibcao,fields=['z'])
      IBCAO.rename_field('z','topo')

      xx=tripolar_n_grid.x_T_bounds.copy()
      yy=tripolar_n_grid.y_T_bounds.copy()

      xx[xx>180.]=xx[xx>180.]-360.
      xx[xx<-180.]=xx[xx<-180.]+360.
      x2,y2 = m(xx,yy,inverse=False)
   
      cart_grid_ncap = supergrid(config='cartesian',axis_units='none',xdat=x2,ydat=y2)

      fnam = 'ncap_topog.nc'
      
      R=IBCAO.subtile('topo',target=cart_grid_ncap)
      R.write_nc(fnam,['mean','max','min','std','count'])   

#### Begin Mercator
   
if do_mercator:

   ingrid=quadmesh('GEBCO_08_v1.nc',var='depth',simple_grid=True,cyclic=True)
   merc_reg=ingrid.geo_region(y=(mercator.y.min()-1.0,mercator.y.max()+1.0))
   TOPO=state('GEBCO_08_v1.nc',grid=ingrid,geo_region=merc_reg,fields=['depth'])
   TOPO.rename_field('depth','topo')
   TOPO.var_dict['topo']['Ztype']='Fixed'
      
   fnam = 'mercator_topog_gebco.nc'
   
   R=TOPO.subtile('topo',target=mercator_grid)
   R.write_nc(fnam,['mean','max','min','std','count'])



      




if do_scap:   


   if use_gebco:

      ingrid=quadmesh('GEBCO_08_v1.nc',var='depth',simple_grid=True,cyclic=True)
      sp_reg=ingrid.geo_region(y=(-90.0,antarctic_sph.y.max()+1.0))

      TOPO=state('GEBCO_08_v1.nc',grid=ingrid,geo_region=sp_reg,fields=['depth'])
      TOPO.rename_field('depth','topo')
      TOPO.var_dict['topo']['Ztype']='Fixed'
      
      fnam = 'scap_topog_gebco.nc'
      fnam2 = 'so_topog_gebco.nc'      
   
      R=TOPO.subtile('topo',target=antarctic_cap_grid)
      R.write_nc(fnam,['mean','max','min','std','count'])

      R=TOPO.subtile('topo',target=antarctic_sph_grid)
      R.write_nc(fnam2,['mean','max','min','std','count'])      

      
   else:
      
      wd=6667000.0
      ht=6667000.0
      m = Basemap(projection='stere',width=wd,height=ht,lon_0=0.0,lat_ts=-71.,lat_0=-90.,resolution='l')   

      f=netCDF4.Dataset('bedmap2.nc')
      x1=sq(f.variables['x'][:])*1000 + 3333000.0
      y1=x1
      nx1=len(x1)
      ny1=len(y1)
      wd=6667000.0
      ht=6667000.0
      x1,y1=np.meshgrid(x1,y1)      
      grid_bedmap = quadmesh(lon=x1,lat=y1,is_latlon=False,is_cartesian=True,simple_grid=True)

      if use_ice_sheet_mask:
         TOPO=state('bedmap2.nc',grid=grid_bedmap,fields=['elev_bed','mask_ice','elev_surf','height_gl04c_wgs84'])
         TOPO.elev_bed = TOPO.elev_bed - TOPO.height_gl04c_wgs84          
         TOPO.elev_bed[TOPO.elev_surf>=1.0]=1.0
         TOPO.elev_bed[TOPO.mask_ice>0.0]=1.0
         TOPO.elev_bed=np.ma.masked_where(np.isnan(TOPO.elev_bed),TOPO.elev_bed)
         TOPO.rename_field('elev_bed','topo')
      else:
         TOPO=state('bedmap2.nc',grid=grid_bedmap,fields=['elev_bed','height_gl04c_wgs84'])
         TOPO.elev_bed = TOPO.elev_bed - TOPO.height_gl04c_wgs84          
         TOPO.elev_bed=np.ma.masked_where(np.isnan(TOPO.elev_bed),TOPO.elev_bed)
         TOPO.rename_field('elev_bed','topo')      
      

      fnam = 'scap_topog_bedmap2.nc'
      fnam2 = 'so_topog_bedmap2.nc'         
         
      xx=antarctic_sph_grid.x_T_bounds.copy()
      yy=antarctic_sph_grid.y_T_bounds.copy()

      xx[xx>180.]=xx[xx>180.]-360.
      xx[xx<-180.]=xx[xx<-180.]+360.

      x2,y2 = m(xx,yy,inverse=False)

      cart_grid_so = supergrid(config='cartesian',axis_units='none',xdat=x2,ydat=y2)

      R=TOPO.subtile('topo',target=cart_grid_so)
         
      R.write_nc(fnam2,['mean','max','min','std','count'])   


      xx=antarctic_cap_grid.x_T_bounds.copy()
      yy=antarctic_cap_grid.y_T_bounds.copy()

      xx[xx>180.]=xx[xx>180.]-360.
      xx[xx<-180.]=xx[xx<-180.]+360.
      x2,y2 = m(xx,yy,inverse=False)


      cart_grid_cap = supergrid(config='cartesian',axis_units='none',xdat=x2,ydat=y2)
      
      R=TOPO.subtile('topo',target=cart_grid_cap)
      R.write_nc(fnam,['mean','max','min','std','count'])


