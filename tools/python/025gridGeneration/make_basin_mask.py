#!/usr/bin/env python

from midas.rectgrid import *
from midas.rectgrid_gen import *
import netCDF4
import numpy

def ice9it(i, j, depth, minD=0.):
  """
  Recursive implementation of "ice 9".
  Returns 1 where depth>minD and is connected to depth[j,i], 0 otherwise.
  """
  wetMask = 0*depth      

  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] <= minD: continue
    wetMask[j,i] = 1
  
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) ) # Periodic beyond i=0

    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (j,0) ) # Periodic beyond i=ni-1
    
    if j>0: stack.add((j-1,i))
        
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) ) # Tri-polar fold beyond j=nj-1

  return wetMask

def ice9(x, y, depth, xy0):
  ji = nearestJI(x, y, xy0)
  return ice9it(ji[1], ji[0], depth)

def nearestJI(x, y, (x0, y0)):
  """
  Find (j,i) of cell with center nearest to (x0,y0).
  """
  return numpy.unravel_index( ((x-x0)**2 + (y-y0)**2).argmin() , x.shape)

def southOf(x, y, xy0, xy1):
  """
  Returns 1 for point south/east of the line that passes through xy0-xy1, 0 otherwise.
  """
  x0 = xy0[0]; y0 = xy0[1]; x1 = xy1[0]; y1 = xy1[1]
  dx = x1 - x0; dy = y1 - y0
  Y = (x-x0)*dy - (y-y0)*dx
  Y[Y>=0] = 1; Y[Y<=0] = 0
  return Y

# Rewrite
print 'Reading grid ...',
x = netCDF4.Dataset('ocean_hgrid.nc').variables['x'][1::2,1::2] # Cell centers
y = netCDF4.Dataset('ocean_hgrid.nc').variables['y'][1::2,1::2] # Cell centers
print 'reading topography ...',
depth = netCDF4.Dataset('ocean_topog.nc').variables['depth'][:]
print 'done.'

print 'Generating global wet mask ...',
wet = ice9(x, y, depth, (0,-35)) # All ocean points seeded from South Atlantic
print 'done.'

code = 0*wet

print 'Finding Cape of Good Hope ...',
tmp = 1 - wet; tmp[x<-30] = 0
tmp = ice9(x, y, tmp, (20,-30.))
yCGH = (tmp*y).min()
print 'done.', yCGH

print 'Finding Melbourne ...',
tmp = 1 - wet; tmp[x>-180] = 0
tmp = ice9(x, y, tmp, (-220,-25.))
yMel = (tmp*y).min()
print 'done.', yMel

print 'Processing Persian Gulf ...'
tmp = wet*( 1-southOf(x, y, (55.,23.), (56.5,27.)) )
tmp = ice9(x, y, tmp, (53.,25.))
code[tmp>0] = 11
wet = wet - tmp # Removed named points

print 'Processing Red Sea ...'
tmp = wet*( 1-southOf(x, y, (40.,11.), (45.,13.)) )
tmp = ice9(x, y, tmp, (40.,18.))
code[tmp>0] = 10
wet = wet - tmp # Removed named points

print 'Processing Black Sea ...'
tmp = wet*( 1-southOf(x, y, (26.,42.), (32.,40.)) )
tmp = ice9(x, y, tmp, (32.,43.))
code[tmp>0] = 7
wet = wet - tmp # Removed named points

print 'Processing Mediterranean ...'
tmp = wet*( southOf(x, y, (-5.7,35.5), (-5.7,36.5)) )
tmp = ice9(x, y, tmp, (4.,38.))
code[tmp>0] = 6
wet = wet - tmp # Removed named points

print 'Processing Baltic ...'
tmp = wet*( southOf(x, y, (8.6,56.), (8.6,60.)) )
tmp = ice9(x, y, tmp, (10.,58.))
code[tmp>0] = 9
wet = wet - tmp # Removed named points

print 'Processing Hudson Bay ...'
tmp = wet*( 
           ( 1-(1-southOf(x, y, (-95.,66.), (-83.5,67.5)))
              *(1-southOf(x, y, (-83.5,67.5), (-84.,71.))) 
           )*( 1-southOf(x, y, (-70.,58.), (-70.,65.)) ) )
tmp = ice9(x, y, tmp, (-85.,60.))
code[tmp>0] = 8
wet = wet - tmp # Removed named points

print 'Processing Arctic ...'
tmp = wet*( 
          (1-southOf(x, y, (-171.,66.), (-166.,65.5))) * (1-southOf(x, y, (-64.,66.4), (-50.,68.5))) # Lab Sea
     +    southOf(x, y, (-50.,0.), (-50.,90.)) * (1- southOf(x, y, (0.,65.5), (360.,65.5))  ) # Denmark Strait
     +    southOf(x, y, (-18.,0.), (-18.,65.)) * (1- southOf(x, y, (0.,64.9), (360.,64.9))  ) # Iceland-Sweden
     +    southOf(x, y, (20.,0.), (20.,90.)) # Barents Sea
     +    (1-southOf(x, y, (-280.,55.), (-200.,65.)))
          )
tmp = ice9(x, y, tmp, (0.,85.))
code[tmp>0] = 4
wet = wet - tmp # Removed named points

print 'Processing Pacific ...'
tmp = wet*( (1-southOf(x, y, (0.,yMel), (360.,yMel)))
           -southOf(x, y, (-257,1), (-257,0))*southOf(x, y, (0,3), (1,3))
           -southOf(x, y, (-254.25,1), (-254.25,0))*southOf(x, y, (0,-5), (1,-5)) 
           -southOf(x, y, (-243.7,1), (-243.7,0))*southOf(x, y, (0,-8.4), (1,-8.4)) 
           -southOf(x, y, (-234.5,1), (-234.5,0))*southOf(x, y, (0,-8.9), (1,-8.9)) 
          )
tmp = ice9(x, y, tmp, (-150.,0.))
code[tmp>0] = 3
wet = wet - tmp # Removed named points

print 'Processing Atlantic ...'
tmp = wet*(1-southOf(x, y, (0.,yCGH), (360.,yCGH)))
tmp = ice9(x, y, tmp, (-20.,0.))
code[tmp>0] = 2
wet = wet - tmp # Removed named points

print 'Processing Indian ...'
tmp = wet*(1-southOf(x, y, (0.,yCGH), (360.,yCGH)))
tmp = ice9(x, y, tmp, (55.,0.))
code[tmp>0] = 5
wet = wet - tmp # Removed named points

print 'Processing Southern Ocean ...'
tmp = ice9(x, y, wet, (0.,-55.))
code[tmp>0] = 1
wet = wet - tmp # Removed named points

code[wet>0] = -9
(j,i) = numpy.unravel_index( wet.argmax(), x.shape)
if j:
  print 'There are leftover points unassigned to a basin code'
  while j:
    print x[j,i],y[j,i],[j,i]
    wet[j,i]=0
    (j,i) = numpy.unravel_index( wet.argmax(), x.shape)
else: print 'All points assigned a basin code'

sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)
grid.D=netCDF4.Dataset('ocean_topog.nc').variables['depth'][:]
grid.wet=numpy.zeros(grid.D.shape)
grid.wet[grid.D>0.]=1.0
grid.lath=grid.y_T[:,grid.im/4]  # should not be needed
grid.latq=grid.y_T_bounds[:,grid.im/4+1] # ditto

mask=code

S=state(grid=grid)
var_dict={}
var_dict['X']='Longitude'
var_dict['xax_data']= grid.lonh
var_dict['xunits']= 'degrees_E'
var_dict['Y']='Latitude'
var_dict['yax_data']= grid.lath
var_dict['yunits']= 'degrees_N'
var_dict['Z']=None
var_dict['T']=None
var_dict['_FillValue']=None
var_dict['missing_value']=None
var_dict['flag_values']='1,2,3,4,5,6,7,8,9,10,11'
var_dict['flag_meanings']='1:Southern Ocean, 2:Atlantic Ocean, 3:Pacific Ocean, 4:Arctic Ocean, 5:Indian Ocean, 6:Mediterranean Sea, 7:Black Sea, 8:Hudson Bay, 9:Baltic Sea, 10:Red Sea, 11:Persian Gulf'

S.add_field_from_array(mask,'basin',var_dict=var_dict)

S.write_nc('basin_codes.nc',['basin'])
