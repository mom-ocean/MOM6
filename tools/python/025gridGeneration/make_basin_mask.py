#!/usr/bin/env python

from midas import *

def ice9it(i,j,depth,wetMask=None,minD=None):
  # Iterative implementation of "ice 9"
  if wetMask is None:
      wetMask = 0*depth      
  if minD is None:
      minD=0.

  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] <= minD: continue
    wetMask[j,i] = 1
  
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) )

    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (j,0) )
    
    if j>0: stack.add((j-1,i))
        
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) )

  return wetMask


sgrid=supergrid(file='ocean_hgrid.nc',cyclic_x=True,tripolar_n=True)
grid=rectgrid(supergrid=sgrid)
grid.D=nc.Dataset('ocean_topog.nc').variables['depth'][:]
grid.wet=np.zeros(grid.D.shape)
grid.wet[grid.D>0.]=1.0
grid.lath=grid.y_T[:,grid.im/4]  # should not be needed
grid.latq=grid.y_T_bounds[:,grid.im/4+1] # ditto


atl_seed = grid.find_geo_bounds(x=(-45.,-45.),y=(25.,25.))
pac_seed = grid.find_geo_bounds(x=(-180.,-180.),y=(0.,0.))
ind_seed = grid.find_geo_bounds(x=(-290.,-290.),y=(-10.,-10.))
arctic_seed = grid.find_geo_bounds(x=(30.,30.),y=(77.,77.))
so_seed = grid.find_geo_bounds(x=(-180.,-180.),y=(-50.,-50.))
hudson_seed = grid.find_geo_bounds(x=(-85.,-85.),y=(60.,60.))
baltic_seed = grid.find_geo_bounds(x=(19.,19.),y=(56.,56.))
med_seed = grid.find_geo_bounds(x=(19.,19.),y=(35.,35.))
black_sea_seed = grid.find_geo_bounds(x=(31.,31.),y=(43.,43.))
red_seed = grid.find_geo_bounds(x=(40.,40.),y=(17.5,17.5))
arabian_seed = grid.find_geo_bounds(x=(53.,53.),y=(25.5,25.5))

atl_30s = grid.find_geo_bounds(x=(-55.,21.),y=(-34.5,-34.5))
ind_30s = grid.find_geo_bounds(x=(21.,117.),y=(-34.5,-34.5))
pac_30s = grid.find_geo_bounds(x=(-210.,-65.),y=(-34.5,-34.5))
bosphorus = grid.find_geo_bounds(x=(29.0,30.0),y=(41.05,41.05))
gibraltar = grid.find_geo_bounds(x=(-5.1,-5.1),y=(35.5,36.5))
bering_st = grid.find_geo_bounds(x=(-173.5,-165.5),y=(66.33,66.33))
davis_st = grid.find_geo_bounds(x=(-70.,-45.),y=(65.,65.))
hudson_bay = grid.find_geo_bounds(x=(-72.5,-72.5),y=(59.,65.))
hudson_bay_n = grid.find_geo_bounds(x=(-82.,-82.),y=(70.5,75.5))
denmark_st = grid.find_geo_bounds(x=(-42.5,-23.),y=(65.7,65.7))
iceland_norway = grid.find_geo_bounds(x=(-20.,15.),y=(65.2,65.2))
baltic = grid.find_geo_bounds(x=(9.,9.),y=(56.9,60.))
itf1 = grid.find_geo_bounds(x=(-235.,-235.),y=(-17.5,7.5))
itf2 = grid.find_geo_bounds(x=(-246.,-235.),y=(-8.6,-8.6))
malacca = grid.find_geo_bounds(x=(-257.2,-257.2),y=(0.5,5.0))
jakarta = grid.find_geo_bounds(x=(-254.6,-254.6),y=(-6.7,-5.7))
rs_entrance = grid.find_geo_bounds(x=(43.9,43.9),y=(12.3,13.0))
as_entrance = grid.find_geo_bounds(x=(56.6,56.6),y=(25.5,27.5))


mask = np.ones(grid.D.shape)

D=grid.D*mask
D[gibraltar[2]:gibraltar[3],gibraltar[0]]=0.  
D[bosphorus[2],bosphorus[0]:bosphorus[1]]=0.  
D[davis_st[2],davis_st[0]:davis_st[1]]=0.
D[atl_30s[2],atl_30s[0]:atl_30s[1]]=0.
D[pac_30s[2],pac_30s[0]:pac_30s[1]]=0.
D[ind_30s[2],ind_30s[0]:]=0.
D[ind_30s[2],0:ind_30s[1]]=0.
D[denmark_st[2],denmark_st[0]:denmark_st[1]]=0.
D[iceland_norway[2],iceland_norway[0]:iceland_norway[1]]=0.
D[hudson_bay[2]:hudson_bay[3],hudson_bay[0]]=0.
D[hudson_bay_n[2]:hudson_bay_n[3],hudson_bay_n[0]]=0.
D[bering_st[2],bering_st[0]:bering_st[1]]=0.
D[baltic[2]:baltic[3],baltic[0]]=0.
D[itf1[2]:itf1[3],itf1[0]]=0.
D[itf2[2],itf2[0]:itf2[1]]=0.   
D[malacca[2]:malacca[3],malacca[0]]=0.
D[jakarta[2]:jakarta[3],jakarta[0]]=0.
D[as_entrance[2]:as_entrance[3],as_entrance[0]]=0.
D[rs_entrance[2]:rs_entrance[3],rs_entrance[0]]=0.

print 'Making Atlantic mask ...'
atl_mask=ice9it(atl_seed[0],atl_seed[2],D,minD=0.)
atl_mask = np.ma.masked_where(atl_mask==0.,atl_mask)
atl_mask[atl_30s[2],atl_30s[0]:atl_30s[1]]=1.
atl_mask[denmark_st[2],denmark_st[0]:denmark_st[1]]=1.
atl_mask[iceland_norway[2],iceland_norway[0]:iceland_norway[1]]=1.
atl_mask[hudson_bay[2]:hudson_bay[3],hudson_bay[0]]=1.
atl_mask[gibraltar[2]:gibraltar[3],gibraltar[0]]=1.  
atl_mask[davis_st[2],davis_st[0]:davis_st[1]]=1.
atl_mask[baltic[2]:baltic[3],baltic[0]]=1.

print 'Making Pacific mask ...'
pac_mask=ice9it(pac_seed[0],pac_seed[2],D,minD=0.)
pac_mask = np.ma.masked_where(pac_mask==0.,pac_mask)
pac_mask[pac_30s[2],pac_30s[0]:pac_30s[1]]=1.
pac_mask[bering_st[2],bering_st[0]:bering_st[1]]=1.
pac_mask[itf1[2]:itf1[3],itf1[0]]=1.
pac_mask[itf2[2],itf2[0]:itf2[1]]=1.

print 'Making Indian mask ...'
ind_mask=ice9it(ind_seed[0],ind_seed[2],D,minD=0.)
ind_mask = np.ma.masked_where(ind_mask==0.,ind_mask)
ind_mask[as_entrance[2]:as_entrance[3],as_entrance[0]]=1.
ind_mask[rs_entrance[2]:rs_entrance[3],rs_entrance[0]]=1.
ind_mask[ind_30s[2],ind_30s[0]:]=1.
ind_mask[ind_30s[2],0:ind_30s[1]]=1.
ind_mask[malacca[2]:malacca[3],malacca[0]]=1.
ind_mask[jakarta[2]:jakarta[3],jakarta[0]]=1.

print 'Making Southern-Ocean mask ...'
so_mask=ice9it(so_seed[0],so_seed[2],D,minD=0.)
so_mask = np.ma.masked_where(so_mask==0.,so_mask)

print 'Making Arctic mask ...'
arctic_mask=ice9it(arctic_seed[0],arctic_seed[2],D,minD=0.)
arctic_mask = np.ma.masked_where(arctic_mask==0.,arctic_mask)
arctic_mask[hudson_bay_n[2]:hudson_bay_n[3],hudson_bay_n[0]]=1.

print 'Making Red Sea mask ...'
rs_mask=ice9it(red_seed[0],red_seed[2],D,minD=0.)
rs_mask = np.ma.masked_where(rs_mask==0.,rs_mask)

print 'Making Persian Gulf mask ...'
as_mask=ice9it(arabian_seed[0],arabian_seed[2],D,minD=0.)
as_mask = np.ma.masked_where(as_mask==0.,as_mask)

print 'Making baltic mask ...'
baltic_mask=ice9it(baltic_seed[0],baltic_seed[2],D,minD=0.)
baltic_mask = np.ma.masked_where(baltic_mask==0.,baltic_mask)

print 'Making Mediterranean mask ...'
med_mask=ice9it(med_seed[0],med_seed[2],D,minD=0.)
med_mask = np.ma.masked_where(med_mask==0.,med_mask)

print 'Making Black Sea mask ...'
black_sea_mask=ice9it(black_sea_seed[0],black_sea_seed[2],D,minD=0.)
black_sea_mask = np.ma.masked_where(black_sea_mask==0.,black_sea_mask)
black_sea_mask[bosphorus[2],bosphorus[0]:bosphorus[1]]=1.  

print 'Making Hudson mask ...'
hudson_mask=ice9it(hudson_seed[0],hudson_seed[2],D,minD=0.)
hudson_mask = np.ma.masked_where(hudson_mask==0.,hudson_mask)

mask[so_mask==1.]=1
mask[atl_mask==1.]=2
mask[pac_mask==1.]=3
mask[arctic_mask==1.]=4
mask[ind_mask==1.]=5
mask[med_mask==1.]=6
mask[black_sea_mask==1.]=7
mask[hudson_mask==1.]=8
mask[baltic_mask==1.]=9
mask[rs_mask==1.]=10
mask[as_mask==1.]=11

mask[grid.wet==0.]=0.


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
var_dict['flag_meanings']='Southern Ocean,Atlantic Ocean, Pacific Ocean, Arctic Ocean, Indian Ocean, Mediterranean Sea, Black Sea, Hudson Bay, Baltic Sea, Red Sea, Persian Gulf'

S.add_field_from_array(mask,'basin',var_dict=var_dict)

S.write_nc('basin_mask.nc',['basin'])
