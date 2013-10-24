#!/usr/bin/env python

def error(msg,code=9):
  print 'Error: ' + msg
  exit(code)


# Imports
try: import argparse
except: error('This version of python is not new enough. python 2.7 or newer is required.')
try: from netCDF4 import Dataset
except: error('Unable to import netCDF4 module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_netcdf4')
try: import numpy as np
except: error('Unable to import numpy module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_numpy')
try: import matplotlib.pyplot as plt
except: error('Unable to import matplotlib.pyplot module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_matplotlib')
from matplotlib.widgets import Button, RadioButtons
from matplotlib.colors import LinearSegmentedColormap
import shutil as sh


def main():

  # Command line arguments
  parser = argparse.ArgumentParser(description=
       '''Point-wise editting of topography.
          Button 1 assigns the prescribed level to the cell at the mouse pointer.
          Adjust the prescribed value with buttons on the bottom.
          Double click button 1 assigns the highest of the nearest ocean points.
          Right click on a cell resets to the original value.
          Scroll wheel zooms in and out.
          Move the "data window" around with the North, South, East and West buttons.
          Closing the window writes the file to the output file if one is specified with --output.
        ''',
       epilog='Written by A.Adcroft, 2013.')
  parser.add_argument('filename', type=str,
                      help='netcdf file to read.')
  parser.add_argument('variable', type=str,
                      nargs='?', default='depth',
                      help='Name of variable to edit. Defaults to "depth".')
  parser.add_argument('--output', type=str,
                      nargs='?', default=' ',
                      help='Write an output file. If no output file is specified, creates the file with the "edit_" prepended to the name  of the input file.')

  optCmdLineArgs = parser.parse_args()

  createGUI(optCmdLineArgs.filename, optCmdLineArgs.variable, optCmdLineArgs.output)


def createGUI(fileName, variable, outFile):

  # Open netcdf file
  try: rg=Dataset( fileName, 'r' );
  except: error('There was a problem opening "'+fileName+'".')

  rgVar = rg.variables[variable] # handle to the variable
  dims = rgVar.dimensions # tuple of dimensions
  depth = rgVar[:] # Read the data
  #depth = depth[0:600,0:600]
  (nj,ni) = depth.shape
  print 'Range of input depths: min=',np.amin(depth),'max=',np.amax(depth)

  try:
    sg=Dataset( 'supergrid.nc', 'r' );
    lon = sg.variables['x'][:]; lon = lon[0:2*nj+1:2,0:2*ni+1:2]
    lat = sg.variables['y'][:]; lat = lat[0:2*nj+1:2,0:2*ni+1:2]
  except:
    lon, lat = np.meshgrid( np.arange(ni+1), np.arange(nj+1) )
  fullData = Topography(lon, lat, depth)

  class Container:
    def __init__(self):
      self.view = None
      self.edits = None
      self.data = None
      self.quadMesh = None
      self.ax = None
      self.syms = None
      cdict = {'red': ((0.0, 0.0, 0.0), (0.5, 0.7, 0.0), (1.0, 0.9, 0.0)),
           'green': ((0.0, 0.0, 0.0), (0.5, 0.7, 0.2), (1.0, 1.0, 0.0)),
            'blue': ((0.0, 0.0, 0.2), (0.5, 1.0, 0.0), (1.0, 0.9, 0.0))}
      self.cmap = LinearSegmentedColormap('my_colormap',cdict,256)
      self.clim = 6000
      self.climLabel = None
  All = Container()
  All.view = View(ni,nj)
  All.edits = Edits()

  # Read edit data, if it exists
  if 'iEdit' in rg.variables:
    jEdit = rg.variables['iEdit'][:]; iEdit = rg.variables['jEdit'][:]
    zEdit = rg.variables['zEdit'][:]
    for l,i in enumerate(iEdit):
      All.edits.setVal( fullData.height[iEdit[l],jEdit[l]] )
      fullData.height[iEdit[l],jEdit[l]] = zEdit[l] # Restore data
      All.edits.add( iEdit[l],jEdit[l] )
  All.data = fullData.cloneWindow( (All.view.i0,All.view.j0), (All.view.iw,All.view.jw) )
  if All.edits.ijz: All.data.applyEdits(fullData, All.edits.ijz)

  # A mask based solely on value of depth
  #notLand = np.where( depth<0, 1, 0)
  #wet = ice9it(600,270,depth)

  All.quadMesh = plt.pcolormesh(All.data.longitude,All.data.latitude,All.data.height,cmap=All.cmap,vmin=-All.clim,vmax=All.clim)
  All.syms = All.edits.plot(fullData)
  dir(All.syms)
  All.ax=plt.gca(); All.ax.set_xlim( All.data.xlim ); All.ax.set_ylim( All.data.ylim )
  All.climLabel = plt.figtext(.97,.97, 'XXXXX', ha='right', va='top')
  All.climLabel.set_text('clim = $\pm$%i'%(All.clim))
  All.edits.label = plt.figtext(.97,.03, 'XXXXX', ha='right', va='bottom')
  All.edits.label.set_text('New depth = %i'%(All.edits.get()))
  lowerButtons = Buttons()
  def resetDto0(event): All.edits.setVal(0)
  lowerButtons.add('Set 0', resetDto0)
  def resetDto100(event): All.edits.addToVal(100)
  lowerButtons.add('+100', resetDto100)
  def resetDto100(event): All.edits.addToVal(30)
  lowerButtons.add('+30', resetDto100)
  def resetDto100(event): All.edits.addToVal(10)
  lowerButtons.add('+10', resetDto100)
  def resetDto100(event): All.edits.addToVal(3)
  lowerButtons.add('+3', resetDto100)
  def resetDto100(event): All.edits.addToVal(1)
  lowerButtons.add('+1', resetDto100)
  def resetDto100(event): All.edits.addToVal(-1)
  lowerButtons.add('-1', resetDto100)
  def resetDto100(event): All.edits.addToVal(-3)
  lowerButtons.add('-3', resetDto100)
  def resetDto100(event): All.edits.addToVal(-10)
  lowerButtons.add('-10', resetDto100)
  def resetDto100(event): All.edits.addToVal(-30)
  lowerButtons.add('-30', resetDto100)
  def resetDto100(event): All.edits.addToVal(-100)
  lowerButtons.add('-100', resetDto100)
  def resetDto100(event): All.edits.addToVal(-500)
  lowerButtons.add('-500', resetDto100)
  def undoLast(event):
    All.edits.pop()
    All.data = fullData.cloneWindow( (All.view.i0,All.view.j0), (All.view.iw,All.view.jw) )
    All.data.applyEdits(fullData, All.edits.ijz)
    All.quadMesh.set_array(All.data.height.ravel())
    All.edits.updatePlot(fullData,All.syms)
    plt.draw()
  lowerButtons.add('Undo', undoLast)
  upperButtons = Buttons(bottom=1-.0615)
  def colorScale(event):
    Levs = [50, 200, 1000, 6000]
    i = Levs.index(All.clim)
    if event=='+clim': i = min(i+1, len(Levs)-1)
    elif event==' -clim': i = max(i-1, 0)
    All.clim = Levs[i]
    #All.quadMesh = plt.pcolormesh(All.data.longitude,All.data.latitude,All.data.height,cmap=All.cmap,vmin=-All.clim,vmax=All.clim)
    #All.ax.set_xlim( All.data.xlim ); All.ax.set_ylim( All.data.ylim )
    All.quadMesh.set_clim(vmin=-All.clim, vmax=All.clim)
    All.climLabel.set_text('clim = $\pm$%i'%(All.clim))
    plt.draw()
  def moveVisData(di,dj):
    All.view.move(di,dj)
    All.data = fullData.cloneWindow( (All.view.i0,All.view.j0), (All.view.iw,All.view.jw) )
    All.data.applyEdits(fullData, All.edits.ijz)
    plt.sca(All.ax); plt.cla()
    All.quadMesh = plt.pcolormesh(All.data.longitude,All.data.latitude,All.data.height,cmap=All.cmap,vmin=-All.clim,vmax=All.clim)
    All.ax.set_xlim( All.data.xlim ); All.ax.set_ylim( All.data.ylim )
    All.syms = All.edits.plot(fullData)
    plt.draw()
  def moveWindowLeft(event): moveVisData(-1,0)
  upperButtons.add('West', moveWindowLeft)
  def moveWindowRight(event): moveVisData(1,0)
  upperButtons.add('East', moveWindowRight);
  def moveWindowDown(event): moveVisData(0,-1)
  upperButtons.add('South', moveWindowDown)
  def moveWindowUp(event): moveVisData(0,1)
  upperButtons.add('North', moveWindowUp)
  climButtons = Buttons(bottom=1-.0615,left=0.65)
  def incrCScale(event): colorScale('+clim')
  climButtons.add('Incr', incrCScale)
  def incrCScale(event): colorScale(' -clim')
  climButtons.add('Decr', incrCScale)
  plt.sca(All.ax)
  def onClick(event): # Mouse button click
    if event.inaxes==All.ax and event.button==1 and event.xdata:
      (i,j) = findPointInMesh(fullData.longitude, fullData.latitude, event.xdata, event.ydata)
      if not i==None:
        (I,J) = findPointInMesh(All.data.longitude, All.data.latitude, event.xdata, event.ydata)
        if event.dblclick:
          nVal = -99999
          if All.data.height[I+1,J]<0: nVal = max(nVal, All.data.height[I+1,J])
          if All.data.height[I-1,J]<0: nVal = max(nVal, All.data.height[I-1,J])
          if All.data.height[I,J+1]<0: nVal = max(nVal, All.data.height[I,J+1])
          if All.data.height[I,J-1]<0: nVal = max(nVal, All.data.height[I,J-1])
          if nVal==-99999: return
          All.edits.add(i,j,nVal)
          All.data.height[I,J] = nVal
        else:
          All.edits.add(i,j)
          All.data.height[I,J] = All.edits.get()
        All.quadMesh.set_array(All.data.height.ravel())
        All.edits.updatePlot(fullData,All.syms)
        plt.draw()
    elif event.inaxes==All.ax and event.button==3 and event.xdata:
      (i,j) = findPointInMesh(fullData.longitude, fullData.latitude, event.xdata, event.ydata)
      if not i==None:
        All.edits.delete(i,j)
        All.data = fullData.cloneWindow( (All.view.i0,All.view.j0), (All.view.iw,All.view.jw) )
        All.data.applyEdits(fullData, All.edits.ijz)
        All.quadMesh.set_array(All.data.height.ravel())
        All.edits.updatePlot(fullData,All.syms)
        plt.draw()
    elif event.inaxes==All.ax and event.button==2 and event.xdata: zoom(event) # Re-center
  plt.gcf().canvas.mpl_connect('button_press_event', onClick)
  def zoom(event): # Scroll wheel up/down
    if event.button == 'up': scale_factor = 1/1.5 # deal with zoom in
    elif event.button == 'down': scale_factor = 1.5 # deal with zoom out
    else: scale_factor = 1.0
    new_xlim, new_ylim = newLims( \
           All.ax.get_xlim(), All.ax.get_ylim(), (event.xdata,event.ydata), \
           All.data.xlim, All.data.ylim, scale_factor)
    if not new_xlim: return # No changein limits
    All.ax.set_xlim(new_xlim[0], new_xlim[1]); All.ax.set_ylim(new_ylim[0], new_ylim[1])
    plt.draw() # force re-draw
  plt.gcf().canvas.mpl_connect('scroll_event', zoom)
  def statusMesg(x,y):
    j,i = findPointInMesh(fullData.longitude, fullData.latitude, x, y)
    if not i==None: return 'lon,lat=%.2f,%.2f  depth(%i,%i)=%.2f'%(x,y,i,j,fullData.height[j,i])
    else: return 'lon,lat=%.3f,%.3f'%(x,y)
  All.ax.format_coord = statusMesg
  plt.show()
  All.edits.list()
  if not outFile: outFile = 'edit_'+fileName
  if not outFile==' ':
    print 'Creating new file "'+outFile+'"'
    # Create new netcdf file
    if not fileName==outFile: sh.copyfile(fileName,outFile)
    try: rg=Dataset( outFile, 'r+' );
    except: error('There was a problem opening "'+outFile+'".')
    rgVar = rg.variables[variable] # handle to the variable
    dims = rgVar.dimensions # tuple of dimensions
    rgVar[:] = fullData.height[:,:] # Write the data
    if All.edits.ijz:
      print 'Applying %i edits'%(len(All.edits.ijz))
      if 'nEdits' in rg.dimensions:
        numEdits = rg.dimensions['nEdits']
      else: numEdits = rg.createDimension('nEdits', 0)#len(All.edits.ijz))
      if 'iEdit' in rg.variables: iEd = rg.variables['iEdit']
      else:
        iEd = rg.createVariable('iEdit','i4',('nEdits',))
        iEd.long_name = 'i-index of edited data'
      if 'jEdit' in rg.variables: jEd = rg.variables['jEdit']
      else:
        jEd = rg.createVariable('jEdit','i4',('nEdits',))
        jEd.long_name = 'j-index of edited data'
      if 'zEdit' in rg.variables: zEd = rg.variables['zEdit']
      else:
        zEd = rg.createVariable('zEdit','f4',('nEdits',))
        zEd.long_name = 'Original value of edited data'
        zEd.units = rgVar.units
      for l,(i,j,z) in enumerate(All.edits.ijz):
        iEd[l] = j; jEd[l] = i; zEd[l] = rgVar[i,j]; rgVar[i,j] = z
    rg.close()


def ice9it(i,j,depth):
  # Iterative implementation of "ice 9"
  wetMask = 0*depth
  (ni,nj) = wetMask.shape
  stack = set()
  stack.add( (i,j) )
  while stack:
    (i,j) = stack.pop()
    if wetMask[i,j] or depth[i,j] >= 0: continue
    wetMask[i,j] = 1
    if i>0: stack.add( (i-1,j) )
    else: stack.add( (ni-1,j) )
    if i<ni-1: stack.add( (i+1,j) )
    else: stack.add( (0,j) )
    if j>0: stack.add( (i,j-1) )
    if j<nj-1: stack.add( (i,j+1) )
  return wetMask


def findPointInMesh(meshX, meshY, pointX, pointY):
  def sign(x):
    if x>0: return 1.0
    elif x<0: return -1.0
    else: return 0.
  def crossProd(u0,v0,u1,v1):
    return sign( u0*v1 - u1*v0 )
  def isPointInConvexPolygon(pX, pY, p):
    u0 = pX[0]-pX[-1]; v0 = pY[0]-pY[-1]
    u1 = pX[-1] - p[0]; v1 = pY[-1] - p[1]
    firstSign = crossProd(u0,v0,u1,v1)
    for n in range(len(pX)-1):
      u0 = pX[n+1]-pX[n]; v0 = pY[n+1]-pY[n]
      u1 = pX[n] - p[0]; v1 = pY[n] - p[1]
      if crossProd(u0,v0,u1,v1)*firstSign<0: return False
    return True 
  def recurIJ(mX, mY, p, ij00, ij22):
    # Unpack indices
    i0 = ij00[0]; i2 = ij22[0]; j0 = ij00[1]; j2 = ij22[1];
    # Test bounding box first (bounding box is larger than polygon)
    xmin=min( np.amin(mX[i0,j0:j2]), np.amin(mX[i2,j0:j2]), np.amin(mX[i0:i2,j0]), np.amin(mX[i0:i2,j2]) )
    xmax=max( np.amax(mX[i0,j0:j2]), np.amax(mX[i2,j0:j2]), np.amax(mX[i0:i2,j0]), np.amax(mX[i0:i2,j2]) )
    ymin=min( np.amin(mY[i0,j0:j2]), np.amin(mY[i2,j0:j2]), np.amin(mY[i0:i2,j0]), np.amin(mY[i0:i2,j2]) )
    ymax=max( np.amax(mY[i0,j0:j2]), np.amax(mY[i2,j0:j2]), np.amax(mY[i0:i2,j0]), np.amax(mY[i0:i2,j2]) )
    if p[0]<xmin or p[0]>xmax or p[1]<ymin or p[1]>ymax: return None, None
    if i2>i0+1:
      i1=int(0.5*(i0+i2))
      if j2>j0+1: # Four quadrants to test
        j1=int(0.5*(j0+j2))
        iAns, jAns = recurIJ(mX, mY, p, (i0,j0), (i1,j1))
        if iAns==None: iAns, jAns = recurIJ(mX, mY, p, (i1,j1), (i2,j2))
        if iAns==None: iAns, jAns = recurIJ(mX, mY, p, (i0,j1), (i1,j2))
        if iAns==None: iAns, jAns = recurIJ(mX, mY, p, (i1,j0), (i2,j1))
      else: # Two halves, east/west, to test
        j1=int(0.5*(j0+j2))
        iAns, jAns = recurIJ(mX, mY, p, (i0,j0), (i1,j2))
        if iAns==None: iAns, jAns = recurIJ(mX, mY, p, (i1,j0), (i2,j2))
    else:
      if j2>j0+1: # Two halves, north/south, to test
        j1=int(0.5*(j0+j2))
        iAns, jAns = recurIJ(mX, mY, p, (i0,j0), (i2,j1))
        if iAns==None: iAns, jAns = recurIJ(mX, mY, p, (i0,j1), (i2,j2))
      else: # Only one cell left (based on the bounding box)
        if not isPointInConvexPolygon( \
              [mX[i0,j0],mX[i0+1,j0],mX[i0+1,j0+1],mX[i0,j0+1]], \
              [mY[i0,j0],mY[i0+1,j0],mY[i0+1,j0+1],mY[i0,j0+1]], \
              p): return None, None
        return i0,j0
    return iAns, jAns
  (ni,nj) = meshX.shape; ij00 = [0, 0]; ij22 = [ni-1, nj-1]
  return recurIJ(meshX, meshY, (pointX, pointY), ij00, ij22)


# Calculate a new window by scaling the current window, centering
# on the cursor if possible.
def newLims(cur_xlim, cur_ylim, cursor, xlim, ylim, scale_factor):
  cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
  cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
  xdata = cursor[0]; ydata = cursor[1]
  new_xrange = cur_xrange*scale_factor; new_yrange = cur_yrange*scale_factor
  xdata = min( max( xdata, xlim[0]+new_xrange ), xlim[1]-new_xrange )
  ydata = min( max( ydata, ylim[0]+new_yrange ), ylim[1]-new_yrange )
  xL = max( xlim[0], xdata - new_xrange ); xR = min( xlim[1], xdata + new_xrange )
  yL = max( ylim[0], ydata - new_yrange ); yR = min( ylim[1], ydata + new_yrange )
  if xL==cur_xlim[0] and xR==cur_xlim[1] and \
     yL==cur_ylim[0] and yR==cur_ylim[1]: return None, None
  return (xL, xR), (yL, yR)


# Class to handle adding buttons to GUI
class Buttons:
  scale = 0.014; space = .01
  def __init__(self, bottom=.015, left=.015):
    self.leftEdge = left
    self.bottomEdge = bottom
    self.height = .05
    self.list = []
  def add(self,label,fn): # fn is callback
    width = self.scale*len(label)
    np = [ self.leftEdge, self.bottomEdge, width, self.height ]
    self.leftEdge = self.leftEdge + width + self.space
    button = Button(plt.axes(np),label); button.on_clicked( fn )
    self.list.append( button )


# Class to contain edits
class Edits:
  def __init__(self):
    self.newDepth = 0
    self.ijz = []
    self.label = None # Handle to text box
  def setVal(self, newVal):
    self.newDepth = newVal
    if self.label: self.label.set_text('New depth = %i'%(self.newDepth))
  def addToVal(self, increment):
    self.newDepth += increment
    if self.label: self.label.set_text('New depth = %i'%(self.newDepth))
    plt.draw()
  def get(self): return self.newDepth
  def delete(self,i,j):
    for I,J,D in self.ijz:
      if (i,j)==(I,J): self.ijz.remove((I,J,D))
  def add(self,i,j,nVal=None):
    self.delete(i,j)
    if not nVal==None: self.ijz.append( (i,j,nVal) )
    else: self.ijz.append( (i,j,self.newDepth) )
  def pop(self):
    if self.ijz: self.ijz.pop()
  def list(self):
    for a in self.ijz: print a
  def plot(self,topo):
    x = []; y= []
    for i,j,z in self.ijz:
      tx,ty = topo.cellCoord(j,i)
      if tx:
        x.append(tx); y.append(ty)
    if x:
      h, = plt.plot(x, y, 'ro', hold=True)
      return h
    else: return None
  def updatePlot(self,topo,h):
    x = []; y= []
    for i,j,z in self.ijz:
      tx,ty = topo.cellCoord(j,i)
      if tx:
        x.append(tx); y.append(ty)
    if x:
      h.set_xdata(x); h.set_ydata(y)


# Class to contain data
class Topography:
  def __init__(self, lon, lat, height):
    self.longitude = lon
    self.latitude = lat
    self.height = np.copy( height )
    self.xlim = ( np.min(lon), np.max(lon) )
    self.ylim = ( np.min(lat), np.max(lat) )
  def cloneWindow(self, (i0,j0), (iw,jw)):
    i1 = i0 + iw; j1 = j0 + jw
    return Topography( self.longitude[j0:j1+1,i0:i1+1], \
                       self.latitude[j0:j1+1,i0:i1+1], \
                       self.height[j0:j1,i0:i1] )
  def applyEdits(self, origData, ijz):
    for i,j,z in ijz:
      x = (origData.longitude[i,j] + origData.longitude[i+1,j+1])/2.
      y = (origData.latitude[i,j] + origData.latitude[i+1,j+1])/2.
      (I,J) = findPointInMesh(self.longitude, self.latitude, x, y)
      if not I==None: self.height[I,J] = z
  def cellCoord(self,j,i):
    #ni, nj = self.longitude.shape
    #if i<0 or j<0 or i>=ni-1 or j>=nj-1: return None, None
    x = (self.longitude[i,j] + self.longitude[i+1,j+1])/2.
    y = (self.latitude[i,j] + self.latitude[i+1,j+1])/2.
    return x,y

# CLass to record the editing window
class View:
  def __init__(self, ni, nj):
    self.ni = ni
    self.nj = nj
    self.i0 = 0
    self.j0 = 0
    self.iw = min(512,ni)
    self.jw = min(512,nj)
  def move(self, di, dj):
    self.i0 = min( max(0, self.i0+int(di*self.iw/2.)), self.ni-self.iw)
    self.j0 = min( max(0, self.j0+int(dj*self.jw/2.)), self.nj-self.jw)
  def geti(self): return (self.i0,self.i0+self.iw)
  def getj(self): return (self.j0,self.j0+self.jw)
# Invoke main()
if __name__ == '__main__': main()

