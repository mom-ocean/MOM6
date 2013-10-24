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
import shutil as sh


def main():
  global debug # Declared global in order to set it

  # Command line arguments
  parser = argparse.ArgumentParser(description=
       'Applies an "ice 9" algorithm to remove detached water from the topography. Also sets land elevation to 0.',
       epilog='Written by A.Adcroft, 2013.')
  parser.add_argument('filename', type=str,
                      help='netcdf file to read.')
  parser.add_argument('variable', type=str,
                      nargs='?', default='depth',
                      help='Name of variable to plot.')
  parser.add_argument('--output', type=str,
                      nargs='?', default=' ',
                      help='name of the output file. If not specified, "iced_" is prepended to the name of the input file.')
  parser.add_argument('-d','--debug', action='store_true',
                      help='Turn on debugging information.')

  optCmdLineArgs = parser.parse_args()

  debug = False;
  if optCmdLineArgs.debug: debug = True

  nFileName = optCmdLineArgs.output
  if nFileName == ' ': nFileName = 'iced_'+optCmdLineArgs.filename
  applyIce9(optCmdLineArgs.filename, nFileName, optCmdLineArgs.variable, 0., -40.)

def applyIce9(fileName, nFileName, variable, x0, y0):

  if not fileName==nFileName: sh.copyfile(fileName,nFileName)

  # Open netcdf file
  try: rg=Dataset( nFileName, 'a' );
  except: error('There was a problem opening "'+nFileName+'".')
  dir(rg)

  rgVar = rg.variables[variable] # handle to the variable
  dims = rgVar.dimensions # tuple of dimensions
  depth = rgVar[:] # Read the data
  #rg.renameVariable('depth','original_depth')
  print 'Range of input depths: min=',np.amin(depth),'max=',np.amax(depth)

  # A mask based solely on value of depth
  #notLand = np.where( depth<0, 1, 0)
  notLand = ice9it(600,270,depth)
  wet = rg.createVariable('wet','f4',dims)
  wet.long_name = 'Wet/dry mask'
  wet.description = 'Values: 1=Ocean, 0=Land'
  wet[:] = notLand # + (1-notLand)*0.3*np.where( depth<0, 1, 0)

  mwet = rg.createVariable('modified_mask','f4',dims)
  mwet.long_name = 'Modified mask'
  mwet.description = 'Values: 1=Ocean, 0=Land, 0.3 indicates water points removed by "Ice 9" algorithm'
  mwet[:] = notLand + (1-notLand)*0.3*np.where( depth<0, 1, 0)

  #zero = np.where( depth<0 & notLand==0, 1, 0)
  #nDepth = rg.createVariable('depth','f4',dims)
  rgVar[:] = depth*notLand
  rg.sync()

def ice9it(i,j,depth):
  # Iterative implementation of "ice 9"
  wetMask = 0*depth
  (nj,ni) = wetMask.shape
  print ni,nj
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] >= 0: continue
    wetMask[j,i] = 1
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) )
    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (0,j) )
    if j>0: stack.add( (j-1,i) )
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) )
  return wetMask

# Invoke main()
if __name__ == '__main__': main()

