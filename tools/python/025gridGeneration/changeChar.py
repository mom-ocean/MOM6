#!/usr/bin/env python

def error(msg,code=9):
  print 'Error: ' + msg
  exit(code)

# Imports
try: import argparse
except: error('This version of python is not new enough. python 2.7 or newer is required.')
try: from netCDF4 import Dataset, stringtochar
except: error('Unable to import netCDF4 module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_netcdf4')
try: import numpy as np
except: error('Unable to import numpy module. Check your PYTHONPATH.\n'
          +'Perhaps try:\n   module load python_numpy')

def main():

  # Command line arguments
  parser = argparse.ArgumentParser(description=
       'Changes the value of a named char variable in a netcdf file.',
       epilog='Written by A.Adcroft, 2013.')
  parser.add_argument('filename', type=str,
                      help='netcdf file to modify.')
  parser.add_argument('variable', type=str,
                      help='Name of char variable to change.')
  parser.add_argument('value', type=str,
                      help='Contents to change string to.')

  optCmdLineArgs = parser.parse_args()

  rg = Dataset(optCmdLineArgs.filename, 'a' );
  if optCmdLineArgs.variable in rg.variables:
    var = rg.variables[optCmdLineArgs.variable]
    dat = np.empty(1,'S'+repr(len(var)))
    dat[0] = optCmdLineArgs.value
    dc = stringtochar(dat)
    var[:] = dc
  else: error('"'+optCmdLineArgs.variable+'" was not found in "'+optCmdLineArgs.filename+'".')
  rg.close()
  print 'File "%s" updated.'%(optCmdLineArgs.filename)

# Invoke main()
if __name__ == '__main__': main()

