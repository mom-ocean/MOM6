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

def main():

  # Command line arguments
  parser = argparse.ArgumentParser(description=
       'Adds the named dimension to a netcdf file.',
       epilog='Written by A.Adcroft, 2013.')
  parser.add_argument('filename', type=str,
                      help='netcdf file to modify.')
  parser.add_argument('dimension', type=str,
                      nargs='?', default='ntiles',
                      help='Name of dimension to add.')
  parser.add_argument('value', type=int,
                      nargs='?', default=1,
                      help='Value to set dimension length to.')

  optCmdLineArgs = parser.parse_args()

  rg = Dataset(optCmdLineArgs.filename, 'a' );
  rg.createDimension(optCmdLineArgs.dimension,optCmdLineArgs.value)
  rg.close()
  print 'File "%s" updated.'%(optCmdLineArgs.filename)

# Invoke main()
if __name__ == '__main__': main()

