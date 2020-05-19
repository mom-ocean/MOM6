#!/usr/bin/env python

import argparse
import os
import re
import sys

def parseCommandLine():
  """
  Parse the command line positional and optional arguments.
  This is the highest level procedure invoked from the very end of the script.
  """

  # Arguments
  parser = argparse.ArgumentParser(description='''trailer.py checks Fortran files for trailing white space.''',
      epilog='Written by A.Adcroft, 2017.')
  parser.add_argument('files_or_dirs', type=str, nargs='+',
      metavar='FILE|DIR',
      help='''Fortran files or director in which to search for Fortran files (with .f, .f90, .F90 suffixes).''')
  parser.add_argument('-e','--exclude_dir', type=str, action='append',
      metavar='DIR',
      help='''Exclude directories from search that end in DIR.''')
  parser.add_argument('-l','--line_length', type=int, default=512,
      help='''Maximum allowed length of a line.''')
  parser.add_argument('-s','--source_line_length', type=int, default=132,
      help='''Maximum allowed length of a source line excluding comments.''')
  parser.add_argument('-d','--debug', action='store_true',
      help='turn on debugging information.')
  args = parser.parse_args()

  global debug
  debug = args.debug

  main(args)

def main(args):
  '''
  Does the actual work
  '''
  if (debug): print(args)

  # Process files_or_dirs argument into list of files
  all_files = []
  for a in args.files_or_dirs:
    if os.path.isfile(a): all_files.append(a)
    elif os.path.isdir(a):
      for d,s,files in os.walk(a):
        ignore = False
        if args.exclude_dir is not None:
          for e in args.exclude_dir:
            if e+'/' in d+'/': ignore = True
        if not ignore:
          for f in files:
            _,ext = os.path.splitext(f)
            if ext in ('.f','.F','.f90','.F90'): all_files.append( os.path.join(d,f) )
    else: raise Exception('Argument '+a+' is not a file or directory! Stopping.')
  if (debug): print('Found: ',all_files)

  # For each file, check for trailing white space
  fail = False
  for filename in all_files:
    this = scan_file(filename, line_length=args.line_length, source_line_length=args.source_line_length)
    fail = fail or this
  if fail: sys.exit(1)

def scan_file(filename, line_length=512, source_line_length=132):
  '''Scans file for trailing white space'''
  def msg(filename,lineno,mesg,line=None):
    if line is None: print('%s, line %i: %s'%(filename,lineno,mesg))
    else: print('%s, line %i: %s "%s"'%(filename,lineno,mesg,line))
  white_space_detected = False
  tabs_space_detected = False
  long_line_detected = False
  with open(filename) as file:
    trailing_space = re.compile(r'.* +$')
    tabs = re.compile(r'.*\t.*')
    lineno = 0
    for line in file.readlines():
      lineno += 1
      line = line.replace('\n','')
      srcline = line.split('!', 1)[0] # Discard comments
      if trailing_space.match(line) is not None:
        if debug: print(filename,lineno,line,trailing_space.match(line))
        if len(line.strip())>0: msg(filename,lineno,'Trailing space detected',line)
        else: msg(filename,lineno,'Blank line contains spaces')
        white_space_detected = True
      if tabs.match(line) is not None:
        if len(line.strip())>0: msg(filename,lineno,'Tab detected',line)
        else: msg(filename,lineno,'Blank line contains tabs')
        tabs_space_detected = True
      if len(line)>line_length:
        if len(line.strip())>0: msg(filename,lineno,'Line length exceeded',line)
        else: msg(filename,lineno,'Blank line exceeds line length limit')
        long_line_detected = True
      if len(srcline)>source_line_length:
        msg(filename,lineno,'Non-comment line length exceeded',line)
  return white_space_detected or tabs_space_detected or long_line_detected

# Invoke parseCommandLine(), the top-level procedure
if __name__ == '__main__': parseCommandLine()
