#!/usr/bin/env python
import argparse
import collections
import json
import os
import sys

import f90nml

record_type = collections.defaultdict(lambda: float)
for rec in ('grain', 'pemin', 'pemax',):
    record_type[rec] = int


def main():
    desc = 'Parse MOM6 model stdout and return clock data in JSON format.'

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--format', '-f', action='store_true')
    parser.add_argument('--dir', '-d')
    parser.add_argument('log')
    args = parser.parse_args()

    config = {}

    if args.dir:
        # Gather model configuration
        input_nml = os.path.join(args.dir, 'input.nml')
        nml = f90nml.read(input_nml)
        config['input.nml'] = nml.todict()

        parameter_filenames = [
            ('params', 'MOM_parameter_doc.all'),
            ('layout', 'MOM_parameter_doc.layout'),
            ('debug', 'MOM_parameter_doc.debugging'),
        ]
        for key, fname in parameter_filenames:
            config[key] = {}
            with open(os.path.join(args.dir, fname)) as param_file:
                params = parse_mom6_param(param_file)
                config[key].update(params)

    # Get log path
    if os.path.isfile(args.log):
        log_path = args.log
    elif os.path.isfile(os.path.join(args.dir, args.log)):
        log_path = os.path.join(args.dir, args.log)
    else:
        sys.exit('stdout log not found.')

    # Parse timings
    with open(log_path) as log:
        clocks = parse_clocks(log)

    config['clocks'] = clocks

    if args.format:
        print(json.dumps(config, indent=4))
    else:
        print(json.dumps(config))


def parse_mom6_param(param_file):
    params = {}
    for line in param_file:
        param_stmt = line.split('!')[0].strip()
        if param_stmt:
            key, val = [s.strip() for s in param_stmt.split('=')]

            # TODO: Convert to equivalent Python types
            if val in ('True', 'False'):
                params[key] = bool(val)
            else:
                params[key] = val

    return params


def parse_clocks(log):
    clock_start_msg = 'Tabulating mpp_clock statistics across'
    clock_end_msg = 'MPP_STACK high water mark='

    fields = []
    for line in log:
        if line.startswith(clock_start_msg):
            npes = line.lstrip(clock_start_msg).split()[0]

            # Get records
            fields = []
            line = next(log)

            # Skip blank lines
            while line.isspace():
                line = next(log)

            fields = line.split()

            # Exit this loop, begin clock parsing
            break

    clocks = {}
    for line in log:
        # Treat MPP_STACK usage as end of clock reports
        if line.lstrip().startswith(clock_end_msg):
            break

        record = line.split()[-len(fields):]

        clk = line.split(record[0])[0].strip()
        clocks[clk] = {}

        for fld, rec in zip(fields, record):
            rtype = record_type[fld]
            clocks[clk][fld] = rtype(rec)

    return clocks


if __name__ == '__main__':
    main()
