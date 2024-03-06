#!/usr/bin/env python3

from __future__ import print_function

import argparse
import json
import math

scale = 1e6  # micro-seconds (should make this dynamic)


def display_timing_file(file, show_all):
    """Parse a JSON file of timing results and pretty-print the results"""

    with open(file) as json_file:
        timing_dict = json.load(json_file)

    print("(Times measured in %5.0e seconds)" % (1./scale))
    print("  Min time Module & function")
    for sub in timing_dict.keys():
        tmin = timing_dict[sub]['min'] * scale
        print("%10.4e %s" % (tmin, sub))

        if show_all:
            tmean = timing_dict[sub]['mean'] * scale
            tmax = timing_dict[sub]['max'] * scale
            tstd = timing_dict[sub]['std'] * scale
            nsamp = timing_dict[sub]['n_samples']
            tsstd = tstd / math.sqrt(nsamp)
            print("           (" +
                  "mean = %10.4e " % (tmean) +
                  "±%7.1e, " % (tsstd) +
                  "max = %10.4e, " % (tmax) +
                  "std = %8.2e, " % (tstd) +
                  "# = %d)" % (nsamp))


def compare_timing_files(file, ref, show_all, significance_threshold):
    """Read and compare two JSON files of timing results"""

    with open(file) as json_file:
        timing_dict = json.load(json_file)

    with open(ref) as json_file:
        ref_dict = json.load(json_file)

    print("(Times measured in %5.0e seconds)" % (1./scale))
    print("  Delta (%)  Module & function")
    for sub in {**ref_dict, **timing_dict}.keys():
        T1 = ref_dict.get(sub)
        T2 = timing_dict.get(sub)
        if T1 is not None:
            # stats for reference (old)
            tmin1 = T1['min'] * scale
            tmean1 = T1['mean'] * scale
        if T2 is not None:
            # stats for reference (old)
            tmin2 = T2['min'] * scale
            tmean2 = T2['mean'] * scale
        if (T1 is not None) and (T2 is not None):
            # change in actual minimum as percentage of old
            dt = (tmin2 - tmin1) * 100 / tmin1
            if dt < -significance_threshold:
                color = '\033[92m'
            elif dt > significance_threshold:
                color = '\033[91m'
            else:
                color = ''
            print("%s%+10.4f%%\033[0m  %s" % (color, dt, sub))
        else:
            if T2 is None:
                print("   removed   %s" % (sub))
            else:
                print("     added   %s" % (sub))

        if show_all:
            if T2 is None:
                print("               --")
            else:
                tmax2 = T2['max'] * scale
                tstd2 = T2['std'] * scale
                n2 = T2['n_samples']
                tsstd2 = tstd2 / math.sqrt(n2)
                print("               %10.4e (" % (tmin2) +
                      "mean = %10.4e " % (tmean2) +
                      "±%7.1e, " % (tsstd2) +
                      "max=%10.4e, " % (tmax2) +
                      "std=%8.2e, " % (tstd2) +
                      "# = %d)" % (n2))
            if T1 is None:
                print("               --")
            else:
                tmax1 = T1['max'] * scale
                tstd1 = T1['std'] * scale
                n1 = T1['n_samples']
                tsstd1 = tstd1 / math.sqrt(n1)
                print("               %10.4e (" % (tmin1) +
                      "mean = %10.4e " % (tmean1) +
                      "±%7.1e, " % (tsstd1) +
                      "max=%10.4e, " % (tmax1) +
                      "std=%8.2e, " % (tstd1) +
                      "# = %d)" % (n1))


# Parse arguments
parser = argparse.ArgumentParser(
    description="Beautify timing output from MOM6 timing tests."
)
parser.add_argument(
    'file',
    help="File to process."
)
parser.add_argument(
    '-a', '--all',
    action='store_true',
    help="Display all metrics rather than just the minimum time."
)
parser.add_argument(
    '-t', '--threshold',
    default=6.0, type=float,
    help="Significance threshold to flag (percentage)."
)
parser.add_argument(
    '-r', '--reference',
    help="Reference file to compare against."
)
args = parser.parse_args()

# Do the thing
if args.reference is None:
    display_timing_file(args.file, args.all)
else:
    compare_timing_files(args.file, args.reference, args.all, args.threshold)
