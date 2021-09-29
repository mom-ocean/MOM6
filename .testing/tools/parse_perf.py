#!/usr/bin/env python
import argparse
import collections
import json
import os
import shlex
import subprocess
import sys


def main():
    desc = 'Parse perf.data and return in JSON format.'

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--format', '-f', action='store_true')
    parser.add_argument('data')
    args = parser.parse_args()

    profile = parse_perf_report(args.data)

    if args.format:
        print(json.dumps(profile, indent=4))
    else:
        print(json.dumps(profile))


def parse_perf_report(perf_data_path):
    profile = {}

    cmd = shlex.split(
        'perf report -s symbol,period -i {}'.format(perf_data_path)
    )
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True) as proc:
        event_name = None
        for line in proc.stdout:
            # Skip blank lines:
            if not line or line.isspace():
                continue

            # Set the current event
            if line.startswith('# Samples: '):
                event_name = line.split()[-1].strip("'")

                # Remove perf modifiers for now
                event_name = event_name.rsplit(':', 1)[0]

                profile[event_name] = {}
                profile[event_name]['symbol'] = {}

            # Get total count
            elif line.startswith('# Event count '):
                event_count = int(line.split()[-1])
                profile[event_name]['count'] = event_count

            # skip all other 'comment' lines
            elif line.startswith('#'):
                continue

            # get per-symbol count
            else:
                tokens = line.split()
                symbol = tokens[2]
                period = int(tokens[3])

                profile[event_name]['symbol'][symbol] = period

    return profile


if __name__ == '__main__':
    main()
