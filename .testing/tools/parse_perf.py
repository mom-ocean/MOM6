#!/usr/bin/env python
import argparse
import collections
import json
import os
import re
import shlex
import subprocess
import sys

perf_scanner = re.Scanner([
  (r'<', lambda scanner, token: token),
  (r'>', lambda scanner, token: token),
  (r'\(', lambda scanner, token: token),
  (r'\)', lambda scanner, token: token),
  (r'[ \t]+', lambda scanner, token: token),
  (r'[^<>() \t]+', lambda scanner, token: token),
])


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
                tokens, remainder = perf_scanner.scan(line)
                if remainder:
                    print('Line could not be tokenized', file=sys.stderr)
                    print(' line:', repr(line), file=sys.stderr)
                    print(' tokens:', tokens, file=sys.stderr)
                    print(' remainder:', remainder, file=sys.stderr)
                    sys.exit(os.EX_DATAERR)

                # Construct record from tokens
                # (NOTE: Not a proper grammar, just dumb bracket counting)
                record = []
                bracks = 0
                parens = 0

                for tok in tokens:
                    if tok == '<':
                        bracks += 1

                    if tok == '(':
                        parens += 1

                    rec = record[-1] if record else None

                    inside_bracket = rec and (bracks > 0 or parens > 0)
                    lead_rec = tok in '<(' and rec and not rec.isspace()
                    tail_rec = not tok.isspace() and rec and rec[-1] in '>)'

                    if inside_bracket or lead_rec or tail_rec:
                        record[-1] += tok
                    else:
                        record.append(tok)

                    if tok == '>':
                        bracks -= 1
                    if tok == ')':
                        parens -= 1

                # Strip any whitespace tokens
                record = [rec for rec in record if not rec.isspace()]

                try:
                    symbol = record[2]
                    period = int(record[3])
                except:
                    print("parse_perf.py: Error extracting symbol count",
                          file=sys.stderr)
                    print("line:", repr(line), file=sys.stderr)
                    print("tokens:", tokens, file=sys.stderr)
                    print("record:", record, file=sys.stderr)
                    raise

                profile[event_name]['symbol'][symbol] = period

    return profile


if __name__ == '__main__':
    main()
