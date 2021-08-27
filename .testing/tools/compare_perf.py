#!/usr/bin/env python
import argparse
import json

# Ignore timers below this threshold (in seconds)
DEFAULT_THRESHOLD = 0.05

# Thresholds for reporting
DT_WARN = 0.10  # Slowdown warning
DT_FAIL = 0.25  # Slowdown abort

ANSI_RED = '\033[31m'
ANSI_GREEN = '\033[32m'
ANSI_YELLOW = '\033[33m'
ANSI_RESET = '\033[0m'


def main():
    desc = (
        'Compare two FMS clock output files and report any differences within '
        'a defined threshold.'
    )

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('expt')
    parser.add_argument('ref')
    parser.add_argument('--threshold')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    threshold = float(args.threshold) if args.threshold else DEFAULT_THRESHOLD
    verbose = args.verbose

    clock_cmp = {}

    print('{:35s}{:8s}  {:8s}'.format('', 'Profile', 'Reference'))
    print()

    with open(args.expt) as profile_expt, open(args.ref) as profile_ref:
        perf_expt = json.load(profile_expt)
        perf_ref = json.load(profile_ref)

        events = [ev for ev in perf_expt if ev in perf_ref]

        for event in events:
            # For now, only report the times
            if event not in ('task-clock', 'cpu-clock'):
                continue

            count_expt = perf_expt[event]['count']
            count_ref = perf_ref[event]['count']

            symbols_expt = perf_expt[event]['symbol']
            symbols_ref = perf_ref[event]['symbol']

            symbols = [
                s for s in symbols_expt
                if s in symbols_ref
                and not s.startswith('0x')
            ]

            for symbol in symbols:
                t_expt = float(symbols_expt[symbol]) / 1e9
                t_ref = float(symbols_ref[symbol]) / 1e9

                # Compare the relative runtimes
                if all(t > threshold for t in (t_expt, t_ref)):
                    dclk = (t_expt - t_ref) / t_ref
                else:
                    dclk = 0.

                # Skip trivially low clocks
                if all(t < threshold for t in (t_expt, t_ref)) and not verbose:
                    continue

                # Report the time differences
                ansi_color = ANSI_RESET

                if abs(t_expt - t_ref) > threshold:
                    if dclk > DT_FAIL:
                        ansi_color = ANSI_RED
                    elif dclk > DT_WARN:
                        ansi_color = ANSI_YELLOW
                    elif dclk < -DT_WARN:
                        ansi_color = ANSI_GREEN

                # Remove module name
                sname = symbol.split('_MOD_', 1)[-1]

                # Strip version from glibc calls
                sname = sname.split('@')[0]

                # Remove GCC optimization renaming
                sname = sname.replace('.constprop.0', '')

                if len(sname) > 32:
                    sname = sname[:29] + '...'

                print('{}{}: {:7.3f}s, {:7.3f}s ({:.1f}%){}'.format(
                    ansi_color,
                    ' ' * (32 - len(sname)) + sname,
                    t_expt,
                    t_ref,
                    100. * dclk,
                    ANSI_RESET,
                ))


if __name__ == '__main__':
    main()
