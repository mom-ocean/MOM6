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

    with open(args.expt) as log_expt, open(args.ref) as log_ref:
        clocks_expt = json.load(log_expt)['clocks']
        clocks_ref = json.load(log_ref)['clocks']

        # Gather timers which appear in both clocks
        clock_tags = [clk for clk in clocks_expt if clk in clocks_ref]

        for clk in clock_tags:
            clock_cmp[clk] = {}

            # For now, we only comparge tavg, the rank-averaged timing
            rec = 'tavg'

            t_expt = clocks_expt[clk][rec]
            t_ref = clocks_ref[clk][rec]

            # Compare the relative runtimes
            if all(t > threshold for t in (t_expt, t_ref)):
                dclk = (t_expt - t_ref) / t_ref
            else:
                dclk = 0.
            clock_cmp[clk][rec] = dclk

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

            print('{}{}: {:7.3f}s, {:7.3f}s ({:.1f}%){}'.format(
                ansi_color,
                ' ' * (32 - len(clk)) + clk,
                t_expt,
                t_ref,
                100. * dclk,
                ANSI_RESET,
            ))


if __name__ == '__main__':
    main()
