#!/usr/bin/env python3
"""Generate leave-one-out bases from existing OFFLINE rotated Ritz modes.

This is NOT an online eigensolver. It does not read or modify A.mtx.
The six files are created exclusively and never overwrite inputs.
"""
import argparse
import os
import sys

from mfamg_mode_audit import read_modes
from mfamg_mode_ritz_selection import write_array


def generate(baseline, rotated, negative_control, output_dir):
    n, base = read_modes(baseline)
    nr, full = read_modes(rotated)
    nn, negative = read_modes(negative_control)
    if n != nr or n != nn or len(full) != len(base) + 6:
        raise ValueError('expected exactly six Ritz directions after baseline')
    if len(negative) != len(full) - 1:
        raise ValueError('negative control must have one fewer column')
    if full[:len(base)] != base:
        raise ValueError('rotated full basis changed the original modes')
    # In this snapshot the positive Ritz direction is last. Establish
    # exact equivalence of omit_06 to the prior negative-only control.
    if negative != full[:-1]:
        raise ValueError('negative control differs from full basis minus last column')
    if not os.path.isdir(output_dir):
        raise ValueError('output directory must exist')
    destinations = [os.path.join(output_dir, 'modes_omit_ritz_%02d.mtx' % i)
                    for i in range(1, 7)]
    for path in destinations:
        if os.path.lexists(path):
            raise ValueError('output exists: %s' % path)
    offset = len(base)
    for i, path in enumerate(destinations):
        write_array(path, full[:offset + i] + full[offset + i + 1:])
        print('RITZ_LOO_OMIT_%02d=%s' % (i + 1, path))
    print('RITZ_LOO_N=%d' % n)
    print('RITZ_LOO_BASE_COUNT=%d' % len(base))
    print('RITZ_LOO_EXTRA_COUNT=6')
    print('RITZ_LOO_VARIANTS=6')
    print('RITZ_LOO_EACH_TOTAL_MODES=%d' % (len(full) - 1))
    print('RITZ_LOO_OMIT_06_MATCHES_PREVIOUS_NEGATIVE_CONTROL=YES')
    print('RITZ_LOO_REPLAY_NEW_CASES=01,02,03,04,05')
    print('RITZ_LOO_NOTE=offline_deletion_not_online_discovery')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline', required=True)
    p.add_argument('--rotated', required=True)
    p.add_argument('--negative-control', required=True)
    p.add_argument('--output-dir', required=True)
    args = p.parse_args()
    try:
        generate(args.baseline, args.rotated, args.negative_control,
                 args.output_dir)
    except (ValueError, OSError, OverflowError) as exc:
        print('RITZ_LOO_ERROR=%s' % exc, file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
