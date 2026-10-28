#!/usr/bin/env python3
"""
cluster/06_depth/glmm/01_compute_percentile_cutoffs.py
======================================================

Step 1 of the callable-depth pipeline (cluster/README.md step 8; Methods
"Callable-position read-depth analysis"). Input: pools27_v4.sanitized.mpileup.gz
(_toplevel/mpileup_v4_merge.sh). Output: pool_chrom_cutoffs_v4.tsv.

Streams through the v4 mpileup once. For every pool x chromosome combination,
estimates the 99th percentile of read depth using t-digest (small-memory
streaming quantile approximation, sub-percent accuracy on the 99th percentile).

Output is a small TSV: one row per (pool, chrom), columns pool, chrom,
max_depth_cutoff. This becomes the lookup table for step 2's max-cov filter.

Why per-chromosome rather than genome-wide? Chromosomes can differ in GC
content, repeat density, and baseline mappability; a genome-wide percentile
biases against repeat-rich chromosomes. PoolSNP's `max-cov 0.99` filter is
also computed per-chromosome, so this matches PoolSNP's definition of
"callable position".

Why t-digest rather than exact? Exact would require collecting every depth
observation per (pool, chrom) in memory. For chr2 (~204 Mb) x 27 pools
that's ~22 GB of int32 values held at once. t-digest gives sub-percent
accuracy at the 99th percentile with negligible memory.

Usage
-----
    python 01_compute_percentile_cutoffs.py \
        --mpileup pools27_v4.sanitized.mpileup.gz \
        --locations locations_v4.txt \
        --out pool_chrom_cutoffs_v4.tsv \
        [--percentile 0.99]

The mpileup must be in standard samtools format with the column layout
chrom, pos, ref, then per-pool triples of (depth, bases, quals). The
locations file is the one-per-line pool-name list matching mpileup column
order (same file used by the mpileup driver itself).

Performance
-----------
Streams the mpileup in a single pass. Wall-clock is I/O-bound, a few hours
for the 27-pool pileup.

Output format
-------------
Tab-separated, header-bearing. Example rows:

    pool        chrom                       max_depth_cutoff
    AGCB        ENA_OZ004311_OZ004311.1     87
    AGCB        ENA_OZ004312_OZ004312.1     91
    ...
    CAL1        ENA_OZ004313_OZ004313.1     142
"""

import argparse
import gzip
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from tdigest import TDigest
    HAVE_TDIGEST = True
except ImportError:
    HAVE_TDIGEST = False


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mpileup', required=True,
                   help='Path to v4 mpileup (.gz or plain text).')
    p.add_argument('--locations', required=True,
                   help='Pool-name list, one per line, in mpileup column order.')
    p.add_argument('--out', required=True,
                   help='Output TSV path.')
    p.add_argument('--percentile', type=float, default=0.99,
                   help='Percentile to compute (default 0.99, matching PoolSNP max-cov).')
    p.add_argument('--progress-every', type=int, default=10_000_000,
                   help='Log progress every N positions (default 10M).')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def open_mpileup(path):
    """Return a text-mode iterator over the mpileup, handling .gz transparently."""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def load_pool_names(locations_path):
    """One pool name per line, in mpileup column order."""
    names = []
    with open(locations_path, 'rt') as fh:
        for line in fh:
            name = line.strip()
            if name:
                names.append(name)
    if not names:
        raise ValueError(f'No pool names found in {locations_path}')
    return names


# ---------------------------------------------------------------------------
# Main streaming logic
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not HAVE_TDIGEST:
        sys.stderr.write(
            'ERROR: tdigest module not available. Install with:\n'
            '    pip install tdigest\n'
            'or activate the conda env that contains it.\n'
        )
        sys.exit(1)

    pool_names = load_pool_names(args.locations)
    n_pools = len(pool_names)
    sys.stderr.write(f'[INFO] loaded {n_pools} pool names from {args.locations}\n')

    # digests[(pool, chrom)] -> TDigest accumulating depth values for that pair.
    # defaultdict so new chromosomes appear automatically.
    digests = defaultdict(TDigest)

    # Per-chrom position counter (for the progress log).
    pos_per_chrom = defaultdict(int)

    t_start = time.time()
    n_lines = 0

    with open_mpileup(args.mpileup) as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            # Standard mpileup: chrom, pos, ref, then triples of (depth, bases, quals)
            # per pool. Total columns: 3 + 3 * n_pools.
            if len(fields) < 3 + 3 * n_pools:
                # Malformed line. Skip with a warning the first few times.
                if n_lines < 5:
                    sys.stderr.write(f'[WARN] short line ({len(fields)} cols): {line[:120]}\n')
                continue

            chrom = fields[0]

            # Per-pool depth lives in column 3 + 3*i (0-indexed: pool 0's depth is
            # field index 3, pool 1's depth is field index 6, etc.).
            for i, pool in enumerate(pool_names):
                depth_str = fields[3 + 3 * i]
                # samtools mpileup writes "0" (not empty) for zero-depth, but be
                # defensive about empty strings.
                if depth_str and depth_str != '0':
                    try:
                        depth = int(depth_str)
                        if depth > 0:
                            digests[(pool, chrom)].update(depth)
                    except ValueError:
                        # Unparseable depth field; rare but possible if the mpileup
                        # is corrupt. Don't let one bad line abort.
                        pass

            n_lines += 1
            pos_per_chrom[chrom] += 1

            if n_lines % args.progress_every == 0:
                elapsed = time.time() - t_start
                rate = n_lines / elapsed if elapsed > 0 else 0
                sys.stderr.write(
                    f'[INFO] processed {n_lines:>14,} positions in {elapsed:>6.0f}s '
                    f'({rate:>9,.0f} pos/s); current chrom={chrom}\n'
                )

    elapsed = time.time() - t_start
    sys.stderr.write(
        f'[INFO] done: {n_lines:,} positions across {len(pos_per_chrom)} chromosomes '
        f'in {elapsed:.0f}s\n'
    )

    # ---- Write the cutoff table -------------------------------------------
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, 'wt') as fh:
        fh.write('pool\tchrom\tmax_depth_cutoff\n')
        for (pool, chrom), td in sorted(digests.items()):
            cutoff = td.percentile(args.percentile * 100)
            # TDigest returns float; round to nearest int for use as an integer
            # depth threshold (consistent with PoolSNP's max-cov semantics).
            cutoff_int = int(round(cutoff))
            fh.write(f'{pool}\t{chrom}\t{cutoff_int}\n')

    n_rows = sum(1 for _ in digests)
    sys.stderr.write(f'[INFO] wrote {n_rows} (pool, chrom) cutoffs to {out_path}\n')


if __name__ == '__main__':
    main()
