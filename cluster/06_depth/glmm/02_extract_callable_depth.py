#!/usr/bin/env python3
"""
cluster/06_depth/glmm/02_extract_callable_depth.py
==================================================

Step 2 of the callable-depth pipeline (cluster/README.md step 8; Methods
"Callable-position read-depth analysis"). Inputs: pools27_v4.sanitized.mpileup.gz,
pool_chrom_cutoffs_v4.tsv (step 1). Output: callable_depth_v4.tsv.gz.

Streams through the v4 mpileup a second time. At each position, applies
the PoolSNP-equivalent position-level filters:

  - min-cov 10:   every pool has depth >= 10
  - max-cov 0.99: no pool's depth exceeds its (pool, chrom) 99th percentile
                  cutoff from step 1

If both pass, writes the chromosome, position, and the per-pool depths for
all pools to the output. Positions failing either filter are skipped.

The output is the callable position x pool depth matrix. Subsequent steps
intersect it with the CDS annotation (step 3) and aggregate to gene x pool
depth sums (step 4).

Usage
-----
    python 02_extract_callable_depth.py \
        --mpileup pools27_v4.sanitized.mpileup.gz \
        --locations locations_v4.txt \
        --cutoffs pool_chrom_cutoffs_v4.tsv \
        --min-cov 10 \
        --out callable_depth_v4.tsv.gz

Output format
-------------
Tab-separated, gzipped, header-bearing. Columns: chrom, pos, then one column
per pool in mpileup column order. Each pool column holds the integer depth
at that callable position. Example:

    chrom                       pos       AGCB    Beltsville  ...  CAL1
    ENA_OZ004311_OZ004311.1     1234      18      22          ...  31
    ENA_OZ004311_OZ004311.1     1235      19      23          ...  30
    ...

Validation
----------
Stats are written to stderr: total positions seen, n passing each filter,
n passing both. As a sanity check, the SNP set in the v4 VCF should be a
subset of the callable set written here. Spot-check by intersecting the
output positions on chr1 with VCF positions on chr1; every VCF position
should be present.
"""

import argparse
import gzip
import sys
import time
from pathlib import Path


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
    p.add_argument('--cutoffs', required=True,
                   help='Output of step 1: TSV with columns pool, chrom, max_depth_cutoff.')
    p.add_argument('--min-cov', type=int, default=10,
                   help='Minimum depth required at every pool (default 10, matching PoolSNP min-cov).')
    p.add_argument('--out', required=True,
                   help='Output gzipped TSV path.')
    p.add_argument('--progress-every', type=int, default=10_000_000,
                   help='Log progress every N positions (default 10M).')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def open_mpileup(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def load_pool_names(locations_path):
    names = []
    with open(locations_path, 'rt') as fh:
        for line in fh:
            name = line.strip()
            if name:
                names.append(name)
    if not names:
        raise ValueError(f'No pool names found in {locations_path}')
    return names


def load_cutoffs(cutoffs_path, pool_names):
    """
    Build a lookup dict: cutoffs[(pool, chrom)] -> int max_depth_cutoff.

    The cutoffs TSV from step 1 may include (pool, chrom) combinations the
    user filters out later (e.g., unplaced contigs); we keep them all. A
    KeyError at lookup time means we encountered an mpileup position whose
    (pool, chrom) wasn't represented in step 1 — that's a pipeline error
    worth surfacing loudly, not silently skipping.
    """
    cutoffs = {}
    with open(cutoffs_path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header != ['pool', 'chrom', 'max_depth_cutoff']:
            raise ValueError(
                f'Unexpected cutoffs header in {cutoffs_path}: {header}'
            )
        for line in fh:
            pool, chrom, cutoff_str = line.rstrip('\n').split('\t')
            cutoffs[(pool, chrom)] = int(cutoff_str)

    # Sanity check: every pool in the locations file must appear in the cutoff
    # table (at least once, on some chrom).
    pools_in_cutoffs = {p for p, _ in cutoffs.keys()}
    missing = set(pool_names) - pools_in_cutoffs
    if missing:
        raise ValueError(
            f'Pools in locations file but missing from cutoffs table: {sorted(missing)}'
        )
    return cutoffs


# ---------------------------------------------------------------------------
# Main streaming logic
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    pool_names = load_pool_names(args.locations)
    n_pools = len(pool_names)
    sys.stderr.write(f'[INFO] loaded {n_pools} pool names from {args.locations}\n')

    cutoffs = load_cutoffs(args.cutoffs, pool_names)
    sys.stderr.write(
        f'[INFO] loaded {len(cutoffs)} (pool, chrom) cutoffs from {args.cutoffs}\n'
    )

    # Counters for the validation summary at the end.
    n_positions = 0
    n_pass_min = 0
    n_pass_both = 0
    n_callable_per_chrom = {}

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    t_start = time.time()

    # Pre-fetch the column indices for each pool's depth field for speed.
    # Pool i's depth lives at field index 3 + 3*i.
    depth_indices = [3 + 3 * i for i in range(n_pools)]

    with open_mpileup(args.mpileup) as fh_in, \
         gzip.open(out_path, 'wt') as fh_out:

        # Header
        fh_out.write('chrom\tpos\t' + '\t'.join(pool_names) + '\n')

        # Cache the last chrom to avoid repeated dict lookups for the cutoff
        # array. When chrom changes, refresh the per-pool cutoff vector.
        current_chrom = None
        chrom_cutoffs = None  # list of int, indexed same as pool_names

        for line in fh_in:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3 + 3 * n_pools:
                continue  # short line; already warned in step 1

            chrom = fields[0]
            pos = fields[1]

            n_positions += 1

            if chrom != current_chrom:
                current_chrom = chrom
                try:
                    chrom_cutoffs = [cutoffs[(p, chrom)] for p in pool_names]
                except KeyError as e:
                    sys.stderr.write(
                        f'[ERROR] cutoff missing for {e.args[0]} on chrom {chrom}. '
                        f'This means step 1 did not produce a cutoff for that '
                        f'(pool, chrom) pair. Re-run step 1.\n'
                    )
                    sys.exit(1)

            # Parse the per-pool depths.
            try:
                depths = [int(fields[idx]) for idx in depth_indices]
            except ValueError:
                # Unparseable depth field; skip this position.
                continue

            # Filter 1: min-cov. Every pool must have depth >= min_cov.
            if min(depths) < args.min_cov:
                continue
            n_pass_min += 1

            # Filter 2: max-cov. No pool's depth may exceed its (pool, chrom) cutoff.
            over_cap = False
            for d, cap in zip(depths, chrom_cutoffs):
                if d > cap:
                    over_cap = True
                    break
            if over_cap:
                continue

            n_pass_both += 1
            n_callable_per_chrom[chrom] = n_callable_per_chrom.get(chrom, 0) + 1

            # Write the callable row.
            fh_out.write(f'{chrom}\t{pos}\t' + '\t'.join(str(d) for d in depths) + '\n')

            if n_positions % args.progress_every == 0:
                elapsed = time.time() - t_start
                rate = n_positions / elapsed if elapsed > 0 else 0
                sys.stderr.write(
                    f'[INFO] {n_positions:>14,} positions seen, '
                    f'{n_pass_both:>14,} callable ({100*n_pass_both/n_positions:>5.1f}%); '
                    f'{rate:>9,.0f} pos/s; current chrom={chrom}\n'
                )

    elapsed = time.time() - t_start

    # ---- Final summary -----------------------------------------------------
    sys.stderr.write('\n')
    sys.stderr.write('=' * 72 + '\n')
    sys.stderr.write(f'[INFO] total positions seen:           {n_positions:>14,}\n')
    sys.stderr.write(f'[INFO] passed min-cov ({args.min_cov:>3}):              {n_pass_min:>14,}'
                     f'  ({100*n_pass_min/max(n_positions,1):>5.1f}%)\n')
    sys.stderr.write(f'[INFO] passed min-cov AND max-cov:    {n_pass_both:>14,}'
                     f'  ({100*n_pass_both/max(n_positions,1):>5.1f}%)\n')
    sys.stderr.write('[INFO] callable positions per chrom:\n')
    for chrom in sorted(n_callable_per_chrom):
        sys.stderr.write(f'         {chrom}: {n_callable_per_chrom[chrom]:>14,}\n')
    sys.stderr.write(f'[INFO] elapsed: {elapsed:.0f}s\n')
    sys.stderr.write(f'[INFO] output written to {out_path}\n')


if __name__ == '__main__':
    main()
