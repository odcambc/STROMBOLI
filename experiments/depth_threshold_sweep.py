#!/usr/bin/env python3
"""Find a defensible minimum-depth cutoff for keeping a barcode cluster.

Policy under test: discard clusters below a depth threshold and call the rest with
SINGLE+QC (single-map + allele-fraction filter). This sweeps depth finely to show
(a) where the false-positive rate collapses and recall plateaus, and
(b) the purity/recall of the RETAINED set for each candidate cutoff.

Reuses the read generation + variant calling from run_mapping_comparison.py.
Run inside the STROMBOLI conda env:
    python experiments/depth_threshold_sweep.py
"""
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

from run_mapping_comparison import (
    make_jobs, call_cluster, QC_MIN_AF, QC_MIN_ALT_READS,
)

DEPTHS = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60]
N_PER_DEPTH = 30
CANDIDATE_CUTOFFS = [2, 4, 6, 8, 10, 12, 15, 20]


def mean(rows, method, key):
    return sum(r["scores"][method][key] for r in rows) / len(rows) if rows else 0.0


def pct(rows, method, key):
    return 100 * mean(rows, method, key)


def main():
    jobs, insert, hotspots = make_jobs(DEPTHS, N_PER_DEPTH)
    print(f"barcodes: {len(jobs)} | depths: {DEPTHS}")
    print(f"QC: AF>={QC_MIN_AF}, alt-reads>={QC_MIN_ALT_READS} | "
          f"calling survivors with SINGLE+QC\n")

    with ThreadPoolExecutor(max_workers=6) as ex:
        results = list(ex.map(call_cluster, jobs))

    by_depth = defaultdict(list)
    for r in results:
        by_depth[r["depth"]].append(r)

    # (a) Fine per-depth curve: where does the FP rate collapse / recall plateau?
    print("Per-depth (SINGLE+QC):              DOUBLE (reference):")
    print(f"{'depth':>5} {'n':>4} {'recall%':>7} {'clean%':>6} {'fpIND':>6} | "
          f"{'recall%':>7} {'clean%':>6} {'fpIND':>6}")
    print("-" * 64)
    for depth in DEPTHS:
        rows = by_depth[depth]
        print(f"{depth:>5} {len(rows):>4} "
              f"{pct(rows, 'single_qc', 'recall'):>7.0f} "
              f"{pct(rows, 'single_qc', 'exact'):>6.0f} "
              f"{mean(rows, 'single_qc', 'fp_indel'):>6.2f} | "
              f"{pct(rows, 'double', 'recall'):>7.0f} "
              f"{pct(rows, 'double', 'exact'):>6.0f} "
              f"{mean(rows, 'double', 'fp_indel'):>6.2f}")

    # (b) Retained-set purity/recall for each candidate min-depth cutoff.
    total = len(results)
    print("\nKeep clusters with depth >= D, call with SINGLE+QC:")
    print(f"{'min D':>6} {'kept':>6} {'kept%':>6} {'recall%':>7} {'clean%':>6} "
          f"{'fpSNP':>6} {'fpIND':>6}")
    print("-" * 52)
    for cutoff in CANDIDATE_CUTOFFS:
        kept = [r for r in results if r["depth"] >= cutoff]
        if not kept:
            continue
        print(f"{cutoff:>6} {len(kept):>6} {100*len(kept)/total:>5.0f}% "
              f"{pct(kept, 'single_qc', 'recall'):>7.0f} "
              f"{pct(kept, 'single_qc', 'exact'):>6.0f} "
              f"{mean(kept, 'single_qc', 'fp_snp'):>6.2f} "
              f"{mean(kept, 'single_qc', 'fp_indel'):>6.2f}")

    print("\nkept% assumes a uniform depth distribution across the swept depths; "
          "real\nyield depends on your library's actual cluster-size distribution.")
    print("recall% = true SNP recovered; clean% = exactly the true SNP and nothing "
          "else;\nfpSNP/fpIND = mean false-positive calls per retained cluster.")


if __name__ == "__main__":
    main()
