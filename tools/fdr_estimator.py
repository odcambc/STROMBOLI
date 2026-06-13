#!/usr/bin/env python3
"""Estimate per-cluster variant false-discovery and barcode-collision rates, to
guide STROMBOLI threshold selection (min_cluster_size, qc_min_af, barcode_distance).

The model is intentionally simple and analytic (no external deps). It has two parts.

1. PER-CLUSTER VARIANT FALSE POSITIVES
   A cluster of N reads over a callable region of length L is called with an
   allele-fraction rule: an ALT is kept if it is supported by >= `min_alt` reads
   AND by a fraction >= `tau` of the reads (tau = qc_min_af for single_qc, or the
   consensus call-fraction for double). A position yields a FALSE positive when, by
   chance, enough reads carry the same erroneous allele to cross that bar.

   Two error sources are modelled:
     * RANDOM substitutions at uniform per-base rate `eps` (each error goes to one
       of 3 bases). The count of a specific wrong base ~ Binomial(N, eps/3).
         E[FP_random per cluster] ~ L * 3 * P(Binom(N, eps/3) >= kmin)
     * SYSTEMATIC sites (e.g. homopolymers): `n_sys` positions where a correlated
       error appears in each read with probability `f_sys` (same direction across
       reads -- this is what ONT does and what naive pileups miscall).
         E[FP_systematic per cluster] ~ n_sys * P(Binom(N, f_sys) >= kmin)
     where kmin = max(min_alt, ceil(tau * N)).

   Key consequences the tool makes quantitative:
     - Random-substitution FDR is killed by a high `tau`, almost regardless of depth.
     - Systematic FDR is only controlled when tau > f_sys; depth then sharpens the
       fraction estimate. If tau <= f_sys, NO depth fixes it -- raise tau instead.
   This is exactly why the depth floor + AF threshold are the levers, and the tool
   recommends a minimum depth for a target FDR.

2. BARCODE COLLISIONS
   With M distinct barcodes of length B over a 4-letter alphabet and clustering
   distance d, two distinct barcodes can be merged when they fall within distance d.
     P(two random B-mers within Hamming d) = sum_{i=0..d} C(B,i) 3^i / 4^B
     E[colliding pairs] ~ C(M,2) * P_collide
   This guides barcode length / clustering distance / library size.

Usage:
    python tools/fdr_estimator.py                       # report with amplicon defaults
    python tools/fdr_estimator.py --error-rate 0.04 --tau 0.8 --target-fdr 0.02
    python tools/fdr_estimator.py --validate experiments/mapping_comparison_results.csv
"""
import argparse
import csv
import math


def binom_tail_ge(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1))


def kmin_for(n, tau, min_alt):
    return max(min_alt, math.ceil(tau * n))


def expected_fp_per_cluster(L, N, eps, tau, min_alt, n_sys=0, f_sys=0.0):
    """Expected number of false-positive variant calls in one cluster of depth N."""
    if N < min_alt:
        return 0.0  # can't make a call at all
    kmin = kmin_for(N, tau, min_alt)
    random_fp = L * 3 * binom_tail_ge(kmin, N, eps / 3)
    systematic_fp = n_sys * binom_tail_ge(kmin, N, f_sys) if n_sys else 0.0
    return random_fp + systematic_fp


def recommend_min_depth(L, eps, tau, min_alt, target_fp, n_sys=0, f_sys=0.0,
                        n_max=1000):
    """Smallest cluster depth N with E[FP per cluster] <= target_fp, or None."""
    for N in range(max(min_alt, 2), n_max + 1):
        if expected_fp_per_cluster(L, N, eps, tau, min_alt, n_sys, f_sys) <= target_fp:
            return N
    return None


def barcode_collision(B, M, d):
    """Return (p_pairwise_within_d, expected_colliding_pairs, p_any_collision)."""
    p = sum(math.comb(B, i) * 3**i for i in range(0, d + 1)) / 4.0**B
    pairs = math.comb(M, 2) * p
    p_any = 1 - math.exp(-pairs) if pairs < 700 else 1.0
    return p, pairs, p_any


def report(args):
    L, eps, tau, min_alt = args.length, args.error_rate, args.tau, args.min_alt
    n_sys, f_sys = args.n_systematic, args.f_systematic

    print("=" * 72)
    print("STROMBOLI FDR estimator")
    print("=" * 72)
    print(f"callable length L = {L} bp | per-base error eps = {eps} | "
          f"AF threshold tau = {tau} | min_alt = {min_alt}")
    print(f"systematic sites n_sys = {n_sys} at fraction f_sys = {f_sys}")
    print()

    print("Expected FALSE-POSITIVE variant calls per cluster vs depth:")
    print(f"{'depth':>6} {'random':>12} {'systematic':>12} {'total':>12}")
    for N in [2, 3, 5, 8, 10, 15, 20, 30, 50, 100, 200]:
        kmin = kmin_for(N, tau, min_alt)
        rand = L * 3 * binom_tail_ge(kmin, N, eps / 3)
        syst = n_sys * binom_tail_ge(kmin, N, f_sys) if n_sys else 0.0
        print(f"{N:>6} {rand:>12.3g} {syst:>12.3g} {rand + syst:>12.3g}")
    print()

    print("Recommended minimum cluster depth (min_cluster_size) for a target "
          "FP/cluster:")
    print(f"{'target FP/cluster':>18} {'min depth':>10}")
    for target in [0.1, 0.05, 0.02, 0.01, 0.005]:
        N = recommend_min_depth(L, eps, tau, min_alt, target, n_sys, f_sys)
        print(f"{target:>18} {('>1000' if N is None else N):>10}")
    if f_sys >= tau and n_sys:
        print(f"  NOTE: f_sys ({f_sys}) >= tau ({tau}): systematic errors cannot be "
              f"removed by depth.\n        Raise qc_min_af above {f_sys} to suppress "
              f"them.")
    print()

    p, pairs, p_any = barcode_collision(args.barcode_length, args.n_barcodes,
                                        args.barcode_distance)
    print("Barcode collisions:")
    print(f"  B = {args.barcode_length} | library M = {args.n_barcodes} | "
          f"clustering distance d = {args.barcode_distance}")
    print(f"  P(two random barcodes within distance d) = {p:.3g}")
    print(f"  expected colliding pairs                  = {pairs:.3g}")
    print(f"  P(>=1 collision in the library)           = {p_any:.3g}")
    per_read_err = 1 - (1 - eps) ** args.barcode_length
    print(f"  P(a barcode read carries >=1 error)       = {per_read_err:.3g}  "
          f"(why clustering is needed)")
    print("=" * 72)


def validate(path, args):
    """Compare model predictions to an empirical sweep (run_mapping_comparison.py)."""
    rows = list(csv.DictReader(open(path)))
    by_depth = {}
    for r in rows:
        d = int(r["depth"])
        by_depth.setdefault(d, []).append(r)
    print(f"Validation against {path} (SINGLE+QC empirical vs model):")
    print(f"  model params: L={args.length}, eps={args.error_rate}, tau={args.tau}, "
          f"n_sys={args.n_systematic}, f_sys={args.f_systematic}")
    print(f"{'depth':>6} {'n':>4} {'emp FP/clust':>13} {'model FP/clust':>15}")
    for depth in sorted(by_depth):
        rs = by_depth[depth]
        emp = sum(int(r["single_qc_fp_snp"]) + int(r["single_qc_fp_indel"])
                  for r in rs) / len(rs)
        model = expected_fp_per_cluster(
            args.length, depth, args.error_rate, args.tau, args.min_alt,
            args.n_systematic, args.f_systematic)
        print(f"{depth:>6} {len(rs):>4} {emp:>13.3g} {model:>15.3g}")
    print("\nThe uniform model captures the random/independent component; residual "
          "empirical\nFP at low depth is the systematic (homopolymer-indel) term -- "
          "tune n_sys/f_sys to it.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Defaults match the shipped amplicon (gp17 ORF 199-3237 = 3039 bp, 20 bp barcode).
    ap.add_argument("--length", type=int, default=3039, help="callable length (bp)")
    ap.add_argument("--error-rate", type=float, default=0.03,
                    help="uniform per-base substitution error rate")
    ap.add_argument("--tau", type=float, default=0.85,
                    help="allele-fraction calling threshold (qc_min_af)")
    ap.add_argument("--min-alt", type=int, default=2,
                    help="minimum ALT-supporting reads (qc_min_alt_reads)")
    ap.add_argument("--n-systematic", type=int, default=12,
                    help="number of systematic (e.g. homopolymer) error sites")
    ap.add_argument("--f-systematic", type=float, default=0.5,
                    help="error fraction at a systematic site (calibrate to data; "
                         "values approaching tau are the worst case and need depth)")
    ap.add_argument("--barcode-length", type=int, default=20)
    ap.add_argument("--n-barcodes", type=int, default=10000,
                    help="number of distinct barcodes in the library")
    ap.add_argument("--barcode-distance", type=int, default=5,
                    help="starcode clustering distance")
    ap.add_argument("--validate", metavar="CSV",
                    help="compare predictions to an empirical sweep CSV")
    args = ap.parse_args()

    if args.validate:
        validate(args.validate, args)
    else:
        report(args)


if __name__ == "__main__":
    main()
