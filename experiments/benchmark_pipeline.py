#!/usr/bin/env python3
"""Benchmark STROMBOLI pipeline performance on the real nanopore subsamples.

Measures, on config/test.yaml (2 real ~10k-read samples):
  1. End-to-end wall-clock + per-stage CPU-time breakdown, for both calling modes.
  2. Cores scaling (validates that threads:1 lets the per-barcode scatter parallelize).
  3. Scatter scaling: wall-clock vs number of barcodes (swept via min_cluster_size),
     fit to wall = fixed_cost + per_barcode_cost * n_barcodes.

Run inside the STROMBOLI conda env:
    python experiments/benchmark_pipeline.py
"""
import glob
import json
import os
import shutil
import subprocess
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "results")


def ncpu():
    try:
        return int(subprocess.run(["sysctl", "-n", "hw.ncpu"],
                                  capture_output=True, text=True).stdout.strip())
    except Exception:
        return os.cpu_count() or 4


def run(mode, cores, min_size=None, until=None):
    """Clear results, run the pipeline (optionally only --until a rule), return
    (wall_seconds, n_barcodes)."""
    shutil.rmtree(RESULTS, ignore_errors=True)
    cmd = ["snakemake", "-s", "workflow/Snakefile", "--configfile", "config/test.yaml",
           "--config", f"calling_mode={mode}"]
    if min_size is not None:
        cmd += [f"min_cluster_size={min_size}"]
    cmd += ["--cores", str(cores), "--quiet", "rules"]
    if until:
        cmd += ["--until", until]
    t = time.time()
    subprocess.run(cmd, cwd=REPO, check=True, capture_output=True, text=True)
    wall = time.time() - t
    n_bc = len(glob.glob(os.path.join(RESULTS, "clusters", "barcodes", "*", "*.fastq")))
    return wall, n_bc


def linfit(xs, ys):
    """Least-squares slope/intercept."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx if sxx else 0.0
    return slope, my - slope * mx


def main():
    cores = min(8, ncpu())
    print(f"machine cores = {ncpu()} | benchmarking with --cores {cores}\n")
    stats_dir = os.path.join(REPO, "experiments")

    # 1. End-to-end, both modes, full cores, plus a fixed-front vs scatter split.
    print("=" * 64)
    print("1. END-TO-END (real data, 2 samples, min_cluster_size=2)")
    print("=" * 64)
    mode_results = {}
    for mode in ("double", "single_qc"):
        wall, n_bc = run(mode, cores, min_size=2)
        mode_results[mode] = (wall, n_bc)
        print(f"  {mode:>10}: {wall:6.1f}s wall | {n_bc} barcodes scattered")
    d_wall = mode_results["double"][0]
    s_wall = mode_results["single_qc"][0]
    print(f"  single_qc speedup: {d_wall / s_wall:.2f}x")

    # Fixed front cost = cutadapt + starcode + clustering (up to the checkpoint).
    front, n_bc = run("single_qc", cores, min_size=2, until="make_cluster_fastas")
    print(f"\n  fixed front (cutadapt+starcode+cluster fastas): {front:.1f}s")
    print(f"  per-barcode scatter+aggregation (single_qc):    "
          f"{s_wall - front:.1f}s for {n_bc} barcodes")

    # 2. Cores scaling (single_qc).
    print("\n" + "=" * 64)
    print("2. CORES SCALING (single_qc, min_cluster_size=2)")
    print("=" * 64)
    core_levels = sorted(set([1, max(2, cores // 2), cores]))
    base = None
    print(f"  {'cores':>6} {'wall(s)':>8} {'speedup':>8}")
    for c in core_levels:
        wall, n_bc = run("single_qc", c, min_size=2)
        if base is None:
            base = wall
        print(f"  {c:>6} {wall:>8.1f} {base / wall:>7.2f}x")

    # 3. Scatter scaling: vary number of barcodes via min_cluster_size.
    print("\n" + "=" * 64)
    print("3. SCATTER SCALING (single_qc, full cores; barcodes via min_cluster_size)")
    print("=" * 64)
    xs, ys = [], []
    print(f"  {'min_size':>9} {'barcodes':>9} {'wall(s)':>8}")
    for ms in (50, 20, 10, 5, 2):
        wall, n_bc = run("single_qc", cores, min_size=ms)
        xs.append(n_bc)
        ys.append(wall)
        print(f"  {ms:>9} {n_bc:>9} {wall:>8.1f}")
    slope, intercept = linfit(xs, ys)
    print(f"\n  fit: wall ~= {intercept:.1f}s fixed + {1000 * slope:.2f} ms/barcode")
    print(f"  (fixed cost = cutadapt+starcode+aggregation; slope = per-barcode "
          f"scatter cost)")


if __name__ == "__main__":
    main()
