#!/usr/bin/env python3
"""Prototype: vectorized per-barcode variant calling via multi-sample mpileup.

STATUS: proven but NOT adopted (as of this exploration). Measured ~7x faster than
the fused per-barcode pipeline (single-threaded, ~190x fewer process spawns) and
recovered the synthetic truth exactly. Deferred because it collapses the per-barcode
evidence trail that keeps errors inspectable -- the FDR estimator, per-cluster
callable fraction, and WT/no-call logic all reason per barcode. Revisit only if
library scale (>~500k barcodes) demands it; would first need FP/recall-equivalence
validation on hard data and batching into sample-groups to bound the multi-sample
matrix. See the roadmap in README.md.


Instead of running the toolchain once per barcode (~7 process spawns x N barcodes),
this tags every read with a read-group whose SAMPLE is its barcode, maps everything
ONCE, and runs a single multi-sample `mpileup | call | norm | csq`. Each barcode is a
sample column; per-barcode calls fall out of one query with an ALT allele-fraction
filter. Total process spawns become O(1) in the number of barcodes.

Given a directory of per-barcode cluster FASTQs (produced by the pipeline's
make_cluster_fastas), it reports timing and -- if a truth CSV is given -- recovery.

    python experiments/prototype_multisample.py --fastq-dir results/clusters/barcodes/synthetic \
        --truth tests/synthetic_truth.csv

Run inside the STROMBOLI conda env.
"""
import argparse
import glob
import os
import subprocess
import tempfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(REPO, "references", "amplicon_ref.fasta")
GFF = os.path.join(REPO, "references", "amplicon_ref.gff")

# awk that injects an @RG line per barcode and tags each read RG:Z:<barcode>,
# where the barcode is the read name minus its trailing _<index>.
TAG_AWK = r"""
/^@/ { print; next }
!hdr { n=split(bcs,a,","); for(i=1;i<=n;i++) print "@RG\tID:"a[i]"\tSM:"a[i]; hdr=1 }
{ bc=$1; sub(/_[0-9]+$/,"",bc); print $0"\tRG:Z:"bc }
"""


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash",
                          capture_output=True, text=True).stdout


def build_tagged_bam(fastq_dir, workdir):
    """Concatenate all barcode FASTQs, map once, tag RG=barcode, sort -> bam."""
    fastqs = sorted(glob.glob(os.path.join(fastq_dir, "*.fastq")))
    barcodes = [os.path.splitext(os.path.basename(f))[0] for f in fastqs]
    combined = os.path.join(workdir, "all.fastq")
    with open(combined, "w") as out:
        for f in fastqs:
            with open(f) as h:
                out.write(h.read())
    n_reads = sum(1 for _ in open(combined)) // 4
    bcs = ",".join(barcodes)
    bam = os.path.join(workdir, "all.bam")
    awk = os.path.join(workdir, "tag.awk")
    with open(awk, "w") as h:
        h.write(TAG_AWK)
    t = time.time()
    sh(f"minimap2 -ax map-ont {REF} {combined} 2>/dev/null "
       f"| awk -v bcs='{bcs}' -f {awk} "
       f"| samtools sort -o {bam} 2>/dev/null")
    sh(f"samtools index {bam}")
    return bam, barcodes, n_reads, time.time() - t


def multisample_call(bam, workdir, af, alt_min):
    """One multi-sample mpileup|call|norm|csq, then per-sample AF-filtered query."""
    annotated = os.path.join(workdir, "anno.bcf")
    t = time.time()
    sh(f"bcftools mpileup -d 5000 -a AD -Ou -f {REF} {bam} 2>/dev/null "
       f"| bcftools call -mv --ploidy 1 -Ou 2>/dev/null "
       f"| bcftools norm -m- -f {REF} -Ou 2>/dev/null "
       f"| bcftools csq -f {REF} -g {GFF} -Ob -o {annotated} --verbose 0 - 2>/dev/null")
    query = sh_out(
        f"bcftools query -f "
        f"'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE\\t%GT\\t%AD]\\n' {annotated}")
    call_time = time.time() - t

    # Per-barcode variant set, applying the single_qc allele-fraction filter.
    per_barcode = {}
    for line in query.splitlines():
        f = line.split("\t")
        chrom, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
        rest = f[4:]
        for i in range(0, len(rest) - 2, 3):
            sample, gt, ad = rest[i], rest[i + 1], rest[i + 2]
            if gt not in ("1", "1/1", "1|1"):
                continue
            parts = ad.split(",")
            if len(parts) < 2:
                continue
            try:
                refd, altd = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            tot = refd + altd
            if altd >= alt_min and tot > 0 and altd / tot >= af:
                per_barcode.setdefault(sample, set()).add((pos, ref.upper(), alt.upper()))
    return per_barcode, call_time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq-dir", required=True)
    ap.add_argument("--truth")
    ap.add_argument("--af", type=float, default=0.85)
    ap.add_argument("--alt-min", type=int, default=2)
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as wd:
        bam, barcodes, n_reads, map_t = build_tagged_bam(args.fastq_dir, wd)
        per_barcode, call_t = multisample_call(bam, wd, args.af, args.alt_min)

    n_bc = len(barcodes)
    print(f"barcodes: {n_bc} | reads: {n_reads}")
    print(f"map+tag+sort : {map_t:6.2f}s   (1 minimap2 + 1 awk + 1 sort, all barcodes)")
    print(f"multisample  : {call_t:6.2f}s   (1 mpileup|call|norm|csq + 1 query)")
    print(f"TOTAL        : {map_t + call_t:6.2f}s   for {n_bc} barcodes  "
          f"({1000 * (map_t + call_t) / max(n_bc, 1):.2f} ms/barcode)")
    print(f"process spawns ~= 9 total (vs ~7 x {n_bc} = {7 * n_bc} per-barcode)")

    if args.truth:
        import csv
        truth = {r["barcode"]: (int(r["pos"]), r["alt"].upper())
                 for r in csv.DictReader(open(args.truth))}
        ok = 0
        for bc, (pos, alt) in truth.items():
            hit = any(v[0] == pos and v[2] == alt for v in per_barcode.get(bc, set()))
            ok += hit
            print(f"  {'OK ' if hit else 'MISS'} {bc}  expect {pos}{alt}  "
                  f"got {sorted(per_barcode.get(bc, set()))}")
        print(f"recovered {ok}/{len(truth)} truth SNPs")


if __name__ == "__main__":
    main()
