#!/usr/bin/env python3
"""Evaluate WT calling under the converged design: single_qc + a cluster depth floor.

In single_qc mode the first-pass mpileup already sees every position, so affirmative
WT calling is a flag change: a position is WT only if it is covered at depth >= a
per-position floor AND matches the reference; below the floor it is a NO-CALL, not
reference. The alternative (Alt 2) is to IMPUTE -- if the cluster clears the
cluster-size floor, call every non-variant position WT regardless of local coverage.

These differ only when a cluster that clears the cluster floor still has per-position
coverage GAPS -- which happens when reads are truncated. This script sweeps the
fraction of full-length reads and measures, for variant clusters whose SNP sits in a
truncation-exposed (5') region:
  * detected            : variant correctly called
  * affirmative false-WT: SNP covered (depth>=floor) but missed -> wrongly WT
  * affirmative no-call : SNP under-covered -> correctly NOT called WT
  * impute   false-WT   : SNP not called for ANY reason -> imputed WT
The gap between impute false-WT and affirmative false-WT is the error affirmative
calling prevents. For WT clusters we also report callable fraction and spurious calls.

Run inside the STROMBOLI conda env:
    python experiments/evaluate_wt_calling.py
"""
import os
import random
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

from run_mapping_comparison import (
    read_reference, run_lengths, other_base, sh, sh_out, REFERENCE, FLANK1,
    P_SUB, P_HP_INDEL, HOMOPOLYMER_MIN, Q_GOOD, Q_BAD,
)

SEED = 77
DEPTH = 20                # reads per cluster (clears the cluster floor)
PER_POS_FLOOR = 5         # per-position depth needed to make an affirmative call
CLUSTER_FLOOR = 10        # min_cluster_size (all clusters here clear it)
P_FULL_SWEEP = [1.0, 0.5, 0.25, 0.1]   # fraction of full-length (untruncated) reads
TRUNC_START = (1300, 2200)             # 5'-truncated reads start in this window
N_VAR_5P = 40             # variant clusters with SNP in the truncation-exposed 5' region
N_VAR_3P = 20             # variant clusters with SNP in the always-covered 3' region
N_WT = 30                 # wild-type clusters
REGION_5P = (250, 1100)
REGION_3P = (2500, 3200)
QC_MIN_AF, QC_MIN_ALT = 0.85, 2

# Filter expression matching the single_qc pipeline rule.
CALL = (
    "bcftools mpileup -d 5000 -a AD -Ou -f {ref} {bam} 2>/dev/null "
    "| bcftools call -vm --ploidy 1 -Ov 2>/dev/null "
    "| bcftools norm -m- -f {ref} -Ov 2>/dev/null "
    "| bcftools view -i 'FMT/AD[0:1]>={ar} && "
    "FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])>={af}' -Ov 2>/dev/null"
)


def emit_read(template, start, rng):
    """Errors (subs + homopolymer indels) over template[start:]; returns seq, qual."""
    runlen = run_lengths(template)
    seq, qual = [], []
    i = start
    n = len(template)
    while i < n:
        base = template[i]
        if runlen[i] >= HOMOPOLYMER_MIN and rng.random() < P_HP_INDEL:
            i += 1
            continue
        if rng.random() < P_SUB:
            seq.append(other_base(rng, base)); qual.append(chr(Q_BAD + 33))
        else:
            seq.append(base); qual.append(chr(Q_GOOD + 33))
        if runlen[i] >= HOMOPOLYMER_MIN and rng.random() < P_HP_INDEL:
            seq.append(base); qual.append(chr(Q_BAD + 33))
        i += 1
    return "".join(seq), "".join(qual)


def parse_variant_positions(vcf_text):
    return {int(l.split("\t")[1]) for l in vcf_text.splitlines()
            if l and not l.startswith("#")}


def parse_depth(text):
    d = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) == 3:
            d[int(f[1])] = int(f[2])
    return d


def run_cluster(args):
    idx, kind, snp_pos, insert, p_full = args
    template = list(insert)
    snp_alt = None
    if snp_pos is not None:
        snp_alt = other_base(random.Random(f"alt-{idx}"), insert[snp_pos - 1])
        template[snp_pos - 1] = snp_alt
    template = "".join(template)
    rng = random.Random(f"{SEED}-{idx}-{p_full}")

    with tempfile.TemporaryDirectory() as d:
        fq = os.path.join(d, "r.fastq")
        with open(fq, "w") as f:
            for r in range(DEPTH):
                start = 0 if rng.random() < p_full else rng.randint(*TRUNC_START)
                seq, qual = emit_read(template, start, rng)
                f.write(f"@r{r}\n{seq}\n+\n{qual}\n")
        bam = os.path.join(d, "c.bam")
        sh(f"minimap2 -ax map-ont {REFERENCE} {fq} 2>/dev/null "
           f"| samtools sort -o {bam} 2>/dev/null")
        sh(f"samtools index {bam}")
        depth = parse_depth(sh_out(f"samtools depth -a {bam}"))
        variants = parse_variant_positions(
            sh_out(CALL.format(ref=REFERENCE, bam=bam, ar=QC_MIN_ALT, af=QC_MIN_AF)))

    return {"idx": idx, "kind": kind, "snp_pos": snp_pos, "p_full": p_full,
            "depth": depth, "variants": variants}


def main():
    rng = random.Random(SEED)
    insert = read_reference(REFERENCE)[:read_reference(REFERENCE).index(FLANK1)]
    orf = list(range(REGION_5P[0], REGION_3P[1] + 1))

    # Build the cluster spec once; reuse across truncation levels.
    specs = []
    idx = 0
    for _ in range(N_VAR_5P):
        specs.append((idx, "var5p", rng.randint(*REGION_5P))); idx += 1
    for _ in range(N_VAR_3P):
        specs.append((idx, "var3p", rng.randint(*REGION_3P))); idx += 1
    for _ in range(N_WT):
        specs.append((idx, "wt", None)); idx += 1

    print(f"insert {len(insert)} bp | depth {DEPTH} | per-pos floor {PER_POS_FLOOR} "
          f"| AF>={QC_MIN_AF} | clusters/level {len(specs)}")
    print(f"5' SNP region {REGION_5P} is exposed by 5' truncation; 3' region "
          f"{REGION_3P} stays covered.\n")

    print(f"{'p_full':>7} {'grp':>6} {'n':>4} | {'detect%':>7} "
          f"{'aff.FWT%':>8} {'aff.NC%':>7} | {'imp.FWT%':>8} | {'callable%':>9}")
    print("-" * 74)

    for p_full in P_FULL_SWEEP:
        jobs = [(i, k, s, insert, p_full) for (i, k, s) in specs]
        with ThreadPoolExecutor(max_workers=6) as ex:
            res = list(ex.map(run_cluster, jobs))

        groups = defaultdict(list)
        for r in res:
            groups[r["kind"]].append(r)

        # Variant groups: classify the SNP outcome.
        for grp in ("var5p", "var3p"):
            rows = groups[grp]
            n = len(rows)
            detect = aff_fwt = aff_nc = 0
            for r in rows:
                pos = r["snp_pos"]
                covered = r["depth"].get(pos, 0) >= PER_POS_FLOOR
                called = pos in r["variants"]
                if called:
                    detect += 1
                elif covered:
                    aff_fwt += 1            # covered but missed -> affirmative false-WT
                else:
                    aff_nc += 1             # under-covered -> affirmative no-call
            imp_fwt = aff_fwt + aff_nc      # impute calls WT whenever not detected
            callable_frac = 100 * sum(
                sum(1 for p in orf if r["depth"].get(p, 0) >= PER_POS_FLOOR) / len(orf)
                for r in rows) / n
            print(f"{p_full:>7} {grp:>6} {n:>4} | {100*detect/n:>7.0f} "
                  f"{100*aff_fwt/n:>8.0f} {100*aff_nc/n:>7.0f} | "
                  f"{100*imp_fwt/n:>8.0f} | {callable_frac:>9.0f}")

        # WT group: spurious variant calls and callable fraction.
        wt = groups["wt"]
        n = len(wt)
        spurious = sum(1 for r in wt if r["variants"])
        callable_frac = 100 * sum(
            sum(1 for p in orf if r["depth"].get(p, 0) >= PER_POS_FLOOR) / len(orf)
            for r in wt) / n
        fully = sum(1 for r in wt
                    if all(r["depth"].get(p, 0) >= PER_POS_FLOOR for p in orf)
                    and not r["variants"])
        print(f"{p_full:>7} {'wt':>6} {n:>4} | {'-':>7} {'-':>8} {'-':>7} | "
              f"{'spur=' + str(spurious):>8} | {callable_frac:>9.0f}  "
              f"(fully-confirmed WT: {100*fully/n:.0f}%)")
        print("-" * 74)

    print("\ndetect% = variant called; aff.FWT% = affirmative false-WT (SNP covered "
          "but missed);\naff.NC% = affirmative no-call (SNP under-covered, correctly "
          "not WT);\nimp.FWT% = impute false-WT (= aff.FWT + aff.NC). callable% = mean "
          "fraction of\nORF positions at depth >= floor. WT 'spur' = clusters with a "
          "spurious variant call.")


if __name__ == "__main__":
    main()
