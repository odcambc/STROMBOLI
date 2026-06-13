#!/usr/bin/env python3
"""Single-map vs double-map vs single-map+depth-QC on synthetic ONT clusters.

STROMBOLI calls variants per barcode cluster via
    reads -> consensus FASTA -> RE-MAP consensus -> mpileup        (DOUBLE map)
We compare three strategies, holding everything else constant:
  * SINGLE      : reads -> mpileup directly on the cluster pileup
  * DOUBLE      : the current pipeline (consensus then re-map)
  * SINGLE+QC   : SINGLE, then filter each call by ALT allele-fraction and minimum
                  supporting reads -- i.e. use the depth information SINGLE keeps
                  (and DOUBLE discards) to reject low-frequency calls.

For each synthetic barcode cluster we generate insert reads carrying a known SNP
plus a realistic ONT error model (random substitutions, homopolymer-length random
indels, and a few SYSTEMATIC read-correlated homopolymer indel hotspots), run the
three strategies with the SAME flags STROMBOLI uses, normalize the VCFs, and score
against truth (the single SNP is the only correct call; anything else is a false
positive).

Run inside the STROMBOLI conda env:
    python experiments/run_mapping_comparison.py
"""
import csv
import os
import random
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

# ----------------------------- configuration -----------------------------------
SEED = 1234
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE = os.path.join(REPO_ROOT, "references", "amplicon_ref.fasta")
RESULTS_CSV = os.path.join(REPO_ROOT, "experiments", "mapping_comparison_results.csv")
FLANK1 = "GCAGTCTGGTGTATGCCTAC"  # first constant flank; insert is everything before it

DEPTHS = [3, 5, 8, 15, 30, 60]
N_PER_DEPTH = 18  # -> 108 barcodes total

# Error model.
P_SUB = 0.02
HOMOPOLYMER_MIN = 3
P_HP_INDEL = 0.06
N_SYSTEMATIC_HOTSPOTS = 12
Q_GOOD, Q_BAD = 24, 8

# Depth-QC thresholds for SINGLE+QC. A true haploid SNP sits in ~all reads (AF~1.0);
# the systematic indel hotspots are injected at <=0.8 frequency, so an AF floor
# above that should remove them while keeping the real variant.
QC_MIN_AF = 0.85
QC_MIN_ALT_READS = 2

BASES = "ACGT"
METHODS = ["single", "double", "single_qc"]

# Pipeline flags from workflow/rules/all.smk, plus -a AD so we get per-allele depth.
MPILEUP_CALL = (
    "bcftools mpileup -d 5000 -a AD -Ou -f {ref} {bam} 2>/dev/null "
    "| bcftools call -vm --ploidy 1 -Ov 2>/dev/null "
    "| bcftools norm -m- -f {ref} -Ov 2>/dev/null"
)


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(
        cmd, shell=True, check=True, executable="/bin/bash",
        capture_output=True, text=True,
    ).stdout


def read_reference(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def run_lengths(seq):
    n = len(seq)
    out = [1] * n
    i = 0
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        for k in range(i, j):
            out[k] = j - i
        i = j
    return out


def other_base(rng, base):
    return rng.choice([b for b in BASES if b != base])


def generate_read(template, runlen, hotspots, rng):
    seq, qual = [], []
    i, n = 0, len(template)
    while i < n:
        base = template[i]
        if i in hotspots:
            op, frac = hotspots[i]
            if rng.random() < frac:
                if op == "del":
                    i += 1
                    continue
                seq.append(base)
                qual.append(chr(Q_BAD + 33))
        if runlen[i] >= HOMOPOLYMER_MIN and rng.random() < P_HP_INDEL:
            i += 1
            continue
        if rng.random() < P_SUB:
            seq.append(other_base(rng, base))
            qual.append(chr(Q_BAD + 33))
        else:
            seq.append(base)
            qual.append(chr(Q_GOOD + 33))
        if runlen[i] >= HOMOPOLYMER_MIN and rng.random() < P_HP_INDEL:
            seq.append(base)
            qual.append(chr(Q_BAD + 33))
        i += 1
    return "".join(seq), "".join(qual)


def parse_records(text):
    """Return list of (pos, ref, alt, af, alt_reads) from a normalized VCF.

    af / alt_reads come from FORMAT/AD ([ref_depth, alt_depth]); when AD is absent
    we pass the call through unfiltered (af=1.0) so QC never silently drops a call
    for a parsing reason."""
    records = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        pos, ref, alt = int(f[1]), f[3].upper(), f[4].upper()
        af, alt_reads = 1.0, 10 ** 6
        if len(f) >= 10:
            fmt = f[8].split(":")
            sample = f[9].split(":")
            if "AD" in fmt:
                ad = sample[fmt.index("AD")].split(",")
                try:
                    vals = [int(x) for x in ad]
                    total = sum(vals)
                    if total > 0 and len(vals) >= 2:
                        alt_reads = vals[1]
                        af = vals[1] / total
                except ValueError:
                    pass
        records.append((pos, ref, alt, af, alt_reads))
    return records


def is_indel(v):
    return len(v[1]) != len(v[2])


def call_cluster(args):
    idx, barcode, depth, snp, near_hotspot, insert, runlen, hotspots = args
    truth = snp  # (pos, ref, alt)

    template = list(insert)
    template[snp[0] - 1] = snp[2]
    template = "".join(template)
    tmpl_runlen = run_lengths(template)
    rng = random.Random(f"{SEED}-{idx}")

    with tempfile.TemporaryDirectory() as d:
        reads_fq = os.path.join(d, "reads.fastq")
        with open(reads_fq, "w") as fq:
            for r in range(depth):
                seq, qual = generate_read(template, tmpl_runlen, hotspots, rng)
                fq.write(f"@{barcode}_{r}\n{seq}\n+\n{qual}\n")

        cluster_bam = os.path.join(d, "cluster.bam")
        sh(f"minimap2 -ax map-ont {REFERENCE} {reads_fq} 2>/dev/null "
           f"| samtools sort -o {cluster_bam} 2>/dev/null")
        sh(f"samtools index {cluster_bam}")

        single_recs = parse_records(
            sh_out(MPILEUP_CALL.format(ref=REFERENCE, bam=cluster_bam)))

        cons_fa = os.path.join(d, "cons.fasta")
        sh(f"samtools consensus -q -c 0.75 -m simple -f fasta {cluster_bam} "
           f"-o {cons_fa} 2>/dev/null")
        cons_bam = os.path.join(d, "cons.bam")
        sh(f"minimap2 -ax map-ont {REFERENCE} {cons_fa} 2>/dev/null "
           f"| samtools sort -o {cons_bam} 2>/dev/null")
        sh(f"samtools index {cons_bam}")
        double_recs = parse_records(
            sh_out(MPILEUP_CALL.format(ref=REFERENCE, bam=cons_bam)))

    def to_set(recs):
        return {(p, r, a) for (p, r, a, af, ar) in recs}

    called = {
        "single": to_set(single_recs),
        "double": to_set(double_recs),
        "single_qc": {
            (p, r, a) for (p, r, a, af, ar) in single_recs
            if af >= QC_MIN_AF and ar >= QC_MIN_ALT_READS
        },
    }

    def score(s):
        fp = s - {truth}
        return {
            "recall": int(truth in s),
            "fp_snp": sum(1 for v in fp if not is_indel(v)),
            "fp_indel": sum(1 for v in fp if is_indel(v)),
            "exact": int(s == {truth}),
        }

    return {"barcode": barcode, "depth": depth, "near_hotspot": near_hotspot,
            "scores": {m: score(called[m]) for m in METHODS}}


def build_insert_and_hotspots():
    """Return (insert, runlen, hotspots, orf_lo, orf_hi) with a fixed seed so the
    error landscape is identical across experiments that import this module."""
    rng = random.Random(SEED)
    reference = read_reference(REFERENCE)
    insert = reference[:reference.index(FLANK1)]
    runlen = run_lengths(insert)
    orf_lo, orf_hi = 220, min(3200, len(insert) - 5)
    hp_sites = [i for i in range(orf_lo, orf_hi) if runlen[i] >= 4]
    rng.shuffle(hp_sites)
    hotspots = {}
    for pos in hp_sites[:N_SYSTEMATIC_HOTSPOTS]:
        hotspots[pos] = (rng.choice(["del", "ins"]),
                         rng.choice([0.3, 0.4, 0.5, 0.6, 0.7, 0.8]))
    return insert, runlen, hotspots, orf_lo, orf_hi


def make_jobs(depths, n_per_depth):
    """Build cluster jobs (one per barcode) for the given depth schedule."""
    insert, runlen, hotspots, orf_lo, orf_hi = build_insert_and_hotspots()
    hotspot_positions = sorted(hotspots)
    rng = random.Random(f"jobs-{SEED}")
    jobs, idx = [], 0
    for depth in depths:
        for _ in range(n_per_depth):
            barcode = "".join(rng.choice(BASES) for _ in range(20))
            near = (idx % 3 == 0) and bool(hotspot_positions)
            if near:
                hp = hotspot_positions[idx % len(hotspot_positions)]
                pos = max(orf_lo, min(orf_hi, hp + rng.choice([-3, -2, 2, 3])))
            else:
                pos = rng.randint(orf_lo, orf_hi)
            ref_base = insert[pos - 1]
            jobs.append((idx, barcode, depth,
                         (pos, ref_base, other_base(rng, ref_base)),
                         near, insert, runlen, hotspots))
            idx += 1
    return jobs, insert, hotspots


def main():
    jobs, insert, hotspots = make_jobs(DEPTHS, N_PER_DEPTH)

    print(f"Insert length: {len(insert)} bp | barcodes: {len(jobs)} | "
          f"systematic hotspots: {len(hotspots)} | "
          f"QC: AF>={QC_MIN_AF}, alt-reads>={QC_MIN_ALT_READS}")
    print("Calling variants three ways per cluster ...")

    with ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(call_cluster, jobs))

    # Per-cluster CSV.
    with open(RESULTS_CSV, "w", newline="") as f:
        w = csv.writer(f)
        header = ["barcode", "depth", "near_hotspot"]
        for m in METHODS:
            header += [f"{m}_recall", f"{m}_fp_snp", f"{m}_fp_indel", f"{m}_exact"]
        w.writerow(header)
        for r in results:
            row = [r["barcode"], r["depth"], int(r["near_hotspot"])]
            for m in METHODS:
                s = r["scores"][m]
                row += [s["recall"], s["fp_snp"], s["fp_indel"], s["exact"]]
            w.writerow(row)

    # Aggregate by (depth, method).
    agg = defaultdict(lambda: defaultdict(float))
    counts = defaultdict(int)
    for r in results:
        counts[r["depth"]] += 1
        for m in METHODS:
            for k, v in r["scores"][m].items():
                agg[(r["depth"], m)][k] += v

    labels = {"single": "SINGLE", "double": "DOUBLE", "single_qc": "SINGLE+QC"}
    print()
    header = f"{'depth':>6} {'n':>4} {'method':>10} | {'recall%':>7} {'clean%':>6} {'fpSNP':>6} {'fpIND':>6}"
    print(header)
    print("-" * len(header))
    totals = defaultdict(lambda: defaultdict(float))
    total_n = 0
    for depth in DEPTHS:
        n = counts[depth]
        total_n += n
        for m in METHODS:
            a = agg[(depth, m)]
            print(f"{depth:>6} {n:>4} {labels[m]:>10} | "
                  f"{100*a['recall']/n:7.0f} {100*a['exact']/n:6.0f} "
                  f"{a['fp_snp']/n:6.2f} {a['fp_indel']/n:6.2f}")
            for k, v in a.items():
                totals[m][k] += v
        print("-" * len(header))
    for m in METHODS:
        a = totals[m]
        print(f"{'ALL':>6} {total_n:>4} {labels[m]:>10} | "
              f"{100*a['recall']/total_n:7.0f} {100*a['exact']/total_n:6.0f} "
              f"{a['fp_snp']/total_n:6.2f} {a['fp_indel']/total_n:6.2f}")

    print("\nSNP recall by proximity to a systematic indel hotspot:")
    print(f"{'subset':>22} {'n':>4} " + " ".join(f"{labels[m]:>10}" for m in METHODS))
    for label, want in (("near hotspot (<=3bp)", True), ("away from hotspot", False)):
        subset = [r for r in results if r["near_hotspot"] == want]
        n = len(subset)
        cells = " ".join(
            f"{100*sum(r['scores'][m]['recall'] for r in subset)/n:9.0f}%"
            for m in METHODS)
        print(f"{label:>22} {n:>4} {cells}")

    print("\nrecall% = true SNP recovered; clean% = called EXACTLY the true SNP and "
          "nothing else;\nfpSNP/fpIND = mean false-positive SNP / indel calls per "
          "cluster.")
    print(f"Per-cluster results written to {RESULTS_CSV}")


if __name__ == "__main__":
    main()
