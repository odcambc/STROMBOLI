#!/usr/bin/env python3
"""Generate a small, deterministic synthetic dataset for STROMBOLI.

Each synthetic barcode is tied to exactly one known SNP in the gp17 ORF. For every
barcode we emit a number of full-amplicon nanopore-like reads in which:

  * the reference's 20xN barcode window is replaced by the barcode,
  * the assigned SNP is applied to the ORF,
  * low-rate substitution errors are sprinkled over the insert (but NOT over the
    flanking constant regions or the barcode, so cutadapt detection and starcode
    clustering stay exact and the ground truth is unambiguous).

Outputs:
  data/synthetic.fastq.gz        - reads to run the pipeline on
  tests/synthetic_truth.csv      - barcode,pos,ref,alt ground truth (1-based pos)

Run from the repository root:  python tests/generate_synthetic_data.py
"""
import csv
import gzip
import os
import random
import re

# Deterministic output.
SEED = 20240607
N_BARCODES = 6
READS_PER_BARCODE = 40
INSERT_ERROR_RATE = 0.015  # per-base substitution rate in the insert
HIGH_Q = 30  # Phred for confident bases
LOW_Q = 10  # Phred for error/low-confidence bases

# SNP positions (1-based, within the gp17 ORF 198-3237), spread along the ORF.
SNP_POSITIONS = [250, 500, 800, 1100, 1700, 2300]

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFERENCE = os.path.join(REPO_ROOT, "references", "amplicon_ref.fasta")
OUT_FASTQ = os.path.join(REPO_ROOT, "data", "synthetic.fastq.gz")
OUT_TRUTH = os.path.join(REPO_ROOT, "tests", "synthetic_truth.csv")

BASES = "ACGT"


def read_reference(path):
    """Return the single reference sequence (uppercased) from a FASTA file."""
    seq = []
    with open(path, "r", encoding="UTF-8") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def random_barcode(rng, length=20):
    return "".join(rng.choice(BASES) for _ in range(length))


def other_base(rng, base):
    return rng.choice([b for b in BASES if b != base])


def apply_insert_errors(rng, seq, protect_span):
    """Substitute bases at INSERT_ERROR_RATE, skipping the protected span
    (the flank+barcode+flank cassette). Returns (sequence, quality) strings."""
    start, end = protect_span
    out = []
    qual = []
    for i, base in enumerate(seq):
        protected = start <= i < end
        if not protected and rng.random() < INSERT_ERROR_RATE:
            out.append(other_base(rng, base))
            qual.append(chr(LOW_Q + 33))
        else:
            out.append(base)
            qual.append(chr(HIGH_Q + 33))
    return "".join(out), "".join(qual)


def main():
    rng = random.Random(SEED)
    reference = read_reference(REFERENCE)

    # Locate the barcode window (run of Ns) and the surrounding cassette so we can
    # protect flanks + barcode from synthetic errors.
    match = re.search(r"N+", reference)
    if not match:
        raise SystemExit("No N barcode window found in reference.")
    n_start, n_end = match.start(), match.end()
    # Protect a margin (the constant flanks) on either side of the barcode.
    protect_span = (n_start - 25, n_end + 25)

    # Build barcode -> SNP truth.
    truth = []
    for i in range(N_BARCODES):
        barcode = random_barcode(rng)
        pos = SNP_POSITIONS[i % len(SNP_POSITIONS)]
        ref_base = reference[pos - 1]
        alt_base = other_base(rng, ref_base)
        truth.append((barcode, pos, ref_base, alt_base))

    os.makedirs(os.path.dirname(OUT_FASTQ), exist_ok=True)
    os.makedirs(os.path.dirname(OUT_TRUTH), exist_ok=True)

    read_no = 0
    with gzip.open(OUT_FASTQ, "wt", encoding="UTF-8") as fq:
        for barcode, pos, ref_base, alt_base in truth:
            # Reference with this barcode inserted and the SNP applied.
            template = list(reference)
            template[n_start:n_end] = list(barcode)
            template[pos - 1] = alt_base
            template = "".join(template)
            for _ in range(READS_PER_BARCODE):
                seq, qual = apply_insert_errors(rng, template, protect_span)
                fq.write(
                    "@synthetic_read_{}_{}\n{}\n+\n{}\n".format(
                        read_no, barcode, seq, qual
                    )
                )
                read_no += 1

    with open(OUT_TRUTH, "w", encoding="UTF-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["barcode", "pos", "ref", "alt"])
        for barcode, pos, ref_base, alt_base in truth:
            writer.writerow([barcode, pos, ref_base, alt_base])

    print("Wrote {} reads for {} barcodes to {}".format(read_no, len(truth), OUT_FASTQ))
    print("Wrote ground truth to {}".format(OUT_TRUTH))


if __name__ == "__main__":
    main()
