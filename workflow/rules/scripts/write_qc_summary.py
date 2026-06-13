"""Aggregate a sample's pipeline outputs into one QC summary JSON.

Writes results/qc/{sample}.stromboli_qc.json, consumed by the multiqc-stromboli plugin.
SCHEMA_VERSION must match the plugin's; bump both together when the format changes.
"""
import csv
import json
import statistics

SCHEMA_VERSION = 1

CLUSTER_BINS = ["1", "2", "3-5", "6-10", "11+"]
AF_BINS = ["0.2-0.4", "0.4-0.6", "0.6-0.85", "0.85-1.0"]


def _bin_cluster_size(n):
    if n <= 1:
        return "1"
    if n == 2:
        return "2"
    if n <= 5:
        return "3-5"
    if n <= 10:
        return "6-10"
    return "11+"


def _bin_af(af):
    if af < 0.4:
        return "0.2-0.4"
    if af < 0.6:
        return "0.4-0.6"
    if af < 0.85:
        return "0.6-0.85"
    return "0.85-1.0"


def _read_cutadapt_counts(path):
    """(reads_total, reads_with_barcode) from a cutadapt --json report."""
    with open(path, encoding="UTF-8") as f:
        rc = json.load(f).get("read_counts", {})
    # With --discard-untrimmed, 'output' is the reads that kept a detected barcode.
    return rc.get("input"), rc.get("output")


def build_qc_summary(
    sample, cutadapt_json, cluster_qc, consensus_variants, variants, flagged,
    min_cluster_size,
):
    reads_total, reads_with_barcode = _read_cutadapt_counts(cutadapt_json)

    sizes = []
    with open(cluster_qc, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sizes.append(int(row["n_reads"]))
    passing = [s for s in sizes if s >= min_cluster_size]
    cluster_hist = {b: 0 for b in CLUSTER_BINS}
    for s in sizes:
        cluster_hist[_bin_cluster_size(s)] += 1

    af_hist = {b: 0 for b in AF_BINS}
    with open(consensus_variants, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            v = row.get("af", "")
            if v not in ("", "None", None):
                try:
                    af_hist[_bin_af(float(v))] += 1
                except ValueError:
                    pass

    with open(variants, encoding="UTF-8") as f:
        n_variants = sum(1 for _ in csv.DictReader(f, delimiter="\t"))

    n_merged = n_mixed = 0
    with open(flagged, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            reason = row.get("reason", "")
            n_merged += "merged" in reason
            n_mixed += "mixed" in reason

    return {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "reads_total": reads_total,
        "reads_with_barcode": reads_with_barcode,
        "n_clusters": len(sizes),
        "n_clusters_passing": len(passing),
        "median_cluster_size": int(statistics.median(passing)) if passing else 0,
        "n_variants": n_variants,
        "n_flagged_merged": n_merged,
        "n_flagged_mixed": n_mixed,
        "cluster_size_histogram": cluster_hist,
        "allele_fraction_histogram": af_hist,
    }


if "snakemake" in globals():
    summary = build_qc_summary(
        snakemake.wildcards["sample"],
        snakemake.input["cutadapt"],
        snakemake.input["cluster_qc"],
        snakemake.input["consensus_variants"],
        snakemake.input["variants"],
        snakemake.input["flagged"],
        snakemake.params["min_cluster_size"],
    )
    with open(snakemake.output[0], "w", encoding="UTF-8") as f:
        json.dump(summary, f, indent=2)
