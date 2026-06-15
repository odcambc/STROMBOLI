"""Aggregate a sample's pipeline outputs into one QC summary JSON.

Writes results/qc/{sample}.stromboli_qc.json, consumed by the multiqc-stromboli plugin.
SCHEMA_VERSION must match the plugin's; bump both together when the format changes.
"""
import csv
import json
import math
import re
import statistics
from collections import Counter

# v4: barcode/variant coverage analytics, all from existing tables — ORF positional
# coverage + evenness (Gini) from amino_acid_change codon positions (needs the orf length,
# passed as a param), variant complexity (SNV/MNV/insertion/deletion) from INDEL + the '+'
# in dna_change, barcodes-per-variant redundancy from variants.tsv, and barcode composition
# (length, GC) from the cluster_qc barcode column.
# v3: variant-consequence classes, variants-per-barcode spread, cluster purity, scalars.
# v2: exact per-size cluster counts; the plugin chooses bins across the whole cohort.
SCHEMA_VERSION = 4

AF_BINS = ["0.2-0.4", "0.4-0.6", "0.6-0.85", "0.85-1.0"]
# second_member_fraction is a purity measure in [0,1]: 0 means a clean single-member
# cluster, higher means a competing barcode shares the cluster (merge pressure). Bounded,
# so pre-binned like AF. "0" is its own bin because clean clusters dominate a good library.
PURITY_BINS = ["0", "0-0.1", "0.1-0.25", "0.25-0.5", "0.5+"]
GC_BINS = ["<0.3", "0.3-0.45", "0.45-0.55", "0.55-0.7", ">0.7"]
# The ORF coverage profile is down-sampled to this many windows so the JSON stays small
# regardless of ORF length; the evenness Gini is still computed at full per-codon resolution.
N_COVERAGE_WINDOWS = 50

_LEADING_INT = re.compile(r"\s*(\d+)")


def _bin_af(af):
    if af < 0.4:
        return "0.2-0.4"
    if af < 0.6:
        return "0.4-0.6"
    if af < 0.85:
        return "0.6-0.85"
    return "0.85-1.0"


def _bin_purity(p):
    if p <= 0:
        return "0"
    if p < 0.1:
        return "0-0.1"
    if p < 0.25:
        return "0.1-0.25"
    if p < 0.5:
        return "0.25-0.5"
    return "0.5+"


def _bin_gc(gc):
    if gc < 0.3:
        return "<0.3"
    if gc < 0.45:
        return "0.3-0.45"
    if gc < 0.55:
        return "0.45-0.55"
    if gc < 0.7:
        return "0.55-0.7"
    return ">0.7"


def _consequence_class(raw):
    """Normalize a bcftools csq 'consequence' field to a coarse class.

    csq prefixes compound/duplicate predictions with '*'/'@' and may join several
    consequences with '&' (first is most severe); a blank field means the variant fell
    outside the annotated transcript (barcode cassette / flanks) -> 'noncoding'.
    """
    c = (raw or "").strip().lstrip("*@")
    if not c:
        return "noncoding"
    return c.split("&")[0]


def _variant_type(indel, ref, alt, dna_change):
    """Coarse structural class of a variant record.

    INDEL set -> insertion/deletion by REF/ALT length; otherwise a '+' in dna_change
    (e.g. '643A>C+644A>T') marks a multi-nucleotide change within one codon (MNV); the
    rest are single-nucleotide variants (SNV).
    """
    if (indel or "").strip():
        ref, alt = ref or "", alt or ""
        if len(alt) > len(ref):
            return "insertion"
        if len(alt) < len(ref):
            return "deletion"
        return "indel_other"
    if "+" in (dna_change or ""):
        return "mnv"
    return "snv"


def _codon_position(amino_acid_change):
    """Leading integer of a csq amino_acid_change ('83V>83I' -> 83); None if absent."""
    m = _LEADING_INT.match(amino_acid_change or "")
    return int(m.group(1)) if m else None


def _orf_codons(orf):
    """Number of codons in the ORF from a 'start-stop' nucleotide range; 0 if unparseable."""
    try:
        start, stop = str(orf).split("-")
        return max(0, (int(stop) - int(start)) // 3)
    except (ValueError, AttributeError):
        return 0


def _gini(values):
    """Gini coefficient of non-negative counts. 0 = perfectly even, ->1 = concentrated."""
    vals = sorted(values)
    n = len(vals)
    total = sum(vals)
    if n == 0 or total == 0:
        return 0.0
    cum = sum(i * v for i, v in enumerate(vals, 1))
    return round((2 * cum) / (n * total) - (n + 1) / n, 4)


def _coverage_profile(codon_counts, orf_codons, n_windows=N_COVERAGE_WINDOWS):
    """Down-sample per-codon counts into <=n_windows ordered windows along the ORF.

    Returns {window_start_codon: n_variants} with every window present (zeros kept) so
    uncovered stretches show as gaps in the profile plot.
    """
    if orf_codons <= 0:
        return {}
    size = max(1, math.ceil(orf_codons / n_windows))
    profile = {str(s): 0 for s in range(1, orf_codons + 1, size)}
    last_key = str(((orf_codons - 1) // size) * size + 1)
    for pos, n in codon_counts.items():
        key = str(((pos - 1) // size) * size + 1)
        profile[key if key in profile else last_key] += n
    return profile


def _read_cutadapt_counts(path):
    """(reads_total, reads_with_barcode) from a cutadapt --json report."""
    with open(path, encoding="UTF-8") as f:
        rc = json.load(f).get("read_counts", {})
    # With --discard-untrimmed, 'output' is the reads that kept a detected barcode.
    return rc.get("input"), rc.get("output")


def _to_int(value):
    """Parse an int from a possibly-blank TSV cell; blank/non-numeric -> 0."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def build_qc_summary(
    sample, cutadapt_json, cluster_qc, consensus_variants, variants, flagged,
    min_cluster_size, orf,
):
    reads_total, reads_with_barcode = _read_cutadapt_counts(cutadapt_json)

    sizes = []
    purity_hist = {b: 0 for b in PURITY_BINS}
    n_impure_clusters = 0
    barcodes = []
    with open(cluster_qc, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sizes.append(int(row["n_reads"]))
            try:
                smf = float(row.get("second_member_fraction", 0) or 0)
            except ValueError:
                smf = 0.0
            purity_hist[_bin_purity(smf)] += 1
            n_impure_clusters += smf > 0
            bc = (row.get("barcode") or "").strip()
            if bc:
                barcodes.append(bc)
    passing = [s for s in sizes if s >= min_cluster_size]
    # Exact size -> n_clusters, ordered by size; the plugin bins this across all samples.
    cluster_size_counts = {str(size): n for size, n in sorted(Counter(sizes).items())}

    # Barcode composition (axis C): length distribution, GC histogram, mean GC.
    barcode_length_counts = {
        str(k): n for k, n in sorted(Counter(len(b) for b in barcodes).items())
    }
    gc_hist = {b: 0 for b in GC_BINS}
    gc_values = []
    for b in barcodes:
        gc = sum(c in "GCgc" for c in b) / len(b)
        gc_hist[_bin_gc(gc)] += 1
        gc_values.append(gc)
    mean_barcode_gc = round(sum(gc_values) / len(gc_values), 4) if gc_values else 0.0

    # Single pass over the per-barcode variant records: allele fraction, consequence
    # classes, mutated positions/codons (coverage), variant complexity, variants-per-barcode.
    af_hist = {b: 0 for b in AF_BINS}
    consequences = Counter()
    variant_types = Counter()
    positions = set()
    codon_counts = Counter()
    n_indels = 0
    variants_per_barcode = Counter()
    with open(consensus_variants, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            v = row.get("af", "")
            if v not in ("", "None", None):
                try:
                    af_hist[_bin_af(float(v))] += 1
                except ValueError:
                    pass
            consequences[_consequence_class(row.get("consequence"))] += 1
            variant_types[_variant_type(
                row.get("INDEL"), row.get("REF"), row.get("ALT"), row.get("dna_change")
            )] += 1
            pos = row.get("POS")
            if pos:
                positions.add(pos)
            codon = _codon_position(row.get("amino_acid_change"))
            if codon:
                codon_counts[codon] += 1
            if (row.get("INDEL") or "").strip():
                n_indels += 1
            bc = row.get("barcode")
            if bc:
                variants_per_barcode[bc] += 1
    # Exact {n_variants_in_barcode: n_barcodes} over variant-bearing barcodes; plugin bins.
    variants_per_barcode_counts = {
        str(k): n for k, n in sorted(Counter(variants_per_barcode.values()).items())
    }

    # ORF coverage evenness (axis A): codon positions come from amino_acid_change; the Gini
    # is over the full ORF (uncovered codons counted as 0), the profile is down-sampled.
    orf_codons = _orf_codons(orf)
    n_codons_covered = len(codon_counts)
    frac_orf_covered = round(n_codons_covered / orf_codons, 4) if orf_codons else 0.0
    coverage_gini = (
        _gini([codon_counts.get(c, 0) for c in range(1, orf_codons + 1)])
        if orf_codons else 0.0
    )
    orf_coverage_profile = _coverage_profile(codon_counts, orf_codons)

    # Barcodes-per-variant redundancy (axis B): group the final mapping by variant identity
    # (dna_change, falling back to POS:REF>ALT) and count how many barcodes carry each.
    variant_barcodes = Counter()
    n_variants = 0
    with open(variants, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            n_variants += 1
            vid = (row.get("dna_change") or "").strip() or "{}:{}>{}".format(
                row.get("POS", ""), row.get("REF", ""), row.get("ALT", "")
            )
            variant_barcodes[vid] += 1
    barcodes_per_variant_counts = {
        str(k): n for k, n in sorted(Counter(variant_barcodes.values()).items())
    }

    n_merged = n_mixed = 0
    n_flagged_confident = n_flagged_ambiguous = 0
    with open(flagged, encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            reason = row.get("reason", "")
            n_merged += "merged" in reason
            n_mixed += "mixed" in reason
            n_flagged_confident += _to_int(row.get("n_confident"))
            n_flagged_ambiguous += _to_int(row.get("n_ambiguous"))

    return {
        "schema_version": SCHEMA_VERSION,
        "sample": sample,
        "reads_total": reads_total,
        "reads_with_barcode": reads_with_barcode,
        "n_clusters": len(sizes),
        "n_clusters_passing": len(passing),
        "median_cluster_size": int(statistics.median(passing)) if passing else 0,
        "n_variants": n_variants,
        "n_distinct_variants": len(variant_barcodes),
        "n_positions_mutated": len(positions),
        "n_codons_covered": n_codons_covered,
        "orf_codons": orf_codons,
        "frac_orf_covered": frac_orf_covered,
        "coverage_gini": coverage_gini,
        "n_indels": n_indels,
        "mean_barcode_gc": mean_barcode_gc,
        "n_flagged_merged": n_merged,
        "n_flagged_mixed": n_mixed,
        "n_impure_clusters": n_impure_clusters,
        "n_flagged_confident": n_flagged_confident,
        "n_flagged_ambiguous": n_flagged_ambiguous,
        "cluster_size_counts": cluster_size_counts,
        "allele_fraction_histogram": af_hist,
        "variant_consequences": dict(consequences),
        "variant_types": dict(variant_types),
        "variants_per_barcode_counts": variants_per_barcode_counts,
        "barcodes_per_variant_counts": barcodes_per_variant_counts,
        "cluster_purity_histogram": purity_hist,
        "orf_coverage_profile": orf_coverage_profile,
        "barcode_length_counts": barcode_length_counts,
        "barcode_gc_histogram": gc_hist,
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
        snakemake.params["orf"],
    )
    with open(snakemake.output[0], "w", encoding="UTF-8") as f:
        json.dump(summary, f, indent=2)
