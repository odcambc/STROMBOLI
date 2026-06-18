import csv
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

clusters_names = ["barcode", "n", "sequences"]


def load_clusters(clusters_path):
    """Map each cluster's canonical barcode to all of its member barcodes.

    starcode --print-clusters output is TAB-delimited:
        canonical_barcode\tcount\tmember1,member2,...
    (the canonical barcode is included among the members).
    """
    clusters = {}
    with open(clusters_path, "r", encoding="UTF-8") as f:
        reader = csv.DictReader(
            f, delimiter="\t", fieldnames=clusters_names, quoting=csv.QUOTE_NONE
        )
        for row in reader:
            clusters[row["barcode"]] = row["sequences"].split(",")
    return clusters


def load_second_fraction(qc_path):
    """Map canonical barcode -> fraction of reads held by its 2nd-largest member."""
    frac = {}
    with open(qc_path, "r", encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            frac[row["barcode"]] = float(row["second_member_fraction"])
    return frac


def load_cluster_depth(qc_path):
    """Map canonical barcode -> cluster read depth (n_reads). This is the number of
    reads backing each variant call; downstream filtering to a target FDR keys on it,
    since per-barcode calls below ~10 reads are noise-dominated regardless of mode."""
    depth = {}
    with open(qc_path, "r", encoding="UTF-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            depth[row["barcode"]] = int(row["n_reads"])
    return depth


def _af(row):
    """Parse the ALT allele fraction; None when absent (double mode)."""
    v = row.get("af", "")
    if v in ("", "None", None):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def classify_and_join(
    clusters_path, variants_path, qc_path, out_main, out_flagged,
    clash_merge_fraction, clash_mixed_af, qc_min_af, exclude_clashes,
):
    """Build the final barcode->variant mapping, flagging and (optionally) excluding
    clashes:
      * merged : the cluster is a merge of distinct barcodes (its 2nd member holds
                 >= clash_merge_fraction of the reads).
      * mixed  : the barcode has a variant at intermediate allele fraction
                 (clash_mixed_af <= AF < qc_min_af) -- a likely collision.
    Confident variants of clean barcodes go to out_main (with the canonical barcode
    expanded to all members); flagged barcodes go to out_flagged with the reason.

    Limitation: a near-perfectly balanced collision (each variant ~50%) can be called
    reference at every site by the haploid caller, producing no variant to flag. Such
    a barcode is simply absent from the mapping (the safe outcome), not flagged. Robust
    detection of that case needs a per-position purity scan -- the same machinery as
    affirmative WT calling (see the roadmap).
    """
    clusters = load_clusters(clusters_path)
    second_frac = load_second_fraction(qc_path)
    cluster_depth = load_cluster_depth(qc_path)

    rows_by_bc = defaultdict(list)
    with open(variants_path, "r", encoding="UTF-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        var_fields = list(reader.fieldnames or [])
        for row in reader:
            rows_by_bc[row["barcode"]].append(row)

    main_fields = ["all_barcodes", "cluster_depth"] + var_fields
    flagged_fields = ["barcode", "reason", "second_member_fraction",
                      "n_confident", "n_ambiguous"]
    n_kept = n_flagged = 0

    with open(out_main, "w", encoding="UTF-8", newline="") as fm, open(
        out_flagged, "w", encoding="UTF-8", newline=""
    ) as ff:
        main = csv.DictWriter(fm, fieldnames=main_fields, delimiter="\t")
        main.writeheader()
        flagged = csv.DictWriter(ff, fieldnames=flagged_fields, delimiter="\t")
        flagged.writeheader()

        seen = set()
        for barcode, rows in rows_by_bc.items():
            seen.add(barcode)
            confident = [r for r in rows if _af(r) is None or _af(r) >= qc_min_af]
            ambiguous = [r for r in rows
                         if _af(r) is not None and clash_mixed_af <= _af(r) < qc_min_af]
            merged = second_frac.get(barcode, 0.0) >= clash_merge_fraction
            mixed = len(ambiguous) > 0

            if merged or mixed:
                n_flagged += 1
                reason = "+".join(
                    r for r, on in (("merged", merged), ("mixed", mixed)) if on)
                flagged.writerow({
                    "barcode": barcode, "reason": reason,
                    "second_member_fraction": round(second_frac.get(barcode, 0.0), 4),
                    "n_confident": len(confident), "n_ambiguous": len(ambiguous),
                })
                if exclude_clashes:
                    continue

            members = clusters.get(barcode, [barcode])
            depth = cluster_depth.get(barcode, "")
            for r in confident:
                main.writerow(
                    {"all_barcodes": ",".join(members), "cluster_depth": depth, **r})
            n_kept += 1

        # Also record merged clusters that produced no variant rows.
        for barcode, frac in second_frac.items():
            if barcode not in seen and frac >= clash_merge_fraction:
                n_flagged += 1
                flagged.writerow({
                    "barcode": barcode, "reason": "merged",
                    "second_member_fraction": round(frac, 4),
                    "n_confident": 0, "n_ambiguous": 0,
                })

    return n_kept, n_flagged


if "snakemake" in globals():
    p = snakemake.params
    k, fl = classify_and_join(
        snakemake.input["clusters"],
        snakemake.input["variants"],
        snakemake.input["qc"],
        snakemake.output["mapping"],
        snakemake.output["flagged"],
        clash_merge_fraction=p["clash_merge_fraction"],
        clash_mixed_af=p["clash_mixed_af"],
        qc_min_af=p["qc_min_af"],
        exclude_clashes=p["exclude_clashes"],
    )
    with open(snakemake.log[0], "w") as log:
        log.write("kept {} barcodes, flagged {} as clashes\n".format(k, fl))
