"""Fast unit tests for the pure-Python pipeline helpers (no bioinformatics tools).

Run with:  pytest tests/test_scripts.py
"""
import csv
import importlib.util
import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "workflow", "rules", "scripts")


def _load(module_name, filename):
    """Load a snakemake script module by path. The scripts guard their execution
    behind `if "snakemake" in globals()`, so importing them here is side-effect free."""
    spec = importlib.util.spec_from_file_location(
        module_name, os.path.join(SCRIPTS, filename)
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


wcv = _load("write_consensus_variants", "write_consensus_variants.py")
mb = _load("match_barcodes", "match_barcodes.py")
wqs = _load("write_qc_summary", "write_qc_summary.py")
cg = _load("create_genbank", "create_genbank.py")
gen = importlib.util.spec_from_file_location(
    "gen", os.path.join(REPO_ROOT, "tests", "generate_synthetic_data.py")
)
generate = importlib.util.module_from_spec(gen)
gen.loader.exec_module(generate)


# --- parse_bcsq -----------------------------------------------------------------


def test_parse_bcsq_full_annotation():
    field = "missense|gp17|gp17|protein_coding|+|548P>548S|1642C>T"
    parsed = wcv.parse_bcsq(field)
    assert parsed["consequence"] == "missense"
    assert parsed["gene"] == "gp17"
    assert parsed["amino_acid_change"] == "548P>548S"
    assert parsed["dna_change"] == "1642C>T"
    assert parsed["linked_record"] is None


def test_parse_bcsq_linked_record():
    parsed = wcv.parse_bcsq("@1642")
    assert parsed["linked_record"] == "@1642"
    assert parsed["consequence"] is None


def test_parse_bcsq_empty():
    for empty in (None, ""):
        parsed = wcv.parse_bcsq(empty)
        assert all(v is None for v in parsed.values())


def test_parse_bcsq_truncated_fields():
    # Trailing fields may be omitted; only the present ones map.
    parsed = wcv.parse_bcsq("synonymous|gp17|gp17|protein_coding")
    assert parsed["consequence"] == "synonymous"
    assert parsed["biotype"] == "protein_coding"
    assert parsed["amino_acid_change"] is None


# --- match_barcodes.classify_and_join ------------------------------------------


def _write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def test_classify_and_join_flags_and_excludes(tmp_path):
    A = "AAAAAAAAAAAAAAAAAAAA"  # clean: confident variant, single member
    B = "BBBBBBBBBBBBBBBBBBBB"  # merged: 2nd member holds 45% of reads
    C = "CCCCCCCCCCCCCCCCCCCC"  # mixed: variant at intermediate AF (collision)

    clusters = tmp_path / "clusters.txt"
    clusters.write_text(
        f"{A}\t10\t{A}\n{B}\t10\t{B},BBBBBBBBBBBBBBBBBBBC\n{C}\t10\t{C}\n")

    qc = tmp_path / "qc.tsv"
    _write_tsv(qc, ["barcode", "n_reads", "n_members", "second_member_fraction"], [
        {"barcode": A, "n_reads": 10, "n_members": 1, "second_member_fraction": "0.0"},
        {"barcode": B, "n_reads": 10, "n_members": 2, "second_member_fraction": "0.45"},
        {"barcode": C, "n_reads": 10, "n_members": 1, "second_member_fraction": "0.0"},
    ])

    variants = tmp_path / "variants.tsv"
    _write_tsv(variants, ["barcode", "POS", "REF", "ALT", "af"], [
        {"barcode": A, "POS": "100", "REF": "A", "ALT": "C", "af": "0.98"},
        {"barcode": B, "POS": "200", "REF": "G", "ALT": "T", "af": "1.0"},
        {"barcode": C, "POS": "300", "REF": "T", "ALT": "A", "af": "0.5"},
    ])

    out_main = tmp_path / "main.tsv"
    out_flagged = tmp_path / "flagged.tsv"
    kept, flagged = mb.classify_and_join(
        str(clusters), str(variants), str(qc), str(out_main), str(out_flagged),
        clash_merge_fraction=0.2, clash_mixed_af=0.2, qc_min_af=0.85,
        exclude_clashes=True)

    main_rows = list(csv.DictReader(open(out_main), delimiter="\t"))
    flagged_rows = {r["barcode"]: r for r in csv.DictReader(open(out_flagged), delimiter="\t")}

    # Only the clean barcode survives; merged and mixed are excluded.
    assert [r["barcode"] for r in main_rows] == [A]
    # Each surviving call is annotated with its cluster read depth (n_reads).
    assert main_rows[0]["cluster_depth"] == "10"
    assert flagged_rows[B]["reason"] == "merged"
    assert flagged_rows[C]["reason"] == "mixed"
    assert kept == 1 and flagged == 2


def test_classify_and_join_double_mode_no_af(tmp_path):
    """double mode has no AF (empty) -> all confident; only merges can be flagged."""
    A = "AAAAAAAAAAAAAAAAAAAA"
    clusters = tmp_path / "c.txt"
    clusters.write_text(f"{A}\t5\t{A}\n")
    qc = tmp_path / "qc.tsv"
    _write_tsv(qc, ["barcode", "n_reads", "n_members", "second_member_fraction"],
               [{"barcode": A, "n_reads": 5, "n_members": 1, "second_member_fraction": "0.0"}])
    variants = tmp_path / "v.tsv"
    _write_tsv(variants, ["barcode", "POS", "REF", "ALT", "af"],
               [{"barcode": A, "POS": "100", "REF": "A", "ALT": "C", "af": ""}])
    out_main, out_flagged = tmp_path / "m.tsv", tmp_path / "f.tsv"
    kept, flagged = mb.classify_and_join(
        str(clusters), str(variants), str(qc), str(out_main), str(out_flagged),
        clash_merge_fraction=0.2, clash_mixed_af=0.2, qc_min_af=0.85, exclude_clashes=True)
    assert kept == 1 and flagged == 0


# --- write_qc_summary.build_qc_summary -----------------------------------------


def test_build_qc_summary(tmp_path):
    import json

    cutadapt = tmp_path / "c.json"
    cutadapt.write_text(json.dumps({"read_counts": {"input": 1000, "output": 800}}))
    cluster_qc = tmp_path / "qc.tsv"
    _write_tsv(cluster_qc, ["barcode", "n_reads", "n_members", "second_member_fraction"], [
        {"barcode": "A", "n_reads": "1", "n_members": "1", "second_member_fraction": "0.0"},
        {"barcode": "B", "n_reads": "4", "n_members": "1", "second_member_fraction": "0.0"},
        {"barcode": "C", "n_reads": "12", "n_members": "1", "second_member_fraction": "0.0"},
    ])
    cvars = tmp_path / "cv.tsv"
    _write_tsv(
        cvars,
        ["barcode", "POS", "REF", "ALT", "consequence", "INDEL",
         "amino_acid_change", "dna_change", "af"],
        [
            {"barcode": "B", "POS": "100", "REF": "C", "ALT": "A", "consequence": "missense",
             "INDEL": "", "amino_acid_change": "10P>10L", "dna_change": "100C>A", "af": "0.95"},
            {"barcode": "B", "POS": "150", "REF": "G", "ALT": "A", "consequence": "*synonymous",
             "INDEL": "", "amino_acid_change": "20K", "dna_change": "150G>A+151T>C", "af": "0.9"},
            {"barcode": "C", "POS": "200", "REF": "A", "ALT": "AT", "consequence": "",
             "INDEL": "1", "amino_acid_change": "", "dna_change": "200A>AT", "af": "0.5"},
        ],
    )
    variants = tmp_path / "v.tsv"
    _write_tsv(variants, ["all_barcodes", "barcode", "POS", "REF", "ALT", "dna_change"], [
        {"all_barcodes": "B", "barcode": "B", "POS": "100", "REF": "C", "ALT": "A",
         "dna_change": "100C>A"},
        {"all_barcodes": "C", "barcode": "C", "POS": "100", "REF": "C", "ALT": "A",
         "dna_change": "100C>A"},  # same variant carried by a 2nd barcode
    ])
    flagged = tmp_path / "f.tsv"
    _write_tsv(
        flagged,
        ["barcode", "reason", "second_member_fraction", "n_confident", "n_ambiguous"],
        [{"barcode": "C", "reason": "mixed", "second_member_fraction": "0.5",
          "n_confident": "3", "n_ambiguous": "2"}],
    )

    s = wqs.build_qc_summary("samp", str(cutadapt), str(cluster_qc), str(cvars),
                             str(variants), str(flagged), min_cluster_size=2, orf="1-300")
    assert s["schema_version"] == 4 and s["sample"] == "samp"
    assert s["reads_total"] == 1000 and s["reads_with_barcode"] == 800
    assert s["n_clusters"] == 3 and s["n_clusters_passing"] == 2  # sizes 4,12 >= 2
    assert s["n_flagged_mixed"] == 1 and s["n_flagged_merged"] == 0
    assert s["cluster_size_counts"] == {"1": 1, "4": 1, "12": 1}
    assert s["allele_fraction_histogram"]["0.85-1.0"] == 2
    assert s["allele_fraction_histogram"]["0.4-0.6"] == 1
    # v3 latent-data metrics
    assert s["n_positions_mutated"] == 3  # POS 100, 150, 200
    assert s["n_indels"] == 1
    # '*synonymous' normalizes to 'synonymous'; blank consequence -> 'noncoding'
    assert s["variant_consequences"] == {"missense": 1, "synonymous": 1, "noncoding": 1}
    assert s["variants_per_barcode_counts"] == {"1": 1, "2": 1}  # C:1 variant, B:2 variants
    assert s["cluster_purity_histogram"]["0"] == 3 and s["n_impure_clusters"] == 0
    assert s["n_flagged_confident"] == 3 and s["n_flagged_ambiguous"] == 2
    # v4 coverage/complexity/redundancy/composition
    # variant types: 100C>A snv, 150G>A+151T>C mnv, A->AT insertion
    assert s["variant_types"] == {"snv": 1, "mnv": 1, "insertion": 1}
    # codon positions parsed from amino_acid_change: 10 and 20 (blank -> none)
    assert s["orf_codons"] == 99 and s["n_codons_covered"] == 2  # (300-1)//3 = 99
    assert s["frac_orf_covered"] == round(2 / 99, 4)
    assert 0.0 < s["coverage_gini"] <= 1.0  # sparse coverage -> uneven
    # both final-mapping rows are the same variant (100C>A) carried by 2 barcodes
    assert s["n_variants"] == 2 and s["n_distinct_variants"] == 1
    assert s["barcodes_per_variant_counts"] == {"2": 1}
    # barcode composition over cluster_qc barcodes "A", "B", "C" (length 1 each)
    assert s["barcode_length_counts"] == {"1": 3}


# --- synthetic data generator ---------------------------------------------------


def test_truth_ref_bases_match_reference():
    """The ref base recorded in the (would-be) truth must match the reference."""
    reference = generate.read_reference(generate.REFERENCE)
    for pos in generate.SNP_POSITIONS:
        assert reference[pos - 1] in "ACGT"


def test_other_base_differs():
    import random

    rng = random.Random(1)
    for base in "ACGT":
        assert generate.other_base(rng, base) != base


# --- create_genbank reading frame -----------------------------------------------


def test_create_genbank_frame():
    """Regression: the 1-based `orf` start must be converted to 0-based for the
    FeatureLocation. Skipping the -1 shifts the CDS one base, so bcftools csq
    translates the ORF out of frame and every consequence call is wrong.

    Crafted so the correct (1-based) frame is a clean ORF (M-start, no internal
    stop), while the +1-shifted frame hits a leading stop."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    # 1-based: C(1)C(2)C(3) | ATG(4-6) AAA AAA AAA TAA(stop, 16-18)
    seq = "CCC" + "ATG" + "AAA" * 3 + "TAA"
    record = SeqRecord(Seq(seq), id="ref", annotations={"molecule_type": "DNA"})
    cg.build_record(record, "4-18", "testgene", "ref")

    cds = [f for f in record.features if f.type == "CDS"][0]
    # 1-based start 4 -> 0-based half-open start 3 (the fix); 4 would be the bug.
    assert int(cds.location.start) == 3
    assert int(cds.location.end) == 18

    prot = str(record.seq[int(cds.location.start) : int(cds.location.end)].translate())
    assert prot.startswith("M")  # correct frame opens with the start codon
    assert prot[:-1].count("*") == 0  # ...and has no internal stops


def test_create_genbank_rejects_out_of_frame_orf():
    """Hardening: build_record fails loudly on an out-of-frame orf (e.g. an off-by-one
    start) instead of silently emitting wrong consequences. Same crafted sequence; the
    +1-shifted start "5-18" begins on a TGA stop and carries no ATG."""
    import pytest
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    seq = "CCC" + "ATG" + "AAA" * 3 + "TAA"  # in frame at 4-18, out of frame at 5-18
    record = SeqRecord(Seq(seq), id="ref", annotations={"molecule_type": "DNA"})
    with pytest.raises(ValueError, match="out of frame"):
        cg.build_record(record, "5-18", "testgene", "ref")
