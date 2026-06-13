"""End-to-end integration test on the synthetic dataset.

This runs the whole Snakemake pipeline and asserts that the final
results/synthetic.variants.tsv recovers the known barcode -> SNP truth table.

It is skipped automatically unless the full toolchain (snakemake + minimap2 +
samtools + bcftools + starcode + cutadapt) is available on PATH, e.g. inside the
conda env. Run explicitly with:  pytest -m integration tests/

Marked `integration` (and `slow`) so it can be excluded from quick runs:
    pytest -m "not integration"
"""
import csv
import importlib.util
import os
import shutil
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REQUIRED_TOOLS = [
    "snakemake",
    "minimap2",
    "samtools",
    "bcftools",
    "starcode",
    "cutadapt",
]

pytestmark = [pytest.mark.integration, pytest.mark.slow]


def _tools_available():
    return all(shutil.which(t) for t in REQUIRED_TOOLS)


def _load_truth(path):
    with open(path) as f:
        return list(csv.DictReader(f))


@pytest.mark.skipif(not _tools_available(), reason="full toolchain not installed")
@pytest.mark.parametrize("calling_mode", ["double", "single_qc"])
def test_pipeline_recovers_synthetic_truth(calling_mode):
    # 1. Generate the deterministic synthetic dataset.
    spec = importlib.util.spec_from_file_location(
        "gen", os.path.join(REPO_ROOT, "tests", "generate_synthetic_data.py")
    )
    generate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generate)
    generate.main()

    truth = _load_truth(os.path.join(REPO_ROOT, "tests", "synthetic_truth.csv"))

    # Clear prior outputs so each mode rebuilds cleanly (the two modes share output
    # paths). We only remove results/ -- never references/, which ship pre-built.
    shutil.rmtree(os.path.join(REPO_ROOT, "results"), ignore_errors=True)

    # 2. Run the pipeline on the synthetic config in the given calling mode. Both
    #    modes must recover the clean-SNP truth set on this (low-error) data.
    subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            "config/synthetic.yaml",
            "--config",
            f"calling_mode={calling_mode}",
            "--cores",
            "4",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    # 3. Load the final barcode -> variant table.
    out_path = os.path.join(REPO_ROOT, "results", "synthetic.variants.tsv")
    assert os.path.exists(out_path), "pipeline did not produce the variants table"
    with open(out_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # 4. Each truth SNP must be recovered for its barcode.
    misses = []
    for entry in truth:
        barcode = entry["barcode"]
        want_pos, want_alt = entry["pos"], entry["alt"].upper()
        found = any(
            barcode in (row.get("all_barcodes") or "")
            and str(row.get("POS")) == want_pos
            and (row.get("ALT") or "").upper() == want_alt
            for row in rows
        )
        if not found:
            misses.append((barcode, want_pos, want_alt))

    assert not misses, "Unrecovered truth SNPs: {}".format(misses)
