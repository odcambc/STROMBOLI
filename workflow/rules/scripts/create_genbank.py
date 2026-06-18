from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os


def validate_orf_frame(seq, cds_start, orf_end, orf):
    """Guard against an out-of-frame ORF and fail loudly. Translates the complete-codon
    span of the CDS (case-insensitively, so soft-masked references are fine) and requires
    a clean reading frame: an ATG/M start and no internal stop codons.

    This catches an off-by-one in either the orf config (e.g. "199-3237" vs "198-3237")
    or the 1-based->0-based conversion. Without it, a frame error is SILENT: barcode/POS
    recovery still works, but bcftools csq translates the gene out of frame and every
    amino-acid consequence is wrong. The trailing partial codon (when the span isn't a
    multiple of 3) is dropped before checking, so a loose end bound doesn't false-alarm."""
    span = orf_end - cds_start
    coding = Seq(str(seq[cds_start:cds_start + span - span % 3]))
    prot = str(coding.translate())
    problems = []
    if not prot.startswith("M"):
        problems.append("CDS does not begin with a start codon (ATG/M)")
    internal_stops = prot.rstrip("*").count("*")
    if internal_stops:
        problems.append(f"{internal_stops} internal stop codon(s)")
    if problems:
        raise ValueError(
            f"ORF {orf!r} appears out of frame on the reference: "
            + "; ".join(problems)
            + ". Check the orf bounds (1-based, inclusive) and the reading frame."
        )


def build_record(record, orf, gene_name, reference_name):
    """Annotate `record` (a SeqRecord) with gene/CDS/transcript/exon features for the
    ORF, returning it. `orf` is the "start-end" config string, 1-based inclusive.

    The `orf` param is 1-based inclusive (e.g. "198-3237"), but Biopython's
    FeatureLocation is 0-based half-open. We convert the start by -1; the half-open
    end equals the 1-based inclusive end, so orf_end is used unchanged. On GFF/GenBank
    export Biopython re-adds 1 to the start, recovering the 1-based CDS start. Omitting
    this -1 shifts the reading frame by one base, so bcftools csq translates the ORF
    out of frame and every consequence call is wrong (see test_create_genbank_frame).
    The reading frame is validated up front so a bad orf fails loudly at reference build.
    """
    orf_start = int(orf.split("-")[0])
    orf_end = int(orf.split("-")[1])
    cds_start = orf_start - 1
    validate_orf_frame(record.seq, cds_start, orf_end, orf)

    record.features = [
        SeqFeature(
            FeatureLocation(cds_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="gene",
        ),
        SeqFeature(
            FeatureLocation(cds_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="CDS",
        ),
        SeqFeature(
            FeatureLocation(cds_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="transcript",
        ),
        SeqFeature(
            FeatureLocation(cds_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="exon",
        ),
    ]
    # Note: bcftools csq only needs gene/transcript/CDS/exon. Earlier versions added
    # three_prime_UTR / five_prime_UTR features, but they were mis-placed (swapped, and
    # inside the CDS) — removed rather than corrected, since csq does not use them.
    return record


if "snakemake" in globals():
    input_fasta = snakemake.input[0]
    output_genbank = snakemake.output[0]
    record = SeqIO.read(input_fasta, "fasta")
    record.annotations = {"molecule_type": "DNA"}
    record = build_record(
        record,
        snakemake.params["orf"],
        snakemake.params["gene_name"],
        snakemake.wildcards["reference_name"],
    )
    SeqIO.write(record, output_genbank, "genbank")
