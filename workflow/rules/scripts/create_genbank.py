from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os


def build_record(record, orf, gene_name, reference_name):
    """Annotate `record` (a SeqRecord) with gene/CDS/transcript/exon features for the
    ORF, returning it. `orf` is the "start-end" config string, 1-based inclusive.

    The `orf` param is 1-based inclusive (e.g. "198-3237"), but Biopython's
    FeatureLocation is 0-based half-open. We convert the start by -1; the half-open
    end equals the 1-based inclusive end, so orf_end is used unchanged. On GFF/GenBank
    export Biopython re-adds 1 to the start, recovering the 1-based CDS start. Omitting
    this -1 shifts the reading frame by one base, so bcftools csq translates the ORF
    out of frame and every consequence call is wrong (see test_create_genbank_frame).
    """
    orf_start = int(orf.split("-")[0])
    orf_end = int(orf.split("-")[1])
    cds_start = orf_start - 1

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
