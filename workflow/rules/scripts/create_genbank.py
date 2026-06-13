from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

input_fasta = snakemake.input[0]
output_genbank = snakemake.output[0]
orf = snakemake.params["orf"]
gene_name = snakemake.params["gene_name"]
reference_name = snakemake.wildcards["reference_name"]

record = SeqIO.read(input_fasta, "fasta")
record.annotations = {"molecule_type": "DNA"}

orf_start = int(orf.split("-")[0])
orf_end = int(orf.split("-")[1])


features = []
features.extend(
    [
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="gene",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="CDS",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="transcript",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": gene_name},
            type="exon",
        ),
    ]
)
# Note: bcftools csq only needs gene/transcript/CDS/exon. Earlier versions added
# three_prime_UTR / five_prime_UTR features, but they were mis-placed (swapped, and
# inside the CDS) — removed rather than corrected, since csq does not use them.

record.features = features

SeqIO.write(record, output_genbank, "genbank")
