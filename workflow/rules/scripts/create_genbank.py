from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os

input_fasta = snakemake.input[0]
output_genbank = snakemake.output[0]
orf = snakemake.params["orf"]
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
            qualifiers={"locus_tag": "gp17"},
            type="gene",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": "gp17"},
            type="CDS",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": "gp17"},
            type="transcript",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": "gp17"},
            type="exon",
        ),
        SeqFeature(
            FeatureLocation(orf_start, orf_start + 9, strand=1, ref=reference_name),
            qualifiers={"locus_tag": "gp17"},
            type="three_prime_UTR",
        ),
        SeqFeature(
            FeatureLocation(orf_end - 9, orf_end, strand=1, ref=reference_name),
            qualifiers={"locus_tag": "gp17"},
            type="five_prime_UTR",
        ),
    ]
)

record.features = features

SeqIO.write(record, output_genbank, "genbank")


"""gp17_orf	ignored_field	gene	1	3195	.	+	.	ID=gene:gp17;biotype=protein_coding;Name=gp17
gp17_orf	ignored_field	transcript	1	3195	.	+	.	ID=transcript:gp17;Parent=gene:gp17;biotype=protein_coding
gp17_orf	ignored_field	exon	1	3195	.	+	.	Parent=transcript:gp17
gp17_orf	ignored_field	three_prime_UTR	1	9	.	+	.	Parent=transcript:gp17
gp17_orf	ignored_field	CDS	1	3195	.	+	0	Parent=transcript:gp17
gp17_orf	ignored_field	five_prime_UTR	3000	3195	.	+	0	Parent=transcript:gp17
"""
