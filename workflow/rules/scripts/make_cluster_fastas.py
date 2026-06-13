import csv
import sys
import os

csv.field_size_limit(sys.maxsize)

cluster_file = snakemake.input[0]
output_dir = os.path.abspath(snakemake.output[0])
min_size = snakemake.params["min_cluster_size"]

# Entries are "<sequence>|<quality>" (see write_sequences.py).
ENTRY_SEP = "|"

clusters = {}

# The clusters file is TAB-delimited: barcode, n_reads, then one entry per read.
with open(cluster_file, "r", encoding="UTF-8") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if int(row[1]) >= min_size:
            clusters[row[0]] = row[2:]

os.makedirs(output_dir, exist_ok=True)

# Write a valid FASTQ record per read so downstream consensus calling can use the
# base qualities (samtools consensus -q).
for cluster, entries in clusters.items():
    with open(
        os.path.join(output_dir, cluster + ".fastq"), "w", encoding="UTF-8"
    ) as f:
        for count, entry in enumerate(entries):
            seq, _, qual = entry.partition(ENTRY_SEP)
            f.write("@{}_{}\n{}\n+\n{}\n".format(cluster, count, seq, qual))
