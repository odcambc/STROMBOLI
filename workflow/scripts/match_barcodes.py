import csv
import sys
import os
csv.field_size_limit(sys.maxsize)

clusters_names = [
    "barcode",
    "n",
    "sequences"
]

clusters = {}

barcode_cluster_file = snakemake.input["clusters"]
barcode_variant_file = snakemake.input["variants"]
all_barcodes_out_file = snakemake.output[0]

# Open barcode cluster file and parse
with open(barcode_cluster_file, "r", encoding='UTF-8') as f:
    reader = csv.DictReader(f, delimiter="\t",
                            fieldnames=clusters_names,
                            quoting=csv.QUOTE_NONE)
    for row in reader:
        clusters[row["barcode"]] = row["sequences"].split(",")

# Read the variant consequences file and parse
with open(barcode_variant_file, "r", encoding='UTF-8') as f:
    with open(all_barcodes_out_file, "w", encoding='UTF-8') as out:
        reader = csv.reader(f, delimiter="\t")
        writer = csv.writer(out, delimiter="\t")
        for row in reader:
            if row[1] in clusters.keys():
                row[1] = ",".join([row[1]] + clusters[row[1]])
                writer.writerow(row)
            else:
                print("Not found: ", row[1], clusters[row[1]])