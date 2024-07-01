import csv
import sys
import os
csv.field_size_limit(sys.maxsize)

cluster_file = snakemake.input[0]

output_dir = os.path.dirname(snakemake.output[0])

clusters = {}
min_size = snakemake.params["min_size"]

with open(cluster_file, 'r', encoding="UTF-8") as f:
    reader = csv.reader(f)
    for row in reader:
        if int(row[1]) >= min_size:
            clusters[row[0]] = row[2:]

for cluster in clusters:
    with open(output_dir + "/" + cluster + ".fasta", 'w+', encoding="UTF-8") as f:
        for count, sequence in enumerate(clusters[cluster]):
            f.write(">" + cluster + "_" + str(count) + "\n" + sequence + "\n")