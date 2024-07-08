import csv
import sys
import logging
import tqdm

csv.field_size_limit(sys.maxsize)

names = [
    "read",
    "n_errors",
    "start",
    "end",
    "seq_l",
    "seq_match",
    "seq_r",
    "name",
    "q_l",
    "q_match",
    "q_r",
]

# TODO: this is a very infefficent script. It needs to be made more memory and thread efficient.

clusters_names = ["barcode", "n", "sequences"]

reads = {}
clusters = {}
barcodes = {}

input_info_file = snakemake.input["info"]
barcode_clusters_file = snakemake.input["barcode_clusters"]

output_clusters_file = snakemake.output[0]
log_file = snakemake.log[0]

# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)


logging.debug("Reading in input info file: {}".format(input_info_file))

# Open tsv and parse
with open(input_info_file, "r", encoding="UTF-8") as f:
    reader = csv.DictReader(f, delimiter="\t", fieldnames=names, quoting=csv.QUOTE_NONE)
    for row in tqdm.tqdm(reader):
        if row["read"] not in reads.keys():
            reads[row["read"]] = {"barcode": "", "seq": ""}
        if row["name"] == "1;1":
            reads[row["read"]]["seq"] = row["seq_l"]
        elif row["name"] == "1;2":
            reads[row["read"]]["barcode"] = row["seq_l"]

logging.debug("Read in {} reads".format(len(reads)))

# Make a dictionary of barcodes
for read, values in tqdm.tqdm(reads.items()):
    if values["barcode"] != "":
        if values["barcode"] not in barcodes.keys():
            barcodes[values["barcode"]] = [values["seq"]]
        else:
            barcodes[values["barcode"]].append(values["seq"])

logging.debug("Read in {} barcodes".format(len(barcodes)))

# Read in barcode clusters
with open(barcode_clusters_file, "r", encoding="UTF-8") as f:
    reader = csv.DictReader(
        f, delimiter="\t", fieldnames=clusters_names, quoting=csv.QUOTE_NONE
    )
    for row in tqdm.tqdm(reader):
        clusters[row["barcode"]] = {
            "n": row["n"],
            "sequences": row["sequences"].split(","),
        }

logging.debug("Read in {} barcode clusters".format(len(clusters)))
logging.debug("Writing out sequences to {}".format(output_clusters_file))

# Write out sequences for each barcode cluster
with open(output_clusters_file, "w", encoding="UTF-8") as f:
    for barcode, values in tqdm.tqdm(clusters.items()):
        f.write(barcode + ",")
        sequences = []
        for seq in values["sequences"]:
            sequences += barcodes[seq]
        n_seq = len(sequences)
        f.write(str(n_seq) + ",")
        f.write(",".join(sequences) + "\n")
