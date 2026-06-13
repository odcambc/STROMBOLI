import csv
import sys
import logging
import tqdm

csv.field_size_limit(sys.maxsize)

# Column names of the cutadapt --info-file output.
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

clusters_names = ["barcode", "n", "sequences"]

reads = {}
clusters = {}
barcodes = {}

input_info_file = snakemake.input["info"]
barcode_clusters_file = snakemake.input["barcode_clusters"]

output_clusters_file = snakemake.output["clusters"]
output_qc_file = snakemake.output["qc"]
log_file = snakemake.log[0]

# Each cluster entry is written as "<sequence>|<quality>" in a TAB-delimited file.
# Tab never occurs in a Phred quality string (ASCII 9 < the Phred range), and "|"
# (ASCII 124 == Q91) is above any quality nanopore produces, so neither separator
# collides with the data. (Comma, the old delimiter, IS a valid quality char.)
ENTRY_SEP = "|"


# Set up logging
if log_file:
    logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.DEBUG)


logging.debug("Reading in input info file: {}".format(input_info_file))

# Open the cutadapt info tsv and reconstruct, per read, the insert sequence and
# its qualities (the "1;1" linked-adapter row) and the barcode (the "1;2" row).
with open(input_info_file, "r", encoding="UTF-8") as f:
    reader = csv.DictReader(f, delimiter="\t", fieldnames=names, quoting=csv.QUOTE_NONE)
    for row in tqdm.tqdm(reader):
        if row["read"] not in reads:
            reads[row["read"]] = {"barcode": "", "seq": "", "qual": ""}
        if row["name"] == "1;1":
            reads[row["read"]]["seq"] = row["seq_l"]
            reads[row["read"]]["qual"] = row["q_l"]
        elif row["name"] == "1;2":
            reads[row["read"]]["barcode"] = row["seq_l"]

logging.debug("Read in {} reads".format(len(reads)))

# Group (sequence, quality) pairs by their barcode. Skip reads that are missing
# either the barcode or the insert sequence.
for read, values in tqdm.tqdm(reads.items()):
    if values["barcode"] != "" and values["seq"] != "":
        entry = values["seq"] + ENTRY_SEP + values["qual"]
        barcodes.setdefault(values["barcode"], []).append(entry)

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

# Write one tab-delimited row per cluster:
#   <canonical barcode>\t<n reads>\t<seq|qual>\t<seq|qual>...
# Fields never contain a tab (seq is ACGTN; quality chars are all below ASCII 9's
# Phred range), so a plain tab-join is unambiguous and avoids csv quoting.
#
# Alongside, write per-cluster QC: how the reads are split across the cluster's member
# barcodes. An error cloud is one dominant member; a MERGE of two genuinely distinct
# barcodes shows a second member holding a large read share (see clash_merge_fraction).
with open(output_clusters_file, "w", encoding="UTF-8") as f, open(
    output_qc_file, "w", encoding="UTF-8"
) as qc:
    qc.write("barcode\tn_reads\tn_members\tsecond_member_fraction\n")
    for barcode, values in tqdm.tqdm(clusters.items()):
        member_counts = [len(barcodes.get(m, [])) for m in values["sequences"]]
        entries = []
        for member in values["sequences"]:
            entries.extend(barcodes.get(member, []))
        total = len(entries)
        nonzero = sorted((c for c in member_counts if c > 0), reverse=True)
        second_frac = (nonzero[1] / total) if len(nonzero) >= 2 and total else 0.0
        qc.write(
            "{}\t{}\t{}\t{:.4f}\n".format(barcode, total, len(nonzero), second_frac)
        )
        f.write("\t".join([barcode, str(total)] + entries) + "\n")
