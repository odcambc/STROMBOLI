import csv
import sys

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

clusters_names = ["barcode", "n", "sequences"]

reads = {}

min_count = snakemake.params["min_count"]
input_info_file = snakemake.input[0]
output_barcodes_file = snakemake.output[0]

# Open tsv and parse
with open(input_info_file, "r", encoding="UTF-8") as f:
    reader = csv.DictReader(f, delimiter="\t", fieldnames=names, quoting=csv.QUOTE_NONE)
    try:
        for row in reader:
            if row["read"] not in reads.keys():
                reads[row["read"]] = {"barcode": "", "seq": ""}
            if row["name"] == "1;1":
                reads[row["read"]]["seq"] = row["seq_l"]
            elif row["name"] == "1;2":
                reads[row["read"]]["barcode"] = row["seq_l"]
    except csv.Error:
        print("Error with row: ", row)


with open(output_barcodes_file, "w", encoding="UTF-8") as f:
    for read, values in reads.items():
        try:
            if values["barcode"] != "" and values["seq"] != "":
                if len(values["barcode"]) > min_count:  # Or < min_count??
                    f.write(values["barcode"] + "\n")
        except TypeError:
            print(read, reads[read])
