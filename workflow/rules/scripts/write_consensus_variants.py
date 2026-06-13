import os
import csv

# NOTE: cyvcf2 is imported lazily inside the snakemake-execution block at the
# bottom so this module's pure helpers (parse_bcsq) can be imported in tests
# without the bioinformatics dependencies installed.

# bcftools csq BCSQ field format:
#   Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
# Records downstream of a compound variant instead carry "@<position>" linking
# back to the record that holds the full annotation.

vcf_keys = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
info_keys = [
    "INDEL",
    "IDV",
    "IMF",
    "DP",
    "VDB",
    "RPBZ",
    "MQBZ",
    "BQBZ",
    "MQSBZ",
    "SCBZ",
    "SGB",
    "MQ0F",
    "AC",
    "AN",
    "DP4",
    "MQ",
    "BCSQ",
]
bcsq_keys = [
    "consequence",
    "gene",
    "transcript",
    "biotype",
    "strand",
    "amino_acid_change",
    "dna_change",
    "linked_record",
]


def parse_bcsq(info_str):
    """Parse a bcftools/csq BCSQ INFO string into the bcsq_keys fields.

    Handles the three shapes seen in practice:
      * a full "Consequence|gene|transcript|..." annotation (trailing fields may
        be omitted when empty),
      * a "@<position>" link to another record (only linked_record is set),
      * an empty/None field (all values None).
    A field may hold several comma-separated consequences; the first is used.
    """
    value_dict = {key: None for key in bcsq_keys}
    if not info_str:
        return value_dict

    first = info_str.split(",")[0]
    if first.startswith("@"):
        value_dict["linked_record"] = first
        return value_dict

    for key, val in zip(bcsq_keys, first.split("|")):
        value_dict[key] = val
    return value_dict


def _alt_fraction(variant):
    """ALT allele fraction from FORMAT/AD, or "" when AD is absent (double mode)."""
    try:
        ad = variant.format("AD")
        if ad is None:
            return ""
        ref_ad, alt_ad = int(ad[0][0]), int(ad[0][1])
        if ref_ad < 0 or alt_ad < 0 or ref_ad + alt_ad == 0:
            return ""
        return round(alt_ad / (ref_ad + alt_ad), 4)
    except Exception:
        return ""


def write_variants(barcode_bcf_files, output_path):
    """Aggregate per-barcode csq BCFs into one TSV keyed by canonical barcode."""
    from cyvcf2 import VCF

    with open(output_path, "w", encoding="UTF-8", newline="") as f:
        # `af` (ALT allele fraction, from FORMAT/AD) is empty in double mode, where the
        # call is on a depth-1 consensus; single_qc populates it and downstream uses it
        # to flag mixed (collision) barcodes.
        headers = ["barcode"] + vcf_keys + info_keys + bcsq_keys + ["af"]
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()

        for barcode_bcf in barcode_bcf_files:
            # Files are named "<barcode>_csq.bcf"; recover the canonical barcode.
            barcode_sequence = os.path.basename(barcode_bcf).split("_")[0]

            for variant in VCF(os.path.abspath(barcode_bcf)):
                record_dict = {"barcode": barcode_sequence, "af": _alt_fraction(variant)}
                record_dict["CHROM"] = variant.CHROM
                record_dict["POS"] = variant.POS
                record_dict["ID"] = variant.ID
                record_dict["REF"] = variant.REF
                record_dict["ALT"] = ",".join(variant.ALT)
                record_dict["QUAL"] = variant.QUAL
                record_dict["FILTER"] = variant.FILTER

                info = variant.INFO
                if info:
                    for key in info_keys:
                        field = info.get(key)
                        if key == "BCSQ":
                            record_dict["BCSQ"] = field
                            record_dict.update(parse_bcsq(field))
                        else:
                            record_dict[key] = field
                else:
                    record_dict.update({key: None for key in info_keys})
                    record_dict.update({key: None for key in bcsq_keys})

                writer.writerow(record_dict)


if "snakemake" in globals():
    write_variants(list(snakemake.input), snakemake.output[0])
