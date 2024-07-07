import sys
import os
import csv
from cyvcf2 import VCF


input_barcode_names = snakemake.input
output_consensus_file = snakemake.output[0]

sample = snakemake.wildcards["sample"]

"""with open(output_consensus_file, "w", encoding="UTF-8") as f:
    for barcode in input_barcode_names:
        barcode_sequence = barcode.split("/")[-1].split(".")[0]
        variant_file = os.path.abspath(barcode)
        with open(variant_file, "r", encoding="UTF-8") as v:
            variant_list = v.readlines()
            variant_sequence = "".join([line.strip() for line in variant_list[1:]])
        f.write(">" + barcode_sequence + "\n" + variant_sequence + "\n")
"""

"""
bcftools csv output info fields

##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">
##INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">
##INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">
##INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">
##INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##INFO=<ID=BCSQ,Number=.,Type=String,Description="Haplotype-aware consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change">
##FORMAT=<ID=BCSQ,Number=.,Type=Integer,Description="Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \"bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]'\" to translate.">

Example output:
gp17_orf	1642	.	C	T	30.4183	.	DP=1;SGB=-0.379885;MQ0F=0;AC=1;AN=1;DP4=0,0,1,0;MQ=60;BCSQ=missense|gp17|gp17|protein_coding|+|548P>548S|1642C>T+1644A>T	GT:PL:BCSQ	1:60,0:1
gp17_orf	1644	.	A	T	30.4183	.	DP=1;SGB=-0.379885;MQ0F=0;AC=1;AN=1;DP4=0,0,1,0;MQ=60;BCSQ=@1642	GT:PL:BCSQ	1:60,0:1

Note for BCSQ format from docs: https://samtools.github.io/bcftools/howtos/csq-calling.html
"The last three fields are omitted when empty. Consequences of compound variants which span multiple sites are printed in one record only,
the remaining records link to it by '@position'. The consequence can start with the asterisk '*' prefix indicating a consequence downstream
from a stop. For more details and examples please see the manual page and the split-vep plugin."


"""
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
    value_dict = {key: None for key in bcsq_keys}
    values = info_str.split("|")
    for value in values:
        if value.startswith("@"):
            value_dict["linked_record"] = value
        else:
            value_dict["linked_record"] = None
            # If the consequences are empty, the last three fields are omitted.
            # Otherwise there will be 7 values
            if len(values) == 7:
                value_dict.update(dict(zip(bcsq_keys, values)))
            if len(values) == 10:
                # 3 prime or 5 prime UTR consequences. Skip these.
                value_dict = {key: None for key in bcsq_keys}
            elif len(values) == 4:
                value_dict.update(dict(zip(bcsq_keys[:4], values)))
            else:
                Warning(f"Unexpected number of values in BCSQ field: {values}")
                value_dict = {key: None for key in bcsq_keys}
                # raise ValueError(f"Unexpected number of values in BCSQ field: {values}")
    return value_dict


with open(output_consensus_file, "w", encoding="UTF-8") as f:
    headers = ["barcode"] + vcf_keys + info_keys + bcsq_keys
    consensus_tsv_writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
    consensus_tsv_writer.writeheader()
    for barcode_bcf in input_barcode_names:
        barcode_sequence = barcode_bcf.split("/")[-1].split("_")[0]

        # Input sequences provided as "results/consensus/{{sample}}/{barcode}_csq.bcf", barcode=barcodes)
        variant_file = os.path.abspath(barcode_bcf)
        for variant in VCF(variant_file):
            # Parse BCF variant fields
            record_dict = {"barcode": barcode_sequence}
            record_dict["CHROM"] = variant.CHROM
            record_dict["POS"] = variant.POS
            record_dict["ID"] = variant.ID
            record_dict["REF"] = variant.REF
            record_dict["ALT"] = ",".join(variant.ALT)
            record_dict["QUAL"] = variant.QUAL
            record_dict["FILTER"] = variant.FILTER
            # Parse info fields
            info = variant.INFO
            for key in info_keys:
                field = info.get(key)
                if key == "BCSQ":
                    record_dict["BCSQ"] = field
                    record_dict.update(parse_bcsq(field))
                else:
                    record_dict[key] = field
            consensus_tsv_writer.writerow(record_dict)
