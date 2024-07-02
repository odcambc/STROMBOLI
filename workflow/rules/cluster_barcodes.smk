rule cutadapt:
    """Detect barcode sequences from nanopore sequencing data using cutadapt."""
    input:
        get_file_from_sample,
    output:
        "results/cutadapt/{sample}.barcodes.info.tsv",
        barcodes_fastq="results/cutadapt/{sample}.barcodes.fastq.gz",
    params:
        adapters=expand(
            "-g {flanking_sequences}", flanking_sequences=config["flanking_sequence"]
        ),
        extra="--discard-untrimmed --info-file=results/cutadapt/{sample}.barcodes.info.tsv --json=results/cutadapt/{sample}.cutadapt.json",
    benchmark:
        "benchmarks/cutadapt/{sample}.benchmark.txt"
    log:
        "logs/cutadapt/{sample}.log",
    threads: 8
    shell:
        "cutadapt -j {threads} "
        "{params.extra} {params.adapters} "
        "-o {output.barcodes_fastq} {input}"


rule starcode:
    """Cluster detected barcode sequences using starcode."""
    input:
        "results/cutadapt/{sample}.barcodes.txt",
    output:
        "results/cutadapt/{sample}.barcodes.clusters.txt",
    params:
        extra="--print-clusters",
    log:
        "logs/starcode/{sample}.log",
    threads: 16
    shell:
        "starcode -i {input} "
        "--print-clusters "
        "-t {threads} "
        "-o {output} 2> {log}"


rule filter_awk:
    """Filter output of cutadapt to find matching sequences.
    This rule finds reads where the barcode was detected and outputs
    the upstream sequence along with the identified barcode."""
    input:
        "results/cutadapt/{sample}.barcodes.info.tsv",
    output:
        "results/cutadapt/{sample}.matches.tsv",
    shell:
        "awk '$2 != /-1/ {print $0}' {input} |  awk '$8 ~ /1;2/ {print $1, $5} > {output}"
