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
        "results/cutadapt/{sample}.barcodes.matches.txt",
    output:
        "results/starcode/{sample}.barcodes.clusters.txt",
    params:
        extra="--print-clusters",
        distance=config["barcode_distance"],
    log:
        "logs/starcode/{sample}.log",
    threads: 16
    shell:
        "starcode -i {input} "
        "{params.extra} "
        "-d {params.distance} "
        "-t {threads} "
        "-o {output} "
        "2> {log}"


rule filter_awk:
    """Filter output of cutadapt to find matching sequences.
    This rule finds reads where the barcode was detected and outputs
    the upstream sequence along with the identified barcode."""
    input:
        "results/cutadapt/{sample}.barcodes.info.tsv",
    output:
        "results/cutadapt/{sample}.barcodes.matches.txt",
    shell:
        "awk '$2 !~ /-1/ {{print $0}}' {input} |  "
        "awk '$8 ~ /1;2/ {{print $5}}' > {output}"


rule write_sequences:
    """Match barcodes to clustered sequences."""
    input:
        info="results/cutadapt/{sample}.barcodes.info.tsv",
        barcode_clusters="results/starcode/{sample}.barcodes.clusters.txt",
    output:
        "results/clusters/{sample}.clusters.csv",
    script:
        "scripts/write_sequences.py"


checkpoint make_cluster_fastas:
    """Generate individual fasta files containing clustered sequences."""
    input:
        "results/clusters/{sample}.clusters.csv",
    output:
        directory("results/clusters/barcodes/{sample}/"),
    params:
        min_size=2,
    script:
        "scripts/make_cluster_fastas.py"


rule minimap2_cluster_aggregate:
    """Map clustered sequences to reference using minimap2."""
    input:
        clusters="results/clusters/barcodes/{sample}/{barcode}.fastq",
        reference="references/amplicon_ref.fasta",
    output:
        "results/clusters/{sample}/mapped/{barcode}.sam",
    log:
        "logs/minimap2/{sample}_{barcode}.log",
    threads: 16
    shell:
        "minimap2 -ax map-ont {input.reference} {input.clusters} -o {output} 2> {log}"


def get_sam_input(wildcards):
    checkpoint_output = checkpoints.make_cluster_fastas.get(**wildcards).output[0]
    sample = wildcards.sample
    barcodes = glob_wildcards(
        os.path.join(checkpoint_output, "{barcode}.fastq")
    ).barcode
    return expand(
        "results/clusters/{{sample}}/mapped/{barcode}.sam",
        sample=wildcards.sample,
        barcode=barcodes,
    )


rule make_list_of_files:
    """Create a list of files to be processed."""
    input:
        get_sam_input,
    output:
        "results/{sample}/barcodes_detected.txt",
    shell:
        "ls {input} > {output}"
