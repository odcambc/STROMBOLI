rule cutadapt:
    """Detect barcode sequences from nanopore sequencing data using cutadapt."""
    input:
        "reads/{sample}.fastq.gz"
    output:
        "results/cutadapt/{sample}.barcodes.txt",
        "results/cutadapt/{sample}.barcodes.info.tsv"
    params:
        adapters=expand("-g {flanking_sequences}", flanking_sequences=config["flanking_sequence"])
        extra="--discard-untrimmed --info-file=results/cutadapt/{sample}.info.tsv --json=results/cutadapt/{sample}.cutadapt.json",
    benchmark:
        "benchmarks/cutadapt/{sample}.benchmark.txt"
    log:
        "logs/cutadapt/{sample}.log",
    threads: 8
    wrapper:
        "v3.6.0/bio/cutadapt/se"

rule starcode:
    """Cluster detected barcode sequences using starcode."""
    input:
        "results/cutadapt/{sample}.barcodes.txt"
    output:
        "results/cutadapt/{sample}.barcodes.clusters.txt"
    params:
        extra="--print-clusters"
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
        "results/cutadapt/{sample}.barcodes.info.tsv"
    output:
        "results/cutadapt/{sample}.matches.tsv"
    shell:
        "awk '$2 != /-1/ {print $0}' {input} |  awk '$8 ~ /1;2/ {print $1, $5} > {output}"

rule consensus:
    """Generate consensus sequences from matching barcode sequences using samtools."""
    input:
        "results/mapped/{barcode}.bam"
    output:
        "results/consensus/{barcode}.fasta"
    log:
        "logs/consensus/{sample}.log",
    threads: 16
    shell:
        "samtools consensus -m simple -f fasta {input} > {output} 2> {log}"

rule mpileup:
    """Call variants using bcftools mpileup."""
    input:
        "results/mapped/{barcode}.bam"
    output:
        "results/consensus/{barcode}.bcf"
    shell:
        "bcftools mpileup -d 5000 -Ou -f {reference} {input} | bcftools call -vm -Ob --ploidy 1 -o {output}"

rule consequence:
    """Annotate variants using bcftools csq."""
    input:
        "results/consensus/{barcode}.bcf"
    output:
        "results/consensus/{barcode}.consequences.tsv"
    shell:
        "bcftools csq -f {reference} -g {gff} {input} | bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' > {output}"

rule minimap2:
    """Map consensus sequences to reference using minimap2."""
    input:
        consensus="results/consensus/{barcode}.fasta",
        reference="reference.fasta"
    output:
        "results/mapped/{barcode}.sam"
    log:
        "logs/minimap2/{barcode}.log",
    threads: 16
    shell:
        "minimap2 -ax map-ont {input.reference} {input.consensus} -o {output} 2> {log}"

rule make_bam:
    """Convert sam to bam using samtools."""
    input:
        "results/mapped/{barcode}.sam"
    output:
        "results/mapped/{barcode}.bam"
    shell:
        "samtools sort {input} > {output}"

rule index_bam:
    """Index bam files using samtools."""
    input:
        "results/mapped/{barcode}.bam"
    output:
        "results/mapped/{barcode}.bam.bai"
    shell:
        "samtools index {input}"

rule make_cluster_fastas:
    """Generate individual fasta files containing clustered sequences."""
    input:
        "results/cutadapt/{sample}.barcodes.clusters.txt"
    output:
        "results/clusters/{barcode}.fasta"
    params:
        min_size = 1000
    script:
        "scripts/make_cluster_fastas.py""

rule match_barcodes:
    """Match barcodes to clustered sequences."""
    input:
        clusters = "results/cutadapt/{sample}.barcodes.clusters.txt",
        variants = "results/consensus/{sample}.variant_consequences.tsv"
    output:
        "results/{sample}.variants.tsv"
    script:
        "scripts/match_barcodes.py"


rule write_barcodes:
    """Parse barcode counts and filter for minimum counts."""
    input:
        "results/cutadapt/{sample}.barcodes.info.tsv"
    output:
        "results/cutadapt/{sample}.filter.barcodes.txt"
    params:
        min_count = 40
    script:
        "scripts/write_barcodes.py"

rule write_sequences:
    """Match barcodes to clustered sequences."""
    input:
        info = "results/cutadapt/{sample}.barcodes.info.tsv"
        barcode_clusters = "results/cutadapt/{sample}.barcodes.clusters.txt"
    output:
        "results/clusters/{sample}.clusters.csv"
    script:
        "scripts/write_sequences.py"