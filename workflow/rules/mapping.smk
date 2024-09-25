rule minimap2_consensus:
    """Map consensus sequences to reference using minimap2."""
    input:
        consensus="results/consensus/{barcode}.fastq",
        reference="reference.fasta",
    output:
        "results/mapped/{barcode}.sam",
    log:
        "logs/minimap2/{barcode}.log",
    threads: 16
    shell:
        "minimap2 -ax map-ont {input.reference} {input.consensus} -o {output} 2> {log}"


rule make_bam:
    """Convert sam to bam using samtools."""
    input:
        "results/mapped/{barcode}.sam",
    output:
        "results/mapped/{barcode}.bam",
    shell:
        "samtools sort {input} > {output}"


rule index_bam:
    """Index bam files using samtools."""
    input:
        "results/mapped/{barcode}.bam",
    output:
        "results/mapped/{barcode}.bam.bai",
    shell:
        "samtools index {input}"


rule minimap2_cluster:
    """Map clustered sequences to reference using minimap2."""
    input:
        cluster="results/clusters/{sample}/{barcode}.fastq",
        reference="reference.fasta",
    output:
        "results/clusters/{sample}/mapped/{barcode}.sam2",
    log:
        "logs/minimap2/{sample}_{barcode}.log",
    threads: 16
    shell:
        "minimap2 -ax map-ont {input.reference} {input.cluster} -o {output} 2> {log}"


rule make_bam_cluster:
    """Convert sam to bam using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.sam",
    output:
        "results/clusters/{sample}/mapped/{barcode}.bam",
    shell:
        "samtools sort {input} > {output}"


rule index_bam_cluster:
    """Index bam files using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.bam",
    output:
        "results/clusters/{sample}/mapped/{barcode}.bam.bai",
    shell:
        "samtools index {input}"
