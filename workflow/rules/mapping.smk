rule minimap2:
    """Map consensus sequences to reference using minimap2."""
    input:
        consensus="results/consensus/{barcode}.fasta",
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
