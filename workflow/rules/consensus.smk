rule consequence:
    """Annotate variants using bcftools csq."""
    input:
        "results/consensus/{barcode}.bcf",
    output:
        "results/consensus/{barcode}.consequences.tsv",
    shell:
        "bcftools csq -f {reference} -g {gff} {input} | bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' > {output}"


rule mpileup:
    """Call variants using bcftools mpileup."""
    input:
        "results/mapped/{barcode}.bam",
    output:
        "results/consensus/{barcode}.bcf",
    shell:
        "bcftools mpileup -d 5000 -Ou -f {reference} {input} | bcftools call -vm -Ob --ploidy 1 -o {output}"


rule consensus:
    """Generate consensus sequences from matching barcode sequences using samtools."""
    input:
        "results/mapped/{barcode}.bam",
    output:
        "results/consensus/{barcode}.fasta",
    log:
        "logs/consensus/{barcode}.log",
    threads: 16
    shell:
        "samtools consensus -m simple -f fasta {input} > {output} 2> {log}"


rule match_barcodes:
    """Match barcodes to clustered sequences."""
    input:
        clusters="results/cutadapt/{sample}.barcodes.clusters.txt",
        variants="results/consensus/{sample}.variant_consequences.tsv",
    output:
        "results/{sample}.variants.tsv",
    script:
        "scripts/match_barcodes.py"


rule write_barcodes:
    """Parse barcode counts and filter for minimum counts."""
    input:
        "results/cutadapt/{sample}.barcodes.info.tsv",
    output:
        "results/cutadapt/{sample}.filter.barcodes.txt",
    params:
        min_count=40,
    script:
        "scripts/write_barcodes.py"


rule make_cluster_fastas:
    """Generate individual fasta files containing clustered sequences."""
    input:
        "results/cutadapt/{sample}.barcodes.clusters.txt",
    output:
        "results/clusters/{sample}.fasta",
    params:
        min_size=1000,
    script:
        "scripts/make_cluster_fastas.py"


rule write_sequences:
    """Match barcodes to clustered sequences."""
    input:
        info="results/cutadapt/{sample}.barcodes.info.tsv",
        barcode_clusters="results/cutadapt/{sample}.barcodes.clusters.txt",
    output:
        "results/clusters/{sample}.clusters.csv",
    script:
        "scripts/write_sequences.py"
