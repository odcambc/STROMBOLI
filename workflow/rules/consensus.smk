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


rule mpileup:
    """Call variants using bcftools mpileup."""
    input:
        "results/mapped/{barcode}.bam",
    output:
        "results/mapped/{barcode}.bcf",
    shell:
        "bcftools mpileup -d 5000 -Ou -f {reference} {input} | bcftools call -vm -Ob --ploidy 1 -o {output}"


rule consequence:
    """Annotate variants using bcftools csq."""
    input:
        "results/mapped/{barcode}.bcf",
    output:
        "results/consensus/{sample}.barcode_consequences.tsv",
    shell:
        "bcftools csq -f {reference} -g {gff} {input} | bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' > {output}"


# 5220  2024-02-19 20:21  ls | grep \.bcf$ | sed 's/\.bcf//' | xargs -I {} sh -c "bcftools csq -f ../../gp17_orf.fasta -g ../../gp17_orf.gff {}.bcf | bcftools query -f'[%CHROM\t%POS\t%SAMPLE\t%TBCSQ\n]' >> variant.consequences.23.tsv"


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


rule match_barcodes:
    """Match barcodes to clustered sequences."""
    input:
        clusters="results/starcode/{sample}.barcodes.clusters.txt",
        variants="results/consensus/{sample}.consequences.tsv",
    output:
        "results/{sample}.variants.tsv",
    script:
        "scripts/match_barcodes.py"
