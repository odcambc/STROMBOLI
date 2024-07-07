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


checkpoint make_cluster_fastas:
    """Generate individual fasta files containing clustered sequences."""
    input:
        "results/clusters/{sample}.clusters.csv",
    output:
        directory("results/clusters/barcodes/{sample}/"),
    params:
        min_cluster_size=config["min_cluster_size"],
    script:
        "scripts/make_cluster_fastas.py"


rule minimap2_map_clusters:
    """Map clustered sequences to reference using minimap2."""
    input:
        clusters="results/clusters/barcodes/{sample}/{barcode}.fastq",
        reference=reference_file,
    output:
        "results/clusters/{sample}/mapped/{barcode}.sam",
    log:
        "logs/minimap2/{sample}/{barcode}.log",
    threads: 16
    shell:
        "minimap2 -ax map-ont {input.reference} {input.clusters} -o {output} 2> {log}"


rule make_bam:
    """Convert sam to bam using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.sam",
    output:
        "results/clusters/{sample}/mapped/{barcode}.bam",
    shell:
        "samtools sort {input} > {output}"


rule index_bam:
    """Index bam files using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.bam",
    output:
        "results/clusters/{sample}/mapped/{barcode}.bam.bai",
    shell:
        "samtools index {input}"


rule consensus:
    """Generate consensus sequences from matching barcode sequences using samtools."""
    input:
        "results/clusters/{sample}/mapped/{barcode}.bam",
    output:
        "results/consensus/{sample}/{barcode}.fasta",
    params:
        use_qual=use_qual,
        call_fract=call_fract,
    log:
        "logs/consensus/{sample}/{barcode}.log",
    threads: 16
    shell:
        "samtools consensus {params.use_qual} {params.call_fract} "
        "-m simple -f fasta {input} > {output} 2> {log}"


rule minimap2_map_consensus_fasta:
    """Map consensus fasta file to reference with minimap2."""
    input:
        consensus="results/consensus/{sample}/{barcode}.fasta",
        reference="references/gp17_orf.fasta",
    output:
        "results/consensus/{sample}/{barcode}_consensus.sam",
    log:
        "logs/minimap2/{sample}/{barcode}_consensus.log",
    threads: 16
    shell:
        "minimap2 -ax map-ont {input.reference} {input.consensus} -o {output} 2> {log}"


rule make_bam_consensus:
    """Convert sam to bam using samtools."""
    input:
        "results/consensus/{sample}/{barcode}_consensus.sam",
    output:
        "results/consensus/{sample}/{barcode}_consensus.bam",
    shell:
        "samtools sort {input} > {output}"


rule index_bam_consensus:
    """Index bam files using samtools."""
    input:
        "results/consensus/{sample}/{barcode}_consensus.bam",
    output:
        "results/consensus/{sample}/{barcode}_consensus.bam.bai",
    shell:
        "samtools index {input}"


rule mpileup:
    """Call variants using bcftools mpileup."""
    input:
        bam="results/consensus/{sample}/{barcode}_consensus.bam",
        reference="references/gp17_orf.fasta",
    output:
        "results/consensus/{sample}/{barcode}_consensus.bcf",
    shell:
        "bcftools mpileup -d 5000 -Ou -f {input.reference} {input.bam} | bcftools call -vm -Ob --ploidy 1 -o {output}"


rule consequence:
    """Annotate variants using bcftools csq.
    Format: 
    Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change
    """
    input:
        bcf="results/consensus/{sample}/{barcode}_consensus.bcf",
        reference="references/gp17_orf.fasta",
    output:
        "results/consensus/{sample}/{barcode}_csq.bcf",
    params:
        gff="references/gp17_orf.gff",
    shell:
        "bcftools csq -f {input.reference} -g {params.gff} {input.bcf} --verbose 2 > {output}"


def get_barcode_names(wildcards):
    checkpoint_output = checkpoints.make_cluster_fastas.get(**wildcards).output[0]
    barcodes = glob_wildcards(
        os.path.join(checkpoint_output, "{barcode}.fastq")
    ).barcode

    return expand("results/consensus/{{sample}}/{barcode}_csq.bcf", barcode=barcodes)


rule write_consensus_variants:
    """Write all consequences to single file."""
    input:
        get_barcode_names,
    output:
        "results/consensus/{sample}_consensus_variants.tsv",
    script:
        "scripts/write_consensus_variants.py"


# ----- To fix or write -----
rule create_genbank:
    """Create a genbank file from a reference fasta."""
    input:
        reference_file,
    output:
        "references/{reference_name}.gb",
    script:
        "scripts/create_genbank.py"


rule create_gff3:
    """Create a gff3 file from a genbank file."""
    input:
        "references/{reference_name}.gb",
    output:
        "references/{reference_name}.gff",
    script:
        "scripts/create_gff3.py"


rule match_barcodes:
    """Match barcodes to clustered sequences."""
    input:
        clusters="results/starcode/{sample}.barcodes.clusters.txt",
        variants="results/consensus/{sample}.consequences.tsv",
    output:
        "results/{sample}.variants.tsv",
    script:
        "scripts/match_barcodes.py"


rule write_sequences:
    """Match barcodes to clustered sequences."""
    input:
        info="results/cutadapt/{sample}.barcodes.info.tsv",
        barcode_clusters="results/starcode/{sample}.barcodes.clusters.txt",
    output:
        "results/clusters/{sample}.clusters.csv",
    script:
        "scripts/write_sequences.py"
