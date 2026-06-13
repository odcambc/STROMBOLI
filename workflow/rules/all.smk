# Barcodes are nucleotide strings; samples may contain dots/underscores/dashes.
# Constraining the wildcards keeps path patterns unambiguous and makes the
# checkpoint's glob_wildcards (below) safe against stray files in the cluster dir.
wildcard_constraints:
    barcode="[ACGT]+",
    sample="[A-Za-z0-9_.-]+",


rule cutadapt:
    """Detect barcode sequences from nanopore sequencing data using cutadapt."""
    input:
        get_file_from_sample,
    output:
        temp("results/cutadapt/{sample}.barcodes.info.tsv"),
        barcodes_fastq="results/cutadapt/{sample}.barcodes.fastq.gz",
        json="results/cutadapt/{sample}.cutadapt.json",
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
        "-o {output.barcodes_fastq} {input} 1> {log}"


rule filter_awk:
    """Filter output of cutadapt to find matching sequences.
    This rule finds reads where the barcode was detected (the second linked
    adapter, name "1;2") and outputs the identified barcode sequence."""
    input:
        "results/cutadapt/{sample}.barcodes.info.tsv",
    output:
        "results/cutadapt/{sample}.barcodes.matches.txt",
    params:
        max_barcode_length=config["max_barcode_length"],
    shell:
        "awk '($2 !~ /-1/ && $8 ~ /1;2/ && length($5) < {params.max_barcode_length}) {{print $5}}' {input} > {output}"


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


rule write_sequences:
    """Group the insert reads belonging to each barcode cluster.
    Carries the per-read qualities (q_l) alongside the sequence so that real
    FASTQ records can be written for the consensus step."""
    input:
        info="results/cutadapt/{sample}.barcodes.info.tsv",
        barcode_clusters="results/starcode/{sample}.barcodes.clusters.txt",
    output:
        clusters="results/clusters/{sample}.clusters.tsv",
        qc="results/clusters/{sample}.cluster_qc.tsv",
    log:
        "logs/starcode/{sample}_write_sequences.log",
    script:
        "scripts/write_sequences.py"


checkpoint make_cluster_fastas:
    """Generate individual FASTQ files containing the reads for each cluster."""
    input:
        "results/clusters/{sample}.clusters.tsv",
    output:
        directory("results/clusters/barcodes/{sample}/"),
    params:
        min_cluster_size=config["min_cluster_size"],
    script:
        "scripts/make_cluster_fastas.py"


# ----- Per-barcode scatter (fused) -----
# Each barcode's entire chain runs as ONE rule in a single shell: map, (consensus for
# double), call, filter, and annotate are piped together. Versus a rule per step this
# removes per-job overhead and the intermediate Snakemake temp BAM/BCF writes; the only
# materialized intermediate is the sorted BAM that bcftools mpileup requires as a file.
# The tools still run, so the saving is orchestration + I/O, not process count. Single-
# threaded (tiny inputs) and grouped so the scatter parallelizes without scheduler churn.
#
# calling_mode (config):
#   "single_qc" : map -> mpileup on the cluster pileup -> filter by ALT allele fraction
#                 -> annotate. Keeps read depth; pair with a min_cluster_size floor; see
#                 tools/fdr_estimator.py for choosing the floor and AF threshold.
#   "double"    : map -> consensus -> re-map -> mpileup (depth-1) -> annotate. The
#                 conservative, depth-agnostic default.
# Both emit results/consensus/{sample}/{barcode}_csq.bcf.
_GFF = expand("references/{reference_name}.gff", reference_name=reference_name)


if config["calling_mode"] == "single_qc":

    rule call_barcode:
        """SINGLE+QC, fused: map|sort > bam; mpileup|call|norm|view(AF)|csq."""
        input:
            clusters="results/clusters/barcodes/{sample}/{barcode}.fastq",
            reference=reference_file,
            gff=_GFF,
        output:
            temp("results/consensus/{sample}/{barcode}_csq.bcf"),
        params:
            max_depth=config["mpileup_max_depth"],
            # Lower bound is clash_mixed_af (not qc_min_af): keep intermediate-AF
            # variants so write_consensus_variants can flag mixed barcodes. The
            # qc_min_af confident-call threshold is applied downstream in match_barcodes.
            filter_expr=(
                "FMT/AD[0:1] >= {ar} && "
                "FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1]) >= {af}"
            ).format(ar=config["qc_min_alt_reads"], af=config["clash_mixed_af"]),
        log:
            "logs/consensus/{sample}/{barcode}.log",
        group:
            "barcode"
        threads: 1
        shell:
            "bam={output}.sort.bam; trap 'rm -f $bam' EXIT; "
            "minimap2 -ax map-ont {input.reference} {input.clusters} 2> {log} "
            "| samtools sort -o $bam 2>> {log}; "
            "bcftools mpileup -d {params.max_depth} -a AD -Ou -f {input.reference} $bam 2>> {log} "
            "| bcftools call -vm --ploidy 1 -Ov 2>> {log} "
            "| bcftools norm -m- -f {input.reference} -Ov 2>> {log} "
            "| bcftools view -i '{params.filter_expr}' -Ou 2>> {log} "
            "| bcftools csq -f {input.reference} -g {input.gff} -Ob -o {output} --verbose 2 - 2>> {log}"

else:  # calling_mode == "double"

    rule call_barcode:
        """DOUBLE, fused: map|sort; consensus; re-map|sort; mpileup|call|csq."""
        input:
            clusters="results/clusters/barcodes/{sample}/{barcode}.fastq",
            reference=reference_file,
            gff=_GFF,
        output:
            temp("results/consensus/{sample}/{barcode}_csq.bcf"),
        params:
            use_qual=use_qual,
            call_fract=call_fract,
            max_depth=config["mpileup_max_depth"],
        log:
            "logs/consensus/{sample}/{barcode}.log",
        group:
            "barcode"
        threads: 1
        shell:
            "cbam={output}.cl.bam; cons={output}.cons.fa; sbam={output}.cons.bam; "
            "trap 'rm -f $cbam $cons $sbam' EXIT; "
            "minimap2 -ax map-ont {input.reference} {input.clusters} 2> {log} "
            "| samtools sort -o $cbam 2>> {log}; "
            "samtools consensus {params.use_qual} {params.call_fract} -m simple -f fasta $cbam > $cons 2>> {log}; "
            "minimap2 -ax map-ont {input.reference} $cons 2>> {log} "
            "| samtools sort -o $sbam 2>> {log}; "
            "bcftools mpileup -d {params.max_depth} -Ou -f {input.reference} $sbam 2>> {log} "
            "| bcftools call -vm --ploidy 1 -Ou 2>> {log} "
            "| bcftools csq -f {input.reference} -g {input.gff} -Ob -o {output} --verbose 2 - 2>> {log}"


def get_barcode_names(wildcards):
    """Gather the per-barcode csq outputs produced after the checkpoint runs."""
    checkpoint_output = checkpoints.make_cluster_fastas.get(**wildcards).output[0]
    barcodes = glob_wildcards(
        os.path.join(checkpoint_output, "{barcode}.fastq")
    ).barcode

    return expand("results/consensus/{{sample}}/{barcode}_csq.bcf", barcode=barcodes)


rule write_consensus_variants:
    """Aggregate all per-barcode consequences for a sample into one table."""
    input:
        get_barcode_names,
    output:
        "results/consensus/{sample}_consensus_variants.tsv",
    log:
        "logs/consensus/{sample}_write_consensus.log",
    script:
        "scripts/write_consensus_variants.py"


rule match_barcodes:
    """Final barcode -> variant mapping: expand each cluster's canonical barcode to
    all members, flag merged/mixed clashes, and (optionally) exclude them. Flagged
    barcodes and the reason are written to results/{sample}.flagged.tsv."""
    input:
        clusters="results/starcode/{sample}.barcodes.clusters.txt",
        variants="results/consensus/{sample}_consensus_variants.tsv",
        qc="results/clusters/{sample}.cluster_qc.tsv",
    output:
        mapping="results/{sample}.variants.tsv",
        flagged="results/{sample}.flagged.tsv",
    params:
        clash_merge_fraction=config["clash_merge_fraction"],
        clash_mixed_af=config["clash_mixed_af"],
        qc_min_af=config["qc_min_af"],
        exclude_clashes=config["exclude_clashes"],
    log:
        "logs/consensus/{sample}_match_barcodes.log",
    script:
        "scripts/match_barcodes.py"


rule qc_summary:
    """Aggregate the sample's outputs into one QC summary JSON for MultiQC
    (consumed by the multiqc-stromboli plugin)."""
    input:
        cutadapt="results/cutadapt/{sample}.cutadapt.json",
        cluster_qc="results/clusters/{sample}.cluster_qc.tsv",
        consensus_variants="results/consensus/{sample}_consensus_variants.tsv",
        variants="results/{sample}.variants.tsv",
        flagged="results/{sample}.flagged.tsv",
    output:
        "results/qc/{sample}.stromboli_qc.json",
    params:
        min_cluster_size=config["min_cluster_size"],
    log:
        "logs/qc/{sample}_qc_summary.log",
    script:
        "scripts/write_qc_summary.py"


# ----- Reference preparation -----
# These only run when the .gff / .gb are missing; the repo ships amplicon_ref.gff.


rule create_gff3:
    """Create a gff3 file from a genbank file.
    convert_genbank_to_gff3.py appends a ##FASTA section, which bcftools csq cannot
    parse, so we strip everything from the FASTA block onward."""
    input:
        "references/{reference_name}.gb",
    output:
        "references/{reference_name}.gff",
    log:
        "logs/ref/{reference_name}_create_gff3.log",
    shell:
        "convert_genbank_to_gff3.py -i {input} -o {output}.tmp 2> {log} && "
        "awk '/^##FASTA/{{exit}} /^>/{{exit}} {{print}}' {output}.tmp > {output} && "
        "rm -f {output}.tmp"


rule create_genbank:
    """Create a genbank file from a reference fasta."""
    input:
        reference_file,
    output:
        "references/{reference_name}.gb",
    params:
        orf=config["orf"],
        gene_name=lambda wc: config.get("gene_name", wc.reference_name),
    log:
        "logs/ref/{reference_name}_create_gb.log",
    script:
        "scripts/create_genbank.py"
