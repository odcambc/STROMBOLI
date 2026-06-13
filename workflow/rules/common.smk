def get_file_from_sample(wildcards):
    """Maps from the sequencing output file names to the sample names defined in the experiment CSV.
    This tries to determine the file name from the sample name, which may or may not be the prefix.
    """
    filename = experiments.loc[
        experiments["sample"] == wildcards.sample, "file"
    ].squeeze()

    if filename is None or type(filename) is pd.Series and filename.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in experiment file.")

    full_name = config["data_dir"] + "/" + filename

    if (
        filename.endswith(".fastq.gz")
        or filename.endswith(".fastq")
        or filename.endswith(".fq")
        or filename.endswith(".fq.gz")
        or filename.endswith(".sam")
        or filename.endswith(".bam")
    ):
        pass
    else:
        full_name = full_name + ".fastq.gz"
    return full_name


def get_ref(wildcards):
    """Removes file suffix from reference fasta file."""
    prefix = config["reference"].split(".fasta")[0]
    return prefix


def get_input(wildcards):
    """Final targets for `rule all`: per-sample barcode->variant table and QC summary."""
    samples_list = list(samples)
    return expand(
        "results/{sample}.variants.tsv", sample=samples_list
    ) + expand("results/qc/{sample}.stromboli_qc.json", sample=samples_list)


# Validate config and experiment files
validate(config, "../schemas/config.schema.yaml")

experiments = (
    pd.read_csv(config["experiment_file"], header=0)
    .dropna(how="all")
    .set_index("sample", drop=False, verify_integrity=True)
)

validate(experiments, "../schemas/experiments.schema.yaml")

# Set variables from config and experiment files
experiment = config["experiment"]
samples = experiments["sample"]
reference_name = get_ref(config["reference"])

reference_file = os.path.join(os.path.abspath(config["ref_dir"]), config["reference"])

if config["use_qual"]:
    use_qual = "-q "
else:
    use_qual = ""

if config["call_fract"]:
    call_fract = "-c " + str(config["call_fract"])
else:
    call_fract = ""
