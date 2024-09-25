chk_experiments = ["experiment1", "experiment2"]
chk_samples = ["sample1", "sample2"]


def get_input_test(wildcards):
    """Generate the input files for the dummy rule all.
    This is necessary to allow optional pipeline outputs."""

    input_list = []

    experiments = ["experiment1", "experiment2"]
    samples = ["sample1", "sample2"]

    input_list.extend(
        expand(
            "my_directory2/{experiment_name}/{samples}.txt",
            experiment_name=chk_experiments,
            samples=chk_samples,
        )
    )
    print(input_list)

    return input_list


# a target rule to define the desired final output
rule all:
    input:
        get_input_test,


# the checkpoint that shall trigger re-evaluation of the DAG
# an number of file is created in a defined directory
checkpoint somestep:
    output:
        directory("my_directory/"),
    run:
        shell("mkdir my_directory")
        shell("cd my_directory")
        for f in chk_samples:
            for e in chk_experiments:
                shell("touch my_directory/{e}_{f}.txt")


# input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.somestep.get(**wildcards).output[0]
    print(checkpoint_output)
    output_files = glob_wildcards(
        os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
    )
    print(
        glob_wildcards(
            os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
        )
    )
    print(
        set(
            glob_wildcards(
                os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
            ).experiment_name
        )
    )
    print(
        expand(
            "my_directory/{experiment_name}_{samples}.txt",
            experiment_name=set(
                glob_wildcards(
                    os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
                ).experiment_name
            ),
            samples=set(
                glob_wildcards(
                    os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
                ).samples
            ),
        )
    )
    return expand(
        "my_directory/{experiment_name}_{samples}.txt",
        experiment_name=glob_wildcards(
            os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
        ).experiment_name,
        samples=glob_wildcards(
            os.path.join(checkpoint_output, "{experiment_name}_{samples}.txt")
        ).samples,
    )


rule aggregate:
    input:
        aggregate_input,
    output:
        "my_directory2/{experiment_name}/{samples}.txt",
    shell:
        "cat {input} > {output}"
