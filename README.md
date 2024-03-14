# Snakemake pipeline for nanopore-based barcode variant mapping

## Quick start

```bash
git clone https://github.com/odcambc/fahrcode
cd fahrcode
conda env create --file fahrcode_env.yaml
conda activate fahrcode
```
### Dependencies

#### Via conda (recommended)
The simplest way to handle dependencies is with [Conda](https://conda.io/docs/) and the provided environment file.

```bash
conda env create --file fahrcode_env.yaml
```

#### Manually

The following are the dependencies required to run the pipeline:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Samtools](http://www.htslib.org/)
* [BCFtools](https://samtools.github.io/bcftools/)
* [minimap2](https://github.com/lh3/minimap2)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [starcode](https://github.com/gui11aume/starcode)

## License

This is licensed under the MIT license. See the LICENSE file for details.

## Contributing

Contributions and feedback are welcome. Please submit an issue or pull request.

## Getting help

For any issues, please open an issue on the GitHub repository. For
questions or feedback, [email Chris](https://www.wcoyotelab.com/members/).
