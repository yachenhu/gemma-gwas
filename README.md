# gwas-gemma

A Snakemake-based pipeline for genome-wide association studies (GWAS) using GEMMA (Genome-wide Efficient Mixed Model Association).

## Overview

This pipeline performs comprehensive genome-wide association analysis using GEMMA, implementing a mixed-model approach to account for population structure and relatedness among samples. The pipeline includes data preprocessing, statistical analysis, annotation, and visualization components.

## Features

- **Genome-wide association analysis** using GEMMA's efficient mixed model algorithms
- **Population stratification correction** using principal components
- **Covariate adjustment** for confounding factors like sex
- **VEP-based annotation** of genetic variants
- **Manhattan and Q-Q plot generation** for result visualization
- **GRM (Genetic Relationship Matrix)** calculation for relatedness modeling

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)

## Dependencies

The pipeline uses several bioinformatics and statistical tools managed through Snakemake environments:
- PLINK - data preprocessing and management
- GEMMA - mixed-model association testing
- VEP (Variant Effect Predictor) - variant annotation
- BCFTools - VCF processing
- R packages for plotting and analysis

## Input Data

The pipeline requires the following input files in the `data/` directory:

- `all.imputed.vcf.gz` - Genotype data in VCF format (configured in `config.yaml`). Sample names must already be in their final form.
- `metadata.tsv` - Sample metadata table with columns for sample identifiers, sex, sire, dam, and traits

Note: The example data included in this repository may not include the VCF file due to size constraints. You will need to provide your own genotype data in VCF format.

## Configuration

Pipeline parameters are configured in `config/config.yaml`. Key parameters include:

- `vcf` - Path to input VCF file
- `metadata` - Sample metadata table
- `ref` - Reference genome information (species, build, release)
- `sample` - Column names for sex, sire, dam in the metadata table
- `trait` - Name of the primary trait to analyze
- `covariates` - List of covariates to include in the model
- `pcs` - Number of principal components to use for population stratification
- `params` - Additional parameters for PLINK, GEMMA, VEP, and BCFTools

## Output

The pipeline generates the following output:

- `results/figures/manhattan.html` - Interactive Manhattan plot
- `results/figures/qq.html` - Quantile-quantile plot
- Various intermediate files in `results/analyzed/` and `results/annotated/`

## Usage

To run the pipeline:

```bash
# Install conda environment with snakemake
conda install -c conda-forge -c bioconda snakemake

# Execute the pipeline (run from the project root, not workflow/)
snakemake -s workflow/Snakefile --use-conda --cores N
```

Replace `N` with the number of CPU cores to use.

## Pipeline Workflow

The pipeline consists of several modules:

1. **Common** - Shared utility functions
2. **Reference** - Reference genome handling
3. **Processing** - Data preprocessing and QC
4. **Statistics** - GEMMA-based association testing
5. **Annotation** - Variant annotation using VEP
6. **Plotting** - Manhattan and Q-Q plot generation

## License

This project is licensed under the terms found in the LICENSE file.
