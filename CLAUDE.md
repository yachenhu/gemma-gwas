# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run the full pipeline (from project root, not workflow/)
snakemake -s workflow/Snakefile --use-conda --cores N

# Dry-run to see what will execute
snakemake -s workflow/Snakefile --use-conda -n

# Run a specific rule or output file
snakemake -s workflow/Snakefile --use-conda --cores N results/statistics/lmm.assoc.txt.gz
```

Replace `N` with the number of CPU cores. Conda environments are defined per-rule in `workflow/envs/` and activated automatically via `--use-conda`.

## Architecture

This is a Snakemake GWAS pipeline using GEMMA's mixed-model association. The main entry point is `workflow/Snakefile`, which loads config from `config/config.yaml` and includes rule modules in dependency order:

1. **`rules/ref.smk`** — Downloads reference genome, annotation, known variants, and VEP cache/plugins. Uses Snakemake wrapper v4.3.0 for most steps.
2. **`rules/processing.smk`** — Filters VCF variants by quality presets, annotates with known variant IDs, converts to PLINK binary format. Defines `FILTER_PRESETS` (relaxed/standard/strict) mapping to bcftools filter expressions.
3. **`rules/statistics.smk`** — Computes PCs via PLINK, builds phenotype and covariate files, estimates the genetic relationship matrix (GRM), and fits the linear mixed model with GEMMA.
4. **`rules/annotation.smk`** — Extracts top hits (`p_wald < 0.01`), runs VEP annotation, and parses VEP CSQ fields into a TSV via bcftools `+split-vep`.
5. **`rules/plotting.smk`** — Joins GWAS summary stats with annotations, generates interactive Manhattan and Q-Q plots using R's manhattanly.

## Key conventions

- **Run from project root**: The Snakefile is at `workflow/Snakefile` but config paths are relative to the project root. Always invoke `snakemake -s workflow/Snakefile` from the repo root.
- **Config-driven**: All parameters (trait, covariates, PCs, filter presets, tool flags) live in `config/config.yaml`. Scripts access them via the `snakemake` object.
- **Python scripts** (`workflow/scripts/*.py`): Receive `snakemake` as a global. Use `snakemake.input`, `snakemake.output`, `snakemake.params`, `snakemake.threads`, `snakemake.log_fmt_shell()`. Dependencies are in `pandas.yaml`, `csvtk.yaml`, `bcftools.yaml`, `pigz.yaml`.
- **R scripts**: Receive `snakemake` as an S4 object (e.g. `snakemake@input[[1]]`, `snakemake@params$top`).
- **PLINK parameters**: The global `--allow-extra-chr --chr-set 95` in config is for chicken (GRCg7b) genome. Chromosome count may need adjusting for other species.
- **GEMMA**: Uses `-gk 1` (centered GRM) and `-lmm 1` (Wald test). GWAS p-values are in column `p_wald`.

## Input data

- `data/all.imputed.vcf.gz` — genotype VCF (sample names must be final; no renaming step)
- `data/metadata.tsv` — tab-separated, first column is sample ID, with columns for sex, sire, dam, and trait(s)

## Output

- `results/figures/manhattan.html`, `results/figures/qq.html` — interactive plots (targets of `rule all`)
- Intermediate files under `results/processed/`, `results/statistics/`, `results/annotated/`
- Reference data cached under `resources/`
