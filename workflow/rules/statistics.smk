from os.path import dirname, basename


rule get_ids:
    input:
        "results/processed/genotypes.fam"
    output:
        temp("results/statistics/data/id.txt")
    shell:
        "cut -d' ' -f2 {input} > {output}"


rule make_phenotype_file:
    input:
        config["metadata"],
        id="results/statistics/data/id.txt",
    output:
        "results/statistics/data/phenotype.txt"
    params:
        trait=config["trait"],
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/make_phenotype_file.log"
    script:
        "../scripts/make_phenotype_file.py"


if config["pcs"]["enabled"]:
    rule get_pcs:
        input:
            bed="results/processed/genotypes.bed",
            bim="results/processed/genotypes.bim",
            fam="results/processed/genotypes.fam",
        output:
            "results/statistics/pca/pcs.eigenvec",
            "results/statistics/pca/pcs.eigenval",
        params:
            bfile=lambda wildcards, input: f"--bfile {input.bed.removesuffix('.bed')}",
            out=lambda wildcards, output: f"--out {output[0].removesuffix('.eigenvec')}",
            pca="--pca {}".format(config["pcs"]["count"]),
            extra=config["params"]["plink"]
        shadow:
            "minimal"
        threads:
            32
        conda:
            "../envs/plink.yaml"
        log:
            "logs/get_pcs.log"
        shell:
            "plink --threads {threads} {params.pca} {params.extra} {params.bfile} "
            "{params.out} >{log} 2>&1"


USE_COVARIATES = config["covariates"]["enabled"] or config["pcs"]["enabled"]


def _covariate_inputs(wildcards):
    d = {
        "metadata": config["metadata"],
        "id": rules.get_ids.output[0],
    }
    if config["pcs"]["enabled"]:
        d["pca"] = rules.get_pcs.output[0]
    return d


if USE_COVARIATES:
    rule make_covariates_file:
        input:
            unpack(_covariate_inputs),
        output:
            "results/statistics/data/covariates.txt"
        params:
            covar=config["covariates"]["columns"]
        conda:
            "../envs/pandas.yaml"
        log:
            "logs/make_covariates_file.log"
        script:
            "../scripts/make_covariates_file.py"


rule estimate_relatedness_matrix:
    input:
        bed="results/processed/genotypes.bed",
        phe="results/statistics/data/phenotype.txt",
    output:
        temp("results/statistics/grm.cXX.txt"),
        log="results/statistics/grm.log.txt"
    params:
        bfile=lambda wildcards, input: f"-bfile {input["bed"].removesuffix('.bed')}",
        outdir=lambda wildcards, output: f"-outdir {dirname(output[0])}",
        o=lambda wildcards, output: f"-o {basename(output[0]).removesuffix('.cXX.txt')}",
        extra=config["params"]["gemma"]["grm"]
    log:
        "logs/estimate_relatedness_matrix.log"
    benchmark:
        "benchmarks/estimate_relatedness_matrix.log"
    conda:
        "../envs/gemma.yaml"
    shell:
        "gemma {params.extra} {params.bfile} {params.outdir} {params.o} "
        "-p {input.phe} >{log} 2>&1"


def _lmm_inputs(wildcards):
    d = {
        "bed": "results/processed/genotypes.bed",
        "cxx": "results/statistics/grm.cXX.txt",
        "phe": "results/statistics/data/phenotype.txt",
    }
    if USE_COVARIATES:
        d["cov"] = "results/statistics/data/covariates.txt"
    return d


rule fit_inear_mixed_model:
    input:
        unpack(_lmm_inputs),
    output:
        temp("results/statistics/lmm.assoc.txt"),
        log="results/statistics/lmm.log.txt",
    params:
        bfile=lambda wildcards, input: f"-bfile {input["bed"].removesuffix('.bed')}",
        outdir=lambda wildcards, output: f"-outdir {dirname(output[0])}",
        o=lambda wildcards, output: f"-o {basename(output[0]).removesuffix('.assoc.txt')}",
        extra=config["params"]["gemma"]["lmm"],
        cov=lambda wildcards, input: f"-c {input.cov}" if input.get("cov") else "",
    log:
        "logs/fit_linear_mixed_model.log"
    benchmark:
        "benchmarks/fit_linear_mixed_model.log"
    conda:
        "../envs/gemma.yaml"
    shell:
        "gemma {params.extra} {params.bfile} {params.outdir} {params.o} "
        "-k {input.cxx} -p {input.phe} {params.cov} >{log} 2>&1"


rule compress_assoc:
    input:
        "results/statistics/lmm.assoc.txt"
    output:
        "results/statistics/lmm.assoc.txt.gz"
    threads:
        64
    conda:
        "../envs/pigz.yaml"
    shell:
        "pigz -p {threads} -c {input} > {output}"
