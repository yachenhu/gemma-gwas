ANN_CONFIG = config["ref"].get("annotation", "ensembl")

if ANN_CONFIG == "ncbi":
    rule get_ncbi_annotation:
        output:
            "resources/ncbi_annotation.gtf.gz",
        params:
            url=config["ref"].get("ncbi_gtf_url", ""),
        log:
            "logs/get_ncbi_annotation.log"
        cache: True
        shell:
            "wget -q --show-progress -O {output} {params.url} 2>{log}"

    _gtf_input  = "resources/ncbi_annotation.gtf.gz"
    _gtf_gene_key = "gene"
else:
    _gtf_input   = "resources/annotation.gtf"
    _gtf_gene_key = "gene_name"


rule plot_manhattan_qq:
    input:
        "results/statistics/lmm.assoc.txt.gz",
    output:
        "results/figures/manhattan.png",
        "results/figures/qq.png",
    params:
        top=100,
    log:
        "logs/plot_manhattan_qq.log"
    conda:
        "../envs/gwaslab.yaml"
    shell:
        "python workflow/scripts/plot_manhattan_qq.py "
        "{input} {output[0]} {output[1]} {params.top} >{log} 2>&1"


rule plot_regional:
    input:
        gwas="results/statistics/lmm.assoc.txt.gz",
        vcf="results/processed/genotypes.filtered.vcf.gz",
        vcf_idx="results/processed/genotypes.filtered.vcf.gz.tbi",
        gtf=_gtf_input,
    output:
        touch("results/figures/regional/.done"),
    params:
        n=5,
        window=500000,
        gene_key=_gtf_gene_key,
    log:
        "logs/plot_regional.log"
    conda:
        "../envs/gwaslab.yaml"
    shell:
        "mkdir -p results/figures/regional && "
        "python workflow/scripts/plot_regional.py "
        "{input.gwas} {input.vcf} {input.gtf} {params.gene_key} "
        "{ANN_CONFIG} "
        "{params.n} {params.window} "
        "results/figures/regional >{log} 2>&1"
