rule get_top_hits:
    input:
        "results/statistics/lmm.assoc.txt.gz",
    output:
        "results/annotated/top_hits.lst",
    params:
        snp="rs",
        p="p_wald",
        cutoff=1e-2,
    conda:
        "../envs/csvtk.yaml"
    shell:
        "set +o pipefail; "
        "csvtk -t sort -k {params.p}:n {input} "
        "| csvtk -t filter2 -f ' ${params.p} < {params.cutoff}' "
        "| csvtk -t cut -f {params.snp} "
        "| sed 1d "
        "> {output} "


rule extract_variants:
    input:
        "results/processed/all.updated.vcf.gz",
        ids="results/annotated/top_hits.lst"
    output:
        temp("results/annotated/top_hits.vcf.gz")
    params:
        extra=lambda wildcards, input: f"--include 'ID=@{input.ids}'"
    threads:
        16
    log:
        "logs/get_variants.py"
    wrapper:
        "v7.2.0/bio/bcftools/view"


rule vep_annotate:
    input:
        calls="results/annotated/top_hits.vcf.gz",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
        fasta="resources/genome.fasta",
        fai="resources/genome.fasta.fai"
    output:
        calls="results/annotated/top_hits.anno.vcf.gz",
        stats="results/annotated/top_hits.anno.stats.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["params"]["vep"]["plugins"],
        extra=config["params"]["vep"]["extra"],
        view_extra="-f .,PASS"
    log:
        "logs/vep/annotate.log",
    benchmark:
        "benchmarks/vep/annotate.json",
    threads: 
        64
    wrapper:
        "v7.2.0/bio/vep/annotate"


rule split_vep:
    input:
        vcf=rules.vep_annotate.output[0],
    output:
        "results/annotated/top_hits.anno.parsed.tsv"
    params:
        attributes=config["params"]["bcftools"]["splitvep"]["attributes"],
        extra=config["params"]["bcftools"]["splitvep"]["extra"]
    log:
        "logs/split_vep.log"
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/split_vep.py"