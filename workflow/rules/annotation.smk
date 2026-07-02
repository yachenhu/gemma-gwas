rule extract_causal_snps:
    input:
        vcf="results/processed/genotypes.annotated.vcf.gz",
        ids="results/finemapping/causal_snps.txt",
    output:
        temp("results/annotated/causal.vcf.gz")
    params:
        extra=lambda wildcards, input: f"--include 'ID=@{input.ids}'"
    threads:
        16
    log:
        "logs/extract_causal_snps.log"
    wrapper:
        "v7.2.0/bio/bcftools/view"


rule vep_annotate:
    input:
        calls="results/annotated/causal.vcf.gz",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins",
        fasta="resources/genome.fasta",
        fai="resources/genome.fasta.fai"
    output:
        calls="results/annotated/causal.vep.vcf.gz",
        stats="results/annotated/causal.vep.html",
    params:
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
        "results/annotated/causal.vep.tsv"
    params:
        attributes=config["params"]["bcftools"]["splitvep"]["attributes"],
        extra=config["params"]["bcftools"]["splitvep"]["extra"]
    log:
        "logs/split_vep.log"
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/split_vep.py"
