FILTER_PRESETS = {
    "relaxed":  "INFO_SCORE >= 0.6 && MAF >= 0.05 && HWE >= 1e-6 && F_MISSING <= 0.1",
    "standard": "INFO_SCORE >= 0.8 && MAF >= 0.05 && HWE >= 1e-6 && F_MISSING <= 0.1",
    "strict":   "INFO_SCORE >= 0.9 && MAF >= 0.05 && HWE >= 1e-4 && F_MISSING <= 0.1",
}

rule filter_variants:
    input:
        config["vcf"]
    output:
        "results/processed/genotypes.filtered.vcf.gz"
    params:
        enabled=config["params"]["filter"].get("enabled", True),
        preset=config["params"]["filter"].get("preset", "standard"),
        custom=config["params"]["filter"].get("custom", ""),
        subset=config.get("chr_subset", ""),
    threads:
        16
    log:
        "logs/filter_variants.log"
    conda:
        "../envs/bcftools.yaml"
    run:
        from snakemake.shell import shell
        if params.enabled:
            expr = params.custom if params.preset == "custom" else FILTER_PRESETS[params.preset]
            extra = "-v snps -m2 -M2"
            if params.subset:
                shell(
                    "bcftools view -r '{params.subset}' {input} "
                    "| bcftools view {extra} "
                    "| bcftools filter --threads {threads} -i '{expr}' "
                    "| bgzip -@ {threads} -c > {output} 2> {log}"
                )
            else:
                shell(
                    "bcftools view {extra} {input} "
                    "| bcftools filter --threads {threads} -i '{expr}' "
                    "| bgzip -@ {threads} -c > {output} 2> {log}"
                )
        else:
            if params.subset:
                shell(
                    "bcftools view -r '{params.subset}' {input} "
                    "| bgzip -@ {threads} -c > {output} 2> {log}"
                )
            else:
                shell("cp {input} {output}")


rule index_filtered_vcf:
    input:
        "results/processed/genotypes.filtered.vcf.gz"
    output:
        "results/processed/genotypes.filtered.vcf.gz.tbi"
    threads:
        16
    log:
        "logs/index_filtered_vcf.log"
    wrapper:
        "v7.2.0/bio/bcftools/index"


rule update_ids:
    input:
        vcf="results/processed/genotypes.filtered.vcf.gz",
        idx="results/processed/genotypes.filtered.vcf.gz.tbi",
        annotations="resources/variation.noiupac.vcf.gz",
        anno_idx="resources/variation.noiupac.vcf.gz.tbi",
    output:
        "results/processed/genotypes.annotated.vcf.gz",
    threads:
        16
    log:
        "logs/update_ids.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools annotate --threads {threads} "
        "-a {input.annotations} -c ID "
        "--set-id +'%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' "
        "-Oz -o {output} {input.vcf} 2>>{log}"


rule vcf2plink:
    input:
        "results/processed/genotypes.annotated.vcf.gz"
    output:
        multiext("results/processed/genotypes", ".bed", ".bim", ".fam")
    params:
        extra=config["params"]["plink"]
    shadow:
        "minimal"
    log:
        "logs/vcf2plink.log"
    conda:
        "../envs/plink.yaml"
    script:
        "../scripts/plink_data.py"
