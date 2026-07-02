rule haplotype_analysis:
    input:
        gwas="results/statistics/lmm.assoc.txt.gz",
        bed="results/processed/genotypes.bed",
        bim="results/processed/genotypes.bim",
        fam="results/processed/genotypes.fam",
        pipcs="results/finemapping/all_pipcs.csv",
        causal="results/finemapping/causal_snps.txt",
    output:
        "results/haplotype/causal_ld_heatmap.png",
        "results/haplotype/genotype_effect.png",
    log:
        "logs/haplotype_analysis.log"
    conda:
        "../envs/susie.yaml"
    shell:
        "mkdir -p results/haplotype && "
        "python workflow/scripts/plot_haplotype.py "
        "results/haplotype results/processed/genotypes {input.causal} $(which plink) >{log} 2>&1"
