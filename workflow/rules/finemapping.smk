if config["params"]["susie"]["enabled"]:

    rule run_susie:
        input:
            gwas="results/statistics/lmm.assoc.txt.gz",
            bed="results/processed/genotypes.bed",
            bim="results/processed/genotypes.bim",
            fam="results/processed/genotypes.fam",
        output:
            touch("results/finemapping/.done"),
        params:
            n=449,
            window=config["params"]["susie"]["window_kb"],
            L=config["params"]["susie"]["L"],
            n_loci=config["params"]["susie"]["n_loci"],
        log:
            "logs/run_susie.log"
        conda:
            "../envs/susie.yaml"
        shell:
            "python workflow/scripts/run_susie.py "
            "{input.gwas} results/processed/genotypes "
            "results/finemapping $(which plink) $(which Rscript) "
            "{params.n} {params.window} {params.L} {params.n_loci} "
            ">{log} 2>&1"
