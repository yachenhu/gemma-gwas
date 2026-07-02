from snakemake.shell import shell

bed_prefix = snakemake.output[0].removesuffix(".bed")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "plink --threads {snakemake.threads} "
    "--vcf {snakemake.input[0]} --const-fid 0 "
    "--make-bed --out {bed_prefix} {extra} {log}"
)
