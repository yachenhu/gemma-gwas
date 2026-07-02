from snakemake.shell import shell

attributes = snakemake.params.get("attributes", [])
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

header = "Variant\t" + "\t".join(attributes)
format_expr = r"%ID\t" + r"\t".join([f"%{x}" for x in attributes])

shell(
    "(bcftools +split-vep {snakemake.input.vcf} "
    "--format {format_expr:q} "
    "{extra} "
    "| sed 1i{header:q} "
    ">{snakemake.output[0]}) {log}"
)
