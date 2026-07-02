import pandas as pd

df = (pd.read_table(snakemake.input[0], dtype=str)
      .set_index(snakemake.params.columns["SNP"], drop=False))
anno = pd.read_table(snakemake.input["anno"], dtype=str, index_col=[0], na_values=".")

df = (
    df.join(anno, how="right")
      .rename(columns={value: key for key, value in snakemake.params.columns.items()})
      .loc[:, snakemake.params.columns.keys()]
      )

df.to_csv(snakemake.output[0], sep="\t", index=False)
