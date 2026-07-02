import pandas as pd

with open(snakemake.input.id, "r") as file:
    index = [line.strip() for line in file]

df = pd.read_csv(snakemake.input[0], sep="\t",
                 index_col=0, dtype={0: str})
df = df[df.index.notnull()]
df = df.reindex(index)

data = pd.to_numeric(df.loc[:, snakemake.params.trait], errors="raise", downcast="integer")
data.to_csv(snakemake.output[0], index=False, header=False, na_rep="NA")
