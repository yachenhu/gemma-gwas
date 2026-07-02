import pandas as pd

def encode(data):
    if data.dtype.kind in "biufc":
        return pd.to_numeric(data, errors="raise", downcast="integer")
    data = pd.get_dummies(data, dtype=int, drop_first=True, dummy_na=True)
    mask = data.iloc[:, -1] == 1
    data.loc[mask, :] = float("nan")
    data = data.iloc[:, :-1]
    return data


with open(snakemake.input.id, "r") as file:
    index = [line.strip() for line in file]

df = pd.DataFrame(index=index)
df["intercept"] = 1

covariates = snakemake.params.get("covar", [])
for covariate in covariates:
    data = pd.read_csv(snakemake.input["metadata"], sep="\t",
                       index_col=0, dtype={0: str})[covariate]
    df = df.join(encode(data))

if "pca" in snakemake.input.keys():
    pca = pd.read_csv(snakemake.input.pca, sep=" ", header=None,
                      index_col=[0, 1], dtype={0: str, 1: str})
    pca.index = pca.index.get_level_values(1)
    df = df.join(pca)

df.to_csv(snakemake.output[0], sep="\t", index=False, header=False, na_rep="NA")
