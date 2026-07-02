import pandas as pd

def encode_sex(x):
    try:
        first = x.capitalize()[0]
    except Exception:
        return 0
    if first == "M":
        return 1
    elif first == "F":
        return 2
    return 0


with open(snakemake.input[0], "r") as file:
    index = [line.strip() for line in file]

df = pd.read_table(snakemake.input["metadata"], dtype=str)
df = df.set_index(df.columns[0])
df = df.reindex(index)

table = pd.DataFrame(index=index)
table["FID"] = 0
table["IID"] = index
table["PAT"] = df[snakemake.params.get("sire")]
table["MAT"] = df[snakemake.params.get("dam")]
table["SEX"] = df[snakemake.params.get("sex")].apply(encode_sex)
table["PHENO"] = -9

table = table.fillna(0)
table[["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv(
    snakemake.output[0], sep="\t", header=False, index=False
)
