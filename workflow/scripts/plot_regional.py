#!/usr/bin/env python3
"""Regional association plots using GWASLab."""
import sys
import os

# Patch GWASLab 3.6.9 bug: stray region_ref=None in matplotlib text() call
import matplotlib.axes
_orig_axes_text = matplotlib.axes.Axes.text
def _safe_text(self, x, y, s, *args, **kwargs):
    kwargs.pop("region_ref", None)
    kwargs.pop("region_ref_total_n", None)
    return _orig_axes_text(self, x, y, s, *args, **kwargs)
matplotlib.axes.Axes.text = _safe_text

import gwaslab as gl

gwas_file = sys.argv[1]
vcf_file  = sys.argv[2]
gtf_file  = sys.argv[3]
gtf_gene_key = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] else "gene_name"
gtf_source   = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] else "ensembl"
n_top     = int(sys.argv[6])
window    = int(sys.argv[7])
out_dir   = sys.argv[8]

os.makedirs(out_dir, exist_ok=True)

# Build chr mapping for NCBI RefSeq chromosome names
_kwargs = {}
if gtf_source == "ncbi":
    _gtf_chr_dict = {}
    for _i in range(1, 42):
        _gtf_chr_dict[str(_i)] = f"NC_0{52531 + _i:05d}.1"
    _gtf_chr_dict["W"] = "NC_052573.1"
    _gtf_chr_dict["Z"] = "NC_052574.1"
    _gtf_chr_dict["MT"] = "NC_001323.1"
    for _k, _v in list(_gtf_chr_dict.items()):
        try: _gtf_chr_dict[int(_k)] = _v
        except: pass
    _kwargs["gtf_chr_dict"] = _gtf_chr_dict

print(f"Loading GWAS: {gwas_file}")
sumstats = gl.Sumstats(
    gwas_file,
    snpid="rs",
    chrom="chr",
    pos="ps",
    p="p_wald",
    ea="allele1",
    nea="allele0",
    eaf="af",
    build="bGalGal1.mat.broiler.GRCg7b",
    fmt=None,
)

d = sumstats.data.sort_values("P")
print(f"Loaded {len(d)} variants. Top p-value: {d['P'].iloc[0]:.2e}")

# Find top independent SNPs (>1Mb apart)
# GWASLab uses Int64 for CHR — keep as-is for internal comparison
leads = []
for _, row in d.iterrows():
    if len(leads) >= n_top:
        break
    chr_i = row["CHR"]
    pos_i = int(row["POS"])
    too_close = any(
        l["CHR"] == chr_i and abs(int(l["POS"]) - pos_i) < window
        for l in leads
    )
    if not too_close:
        leads.append(row)
        print(f"  Lead #{len(leads)}: chr{chr_i}:{pos_i}  {row['SNPID']}  p={row['P']:.2e}")

# Plot each region — pass CHR as int (GWASLab uses Int64 for CHR, vcf_chr_dict has int keys)
for i, lead in enumerate(leads):
    chr_int = int(lead["CHR"])
    pos_center = int(lead["POS"])
    start = pos_center - window
    end = pos_center + window

    out_file = os.path.join(out_dir, f"regional_chr{chr_int}_{pos_center}.png")
    print(f"\nPlotting chr{chr_int}:{start}-{end} -> {out_file}")

    try:
        sumstats.plot_mqq(
            mode="r",
            region=(chr_int, start, end),
            vcf_path=vcf_file,
            gtf_path=gtf_file,
            gtf_gene_name=gtf_gene_key,
            region_title=f"chr{chr_int}:{start:,} - {end:,}",
            region_lead_grid=True,
            region_hspace=0.02,
            region_step=21,
            build="bGalGal1.mat.broiler.GRCg7b",
            save=out_file,
            save_args={"dpi": 200, "facecolor": "white"},
            **_kwargs,
        )
        print(f"  -> Saved")
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

print("\nDone.")
