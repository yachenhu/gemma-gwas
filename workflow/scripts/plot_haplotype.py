#!/usr/bin/env python3
"""Haplotype and LD analysis for causal SNPs."""
import sys, os, subprocess as sp
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

out_dir = sys.argv[1]
bfile   = sys.argv[2]  # PLINK bfile for full region
snps_file = sys.argv[3]
plink   = sys.argv[4]

os.makedirs(out_dir, exist_ok=True)

if os.path.getsize(snps_file) == 0:
    print("No causal SNPs found — skipping haplotype analysis.")
    for fname in ["causal_ld_heatmap.png", "genotype_effect.png"]:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.text(0.5, 0.5, "No causal SNPs (PIP > 0.01) identified",
                ha='center', va='center', fontsize=14, transform=ax.transAxes)
        fig.savefig(os.path.join(out_dir, fname), dpi=150, facecolor='white', bbox_inches='tight')
        plt.close()
    sys.exit(0)

# ---- 1. Load genotype data for causal SNPs ----
# Extract genotypes in recoded format (0/1/2)
sp.run([plink, "--bfile", bfile, "--chr", "1",
        "--extract", snps_file, "--recode", "A",
        "--out", f"{out_dir}/causal_geno"],
       stdout=sp.DEVNULL, stderr=sp.DEVNULL)

geno = pd.read_csv(f"{out_dir}/causal_geno.raw", sep=r"\s+")
# PLINK --recode A: columns are FID IID PAT MAT SEX PHENOTYPE rsXXX_A1 ...
snp_cols = [c for c in geno.columns if c.split('_')[0] in
            pd.read_csv(snps_file, header=None)[0].tolist()]
snp_names_clean = [c.split('_')[0] for c in snp_cols]
# Remove duplicate columns (PLINK can have MAT/PAT columns too)
geno_clean = geno[snp_cols].copy()
geno_clean.columns = snp_names_clean

# ---- 2. Load phenotype ----
phe = pd.read_csv("results/statistics/data/phenotype.txt", sep="\t", header=None,
                  names=['pheno'], dtype={0: float})
pheno_values = phe['pheno'].values.astype(float)
print(f"  Phenotype: {len(pheno_values)} samples, "
      f"pheno=0: {(pheno_values==0).sum()}, pheno=1: {(pheno_values==1).sum()}")

# ---- 3. Genotype-phenotype analysis ----
print("=== Genotype-phenotype per SNP ===")
results = []
snp_names = list(geno_clean.columns)
for snp in snp_names:
    g = geno_clean[snp].dropna().astype(int)
    phe_g = pheno_values[:len(g)]

    g0 = phe_g[g == 0]
    g1 = phe_g[g == 1]
    g2 = phe_g[g == 2]

    mean0 = np.mean(g0) if len(g0) > 0 else np.nan
    mean1 = np.mean(g1) if len(g1) > 0 else np.nan
    mean2 = np.mean(g2) if len(g2) > 0 else np.nan

    # Simple chi2 test
    from scipy.stats import chi2_contingency
    table = np.zeros((3, 2))
    for ai, a in enumerate([0, 1, 2]):
        for bi, b in enumerate([0, 1]):
            table[ai, bi] = np.sum((g == a) & (phe_g == b))
    try:
        chi2, p, _, _ = chi2_contingency(table)
    except:
        chi2, p = np.nan, np.nan

    results.append({'SNP': snp, 'mean_0': mean0, 'mean_1': mean1, 'mean_2': mean2, 'p': p})
    print(f"  {snp:20s}  g0={mean0:.3f}(n={len(g0)})  g1={mean1:.3f}(n={len(g1)})  g2={mean2:.3f}(n={len(g2)})  p={p:.2e}")

# ---- 4. LD Heatmap ----
print("\n=== LD Heatmap ===")
geno_filled = geno_clean.fillna(geno_clean.median()).astype(int)
ld_matrix = np.zeros((len(snp_names), len(snp_names)))
for i, s1 in enumerate(snp_names):
    for j, s2 in enumerate(snp_names):
        ld_matrix[i, j] = np.corrcoef(geno_filled[s1], geno_filled[s2])[0, 1] ** 2

fig, ax = plt.subplots(figsize=(12, 10))
im = ax.imshow(ld_matrix, cmap='RdBu_r', vmin=0, vmax=1)
ax.set_xticks(range(len(snp_names)))
ax.set_yticks(range(len(snp_names)))
ax.set_xticklabels([s[:15] for s in snp_names], rotation=45, ha='right', fontsize=8)
ax.set_yticklabels([s[:15] for s in snp_names], fontsize=8)
for i in range(len(snp_names)):
    for j in range(len(snp_names)):
        text = ax.text(j, i, f'{ld_matrix[i,j]:.2f}', ha='center', va='center',
                      color='white' if ld_matrix[i,j] > 0.5 else 'black', fontsize=7)
ax.set_title(f'Pairwise LD (r²) among {len(snp_names)} causal SNPs\n'
             f'(4.58-5.69 Mb, chr1)', fontsize=12)
plt.colorbar(im, ax=ax, shrink=0.8)
plt.tight_layout()
ld_file = os.path.join(out_dir, "causal_ld_heatmap.png")
fig.savefig(ld_file, dpi=200, facecolor='white', bbox_inches='tight')
plt.close()
print(f"  Saved {ld_file}")

# ---- 5. Genotype effect bar plot ----
print("\n=== Genotype Effect Plot ===")
df_r = pd.DataFrame(results)
df_r = df_r.sort_values('p')

fig, axes = plt.subplots(2, 1, figsize=(14, 10), gridspec_kw={'height_ratios': [2, 1]})

# Bar plot: mean phenotype per genotype
x = np.arange(len(df_r))
width = 0.25
ax = axes[0]
for gi, (label, offset) in enumerate([('g0 (hom ref)', -width), ('g1 (het)', 0), ('g2 (hom alt)', width)]):
    vals = [df_r[f'mean_{gi}'].iloc[i] for i in range(len(df_r))]
    ax.bar(x + offset, vals, width, label=label, alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels([s[:15] for s in df_r['SNP']], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('Mean phenotype (feather_color)')
ax.set_title('Genotype-Phenotype Association for Causal SNPs')
ax.legend(fontsize=8)

# -log10 p-value
ax = axes[1]
colors = ['#E41A1C' if p < 0.01 else '#377EB8' for p in df_r['p']]
ax.bar(x, -np.log10(df_r['p']), color=colors, alpha=0.8)
ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.5)
ax.axhline(-np.log10(0.01), color='red', linestyle='--', alpha=0.5)
ax.text(len(df_r)-1, -np.log10(0.05), 'p=0.05', ha='right', va='bottom', fontsize=7, color='grey')
ax.text(len(df_r)-1, -np.log10(0.01), 'p=0.01', ha='right', va='bottom', fontsize=7, color='red')
ax.set_xticks(x)
ax.set_xticklabels([s[:15] for s in df_r['SNP']], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('-log10(p)')
ax.set_xlabel('Causal SNP')

plt.tight_layout()
eff_file = os.path.join(out_dir, "genotype_effect.png")
fig.savefig(eff_file, dpi=200, facecolor='white', bbox_inches='tight')
plt.close()
print(f"  Saved {eff_file}")

# ---- 6. Allele effect direction plot ----
print("\n=== Allele Effect Summary ===")
gwas = pd.read_table("results/statistics/lmm.assoc.txt.gz", sep="\t")
gwas = gwas[gwas['rs'].isin(snp_names)]
for _, r in gwas.iterrows():
    direction = 'RISK (+)' if r['beta'] > 0 else 'PROTECTIVE (-)'
    print(f"  {r['rs']:20s}  {r['allele1']}>{r['allele0']}  "
          f"β={r['beta']:+.4f}  AF={r['af']:.3f}  p={r['p_wald']:.1e}  {direction}")

# ---- 7. Combined haplotype analysis (top multi-SNP combinations) ----
print("\n=== Multi-SNP Genotype Combinations (proxy haplotypes) ===")
# Dynamically group causal SNPs into loci by position proximity (250 kb window)
gwas_causal = gwas[gwas['rs'].isin(snp_names)].sort_values('ps')
loci = {}
window_bp = 250_000
for _, row in gwas_causal.iterrows():
    pos = row['ps']
    placed = False
    for name, (locus_chr, locus_start, locus_end, snps) in loci.items():
        if abs(pos - locus_start) < window_bp or abs(pos - locus_end) < window_bp:
            loci[name] = (locus_chr, min(locus_start, pos), max(locus_end, pos), snps + [row['rs']])
            placed = True
            break
    if not placed:
        loc_name = f"Locus_{len(loci)+1}_{pos/1e6:.2f}Mb"
        loci[loc_name] = (row['chr'], pos, pos, [row['rs']])

for loc_name, (chrom, start, end, loc_snps) in loci.items():
    if len(loc_snps) < 2:
        continue
    valid = [s for s in loc_snps if s in geno_clean.columns]
    if len(valid) < 2:
        continue
    print(f"\n  --- {loc_name} ({len(valid)} SNPs, chr{chrom}:{start}-{end}) ---")
    loc_geno = geno_clean[valid].fillna(1).astype(int)
    combo = loc_geno.apply(lambda row: ''.join(row.astype(str)), axis=1)
    counts = combo.value_counts()

    for hap, count in counts.head(8).items():
        phe_mean = pheno_values[combo == hap].mean() if count > 0 else np.nan
        freq = count / len(combo)
        print(f"    {hap}  n={count:4d}  freq={freq:.3f}  phe_mean={phe_mean:.3f}")

print("\nDone.")
