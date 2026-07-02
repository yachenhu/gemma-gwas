#!/usr/bin/env python3
"""SuSiE RSS fine-mapping for GWAS loci."""
import subprocess as sp
import os, sys
import pandas as pd
import numpy as np
import gc

gwas_file  = sys.argv[1]
bfile      = sys.argv[2]
out_dir    = sys.argv[3]
plink      = sys.argv[4]
rscript    = sys.argv[5]
n_sample   = int(sys.argv[6])
window_kb  = int(sys.argv[7])
L          = int(sys.argv[8])
n_loci     = int(sys.argv[9]) if len(sys.argv) > 9 else 5

os.makedirs(out_dir, exist_ok=True)

# ---- Load GWAS data ----
print(f"Loading {gwas_file}")
gwas = pd.read_table(gwas_file, sep="\t")
gwas = gwas[gwas['p_wald'].notna() & (gwas['p_wald'] > 0)]
gwas['CHR']   = gwas['chr'].astype(str)
gwas['POS']   = gwas['ps'].astype(int)
gwas['SNPID'] = gwas['rs']
gwas['P']     = gwas['p_wald'].astype(float)
gwas['BETA']  = gwas['beta'].astype(float)
gwas['SE']    = gwas['se'].astype(float)
gwas['EA']    = gwas['allele1']
gwas['NEA']   = gwas['allele0']
print(f"  {len(gwas)} variants")

# ---- Find lead SNPs ----
print("Finding lead SNPs...")
sig = gwas[gwas['P'] < 5e-8].sort_values('P')
leads = []
window_bp = window_kb * 1000
for _, row in sig.iterrows():
    if any(abs(row['POS'] - l['POS']) < window_bp and row['CHR'] == l['CHR']
           for l in leads):
        continue
    leads.append(row)

leads = leads[:n_loci]
print(f"  Top {len(leads)} loci selected")

# ---- For each locus, compute LD and run SuSiE ----
all_pips = []
for i, lead in enumerate(leads):
    chrom  = lead['CHR']
    pos    = lead['POS']
    snpid  = lead['SNPID']
    start  = max(0, pos - window_bp)
    end    = pos + window_bp
    name   = f"locus_{i+1}_{snpid}"
    print(f"\n=== {name}: chr{chrom}:{start}-{end} ===")

    locus = gwas[(gwas['CHR'] == chrom) &
                 (gwas['POS'] >= start) &
                 (gwas['POS'] <= end)].copy().sort_values('POS')
    if len(locus) < 10:
        print("  Too few variants, skipping")
        continue
    print(f"  {len(locus)} variants in region")

    # Write sumstats and snplist
    ss_file  = os.path.join(out_dir, f"{name}.sumstats")
    snp_file = os.path.join(out_dir, f"{name}.snplist")
    locus[['SNPID','CHR','POS','EA','NEA','BETA','SE','P']].to_csv(
        ss_file, sep="\t", index=False)
    locus[['SNPID']].to_csv(snp_file, index=False, header=False)

    # PLINK LD matrix
    ld_file = os.path.join(out_dir, f"{name}.ld")
    if not os.path.exists(ld_file):
        sp.run([plink, "--bfile", bfile, "--chr", chrom,
                "--from-bp", str(start), "--to-bp", str(end),
                "--extract", snp_file, "--r2", "square",
                "--allow-extra-chr", "--chr-set", "95",
                "--keep-allele-order", "--out", f"{out_dir}/{name}"],
               stdout=sp.DEVNULL, stderr=sp.DEVNULL)

    if not os.path.exists(ld_file):
        print("  LD computation failed")
        continue

    # Align LD matrix with sumstats
    ld_data = np.loadtxt(ld_file)
    snplist_ids = pd.read_csv(snp_file, header=None)[0].tolist()
    n_ld = ld_data.shape[0]
    # PLINK may return slightly different count (SNP ID duplicates in bim)
    if n_ld != len(snplist_ids):
        print(f"  LD size {n_ld} vs snplist {len(snplist_ids)}, aligning...")
    ld_snps = snplist_ids[:n_ld] if n_ld <= len(snplist_ids) else snplist_ids
    if n_ld > len(snplist_ids):
        ld_snps += [f"UNK{j}" for j in range(n_ld - len(snplist_ids))]

    common = [s for s in locus['SNPID'] if s in ld_snps]
    if len(common) < 10:
        print(f"  Too few overlapping SNPs: {len(common)}")
        continue
    idx = np.array([ld_snps.index(s) for s in common])
    R = ld_data[np.ix_(idx, idx)]
    n_use = len(common)
    print(f"  LD matrix: {n_use}x{n_use}")

    # Write matched sumstats and R matrix
    matched_ss = os.path.join(out_dir, f"{name}_matched.sumstats")
    matched_locus = locus[locus['SNPID'].isin(common)].sort_values('POS')
    matched_locus[['SNPID','CHR','POS','EA','NEA','BETA','SE','P']].to_csv(
        matched_ss, sep="\t", index=False)
    r_mat_file = os.path.join(out_dir, f"{name}_R.csv")
    np.savetxt(r_mat_file, R, delimiter=",")
    del ld_data, R, idx
    gc.collect()

    # Write and run R script
    r_script = f'''
library(susieR)
sumstats <- read.csv("{matched_ss}", sep="\\t")
R <- as.matrix(read.csv("{r_mat_file}", header=FALSE))
z <- sumstats$BETA / sumstats$SE
n <- {n_sample}
fitted <- susie_rss(z = z, R = R, n = n, L = {L}, max_iter = 100,
                    min_abs_corr = 0.5, refine = FALSE)
result <- summary(fitted)$vars
if (nrow(result) > 0) {{
    result$SNPID <- sumstats$SNPID[result$variable]
    result$LOCUS <- "{name}"
    write.csv(result, "{out_dir}/{name}.pipcs", row.names = FALSE)
    cat("  Credible sets:", nrow(result), "\\n")
}} else {{
    cat("  No credible sets found\\n")
}}
'''
    r_file = os.path.join(out_dir, f"{name}.R")
    with open(r_file, 'w') as f:
        f.write(r_script)

    print("  Running SuSiE RSS...")
    result = sp.run([rscript, r_file], capture_output=True, text=True, timeout=600)
    if result.stdout:
        print(result.stdout.strip())
    for line in result.stderr.split('\n'):
        if 'Error' in line or 'Warning' in line:
            print(f"  R: {line.strip()}")

    # Read results
    pip_file = os.path.join(out_dir, f"{name}.pipcs")
    if os.path.exists(pip_file):
        pips = pd.read_csv(pip_file)
        all_pips.append(pips)
        top = pips[pips['variable_prob'] > 0.01].sort_values('variable_prob', ascending=False)
        for _, r in top.iterrows():
            print(f"    {r['SNPID']:20s} PIP={r['variable_prob']:.4f}  CS={r['cs']}")

    # Clean up large intermediate files
    for f in [ld_file, r_mat_file, matched_ss]:
        if os.path.exists(f):
            os.remove(f)

# ---- Save combined results ----
out_csv = os.path.join(out_dir, "all_pipcs.csv")
snp_file = os.path.join(out_dir, "causal_snps.txt")
if all_pips:
    combined = pd.concat(all_pips, ignore_index=True)
    combined.to_csv(out_csv, index=False)
    print(f"\nResults: {len(combined)} variants, saved to {out_csv}")
    top = combined[combined['variable_prob'] > 0.01].sort_values('variable_prob', ascending=False)
    if len(top) > 0:
        print("High-confidence causal variants (PIP>0.01):")
        causal_snps = []
        for _, r in top.iterrows():
            print(f"  {r['SNPID']:20s} {r.get('LOCUS',''):30s} PIP={r['variable_prob']:.4f}  CS={r['cs']}")
            causal_snps.append(r['SNPID'])
        with open(snp_file, 'w') as f:
            f.write('\n'.join(causal_snps) + '\n')
        print(f"  Causal SNP list saved to {snp_file}")
    else:
        with open(snp_file, 'w') as f:
            f.write('')
else:
    print("No results generated.")
    pd.DataFrame(columns=['variable','SNPID','variable_prob','cs','LOCUS']).to_csv(out_csv, index=False)
    with open(snp_file, 'w') as f:
        f.write('')
