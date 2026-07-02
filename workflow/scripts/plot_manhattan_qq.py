#!/usr/bin/env python3
"""Manhattan + QQ plot using GWASLab."""
import sys
import gwaslab as gl

gwas_file   = sys.argv[1]
out_m       = sys.argv[2]
out_q       = sys.argv[3]
top_n       = int(sys.argv[4]) if len(sys.argv) > 4 else 100

ss = gl.Sumstats(gwas_file,
                 snpid="rs", chrom="chr", pos="ps", p="p_wald",
                 ea="allele1", nea="allele0", eaf="af",
                 build="bGalGal1.mat.broiler.GRCg7b", fmt=None)

ss.plot_mqq(mode="m",
            save=out_m,
            save_args={"dpi": 200, "facecolor": "white"},
            build="bGalGal1.mat.broiler.GRCg7b")

ss.plot_mqq(mode="qq",
            save=out_q,
            save_args={"dpi": 200, "facecolor": "white"},
            build="bGalGal1.mat.broiler.GRCg7b")
print("Done")
