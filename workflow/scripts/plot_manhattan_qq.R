library(dplyr)
library(manhattanly)
library(htmlwidgets)

table <- read.csv2(snakemake@input[[1]], sep = "\t")

table$BP <- as.numeric(table$BP)
table$P <- as.numeric(table$P)

top <- snakemake@params$top
highlight <- NULL
if (!is.null(top)) {
    sorted <- table[order(table$P), ]
    highlight <- head(sorted$SNP, top)
}

fig <- manhattanly(table, snp = "SNP", gene = "GENE", highlight = highlight)
saveWidget(fig, snakemake@output[[1]])

fig <- qqly(table)
saveWidget(fig, snakemake@output[[2]])
