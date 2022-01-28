#!/usr/bin/env Rscript
library(dplyr)
load('../sig_genes.rda')
top20 = arrange(sig_genes[1:20,],pval)
top = arrange(sig_genes,pval)
head(top20)
write.table(top20, file = "top20_genes.tsv", row.names=FALSE, sep="\t")
write.table(top, file = "top_genes.tsv", row.names=FALSE, sep="\t")
