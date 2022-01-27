#!/usr/bin/env Rscript
library(clusterProfiler)
library(GSEABase)
library(org.Mm.eg.db)
top_genes <- read.table(file = 'top20_genes.tsv', sep = '\t', header = TRUE)
head(top_genes)
filename<- "c7.all.v7.1.entrez.gmt"
gmtfile <- system.file(filename)

t <- c(top_genes$gene_name)

et <- bitr(t, fromType="SYMBOL", toType=(c("ENTREZID","GENENAME")), OrgDb="org.Mm.eg.db")
head(et)
#c6 <- read.gmt(gmtfile)
c6 <- read.gmt(filename)
yourEntrezIdList<- c(et$ENTREZID) #USE YOUR OWN ENTREZ GENEID LIST HERE 
ImmunSigEnrich <- enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.01) 
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenes.csv") 
goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Mm.eg.db, ont= "ALL",pAdjustMethod= "BH",pvalueCutoff = 0.01,readable= TRUE) 
write.csv(goEnrich,"MyGORelatedGenes.csv")
keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "mmu",pAdjustMethod="BH", pvalueCutoff = 0.01)
write.csv(keggEnrich,"MyKEGGRelatedGenes.csv")
