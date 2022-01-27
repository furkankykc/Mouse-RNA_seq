#!/usr/bin/env Rscript
# Load libraries needed for this analysis
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(ggplot2)
library(gplots)

# Define a path for the output PDF to be written
outfile="Tutorial_Part2_ballgown_output.pdf"

# Load phenotype data
pheno_data = read.csv("SLN_vs_NCT.csv")

# Display the phenotype data
pheno_data

# Load the ballgown object from file
load('bg.rda')
# Load founded significant genes from Part 1
load('sig_genes.rda')
# The load command, loads an R object from a file into memory in our R session. 
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg

# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")
bg_table = texpr(bg, 'all')
row.names(fpkm) <-bg_table$gene_id
bg_gene_names = unique(bg_table[, 9:10])
gene_expression = as.data.frame(gexpr(bg))
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
sig_genes <-arrange(results_genes[sig,],qval)
save(sig_genes, file='sig_genes.rda')
sig=which(results_genes$pval<0.05 & abs(results_genes$fc)>2)
top20_genes <- head(arrange(results_genes[sig,],qval),20)
most_sig_gene_index <- which(geneIDs(bg)==top20_genes$id[1])[1]

# top_genes<- read.csv("SLN_vs_NCT_gene_results_sig.tsv",sep="\t")
# treshold = 3
# top_genes$fc <- log2(top_genes$fc+1)
# sig=which(top_genes$pval<0.05 & abs(top_genes$fc)>treshold)
# top_genes_filtered <- top_genes[sig,]
# top_genes_filtered <- arrange(top_genes_filtered[sig,],qval)
# most_sig_gene_index <- which(geneIDs(bg)==top_genes_filtered$id[1])[1]
# dim(top_genes_filtered)
# head(top_genes_filtered)
# dim(top_genes)
# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(factor(pheno_data$type))+1,las=2,ylab='log2(FPKM+1)')


# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[most_sig_gene_index]

# Display the gene name for a single row of data 
ballgown::geneNames(bg)[most_sig_gene_index]

# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
#plot(fpkm[most_sig_gene_index,] ~ as.numeric(factor(pheno_data$type)), border=c(2,3), main=paste(ballgown::geneNames(bg)[most_sig_gene_index],' : ', ballgown::transcriptNames(bg)[most_sig_gene_index]),pch=19, xlab="Type", ylab='log2(FPKM+1)')
plot(fpkm[most_sig_gene_index,] ~ as.factor(pheno_data$type), border=c(2,3), main=paste(ballgown::geneNames(bg)[most_sig_gene_index],' : ', ballgown::transcriptNames(bg)[most_sig_gene_index]),pch=19, xlab="Type", ylab='log2(FPKM+1)')

# Add the FPKM values for each sample onto the plot 
points(fpkm[most_sig_gene_index,] ~ jitter(as.numeric(factor(pheno_data$type))), col=as.numeric(factor(pheno_data$type))+1, pch=16)

# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[most_sig_gene_index], bg, main=c('Gene in sample NCT_Rep1'), sample=c('NCT_Rep1'))
plotTranscripts(ballgown::geneIDs(bg)[most_sig_gene_index], bg, main=c('Gene in sample NCT_Rep2'), sample=c('NCT_Rep2'))
plotTranscripts(ballgown::geneIDs(bg)[most_sig_gene_index], bg, main=c('Gene in sample SLN_Rep1'), sample=c('SLN_Rep1'))
plotTranscripts(ballgown::geneIDs(bg)[most_sig_gene_index], bg, main=c('Gene in sample SLN_Rep2'), sample=c('SLN_Rep2'))
transcript_gene_table = indexes(bg)$t2g

counts=table(transcript_gene_table[,"g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

#### Plot #2 - the distribution of transcript sizes as a histogram
#In this analysis we supplied StringTie with transcript models so the lengths will be those of known transcripts
#However, if we had used a de novo transcript discovery mode, this step would give us some idea of how well transcripts were being assembled
#If we had a low coverage library, or other problems, we might get short 'transcripts' that are actually only pieces of real transcripts

full_table <- texpr(bg , 'all')
hist(full_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

#### Summarize FPKM values for all 6 replicates
#What are the minimum and maximum FPKM values for a particular library?
min(gene_expression[,"FPKM.SLN_Rep1"])
max(gene_expression[,"FPKM.SLN_Rep1"])

#Set the minimum non-zero FPKM values for use later.
#Do this by grabbing a copy of all data values, coverting 0's to NA, and calculating the minimum or all non NA values
#zz = fpkm_matrix[,data_columns]
#zz[zz==0] = NA
#min_nonzero = min(zz, na.rm=TRUE)
#min_nonzero

#Alternatively just set min value to 1
min_nonzero=1

# Set the columns for finding FPKM and create shorter names for figures
data_columns=c(1:4)
short_names=c("SLN_1","SLN_2","NCT_1","NCT_2")

#### Plot #3 - View the range of values and general distribution of FPKM values for all 4 libraries
#Create boxplots for this purpose
#Display on a log2 scale and add the minimum non-zero value to avoid log2(0)
data_colors=c("tomato1","tomato2","royalblue1","royalblue2")

boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 6 libraries")
#Note that the bold horizontal line on each boxplot is the median

#### Plot #4 - plot a pair of replicates to assess reproducibility of technical replicates
#Tranform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0)
x = gene_expression[,"FPKM.SLN_Rep1"]
y = gene_expression[,"FPKM.SLN_Rep2"]
plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="FPKM (SLN, Replicate 1)", ylab="FPKM (SLN, Replicate 2)", main="Comparison of expression values for a pair of replicates")

#Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)

#Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")

#### Plot #5 - Scatter plots with a large number of data points can be misleading ... regenerate this figure as a density scatter plot
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="FPKM (SLN, Replicate 1)", ylab="FPKM (SLN, Replicate 2)", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)

#### Plot all sets of replicates on a single plot
#Create an function that generates an R plot.  This function will take as input the two libraries to be compared and a plot name and color
plotCor = function(lib1, lib2, name, color){
	x=gene_expression[,lib1]
	y=gene_expression[,lib2]
	zero_count = length(which(x==0)) + length(which(y==0))
	plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col=color, cex=0.25, xlab=lib1, ylab=lib2, main=name)
	abline(a=0,b=1)
	rs=cor(x,y, method="pearson")^2
	legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
	legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}
#Open a plotting page with room for two plots on one page
par(mfrow=c(1,2))

#Plot #6 - Now make a call to our custom function created above, once for each library comparison
plotCor("FPKM.SLN_Rep1", "FPKM.NCT_Rep1", "SLN_1 vs NCT_1", "tomato2")
plotCor("FPKM.SLN_Rep2", "FPKM.NCT_Rep2", "SLN_2 vs NCT_2", "royalblue2")


##### One problem with these plots is that there are so many data points on top of each other, that information is being lost
#Regenerate these plots using a density scatter plot
plotCor2 = function(lib1, lib2, name, color){
	x=gene_expression[,lib1]
	y=gene_expression[,lib2]
	zero_count = length(which(x==0)) + length(which(y==0))
	colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab=lib1, ylab=lib2, main=name, colramp=colors, nbin=275)
	abline(a=0,b=1)
	rs=cor(x,y, method="pearson")^2
	legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
	legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}

#### Plot #7 - Now make a call to our custom function created above, once for each library comparison
par(mfrow=c(1,2))
plotCor2("FPKM.SLN_Rep1", "FPKM.NCT_Rep1", "SLN_1 vs NCT_1", "tomato2")
plotCor2("FPKM.SLN_Rep2", "FPKM.NCT_Rep2", "SLN_2 vs NCT_2", "royalblue2")


#### Compare the correlation 'distance' between all replicates
#Do we see the expected pattern for all eight libraries (i.e. replicates most similar, then tumor vs. normal)?

#Calculate the FPKM sum for all 6 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

#Identify the genes with a grand sum FPKM of at least 5 - we will filter out the genes with very low expression across the board
i = which(gene_expression[,"sum"] > 5)

#Calculate the correlation between all pairs of data
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")

#Print out these correlation values
r

#### Plot #8 - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries
#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together
#What pattern do we expect to see, given the types of libraries we have (technical replicates, biologal replicates, tumor/normal)?
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)")
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

# Calculate the differential expression results including significance
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes,bg_gene_names,by.x=c("id"),by.y=c("gene_id"))

#### Plot #9 - View the distribution of differential expression values as a histogram
#Display only those that are significant according to Ballgown

sig=which(results_genes$pval<0.05)
results_genes[,"de"] = log2(results_genes[,"fc"])
hist(results_genes[sig,"de"], breaks=50, col="seagreen", xlab="log2(Fold change) SLN vs NCT", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)

#### Plot #10 - Display the grand expression values from SLN and NCT and mark those that are significantly differentially expressed
gene_expression[,"SLN"]=apply(gene_expression[,c(1:2)], 1, mean)
gene_expression[,"NCT"]=apply(gene_expression[,c(3:5)], 1, mean)

x=log2(gene_expression[,"SLN"]+min_nonzero)
y=log2(gene_expression[,"NCT"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="SLN FPKM (log2)", ylab="NCT FPKM (log2)", main="SLN vs NCT FPKMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
topn = order(abs(results_genes[sig,"fc"]), decreasing=TRUE)[1:25]
topn = order(results_genes[sig,"qval"])[1:25]
text(x[topn], y[topn], results_genes[topn,"gene_name"], col="black", cex=0.75, srt=45)


#### Write a simple table of differentially expressed transcripts to an output file
#Each should be significant with a log2 fold-change >= 2
sigpi = which(results_genes[,"pval"]<0.05)
sigp = results_genes[sigpi,]
sigde = which(abs(sigp[,"de"]) >= 2)
sig_tn_de = sigp[sigde,]

#Order the output by or p-value and then break ties using fold-change
o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"de"]), decreasing=FALSE)

output = sig_tn_de[o,c("gene_name","id","fc","pval","qval","de")]
write.table(output, file="SigDE_supplementary_R.txt", sep="\t", row.names=FALSE, quote=FALSE)

#View selected columns of the first 25 lines of output
output[1:25,c(1,4,5)]

#You can open the file "SigDE.txt" in Excel, Calc, etc.
#It should have been written to the current working directory that you set at the beginning of the R tutorial
dir()


#### Plot #11 - Create a heatmap to vizualize expression differences between the eight samples
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

main_title="sig DE Transcripts"
par(cex.main=0.8)
sig_genes_de=sig_tn_de[,"id"]
sig_gene_names_de=sig_tn_de[,"gene_name"]

data=log2(as.matrix(gene_expression[as.vector(sig_genes_de),data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sig_gene_names_de,col=rev(heat.colors(75)))

#### Plot #12 - Volcano plot

# default all genes to "no change"
results_genes$diffexpressed <- "No"
# if log2Foldchange > 2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$de > 0.6 & results_genes$pval < 0.05] <- "Up"
# if log2Foldchange < -2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$de < -0.6 & results_genes$pval < 0.05] <- "Down"

results_genes$gene_label <- NA
# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label[results_genes$diffexpressed != "No"] <- results_genes$gene_name[results_genes$diffexpressed != "No"]

ggplot(data=results_genes[results_genes$diffexpressed != "No",], aes(x=de, y=-log10(pval), label=gene_label, color = diffexpressed)) +
             xlab("log2Foldchange") +
             scale_color_manual(name = "Differentially expressed", values=c("blue", "red")) +
             geom_point() +
             theme_minimal() +
             geom_text_repel() +
             geom_vline(xintercept=c(-0.6, 0.6), col="red") +
             geom_hline(yintercept=-log10(0.05), col="red") +
             guides(colour = guide_legend(override.aes = list(size=5))) +
             geom_point(data = results_genes[results_genes$diffexpressed == "No",], aes(x=de, y=-log10(pval)), colour = "black")




#plotMeans('TST',bg,groupvar="type",legend=FALSE)

# Close the PDF device where we have been saving our plots
dev.off()

# Exit the R session
quit(save="no")
