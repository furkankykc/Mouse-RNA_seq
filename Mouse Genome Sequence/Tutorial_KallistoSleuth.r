#load sleuth library
suppressMessages({
      library("sleuth")
      library("dplyr")
})

#set input and output dirs
datapath = "/home/furkankykc/workspace/rnaseq/de/mice/sleuth/input"
resultdir = '/home/furkankykc/workspace/rnaseq/de/mice/sleuth/results'
setwd(resultdir)

#create a sample to condition metadata description
sample_id = dir(file.path(datapath))
kal_dirs = file.path(datapath, sample_id)
print(kal_dirs)
sample = c("NCT_Rep1", "NCT_Rep2", "SLN_Rep1", "SLN_Rep2")
condition = c("NCT", "NCT", "SLN", "SLN")
s2c = data.frame(sample,condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

#run sleuth on the data
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

#summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)
gene_for_plot<- arrange(sleuth_significant,pval)[1,1]
#plot an example DE transcript result
p1 = plot_bootstrap(so, gene_for_plot, units = "est_counts", color_by = "condition")
p2 = plot_pca(so, color_by = 'condition')

#Print out the plots created above and store in a single PDF file
pdf(file="SleuthResults.pdf")
print(p1)
print(p2)
dev.off()

#Quit R
quit(save="no")

