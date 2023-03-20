if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#install the library in r/4.1.2
library(DESeq2)


setwd("/global/scratch/hpcXXXX/analysis/RNA_Seq/Counts")

#Load the files
control_files <- list.files(path = ".", pattern = "CR.gene.count", full.names = TRUE)
case_files <- list.files(path = ".", pattern = "Dx.gene.count", full.names = TRUE)
positive_files <- list.files(path = "rela/", pattern = "-Rela.gene.count", full.names = TRUE)

#read in counts data
count_data <- list()

for (file in control_files) {
  sample_name <- gsub(".gene.count", "", basename(file))
  count_data[[sample_name]] <- read.table(file, sep = "\t", header = F, row.names = 1)
}


for (file in case_files) {
  sample_name <- gsub(".gene.count", "", basename(file))
  count_data[[sample_name]] <- read.table(file, sep = "\t", header = F, row.names = 1)
}

for (file in positive_files) {
  sample_name <- gsub(".gene.count", "", basename(file))
  count_data[[sample_name]] <- read.table(file, sep = "\t", header = F, row.names = 1)
}

#create the meta data from the file names
metadata <- data.frame(
  sample = names(count_data),
  type = c(rep("CR", length(control_files)), rep("Dx", length(case_files)), rep("Rela", length(positive_files)))
)



####DeSeq2 dataset object

dds <- DESeqDataSetFromMatrix(countData = do.call(cbind, count_data),
                              colData = metadata,
                              design = ~ type)

#filtering the rows with low gene counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Run the DESeq
dds <- DESeq(dds)
res <- results(dds)
de_genes <- res[which(res$padj < 0.05),]

#plots to visualize the results

jpeg('../volcano_plot.jpeg')
plot(log2FoldChange(dds), -log10(pvalue(dds)), pch=20, cex=0.4, main="Volcano plot")
abline(h = -log10(0.05), col = "red", lty = 2)
abline(v = c(-1,1), col = "blue", lty = 2)
dev.off()
