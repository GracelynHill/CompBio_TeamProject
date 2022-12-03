#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("GEOquery")
library(GEOquery)
gse <- getGEO("GSE206397")
gpl <- getGEO("GPL20301")

samples <- read.table(file.path("GSE206397_countsAll.txt"), header=TRUE)
sra <- read.csv(file.path("SraRunTable.csv"), header=TRUE)


runslist <- gse[["GSE206397_series_matrix.txt.gz"]]@phenoData@data[["title"]]

library(stringr)
str_replace_all(runslist, c("-" = "."))

samples <- subset(samples, select = c("SR3.1", "SR3.2", "SR3.3", "SR3.4", "SR5.1", "SR5.2", "SR5.3", "SR5.4", "SR6.1", "SR6.2", "SR6.3", "SR6.4", "SR7.1", "SR7.2", "SR7.3", "SR1.4", "SR1.5",
                                       "SR1.6", "SR1.7", "SR2.2", "SR2.3", "SR2.4", "SR2.5", "SR4.4", "SR4.5", "SR4.6", "SR4.7"))

sra[is.na(sra)] <- "unknown"

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = samples,
                              colData = sra, design = ~ Genotype)
dds

dds <- DESeq(dds)
res <- results(dds)
res

plotMA(res, ylim=c(-2,2))

resultsNames(dds)

plotCounts(dds, gene=which.min(res$padj), intgroup="Genotype")
