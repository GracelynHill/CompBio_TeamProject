#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("GEOquery")
#BiocManager::install("scran")

setwd("C:/Users/Gracelyn/Desktop/codezone/Team_Project")

#datasets from GEO
library(GEOquery)
gse <- getGEO("GSE206397")
gpl <- getGEO("GPL20301")

#reading in expression by sample & metadata
samples <- read.table(file.path("GSE206397_countsAll.txt"), header=TRUE)
sra <- read.table(file.path("sra_full.txt"), header=TRUE, sep="\t")
sra2 <- read.table(file.path("sra.txt"), header=TRUE, sep="\t")


#correct list of sample names
runslist <- gse[["GSE206397_series_matrix.txt.gz"]]@phenoData@data[["title"]]

#replacing names with the correct ones
library(stringr)
str_replace_all(runslist, c("-" = "."))
samples <- subset(samples, select = c("SR3.1", "SR3.2", "SR3.3", "SR3.4", "SR5.1", "SR5.2", "SR5.3", "SR5.4", "SR6.1", "SR6.2", "SR6.3", "SR6.4", "SR7.1", "SR7.2", "SR7.3", "SR1.4", "SR1.5",
                                      "SR1.6", "SR1.7", "SR2.2", "SR2.3", "SR2.4", "SR2.5", "SR4.4", "SR4.5", "SR4.6", "SR4.7"))

#swapping out NA

sra[is.na(sra)] <- "unknown"

#DESEQ run 1: by treatment
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = samples,
                              colData = sra, design = ~ TREATMENT)
dds

dds <- DESeq(dds)
res <- results(dds)
res

plotMA(res, ylim=c(-2,2))

resultsNames(dds)

plotCounts(dds, gene=which.min(res$padj), intgroup="TREATMENT")

#distance matrix calculation

vsd <- vst(dds, blind=FALSE)
corrs <- cor(assay(vsd), method="spearman")
corr.dists <- as.dist(1 - corrs)

library(RColorBrewer)
library(pheatmap)
diag(corrs) <- NA
colors <- colorRampPalette(brewer.pal(9, "Blues"))(99)
pheatmap(corrs, 
         clustering_distance_rows=corr.dists,
         clustering_distance_cols=corr.dists,
         col=colors)

#writing distance matrix. May need to change NA to 0 in order to run it with omeClust
write.table(corrs, file = "dist_treatment.txt", sep = "\t")

#At this point you'll have everything you need to run omeClust
#go to terminal and enter the following (adjust as necessary):
# omeClust -i dist_treatment.txt --metadata sra_treat.txt -o OmeClust_treat --plot


#re-run of DESEQ analysis, this time by Genotype

dds <- DESeqDataSetFromMatrix(countData = samples,
                              colData = sra, design = ~ Instrument)
dds

dds <- DESeq(dds)
res <- results(dds)
res

plotMA(res, ylim=c(-2,2))

resultsNames(dds)

plotCounts(dds, gene=which.min(res$padj), intgroup="Instrument")


vsd <- vst(dds, blind=FALSE)
corrs <- cor(assay(vsd), method="spearman")
corr.dists <- as.dist(1 - corrs)

library(RColorBrewer)
library(pheatmap)
diag(corrs) <- NA
colors <- colorRampPalette(brewer.pal(9, "Blues"))(99)
pheatmap(corrs, 
         clustering_distance_rows=corr.dists,
         clustering_distance_cols=corr.dists,
         col=colors)

write.table(corrs, file = "dist_geno.txt", sep = "\t")



library(Tweedieverse)
samples_flip <- data.frame(t(samples[-1]))
colnames(samples_flip) <- rownames(samples)

round(sum(samples_flip==0)/(nrow(samples_flip)*ncol(samples_flip))*100, 1)


library(scran)
library(ggplot2)
library(omePath)
# SCRAN normalization
sce<-t(samples_flip)
scale_factor <- scran::calculateSumFactors(sce)

#Tweedieverse expects input rows and columns to be flipped compared to DESEQ, so I flip them

samples_flip$scale_factor<-scale_factor
samples_flip[is.na(samples_flip)] <- 0

#Tweedieverse by genotype

Tweedieverse::Tweedieverse(
  samples_flip,
  sra2,
  output = 'Genotype_Output', 
  base_model = 'CPLM',
  fixed_effects = c('Instrument'),
  reference = c("Instrument,unknown"),
  link = 'identity',
  standardize = TRUE,
  adjust_offset = FALSE,
  plot_scatter = TRUE)

#Tweediverse by Treatment

Tweedieverse::Tweedieverse(
  samples_flip,
  sra2, 
  output = 'Treatment_Output', # Assuming demo_output exists
  base_model = 'CPLM',
  fixed_effects = c('Treatment'),
  reference = c("Treatment,Pre-reinfection"),
  standardize = TRUE,
  adjust_offset = FALSE,
  plot_scatter = TRUE)

#plotting. There were no Tweedieverse significant associations, so I ran them over all results

library(readr)
geno_assoc <- read_tsv("./Genotype_Output/all_results.tsv")
treatment_assoc <- read_tsv("./Treatment_Output/all_results.tsv")

library(ggplot2)

k <- ggplot(geno_assoc, aes(x=value, y=coef))
k <- k + geom_boxplot()
plot(k)

q <- ggplot(treatment_assoc, aes(x=value, y=coef))
q <- q + geom_boxplot()
plot(q)


