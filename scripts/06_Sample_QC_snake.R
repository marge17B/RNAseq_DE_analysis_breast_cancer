sink(snakemake@log[[1]], split=TRUE)

# Load libraries and data
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(pheatmap)

# Load tximport object and sample metadata
samples <- readRDS(snakemake@input[["samples_rds"]])
rld <- readRDS(snakemake@input[["rld"]])

#plotPCA(rld, intgroup="condition")

pca <- plotPCA(rld, intgroup="condition", returnData=TRUE)
pca$sample_name <- samples$sample_name

#PCA plot with sample labels
pca_plot <- ggplot(pca, aes(PC1, PC2, color = condition, label = sample_name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 4) 

#save as png
ggsave(snakemake@output[["PCA_plot"]], pca_plot, width = 6, height = 5, dpi = 150)

#Inspect Cook's distance for SRR7819995

dds_filter <- readRDS(snakemake@input[["dds_filter"]])
cooks <- assays(dds_filter)[["cooks"]]
summary(cooks[, "SRR7819995"])
summary(cooks[, "SRR7819993"])
summary(cooks[, "SRR7819994"])
summary(cooks[, "SRR7819992"])

#heatmap plot

rld_mat <- assay(rld)  
rld_cor <- cor(rld_mat)
colnames(rld_cor) <- samples$sample_name
rownames(rld_cor) <- samples$sample_name

#annotation for sample conditions
annotation <- data.frame(condition = samples$condition)
rownames(annotation) <- samples$sample_name

#save as png
png(snakemake@output[["QC_heatmap"]], width = 800, height = 800, res = 150)
pheatmap(rld_cor, annotation_col = annotation)
dev.off()