## Load libraries and data
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(pheatmap)

# Load tximport object and sample metadata
samples <- readRDS("data/rds/samples.rds")
rld <- readRDS("data/rds/rld.rds")

#plotPCA(rld, intgroup="condition")

pca <- plotPCA(rld, intgroup="condition", returnData=TRUE)
pca$sample_name <- samples$sample_name

#Custom PCA plot with sample labels
pca_plot <- ggplot(pca, aes(PC1, PC2, color = condition, label = sample_name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 4) 

#save as png
ggsave("results/figures_plots/PCA_plot.png", pca_plot, width = 6, height = 5, dpi = 150)

# HEATMAP CREATION 

#matrix
rld_mat <- assay(rld)  
rld_cor <- cor(rld_mat)    ## cor() is a base R function
colnames(rld_cor) <- samples$sample_name
rownames(rld_cor) <- samples$sample_name

#annotation for sample conditions
annotation <- data.frame(condition = samples$condition)
rownames(annotation) <- samples$sample_name

#save
png("results/figures_plots/QC_heatmap.png", width = 800, height = 800, res = 150)
pheatmap(rld_cor, annotation_col = annotation)
dev.off()