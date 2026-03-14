sink(snakemake@log[[1]], split=TRUE)

## Load libraries and data
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(apeglm)
library(ggrepel)
library(RColorBrewer)

# Load previously saved objects
txi <- readRDS(snakemake@input[["txi"]])
samples <- readRDS(snakemake@input[["samples_rds"]])
dds <- readRDS(snakemake@input[["dds_filter"]])
rld <- readRDS(snakemake@input[["rld"]])

# Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds, name="condition_treated_vs_control")

# MA plot

# check LFC range and stats
range(res$log2FoldChange, na.rm = TRUE)
summary(res$log2FoldChange)
#range:-23.72896  25.24171
png(snakemake@output[["MA_plot_all_range"]],
    width = 1200, height = 1200, res = 150)

plotMA(res, ylim = c(-30, 30))
dev.off()

# Shrink Log2 Fold Changes and MA Plot
res_shrink <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")

#Check range and summary
range(res_shrink$log2FoldChange, na.rm = TRUE)
summary(res_shrink$log2FoldChange)

png(snakemake@output[["MA_plot_shrunken"]],
    width = 1200, height = 1200, res = 150)

plotMA(res_shrink, ylim = c(-6, 6))

dev.off()

# Gene Annotation Using GTF

#Read gtf
gtf <- read.table(snakemake@input[["gtf_path"]], header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gtf_gene <- gtf[gtf$V3 == "gene", ]

#Parse attribute column
gtf_parsed <- gtf_gene %>%
  mutate(id = row_number()) %>%                 
  separate_rows(V9, sep="; ") %>%                
  separate(V9, c("name","value"), sep=" ") %>%   
  mutate(value = gsub(";","",value)) %>%         
  spread(name, value)

#View(gtf_parsed)

gtf_gene_names_id <- gtf_parsed %>% select(gene_id, gene_name)

#Convert DESeq2 results to dataframe
res_df <- as.data.frame(res)

res_df <- res_df %>%
  tibble::rownames_to_column(var = "gene_id")

#Merge
res_annotated <- merge(res_df, gtf_gene_names_id, by = "gene_id", all.x = TRUE)

# Save as RDS 
saveRDS(res_annotated,snakemake@output[["res_annotated"]])

# Volcano plot

res_annotated <- res_annotated %>%
  mutate(significant = case_when(padj < 0.05 & log2FoldChange > 0.6  ~ "Upregulated", padj < 0.05 & log2FoldChange < -0.6 ~ "Downregulated", TRUE ~ "Not significant"))

labels_df <- res_annotated %>%
  filter(significant != "Not significant") %>%
  arrange(padj) %>% 
  head(30)  

volcano_plot <- ggplot(res_annotated, aes(x = log2FoldChange, y = -log10(padj), col = significant))  +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +  
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("Downregulated" = "#00AFBB","Not significant" = "grey", "Upregulated" = "#FFDB6D")) +  geom_text_repel(data = labels_df,
                                                                                                                                       aes(label = gene_name), max.overlaps = Inf,
                                                                                                                                       box.padding = 0.4, point.padding = 0.2, segment.size = 0.3) 
ggsave(
  filename = snakemake@output[["volcano_plot"]],
  plot = volcano_plot,
  width = 8,
  height = 7,
  dpi = 300
)

# Differentially Expressed Genes

upregulated <- sum(res_annotated$significant == "Upregulated")
downregulated <- sum(res_annotated$significant == "Downregulated")
total_DE <- upregulated + downregulated

print(paste("Upregulated genes:", upregulated)) 
print(paste("Downregulated genes:", downregulated)) 
print(paste("Total DE genes:", total_DE)) 

# Top Differentially Expressed Genes

res_ordered <- res_annotated[order(res_annotated$padj), ]

#top 50 genes
top50_genes <- res_ordered %>%
  slice(1:50)

write.csv(top50_genes,
          snakemake@output[["top50_DE_genes"]])

# Heatmap for top genes 

top_gene_ids <- res_ordered$gene_id[1:50]
top_gene_names <- res_ordered$gene_name[1:50]
rld_mat <- assay(rld)
mat <- rld_mat[top_gene_ids, ]

#replace id to name
rownames(mat) <- top_gene_names

#Change column names (SRR → control_1 etc.)
samples_names_condition <- samples[colnames(mat), ] 
colnames(mat) <- samples_names_condition$sample_name
annotation_col <- data.frame(condition = samples_names_condition$condition)
rownames(annotation_col) <- samples_names_condition$sample_name

png(snakemake@output[["heatmap_top50"]], width = 1200, height = 1200, res = 150)

pheatmap(mat,
            scale='row',
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            annotation_col = annotation_col,
            show_rownames = TRUE,
            fontsize_row = 6,
            fontsize_col = 10)
dev.off()

# Heatmaps of Top Upregulated and Downregulated Genes

up_genes <- res_ordered %>%
  filter(significant=='Upregulated') %>%
  slice(1:30)

down_genes <- res_ordered %>%
  filter(significant=='Downregulated') %>%
  slice(1:30)

up_ids <- up_genes$gene_id
down_ids <- down_genes$gene_id

# Build matrices
mat_up <- rld_mat[up_ids, ]
mat_down <- rld_mat[down_ids, ]

# Replace IDs with gene names
rownames(mat_up) <- up_genes$gene_name
rownames(mat_down) <- down_genes$gene_name

samples_ordered <- samples[colnames(rld_mat), ]
colnames(mat_up) <- samples_ordered$sample_name
colnames(mat_down) <- samples_ordered$sample_name

png(snakemake@output[["heatmap_upregulated"]], width = 1200, height = 1200, res = 150)
pheatmap(mat_up,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Top Upregulated Genes (Treated vs Control)")
dev.off()

png(snakemake@output[["heatmap_downregulated"]], width = 1200, height = 1200, res = 150)

pheatmap(mat_down,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Top Downregulated Genes (Treated vs Control)")

dev.off()

sink()