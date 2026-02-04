Differential Expression Analysis
================
maria
2025-12-07

## Load libraries and data

``` r
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(apeglm)
library(ggrepel)
library(RColorBrewer)

# Load previously saved objects
txi <- readRDS("../data/rds/txi_salmon.rds")
samples <- readRDS("../data/rds/samples.rds")
dds <- readRDS("../data/rds/dds_filter.rds")
res_annotated <- readRDS("../data/rds/res_annotated.rds")
rld <- readRDS("../data/rds/rld.rds")
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.path = "results/figures_plots/",
  fig.width = 7,
  fig.height = 6
)

dir.create("results/figures_plots", recursive = TRUE, showWarnings = FALSE)


```

## Differential Expression Analysis

Run DESeq to estimate gene expression differences between treated and
control:

``` r
dds <- DESeq(dds)
```

    ## using pre-existing normalization factors

    ## estimating dispersions

    ## found already estimated dispersions, replacing these

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
res <- results(dds, name="condition_treated_vs_control")

summary(res)
```

    ## 
    ## out of 17415 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 2158, 12%
    ## LFC < 0 (down)     : 1912, 11%
    ## outliers [1]       : 30, 0.17%
    ## low counts [2]     : 3372, 19%
    ## (mean count < 19)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Differential expression analysis identified 4070 genes with significant
expression changes between treated and control samples (p-value \< 0.1).
A similar number of genes were upregulated (2158) and downregulated
(1912), indicating that the treatment affects gene expression in both
directions. Low-count genes and outliers were filtered to improve
statistical robustness.

## MA plot

MA plots visualize log2 fold change (M) versus mean expression (A) to
show differential expression:

``` r
#standard MA plot
plotMA(res)
```


``` r
# check LFC range and stats
range(res$log2FoldChange, na.rm = TRUE)
```

    ## [1] -23.72896  25.24171

``` r
summary(res$log2FoldChange)
```

    ##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
    ## -23.728957  -0.206925  -0.008187   0.016079   0.202263  25.241714

``` r
#range:-23.72896  25.24171
plotMA(res, ylim = c(-30, 30))
```


The x-axis shows the mean normalized read counts and the y-axis shows
log2 fold changes between treated and control conditions. Blue points
indicate significantly differentially expressed genes (p-value \< 0.1).
Genes with positive log2 fold changes are upregulated in the treated
condition, whereas genes with negative log2 fold changes are
downregulated. The raw MA plot displays extreme log2 fold changes for
some genes, with values exceeding ±20. These large fold changes are
probably driven by genes with low read counts and are not biologically
reliable, motivating the use of log2 fold change shrinkage.

## Shrink Log2 Fold Changes and MA Plot

``` r
res_shrink <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm")
```

    ## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    ##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    ##     sequence count data: removing the noise and preserving large differences.
    ##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
#Check range and summary
range(res_shrink$log2FoldChange, na.rm = TRUE)
```

    ## [1] -9.082294  7.250765

``` r
summary(res_shrink$log2FoldChange)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## -9.082294 -0.070536 -0.001729  0.012728  0.067262  7.250765

``` r
plotMA(res_shrink, ylim = c(-6, 6))
```

    ## png 
    ##   2

``` r
#plot
plotMA(res_shrink, ylim = c(-6, 6))
```

![Shrunken MA Plot](../results/figures_plots/MA_plot_shrunken.png)

After applying apeglm shrinkage, log2 fold changes are compressed into a
more biologically plausible range (−9 to +7), reducing the influence of
low-count genes.

## Gene Annotation Using GTF

``` r
#Read gtf
gtf <- read.table("../data/references/Homo_sapiens.GRCh38.115.gtf", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gtf_gene <- gtf[gtf$V3 == "gene", ]

#Parse attribute column
gtf_parsed <- gtf_gene %>%
  mutate(id = row_number()) %>%                 
  separate_rows(V9, sep="; ") %>%                
  separate(V9, c("name","value"), sep=" ") %>%   
  mutate(value = gsub(";","",value)) %>%         
  spread(name, value)
head(gtf_parsed)
```

    ## # A tibble: 6 × 14
    ##   V1    V2           V3        V4     V5 V6    V7    V8       id gene_biotype gene_id
    ##   <chr> <chr>        <chr>  <int>  <int> <chr> <chr> <chr> <int> <chr>        <chr>  
    ## 1 1     ensembl_hav… gene  3.07e6 3.44e6 .     +     .         1 protein_cod… ENSG00…
    ## 2 1     havana       gene  5.30e6 5.31e6 .     -     .         2 lncRNA       ENSG00…
    ## 3 1     havana       gene  5.49e6 5.49e6 .     +     .         3 lncRNA       ENSG00…
    ## 4 1     havana       gene  4.18e6 4.18e6 .     -     .         4 processed_p… ENSG00…
    ## 5 1     havana       gene  4.57e6 4.59e6 .     +     .         5 lncRNA       ENSG00…
    ## 6 1     havana       gene  5.09e6 5.09e6 .     -     .         6 lncRNA       ENSG00…
    ## # ℹ 3 more variables: gene_name <chr>, gene_source <chr>, gene_version <chr>

``` r
#View(gtf_parsed)

gtf_gene_names_id <- gtf_parsed %>% select(gene_id, gene_name)

#Convert DESeq2 results to dataframe
res_df <- as.data.frame(res)

res_df <- res_df %>%
  tibble::rownames_to_column(var = "gene_id")

#Merge
res_annotated <- merge(res_df, gtf_gene_names_id, by = "gene_id", all.x = TRUE)

# Save as RDS 
saveRDS(res_annotated,"../data/res_annotated.rds")
```

Gene annotation was performed by merging DESeq2 results with Ensembl
gene names extracted from the GTF file.

## Volcano plot

``` r
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

volcano_plot
```

    ## Warning: Removed 3402 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
ggsave(
  filename = "../results/figures_plots/volcano_plot.png",
  plot = volcano_plot,
  width = 8,
  height = 7,
  dpi = 300
)
```

    ## Warning: Removed 3402 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

The volcano plot shows the differential gene expression between treated
and control samples. Each point represents a gene. The x-axis shows the
log2 fold change, and the y-axis shows −log10(padj). Vertical dashed
lines indicate log2 fold change thresholds (±0.6), and the horizontal
dashed line indicates the significance threshold (padj = 0.05).
Upregulated genes are shown in yellow, downregulated genes in blue, and
non-significant genes in grey. Selected highly significant genes are
labeled.

## Differentially Expressed Genes

``` r
upregulated <- sum(res_annotated$significant == "Upregulated")
downregulated <- sum(res_annotated$significant == "Downregulated")
total_DE <- upregulated + downregulated

print(paste("Upregulated genes:", upregulated)) 
```

    ## [1] "Upregulated genes: 504"

``` r
print(paste("Downregulated genes:", downregulated)) 
```

    ## [1] "Downregulated genes: 230"

``` r
print(paste("Total DE genes:", total_DE)) 
```

    ## [1] "Total DE genes: 734"

## Top Differentially Expressed Genes

``` r
res_ordered <- res_annotated[order(res_annotated$padj), ]

#top 50 genes
top50_genes <- res_ordered %>%
  slice(1:50)

top50_genes
```

    ##            gene_id   baseMean log2FoldChange      lfcSE      stat        pvalue
    ## 1  ENSG00000175334  6504.9549      1.6917866 0.06565098  25.76940 1.954060e-146
    ## 2  ENSG00000163041  7956.5175      1.6189332 0.06406308  25.27092 6.670901e-141
    ## 3  ENSG00000196396  6606.8813      1.1469409 0.05031878  22.79350 5.319453e-115
    ## 4  ENSG00000105976  9462.9209      1.5158488 0.07064756  21.45649 3.971736e-102
    ## 5  ENSG00000128595 22835.1821      1.4744238 0.06964289  21.17120  1.760107e-99
    ## 6  ENSG00000180398 20565.0669      0.8968957 0.04361431  20.56425  5.737786e-94
    ## 7  ENSG00000101384 11781.4601      1.3285257 0.06785554  19.57874  2.347753e-85
    ## 8  ENSG00000119720   856.2523     -4.1970624 0.22027062 -19.05412  6.073198e-81
    ## 9  ENSG00000292366  2709.8369      1.4259820 0.07552146  18.88181  1.609662e-79
    ## 10 ENSG00000143384 23405.6502      0.9800892 0.05274724  18.58086  4.590645e-77
    ## 11 ENSG00000213281  6952.2531      1.1951762 0.06488845  18.41894  9.260767e-76
    ## 12 ENSG00000117632 16741.2150      1.3438558 0.07360063  18.25875  1.762652e-74
    ## 13 ENSG00000145919  2664.9687      1.1776209 0.06676227  17.63902  1.235879e-69
    ## 14 ENSG00000075785  8927.7532      0.9949124 0.05729643  17.36430  1.537655e-67
    ## 15 ENSG00000102024 12234.2969      0.9221593 0.05796483  15.90895  5.493022e-57
    ## 16 ENSG00000148572  1597.9256      1.2725127 0.08115352  15.68031  2.062299e-55
    ## 17 ENSG00000134899  1255.2563      1.2174577 0.07769217  15.67028  2.415311e-55
    ## 18 ENSG00000059804  8560.2893      0.8385370 0.05378268  15.59121  8.353785e-55
    ## 19 ENSG00000164066  1280.9835     -1.5627176 0.10093920 -15.48177  4.606121e-54
    ## 20 ENSG00000112305  2605.5335      0.9832447 0.06363548  15.45120  7.404768e-54
    ## 21 ENSG00000079739  7155.5807      0.8892734 0.05806126  15.31612  5.967375e-53
    ## 22 ENSG00000110536  1245.9647      1.9957578 0.13046223  15.29759  7.933950e-53
    ## 23 ENSG00000182934 13196.3100      0.8598733 0.05704942  15.07243  2.459029e-51
    ## 24 ENSG00000144028 12329.5473     -0.8151907 0.05501807 -14.81678  1.141246e-49
    ## 25 ENSG00000065978 23877.8112      0.7887879 0.05339238  14.77342  2.174012e-49
    ## 26 ENSG00000125450  3212.8729     -0.9115559 0.06191537 -14.72261  4.614884e-49
    ## 27 ENSG00000137142   465.9554      1.8175319 0.12711998  14.29777  2.259438e-46
    ## 28 ENSG00000075643  3720.3736     -1.1974203 0.08517700 -14.05802  6.877969e-45
    ## 29 ENSG00000159140 14745.6109      0.6958688 0.04953232  14.04878  7.836875e-45
    ## 30 ENSG00000140943  5503.1328      1.0683617 0.07649460  13.96650  2.496010e-44
    ## 31 ENSG00000155660 22938.2124      0.7286613 0.05228331  13.93679  3.786066e-44
    ## 32 ENSG00000167193  4714.7107      0.9647368 0.06976877  13.82763  1.736459e-43
    ## 33 ENSG00000121236   980.1952     -1.4201526 0.10329077 -13.74908  5.158412e-43
    ## 34 ENSG00000198768  3189.8227      0.7608360 0.05627853  13.51912  1.206165e-41
    ## 35 ENSG00000111328  2480.2283      0.8558992 0.06385480  13.40383  5.741814e-41
    ## 36 ENSG00000198873   995.9021      1.3081846 0.09768603  13.39173  6.759099e-41
    ## 37 ENSG00000105447  2076.0842      0.9900628 0.07401327  13.37683  8.259749e-41
    ## 38 ENSG00000151376  2393.1527     -0.9016180 0.06765568 -13.32657  1.621914e-40
    ## 39 ENSG00000182768  3492.8517      0.9159036 0.06902967  13.26826  3.536966e-40
    ## 40 ENSG00000120805  4813.9412      0.9050211 0.06828281  13.25401  4.277232e-40
    ## 41 ENSG00000105058  4132.0951      0.9215197 0.06976776  13.20839  7.848716e-40
    ## 42 ENSG00000162692   745.1983      1.5343016 0.11821852  12.97852  1.619767e-38
    ## 43 ENSG00000153310  4141.2350      1.3996542 0.10878177  12.86662  6.936798e-38
    ## 44 ENSG00000080824 60765.7339      0.7967888 0.06207729  12.83543  1.038228e-37
    ## 45 ENSG00000117868 12304.2150      1.1678989 0.09120758  12.80485  1.540325e-37
    ## 46 ENSG00000109046  4815.0253      1.0107929 0.07908314  12.78139  2.082996e-37
    ## 47 ENSG00000116478  5811.1664     -0.7269764 0.05704965 -12.74287  3.415874e-37
    ## 48 ENSG00000165802  3568.1261      0.8162746 0.06445336  12.66458  9.291807e-37
    ## 49 ENSG00000286001   622.3256      1.3938411 0.11044847  12.61983  1.641727e-36
    ## 50 ENSG00000253276  1070.1885      0.9982102 0.07952831  12.55163  3.893903e-36
    ##             padj gene_name   significant
    ## 1  2.738224e-142     BANF1   Upregulated
    ## 2  4.673967e-137     H3-3A   Upregulated
    ## 3  2.484717e-111     PTPN1   Upregulated
    ## 4   1.391398e-98       MET   Upregulated
    ## 5   4.932876e-96      CALU   Upregulated
    ## 6   1.340060e-90     MCFD2   Upregulated
    ## 7   4.699866e-82      JAG1   Upregulated
    ## 8   1.063797e-77     NRDE2 Downregulated
    ## 9   2.506243e-76     VAMP7   Upregulated
    ## 10  6.432871e-74      MCL1   Upregulated
    ## 11  1.179738e-72      NRAS   Upregulated
    ## 12  2.058337e-71     STMN1   Upregulated
    ## 13  1.332183e-66      BOD1   Upregulated
    ## 14  1.539082e-64     RAB7A   Upregulated
    ## 15  5.131581e-54      PLS3   Upregulated
    ## 16  1.806187e-52     NRBF2   Upregulated
    ## 17  1.990926e-52     ERCC5   Upregulated
    ## 18  6.503422e-52    SLC2A3   Upregulated
    ## 19  3.397136e-51      INTU Downregulated
    ## 20  5.188151e-51     SMAP1   Upregulated
    ## 21  3.981944e-50      PGM1   Upregulated
    ## 22  5.053566e-50    PTPMT1   Upregulated
    ## 23  1.498190e-48     SRPRA   Upregulated
    ## 24  6.663450e-47  SNRNP200 Downregulated
    ## 25  1.218577e-46      YBX1   Upregulated
    ## 26  2.487245e-46     NUP85 Downregulated
    ## 27  1.172648e-43   IGFBPL1   Upregulated
    ## 28  3.442178e-42     MOCOS Downregulated
    ## 29  3.786832e-42       SON   Upregulated
    ## 30  1.165886e-41    MBTPS1   Upregulated
    ## 31  1.711424e-41     PDIA4   Upregulated
    ## 32  7.604064e-41       CRK   Upregulated
    ## 33  2.190449e-40     TRIM6 Downregulated
    ## 34  4.971172e-39   APCDD1L   Upregulated
    ## 35  2.298858e-38   CDK2AP1   Upregulated
    ## 36  2.630979e-38      GRK5   Upregulated
    ## 37  3.128212e-38     GRWD1   Upregulated
    ## 38  5.981021e-38       ME3 Downregulated
    ## 39  1.270859e-37      NGRN   Upregulated
    ## 40  1.498421e-37      ARL1   Upregulated
    ## 41  2.682538e-37    FAM32A   Upregulated
    ## 42  5.404237e-36     VCAM1   Upregulated
    ## 43  2.260590e-35     CYRIB   Upregulated
    ## 44  3.306519e-35  HSP90AA1   Upregulated
    ## 45  4.796571e-35     ESYT2   Upregulated
    ## 46  6.345440e-35      WSB1   Upregulated
    ## 47  1.018439e-34     HDAC1 Downregulated
    ## 48  2.712627e-34      NSMF   Upregulated
    ## 49  4.695003e-34      <NA>   Upregulated
    ## 50  1.091305e-33   CCDC71L   Upregulated

``` r
write.csv(top50_genes,
          "../results/tables/top50_DE_genes.csv")
```

## Heatmap for top genes

``` r
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

pheatmap <- pheatmap(mat,
         scale='row',
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10)
```

``` r
png("../results/figures_plots/heatmap_top50",width = 1200, height = 1200, res = 150)
pheatmap
dev.off()
```

    ## png 
    ##   2

Heatmap of the top 50 differentially expressed genes. Genes were
selected based on adjusted p-value \< 0.05 and log2 fold change \> 0.6
(upregulated) or \< -0.6 (downregulated). Gene expression values are
shown as z-scores, where red indicates higher expression and blue
indicates lower expression relative to each gene’s mean. Hierarchical
clustering was applied to both genes and samples.

## Heatmaps of Top Upregulated and Downregulated Genes

``` r
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

pheatmap(mat_up,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Top Upregulated Genes (Treated vs Control)")
```

``` r
pheatmap(mat_down,
          scale = "row",
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          annotation_col = annotation_col,
          show_rownames = TRUE,
          fontsize_row = 6,
          fontsize_col = 10,
          main = "Top Downregulated Genes (Treated vs Control)")
```


``` r
png("../results/figures_plots/heatmap_upregulated.png",
    width = 1200, height = 1200, res = 150)

dev.off()
```

    ## png 
    ##   2
