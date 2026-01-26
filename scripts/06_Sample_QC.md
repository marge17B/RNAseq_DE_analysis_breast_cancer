Sample QC
================
maria
2026-01-21

## Load libraries and data

``` r
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(pheatmap)

# Load tximport object and sample metadata
samples <- readRDS("../data/rds/samples.rds")
rld <- readRDS("../data/rds/rld.rds")
```

## Principal Component Analysis (PCA)

``` r
#plotPCA(rld, intgroup="condition")

pca <- plotPCA(rld, intgroup="condition", returnData=TRUE)
```

    ## using ntop=500 top features by variance

``` r
pca$sample_name <- samples$sample_name

#Inspect PCA data
pca
```

    ##                  PC1        PC2   group condition       name sample_name
    ## SRR7819990 -4.588670  1.9082375 control   control SRR7819990   control_1
    ## SRR7819991 -5.118626  2.9865056 control   control SRR7819991   control_2
    ## SRR7819992 -5.121976 -0.9451402 control   control SRR7819992   control_3
    ## SRR7819993  2.901784 -4.7033680 treated   treated SRR7819993   treated_1
    ## SRR7819994  3.215554 -4.7819655 treated   treated SRR7819994   treated_2
    ## SRR7819995  8.711934  5.5357307 treated   treated SRR7819995   treated_3

``` r
#Custom PCA plot with sample labels
pca_plot <- ggplot(pca, aes(PC1, PC2, color = condition, label = sample_name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 4) 

pca_plot
```

![](06_Sample_QC_files/figure-gfm/pca-1.png)<!-- -->

``` r
#save as png
ggsave("../results/figures_plots/PCA_plot.png", pca_plot, width = 6, height = 5, dpi = 150)
```

PCA reveals clear separation between control and treated samples. One
treated sample (treated_3 / SRR7819995) appears slightly separated from
other treated samples, suggesting a potential outlier that needs further
investigation.

## Inspect Cook’s distance for the treated_3 (SRR7819995)

``` r
dds_filter <- readRDS("../data/rds/dds_filter.rds")
cooks <- assays(dds_filter)[["cooks"]]
summary(cooks[, "SRR7819995"])
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##  0.00000  0.00494  0.02730  0.18535  0.12540 35.63046

``` r
summary(cooks[, "SRR7819993"])
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##  0.00000  0.00684  0.03611  0.16688  0.13768 32.18834

``` r
summary(cooks[, "SRR7819994"])
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ##  0.000000  0.003651  0.022141  0.148821  0.111012 30.227385

``` r
summary(cooks[, "SRR7819992"])
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##  0.00000  0.00895  0.04202  0.15755  0.13679 33.27046

The distribution of Cook’s distances for SRR7819995 was highly
comparable to the other samples. Although this sample appears slightly
separated from the other treated samples in the PCA, all technical
quality metrics, including mapping rate (~91%), and Cook’s distance,
were within normal ranges. Notably, SRR7819995 has slightly lower
library size (~40.7M reads vs ~50–52M in others), as confirmed by read
count summaries and MultiQC reports, which may contribute to this
separation. Thus the sample retained in the analysis.

## Heatmap

``` r
#matrix
rld_mat <- assay(rld)  
rld_cor <- cor(rld_mat)    ## cor() is a base R function
colnames(rld_cor) <- samples$sample_name
rownames(rld_cor) <- samples$sample_name
head(rld_cor) 
```

    ##           control_1 control_2 control_3 treated_1 treated_2 treated_3
    ## control_1 1.0000000 0.9996676 0.9994700 0.9991169 0.9991502 0.9989424
    ## control_2 0.9996676 1.0000000 0.9995061 0.9990914 0.9991178 0.9989193
    ## control_3 0.9994700 0.9995061 1.0000000 0.9991670 0.9991023 0.9988385
    ## treated_1 0.9991169 0.9990914 0.9991670 1.0000000 0.9995848 0.9993356
    ## treated_2 0.9991502 0.9991178 0.9991023 0.9995848 1.0000000 0.9994187
    ## treated_3 0.9989424 0.9989193 0.9988385 0.9993356 0.9994187 1.0000000

``` r
#annotation for sample conditions
annotation <- data.frame(condition = samples$condition)
rownames(annotation) <- samples$sample_name

#save
png("../results/figures_plots/heatmap.png", width = 800, height = 800, res = 150)
pheatmap(rld_cor, annotation_col = annotation)
```

![](06_Sample_QC_files/figure-gfm/heatmap-1.png)<!-- -->

``` r
dev.off()
```

    ## png 
    ##   3

Samples show high correlation (values close to 1) and cluster primarily
by condition (control vs treated). One treated sample (SRR7819995) shows
slightly lower correlation with other treated samples, consistent with
PCA, but overall the heatmap confirms data consistency and reliability
before differential expression analysis.
