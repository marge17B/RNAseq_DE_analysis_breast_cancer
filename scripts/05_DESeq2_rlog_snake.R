sink(snakemake@log[[1]], split=TRUE)
library("readr")
library('DESeq2')

#load data
txi <- readRDS(snakemake@input[["txi"]])
samples <- readRDS(snakemake@input[["samples_rds"]])


# Create DESEq2 object
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition)

#Inspect dataset
head(dds)

#Inspect raw counts
counts(dds)

#Filter low-count genes
keep = rowSums(counts(dds)) >=10 
dds_filter = dds[keep,]

#Inspect filtered counts
counts(dds_filter)

#Run DESeq2
dds_filter = DESeq(dds_filter)

#Save 
saveRDS(dds_filter, snakemake@output[["dds_filter"]])

## rlog transformation
rld <- rlog(dds_filter, blind=TRUE)
rld$sample_name <- samples$sample_name

saveRDS(rld, snakemake@output[["rld"]])