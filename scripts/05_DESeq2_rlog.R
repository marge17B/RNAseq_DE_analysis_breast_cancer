library("readr")
library('DESeq2')

#load data
txi <- readRDS("data/rds/txi_salmon.rds")
samples <- readRDS("data/rds/samples.rds")


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
saveRDS(dds_filter, "data/rds/dds_filter.rds")

## rlog transformation
rld <- rlog(dds_filter, blind=TRUE)
rld$sample_name <- samples$sample_name
rld

saveRDS(rld, "data/rds/rld.rds")