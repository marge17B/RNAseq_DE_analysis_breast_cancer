sink(snakemake@log[[1]], split=TRUE)
library("tximport")
library("readr")

tx2gene <- read.csv(snakemake@input[["tx2gene"]])

# Sample metadata

samples <- data.frame(
  sample = c("SRR7819990","SRR7819991","SRR7819992", "SRR7819993","SRR7819994","SRR7819995"),
  condition = factor(c("control","control","control","treated","treated","treated")), #convert a character vector into a factor
  sample_name = c("control_1","control_2","control_3","treated_1","treated_2","treated_3")
)

rownames(samples) <- samples$sample

saveRDS(samples, snakemake@output[["samples_rds"]])
write.csv(samples, snakemake@output[["samples_metadata"]], row.names = FALSE)

# Import salmon data
files <- file.path(snakemake@input[["salmon_quant"]], samples$sample, "quant.sf")
names(files) <- samples$sample

#create txi
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

#inspect txi
names(txi)           # counts, abundance, length
dim(txi$counts)      # genes x samples
head(txi$counts)     # first rows of raw counts
head(txi$abundance)  # TPM values
head(txi$length)     # gene lengths
head(txi$countsFromAbundance)   

# Save the txi object
saveRDS(txi, snakemake@output[["txi"]])

print('txi succesfully saved!')