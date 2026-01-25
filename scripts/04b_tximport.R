library("tximport")
library("readr")

tx2gene <- read.csv(file= "data/tx2gene.csv")

## Sample metadata

samples <- data.frame(
  sample = c("SRR7819990","SRR7819991","SRR7819992", "SRR7819993","SRR7819994","SRR7819995"),
  condition = factor(c("control","control","control","treated","treated","treated")), #convert a character vector into a factor
  sample_name = c("control_1","control_2","control_3","treated_1","treated_2","treated_3")
)

rownames(samples) <- samples$sample

saveRDS(samples, file = "data/rds/samples.rds")
write.csv(samples, "data/samples_metadata.csv", row.names = FALSE)

## Import salmon data
files <- file.path( "data/salmon_quant", samples$sample, "quant.sf")
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
saveRDS(txi, file = "data/rds/txi_salmon.rds")

print('txi succesfully saved!')