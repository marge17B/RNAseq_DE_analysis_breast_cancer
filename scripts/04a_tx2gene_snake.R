sink(snakemake@log[[1]], split=TRUE)
library(GenomicFeatures)

gtf_path <- snakemake@input[["gtf_path"]]

txdb <- makeTxDbFromGFF(gtf_path)

tx2gene <- select(txdb, keys=keys(txdb, keytype="TXNAME"), columns="GENEID", keytype="TXNAME")
head(tx2gene)

write.csv(tx2gene, snakemake@output[["tx2gene"]], row.names = FALSE)

print("tx2gene created successfully")