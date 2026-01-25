library(GenomicFeatures)

gtf_path <- "data/references/Homo_sapiens.GRCh38.115.gtf"

txdb <- makeTxDbFromGFF(gtf_path)

tx2gene <- select(txdb, keys=keys(txdb, keytype="TXNAME"), columns="GENEID", keytype="TXNAME")
head(tx2gene)

write.csv(tx2gene, "data/tx2gene.csv", row.names = FALSE)
