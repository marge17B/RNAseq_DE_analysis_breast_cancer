#!/bin/bash

for seq in data/*.fastq.gz
do
    sample=$(basename "$seq" .fastq.gz)            #basename "$seq": gets the file name with the suffix. by addinf  .fastq.gz i remve suffix
    echo "Running Salmon quantification on $sample"

    salmon quant \
        -i /home/maria/Documents/python_r/ngs/rna-breast/data/references/Homo_index/ \
        -l A \
        -r "$seq" \
        --validateMappings \
        --gcBias \
        -o salmon_quant/"$sample"

    echo "Finished $sample"
done

echo "All Salmon quantifications finished!"
