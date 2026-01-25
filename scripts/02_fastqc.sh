#!/bin/bash

mkdir -p fastqc_reports
mkdir -p fastqc_logs

for fq in data/*.fastq.gz; 
do
    echo "Running FastQC on $fq"
    fastqc $fq -o fastqc_reports/ > fastqc_logs/"$sample".log 2>&1
done

echo "FastQC finished!"
