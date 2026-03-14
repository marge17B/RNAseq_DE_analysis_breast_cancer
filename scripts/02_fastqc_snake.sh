#!/bin/bash

mkdir -p results/fastqc_reports
mkdir -p results/fastqc_logs

for fq in data/*.fastq.gz; 
do
    echo "Running FastQC on $fq"
    fastqc $fq -o results/fastqc_reports/ > results/fastqc_logs/"$sample".log 2>&1
done

echo "FastQC finished!"