RNA-seq Differential Expression Analysis (Breast Cancer Dataset)
================
maria
2025-12-16

## Project Overview

This project presents a reproducible RNA-seq analysis pipeline to
investigate the transcriptomic effects of NRDE2 depletion in human
breast cancer cells using publicly available data from BioProject
PRJNA490376. NRDE2 is known to be involved in RNA interference, gene
silencing,RNA splicing, and nuclear export; understanding how its loss
alters gene expression may provide insight into genomic stability and
cancer progression. The pipeline processes RNA-seq data from
NRDE2-depleted vs. control breast cancer cells, covering data import,
quality control, transcript quantification with Salmon, normalization,
differential expression analysis using DESeq2, and downstream
visualization. This repository is designed to be fully reproducible and
portfolio-ready, providing a transparent workflow for RNA-seq analysis
from raw data to biological interpretation.

*Important note:* Due to file size limitations, large files such as raw
FASTQ files, reference genomes, Salmon indices, quantification outputs
(quant.sf), and RDS objects are not stored in this repository. All
results and data can be reproduced by downloading the required data and
running the provided scripts.

## Data Description

This project uses publicly available RNA-seq data, obtained from the
European Nucleotide Archive (ENA).

- **BioProject:**
  [PRJNA490376](https://www.ebi.ac.uk/ena/browser/view/PRJNA490376)
- **Organism:** *Homo sapiens*
- **Cell Model:** Breast cancer cells
- **Platform:** Illumina NextSeq 500
- **Experimental Design:**
  - Control: 3 samples
  - si-NRDE2: 3 samples
- **Read Type:** single-end

## Pipeline Overview

### 1. Data acquisition

Raw sequencing files were downloaded from the European Nucleotide
Archive (ENA) in FASTQ format.

**Execution:**

``` bash
scripts/01_download_data.sh
```

Detailed sample metadata, including accession IDs, sample names and
condition labels, is provided in:

``` bash
data/samples_metadata.csv
```

### 2. Quality Control

**Tools:** `FastQC`, `MultiQC`

**Description:** Initial quality assessment of raw read quality was
performed using **FastQC** to evaluate per-base quality scores, GC
content, and potential adapter contamination. **MultiQC** was used to
aggregate individual FastQC reports into a single summary report across
all samples. This step is automated using a Bash script,
`scripts/fastqc.sh`, that processes all FASTQ files in parallel and
generates individual reports and execution logs.

**Execution:**

``` bash
# Run FastQC on all FASTQ files
bash scripts/02_fastqc.sh

# MultiQC summary report
multiqc results/fastqc_reports/ -o results/multiqc_report/
```

**Output Directory:** 

* FastQC HTML reports for each sample: 'results/fastqc_reports/*.html' (not included in repository)
* MultiQC summary report: 'results/multiqc_report/multiqc_report.html' (not included in repository)

### 3. Transcript Quantification

**Tools:** `Salmon`

**Description:** Transcript-level quantification was performed using
**Salmon**. Reads were mapped against the human reference transcriptome
(GRCh38) to estimate abundances (TPM) and counts (NumReads)

**Execution:**

``` bash
bash scripts/03_salmon_quant.sh
```

**Output Directory:**

``` bash
data/salmon_quant  (not included in repository)
```

### 4. Import of Transcript-Level Estimates & Gene-Level Summarization

**Tools:** `tximport`, `GenomicFeatures` (R)

**Description:** This step converts raw transcript quantification into a
gene-level count matrix ready for statistical analysis.

- **Transcript-to-Gene Mapping:** A *tx2gene* mapping file was generated
  from the Ensembl GRCh38.115 GTF to link transcript IDs to gene
  symbols.

- **Summarization:** Using `tximport`, salmon transcript-level abundance
  estimates were imported into R and summarized to the gene level using
  the tx2gene mapping.

**Execution:**

``` bash
Rscript scripts/04a_tx2gene.R
Rscript scripts/04b_tximport.R
```

**Output Directory:**

``` bash
data/tx2gene.csv
data/rds/txi_salmon.rds  (not included in repository)
```

### 5. DESeq2 Object Creation & rlog Transformation

**Tools:** `DESeq2` (R)

**Description:**  
This script creates the DESeq2 object and filters low-count genes
(rowSums < 10). It also performs an rlog (regularized log)
transformation to provide normalized counts for quality control
analysis.

**Execution:**

``` bash
Rscript scripts/05_DESeq2_rlog.R
```

**Outputs:**

``` bash
data/rds/dds_filter.rds   # DESeq2 object for differential expression (input for Step 7, not included in repository).
data/rds/rld.rds          # rlog-transformed counts for PCA & heatmaps (input for Step 6, not included in repository).
```

### 6. QC: PCA & Sample Correlation Heatmap

**Tools:** `DESeq2`, `ggplot2`, `ggrepel`, `pheatmap` (R)

**Description:** This step performs sample-level QC to assess sample
similarity and detect potential outliers. For the full QC report (code, plots and commentary), see: **`results/reports/Sample_QC.md`**

- **Principal Component Analysis (PCA) :** Visualizes the major sources
  of variation between samples.

- **Sample correlation heatmap:** Displays pairwise correlations of gene
  expression between all samples.

**Execution:**

    Rscript -e "rmarkdown::render('scripts/06_Sample_QC.Rmd')"

**Outputs:**

``` bash
results/figures_plots/PCA_plot.png
results/figures_plots/qc_heatmap.png
```

*Notes:*

- PCA highlighted a potential outlier samples (SRR7819995), which
  appears more distant from the other treated samples. To evaluate this,
  several QC metrics were examined. The sample showed a normal mapping
  rate (~91%), consistent with the other samples, but had slightly fewer
  total reads, which could explain its PCA position. The sample is
  therefore considered acceptable to keep in the analysis.
- Correlation heatmap confirms high reproducibility among replicates,
  with samples clustering primarily by treatment group.
  

### 7. Differential Expression Analysis

**Tools:** `DESeq2`, `ggplot2`, `apeglm`, `pheatmap` (R)

**Description:** This step performs differential gene expression
analysis between treated and control samples using DESeq2. Genes with an
adjusted p-value \< 0.05 and \|log2 fold change\| \> 0.6 were considered
significantly differentially expressed. Log2 fold changes were shrunk
using `apeglm` algorithm to reduce noise and stabilize effect size
estimates for genes with low read counts. Gene annotation was performed
by merging DESeq2 results with Ensembl gene names extracted from the GTF
file. For the full DE report (code, plots and commentary), see: **`results/reports/DE_analysis.md`**

**Results:**  
- Total genes analyzed: 17415  
- Differentially expressed genes (FDR \< 0.05, \|log2FC\| \> 0.6): 734  
- Upregulated: 504  
- Downregulated: 230

**Visualizations:**  
- Raw and shrunken MA plots  
- Volcano plot of differential expression results  
- Heatmaps of top 50 DE genes

*Note:* Detailed statistical results, figures, and plots are available
in the markdown report. 'results/07_DE_analysis.Rmd'

**Execution:**

    Rscript -e "rmarkdown::render('scripts/07_DE_analysis.Rmd')"

**Outputs:**

``` bash
results/figures_plots/MA_plot_all_range.png
results/figures_plots/MA_plot_shrunken.png
results/figures_plots/volcano_plot.png
results/figures_plots/top50_DE_genes.csv
results/figures_plots/heatmap_top50
results/figures_plots/heatmap_upregulated.png
results/figures_plots/heatmap_downregulated.png
```

## Repository Structure

    ├── data/                   # sample metadata and RDS objects salmon references
    │   ├── rds/                # RDS objects (not tracked in this repository)
    │   ├── samples_metadata.csv    
    │   ├── references/         # Reference genome, GTF file (not tracked)
    │   └── salmon_quant/       # Salmon transcript quantification outputs (not tracked)
    ├── scripts/                # Bash and R scripts
    ├── results/            
    │   ├── fastqc_reports/     # FastQC reports  (not tracked)
    │   ├── multiqc_reports/    # MultiQC report  (not tracked)
    │   ├── multiqc_quant/      # MultiQC quantification summaries  (not tracked)
    │   ├── figures_plots/      # Final PNGs
    │   └── tables/             # CSVs of DE genes
    └── README.md

## Software & Versions

- FastQC v0.12.1
- MultiQC v1.18
- Salmon v1.10.2
- R v4.3.3
- R packages:
  - DESeq2 v1.42.1
  - tximport v1.30.0
  - ggplot2 v4.0.1
  - GenomicFeatures v1.54.4
  - apeglm v1.24.0
  - pheatmap v1.0.13
  - ggrepel v0.9.6
  - tidyverse v2.0.0

## Reproducibility

The analysis is fully reproducible. All scripts are provided in the
`scripts/` directory and should be executed in the order below.

1.  Download raw FASTQ files from ENA
    (<https://www.ebi.ac.uk/ena/browser/view/PRJNA490376>) using wget.

``` bash
bash scripts/01_download_data.sh
```

2.  Run quality control:

``` bash
bash scripts/02_fastqc.sh
```

3.  Transcript quantification using Salmon:

``` bash
bash scripts/03_salmon_quant.sh
```

4a. Generate tx2gene mapping:

``` bash
Rscript scripts/04a_tx2gene.R
```

4b. Import and summarize counts:

``` bash
Rscript scripts/04b_tximport.R
```

6.  Sample-level QC (PCA & heatmap):

``` bash
Rscript -e "rmarkdown::render('scripts/06_Sample_QC.Rmd')"
```

7.  Differential expression analysis:

``` bash
Rscript -e "rmarkdown::render('scripts/07_DE_analysis.Rmd')"
```
