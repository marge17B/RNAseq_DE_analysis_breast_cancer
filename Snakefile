SAMPLES = [
    "SRR7819990",
    "SRR7819991",
    "SRR7819992",
    "SRR7819993",
    "SRR7819994",
    "SRR7819995"
]

rule all:
    input:
        expand("results/fastqc_reports/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/fastqc_reports/{sample}_fastqc.zip", sample=SAMPLES),
        "results/figures_plots/PCA_plot_snake.png",
        "results/figures_plots/QC_heatmap.png",
        "data/samples_metadata.csv",
        "results/figures_plots/MA_plot_all_range.png",
        "results/figures_plots/MA_plot_shrunken.png",
        "results/figures_plots/volcano_plot.png",
        "results/tables/top50_DE_genes.csv",
        "results/figures_plots/heatmap_top50.png",
        "results/figures_plots/heatmap_upregulated.png",
        "results/figures_plots/heatmap_downregulated.png"

rule download_fastq:
    output:
        expand("data/{sample}.fastq.gz", sample=SAMPLES)
    log:
        "logs/logs_download.log"
    shell:
        "bash scripts/01_download_data_snake.sh &> {log}"
		
rule fastqc:
    input:
        expand("data/{sample}.fastq.gz", sample=SAMPLES)
    output:
        expand("results/fastqc_reports/{sample}_fastqc.html", sample=SAMPLES),
        expand("results/fastqc_reports/{sample}_fastqc.zip", sample=SAMPLES)
    log:
        "logs/logs_fastqc.log"
    shell:
        "bash  scripts/02_fastqc_snake.sh &> {log}"
		
rule quant:
    input:
        fastq=expand("data/{sample}.fastq.gz", sample=SAMPLES),
        index="data/references/Homo_index"
    output:
        expand("data/salmon_quant/{sample}/quant.sf", sample=SAMPLES)
    log:
        "logs/logs_quant.log"
    shell:
        "bash scripts/03_salmon_quant_snake.sh &> {log}"

rule tx2gene:
    input:
        gtf_path="data/references/Homo_sapiens.GRCh38.115.gtf"
    output:
        tx2gene="data/tx2gene.csv"
    log:
        "logs/tx2gene.log"
    script:
        "scripts/04a_tx2gene.R"

rule tximport:
    input:
        tx2gene="data/tx2gene.csv",
        salmon_quant="data/salmon_quant"
    output:
        samples_rds="data/rds/samples.rds",
        samples_metadata="data/samples_metadata.csv",
        txi="data/rds/txi_salmon.rds"
    log:
        "logs/tximport.log"
    script:
        "scripts/04b_tximport_snake.R"

rule DESEq2_rlog:
    input:
        txi="data/rds/txi_salmon.rds",
        samples_rds="data/rds/samples.rds"
    output:
        dds_filter="data/rds/dds_filter.rds",
        rld="data/rds/rld.csv"
    log:
        "logs/DESEq2_rlog.log"
    script:
        "scripts/05_DESeq2_rlog_snake.R"

rule Sample_QC:
    input:
        rld="data/rds/rld.csv",
        dds_filter="data/rds/dds_filter.rds",
        samples_rds="data/rds/samples.rds"
    output:
        PCA_plot="results/figures_plots/PCA_plot.png",
        QC_heatmap="results/figures_plots/QC_heatmap.png"
    log:
        "logs/Sample_QC.log"
    script:
        "scripts/06_Sample_QC_snake.R"

rule DE_analysis:
    input:	
        rld="data/rds/rld.csv",
        txi="data/rds/txi_salmon.rds",
        samples_rds="data/rds/samples.rds",
        dds_filter="data/rds/dds_filter.rds",
        gtf_path="data/references/Homo_sapiens.GRCh38.115.gtf"
    output:
        MA_plot_all_range="results/figures_plots/MA_plot_all_range.png",
        MA_plot_shrunken="results/figures_plots/MA_plot_shrunken.png",
        res_annotated="data/res_annotated.rds",
        volcano_plot="results/figures_plots/volcano_plot.png",
        top50_DE_genes="results/tables/top50_DE_genes.csv",
        heatmap_top50="results/figures_plots/heatmap_top50.png",
        heatmap_upregulated="results/figures_plots/heatmap_upregulated.png",
        heatmap_downregulated="results/figures_plots/heatmap_downregulated.png"
    log:
        "logs/DE_analysis.log"
    script:
        "scripts/07_DE_analysis_snake.R"
