# RNAseq-Functional-Enrichment-Analysis
Automated R workflow for identifying significant differentially expressed genes (DEGs) from DESeq2 results and performing GO (BP, CC, MF) and KEGG functional enrichment analysis with publication-ready bar plots and bubble plots.

## Overview

A reproducible **R workflow** for identifying significant differentially expressed genes (DEGs) from complete **DESeq2** results and performing **Gene Ontology (GO)** and **KEGG pathway enrichment analysis** using **g:Profiler**. The workflow automatically filters significant DEGs, separates upregulated and downregulated genes, performs functional enrichment analysis, and generates publication-quality visualizations and summary tables.

Although developed using rice RNA-seq data, the workflow can be adapted to any organism supported by **g:Profiler** by changing the organism identifier.

## Features

- Reads the complete DESeq2 differential expression results.
- Automatically identifies significant DEGs using:
  - Adjusted p-value (padj) < 0.05
  - |log2FoldChange| > 1
- Separates genes into:
  - All significant DEGs
  - Upregulated DEGs
  - Downregulated DEGs
- Uses all DESeq2-tested genes as the background for enrichment analysis.
- Performs functional enrichment using **g:Profiler**.
- Supports:
  - GO Biological Process (GO:BP)
  - GO Cellular Component (GO:CC)
  - GO Molecular Function (GO:MF)
  - KEGG Pathways
- Generates publication-quality:
  - Bar plots
  - Bubble plots
- Displays GO and KEGG identifiers together with term names.
- Exports enrichment tables in CSV format.
- Saves figures in PNG, PDF and SVG formats.
- Automatically records R session information for reproducibility.


## Repository Structure

```
RNAseq-functional-enrichment-analysis/
│
├── scripts/
│   └── DEG_functional_enrichment_analysis.R
│
├── results/
│  ├── figures/
|  ├── tables/       
│
├── README.md
├── LICENSE

## Input File

The workflow requires the **complete DESeq2 results table**, including both significant and non-significant genes.
 Required columns:

| Column | Description |
|---------|-------------|
| Gene_ID | Gene identifier (e.g., Os01g0884300) |
| log2FoldChange | Log2 fold change |
| padj | Adjusted p-value (False Discovery Rate) |

Example:

```csv
Gene_ID,baseMean,log2FoldChange,padj
Os01g0100100,250.52,2.35,0.00004
Os01g0100200,185.21,-1.82,0.00021
Os01g0100300,98.70,0.45,0.284

## Output

### DEG Tables

results/functional_enrichment/tables/
DEG_Summary.csv
All_DESeq2_Tested_Genes.csv
All_Significant_DEGs.csv
Upregulated_Significant_DEGs.csv
Downregulated_Significant_DEGs.csv

## Functional Enrichment Tables

Separate enrichment tables are generated for:

- GO Biological Process
- GO Cellular Component
- GO Molecular Function
- KEGG Pathways

for:

- All significant DEGs
- Upregulated DEGs
- Downregulated DEGs

Each table contains:

- GO/KEGG ID
- Term name
- Gene count
- Gene ratio
- Background ratio
- Adjusted p-value
- −log10(adjusted p-value)
- Genes associated with each enriched term

## Generated Figures

For each DEG group, the workflow generates:

### GO Biological Process

- Top 10 Bar Plot
- Top 10 Bubble Plot

### GO Cellular Component

- Top 10 Bar Plot
- Top 10 Bubble Plot

### GO Molecular Function

- Top 10 Bar Plot
- Top 10 Bubble Plot

### KEGG Pathways

- Top 10 Bar Plot
- Top 10 Bubble Plot

Example output:

```
All_Significant_DEGs_GO_BP_Top10_Barplot.png
All_Significant_DEGs_GO_BP_Top10_Bubbleplot.png

All_Significant_DEGs_GO_CC_Top10_Barplot.png
All_Significant_DEGs_GO_CC_Top10_Bubbleplot.png

All_Significant_DEGs_GO_MF_Top10_Barplot.png
All_Significant_DEGs_GO_MF_Top10_Bubbleplot.png

All_Significant_DEGs_KEGG_Top10_Barplot.png
All_Significant_DEGs_KEGG_Top10_Bubbleplot.png
```

The same set of figures is generated for upregulated and downregulated DEGs.


## R Packages

The workflow requires:

- gprofiler2
- ggplot2
- dplyr
- readr
- stringr
- forcats
- tibble
- scales
- svglite


## Usage

1. Clone this repository.

2. Place the complete DESeq2 results file in the **data** directory.

3. Run:

source("scripts/DEG_functional_enrichment_analysis.R")

The workflow will automatically:

- Read the DESeq2 results.
- Identify significant DEGs.
- Separate upregulated and downregulated genes.
- Perform GO and KEGG enrichment analysis.
- Generate publication-quality figures.
- Export enrichment tables.
- Save R session information.
