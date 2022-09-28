# Winter moth RNAseq
This folder contains all the scripts needed to reproduce the analysis of RNAseq data from winter moth embryos, published in Molecular Ecology https://doi.org/10.1111/mec.16705 : _Transcriptional regulation underlying the temperature response of embryonic development rate in the winter moth._

All used software versions are reported in the manuscript.

**NB: The raw RNAseq reads can be found on the European Nucleotide Archive (ENA) under accession no. PRJEB55675. All processed data used to run the analysis can be found on Dryad https://doi.org/10.5061/dryad.hx3ffbghd, including the final transcriptome incl. functional annotation, GO annotation table, final gene counts matrix, and phenotypic data.**

&nbsp;

## Authors
Natalie E. van Dis, ORCID ID: 0000-0002-9934-6751

&nbsp;

## Step1: RNAseq processing
### Script: ```1_Snakefile ```
File for ```snakemake``` in Linux to run the bioinformatics pipeline for RNAseq processing including QC, trimming, mapping, and transcript quantification.

The conda working environment listing all used software and software versions can be found in ```_src/env_WM_RNAseq.yml```.
Adapter sequences used for trimming can be found in ```_src/TruSeq3-PE-2.fa```

&nbsp;

## Step2: Functional annotation 
### Script: ```2_functional_annotation.md```
Command line code used to produce the functional annotation of the analyzed winter moth transcriptome.

&nbsp;

## Step3: Prepping the data
### Script: ```3a_annotion.R```
R script to process the functional annotation results from Step2 into an annotation table. Input files not deposited.

### Script: ```3b_preprocessing.R```
R script for data quality checks incl. filtering and PCA analysis.

&nbsp;

## Step4: Statistical analysis
### Script: ```4a_analysis_limma_allweeks.R```
R script for Differential expression analysis using R package ```limma```. See ```_src/env_limma.txt``` for R package versions.

### Script: ```4b_analysis_WGCNA_coexpr.R```
R script for Co-expression analysis using R package ```WGCNA```.

&nbsp;

## Step5: Visualizing the results 
### Script: ```5a_visual_limma_allweeks.R```
R script to visualize Differential expression analysis results with Venn diagrams and heatmaps.

### Script: ```5b_explore_DEGs.R```
R script to visualize patterns of individual differentially expressed genes (DEGs), check overlap between analyses, and produce MA and Volcano plots.

See ```_src/env_visualization.txt``` for R package versions.

&nbsp;

## Step6: Gene Ontology (GO) Overrepresentation analysis
### Script: ```6a_WGCNA-results_topGO.R```
R script for GO overrepresentation analysis on Co-expression analysis results with R package ```topGO```

### Script: ```6b_limma-results_topGO.R```
R script for GO overrepresentation analysis on Differential expression analysis results with R package ```topGO```. 

See ```_src/env_topGO.txt``` for R package versions.

&nbsp;

## Step7: Hierarchical clustering of GO results
### Script: ```7_topGO-results_ViSEAGO.R```
R script for hierarchical clustering of GO results with R package ```ViSEAGO```. See ```_src/env_ViSEAGO.txt``` for R package versions.

&nbsp;

## Step8: Study Descriptives
### Script: ```8_descriptives.R```

&nbsp;

## Reproducing manuscript Figures
### Script: ```figure1_samplesize.R```
### Script: ```figure2_modules.R```
### Script: ```figure3_Heatmaps.R```
### Script: ```figure4_DEGs.R```
### Script: ```figure5_VennDiagram.R```
