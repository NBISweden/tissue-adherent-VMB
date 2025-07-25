# Distinct cervical tissue-adherent and luminal microbiome communities correlate with mucosal host gene expression and protein levels in Kenyan sex workers

## General info

This repository contains the code related to the article "Distinct cervical tissue-adherent and luminal microbiome communities correlate with mucosal host gene expression and protein levels in Kenyan sex workers" [Edfeldt, Kaldhusdal et al (2023)](https://doi.org/10.1186/s40168-023-01502-4)

This project used multiple datasets: 
- Transcriptomics data (bulk mRNA-SEQ)
- Microbiota data (16S) 
- Protein data (bead-based affinity assay)
- Cytokine data

The analysis can be replicated by installing `conda` and running the individual rmarkdown scripts.\
The `WORKFLOW.md` file describes the input and output of each rmarkdown script.\
Raw counts for the transcriptomics data is downloaded from GEO and the raw counts for the microbiota is downloaded from supplement files in the rmarkdown script named `03_normalize_data.Rmd`. Normalized MFI protein values and the cytokine values are downloaded in `Figure_06.Rmd`.

## Data Availability Statement

**Microbiome sequencing data** for this study has been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB50325. The raw counts table with taxonomic annotation can be found in Suppl. Table 1

**Transcriptomics count data** can be accessed in the Gene Expression Omnibus public repository, SuperSeries ID GSE217237. The raw transcriptomic sequencing data cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo\@ki.se), and access will be granted if the request meets the requirements of the data policy.

**Protein data** is previously published in doi: [10.1371/journal.ppat.1010494](https://doi.org/10.1371/journal.ppat.1010494.s017)\
**Cytokine data** Can be found in Suppl. Table 14\
**Sociodemographic and clinical characteristics** Can be found in Suppl. Table 2 

Supplementary material can be downloaded [here](https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-023-01502-4/MediaObjects/40168_2023_1502_MOESM2_ESM.zip)

## Workflow

### Analysis scripts

0.  [00_dada2_script_run1](https://nbisweden.github.io/tissue-adherent-VMB/00_dada2_script_run1)
0.  [00_dada2_script_run2](https://nbisweden.github.io/tissue-adherent-VMB/00_dada2_script_run2)
1.  [01_taxonomy_improved](https://nbisweden.github.io/tissue-adherent-VMB/01_taxonomy_improved)
2.  [02_data_preprocessing](https://nbisweden.github.io/tissue-adherent-VMB/02_data_preprocessing)
3.  [03_normalize_data](https://nbisweden.github.io/tissue-adherent-VMB/03_normalize_data)
4.  [04_clustering](https://nbisweden.github.io/tissue-adherent-VMB/04_clustering)

### Manuscript figures
1. [Figure1.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Figure1)
2. [Figure2.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Figure2)
3. [Figure3.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Figure3)
4. [Figure4-5.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Figure4-5)
5. [Figure6.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Figure6)
6. [Suppl.Figure1.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Suppl.Figure1)
7. [Suppl.Figure2.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Suppl.Figure3)
8. [Suppl.Figure4-5.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Suppl.Figure4-5)
9. [Suppl.Figure6.Rmd](https://nbisweden.github.io/tissue-adherent-VMB/Suppl.Figure6)

## Dependencies

Project is created with:\
\* R version: 3.6.2\
\* RStudio version: 2022.07.1+554

## Setup Conda

1.  Clone the repo
2.  If not already installed download mini conda/conda
3.  In the terminal navigate to the project directory
4.  create a new environment:<br/>
    `conda env create -n tissue_adherent_VMB -f ./resources/environment_tissue_adherent_VMB.yml`
5.  Activate the environment:<br/>
    `conda activate tissue_adherent_VMB`
6.  Open your preferred text editor from terminal by :<br/> 
    `/path/to/rstudio &`
7.  Install niceRplots package from github:<br/> 
    `remotes::install_github("czarnewski/niceRplots",force=T)`

## Repo structure

    project
    │   README.md
    │  
    └───code
    │   │   differential_abundance.R
    │   │   enrichment_function.R 
    │   │   knit_function.R
    │   │   Picrust_pipe.sh
    └───data
    │ 
    └───reports
    │   │   00_Preprocessing.Rmd
    │   │   02_Analysis.Rmd
    │   │   ...
    │   └───manuscript
    │       │   Figure01.Rmd
    │       │   Figure02.Rmd
    │       │   ...
    │       └───Figures
    │           │  Figure 1.jpeg 
    │           │   ...
    └───resources
    │   └───KEGG_GO_database
    │   └───OptiVag_BVAB_database
    │   │   ...
    │ 
