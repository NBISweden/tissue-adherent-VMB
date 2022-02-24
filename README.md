# Spatial_DMPA analysis and figures

## Table of contents
* [General info](#general-info)
* [Dependencies](#dependencies)
* [Data Availability Statment](#data-availability-statment)
* [Repo description](#repo-description)
* [Omics integration](#omics-integration)
* [Setup](#setup)

## General info
Publication: "Distinct human cervical tissue-adherent and luminal microbiome communities correlate with mucosal host gene expression and protein levels"  
doi: [](link)  

This project used multiple omics data:  
- Transcriptomics data (bulk mRNA-SEQ)  
- Microbiome data (16S)  
- Protein data (bead-based affinity assay)  

## Data Availability Statment
The **raw microbiome sequencing data** for this study has been deposited
in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB50325.  

The **processed transcriptomics sequencing data** files can be accessed in the Gene Expression Omnibus public repository, accession ID GSE194276. 

The **raw transcriptomic sequencing data** and sociodemographic and clinical characteristics of the study participants cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo@ki.se), and access will be granted if the request meets the requirements of the data policy. **Protein/cytokine data** is available upon request. 
	
## Dependencies
Project is created with:
* R version: 3.6.2
* RStudio version: 1.1.456
* renv version: 0.15.2

## Repo description
- **src**  
  contains all the analysis script including all preprocessing steps
- **manuscript**  
  reproducible code for figures included in the manuscript
- **data**  
  Contains ...

```
project
│   README.md
│   renv.loc    
└───src
│   │   00_Preprocessing.Rmd
│   │   02_Analysis.Rmd
│   │   ...
│   └───manuscript
│       │   Figure01.Rmd
│       │   Figure02.Rmd
│       │   ...
└───bin
│   │   file021.txt
│   │   file022.txt
└───data
│   │   file021.txt
│   │   file022.txt
│
```

## Omics integration


## Setup
There are two ways to run this project:

Conda + renv

1. Clone the repo
2. If not allready installed download mini conda/conda
3. In the terminal navigate to the project directory
4. create a new enviroment:<br/>
  `conda env create -n Spatial_DMPA -f environment.yml`
5. Activate the enviroment:<br/>
  `conda activate Spatial_DMPA`
6. Open Rstudio:<br/>
  `rstudio& Spatial DMPA`
5. Install all packages specified by the lockfile:<br/>
  `renv::restore()`
  


