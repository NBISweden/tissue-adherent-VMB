Figure 3. Tissue Microbiome Abundances
================



``` r
##################
# LOAD LIBRARIES #
##################
suppressWarnings({suppressMessages({suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(vegan)
  library(RColorBrewer)
  # remotes::install_github("czarnewski/niceRplots",force=T)
  library(niceRplots)
  library(openxlsx)
})  })  })

#############
# LODA DATA #
#############
datasets_all_samples <- readRDS("../../results/03_normalize_data_output/datasets_all_samples.RDS")
metadata <- read.csv("../../data/metadata.csv",row.names = 1, stringsAsFactors = F)
raw_data <- c("ASV_Luminal", "ASV_Tissue") %>% 
  set_names() %>%
  map(., ~read_csv(paste0("../../results/02_data_preprocessing_output/",.x,"_raw_counts.csv")))
TRX <- datasets_all_samples[["Tissue_RNAseq_V3_normalized"]]
sample_use <- intersect(colnames(raw_data$ASV_Tissue), colnames(TRX)[colSums(TRX)!=0])

m <- metadata %>%
  select(ID, Luminal_gr, Tissue_gr) %>%
  arrange(Luminal_gr) %>%
  group_by(Luminal_gr) %>%
  nest() %>%
  mutate(L = map(data, ~arrange(.x, Tissue_gr))) %>%
  unnest(L)

group_annotation <- factor(setNames(m$Luminal_gr,m$ID))
group_annotation <- factor(setNames(metadata$Luminal_gr,metadata$ID))
group_annotation <- factor(setNames(metadata$Tissue_gr,metadata$ID))

#################
# COLOR PALETTS #
#################
lvl <- c("L. iners", "Gardnerella", "L. crispatus/acidophilus","Prevotella", "Atopobium","Sneathia", "L. jensenii", "Megasphaera", "Streptococcus", "Anaerococcus", "BVAB2", "Escherichia/Shigella", "BVAB1", "Dialister", "Mycoplasma", "Bifidobacterium", "Other")
pal <- c( "#0072B2", "#009E73","#D55E00", "#CC79A7", "#E69F00", "#999999")
taxa_pal <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal(8,"Pastel1"),"grey90")
```

<img src="./Figures/Figure 3.png" style="display: block; margin: auto;" />

**Fig. 3 Identification of a distinct ectocervical tissue-adherent
microbiome.** Ectocervical tissue samples were assessed for presence of
a tissue-adherent microbiome. a Bar plots of alpha diversity indices and
taxonomy profiles for each individual tissue sample. Color-coded squares
above the stacked bar plots show bacterial vaginosis (BV, binned
Nugent’s scores) and HIV diagnosis, respectively. Gray: negative,
orange: intermediate, red: positive BV; Gray: HIV seronegative, red: HIV
seropositive. b Total relative abundance in the luminal and tissue
microbiome datasets. All taxa with a total relative abundance \< 0.55
are included in the “other” category. c Microbiome profile shift between
luminal and tissue
