Figure 6. Differential Protein expression
================



``` r
##################
# LOAD LIBRARIES #
##################
suppressWarnings({suppressMessages({suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(openxlsx)
  library(readxl)
  library(rstatix)
  #remotes::install_github("czarnewski/niceRplots",force=T)
  library(niceRplots)
})  })  })

#########
# PATHS #
#########
result_dir <- "./Suppl.Tbl/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

prot_url <- "https://doi.org/10.1371/journal.ppat.1010494.s017"
#cyt_url <- ""
#meta_url <- ""

############
# DOWNLOAD #
############
if(isFALSE(file.exists("../../data/Protein_norm_MFI.xlsx"))){download.file(prot_url, paste0("../../data/", "Protein_norm_MFI.xlsx"), method="auto")}
#if(isFALSE(file.exists("../../data/Cytokine_data.xlsx"))){download.file(cyt_url, paste0("../../data/", "Cytokine_data.xlsx"), method="auto")}
#if(isFALSE(file.exists("../../data/metadata.xlsx"))){download.file(meta_url, paste0("../../data/", "metadata.xlsx"), method="auto")}

#############
# LODA DATA #
#############
# datasets_all_samples <- readRDS("../../results/03_normalize_data_output/datasets_all_samples.RDS")
Protein_V3 <- read_xlsx("../../data/Protein_norm_MFI.xlsx", skip = 2) %>% rename(ID="Patient/Sample ID")
# cytokine_data <- read.csv("../../data/210125_Cytokines_visit3CVL.csv")
metadata <- read.csv("../../data/metadata.csv",row.names = 1, stringsAsFactors = F)

cytokine_data <- read_xlsx(paste0("/Users/vilkal/work/Brolidens_work/Projects/broliden_5325/data/", "Cytokine_data.xlsx"), na = "NA", sheet = "Cytokine_data") 
#################
# COLOR PALETTS #
#################
pal <- c( "#0072B2", "#009E73","#D55E00", "#CC79A7", "#E69F00", "#999999")
taxa_pal <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal(8,"Pastel1"),"grey90")
```

``` r
cytokine_data <- cytokine_data %>% 
  column_to_rownames(var = "ID") %>% 
  as.matrix()

cytokine_data <- t(cytokine_data)
cytokine_data[is.na(cytokine_data)] <- 0.01
cytokine_data <- log2(cytokine_data+1)
cytokine_data_log <- cytokine_data

cyto2 <- as.data.frame(t(cytokine_data))
cyto2$ID<-row.names(cyto2)

#select data and re-format to long format
cytogr.long <- metadata %>%
  select(ID, Luminal_gr, Tissue_gr) %>%
  left_join(.,cyto2, by="ID") %>%
  pivot_longer(c(-Luminal_gr, -Tissue_gr, -ID), names_to = "Cytokines", values_to = "value") 

test.krusk_dunn.fun <- function(df, var, group) {
  var <- enquo(var)
  
  stat.test.krusk <- df %>%
  group_by(!!(var)) %>%
  kruskal_test(., as.formula(paste("value", "~", group))) %>%
  adjust_pvalue(method = "BH") %>%
  select(-.y., -n, -df) %>%
  arrange(p) %>%
  rename(p.value = "p") %>%
  add_significance()

# select the significant proteins from kruskal
sig <- stat.test.krusk %>% filter(p.adj<0.05) %>% pull(!!var)

### Post hoc test
stat.test.dunn <- df %>%
  filter(!!var %in% sig) %>%
  group_by(!!var) %>%
  dunn_test(., as.formula(paste("value", "~", group))) %>%
  #dunn_test(value ~ !!group) %>%
  select(-.y.) %>%
  arrange(p) %>%
  rename(p.value = "p") %>%
  filter(p.adj<0.05)

  sig2 <- stat.test.dunn[stat.test.dunn$p.adj<0.05,]
  
 l <- list(stat.test.krusk, stat.test.dunn)
 n <-  c("Kruskal", "Dunn's") %>% paste0(., "_", group, sep = "")
 l <- set_names(l, n)

return(l)
}

n <- c("Luminal_gr", "Tissue_gr") %>% set_names(c("Luminal", "Tissue"))
stat <- imap(n, ~test.krusk_dunn.fun(cytogr.long, Cytokines, .x))
s <- flatten(stat)
```

``` r
names <- c("15 Luminal", "16 Tissue") %>% set_names()

# save TopTables pairwise
l_top_p <- list(pairwise_dpe_L$top, pairwise_dge_T$top)  %>% set_names(names)

# save TopTables across + cytokines
l_top_a <- list( list(Across = dpe_L$top), list(Across = dpe_T$top))  %>% set_names(names)
l <- list(l_top_a, stat, l_top_p)
p_list <- imap(names, ~flatten(map_depth(l, 1, .x))) 

imap(p_list, ~write.xlsx(.x,file=paste0(result_dir,"Suppl.Tbl.",.y,"_Proteins.xlsx"))) # Excel
write.xlsx(list("L2T2 vs L2T3"=dpe_LT$top),file=paste0(result_dir,"Suppl.Tbl.17"," L2T2_L2T3_Proteins.xlsx"))
```

``` r
# dev.new(height=10, width=6.6929133858, noRStudioGD = TRUE)
figlabels <- letters
####################################
# TOP DIFF. EXPRESSED PROT HEATMAP #
####################################
### A
# Luminal groups
t <- top_dpe.fun(dpe_L$top)
plot_heatmap.fun(t, "Luminal_gr", scale = T, layout = "multi")

### C
# Tissue groups
t <- top_dpe.fun(dpe_T$top)
par(mfg = c(2,1)) # adds the remaning plots to the previous layout
plot_heatmap.fun(t, "Tissue_gr", scale = T)


##############################
# SIG. CYTOKINES VIOLIN PLOT #
##############################
### B
# Luminal groups
g <- c("MIG", "IL1b", "IL12_p70", "MCP1", "IL1a", "IP10") 
#g <- c("MIG", "IL12_p70","IL1a", "IP10", "IL1b", "MCP1" )
lab <- c("MIG", "IL-1\u03b2", "IL-12-p70", "MCP-1", "IL-1\u03b1", "IP-10" )
#lab <- c("IP-10", "IL-1\u03b2", "IL-12-p70", "MCP-1","IL-1\u03b1", "MIG")
Luminal <- sig_stars.fun(s$Kruskal_Luminal_gr, s$`Dunn's_Luminal_gr`, "Luminal_gr", order = g[1:4])
rownames(Luminal$sig_df) <- lab
stars <- Luminal$stars_df 

par( mar=c(1.4,5,1.5,2)) #b,l,t,r 
violist(data=Luminal$sig_df,
          srt=45, y_padding = 2, smooth = .2,
          genes= lab[1:4],
          #cex.axis = 1.1,
          #add_ylims = T, #gives starting value of y-axis and max value of all data points
          col = pal, transparency = 50,
          clustering=c(Luminal$groupings))
segments(x0 = stars$line_x0, x1 = stars$line_x1, y0 = stars$y_, y1 = stars$y_) 
text( stars$x , stars$y_+.05,
        labels = stars$p.adj.signif, srt = 0, xpd = TRUE, cex=1)
title(main = "Luminal gr. Cytokines", line = 0.7, cex.main = 1.1)

add_letter("b")

Luminal <- sig_stars.fun(s$Kruskal_Luminal_gr, s$`Dunn's_Luminal_gr`, "Luminal_gr", order = g[5:6])
rownames(Luminal$sig_df) <- lab
stars <- Luminal$stars_df 
par(mar=c(1.4,5,1.5,2) ) #b,l,t,r
violist(data=Luminal$sig_df,
          srt=45, y_padding = 2, smooth = .2,
          genes= lab[5:6],
          #add_ylims = T, #gives starting value of y-axis and max value of all data points
          col = pal, transparency = 50,
          clustering=c(Luminal$groupings))
segments(x0 = stars$line_x0, x1 = stars$line_x1, y0 = stars$y_, y1 = stars$y_)
text( stars$x , stars$y_+.05,
        labels = stars$p.adj.signif, srt = 0, xpd = TRUE, cex=1)
title(main = "Luminal gr. Cytokines", line = 0.7, cex.main = 1.1)

### D
# Tissue groups
lab <- c("IL-1\u03b1", "IP-10")
Tissue <- sig_stars.fun(s$Kruskal_Tissue_gr, s$`Dunn's_Tissue_gr`, "Tissue_gr")
rownames(Tissue$sig_df) <- lab
stars <- Tissue$stars_df

#par(mar=c(20.4,6,1.5,3)) #b,l,t,r # ps = 20, point size
par(mar=c(1,5,1.5,5)) #b,l,t,r
violist(data=Tissue$sig_df,
          srt=45, y_padding = 2, smooth = .2,
          genes= rownames(Tissue$sig_df),
          col = pal, transparency = 50,
          # add_ylims = T,
          clustering=Tissue$groupings )
segments(x0 = stars$line_x0, x1 = stars$line_x1, y0 = stars$y_, y1 = stars$y_)
text( stars$x , stars$y_+.05,
        labels = stars$p.adj.signif, srt = 0, xpd = TRUE, cex=1)
title(main = "Tissue gr. Cytokines", line = 0.7, cex.main = 1.1)

add_letter("d")
```

<img src="./Figures/Figure 6.png" style="display: block; margin: auto;" />

**Figure 6. Characterization of the host protein profile as stratified
by the luminal and tissue-adherent microbiome study groups.** All study
samples were assessed for significant differences in protein levels
across the study groups. **a**,**b** The five luminal study groups and
**c**,**d** the five corresponding tissue groups. a, c Proteins with a
p-value \> 0.05 were omitted from the heatmap. Proteins with
significantly different levels across the groups were clustered by
hierarchical agglomerative clustering using inverse Pearson’s
correlation as distance measure and Ward’s method (“ward. D2”) for
linkage. Color-coded rows below the heatmaps show clinical diagnosis of
bacterial vaginosis (BV, binned Nugent’s scores) and HIV diagnosis,
respectively. Color coding for BV: Gray: negative, orange: intermediate,
red: positive; and for HIV: Gray: HIV seronegative, red: HIV
seropositive. b, d Violin plots of log2 transformed cytokine levels
across the luminal and tissue study groups (as indicated on the x-axis),
respectively. The asterisk and lines indicate statistically significant
results of Dunn’s test with Benjamini Hochberg’s correction analysis.
Adjusted p-values: \* \< 0.05, \*\* \< 0.01, \*\*\* \< 0.001, \*\*\*\*
\< 0.0001
