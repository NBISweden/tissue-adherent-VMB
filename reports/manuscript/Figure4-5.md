Figure 4. & 5. DGE’s and Enrichment
================



``` r
##################
# LOAD LIBRARIES #
##################
suppressWarnings({suppressMessages({suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(openxlsx)
  library(scales)
  library(fgsea)
  library(RColorBrewer)
  #remotes::install_github("czarnewski/niceRplots",force=T)
  library(niceRplots)
  library(readxl)
  library(enrichR)
  library(igraph)
  library(rafalib)
})  })  })

#########
# PATHS #
#########
result_dir <- "./Suppl.Tbl/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

DEG_files <- c("Suppl.Tbl.04 DEGs_Luminal.xlsx",
               "Suppl.Tbl.08 DEGs_Tissue.xlsx")
enrichR_files <- c("Suppl.Tbl.07 TF_PPI_Luminal.xlsx",
                   "Suppl.Tbl.11 TF_PPI_Tissue.xlsx")

#############
# LODA DATA #
#############
datasets_all_samples <- readRDS("../../results/03_normalize_data_output/datasets_all_samples.RDS")
metadata <- read.csv("../../data/metadata.csv",row.names = 1, stringsAsFactors = F)
sample_use <- metadata$ID
gr_name <- c("Luminal", "Tissue")

source("../../code/enrichment_function.R")
GO_database <- fgsea::gmtPathways("../../resources/KEGG_GO_database/c5.bp.v6.2.symbols.gmt.txt")
KEGG_database <- fgsea::gmtPathways("../../resources/KEGG_GO_database/c2.cp.kegg.v6.2.symbols.gmt.txt")

##################
# DEFINE CUTOFFS #
##################
dge_cutoff <- 0.01 # p-value
enrich_cutoff <- 0.1 # p-value
no_genes_cutoff <- 3
col_name <- paste0("genes_in(pValue<",dge_cutoff,")")

#################
# COLOR PALETTS #
#################
pal <- c( "#0072B2", "#009E73","#D55E00", "#CC79A7", "#E69F00", "#999999")
taxa_pal <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal(8,"Pastel1"),"grey90")
```

## Heatmap DGEs and enrichment

``` r
##################
# DEG's FUNCTION #
##################
#group <- "Tissue_gr"
#group <- "Luminal_gr"
edgR_dge.fun <- function(group, confound) {
  
  groups <- metadata[sample_use, group]
  n_gr <- length(na.omit(unique(groups)))
  
  if(any(is.na(groups))) {
    groups <- replace_na(groups, "X")
  }

  if(confound==FALSE) {
    c <- "_no_conf"
    design <- model.matrix(~groups, data=metadata[sample_use,]) # without confounders
  }else{
    c <- "_conf"
    design <- model.matrix(~groups+HIVstatus+Contraception,
                         data=metadata[sample_use,])
  }
  
  colnames(design) <- sub("groups","",colnames(design))
  colnames(design) <- sub("\\(Intercept\\)","Intercept",colnames(design))
  colnames(design) <- sub(" ","_",colnames(design))
  
  y <- DGEList(counts=TRX_counts[,sample_use], remove.zeros = T)
  y <- calcNormFactors(y,method = "TMM")
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTagwiseDisp(y,design)
  fit <- glmFit(y, design)
  
  lrt <- glmLRT(fit,coef=2:n_gr) # with intercept
  top <- topTags(lrt,adjust.method = "BH",n = "all",sort.by = "p.value")[[1]]
  colnames(top) <- sub(group,"",colnames(top))
  top <- cbind(LogFC.intercept=0,top) # with intercept
  top <-  rownames_to_column(top, var = "Genes")
  
return(list(design=design, y=y, fit=fit, lrt=lrt, top=top, c=c))
}

#######################
# ENRICHMENT FUNCTION #
#######################
# group_name <- "Luminal_gr"
# top <- dge_L$top
# y <- dge_L$y
# c <- dge_L$c
gene_mod_enrich.fun <- function(top, y, c, db, group_name, cuttoff){
  p_val <- pull(top, "PValue")
  top_dge <- (p_val < 0.01) #& (rowSums( abs(top[,grep("logFC",colnames(top))]) >= log2(1.5) ) >= 1)
  top_dge <- top$Genes[ top_dge ]
  
  top_TRX_counts <- edgeR::cpm(y,normalized.lib.sizes = T,log = T)[top_dge, ]
  top_TRX_counts <- t(apply(top_TRX_counts,1,function(x){scale(x,T,T)}))
  top_TRX_counts[top_TRX_counts > 5] <- 5
  colnames(top_TRX_counts) <- colnames(TRX_counts)
  
  
  h <- hclust( as.dist( (1- cor(t(top_TRX_counts)))/2 ), method = "ward.D2")
  cuttoff <- cuttoff[names(cuttoff) == db]
  gene_module <- cutree(h, h=cuttoff)
  
  if(db == "GO"){
  gmt_list <- GO_database}else{
  gmt_list <- KEGG_database}
  res_list <- lapply( unique( gene_module[h$order] ),function(x){
    temp <- compute_enrichment(genes = names(gene_module)[gene_module==x],
                              gmt_list = gmt_list,
                              min_terms_pathway = 10,
                              max_terms_pathway = 300,
                              min_overlap = 3,
                              sort_by_pvalue = T)
    temp <- temp[!grepl("REGULATION",rownames(temp)),]
    return(temp)
  } )
  
  return(list(res_list=res_list, top_TRX_counts=top_TRX_counts, 
              gene_module=gene_module, cuttoff=cuttoff, h=h, c=c))
}
```

``` r
###########################
# PAIRWISE DEG's FUNCTION #
###########################
pairwise_dge.fun <- function(group, fit, design) {
  groups <- metadata[sample_use, group]
  levels <- as.character(sort(unique(groups),na.last = NA))
  
  # without intercept:
  n <- combn(levels,2) 
  g <- map_chr(seq_along(1:ncol(n)), ~paste(n[1,.x], n[2,.x], sep = "-"))
  s <- seq_along(1:ncol(n)) %>% set_names(g)
  #name <- map_chr(seq_along(1:ncol(n)), ~paste(str_extract(n[1,.x], "^g."), str_extract(n[2,.x], "^g."), sep = "_"))
  
  # with intercept:
  g <- sub(paste0("^(",levels[1], ")-(.*)"), '-\\2', g)
  
  CONTRASTS <- makeContrasts( contrasts=g,
                              levels = design )
  
  #fit2 <- contrasts.fit(fit, CONTRASTS)
  #lrt_ <- glmLRT(fit2, coef = 2:5)
  lrt <- map(s, ~ glmLRT(fit, contrast = CONTRASTS[,.x]))
  top <- map(lrt, ~topTags(.x,adjust.method = "BH",n = "all",sort.by = "p.value")[[1]]) %>% 
         map(., ~rownames_to_column(.x, var = "Genes"))
  
  results <- map_df(lrt, ~decideTests(.x, method="separate", p.value = 0.05))
  summary <- bind_cols( "FDR<0.05"= c("Down","NotSig","Up"), summary(results))
  colnames(summary) <- c("FDR<0.05", g)
  top <- c(summary = list(summary), top)
  
return(list(lrt=lrt, top=top, CONTRASTS=CONTRASTS))
#return(list(lrt=lrt, top=top, CONTRASTS=CONTRASTS, results=results))
}
################################
# PAIRWISE ENRICHMENT FUNCTION #
################################
enrich_list <- function(x, gmt_list){
  temp <- compute_enrichment(genes = x,
                            gmt_list = gmt_list,
                            min_terms_pathway = 10,
                            max_terms_pathway = 300,
                            min_overlap = 3,
                            sort_by_pvalue = T)
  temp <- temp[!grepl("REGULATION",rownames(temp)),]
  temp <- rownames_to_column(temp, var = "Terms")
  return(temp)
} 
```

``` r
TRX <- datasets_all_samples[["Tissue_RNAseq_V3_normalized"]]
TRX_counts <- round(2^TRX - 1)
TRX_counts <- TRX_counts [ rowSums(TRX_counts>0)>= 1 , ]
TRX_counts <- TRX_counts[,sample_use]

#########
# DEG's #
#########
dge_L <- edgR_dge.fun("Luminal_gr", confound=T)
dge_T <- edgR_dge.fun("Tissue_gr", confound=T)


###############
# ENRICHMENT #
##############
db <- c("GO", "KEGG") %>% set_names()
enrich_L <- db %>%
  map(., ~gene_mod_enrich.fun(dge_L$top,dge_L$y,dge_L$c,
                              db = .x,
                              group_name = "Luminal_gr", 
                              cuttoff = c(GO=2.5, KEGG=2.5)#3
                              )) 
enrich_T <- db %>%
  map(., ~gene_mod_enrich.fun(dge_T$top,dge_T$y,dge_T$c,
                              db = .x,
                              group_name = "Tissue_gr",
                              cuttoff = c(GO=1.6, KEGG=1.6)
                              ))
```

``` r
##################
# PAIRWISE DEG's #
##################
pairwise_dge_L <- pairwise_dge.fun(group = "Luminal_gr", 
                                   fit=dge_L$fit, 
                                   design=dge_L$design )
pairwise_dge_T <- pairwise_dge.fun(group = "Tissue_gr",
                                   fit=dge_T$fit,
                                   design=dge_T$design)

###############################
# PAIRWISE ENRICHMENT GO/KEGG #
###############################
#info_logFC <- map(top_dfs[[1]], ~table(abs(.x$logFC) > 1))
#info_PValue <- map(top_dfs[[1]], ~table(.$PValue < 0.01))

l_top <- list(Luminal = pairwise_dge_L$top[-1], Tissue = pairwise_dge_T$top[-1])
# up and down regulated genes seperatly
t_split <- modify_depth(l_top, 2, ~ .x %>%
        dplyr::filter(., .$PValue < 0.01) %>%
        #dplyr::filter(., abs(.$logFC) < 1) %>%
        mutate(., reg = ifelse(.$logFC > 0, "up", "down")) %>%
        split(., .$reg)
  )

## GO enrichment 
pairwise_GO <- map_depth(t_split, 3, ~enrich_list(.x$Genes, GO_database)) %>%
  map_depth(., 2, ~bind_rows(.x, .id = "regulation")) 
## KEGG enrichment
pairwise_KEGG <- map_depth(t_split, 3, ~enrich_list(.x$Genes, KEGG_database)) %>%
  map_depth(., 2, ~bind_rows(.x, .id = "regulation")) 
```

``` r
# DEG's across and pairwise
pages <- list(list(Across=dge_L$top), pairwise=pairwise_dge_L$top) %>% flatten()
write.xlsx(pages, file=paste0(result_dir,DEG_files[1]))

pages <- list(list(Across=dge_T$top), pairwise=pairwise_dge_T$top) %>% flatten()
write.xlsx(pages, file=paste0(result_dir,DEG_files[2]))

#Files with filtered p-value < 0.01 
# Have to update the summary table, however the decideTests
# function does not provide support for p-value cutoffs only FDR
# therefore I need to write code to get the correct summary

p_L <- map(pairwise_dge_L$top[-1], ~filter(.x, `PValue`<0.01))
p_T <- map(pairwise_dge_T$top[-1], ~filter(.x, `PValue`<0.01))

pages <- list(list(Across=dge_L$top), pairwise=p_L) %>% flatten()
write.xlsx(pages, file=paste0(result_dir,DEG_files[1]))

pages <- list(list(Across=dge_T$top), pairwise=p_T) %>% flatten()
write.xlsx(pages, file=paste0(result_dir,DEG_files[2]))
```

``` r
# Enrichment Across
df_across <- enframe(list(Luminal=enrich_L, Tissue=enrich_T)) %>%  
  unnest_longer(value) %>%
  mutate(res = map(.$value, "res_list")) %>%
  unnest_longer(res) %>%
  mutate(res = map(res, ~rownames_to_column(.x, var = "Terms"))) %>%
  mutate(res = map(res, ~filter(.x, pvalue <= 0.01))) %>%
  group_by(name, value_id) %>%
  mutate(mod = row_number()) %>% ungroup() %>%
  mutate(res = set_names(.$res, paste0(.$value_id,"_gene_module_",.$mod))) %>%
  group_split(name) %>% set_names(., map(., ~unique(.x$name)))

write.xlsx(df_across[[1]]$res, file=paste0(result_dir,"Suppl.Tbl.05 Enrichment_Across_Luminal",".xlsx"))
write.xlsx(df_across[[2]]$res, file=paste0(result_dir,"Suppl.Tbl.09 Enrichment_Across_Tissue",".xlsx"))

# Enrichment pairwise
df_pairwise <- enframe(list(GO=pairwise_GO, KEGG=pairwise_KEGG)) %>%
  unnest_longer(value) %>%
  unnest_longer(value, names_repair ="unique") %>%
  mutate(value = map(value, ~filter(.x, FDR <= 0.05))) %>%
  mutate(value = set_names(.$value, paste0(.$name,"_",.$value_id...3))) %>%
  group_split(value_id...4) %>% set_names(., map(., ~unique(.x$value_id...4)))
  
write.xlsx(df_pairwise[[1]]$value, file=paste0(result_dir,"Suppl.Tbl.06 Enrichment_Pairwise_Luminal",".xlsx"))
write.xlsx(df_pairwise[[2]]$value, file=paste0(result_dir,"Suppl.Tbl.10 Enrichment_Pairwise_Tissue",".xlsx"))
```

<img src="./Figures/Figure 4.png" style="display: block; margin: auto;" />

<img src="./Figures/Figure 5.png" style="display: block; margin: auto;" />

**Figure. 4-5. Characterization of the host transcriptome as stratified
by the luminal microbiome study groups.** The luminal samples were
assessed for differential gene expression across the study groups. **a**
Differential gene expression analysis was applied across the five
luminal study groups. Significant DEGs (p-value \< 0.01) were divided
into six modules by hierarchical agglomerative clustering using inverse
Pearson’s correlation as distance measure and Ward’s method (“ward. D2”)
for linkage. Enrichment analysis was performed on each module using both
the KEGG and GO databases. The three most significant terms were
included in the heatmap. **b** Pairwise enrichment analysis of
protein-protein interactions of transcription factors (TF-PPI). Top-10
up- and down regulated transcription factors with p-value \< 0.01 were
included in the bar plots

``` r
metadata <- metadata %>% 
  mutate(L_T_groups = paste0(.$Luminal_gr, .$Tissue_gr)) %>%
  mutate(L_T_groups = ifelse(grepl("L2T2|L2T3", .$L_T_groups), .$L_T_groups, NA))

dge_LT <- edgR_dge.fun("L_T_groups", confound=T)
enrich_LT <- db %>%
  map(., ~gene_mod_enrich.fun(dge_LT$top,dge_LT$y,dge_LT$c,
                              db = .x,
                              group_name = "L_T_groups",
                              cuttoff = c(GO=2, KEGG=2) # c(GO=2, KEGG=2.7)
                              ))

# Plot the L2T2 vs L2T3 heatmap:
#pal <- c("#54ba63","#ad58c6","#d03e7a","#4cb6a9","#bb6fac","#cca442","#8fb53c","#6a7ecd","#63813e","#c46570","#b3743e","#d35238")
#plot_heatmap.fun(enrich_LT, "KEGG", "L_T_groups", "L_T_groups")

t <- dge_LT$top %>% filter(PValue <= 0.01)

df_across <- enframe(list(LtoT=enrich_LT)) %>%  
  unnest_longer(value) %>%
  mutate(res = map(.$value, "res_list")) %>%
  unnest_longer(res) %>%
  mutate(res = map(res, ~rownames_to_column(.x, var = "Terms"))) %>%
  mutate(res = map(res, ~filter(.x, pvalue <= 0.01))) %>%
  group_by(name, value_id) %>%
  mutate(mod = row_number()) %>% ungroup() %>%
  mutate(res = set_names(.$res, paste0(.$value_id,"_gene_module_",.$mod))) %>%
  group_split(name) %>% set_names(., map(., ~unique(.x$name)))


pages <- list(list("DEGs L2T2 vs L2T3"=t), df_across[[1]]$res) %>% flatten()
write.xlsx(pages, file=paste0(result_dir,"Suppl.Tbl.12 L2T2_L2T3_DEGs_and_enrichment.xlsx")) 
```

## igraph network plot

``` r
# 1. Change the groups of interest by assigning condition_use
# 2. Change the number of genes to be included in the network by changing n
#    look at tempD_i and tempU_i to see how many genes was significant
# 3. Adjust the size of the gene nodes and text in the graph by changing; square, circle and txt

# combinations
sheets <- map(DEG_files, ~excel_sheets(paste0(result_dir, .x))[-c(1:2)] ) %>%
  map(., ~set_names(.x))

# Set enrichment conditions
n <- 10  # maximum number of TF's to plot
# sizes:
square <- 10
cricle <- 8
txt <- .7
  
# Create a graph from the list of UP and DOWN regulated genes
i_graph <- function(temp_i, n, square, cricle){
  temp <- lapply( 1:n, function(x){ strsplit(as.character(temp_i$Genes[x]),split = ";")[[1]]  } )
  temp <- map(temp, purrr::discard, is.na) %>% compact() # removes NA
  names(temp) <- temp_i$`Transcription Factor`[1:length(temp)]
  temp2 <- lapply( names(temp) , function(x){ rep(x,length(temp[[x]])) })
  temp <- cbind(unlist(temp2),unlist(temp))
  g <- graph_from_data_frame(temp,directed = F)
  v_shape <- ifelse( names(V(g)) %in% unique(temp[,1]) , "square" , "circle" )
  v_color <- ifelse( names(V(g)) %in% unique(temp[,1]) , "honeydew" , "cornsilk" )
  v_size <- ifelse( names(V(g)) %in% unique(temp[,1]) , square , cricle ) 
  
  return(list(g=g, v_shape=v_shape, v_color=v_color, v_size=v_size))
  }


# Test parameters
# i <- "TF_PPI_Tissue"
# i <- "TF_TRRUST_enrichment"
# gr_name <- "Tissue"
# comp <- sheets

# Save all combination plots in one file
save_pdf.fun <- function(i, gr_name, comp){
  temp_i_name <- grep(i,list.files(paste0(result_dir),full.names=F),value=T)
  p <- map(comp, ~read_xlsx(paste0(result_dir, temp_i_name), sheet = .x)) 
  
  p_split <- p %>%
    imap(., ~.x %>%
          mutate(Regulation = paste0(.y, "_", .$Regulation)) %>%
          split(., .$Regulation)
          ) %>% compact() %>% flatten()
  
  obj_igraph <- imap(p_split, ~i_graph(.x, n, square, cricle))
  
  
  par(mar = c(0, 0, .7, 0))  #bottom, left, top, right
  pdf(file = paste0(result_dir, "Plots/", i, "_network_",".pdf", sep=""))
    p <- imap(obj_igraph,
         ~with(.x, 
                 plot.igraph(g, label=V(g), page = T,
                 main = paste0(str_replace(.y, "_", " ")," regulated genes"),
                 layout=layout_with_kk(g), vertex.label.color="black" ,
                 vertex.size=v_size, vertex.label.cex=txt, vertex.label.font=2, 
                 vertex.shape=v_shape, vertex.color=v_color, marigin = c(-1,-4,0,0),
                 vertex.frame.color="grey60", curved=F)
               
    ))
    dev.off()
}

# if(isFALSE(dir.exists(paste0(result_dir,"")))){dir.create(paste0("../",result.dir),recursive = TRUE)}
# save_pdf.fun("TF_PPI_Luminal", "Luminal", comp = sheets[[1]])
# save_pdf.fun("TF_PPI_Tissue", "Tissue", comp = sheets[[2]])
```
