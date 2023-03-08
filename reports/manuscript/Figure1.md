Figure 1. Luminal Study Groups
================



``` r
##################
# LOAD LIBRARIES #
##################
suppressWarnings({suppressMessages({suppressPackageStartupMessages({
  library(tidyverse)
  library(magick)
  library(uwot)
  library(scales)
  library(igraph)
  library(RColorBrewer)
  #remotes::install_github("czarnewski/niceRplots",force=T)
  library(niceRplots)
  library(openxlsx)
})  })  })

#########
# PATHS #
#########
result_dir <- "./Suppl.Tbl/"
if( isFALSE(dir.exists(result_dir)) ) { dir.create(result_dir,recursive = TRUE) }

#############
# LODA DATA #
#############
fig1a <- magick::image_read("../../resources/Scematic figure.pdf", density = 900 )
datasets_all_samples <- readRDS("../../results/03_normalize_data_output/datasets_all_samples.RDS")
metadata <- read.csv("../../data/metadata.csv",row.names = 1)
group_annotation <- factor(setNames(metadata$Luminal_gr,metadata$ID))

SNN_participant <- read.csv("../../results/04_clustering_output/participant_SNN_graph.csv",row.names = 1)
bac_communities <- read.csv("../../results/04_clustering_output/bacterial_communities.csv",row.names = 1)

TRX <- datasets_all_samples[["Tissue_RNAseq_V3_normalized"]]
sample_use <- colnames(TRX)[colSums(TRX)!=0]

#################
# COLOR PALETTS #
#################
lvl <- c("L. iners", "Gardnerella", "L. crispatus/acidophilus","Prevotella", "Atopobium","Sneathia", "L. jensenii", "Megasphaera", "Streptococcus", "Anaerococcus", "BVAB2", "Escherichia/Shigella", "BVAB1", "Dialister", "Mycoplasma", "Bifidobacterium", "Other")

pal <- c( "#0072B2", "#009E73","#D55E00", "#CC79A7", "#E69F00", "#999999")
taxa_pal <- c(RColorBrewer::brewer.pal(8,"Pastel2"),RColorBrewer::brewer.pal(8,"Pastel1"),"grey90")
# scales::show_col(taxa_pal)
```

``` r
# , fig.height=6.2
# fig.height=7, fig.asp=1,
figlabels <- letters

layout(matrix(c(1,1,1,1,1,1,
                2,2,2,2,2,2,
                3,3,3,3,3,3,
                5,6,7,7,8,8,
                5,4,7,7,8,8,
                5,4,7,7,8,8) , ncol= 6,byrow = T),
                               widths = c(1.87,1.42,.7,.8,.8,.7,.8),
                               heights = lcm(c(4.7,2,4.5,2,.5)) # height in cm
                              # heights = c(1.4,.55,1.4,.61,.75,0)
       ) 
### A
#################################
# COHORT SCHEMATIC ILLUSTRATION #
#################################
par(mar=c(0.1,3.7,0.3,1)) #b,l,t,r
plot(c(0,1),c(0,1),axes=F,type="n",xlab="",ylab="",xaxs="i",yaxs="i")
rasterImage(fig1a, 0, 0, 1, 1)

#add label
add_letter(figlabels[1]); figlabels <- figlabels[-1]


### B
#############################
# DIVERSITY + TAX ABUNDANCE #
#############################

ASV <- list(c("Luminal_gr","ASV_Luminal_normalized"),
            c("Tissue_gr","ASV_Tissue_normalized")
         )
div_list <- list()
for(i in ASV){
  x <- datasets_all_samples[[i[2]]]
  gr <- factor(setNames(metadata[[i[1]]],metadata$ID))

  dim(x)
  x <- x[ rowSums(x>0) >= 3 ,  ]
  #x_ <- x[ rowSums(x>0) ,  ]
  dim(x)
  temp <- t(t(2^x-1)/colSums(2^x-1))

  shann <- vegan::diversity(t(temp),index = c("shannon") )
  simp <- vegan::diversity(t(temp),index = c("simpson") )
  invsimp <- vegan::diversity(t(temp),index = c("invsimpson") )
  divers <- as.matrix(rbind(shann,simp,invsimp) )
  divers[is.infinite(divers)] <- 0
  rownames(divers) <- c("Shannon","Simpson","Inv Simpson")
  
  div <- tibble(ID=names(shann),shann, simp, invsimp ) %>% rename_with(., ~paste0(c("ID","Shannon","Simpson","Inv Simpson")))
  div_list[[length(div_list) + 1]] <- div
  
  if(i[2] == "ASV_Tissue_normalized"){break}
  
  par(mar=c(0,5,1.5,1.5)) #b,l,t,r
  barlist( data = divers[,order(group_annotation)],
         main = "", xlab="", #bg="white", #cex=.6, 
         genes = c("Shannon","Simpson","Inv Simpson"),
         draw_mean_lines=F,col = pal[gr[order(group_annotation)]])
  
  if(i[2] == "ASV_Luminal_normalized"){
    end <- table(pal[factor(group_annotation[order(group_annotation)])])
    end <- end[order(factor(names(end), levels = pal))]
    end <- map_dbl(cumsum(end), ~ (.x * (par("usr")[2])/length(colnames(divers)))-.2 )
    start <- c(0.5, map_dbl(end[1:4], ~.x+0.7))
    
    axis(1, at = c(start, end), label = F, pos =3.3 , xpd=T, col="white",col.ticks="black")
    map2(start, end, ~lines(x=c(.y, .x), y= c(3.3, 3.3), xpd=T, cex=.8))
    t <- map2(start, end, ~(.x+.y)/2)
    map2(t,levels(gr), ~text(x=.x, y=3.8, .y,xpd=T))
  }
  
  #add label
  add_letter(figlabels[1]); figlabels <- figlabels[-1]

  par(mar=c(1.5,5,1,1.5)) #b,l,t,r
  temp[is.nan(temp)] <- 0

  o <- order( rowSums(temp) , decreasing = T)
  sample_o <- names(group_annotation)[ order( group_annotation )]
  temp <- temp[o,sample_o]*100
  temp <- rbind(temp[1:16,], Other = colSums(temp[17:nrow(temp),]) )

  barplot( temp, las=2 , border = NA , yaxs="i", xaxs="i",ylab="abundance",main=gsub("_gr","",i[1]),font.main=1,
           col = taxa_pal[factor(rownames(temp),levels = lvl)] ,names.arg = rep("",ncol(temp)))

  points( (1:ncol(temp))*1.2 -.5, rep(par("usr")[3],ncol(temp))-2, bg="white",
        col= c("tomato","orange","#d7d7d7")[as.numeric(metadata$BV_Diagnosis_v3)][order( group_annotation )],
        pch=15,cex=.7,xpd=T)
  text(par("usr")[2],par("usr")[3]-2,labels = "      BV",cex=.5,xpd=T)
  points( (1:ncol(temp))*1.2 -.5, rep(par("usr")[3],ncol(temp))-7,bg="white",
        col= c("#d7d7d7","tomato")[as.numeric(metadata$HIVstatus)][order( group_annotation )],
        pch=15,cex=.7,xpd=T)
  text(par("usr")[2],par("usr")[3]-8,labels = "       HIV",cex=.5,xpd=T)

  l <- c("BV:","Normal","Intermediate","BV")
  legend(x = 10, y = 111, #"bottom", inset = c(0, -0.5),
         bty = "n", horiz=TRUE, cex=0.7, xpd = T, 
         border = "white", 
         legend = l, x.intersp=0.5, xjust=0,
         text.width = c(0, map_dbl(l, ~strwidth(.x)/2)) ,
         fill=c(NA, "#d7d7d7","orange","tomato"))
  l2 <- c("HIV:","HIV-","HIV+")
  legend(x = 80, y = 111, # bottom position: par("usr")[3]-5
         bty = "n", horiz=TRUE, cex=0.7, xpd = T,
         border = "white",
         legend = l2, x.intersp=0.5, 
         text.width =c(0, map_dbl(l2, ~strwidth(.x)/2)),
         fill=c(NA, "#d7d7d7", "tomato"))
}


###########
# LEGENDS #
###########
par(mar=c(0,0,0,.5)) #b,l,t,r
empty_plot()
legend(x = -0.04, y = 1.04, 
       legend = c("L1 L. crispatus/jensenii", "L2 L. iners", "L3 Gardnerella", "L4 High diverse", "L5 Other"),
       xjust = 0,yjust = 1,title.adj = 0,  ncol =1,
       bty = "n",pch = 22,pt.cex = 2,cex = 1,xpd=T,title = "Study Groups",
       pt.bg = pal)


par(mar=c(0,4.5,0,.5)) #b,l,t,r
empty_plot()
legend(x = par("usr")[c(1)],
         y = par("usr")[c(4)],
         legend = lvl[1:12],xjust = 0,yjust = 1,
         bty = "n",pch = 22,pt.cex = 2,cex = 1,xpd=T,
         pt.bg = taxa_pal)


par(mar=c(0,0,0,1)) #b,l,t,r
empty_plot()
legend(x = par("usr")[c(1)],
         y = par("usr")[c(4)],
         legend = lvl[13:17],xjust = 0,yjust = 1,
         bty = "n",pch = 22,pt.cex = 2,cex = 1,xpd=T,
         pt.bg = taxa_pal[13:17])


### C + D
################
# STUDY GROUPS #
################

### Louvain ###
g <- graph_from_adjacency_matrix(as.matrix(SNN_participant), mode = "undirected",diag = F,weighted = T)
set.seed(1)
l <- layout_with_fr(g,niter=3000,start.temp=30)

par(mar=c(0,0,0,1)) #b,l,t,r
plot( g , vertex.label.cex=0.000001 , vertex.color = pal[factor(group_annotation)] ,
      edge.width=  ( E(g)$weight / max(E(g)$weight)) ,
      vertex.size=10, layout=l,
      edge.color=colorRampPalette(c("grey95","black"))(90) [ round( E(g)$weight / max(E(g)$weight) * 88 )+1 ] )
title(main = "Louvain", line = -1, cex.main = 1)

#add label
add_letter(figlabels[1]); figlabels <- figlabels[-1]

### UMAP ###
all_microbiome <- datasets_all_samples[["ASV_Luminal_normalized"]]
all_microbiome <- all_microbiome[ rowSums(all_microbiome>0) >= 3 ,  ]
cors <- cor(all_microbiome[,], method = "pearson")
adj <- (1-cors)/2
U <- uwot::umap(adj,n_neighbors = 30)
plot( g , vertex.label.cex=0.000001 , vertex.color = pal[factor(group_annotation)] , 
      edge.width=  ( E(g)$weight / max(E(g)$weight)) ,layout=U,
      vertex.size=10,
      edge.color=colorRampPalette(c("grey95","black"))(90) [ round( E(g)$weight / max(E(g)$weight) * 88 )+1 ] )
title(main = "UMAP", line = -1, cex.main = 1)

#add label
add_letter(figlabels[1]); figlabels <- figlabels[-1]
```

<img src="./Figures/Figure 1.png" style="display: block; margin: auto;" />

**Figure 1. Characterization of a highly diversemicrobiome in
cervicovaginal (luminal) samples from Kenyan sex workers.**
Cervicovaginal lavage (luminal) and ectocervical tissue study samples
were assessed by 16S rRNA sequencing, gene expression, and protein
profiling. **a** Schematic drawing depicting the sampling scheme and the
resulting omics datasets. **b** Bar plots of alpha diversity indices and
taxonomy profiles for each individual luminal sample. Color-coded
squares above the stacked bar plots indicate bacterial vaginosis (BV,
binned Nugent’s scores): Gray: negative, orange: intermediate, red:
positive; and HIV diagnosis: Gray: HIV seronegative, red: HIV
seropositive. Two12-SNN graphs were constructed using: **c** Louvain
community detection algorithm, and d Uniform Manifold Approximation
(UMAP). The two graphs were overlayed in color with the predefined
luminal study groups. The undirected edges are included in gray
connecting the nodes

``` r
# Total Relative Abundance 
raw <- c("ASV_Luminal", "ASV_Tissue") %>% set_names()
raw_list <- raw %>%
  map(., ~read_csv(paste0("../../results/02_data_preprocessing_output/",.x,"_raw_counts.csv")))  %>%
  map(., ~column_to_rownames(., var = "X1")) %>%          
  map(., ~dplyr::select(., any_of(sample_use))) 
            
r <- raw_list %>%
  map(., ~mutate(.,across(where(is.numeric), ~ ./sum(.))) ) %>%
  map(., ~(rowSums(.x)/length(colnames(.x)) )*100 ) %>%
  bind_rows(., .id="ID") %>%
    pivot_longer(-ID) %>% 
    pivot_wider(names_from = ID, values_from = value) %>%
    mutate_at(vars(name), ~replace(., ASV_Luminal<0.5514392, "Other")) %>%
    mutate(name = ifelse(.$name=="Not assigned", "Other", .$name)) %>%
    group_by(name) %>%
    summarize(across(where(is.numeric), ~sum(.))) %>%
    arrange(-ASV_Luminal) %>%
    dplyr::rename(Luminal="ASV_Luminal",Tissue="ASV_Tissue") %>%
    rename(Taxa="name") %>%
  mutate(across(where(is.numeric), ~round(.x, digits = 1)))
  
# Alfa Diversity Stats
div_stat <- div_list %>%
  set_names(map(ASV,1)) %>%
  map(., ~summarise_if(.x, is.numeric, list(median=median, IQR=IQR), na.rm = TRUE)) %>%
  map(., ~pivot_longer(.x, everything(), names_sep = "_", names_to = c(".value", "Stats"), )) %>%
  bind_rows(.id = "Sample type") %>%
  mutate("Sample type" = gsub("_gr", "", .$"Sample type")) %>%
  mutate(across(where(is.numeric), ~round(.x, digits = 2)))

# Bacterial Communities
bac_communities <- bac_communities %>% 
  mutate("bacterial communities"= factor(paste0("BC",sprintf("%02d",bac_communities$Bact_com)))) %>%
  select(-Bact_com) %>%
  arrange(`bacterial communities`)

pages <- list("Total Relative Abundance"=r ,"Alfa Diversity"=div_stat, "Bacterial Communities"=bac_communities)
write.xlsx(pages, file=paste0(result_dir,"Suppl.Tbl.01",".xlsx"))
alpa_div <- set_names(div_list, c("Luminal","Tissue")) %>%
  map2(., c("Luminal_gr", "Tissue_gr"), ~left_join(.x, select(metadata, ID, .y), by="ID")) %>%
  map2(., c("Luminal_gr", "Tissue_gr"), ~select(.x, ID, .y, everything()))
write.xlsx(alpa_div, file=paste0(result_dir,"alpha_diversity",".xlsx"))
```
