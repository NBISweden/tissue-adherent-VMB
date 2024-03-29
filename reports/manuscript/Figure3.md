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

``` r
# fig.asp=1, fig.height=7,
### A
#############################
# DIVERSITY + TAX ABUNDANCE #
#############################
figlabels <- letters

layout(matrix(c(1,1,1,1,1,1,
                2,2,2,2,2,2,
                3,3,6,7,4,4,
                3,3,6,5,4,4) , ncol= 6,byrow = T),
                               widths = c(1.2,1,1.83,1.7,.85,1.4),
                               heights = lcm(c(2,4.5,2.5,2.3))
                               # heights = c(.7,1.6,1,1.1,1,.8,.2)
       )

ASV <- list(#c("Luminal_gr","ASV_Luminal_normalized"),
             c("Tissue_gr","ASV_Tissue_normalized")
         )
div_list <- list()
for(i in ASV){
  x <- datasets_all_samples[[i[2]]]
  
  df <- na.omit(metadata[sample_use, c(i[1], "ID")])
  metadata <- metadata[df[,"ID"],]
  gr <- factor(setNames(metadata[[i[1]]],metadata$ID))
  group_annotation <- na.omit(group_annotation)

  par(mar=c(0,5,1.5,1.5)) #b,l,t,r
  temp <- t(t(2^x-1)/colSums(2^x-1))
  temp <- temp[, df[,"ID"]]
  
  shann <- vegan::diversity(t(temp),index = c("shannon") )
  simp <- vegan::diversity(t(temp),index = c("simpson") )
  invsimp <- vegan::diversity(t(temp),index = c("invsimpson") )
  divers <- as.matrix(rbind(shann,simp,invsimp) )
  divers[is.infinite(divers)] <- 0
  rownames(divers) <- c("Shannon","Simpson","Inv Simpson")
  
  barlist( data = divers[,order(group_annotation)],
         main = "", xlab="",
         genes = c("Shannon","Simpson","Inv Simpson"),
         draw_mean_lines=F,col = pal[gr[order(group_annotation)]])
  
  div <- tibble(ID=names(shann),shann, simp, invsimp ) %>% 
    rename_with(., ~paste0(c("ID","Shannon","Simpson","Inv Simpson")))
  div_list[[length(div_list) + 1]] <- div
  
  if(i[2] == "ASV_Tissue_normalized"){
    end <- table(pal[factor(group_annotation[sample_use])])
    end <- end[order(factor(names(end), levels = pal))]
    end <- map_dbl(cumsum(end), ~ (.x * (par("usr")[2])/length(colnames(divers)))-.2 )
    start <- c(0.5, map_dbl(end[1:4], ~.x+0.7))
    
    axis(1, at = c(start, end), label = F, pos =3.3 , xpd=T, col=NA, col.ticks="black")
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
  temp <- rbind(temp[c(lvl[1:16]),], Other = colSums(temp[!(row.names(temp) %in% c(lvl[1:16])), ]))
  
  barplot( temp, las=2 , border = NA , yaxs="i", xaxs="i",ylab="abundance",
           main=gsub("_gr","",i[1]),font.main=1,names.arg = rep("",ncol(temp)),
           col = taxa_pal[factor(rownames(temp),levels = lvl)])
  
  points( (1:ncol(temp))*1.2 -.5, rep(par("usr")[3],ncol(temp))-2, 
          bg="white", pch=15,cex=.7,xpd=T, 
          col= c("tomato","orange","#d7d7d7")[factor(metadata$BV_Diagnosis_v3)][order( group_annotation)])
  points( (1:ncol(temp))*1.2 -.5, rep(par("usr")[3],ncol(temp))-7,bg="white",
        col= c("#d7d7d7","tomato")[factor(metadata$HIVstatus)][order( group_annotation )],
        pch=15,cex=.7,xpd=T)
  text(par("usr")[2],par("usr")[3]-2,labels = "      BV",cex=.5,xpd=T)
  text(par("usr")[2],par("usr")[3]-7,labels = "       HIV",cex=.5,xpd=T) # pos=4,
  
 l <- c("BV:","Normal","Intermediate","BV")
  legend(x = 10, y = 111, #"bottom", inset = c(0, -0.5),
         bty = "n", horiz=TRUE, cex=0.7, xpd = T, 
         border = NA, #pch=15,
         legend = l, x.intersp=0.5, xjust=0,
         text.width = c(0, map_dbl(l, ~strwidth(.x)/2)) ,
         fill=c(NA, "#d7d7d7","orange","tomato"))
  l2 <- c("HIV:","HIV-","HIV+")
  legend(x = 70, y = 111, # bottom position: par("usr")[3]-5
         bty = "n", horiz=TRUE, cex=0.7, xpd = T,
         border = NA,
         legend = l2, x.intersp=0.5, 
         text.width =c(0, map_dbl(l2, ~strwidth(.x)/2)),
         fill=c(NA, "#d7d7d7", "tomato"))
}


### B
############################
# TOTAL RELATIVE ABUNDANCE #
############################
sample_use <- colnames(TRX)[colSums(TRX)!=0]
raw_list <- raw_data %>%
  map(., ~column_to_rownames(., var = "X1")) %>%          
  map(., ~dplyr::select(., any_of(sample_use))) 
            
r <- raw_list %>%
  map(., ~mutate(.,across(where(is.numeric), ~ ./sum(.))) ) %>%
  map(., ~(rowSums(.x)/length(colnames(.x)) )*100 ) %>%
  bind_rows(., .id="ID") %>%
    pivot_longer(-ID) %>% 
    pivot_wider(names_from = ID, values_from = value) %>%
    mutate_at(vars(name), ~replace(., ASV_Luminal<0.55, "Other")) %>%
    mutate(name = ifelse(.$name=="Not assigned", "Other", .$name)) %>%
    group_by(name) %>%
    summarize(across(where(is.numeric), ~sum(.))) %>%
    arrange(-ASV_Luminal) %>%
    dplyr::rename(Luminal="ASV_Luminal",Tissue="ASV_Tissue") %>%
    `row.names<-`(., NULL) %>% 
    column_to_rownames(var = "name") %>%
    as.matrix()

par(mar=c(1.3,5,2,.3)) #b,l,t,r
barplot( r, border = NA ,ylab="abundance", xlab ="", axisnames = F,
         col = taxa_pal[factor(rownames(r),levels = lvl)])
title(main = "Total relative abundance ", line = 0.5, cex.main = .9, font.main= 1)
mtext("Luminal        Tissue ", line = 0, cex=.6, side=1)

#add label
add_letter(figlabels[1]); figlabels <- figlabels[-1]


### C
#############################
# LUMINAL-TISSUE COMPARISON #
#############################
par(mar=c(0,1.1,2,2)) #b,l,t,r
x <- as.character(metadata[,"Luminal_gr"])
x[is.na(x)] <- "NA"
y <- as.character(metadata[,"Tissue_gr"])
y[is.na(y)] <- "NA"
x1 <- x [ (x!="NA") & (y!="NA") ]
y1 <- y [ (x!="NA") & (y!="NA") ]
df <- data.frame(x1,y1)
p <- plot_sankey(df, color_by = 1, 
             plot_labels = T, plot_weights = F, order_1_by = "NULL",
             order_2_by = "disentangled",use_w1 = T,use_w2 = T,
             gap2v = 0,smoothing = .1,pal= pal
             )
text(c(0,1),c(1,1),labels = c( "Luminal" , "Tissue" ),pos=3, xpd=T)

#add label
add_letter(figlabels[1]); figlabels <- figlabels[-1]


###########
# LEGENDS #
###########
par(mar=c(0,0,0,.5)) #b,l,t,r # par(mar=c(1,.3,0,.5)) #b,l,t,r
empty_plot()
legend(x = -0.04, 
       y = 1.04, 
       legend = c("L. crispatus/jensenii", "L. iners", "Gardnerella", "High diverse", "Other") ,
       xjust = 0,yjust = 1,title.adj = 0,  ncol =1,
       bty = "n",pch = 22,pt.cex = 2,cex = 1,xpd=T,title = "Study Groups",
       pt.bg = pal)

par(mar=c(0,.2,1.5,.4)) #b,l,t,r
empty_plot()
legend(x = par("usr")[c(1)],
         y = par("usr")[c(4)],
         legend = lvl[1:12],xjust = 0,yjust = 1,
         bty = "n",pch = 22,pt.cex = 2,cex = 1,xpd=T,
         pt.bg = taxa_pal)

par(mar=c(0,0,1.5,1)) #b,l,t,r
empty_plot()
legend(x = par("usr")[c(1)],
         y = par("usr")[c(4)],
         legend = lvl[13:17],xjust = 0,yjust = 1,
         bty = "n",pch = 22,pt.cex = 2,cex = 1,xpd=T,
         pt.bg = taxa_pal[13:17])
```

<img src="./Figures/Figure 3.jpeg" style="display: block; margin: auto;" />

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
