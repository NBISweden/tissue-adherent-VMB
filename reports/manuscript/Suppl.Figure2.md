Suppl. Figure 2. Differential bacterial abundance across Tissue and
Luminal microbiome
================



``` r
##################
# LOAD LIBRARIES #
##################
suppressWarnings({suppressMessages({suppressPackageStartupMessages({
  #remotes::install_github("czarnewski/niceRplots",force=T)
  library(niceRplots)
  library(openxlsx)
})  })  })

#############
# LODA DATA #
#############
datasets_all_samples <- readRDS("../../results/03_normalize_data_output/datasets_all_samples.RDS")
metadata <- read.csv("../../data/metadata.csv", row.names = 1, stringsAsFactors = F)

# setwd("/Users/vilkal/work/Brolidens_work/Projects/Gabriella_repo/reports/rmarkdown/manuscript")
# datasets_all_samples <- readRDS("../../../results/datasets_all_samples.RDS")
# metadata <- read.csv2("../../../results/metadata_integration.csv",row.names = 1, stringsAsFactors = F)
group_annotation <- factor(setNames(metadata$Luminal_gr,metadata$ID))

#Graph construction  
fct <- .25 # FC threshold for differential expression
pvt <- 0.01 # Pvalue threshold for differential expression
min_pct <- .05 # minimun level of detected bacteria in each sample group

#################
# COLOR PALETTS #
#################
pal <- c( "#0072B2", "#009E73","#D55E00", "#CC79A7", "#E69F00", "#999999")
```

### Dot plot and bar plot

``` r
all_microbiome <- cbind(datasets_all_samples[["ASV_Tissue_normalized"]],
                        datasets_all_samples[["ASV_Luminal_normalized"]])
datasets <- factor(c( rep("Tissue",ncol(datasets_all_samples[["ASV_Tissue_normalized"]])) , rep("Luminal",ncol(datasets_all_samples[["ASV_Luminal_normalized"]])) ))
datasets <- datasets[ colSums(all_microbiome)>0 ]
all_microbiome <- all_microbiome[,colSums(all_microbiome)>0]

NN <- min(table(datasets))


res <- data.frame( matrix(0,nrow = 1,ncol = 6) )
for(i in levels(datasets)){
  for(j in rownames(all_microbiome) ){
    
    a <- all_microbiome[j,datasets == i]
    b <- all_microbiome[j,datasets != i]
    
    # calc % of samples that have a presence of bacteria j 
    perc1 <- sum(all_microbiome[j,datasets == i]>0) / sum(datasets == i) 
    perc2 <- sum(all_microbiome[j,datasets != i]>0) / sum(datasets != i)
    
    temp <- wilcox.test(x=a, y=b)
    
    fc <- log2( (mean(a)+1e-3) / (mean(b)+1e-3) )
    
    res <- rbind(res, setNames(c(j,i,fc,perc1,perc2,unlist(temp)[2] ),
                          c("bacteria","cluster","fc","perc.1","perc.2","pvalue")) )
    colnames(res) <- c("bacteria","cluster","fc","perc.1","perc.2","pvalue")
  }
}
res <- res[-1,]
res$pvalue <- as.numeric(res$pvalue)
res$perc.1 <- as.numeric(res$perc.1)
res$perc.2 <- as.numeric(res$perc.2)
res$perc.diff <- res$perc.1 - res$perc.2

res$fc <- as.numeric(res$fc)
res$FDR <- p.adjust(res$pvalue)
res <- res[order(res$pvalue),]
res <- res[abs(res$fc) > fct,]
res <- res[abs(res$pvalue) < pvt,]
res <- res[ (res$perc.1 > min_pct) | (res$perc.2 > min_pct) ,]
# dim(res)

#write.xlsx(res, file=paste0("../../../results/","DBA_across_Luminal_vs_Tissue",".xlsx"))

ord <- getcluster(data = all_microbiome, genes =  unique(as.character(res$bacteria)), clustering =  datasets)
ord <- ord[unique(as.character(res$bacteria))]

### A & B
####################################
# DIFFERENTIAL BACTERIAL ABUNDANCE #
####################################
figlabels <- letters
par(mfrow=c(1,2),mar=c(2,12,1,3.1))
plot_dots( all_microbiome, genes = names(sort(ord)) , clustering = datasets,
           show_grid = T,cex.main=1,font.main=1,
           cex.row = .8,cex.col = .8,srt = 0)
add_letter("a")
barlist( all_microbiome, names(sort(ord)) , clustering = datasets,
         show_grid = T,cex.main=1,font.main=1,
         cex.axis = .8,srt = 0)
add_letter("b")
```

<img src="./Suppl.Figures/Suppl.Fig.2.jpeg" style="display: block; margin: auto;" />

``` r
dev.off()
```

    ## null device 
    ##           1

**Suppl. Figure 2. Differential bacterial abundance across the luminal
and tissue microbiome datasets.** Differential bacterial abundance was
compared between the luminal and tissue-adherent microbiome data sets.
The results are shown as **a)** dot plots, and **b)** bar plots,
respectively. Bacteria with log2FC above 0.25 and p-value 0.01 (from the
Wilcoxon’s test) were considered significantly different and were sorted
by the highest expression. The colour scale indicates the difference in
total abundance between the datasets as a proportion, where “max” is the
highest abundance of the two datasets, and the other becomes a
proportion of this value. The size of the dots indicates the average
abundance of the given bacteria in the given data set.
