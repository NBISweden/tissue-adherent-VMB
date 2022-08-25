compute_differential_abundance <- function(x){
  
  NN <- min(table(datasets))
  
  res <- data.frame( matrix(0,nrow = 1,ncol = 6) )
  for(i in levels(datasets)){
    for(j in rownames(all_microbiome) ){
      
      set.seed(1)
      a <- c(sample(all_microbiome[j,datasets == i],NN))
      set.seed(1)
      b <- sample(all_microbiome[j,datasets != i],NN)
      
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
  # res <- res[order(res$pvalue),]
  res <- res[abs(res$fc) > fct,]
  res <- res[abs(res$pvalue) < 0.05,]
  res <- res[ (res$perc.1 > min_pct) | (res$perc.2 > min_pct) ,]
  dim(res)
}