compute_enrichment <- function(genes , gmt_list , min_overlap = NULL, min_terms_pathway=NULL, max_terms_pathway=NULL, sort_by_pvalue=F){
  
  allgenes <- unique( c(unlist(gmt_list) , genes) )
  
  res <- t(sapply( gmt_list ,
                   genes=genes,
                   Total=length(allgenes),
                   inselection=length(genes),
                   function(x , genes,Total, inselection){
                     Overlap <- sum(x%in%genes)
                     inpathway <- length(x)
                     return( c(overlap=Overlap,
                               in_input_list=inselection,
                               in_pathway=inpathway,
                               all_gene=Total,
                               percent_overlap=Overlap/inselection,
                               overlap_coeficient=Overlap/min(inpathway,inselection),
                               jaccard_index=Overlap / (inpathway+inselection-Overlap),
                               pvalue=phyper(Overlap-1, 
                                             inpathway, 
                                             Total-inpathway, 
                                             inselection, 
                                             lower.tail= FALSE)) )
                   }))
  res <- as.data.frame(res)
  res$FDR <- p.adjust(res$pvalue,method = "BH" )
  res$genes <- sapply(gmt_list,genes=genes,function(x,genes) paste0(x[x%in%genes],collapse = ";") )
  
  if(sort_by_pvalue){ res <- res[order(res$pvalue),] }
  if(!is.null(min_overlap)){ res <- res[res$overlap >= min_overlap, ] }
  if(!is.null(min_terms_pathway)){ res <- res[res$in_pathway >= min_terms_pathway, ] }
  if(!is.null(max_terms_pathway)){ res <- res[res$in_pathway <= max_terms_pathway, ] }
  
  return(res)
}

# 
# mypar()
# sel <- c("grey80","red")[ (res$in_pathway > 1000) + 1 ]
# plot(res$jaccard_index, res$percent_overlap,cex=.2,col= sel )
# plot(res$jaccard_index, res$overlap_coeficient,cex=.2,col= sel )
# plot(res$percent_overlap, res$overlap_coeficient,cex=.2,col= sel )
# 
# 
# plot(res$jaccard_index, -log10(res$pvalue),cex=.2,col= sel )
# plot(res$percent_overlap, -log10(res$pvalue),cex=.2,col= sel )
# plot(res$overlap_coeficient, -log10(res$pvalue),cex=.2,col= sel )
# plot(res$in_pathway, -log10(res$pvalue),cex=.2,col= sel )
# plot(res$in_pathway, -log10(res$FDR),cex=.2,col= sel )
# 
# plot(res$in_pathway, -log10(res$FDR*res$in_pathway/res$all_gene),cex=.2,col= sel )
# plot(-log10(res$pvalue), -log10(res$FDR),cex=.2,col= sel )
# 
# 
# plot(res$in_pathway, res$jaccard_index,cex=.2,col= sel )
# plot(res$in_pathway, res$overlap_coeficient,cex=.2,col= sel )
# plot(res$in_pathway, res$percent_overlap,cex=.2,col= sel )
# plot(-log10(res$pvalue), -log10(res$pvalue*res$in_pathway ),cex=.2,col= sel )
# 
# lines(c(0,1),c(0,6))
# 
# res[order(res$percent_overlap),]

