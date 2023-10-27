# Coexpression Enrichment Analysis
# (C) Christian Forst (chris@santafe.edu) 2018-2023
# assume gn data frame from co-expression module file
# assume d data frame with DE genes
# assume gn is list of genes

# FisherTest
#input vector:=[popilation, population_hit, list, list_hit]
fisherTest = function(poplisthits, minPopulationHits=5)
{
   if ( (poplisthits[2]<minPopulationHits) | (poplisthits[4]==0) )
      return (1)

   q= poplisthits[4] #list hit
   m= poplisthits[2] #population hit
   n= poplisthits[1]-poplisthits[2] #population non-hit
   k= poplisthits[3] #list total

   myp=phyper(q-1, m, n, k, lower.tail=F)
   signif(myp,3)
}


coexEnrich <- function(gn,d,bg=NULL,fc.thresh=0,pvalCut=NULL) {
  
  if (is.null(bg)) { bg <- unique(unlist(gn)) }
  if (is.null(pvalCut))  { pvalCut=0.05 }
  
  N<-length(unique(bg))
  X<-length(intersect(unique(bg),d))
  out <- NULL
  for (m in 1:length(gn) ) {
    n<-length(unique(gn[[m]]))
    x<-length(intersect(gn[[m]],d))
    message("Processing ", names(gn)[m], ": ", n, "," , x)
    fc <- (x/n)/(X/N)
    genes <- paste(intersect(gn[[m]],d), collapse=" ") # genes from 'x'
    hgd<- fisherTest(c(N,X,n,x))
    if (!is.nan(fc)) {
      if (fc >= fc.thresh) {
        o <- c(names(gn)[m], x, n, X, N, fc, hgd, genes)
      }
    }
    else {
      o <- c(m, NA, NA, NA, NA, NA, NA, NA)
    }
    out <- rbind(out,o)
  }
  p.val <- as.numeric(out[,7]) # p.value  
  all.pws<-length(which(as.numeric(out[,2])!=0)) # x
  corr.p.val<-p.adjust(p.val, method="fdr")
  out<-cbind(out[,-8],corr.p.val,out[,8])
  colnames(out)<-c("Module", "x", "n", "X", "N", "fold.change", "p.value", "corr.p.val", "Genes")
return(out)
}

