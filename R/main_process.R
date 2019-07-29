main_process<-function(interactions,numberOfReadPairs,parallel,cores){
  #
  #
  ##########################################################
  if(parallel){
    
    requireNamespace("parallel")
    if(nrow(interactions)>1e8){
      
      t <- ceiling(nrow(interactions)/1e8)
      dfList <- list()
      for(i in 2:t){
        dfList[[i]] <- interactions[(((i-1)*1e8)+1):min((i*1e8),nrow(interactions)),]
        
      }
      dtList <- lapply(dfList, function(x) as.data.frame(t(cbind(as.numeric(x[["frequencies"]]), 
                                                                 as.numeric(x[["probabilityOfInteraction"]])))))
      
      pvalues=list()
      for(i in 1:length(dtList)){
        pvalues[[i]] <-unlist(parallel::mclapply(dtList[[i]], function(x)
        {
          binom.test(x[1]-1, numberOfReadPairs, x[2], alternative = "greater")$p.value
        }, 
        mc.cores=cores))
        
      }
      
      pvals=unlist(pvalues)
      interactions$pvalue <- pvals
      
    }else{
      
      binomParams <- as.data.frame(t(cbind(
        as.numeric(interactions[["frequencies"]]), 
        as.numeric(interactions[["probabilityOfInteraction"]]))))
      
      interactions$pvalue <- unlist(parallel::mclapply(binomParams, function(x){
        binom.test(x[1]-1, numberOfReadPairs, x[2], alternative = "greater")$p.value
      }, mc.cores=cores))
    }
    
  }
  ##########################################################
  #
  else{
    
    if(nrow(interactions)>1e8){
      
      t <- ceiling(nrow(interactions)/1e8)
      dfList <- list()
      dfList[[1]] <- interactions[1:1e8,]
      for(i in 2:t){
        dfList[[i]] <- interactions[(((i-1)*1e8)+1):min((i*1e8),nrow(interactions)),]
      }
      pvalues=list()
      
      for(i in 1:length(dfList)){
        pvalues[[i]] <-apply(dfList[[i]], 1, function(x)
        {
          binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
        }	
        )
      } 
      pvals=unlist(pvalues)
      interactions$pvalue <- pvals
      
    }else{
      
      interactions$pvalue <- apply(interactions, 1, function(x)
      {
        binom.test(as.numeric(x[["frequencies"]])-1, numberOfReadPairs, as.numeric(x[["probabilityOfInteraction"]]), alternative = "greater")$p.value
      }	
      )
      
    }
    
  }
  ##########################################################
  return(interactions)
}