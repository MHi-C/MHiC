calculate_coverage<-function(interactions,flag = TRUE){
  
  #
  ##########################################################
  if(nrow(interactions)>1e8){
    
    t <- ceiling(nrow(interactions)/1e8)
    IList <- list()
    IList[[1]] <- interactions[1:1e8,]
    
    for(i in 2:t){
      IList[[i]] <- interactions[(((i-1)*1e8)+1):min((i*1e8),nrow(interactions)),]
    }
    dtList <- lapply(IList, data.table)
    
    covAs <- lapply(dtList, function(x) x[,sum(frequencies), by=int1])
    covBs <- lapply(dtList, function(x) x[,sum(frequencies), by=int2])
    
    covAm <- do.call(rbind, covAs)
    covBm <- do.call(rbind, covBs)
    covA <- covAm[,sum(V1),by=int1]
    covB <- covBm[,sum(V1),by=int2]	
    
  }
  ##########################################################
  else{
    
    binned_interactions=data.table(interactions)
    covA <- binned_interactions[,sum(frequencies),by=int1]	
    covB <- binned_interactions[,sum(frequencies),by=int2]
    
  }
  ##########################################################
  covA <- setkey(covA,key='int1')
  setnames(covB, 1,'int1')
  covB <- setkey(covB,key='int1')
  
  cov=merge(covA,covB,all.x=TRUE,all.y=TRUE,by='int1')
  cov$V1.x[is.na(cov$V1.x)]=0
  cov$V1.y[is.na(cov$V1.y)]=0
  cov$coverage=cov$V1.x+cov$V1.y
  coverage=cov$coverage
  names(coverage)=cov$int1
  sumcov <- sum(coverage)
  relative_coverage <- coverage/sumcov
  names(relative_coverage)=names(coverage)
  interactions$coverage_source <- relative_coverage[interactions$int1]
  interactions$coverage_target <- relative_coverage[interactions$int2]
  ##########################################################
  
  if(flag){return(interactions)}
  
  else{return(relative_coverage)}
  
  
}