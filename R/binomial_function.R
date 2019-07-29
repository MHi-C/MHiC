binomial_function<-function(reads_file,tools_name
                            ,cistrans,parallel,cores
                             ,removeDiagonal){
 
  ##########################################################
  
  #
  if(parallel)
  {
    requireNamespace("parallel")
    print("running garbage collector before parallel fork")
    gc()
  }
  ##########################################################
  interactions <-reads_file
  interactions$int1 <-paste(interactions$chr1,interactions$locus1_start
                            ,sep='_')
  interactions$int2 <-paste(interactions$chr2,interactions$locus2_start
                            ,sep='_')
  ##########################################################
  #diagonal removal
  if(removeDiagonal)
  {
    interactions <- interactions[interactions$int1!=interactions$int2,]	
  }
  #
  if(cistrans=='cis'){
    interactions <- interactions[interactions$chr1==interactions$chr2,]	
  }
  if(cistrans=='trans'){
    interactions <- interactions[interactions$chr1!=interactions$chr2,]	
  }
  ##########################################################
  #all read pairs used in binomial 
  
  numberOfReadPairs  <- sum(interactions$frequencies)
  
  all_bins <- unique(c(unique(interactions$int1), unique(interactions$int2)))
  all_bins <- sort(all_bins)
  ##########################################################
  #calculate coverage
  
  interactions<-calculate_coverage(interactions)
  relative_coverage<-calculate_coverage(interactions, FALSE)
  ##########################################################
  #probability correction assuming on average equal probabilities for all interactions
  numberOfAllInteractions <- length(all_bins)^2
  upperhalfBinNumber <- (length(all_bins)^2-length(all_bins))/2
  ##########################################################
  #
  if(cistrans!='all'){
    chromos <- unique(interactions$chr1)
    chrlens <- c()
    for(cr in chromos){ 
      chrlens[cr] <- length(unique(c(unique(interactions$locus1_start[interactions$chr1==cr])
                                     ,unique(interactions$locus2_start[interactions$chr2==cr]))))
    }
    cisBinNumber <-(sum(chrlens^2)-length(all_bins))/2	
    transBinNumber <- upperhalfBinNumber-cisBinNumber
  }
  ##########################################################
  #
  diagonalProb <- sum(relative_coverage^2)
  if(cistrans=='all'){
    probabilityCorrection <- if(removeDiagonal){1/(1-diagonalProb)}else{1}
  }
  if(cistrans=='cis'){
    probabilityCorrection <- upperhalfBinNumber/cisBinNumber
  }
  if(cistrans=='trans'){
    probabilityCorrection <- upperhalfBinNumber/transBinNumber
  }
  ##########################################################
  max1=max(interactions$locus1_start)
  max2=max(interactions$locus2_start)
  max= max(max1,max2)
  min1=min(interactions$locus1_start)
  min2=min(interactions$locus2_start)
  min= min(min1,min2)
  
  interactions$probabilityOfInteraction <- interactions$coverage_source*interactions$coverage_target*2*probabilityCorrection
  z=(max-abs(interactions$locus1_start-interactions$locus2))/(nrow(interactions))
  interactions$probabilityOfInteraction <-interactions$probabilityOfInteraction+((1-sum(interactions$probabilityOfInteraction))/(max - min))*z
  
  interactions$probabilityOfInteraction <-interactions$probabilityOfInteraction+((1-sum(interactions$probabilityOfInteraction))/nrow(interactions))
  interactions$predicted <- interactions$probabilityOfInteraction * numberOfReadPairs
  ##########################################################
  #
  interactions<-main_process(interactions,numberOfReadPairs,parallel,cores)
  ##########################################################
  #
  #observed over expected log ratio
  
  
  interactions$logFoldChange <- log2(interactions$frequencies/interactions$predicted)
  #multiple testing correction separately for matrices with all interactions/only cis/only transs
  
  
  if(cistrans=='all'){
    interactions$qvalue <- if(removeDiagonal){p.adjust(interactions$pvalue
                                               , method = "BH", n=upperhalfBinNumber)}else{p.adjust(interactions$pvalue
                                                                                                    , method = "BH", n=upperhalfBinNumber+length(all_bins))}
  }
  else if(cistrans=='cis'){
    interactions$qvalue <- if(removeDiagonal){p.adjust(interactions$pvalue
                                               , method = "BH", n=cisBinNumber)}else{p.adjust(interactions$pvalue
                                                                                              , method = "BH", n=cisBinNumber+length(all_bins))}	
  }
  else if(cistrans=='trans'){
    interactions$qvalue <- p.adjust(interactions$pvalue, method = "BH"
                                    , n=transBinNumber)	
  }
  ##########################################################
  #
  interactions=interactions[,c('chr1',"locus1_start",'chr2',"locus2_start"
                               ,'coverage_source','coverage_target','probabilityOfInteraction'
                               , 'predicted','frequencies', 'pvalue','qvalue','logFoldChange')]
  
  colnames(interactions)=c('chr1',"locus1",'chr2',"locus2"
                           ,'relCoverage1','relCoverage2','probability', 'expected'
                           ,'readCount', 'pvalue','qvalue','logObservedOverExpected')
  
  
  return(interactions)
  
  
}