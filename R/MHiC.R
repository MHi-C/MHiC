MHiC<-function(reads_file,Digest_file = NULL , outdir ,method="gothic",sample_name
               ,tools_name = NULL,res = 1000000 ,cistrans = "all"
               , parallel=FALSE, cores=NULL, removeDiagonal=TRUE,save = TRUE 
               , remove = FALSE
               ,min_cov = 2
               ,min_len = 0.1
               ,min_gc = 0.3
               ,min_map = 0.8
               , biasfile="none", noOfPasses=1,
               noOfBins=100, mappabilityThreshold=1, distUpThres=-1,
               distLowThres=-1
               ){
  
  if(tools_name == 'HOMER' | tools_name == 'homer'){
  
    stop("This part is Under the update")
    
  }
 
  
  message("MHiC processing started")
  #this tools get Hi-C data from different tools
  #function for get different Hi-C structure
  interactions<-get_Hic_data(reads_file,Digest_file,tools_name,method)
  
  return(interactions)
 
    if( method=="gothic" || method=="GOTHiC"){
    
  interactions<-interaction_res(interactions,res)
  
  }
  

  
  
  if(method=="gothic" || method=="GOTHiC"){
  #apply binomial function for identify significant interaction
  Result<-binomial_function(interactions,tools_name
                            ,cistrans,parallel
                            ,cores,removeDiagonal)
  Result$valid<-0
  Result[-log10(Result$qvalue) >= 1.3,]$valid<-1
  #remove non significant interaction
  finalResult<-Result[-log10(Result$qvalue) >= 1.3,]
  }
  else if(method=="hicnorm"){
    
    
      Result<-HiCnorm(interactions,tools_name
                                ,cistrans
                                ,Digest_file,parallel
                                ,cores,removeDiagonal
                                ,min_cov 
                                ,min_len 
                                ,min_gc 
                                ,min_map 
                                )
      Result$len<-NULL
      Result$gc<-NULL
      Result$map<-NULL
      fm<-abs(Result$locus2-Result$locus1);
      fm<-fm[which(fm!=0)]
     if( min(fm)!=res ){
      
       Result<-addhicnorm(Result,res)
       
     }
      print(" HiCNorm cannot remove invalid interactions")
      Result<-Result[which(!is.infinite(Result$expected)),]
      finalResult<-Result[which(!is.infinite(Result$expected)),]
      
      
  }
  else if(method=="fithic" || method=="FitHiC"){
    
   
  FitHiCm(intersData= interactions$interactions ,fragsData=interactions$digest ,
                   biasfile=biasfile, noOfPasses=noOfPasses,
                   noOfBins=noOfBins, mappabilityThreshold=mappabilityThreshold, 
                   libname=sample_name, distUpThres=distUpThres,
                   distLowThres=distLowThres)
    
    return("Finish")
    
  }
  #save result
  if(remove==TRUE){
    print("Operation Successfully Completed")
  if(save == TRUE){
    
  write.csv(finalResult,file=paste(sample_name, ".csv", sep=""),row.names = FALSE)
    
  }
    return(finalResult)
  }
  else{
    print("Operation Successfully Completed")
    if(save == TRUE){
      
      write.csv(Result,file=paste(sample_name, ".csv", sep=""),row.names = FALSE)
      
    }
    return(Result)
  }
  
  
}
