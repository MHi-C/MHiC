get_HiCUP<-function(reads_file,Digest_file,method){
  if(method=="gothic" || method=="GOTHiC"){
    
  compressed<-grepl("\\.gz$", reads_file)
  filetype<-grepl("\\.table$", reads_file)
  if(filetype == FALSE){
  filetype<-grepl("\\.txt$", reads_file)}
  ##########################################################
  
  #import digest file
  ranges<-hicup_ranges(Digest_file)
  
  ##########################################################
  
  #import  reads file
  if(compressed){
    
    reads<-read.table(gzfile(reads_file)) 
    
  } 
  else{
  if(filetype){
    reads<-read.table(reads_file) 
  }
 
  else{stop('HiCUP output must be in a table format')}
  
  }
  names(reads) <- c("id", "flag", "chr", "locus")
  reads<-reads[order(reads$id, reads$flag, reads$chr, reads$locus),]
  ##########################################################
  
  
  #in the sam file the two ends of the interactions are consecutive rows
  odd <- reads[seq(1,nrow(reads),2), ]
  even <- reads[seq(2,nrow(reads),2), ]
  names(odd) <- c("id1", "flag1", "chr1", "locus1")
  names(even) <- c("id2", "flag2", "chr2", "locus2")	
  
  levels(odd$chr1)<-fixChromosomeNames(levels(odd$chr1))
  levels(even$chr2)<-fixChromosomeNames(levels(even$chr2))
  
  
  if(nrow(odd) != nrow(even)){stop("odd and even rows of data must have same counts")}
  joined <- cbind(odd, even)
  if(!all(joined$id1==joined$id2)){stop(" each odd and even row of data must have same id ")}
  ##########################################################
  
  
  interactions<-joined[,c("chr1", "locus1","chr2", "locus2")]
  interactingLoci<-map_reads_HiCUP(interactions,ranges)
  ##########################################################
  
  return(interactingLoci)
  }
  
  ####################################################################################################################
  
  
  if( method=="fithic" || method=="FitHiC"){
    
    compressed<-grepl("\\.gz$", reads_file)
    filetype<-grepl("\\.table$", reads_file)
    if(filetype == FALSE){
      filetype<-grepl("\\.txt$", reads_file)}
    ##########################################################
    
    #import digest file
    ranges<-hicup_ranges(Digest_file)
    
    
    fragment <- data.table(chr=ranges$chr, mid=(ranges$start + ranges$end) / 2,
                           index=ranges$id)
    names(fragment)<-c("chr","mid","index")
    
    ##########################################################
    
    #import  reads file
    if(compressed){
      
      reads<-read.table(gzfile(reads_file)) 
      
    } 
    else{
      if(filetype){
        reads<-read.table(reads_file) 
      }
      
      else{stop('HiCUP output must be in a table format')}
      
    }
    names(reads) <- c("id", "flag", "chr", "locus")
    reads<-reads[order(reads$id, reads$flag, reads$chr, reads$locus),]
    ##########################################################
    
    
    #in the sam file the two ends of the interactions are consecutive rows
    odd <- reads[seq(1,nrow(reads),2), ]
    even <- reads[seq(2,nrow(reads),2), ]
    names(odd) <- c("id1", "flag1", "chr1", "locus1")
    names(even) <- c("id2", "flag2", "chr2", "locus2")	
    
    levels(odd$chr1)<-fixChromosomeNames(levels(odd$chr1))
    levels(even$chr2)<-fixChromosomeNames(levels(even$chr2))
    
    
    if(nrow(odd) != nrow(even)){stop("odd and even rows of data must have same counts")}
    joined <- cbind(odd, even)
    if(!all(joined$id1==joined$id2)){stop(" each odd and even row of data must have same id ")}
    ##########################################################
    
    
    interactions<-joined[,c("chr1", "locus1","chr2", "locus2")]
    interactingLoci<-map_reads_fithic(interactions,ranges)
    ##########################################################
    ranges<-fragment;
  
    FitHiCstruct<-list(interactions=interactingLoci,digest=ranges)
    
    return(FitHiCstruct)
  }
}