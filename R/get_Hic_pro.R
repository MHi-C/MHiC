get_Hic_pro<-function(reads_file,Digest_file = NULL,method){
  
  if(method=="gothic" || method=="GOTHiC"){
  filetype<-grepl("\\.matrix$", reads_file)
  if(filetype == FALSE){filetype<-grepl("\\.table$", reads_file)}
  if(filetype == FALSE){
  filetype<-grepl("\\.txt$", reads_file)}
  ##########################################################
  #import digest file
  
  ranges<-read.table(Digest_file, header = FALSE)
  
  ranges$chr<-ranges$V1
  ranges$start<-ranges$V2
  ranges$end<-ranges$V3
  ranges$id<-ranges$V4
  ranges<-ranges[,c("chr","start","end","id")]
  levels(ranges$chr) <- fixChromosomeNames(levels(ranges$chr))
  ##########################################################
  #import  reads file
  if(filetype){
    reads<-read.table(reads_file, header = TRUE) 
  }
  else{stop('Hic_pro output must be in a matrix format')}
  
  colnames(reads)<-c("locus1","locus2","interaction_counts")
  ##########################################################
  # map reads 
  fragments1<- reads$locus1
  fragments2<- reads$locus2
  
  loci1 <- ranges[match(fragments1,ranges$id),]
  loci2 <- ranges[match(fragments2,ranges$id),]
  rm(fragments1,fragments2)
  ##########################################################
  interactions<-data.frame(loci1$chr,loci1$start,loci2$chr,loci2$start
                           ,reads$interaction_counts)
  
  names(interactions)<-c("chr1", "locus1_start", "chr2", "locus2_start"
                         ,"frequencies")
  
  return(interactions)
  }
  
  
  
  
  
  ####################################################################################################################
  
  
  if( method=="fithic" || method=="FitHiC"){
    
   
    filetype<-grepl("\\.matrix$", reads_file)
    if(filetype == FALSE){filetype<-grepl("\\.table$", reads_file)}
    if(filetype == FALSE){filetype<-grepl("\\.txt$", reads_file)}
    if(filetype == FALSE){filetype<-grepl("\\.matrix.gz$", reads_file)}
    ##########################################################
    #import digest file
    
    ranges<-read.table(Digest_file, header = FALSE)
   
    
    names(ranges)<-c("chr","start","end","id")

    levels(ranges$chr) <- fixChromosomeNames(levels(ranges$chr))
    
    fragment <- data.frame(chr=ranges$chr, mid=((ranges$start + ranges$end) / 2),
                           index=ranges$id)
    
    names(fragment)<-c("chr","mid","index")
    ranges<-fragment;
  
    ##########################################################
    #import  reads file
    if(filetype){
      reads<-read.table(reads_file, header = FALSE) 
     
    }
    else{stop('Hic_pro output must be in a matrix format')}
    
    
    colnames(reads)<-c("locus1","locus2","interaction_counts")
    ##########################################################
    # map reads 
    fragments1<- reads$locus1
    fragments2<- reads$locus2
    
    loci1 <- ranges[match(fragments1,ranges$index),]
    loci2 <- ranges[match(fragments2,ranges$index),]
  
    rm(fragments1,fragments2)
    ##########################################################
    interactions<-data.frame(loci1$chr,loci1$mid,loci2$chr,loci2$mid
                             ,reads$interaction_counts)
    
    names(interactions)<-c("chr1", "mid1", "chr2", "mid2"
                           ,"hitCount")
    
    FitHiCstruct<-list(interactions=interactions,digest=ranges)
    
    return(FitHiCstruct)
    
    
  }
}