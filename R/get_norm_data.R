get_norm_data<-function(reads_file,Digest_file){
 
  filetype<-grepl("\\.matrix$", reads_file)
  if(filetype == FALSE){filetype<-grepl("\\.table$", reads_file)}
  if(filetype== FALSE){
    filetype<-grepl("\\.txt$", reads_file)}
  
#################################################################### 
  filetype1<-grepl("\\.matrix$", Digest_file)
  if(filetype1 == FALSE){filetype1<-grepl("\\.table$", Digest_file)}
  if(filetype1== FALSE){
    filetype1<-grepl("\\.txt$", Digest_file)}
  if(filetype1== FALSE){
    filetype1<-grepl("\\.bed$", Digest_file)}

####################################################################
  if(filetype){
    reads<-read.table(reads_file, header = FALSE) 
   
  }
  
  #if(filetype1){
  #  Digest<-read.table(Digest_file, header = FALSE) 
    
  #}
 
  names(reads) <- c("chr",'bin1', 'bin2', 'raw')
  #names(Digest) <- c('chr', 'start', 'end',"gc","len","map")
  ####################################################################
  #print(reads)
  
  #fragments1<- reads$binid1
  #fragments2<- reads$binid2
  
  #loci1 <- Digest[match(fragments1,Digest$start),]
  #loci2 <- Digest[match(fragments2,Digest$end),]
  #rm(fragments1,fragments2)
  ##########################################################
  interactions<-data.frame(reads$chr,reads$bin1,reads$chr,reads$bin2
                           ,reads$raw)
  
  
  
  
  
  names(interactions)<-c("chr1", "locus1_start", "chr2", "locus2_start"
                        ,"frequencies")

  
  return(interactions)
}