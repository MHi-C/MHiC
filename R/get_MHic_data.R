get_MHic_data<-function(reads_file,method){
  
  filetype<-grepl("\\.matrix$", reads_file)
  if(filetype == FALSE){filetype<-grepl("\\.table$", reads_file)}
  if(filetype== FALSE){
    filetype<-grepl("\\.txt$", reads_file)}
  
  ##########################################################
  if(filetype){
    reads<-read.table(reads_file, header = FALSE, sep = "\t") 
  }
  
  else{stop('data must be in MHiC format')}
  
  if(method=="gothic"){
  interactions<-data.frame(reads$V1,reads$V2,reads$V3,reads$V4
                           ,reads$V5)
  
  names(interactions)<-c("chr1", "locus1_start", "chr2", "locus2_start"
                         ,"frequencies")
  }
   if(method=="hicnorm"){
     
    interactions<-data.frame(reads$V1,reads$V2,reads$V3,reads$V4
                             ,reads$V5,reads$V6,reads$V7,reads$V8)
    
    names(interactions)<-c("chr1", "locus1_start", "chr2", "locus2_start"
                           ,"frequencies","len","gc","map")

  }
  return(interactions)
  
}