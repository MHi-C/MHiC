get_HOMER<-function(interaction_file){
  #get Homer interaction matrix 
  #interaction matrix is output from analyzeHiC part in homer
  #import  file
  filetype<-grepl("\\.txt$", interaction_file)
  

    if(!filetype){
      reads<-read.matrix(interaction_file, header = FALSE, sep = "\t") 
    }
  #homer tools can make different interaction matrix
  #this tools accept Homer outputs with raw interactions counts between the regions  
  reads<-reads[-1,-1]
  
  names<-reads[,1]
  reads<-reads[,-1]
  reads <- Matrix(reads, sparse = TRUE)
  colnames(reads)<-names
  rownames(reads)<-names
  
  interactions <- summary(reads)
  interactions<-data.frame(Origin      = as.character(rownames(reads)[interactions$i]),
                           Destination = as.character(colnames(reads)[interactions$j]),
                           frequencies = as.numeric(interactions$x))
  interactions$chr1<-as.character(unlist(strsplit(interactions$Origin, "-"))[1])
  interactions$locus1_start<-as.numeric(unlist(strsplit(interactions$Origin, "-"))[2])
  interactions$chr2<-as.character(unlist(strsplit(interactions$Destination, "-"))[1])
  interactions$locus2_start<-as.numeric(unlist(strsplit(interactions$Destination, "-"))[2])
  
  
  interactions<-interactions[,c("chr1", "locus1_start", "chr2", "locus2_start"
                                ,"frequencies")]
  
  names(interactions)<-c("chr1", "locus1_start", "chr2", "locus2_start"
                         ,"frequencies")
  return(interactions)
}