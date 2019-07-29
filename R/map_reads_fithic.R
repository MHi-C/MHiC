map_reads_fithic<-function(reads_file,Digest_file){
  ##########################################################
  #
  ranges<-GRanges(seqnames=Digest_file$chr
                  , ranges=IRanges(start=Digest_file$start, end=Digest_file$end))
  
  ##########################################################
  fragment1<-reads_file[,c("chr1","locus1")]
  fragment1<-overlap(fragment1,Digest_file)
  
  fragment2<-reads_file[,c("chr2","locus2")]
  fragment2<-overlap(fragment2,Digest_file)
  ##########################################################
  interaction<-!is.null(fragment1) & !is.null(fragment2)
  
  
  fragment1 <- fragment1[interaction]
  fragment2 <- fragment2[interaction]
  fragment1_1 <- pmin(fragment1, fragment2)
  fragment2_2 <- pmax(fragment1, fragment2)
  
  ##########################################################
  interactions <- as.data.frame(cbind(as.character(seqnames(ranges[fragment1_1]))
                                      , abs(start(ranges(ranges[fragment1_1]))+
                                        end(ranges(ranges[fragment1_1])))/2
                                      , as.character(seqnames(ranges[fragment2_2]))
                                      , abs(start(ranges(ranges[fragment2_2]))+
                                        end(ranges(ranges[fragment2_2])))/2
                                      , rep(1, times=length(fragment1_1)))
                                , stringsAsFactors=FALSE)
  
  names(interactions)<-c("chr1","mid1","chr2","mid2"
                         ,"hitCount")
  ##########################################################
  interactions$mid1=as.numeric(interactions$mid1)
  interactions$mid2=as.numeric(interactions$mid2)
  interactions$hitCount=as.numeric(interactions$hitCount)
  ##########################################################
  return(interactions)
}