overlap<-function(reads_file,Digest_file){
  
ranges<-GRanges(seqnames=Digest_file$chr
                , ranges=IRanges(start=Digest_file$start, end=Digest_file$end))
reads<-GRanges(seqnames =	Rle(reads_file[,1])
               , ranges = IRanges(reads_file[,2], width = 1),strand = Rle("*", nrow(reads_file)))
fragment<- findOverlaps(reads, ranges, select="first")

  return(fragment)
}