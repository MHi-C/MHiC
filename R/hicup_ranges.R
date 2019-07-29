
hicup_ranges<-function(Digest_file){
  
  #
ranges<-read.table(Digest_file,stringsAsFactors=FALSE ,header = TRUE, skip = 1)

levels(ranges$Chromosome) <- fixChromosomeNames(levels(ranges$Chromosome))

##########################################################
ranges$chr<-ranges$Chromosome
ranges$start<-ranges$Fragment_Start_Position
ranges$end<-ranges$Fragment_End_Position
ranges$id<-rownames(ranges)
ranges<-ranges[,c("chr","start","end","id")]
##########################################################

return(ranges)

}