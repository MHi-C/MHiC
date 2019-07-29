interaction_res<-function(interactions,res){
  
  if(!identical(colnames(interactions)[1:5], c("chr1", "locus1_start"
                                               , "chr2", "locus2_start"
                                               ,"frequencies"))){
    
    stop("expecting columns chr1,locus1_start,chr2,locus2_start
         ,frequencies")
  }
  ##########################################################
  #map data to bin size
  interactions$locus1_start <- (interactions$locus1_start %/% res) * res
  
  interactions$locus2_start <- (interactions$locus2_start %/% res) * res
  
  ##########################################################
  #put smaller locus
  firstSmaller <- subset(interactions, (as.character(chr1) < as.character(chr2)) 
          | (as.character(chr1) == (as.character(chr2)) & locus1_start <= locus2_start))
  secondSmaller <- subset(interactions, !((as.character(chr1) < as.character(chr2)) 
          | (as.character(chr1) == (as.character(chr2)) & locus1_start <= locus2_start)))
  
  cnames <- colnames(interactions)
  cnames[1:5] <- c("chr2", "locus2_start","chr1"
                   , "locus1_start", "frequencies")
  setnames(secondSmaller,colnames(secondSmaller),cnames)
  interactions <- rbind(firstSmaller, secondSmaller)
  interactions <- interactions[order(interactions$chr1, interactions$locus1_start
                                     , interactions$chr2, interactions$locus2_start), ]
  ##########################################################
  #merge interactions which in a same loci
  interactions <-countDuplicates(interactions)	
  ##########################################################
  return(interactions)
}