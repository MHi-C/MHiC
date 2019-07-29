addhicnorm<-function(interactions,res){
  
  
  interactions<-interactions[which(!is.infinite(interactions$expected)),]
  interactions$locus1<-(interactions$locus1 %/% res)*res
  interactions$locus2<-(interactions$locus2 %/% res)*res
  final<-ddply(interactions, .(chr1,locus1,chr2,locus2), numcolwise(sum))
  return(final)
  
}