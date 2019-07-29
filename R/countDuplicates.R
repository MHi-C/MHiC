countDuplicates<-function(interactions,column_names = c("frequencies","frequency"
                                                        , "read_counts","count")){
  
  include<-any(column_names %in% names(interactions))
  ##########################################################
   if(include){
     column_name<- column_names[column_names %in% names(interactions)][1]
     if(is.null(interactions[,column_name]) ){
       stop("frequencies column missing")
     }
   }
  ##########################################################
  else{
    column_name <- column_names[1]
    interactions[, column_name] <- 1
  }
  ##########################################################
  interactions <- interactions[order(interactions$chr1, interactions$locus1_start
                                     , interactions$chr2, interactions$locus2_start), ]
  ##########################################################
  
  duplicates <- c(duplicated(interactions[, -which(names(interactions) == column_name)]))
  
  
  cumulative <- c(0, cumsum(interactions[, column_name]))
  beforeFirstUnique <- cumulative[!c(duplicates, FALSE)] 
  
  result <- interactions[!duplicates, ]
  result[, column_name] = diff(beforeFirstUnique)
  
  return(result)
  
}