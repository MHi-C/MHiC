HiCnorm<-function(reads_file,tools_name
                            ,cistrans= "cis",digestfile ,parallel,cores
                            ,removeDiagonal,min_cov,min_len,min_gc,min_map){
  
  
 
  
  interactions <-reads_file

  interactions$int1 <-paste(interactions$chr1,interactions$locus1_start
                            ,sep='_')
  interactions$int2 <-paste(interactions$chr2,interactions$locus2_start
                            ,sep='_')
  
  ##########################################################
  #diagonal removal
  if(removeDiagonal)
  {
    interactions <- interactions[interactions$int1!=interactions$int2,]	
  }
  #
 
  if(cistrans=='trans'){
    interactions <- interactions[interactions$chr1!=interactions$chr2,]	
    stop("HiCNorm only work on cis interacion")
  }
  else(

      interactions <- interactions[interactions$chr1==interactions$chr2,]	
    
  )
  
  interactions$int1<-NULL
  interactions$int2<-NULL
  #if(gc== FALSE && len== FALSE && map== FALSE){
  #  reads_file$Dis<-abs(as.numeric(reads_file$locus1_start1)-as.numeric(reads_file$locus2_start))
  #  model<-glm(frequencies~Dis, data = reads_file, family = "poisson")
    
    
    
  #}
  #if(gc== FALSE && len== FALSE && map== TRUE){
  #  reads_file$Dis<-abs(as.numeric(reads_file$locus1_start1)-as.numeric(reads_file$locus2_start))
  #  model<-glm(frequencies~Dis, data = reads_file, family = "poisson")
    
    
    
  #}
  
  
  
        
        chr_n <- table(interactions$chr1)
        chr_id <- names(chr_n)
        final<-NULL
        
        
        if(!is.null(digestfile)){
          
          Digest<-read.table(digestfile, header = FALSE, sep = "\t")
          names(Digest) <- c('chr', 'bin1', 'bin2', 'len', 'gc', 'map')
        for (chr in chr_id) {
          tryCatch(
            {
              genomic_features_chr<- Digest[which(Digest$chr==chr),]
              map_chr <- interactions[which(interactions$chr1==chr),]
             
              conf_bins <- which(genomic_features_chr$len/(genomic_features_chr$bin2-genomic_features_chr$bin1)>=min_len & genomic_features_chr$gc>=min_gc & genomic_features_chr$map>=min_map) 
              
              ind1 <- match(map_chr$locus1_start,genomic_features_chr$bin1)
              ind2 <- match(map_chr$locus2_start,genomic_features_chr$bin1)
              
              map_chr$len <- log(genomic_features_chr$len[ind1] * genomic_features_chr$len[ind2])
              map_chr$gc <- log(genomic_features_chr$gc[ind1] * genomic_features_chr$gc[ind2])
              map_chr$map <- log(genomic_features_chr$map[ind1] * genomic_features_chr$map[ind2])
              
              res <- map_chr
              
              map_fit <- map_chr[!is.na(map_chr$len) & map_chr$frequencies >= min_cov & ind1 %in% conf_bins & ind2 %in% conf_bins,]
              map_chr$len <- (map_chr$len - mean(map_fit$len)) / sd(map_fit$len)
              map_chr$gc <- (map_chr$gc - mean(map_fit$gc)) / sd(map_fit$gc)
              
              map_fit$len <- (map_fit$len - mean(map_fit$len)) / sd(map_fit$len)
              map_fit$gc <- (map_fit$gc - mean(map_fit$gc)) / sd(map_fit$gc)
              
             
                fit <- glm(frequencies~len+gc+offset(map), data = map_fit, family = "poisson")
              
              coeff <- round(fit$coeff,4)
              res$nor <- round(map_chr$frequencies / exp(coeff[1] + coeff[2] * map_chr$len + coeff[3] * map_chr$gc + map_chr$map), 4)
              
              final <- rbind(final,res)
            })
        } 
        }
        else{
          
          
          Digest<-data.frame(interactions$chr1
                             ,interactions$locus1_start
                             ,interactions$locus2_start
                             ,interactions$len
                             ,interactions$gc
                             ,interactions$map)
          interactions$gc<-NULL
          interactions$len<-NULL
          interactions$map<-NULL
          names(Digest) <- c('chr', 'bin1', 'bin2', 'len', 'gc', 'map')
          for (chr in chr_id) {
            tryCatch(
              {
                genomic_features_chr<- Digest[which(Digest$chr==chr),]
                map_chr <- interactions[which(interactions$chr1==chr),]
                
                conf_bins <- which(genomic_features_chr$len/(genomic_features_chr$bin2-genomic_features_chr$bin1)>=min_len & genomic_features_chr$gc>=min_gc & genomic_features_chr$map>=min_map) 
                
                ind1 <- match(map_chr$locus1_start,genomic_features_chr$bin1)
                ind2 <- match(map_chr$locus2_start,genomic_features_chr$bin1)
                
                map_chr$len <- log(genomic_features_chr$len[ind1] * genomic_features_chr$len[ind2])
                map_chr$gc <- log(genomic_features_chr$gc[ind1] * genomic_features_chr$gc[ind2])
                map_chr$map <- log(genomic_features_chr$map[ind1] * genomic_features_chr$map[ind2])
                
                res <- map_chr
                
                map_fit <- map_chr[!is.na(map_chr$len) & map_chr$frequencies >= min_cov & ind1 %in% conf_bins & ind2 %in% conf_bins,]
                map_chr$len <- (map_chr$len - mean(map_fit$len)) / sd(map_fit$len)
                map_chr$gc <- (map_chr$gc - mean(map_fit$gc)) / sd(map_fit$gc)
                
                map_fit$len <- (map_fit$len - mean(map_fit$len)) / sd(map_fit$len)
                map_fit$gc <- (map_fit$gc - mean(map_fit$gc)) / sd(map_fit$gc)
                
               
                fit <- glm(frequencies~len+gc+offset(map), data = map_fit, family = "poisson")
                
                coeff <- round(fit$coeff,4)
               
                res$nor <- round(map_chr$frequencies / exp(coeff[1] + coeff[2] * map_chr$len + coeff[3] * map_chr$gc + map_chr$map), 4)
                
                final <- rbind(final,res)
              })
          } 
          
        }
        
        colnames(final)=c('chr1',"locus1",'chr2',"locus2"
                                 ,'readCount', 'len','gc','map','expected')
        
        
        return(final)
       
  
  
}