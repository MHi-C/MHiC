get_Hic_data<-function(reads_file,Digest_file= NULL,tools_name,method){
  
  if(method=="gothic" || method=="GOTHiC" || method=="fithic" || method=="FitHiC"){
  if(is.null(tools_name)){
    
    interactions<-get_MHic_data(reads_file,method)
  }
  ##########################################################
  
  if(tools_name == 'HiC_PRO'|tools_name == 'HIC_PRO' |tools_name == 'hic_pro' |tools_name == 'Hic_Pro'){
    
  
      interactions<-get_Hic_pro(reads_file,Digest_file,method)
      
    
  }
  ##########################################################
  
  else if(tools_name == 'HiCUP' | tools_name == 'HICUP'|tools_name == 'hicup'){
    
      interactions<-get_HiCUP(reads_file,Digest_file,method)
  }
  ##########################################################
  
  else if(tools_name == 'HOMER' | tools_name == 'homer'){
    
    interactions<-get_HOMER(reads_file)
  }
  ##########################################################
  
    return(interactions)
  }
  else if(method=="hicnorm"){
    
    
    
    interactions<-get_norm_data(reads_file,Digest_file)
      
    
    return(interactions)
  }
 
  
}
