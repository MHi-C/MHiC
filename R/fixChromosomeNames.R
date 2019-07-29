fixChromosomeNames <- function(chrnames)
{
  #capital to small
  chrnames <- paste0("chr", chrnames)
  chrnames <- sub("CHR", "chr", chrnames)
  chrnames <- sub("chrchr", "chr", chrnames)
  
 
}