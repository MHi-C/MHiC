FitHiCm<-function(intersData ,fragsData, biasfile="none", noOfPasses=1,
                          noOfBins=100, mappabilityThreshold=1, libname,
                          distUpThres=-1,
                          distLowThres=-1, visual=FALSE, useHiCPro=TRUE){
  outdir <- file.path(getwd(), libname)
  
  distScaling <- 10000.0
  toKb <- 10^-3
  toMb <- 10^-6
  toProb <- 10^5
  
  useBinning <- TRUE 
  useInters <- FALSE 
  
  if (distUpThres == -1) {
    distUpThres <- Inf #no upper bound
  }
  if (distLowThres == -1) {
    distLowThres <- (-Inf) #no lower bound
  }
  

  r1 <- generate_FragPairs(fragsData, distUpThres, distLowThres)
  
  possibleInterAllCount <- r1[["possibleInterAllCount"]]
  possibleIntraAllCount <- r1[["possibleIntraAllCount"]]
  possibleIntraInRangeCount <- r1[["possibleIntraInRangeCount"]]
  
  
  listOfMappableFrags <- r1[["listOfMappableFrags"]]
  
  possiblePairsPerDistance <- r1[["possiblePairsPerDistance"]]
  
  baselineInterChrProb <- r1[["baselineInterChrProb"]]
  baselineIntraChrProb <- r1[["baselineIntraChrProb"]]
  
  biasDic <- NULL
  if (biasfile != "none") {
    biasDic <- read_ICE_biases(biasfile, fragsData, useHiCPro)
  }
  
  
  
  
  r2 <- read_All_Interactions(intersData, biasDic, listOfMappableFrags,
                              possiblePairsPerDistance, distUpThres, distLowThres)
  
  sortedInteractions <- r2[["sortedInteractions"]]
  
  observedInterAllSum <- r2[["observedInterAllSum"]]
  observedIntraAllSum <- r2[["observedIntraAllSum"]]
  observedIntraInRangeSum <- r2[["observedIntraInRangeSum"]]
  
  observedInterAllCount <- r2[["observedInterAllCount"]]
  observedIntraAllCount <- r2[["observedIntraAllCount"]]
  observedIntraInRangeCount <- r2[["observedIntraInRangeCount"]]
  
  minObservedGenomicDist <- r2[["minObservedGenomicDist"]]
  maxObservedGenomicDist <- r2[["maxObservedGenomicDist"]]
  
  possiblePairsPerDistance <- r2[["possiblePairsPerDistance"]]
  
  
  

  tempData <- calculate_Probabilities(sortedInteractions,
                                      rep(0, nrow(sortedInteractions)),
                                      paste(libname, ".fithic_pass1", sep=""), noOfBins, useBinning,
                                      distScaling, observedIntraInRangeSum, outdir, visual, libname, toKb,
                                      toProb)

  
  isOutlier <- fit_Spline(tempData$x, tempData$y, tempData$yerr, intersData,
                          sortedInteractions, biasDic, paste(libname, ".spline_pass1", sep=""),
                          1, outdir, visual, distLowThres, distUpThres, toKb, toProb, useInters,
                          baselineIntraChrProb, baselineInterChrProb, observedInterAllSum,
                          observedIntraAllSum, observedIntraInRangeSum, possibleInterAllCount,
                          possibleIntraAllCount, possibleIntraInRangeCount,
                          maxObservedGenomicDist)
  
  

  
  if (noOfPasses < 1) {
    stop("Number of passes must be greater than 0")
  }
  
  for (i in seq(2, noOfPasses + 1)) {
    tempData <- calculate_Probabilities(sortedInteractions, isOutlier,
                                        paste(libname, ".fithic_pass", i, sep=""), noOfBins, useBinning,
                                        distScaling, observedIntraInRangeSum, outdir, visual, libname,
                                        toKb, toProb)
 
    
    
    isOutlier <- fit_Spline(tempData$x, tempData$y, tempData$yerr,
                            intersData, sortedInteractions, biasDic,
                            paste(libname, ".spline_pass", i, sep=""), i, outdir, visual,
                            distLowThres, distUpThres, toKb, toProb, useInters,
                            baselineIntraChrProb, baselineInterChrProb, observedInterAllSum,
                            observedIntraAllSum, observedIntraInRangeSum, possibleInterAllCount,
                            possibleIntraAllCount, possibleIntraInRangeCount,
                            maxObservedGenomicDist)
  }
  
  
  message("Execution of Fit-Hi-C By MHiC completed successfully. [DONE]")
  
  return
}