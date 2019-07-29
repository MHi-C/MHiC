fit_Spline <- function(x, y, yerr, intersData, sortedInteractions, biasDic,
                       figname, passNo, outdir, visual, distLowThres, distUpThres, toKb, toProb,
                       useInters, baselineIntraChrProb, baselineInterChrProb, observedInterAllSum,
                       observedIntraAllSum, observedIntraInRangeSum, possibleInterAllCount,
                       possibleIntraAllCount, possibleIntraInRangeCount, maxObservedGenomicDist) {
  
  message("Running fit_Spline method ...")
  
  bias <- bias1 <- bias2 <- chr1 <- chr2 <- data <- hitCount <-
    interactionDistance <- isOutlier <- p_val <- NULL
  
  ius <- smooth.spline(x, y)
  
  tempMaxX <- max(x)
  tempMinX <- min(x)
  tempList <- unique(trunc(sortedInteractions$interactionDistance))
  
  splineX <- tempList[tempList >= tempMinX & tempList <= tempMaxX]
  splineY <- predict(ius, splineX)$y
  
  newSplineY <- monoreg(splineX, splineY, type="antitonic")$yf
  
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
  }
  
  if (visual) {
    message("Plotting ", figname, ".png")
    
    png(filename=file.path(outdir, paste(figname, ".png", sep="")),
        width=800, height=600)
    
    par(mfrow=c(2, 1))
    
    if (distLowThres > -1 & distUpThres > -1) {
      plot(splineX * toKb, newSplineY * toProb, type="l", col="green",
           xlim=c(distLowThres * toKb, distUpThres * toKb),
           xlab="Genomic distance (kb)",
           ylab=expression(paste("Contact probability (x", 10^-5, ")",
                                 sep="")))
    } else {
      plot(splineX * toKb, newSplineY * toProb, type="l", col="green",
           xlab="Genomic distance (kb)",
           ylab=expression(paste("Contact probability (x", 10^-5, ")",
                                 sep="")))
    }
    arrows(x * toKb, y * toProb - yerr * toProb, x * toKb,
           y * toProb + yerr * toProb, col="red", length=0.05, angle=90,
           code=3)
    if (useInters) {
      lines(x * toKb, rep(baselineIntraChrProb, length(x)) * toProb,
            col="black")
      lines(x * toKb, rep(baselineInterChrProb, length(x)) * toProb,
            col="blue")
      legend("topright", legend=c(paste("spline-", passNo, sep=""),
                                  "Mean with std. error", "Baseline intra-chromosomal",
                                  "Baseline inter-chromosomal"), pch=c(NA, 91, NA, NA),
             lty=c(1, NA, 1, 1), col=c("green", "red", "black", "blue"))
    } else {
      legend("topright", legend=c(paste("spline-", passNo, sep=""),
                                  "Mean with std. error"), pch=c(NA, 91),
             lty=c(1, NA), col=c("green", "red"))
    }
    
    if (distLowThres > -1 & distUpThres > -1) {
      plot(splineX, newSplineY, type="l", log="xy", col="green",
           xlim=c(distLowThres, distUpThres),
           xlab="Genomic distance (log-scale)",
           ylab="Contact probability (log-scale)")
    } else {
      plot(splineX, newSplineY, type="l", log="xy", col="green",
           xlab="Genomic distance (log-scale)",
           ylab="Contact probability (log-scale)")
    }
    arrows(x, y - yerr, x, y + yerr, col="red", length=0.05, angle=90,
           code=3)
    if (useInters) {
      lines(x, rep(baselineIntraChrProb, length(x)), col="black")
      lines(x, rep(baselineInterChrProb, length(x)), col="blue")
    }
    
    dev.off()
  }
  
  if (!is.null(biasDic)) {
    names(biasDic) <- c("chr1", "mid1", "bias1")
    data <- merge(intersData, biasDic, by=c("chr1", "mid1"), all.x=TRUE)
    
    names(biasDic) <- c("chr2", "mid2", "bias2")
    data <- merge(data, biasDic, by=c("chr2", "mid2"), all.x=TRUE)
    
    data[is.na(data)] <- 1.0
  } else {
    data <- data.table(intersData, bias1=rep(1.0, nrow(intersData)),
                       bias2=rep(1.0, nrow(intersData)))
  }
  
  ### chr1 mid1 chr2 mid2 interactionDistance hitCount bias1 bias2 p_val ###
  data <- data.table(chr1=data$chr1, mid1=data$mid1, chr2=data$chr2,
                     mid2=data$mid2, interactionDistance=abs(data$mid1 - data$mid2),
                     hitCount=data$hitCount, bias1=data$bias1, bias2=data$bias2,
                     p_val=rep(Inf, nrow(data)))
  data[chr1 == chr2 & useInters, p_val := 1 - pbinom(hitCount - 1,
                                                     observedInterAllSum, baselineInterChrProb * bias1 * bias2)]
  data[chr1 == chr2 & interactionDistance > distUpThres,
       p_val := 1 - pbinom(hitCount - 1, observedIntraAllSum, 1)]
  data[chr1 == chr2 & interactionDistance <= distLowThres, p_val := 1]
  data[chr1 == chr2 & in_range_check(interactionDistance,
                                     distLowThres, distUpThres), p_val := 1 - pbinom(hitCount - 1,
                                                                                     observedIntraInRangeSum, newSplineY[pmin(bisect_left(splineX,
                                                                                                                                          pmin(pmax(interactionDistance, tempMinX), tempMaxX)),
                                                                                                                              length(splineX))] * bias1 * bias2)]
  data[chr1 == chr2 & (bias1 < 0 | bias2 < 0), p_val := 1]
  
  p_vals <- data$p_val
  p_vals <- p_vals[p_vals != Inf]
  
  if (useInters) {
    q_vals <- benjamini_hochberg_correction(p_vals,
                                            possibleInterAllCount + possibleIntraAllCount)
  } else {
    q_vals <- benjamini_hochberg_correction(p_vals,
                                            possibleIntraInRangeCount)
  }
  
  message("Writing p-values to file ", figname, ".significances.txt.gz")
  file.create(file.path(outdir,
                        paste(figname, ".significances.txt.gz", sep="")))
  gz <- gzfile(file.path(outdir,
                         paste(figname, ".significances.txt.gz", sep="")))
  outputData <- data.table(chr1=data$chr1, fragmentMid1=data$mid1,
                           chr2=data$chr2, fragmentMid2=data$mid2, contactCount=data$hitCount,
                           p_value=p_vals, q_value=q_vals)
  if (useInters) {
    outputData <- subset(outputData, chr1 != chr2)
  } else {
    outputData <- subset(outputData, chr1 == chr2 &
                           in_range_check(abs(outputData$fragmentMid1 -
                                                outputData$fragmentMid2), distLowThres, distUpThres))
  }
  write.table(outputData, gz, quote=FALSE, sep="\t", row.names=FALSE,
              col.names=TRUE)
  
  outlierThres <- 1 / possibleIntraInRangeCount
  tempData <- data.table(
    interactionDistance=sortedInteractions$interactionDistance,
    hitCount=sortedInteractions$hitCount,
    bias=sortedInteractions$bias,
    p_val=rep(Inf, nrow(sortedInteractions)),
    isOutlier=rep(0, nrow(sortedInteractions)))
  tempData[, p_val := 1 - pbinom(hitCount - 1, observedIntraInRangeSum,
                                 newSplineY[pmin(bisect_left(splineX, pmin(pmax(interactionDistance,
                                                                                tempMinX), tempMaxX)), length(splineX))] * bias)]
  tempData[p_val < outlierThres, isOutlier := 1]
  
  belowData <- subset(tempData, isOutlier == 1)
  aboveData <- subset(tempData, isOutlier == 0)
  
  distsBelow <- belowData$interactionDistance
  distsAbove <- aboveData$interactionDistance
  
  intcountsBelow <- belowData$hitCount
  intcountsAbove <- aboveData$hitCount
  
  belowThresCount <- nrow(belowData)
  aboveThresCount <- nrow(aboveData)
  
  if (visual) {
    message("Plotting results of extracting outliers to file ", figname,
            ".extractOutliers.png")
    
    png(filename=file.path(outdir, paste(figname, ".extractOutliers.png",
                                         sep="")), width=800, height=600)
    
    downsample <- 30
    randIndcsAbove <- sample(seq(1, length(intcountsAbove)),
                             trunc(length(intcountsAbove) / downsample))
    randIndcsAbove <- randIndcsAbove[order(randIndcsAbove)]
    downsample <- 20
    randIndcsBelow <- sample(seq(1, length(intcountsBelow)),
                             trunc(length(intcountsBelow) / downsample))
    randIndcsBelow <- randIndcsBelow[order(randIndcsBelow)]
    
    if (length(intcountsBelow) > 0 &
        (distLowThres > -1 & distUpThres > -1)) {
      plot(distsBelow[randIndcsBelow] * toKb,
           intcountsBelow[randIndcsBelow], pch=20, col="red",
           xlim=c(0, distUpThres * toKb),
           ylim=c(0, min(max(intcountsBelow), 1500)),
           xlab="Genomic distance (kb)", ylab="Contact counts")
    } else if (length(intcountsBelow) > 0) {
      plot(distsBelow[randIndcsBelow] * toKb,
           intcountsBelow[randIndcsBelow], pch=20, col="red",
           ylim=c(0, min(max(intcountsBelow), 1500)),
           xlab="Genomic distance (kb)", ylab="Contact counts")
    } else if (distLowThres > -1 & distUpThres > -1) {
      plot(distsBelow[randIndcsBelow] * toKb,
           intcountsBelow[randIndcsBelow], pch=20, col="red",
           xlim=c(0, distUpThres * toKb), xlab="Genomic distance (kb)",
           ylab="Contact counts")
    } else {
      plot(distsBelow[randIndcsBelow] * toKb,
           intcountsBelow[randIndcsBelow], pch=20, col="red",
           xlab="Genomic distance (kb)", ylab="Contact counts")
    }
    lines(append(splineX, maxObservedGenomicDist) * toKb, append(newSplineY,
                                                                 newSplineY[length(newSplineY)]) * observedIntraInRangeSum,
          col="green")
    legend("topright", legend=c("Outliers (p-value < 1/M)",
                                paste("spline-", passNo, " (x N)", sep="")), pch=c(20, NA),
           lty=c(NA, 1), col=c("red", "green"))
    
    dev.off()
  }
  
  if (visual) {
    message("Plotting q-values to file ", figname, ".qplot.png")
  }
  minFDR <- 0
  maxFDR <- 0.05
  increment <- 0.001
  plot_qvalues(q_vals, minFDR, maxFDR, increment, paste(figname, ".qplot",
                                                        sep=""), outdir, visual)
  
  message("Complete fit_Spline method [OK]")
  
  return(tempData$isOutlier)
}

plot_qvalues <- function(q_values, minFDR, maxFDR, increment, figname, outdir,
                         visual) {
  
  qvalTicks <- seq(minFDR, maxFDR + increment - increment / 2, increment)
  qvalBins <- floor(q_values / increment)
  
  data <- qvalBins[qvalBins < length(qvalTicks)] + 1
  significantTicks <- as.data.frame(table(append(seq(1, length(qvalTicks)),
                                                 data)))$Freq - 1
  
  # make it cumulative
  significantTicks <- cumsum(significantTicks)
  # shift them by 1
  significantTicks <- append(0,
                             significantTicks[1 : length(significantTicks) - 1])
  
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
  }
  
  if (visual) {
    png(filename=file.path(outdir, paste(figname, ".png", sep="")),
        width=800, height=600)
    plot(qvalTicks, significantTicks, type="o", pch=17, col="blue", lty=1,
         xlab="FDR threshold", ylab="Significant contacts")
    dev.off()
  }
}

in_range_check <- function(interactionDistance, distLowThres, distUpThres) {
  return(interactionDistance > distLowThres &
           interactionDistance <= distUpThres)
}

bisect_left <- function(a, x) {
  return(length(a) + 1 - findInterval(-x, sort(-a)))
}

benjamini_hochberg_correction <- function(p_values, num_total_tests) {
  stopifnot(
    is.numeric(p_values),
    is.wholenumber(num_total_tests), length(num_total_tests) == 1L)
  
  order <- order(p_values)
  sorted_pvals <- p_values[order]
  
  bh_values <- ifelse(
    sorted_pvals == 1, 1,
    pmin(sorted_pvals * num_total_tests / seq_along(sorted_pvals), 1))
  
  return(cummax(bh_values)[order(order)])
}

is.wholenumber <- function(x, tol=.Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}