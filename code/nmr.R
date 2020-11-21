# scrip for getting all NMR results for the application

getMortResults = function(seed=123) {
  set.seed(seed)
  
  # first for the risk model
  timeSPDE = system.time(resultsSPDE <- fitSPDEKenyaDat(dat=mort, dataType="mort", 
                                                        significanceCI=.8, 
                                                        nPostSamples=1000, verbose=TRUE, seed=NULL, 
                                                        urbanEffect=TRUE, clusterEffect=TRUE, 
                                                        leaveOutRegionName=NULL, doValidation=FALSE))[3]
  
  # now aggregate
  timeAllAgg = system.time(agg <- modLCPB(uDraws=resultsSPDE$uDraws, resultsSPDE$sigmaEpsilonDraws, 
                                          includeUrban=TRUE, clusterLevel=FALSE, pixelLevel=TRUE, constituencyLevel=TRUE, countyLevel=TRUE, 
                                          regionLevel=TRUE, nationalLevel=TRUE, doModifiedPixelLevel=FALSE, 
                                          onlyDoModifiedPixelLevel=FALSE, 
                                          doLCPb=TRUE, doLCpb=TRUE, doLcpb=TRUE, urbanEffect=resultsSPDE$fixedEffectDraws[2,]))[3]
  
  # save results
  save(resultsSPDE, timeSPDE, timeAllAgg, file="savedOutput/application/mortResultsRisk.RData")
  save(agg, timeSPDE, timeAllAgg, file="savedOutput/application/mortResultsAgg.RData")
  
  invisible(NULL)
}

# Make plots for the neonatal mortality application
makeMortPlots = function() {
  # first load the model predictions
  out = load("savedOutput/application/mortResultsAgg.RData")
  browser()
  # calculate the range of predictions and CI widths
  rangePrevalencePredPixel = c()
  rangePrevalencePredConstituency = c()
  rangePrevalencePredCounty = c()
  rangePrevalencePredProvince = c()
  rangePrevalenceCIWidthPixel = c()
  rangePrevalenceCIWidthConstituency = c()
  rangePrevalenceCIWidthCounty = c()
  rangePrevalenceCIWidthProvince = c()
  
  rangeCountPredPixel = c()
  rangeCountPredConstituency = c()
  rangeCountPredCounty = c()
  rangeCountPredProvince = c()
  rangeCountCIWidthPixel = c()
  rangeCountCIWidthConstituency = c()
  rangeCountCIWidthCounty = c()
  rangeCountCIWidthProvince = c()
  
  rangeRelativePrevalencePredConstituency = c()
  rangeRelativePrevalencePredCounty = c()
  rangeRelativePrevalencePredProvince = c()
  rangeRelativePrevalenceCIWidthConstituency = c()
  rangeRelativePrevalenceCIWidthCounty = c()
  rangeRelativePrevalenceCIWidthProvince = c()
  
  # make the color scales
  areaLevels = c("pixel", "constituency", "county", "province")
  for(i in 1:length(areaLevels)) {
    thisLevel = areaLevels[i]
    
    if(thisLevel == "pixel") {
      rangePrevalencePredPixel = range(c(rangePrevalencePredPixel, 
                                         rowMeans(agg$pixelMatriceslcpb$p, na.rm=TRUE)))
      rangeCountPredPixel = range(c(rangeCountPredPixel, 
                                    rowMeans(agg$pixelMatriceslcpb$Z, na.rm=TRUE)))
      
      # threshold average population denominator per pixel to reduce noise
      prevalenceCIWidthsPixel = apply(agg$pixelMatricesLCPB$p, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})
      countCIWidthsPixel = apply(agg$pixelMatricesLCPB$Z, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})
      nPerPixel = rowSums(agg$pixelMatricesLCPB$N)
      prevalenceCIWidthsPixel[nPerPixel <= 1000*1] = NA
      countCIWidthsPixel[nPerPixel <= 1000*1 | countCIWidthsPixel < 1] = NA
      rangePrevalenceCIWidthPixel = range(prevalenceCIWidthsPixel, na.rm=TRUE)
      rangeCountCIWidthPixel = range(countCIWidthsPixel, na.rm=TRUE)
    } else if(thisLevel == "constituency") {
      rangePrevalencePredConstituency = range(c(rangePrevalencePredConstituency, 
                                                rowMeans(agg$aggregatedResultsLCpb$constituencyMatrices$p, na.rm=TRUE)))
      rangeCountPredConstituency = range(c(rangeCountPredConstituency, 
                                           rowMeans(agg$aggregatedResultsLCpb$constituencyMatrices$Z)), na.rm=TRUE)
      rangeRelativePrevalencePredConstituency = range(c(rangeRelativePrevalencePredConstituency, 
                                                        rowMeans(agg$aggregatedResultsLCPB$constituencyMatrices$pUrban/agg$aggregatedResultsLCPB$constituencyMatrices$pRural, na.rm=TRUE)))
      
      rangePrevalenceCIWidthConstituency = range(c(rangePrevalenceCIWidthConstituency, 
                                                   apply(agg$aggregatedResultsLCPB$constituencyMatrices$p, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
      rangeCountCIWidthConstituency = range(c(rangeCountCIWidthConstituency, 
                                              apply(agg$aggregatedResultsLCPB$constituencyMatrices$Z, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
      rangeRelativePrevalenceCIWidthConstituency = range(c(rangeRelativePrevalenceCIWidthConstituency, 
                                                           apply(agg$aggregatedResultsLCPB$constituencyMatrices$pUrban/agg$aggregatedResultsLCPB$constituencyMatrices$pRural, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
    } else if(thisLevel == "county") {
      rangePrevalencePredCounty = range(c(rangePrevalencePredCounty, 
                                          rowMeans(agg$aggregatedResultsLCpb$countyMatrices$p, na.rm=TRUE)))
      rangeCountPredCounty = range(c(rangeCountPredCounty, 
                                     rowMeans(agg$aggregatedResultsLCpb$countyMatrices$Z, na.rm=TRUE)))
      rangeRelativePrevalencePredCounty = range(c(rangeRelativePrevalencePredCounty, 
                                                  rowMeans(agg$aggregatedResultsLCPB$countyMatrices$pUrban/agg$aggregatedResultsLCPB$countyMatrices$pRural, na.rm=TRUE)))
      
      rangePrevalenceCIWidthCounty = range(c(rangePrevalenceCIWidthCounty, 
                                             apply(agg$aggregatedResultsLCPB$countyMatrices$p, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
      rangeCountCIWidthCounty = range(c(rangeCountCIWidthCounty, 
                                        apply(agg$aggregatedResultsLCPB$countyMatrices$Z, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
      rangeRelativePrevalenceCIWidthCounty = range(c(rangeRelativePrevalenceCIWidthCounty, 
                                                     apply(agg$aggregatedResultsLCPB$countyMatrices$pUrban/agg$aggregatedResultsLCPB$countyMatrices$pRural, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
    } else if(thisLevel == "province") {
      rangePrevalencePredProvince = range(c(rangePrevalencePredProvince, 
                                            rowMeans(agg$aggregatedResultsLCpb$regionMatrices$p, na.rm=TRUE)))
      rangeCountPredProvince = range(c(rangeCountPredProvince, 
                                       rowMeans(agg$aggregatedResultsLCpb$regionMatrices$Z, na.rm=TRUE)))
      rangeRelativePrevalencePredProvince = range(c(rangeRelativePrevalencePredProvince, 
                                                    rowMeans(agg$aggregatedResultsLCPB$regionMatrices$pUrban/agg$aggregatedResultsLCPB$regionMatrices$pRural, na.rm=TRUE)))
      
      rangePrevalenceCIWidthProvince = range(c(rangePrevalenceCIWidthProvince, 
                                             apply(agg$aggregatedResultsLCPB$regionMatrices$p, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
      rangeCountCIWidthProvince = range(c(rangeCountCIWidthProvince, 
                                        apply(agg$aggregatedResultsLCPB$regionMatrices$Z, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
      rangeRelativePrevalenceCIWidthProvince = range(c(rangeRelativePrevalenceCIWidthProvince, 
                                                     apply(agg$aggregatedResultsLCPB$regionMatrices$pUrban/agg$aggregatedResultsLCPB$regionMatrices$pRural, 1, function(x) {diff(quantile(x, probs=c(.1, .9), na.rm=TRUE))})))
    }
  }
  
  ##### make plots
  
  # load shape files for plotting
  # require(maptools)
  # regionMap = readShapePoly("data/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  # regionMap = readOGR("data/mapData/kenya_region_shapefile/kenya_region_shapefile.shp")
  kenyaMap = adm0
  countyMap = adm1
  constituencyMap = adm2
  
  # make color scales
  meanCols=makeRedBlueDivergingColors(64, rev = TRUE)
  sdCols=makeBlueYellowSequentialColors(64)
  popCols=makePurpleYellowSequentialColors(64, rev=TRUE)
  urbCols=makeGreenBlueSequentialColors(64)
  
  getWidth = function(x) {
    diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))
  }
  
  ## 2 x 2 plot of predictions (prevalence)
  
  # plot mean
  pixelMean = rowMeans(agg$pixelMatriceslcpb$p)
  constituencyMean = rowMeans(agg$aggregatedResultsLCpb$constituencyMatrices$p)
  countyMean = rowMeans(agg$aggregatedResultsLCpb$countyMatrices$p)
  provinceMean = rowMeans(agg$aggregatedResultsLCpb$regionMatrices$p)
  widthRange = range(c(rangePrevalencePredPixel, 
                      rangePrevalencePredConstituency, 
                      rangePrevalencePredCounty, 
                      rangePrevalencePredProvince))
  
  png(paste0(figDirectory, "application/prevalenceMean.png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 6.5, 1.1, 2.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){mean(x, na.rm=TRUE)}, 
             zlim=meanRange, nx=160, ny=160, main="", cex.main=3, col=meanCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # meanTicks = pretty(meanRange, n=5)[-1]
  # meanTickLabels = as.character(meanTicks)
  # image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
  #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
  #            legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  plotMapDat(plotVar=provinceMean, new = TRUE, mapDat=regionMap, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             ylab="")
  # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Posterior mean prevalence", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  # plot credible interval widths
  pixelWidth = prevalenceCIWidthsPixel
  constituencyWidth = apply(agg$aggregatedResultsLCPB$constituencyMatrices$p, 1, getWidth)
  countyWidth = apply(agg$aggregatedResultsLCPB$countyMatrices$p, 1, getWidth)
  provinceWidth = apply(agg$aggregatedResultsLCPB$regionMatrices$p, 1, getWidth)
  widthRange = range(c(rangePrevalenceCIWidthConstituency, 
                       rangePrevalenceCIWidthCounty, 
                       rangePrevalenceCIWidthProvince))
  widthRangePixel = rangePrevalenceCIWidthPixel
  
  png(paste0(figDirectory, "application/prevalenceCIWidth.png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelWidth, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(widthRangePixel), nx=160, ny=160, main="", cex.main=3, col=sdCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  widthTicksPixel = pretty(rangePrevalenceCIWidthPixel, n=8)[-c(1, 9)]
  widthTickLabelsPixel = as.character(widthTicksPixel)
  widthTicks = pretty(widthRange, n=5)
  widthTickLabels = as.character(widthTicks)
  image.plot(zlim=range(log(widthRangePixel)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
             col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksPixel), labels=widthTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyWidth, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyWidth, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  plotMapDat(plotVar=provinceWidth, new = TRUE, mapDat=regionMap, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             ylab="")
  # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Posterior prevalence 80% CI Width", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  ## 2 x 2 plot of predictions (counts)
  
  # plot means
  pixelMean = rowMeans(agg$pixelMatricesLCpb$Z)
  constituencyMean = rowMeans(agg$aggregatedResultsLCPB$constituencyMatrices$Z)
  countyMean = rowMeans(agg$aggregatedResultsLCPB$countyMatrices$Z)
  provinceMean = rowMeans(agg$aggregatedResultsLCPB$regionMatrices$Z)
  meanRangePixel = range(pixelMean[pixelMean >= .05])
  meanRangeConstituency = rangeCountPredConstituency
  meanRangeCounty = rangeCountPredCounty
  meanRangeProvince = rangeCountPredProvince
  
  png(paste0(figDirectory, "application/countMean.png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(meanRangePixel), nx=160, ny=160, main="", cex.main=3, col=meanCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  meanTicksPixel = c(1, 10, 50, 100, pretty(meanRangePixel, n=5)[-c(1, 4, 6)])
  meanTickLabelsPixel = as.character(meanTicksPixel)
  meanTicksConstituency = c(100, 500, pretty(meanRangeConstituency, n=8)[-1])
  meanTickLabelsConstituency = as.character(meanTicksConstituency)
  meanTicksCounty = c(500, 1000, pretty(meanRangeCounty, n=8)[-c(1, 6)])
  meanTickLabelsCounty = as.character(meanTicksCounty)
  meanTicksProvince = pretty(meanRangeProvince, n=8)[-1]
  meanTickLabelsProvince = as.character(meanTicksProvince)
  image.plot(zlim=range(log(meanRangePixel)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
             col=meanCols, add = TRUE, axis.args=list(at=log(meanTicksPixel), labels=meanTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.mean=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeConstituency), ticks=meanTicksConstituency, tickLabels=meanTickLabelsConstituency, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.mean=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeCounty), ticks=meanTicksCounty, tickLabels=meanTickLabelsCounty, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.mean=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  plotMapDat(plotVar=provinceMean, new = TRUE, mapDat=regionMap, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeProvince), ticks=meanTicksProvince, tickLabels=meanTickLabelsProvince, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.mean=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             ylab="")
  # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Posterior mean total deaths", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  # plot credible interval means
  pixelWidth = countCIWidthsPixel
  constituencyWidth = apply(agg$aggregatedResultsLCPB$constituencyMatrices$Z, 1, getWidth)
  countyWidth = apply(agg$aggregatedResultsLCPB$countyMatrices$Z, 1, getWidth)
  provinceWidth = apply(agg$aggregatedResultsLCPB$regionMatrices$Z, 1, getWidth)
  widthRangePixel = rangeCountCIWidthPixel # range(pixelWidth)
  widthRangeConstituency = rangeCountCIWidthConstituency
  widthRangeCounty = rangeCountCIWidthCounty
  widthRangeProvince = rangeCountCIWidthProvince
  
  png(paste0(figDirectory, "application/countWidth.png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(widthRangePixel), nx=160, ny=160, main="", cex.main=3, col=sdCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  widthTicksPixel = c(1, 10, 50, 100, pretty(widthRangePixel, n=5)[-c(1)])
  widthTickLabelsPixel = as.character(widthTicksPixel)
  widthTicksConstituency = c(100, pretty(widthRangeConstituency, n=8)[-1])
  widthTickLabelsConstituency = as.character(widthTicksConstituency)
  widthTicksCounty = c(500, pretty(widthRangeCounty, n=5)[-c(1)])
  widthTickLabelsCounty = as.character(widthTicksCounty)
  widthTicksProvince = pretty(widthRangeProvince, n=8)
  widthTickLabelsProvince = as.character(widthTicksProvince)
  image.plot(zlim=range(log(widthRangePixel)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
             col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksPixel), labels=widthTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeConstituency), ticks=widthTicksConstituency, tickLabels=widthTickLabelsConstituency, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeCounty), ticks=widthTicksCounty, tickLabels=widthTickLabelsCounty, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  plotMapDat(plotVar=provinceWidth, new = TRUE, mapDat=regionMap, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeProvince), ticks=widthTicksProvince, tickLabels=widthTickLabelsProvince, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             ylab="")
  # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Total deaths 80% CI width", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  ## 1 x 3 plot of predictions (relative prevalence)
  
  # plot mean
  pixelMean = rowMeans(agg$pixelMatriceslcpb$p)
  constituencyMean = rowMeans(agg$aggregatedResultsLCpb$constituencyMatrices$p)
  countyMean = rowMeans(agg$aggregatedResultsLCpb$countyMatrices$p)
  provinceMean = rowMeans(agg$aggregatedResultsLCpb$regionMatrices$p)
  widthRange = range(c(rangePrevalencePredPixel, 
                       rangePrevalencePredConstituency, 
                       rangePrevalencePredCounty, 
                       rangePrevalencePredProvince))
  
  png(paste0(figDirectory, "application/relativePrevalence.png"), width=1500, height=1000)
  par(mfrow=c(2,3), oma=c( 0,4,4,7), mar=c(6.1, 6.5, 1.1, 2.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){mean(x, na.rm=TRUE)}, 
             zlim=meanRange, nx=160, ny=160, main="", cex.main=3, col=meanCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # meanTicks = pretty(meanRange, n=5)[-1]
  # meanTickLabels = as.character(meanTicks)
  # image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
  #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
  #            legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  plotMapDat(plotVar=provinceMean, new = TRUE, mapDat=regionMap, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             ylab="")
  # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Posterior mean prevalence", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
}



