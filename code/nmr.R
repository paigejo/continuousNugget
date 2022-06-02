# script for getting all NMR results for the application
# 25GB
getMortResults = function(seed=123, useCoarseGrid=FALSE, logisticApproximation=FALSE, nPostSamples=1000) {
  time1 = proc.time()[3]
  set.seed(seed)
  
  easpa = makeDefaultEASPA()
  
  poppsub = poppsubKenya
  
  if(useCoarseGrid) {
    popMat = popGridCoarse
    popMatAdjusted = popGridCoarseAdjusted
  } else {
    popMat = popGrid
    popMatAdjusted = popGridAdjusted
  }
  
  # Fit the risk model:
  time2 = proc.time()[3]
  riskOut = fitSPDEKenyaDat(dat=mort, nPostSamples=nPostSamples, popMat=popMat)
  
  # remove large parts of the INLA model we don't need
  riskOut$mod$all.hyper = NULL
  riskOut$mod$.args = NULL
  
  logitDraws = riskOut$uDraws # uDraws includes fixed effects and is on logit scale
  sigmaEpsilonDraws = riskOut$sigmaEpsilonDraws
  time3 = proc.time()[3]
  
  # run the models
  out = simPopCustom(logitRiskDraws=logitDraws, sigmaEpsilonDraws=sigmaEpsilonDraws, 
                     easpa=easpa, popMat=popMat, 
                     targetPopMat=popMatAdjusted, poppsub=poppsub, 
                     stratifyByUrban=TRUE, subareaLevel=TRUE, gridLevel=TRUE, 
                     doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
                     doGriddedRisk=FALSE, doSmoothRiskLogisticApprox=logisticApproximation, 
                     min1PerSubarea=TRUE)
  
  pixelPop = out$pixelPop
  subareaPop = out$subareaPop$aggregationResults
  areaPop = out$areaPop$aggregationResults
  allTimings = out$allTimings
  processedTimings = out$processedTimings
  aggregationTimings=list(allTimings=allTimings, processedTimings=processedTimings)
  time4 = proc.time()[3]
  
  # get timings
  rawTimes = c(time1, time2, time3, time4)
  totalTimes = diff(rawTimes)
  totalTimes = c(totalTimes, sum(totalTimes))
  names(totalTimes) = c("setup", "SPDEmodel", "aggregationModel", "totalTime")
  
  # save file
  fileName = paste0("savedOutput/application/finalMort", ifelse(useCoarseGrid, "Coarse", ""), ".RData")
  save(riskOut, pixelPop, subareaPop, areaPop, aggregationTimings, rawTimes, totalTimes, file=fileName)
  
  list(riskOut=riskOut, subareaPopAggRes=subareaPop, 
       areaPopAggRes=areaPop, 
       aggregationTimings=aggregationTimings, 
       rawTimes=rawTimes, totalTimes=totalTimes)
}

# Make plots for the neonatal mortality application
makeMortPlots = function(logisticApproximation=FALSE, signif=.95) {
  alpha = 1 - signif
  
  # first load the model predictions
  logisticText = ifelse(!logisticApproximation, "", "logisticApprox")
  out = load(paste0("savedOutput/application/finalMort", logisticText, ".RData"))
  
  # calculate the range of predictions and CI widths
  rangePrevalencePredPixel = c()
  rangePrevalencePredConstituency = c()
  rangePrevalencePredCounty = c()
  rangePrevalencePredProvince = c()
  rangePrevalenceCIWidthPixel = c()
  rangePrevalenceCIWidthConstituency = c()
  rangePrevalenceCIWidthCounty = c()
  rangePrevalenceCIWidthProvince = c()
  rangePrevalenceSDPixel = c()
  rangePrevalenceSDConstituency = c()
  rangePrevalenceSDCounty = c()
  rangePrevalenceSDProvince = c()
  
  rangeCountPredPixel = c()
  rangeCountPredConstituency = c()
  rangeCountPredCounty = c()
  rangeCountPredProvince = c()
  rangeCountCIWidthPixel = c()
  rangeCountCIWidthConstituency = c()
  rangeCountCIWidthCounty = c()
  rangeCountCIWidthProvince = c()
  rangeCountSDPixel = c()
  rangeCountSDConstituency = c()
  rangeCountSDCounty = c()
  rangeCountSDProvince = c()
  
  rangeRelativePrevalencePredConstituency = c()
  rangeRelativePrevalencePredCounty = c()
  rangeRelativePrevalencePredProvince = c()
  rangeRelativePrevalenceCIWidthConstituency = c()
  rangeRelativePrevalenceCIWidthCounty = c()
  rangeRelativePrevalenceCIWidthProvince = c()
  rangeRelativePrevalenceSDConstituency = c()
  rangeRelativePrevalenceSDCounty = c()
  rangeRelativePrevalenceSDProvince = c()
  
  # make the color scales
  
  areaLevels = c("pixel", "constituency", "county") # province results no longer supported
  for(i in 1:length(areaLevels)) {
    thisLevel = areaLevels[i]
    
    if(thisLevel == "pixel") {
      rangePrevalencePredPixel = range(c(rangePrevalencePredPixel, 
                                         rowMeans(pixelPop$pSmoothRisk, na.rm=TRUE)))
      rangeCountPredPixel = range(c(rangeCountPredPixel, 
                                    rowMeans(pixelPop$ZSmoothRisk, na.rm=TRUE)))
      
      # threshold average population denominator per pixel to reduce noise
      prevalenceCIWidthPixel = apply(pixelPop$pFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthPixel = apply(pixelPop$ZFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      nPerPixel = rowSums(pixelPop$NFineScalePrevalence)
      prevalenceCIWidthPixel[nPerPixel <= 1000*1] = NA
      countCIWidthPixel[nPerPixel <= 1000*1 | countCIWidthPixel < 1] = NA
      rangePrevalenceCIWidthPixel = range(prevalenceCIWidthPixel, na.rm=TRUE)
      rangeCountCIWidthPixel = range(countCIWidthPixel, na.rm=TRUE)
      
      prevalenceSDPixel = apply(pixelPop$pFineScalePrevalence, 1, sd, na.rm=TRUE)
      countSDPixel = apply(pixelPop$ZFineScalePrevalence, 1, sd, na.rm=TRUE)
      nPerPixel = rowSums(pixelPop$NFineScalePrevalence)
      prevalenceSDPixel[nPerPixel <= 1000*1] = NA
      countSDPixel[nPerPixel <= 1000*1 | countSDPixel < 1] = NA
      rangePrevalenceSDPixel = range(prevalenceSDPixel, na.rm=TRUE)
      rangeCountSDPixel = range(countSDPixel, na.rm=TRUE)
      
      # do the same but for the lcpb model
      prevalenceCIWidthPixellcpb = apply(pixelPop$pSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthPixellcpb = apply(pixelPop$ZSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      # prevalenceCIWidthPixellcpb[nPerPixel <= 1000*1] = NA
      # countCIWidthPixellcpb[nPerPixel <= 1000*1 | is.na(countCIWidthPixel)] = NA
      rangePrevalenceCIWidthPixellcpb = range(prevalenceCIWidthPixellcpb, na.rm=TRUE)
      rangeCountCIWidthPixellcpb = range(countCIWidthPixellcpb, na.rm=TRUE)
      
      prevalenceSDPixellcpb = apply(pixelPop$pSmoothRisk, 1, sd, na.rm=TRUE)
      countSDPixellcpb = apply(pixelPop$ZSmoothRisk, 1, sd, na.rm=TRUE)
      prevalenceSDPixellcpb[nPerPixel <= 1000*1] = NA
      countSDPixellcpb[nPerPixel <= 1000*1 | is.na(countSDPixel)] = NA
      rangePrevalenceSDPixellcpb = range(prevalenceSDPixellcpb, na.rm=TRUE)
      rangeCountSDPixellcpb = range(countSDPixellcpb, na.rm=TRUE)
      
      # do the same but for the LCPb model
      prevalenceCIWidthPixelLCPb = apply(pixelPop$pFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthPixelLCPb = apply(pixelPop$ZFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      prevalenceCIWidthPixelLCPb[nPerPixel <= 1000*1] = NA
      countCIWidthPixelLCPb[nPerPixel <= 1000*1 | is.na(countCIWidthPixel)] = NA
      rangePrevalenceCIWidthPixelLCPb = range(prevalenceCIWidthPixelLCPb, na.rm=TRUE)
      rangeCountCIWidthPixelLCPb = range(countCIWidthPixelLCPb, na.rm=TRUE)
      
      prevalenceSDPixelLCPb = apply(pixelPop$pFineScaleRisk, 1, sd, na.rm=TRUE)
      countSDPixelLCPb = apply(pixelPop$ZFineScaleRisk, 1, sd, na.rm=TRUE)
      prevalenceSDPixelLCPb[nPerPixel <= 1000*1] = NA
      countSDPixelLCPb[nPerPixel <= 1000*1 | is.na(countSDPixel)] = NA
      rangePrevalenceSDPixelLCPb = range(prevalenceSDPixelLCPb, na.rm=TRUE)
      rangeCountSDPixelLCPb = range(countSDPixelLCPb, na.rm=TRUE)
    } else if(thisLevel == "constituency") {
      urbanConstituencies = poppcon$County == "Nairobi" | poppcon$County == "Mombasa"
      undefinedRelativePrevalenceConstituencies = (poppcon$popUrb == 0) | (poppcon$popRur == 0)
      rangePrevalencePredConstituency = range(c(rangePrevalencePredConstituency, 
                                                rowMeans(subareaPop$pFineScaleRisk, na.rm=TRUE)))
      rangeCountPredConstituency = range(c(rangeCountPredConstituency, 
                                           rowMeans(subareaPop$ZFineScaleRisk)), na.rm=TRUE)
      relativePrevalencePredConstituency = rowMeans(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, na.rm=TRUE)
      relativePrevalencePredConstituency[undefinedRelativePrevalenceConstituencies] = NA
      rangeRelativePrevalencePredConstituency = range(relativePrevalencePredConstituency[is.finite(relativePrevalencePredConstituency)], na.rm=TRUE)
      
      prevalenceCIWidthConstituency = apply(subareaPop$pFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthConstituency = apply(subareaPop$ZFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      rangePrevalenceCIWidthConstituency = range(prevalenceCIWidthConstituency)
      rangeCountCIWidthConstituency = range(countCIWidthConstituency)
      relativePrevalenceCIWidthConstituency = apply(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthConstituency[undefinedRelativePrevalenceConstituencies] = NA
      rangeRelativePrevalenceCIWidthConstituency = range(relativePrevalenceCIWidthConstituency[is.finite(relativePrevalenceCIWidthConstituency)], na.rm=TRUE)
      
      prevalenceSDConstituency = apply(subareaPop$pFineScalePrevalence, 1, sd, na.rm=TRUE)
      countSDConstituency = apply(subareaPop$ZFineScalePrevalence, 1, sd, na.rm=TRUE)
      rangePrevalenceSDConstituency = range(prevalenceSDConstituency)
      rangeCountSDConstituency = range(countSDConstituency)
      relativePrevalenceSDConstituency = apply(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, 1, sd, na.rm=TRUE)
      relativePrevalenceSDConstituency[undefinedRelativePrevalenceConstituencies] = NA
      rangeRelativePrevalenceSDConstituency = range(relativePrevalenceSDConstituency[is.finite(relativePrevalenceSDConstituency)], na.rm=TRUE)
      
      # get credible interval widths for the lcpb model
      prevalenceCIWidthConstituencylcpb = apply(subareaPop$pSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthConstituencylcpb = apply(subareaPop$ZSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthConstituencylcpb = apply(subareaPop$pUrbanSmoothRisk/subareaPop$pRuralSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthConstituencylcpb[undefinedRelativePrevalenceConstituencies] = NA
      
      prevalenceSDConstituencylcpb = apply(subareaPop$pSmoothRisk, 1, sd, na.rm=TRUE)
      countSDConstituencylcpb = apply(subareaPop$ZSmoothRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDConstituencylcpb = apply(subareaPop$pUrbanSmoothRisk/subareaPop$pRuralSmoothRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDConstituencylcpb[undefinedRelativePrevalenceConstituencies] = NA
      
      # do the same for the LCPb model
      prevalenceCIWidthConstituencyLCPb = apply(subareaPop$pFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthConstituencyLCPb = apply(subareaPop$ZFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthConstituencyLCPb = apply(subareaPop$pUrbanFineScaleRisk/subareaPop$pRuralFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthConstituencyLCPb[undefinedRelativePrevalenceConstituencies] = NA
      
      prevalenceSDConstituencyLCPb = apply(subareaPop$pFineScaleRisk, 1, sd, na.rm=TRUE)
      countSDConstituencyLCPb = apply(subareaPop$ZFineScaleRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDConstituencyLCPb = apply(subareaPop$pUrbanFineScaleRisk/subareaPop$pRuralFineScaleRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDConstituencyLCPb[undefinedRelativePrevalenceConstituencies] = NA
    } else if(thisLevel == "county") {
      urbanCounties = sort(poppc$County) == "Nairobi" | sort(poppc$County) == "Mombasa"
      rangePrevalencePredCounty = range(c(rangePrevalencePredCounty, 
                                          rowMeans(areaPop$pFineScaleRisk, na.rm=TRUE)))
      rangeCountPredCounty = range(c(rangeCountPredCounty, 
                                     rowMeans(areaPop$ZFineScaleRisk, na.rm=TRUE)))
      relativePrevalencePredCounty = rowMeans(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence, na.rm=TRUE)
      relativePrevalencePredCounty[urbanCounties] = NA
      rangeRelativePrevalencePredCounty = range(relativePrevalencePredCounty[is.finite(relativePrevalencePredCounty)], na.rm=TRUE)
      
      prevalenceCIWidthCounty = apply(areaPop$pFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthCounty = apply(areaPop$ZFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      rangePrevalenceCIWidthCounty = range(prevalenceCIWidthCounty)
      rangeCountCIWidthCounty = range(countCIWidthCounty)
      relativePrevalenceCIWidthCounty = apply(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthCounty[urbanCounties] = NA
      rangeRelativePrevalenceCIWidthCounty = range(relativePrevalenceCIWidthCounty[is.finite(relativePrevalenceCIWidthCounty)], na.rm=TRUE)
      
      prevalenceSDCounty = apply(areaPop$pFineScalePrevalence, 1, sd, na.rm=TRUE)
      countSDCounty = apply(areaPop$ZFineScalePrevalence, 1, sd, na.rm=TRUE)
      rangePrevalenceSDCounty = range(prevalenceSDCounty)
      rangeCountSDCounty = range(countSDCounty)
      relativePrevalenceSDCounty = apply(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence, 1, sd, na.rm=TRUE)
      relativePrevalenceSDCounty[urbanCounties] = NA
      rangeRelativePrevalenceSDCounty = range(relativePrevalenceSDCounty[is.finite(relativePrevalenceSDCounty)], na.rm=TRUE)
      
      # get credible interval widths for the lcpb model
      prevalenceCIWidthCountylcpb = apply(areaPop$pSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthCountylcpb = apply(areaPop$ZSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthCountylcpb = apply(areaPop$pUrbanSmoothRisk/areaPop$pRuralSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthCountylcpb[urbanCounties] = NA
      
      prevalenceSDCountylcpb = apply(areaPop$pSmoothRisk, 1, sd, na.rm=TRUE)
      countSDCountylcpb = apply(areaPop$ZSmoothRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDCountylcpb = apply(areaPop$pUrbanSmoothRisk/areaPop$pRuralSmoothRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDCountylcpb[urbanCounties] = NA
      
      # get credible interval widths for the LCPb model
      prevalenceCIWidthCountyLCPb = apply(areaPop$pFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthCountyLCPb = apply(areaPop$ZFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthCountyLCPb = apply(areaPop$pUrbanFineScaleRisk/areaPop$pRuralFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthCountyLCPb[urbanCounties] = NA
      
      prevalenceSDCountyLCPb = apply(areaPop$pFineScaleRisk, 1, sd, na.rm=TRUE)
      countSDCountyLCPb = apply(areaPop$ZFineScaleRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDCountyLCPb = apply(areaPop$pUrbanFineScaleRisk/areaPop$pRuralFineScaleRisk, 1, sd, na.rm=TRUE)
      relativePrevalenceSDCountyLCPb[urbanCounties] = NA
    } else if(thisLevel == "province") {
      urbanProvinces = sort(regionMap@data$NAME_1) == "Nairobi"
      rangePrevalencePredProvince = range(c(rangePrevalencePredProvince, 
                                            rowMeans(agg$aggregatedResultsLCPb$regionMatrices$p, na.rm=TRUE)))
      rangeCountPredProvince = range(c(rangeCountPredProvince, 
                                       rowMeans(agg$aggregatedResultsLCPb$regionMatrices$Z, na.rm=TRUE)))
      relativePrevalencePredProvince = rowMeans(agg$aggregatedResultsLCPB$regionMatrices$pUrban/agg$aggregatedResultsLCPB$regionMatrices$pRural, na.rm=TRUE)
      relativePrevalencePredProvince[urbanProvinces] = NA
      rangeRelativePrevalencePredProvince = range(relativePrevalencePredProvince[is.finite(relativePrevalencePredProvince)], na.rm=TRUE)
      
      prevalenceCIWidthProvince = apply(agg$aggregatedResultsLCPB$regionMatrices$p, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthProvince = apply(agg$aggregatedResultsLCPB$regionMatrices$Z, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      rangePrevalenceCIWidthProvince = range(prevalenceCIWidthProvince)
      rangeCountCIWidthProvince = range(countCIWidthProvince)
      relativePrevalenceCIWidthProvince = apply(agg$aggregatedResultsLCPB$regionMatrices$pUrban/agg$aggregatedResultsLCPB$regionMatrices$pRural, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthProvince[urbanProvinces] = NA
      rangeRelativePrevalenceCIWidthProvince = range(relativePrevalenceCIWidthProvince[is.finite(relativePrevalenceCIWidthProvince)], na.rm=TRUE)
      
      prevalenceSDProvince = apply(agg$aggregatedResultsLCPB$regionMatrices$p, 1, sd, na.rm=TRUE)
      countSDProvince = apply(agg$aggregatedResultsLCPB$regionMatrices$Z, 1, sd, na.rm=TRUE)
      rangePrevalenceSDProvince = range(prevalenceSDProvince)
      rangeCountSDProvince = range(countSDProvince)
      relativePrevalenceSDProvince = apply(agg$aggregatedResultsLCPB$regionMatrices$pUrban/agg$aggregatedResultsLCPB$regionMatrices$pRural, 1, sd, na.rm=TRUE)
      relativePrevalenceSDProvince[urbanProvinces] = NA
      rangeRelativePrevalenceSDProvince = range(relativePrevalenceSDProvince[is.finite(relativePrevalenceSDProvince)], na.rm=TRUE)
      
      # now calculate credible interval widths for the lcpb model
      prevalenceCIWidthProvincelcpb = apply(agg$aggregatedResultslcpb$regionMatrices$p, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthProvincelcpb = apply(agg$aggregatedResultslcpb$regionMatrices$Z, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthProvincelcpb = apply(agg$aggregatedResultslcpb$regionMatrices$pUrban/agg$aggregatedResultslcpb$regionMatrices$pRural, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthProvincelcpb[urbanProvinces] = NA
      
      prevalenceSDProvincelcpb = apply(agg$aggregatedResultslcpb$regionMatrices$p, 1, sd, na.rm=TRUE)
      countSDProvincelcpb = apply(agg$aggregatedResultslcpb$regionMatrices$Z, 1, sd, na.rm=TRUE)
      relativePrevalenceSDProvincelcpb = apply(agg$aggregatedResultslcpb$regionMatrices$pUrban/agg$aggregatedResultslcpb$regionMatrices$pRural, 1, sd, na.rm=TRUE)
      relativePrevalenceSDProvincelcpb[urbanProvinces] = NA
      
      # now calculate credible interval widths for the LCPb model
      prevalenceCIWidthProvinceLCPb = apply(agg$aggregatedResultsLCPb$regionMatrices$p, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      countCIWidthProvinceLCPb = apply(agg$aggregatedResultsLCPb$regionMatrices$Z, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthProvinceLCPb = apply(agg$aggregatedResultsLCPb$regionMatrices$pUrban/agg$aggregatedResultsLCPb$regionMatrices$pRural, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))})
      relativePrevalenceCIWidthProvinceLCPb[urbanProvinces] = NA
      
      prevalenceSDProvinceLCPb = apply(agg$aggregatedResultsLCPb$regionMatrices$p, 1, sd, na.rm=TRUE)
      countSDProvinceLCPb = apply(agg$aggregatedResultsLCPb$regionMatrices$Z, 1, sd, na.rm=TRUE)
      relativePrevalenceSDProvinceLCPb = apply(agg$aggregatedResultsLCPb$regionMatrices$pUrban/agg$aggregatedResultsLCPb$regionMatrices$pRural, 1, sd, na.rm=TRUE)
      relativePrevalenceSDProvinceLCPb[urbanProvinces] = NA
    }
  }
  
  ##### make plots
  
  # load shape files for plotting
  # require(maptools)
  # regionMap = readShapePoly("data/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  # regionMap = readOGR("data/mapData/kenya_region_shapefile/kenya_region_shapefile.shp")
  kenyaMap = adm0
  countyMap = adm1compressed
  constituencyMap = adm2compressed
  
  # make color scales
  meanCols=makeRedBlueDivergingColors(64, rev = TRUE)
  sdCols=makeBlueGreenYellowSequentialColors(64)
  popCols=makePurpleYellowSequentialColors(64, rev=TRUE)
  urbCols=makeGreenBlueSequentialColors(64)
  
  getWidth = function(x) {
    diff(quantile(x, prob=c(alpha/2, 1-alpha/2), na.rm=TRUE))
  }
  
  ## 3 x 2 plot of predictions (prevalence mean and CI width)
  
  # plot mean
  pixelMean = rowMeans(pixelPop$pSmoothRisk)
  constituencyMean = rowMeans(subareaPop$pFineScaleRisk)
  countyMean = rowMeans(areaPop$pFineScaleRisk)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPb$regionMatrices$p)
  meanRange = range(pixelMean, constituencyMean, countyMean)
  # widthRange = range(c(rangePrevalenceCIWidthPixel, 
  #                     rangePrevalenceCIWidthConstituency, 
  #                     rangePrevalenceCIWidthCounty))
  
  # browser()
  png(paste0(figDirectory, "application/prevalenceMeanCIWidth", logisticText, ".png"), width=1500, height=1000)
  par(mfrow=c(2,3), oma=c(5,5,5,5), mar=c(3.1, 5.1, 1.1, 7.1))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){logit(mean(x, na.rm=TRUE))}, 
             zlim=logit(meanRange), nx=160, ny=160, main="", cex.main=3, col=meanCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  meanTicksPixel = pretty(rangePrevalencePredPixel, n=8)[-c(1, 9)]
  meanTickLabelsPixel = as.character(meanTicksPixel)
  meanTicks = pretty(meanRange, n=5)
  meanTickLabels = as.character(meanTicks)
  # image.plot(zlim=range(log(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
  #            col=meanCols, add = TRUE, axis.args=list(at=log(meanTicksPixel), labels=meanTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
  #            legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # pixel level
  
  quilt.plot(popGrid$lon, popGrid$lat, prevalenceCIWidthPixellcpb, FUN=function(x){logit(mean(x, na.rm=TRUE))}, 
             nx=160, ny=160, main="", cex.main=3, col=sdCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  widthTicksPixel = pretty(rangePrevalenceCIWidthPixel, n=8)[-c(1, 9)]
  widthTickLabelsPixel = as.character(widthTicksPixel)
  widthRangePixel = rangePrevalenceCIWidthPixel
  image.plot(zlim=range(logit(widthRangePixel)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
             col=sdCols, add = TRUE, axis.args=list(at=logit(widthTicksPixel), labels=widthTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.width=3, smallplot=c(.88,.91,.1,.9))
  
  # constituency level
  widthTicksSubarea = pretty(rangePrevalenceCIWidthConstituency, n=8)
  widthTickLabelsSubarea = as.character(widthTicksSubarea)
  widthRangeSubarea = rangePrevalenceCIWidthConstituency
  plotMapDat(plotVar=prevalenceCIWidthConstituency, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=logit, scaleFunInverse=expit, 
             cols=sdCols, ticks=widthTicksSubarea, tickLabels=widthTickLabelsSubarea, 
             zlim=logit(widthRangeSubarea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # county level
  widthTicksArea = pretty(rangePrevalenceCIWidthCounty, n=8)
  widthTickLabelsArea = as.character(widthTicksArea)
  widthRangeArea = rangePrevalenceCIWidthCounty
  plotMapDat(plotVar=prevalenceCIWidthCounty, new = TRUE, 
             main="", scaleFun=logit, scaleFunInverse=expit, 
             cols=sdCols, ticks=widthTicksArea, tickLabels=widthTickLabelsArea, 
             zlim=logit(widthRangeArea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # add title in the top margin
  mtext(side = 2, paste0("Mean Prevalence"), line = 2.5, cex=2.5, outer=TRUE, at=.76)
  mtext(side = 2, paste0(100*signif, "% CI Width"), line = 2.5, cex=2.5, outer=TRUE, at=.26)
  
  mtext(side = 1, "Longitude", line = 2.5, cex=2, outer=TRUE, at=.167)
  mtext(side = 1, "Longitude", line = 2.5, cex=2, outer=TRUE, at=.5)
  mtext(side = 1, "Longitude", line = 2.5, cex=2, outer=TRUE, at=.833)
  
  mtext(side = 3, "5km Pixel Level", line = 2, cex=2.5, outer=TRUE, at=.167)
  mtext(side = 3, "Constituency Level", line = 2, cex=2.5, outer=TRUE, at=.5)
  mtext(side = 3, "County Level", line = 2, cex=2.5, outer=TRUE, at=.833)
  
  dev.off()
  
  browser()
  # same plot but for burden
  pixelMean = rowMeans(pixelPop$ZSmoothRisk)
  constituencyMean = rowMeans(subareaPop$ZFineScaleRisk)
  countyMean = rowMeans(areaPop$ZFineScaleRisk)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPb$regionMatrices$p)
  meanRange = c(.01, max(c(pixelMean, constituencyMean, countyMean)))
  meanRangePixel = c(.01, max(pixelMean))
  meanRangeConstituency = range(constituencyMean)
  meanRangeCounty = range(countyMean)
  
  png(paste0(figDirectory, "application/burdenMeanCIWidth", logisticText, ".png"), width=1500, height=1000)
  par(mfrow=c(2,3), oma=c(5,5,5,5), mar=c(3.1, 5.1, 1.1, 7.1))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(meanRangePixel), nx=160, ny=160, main="", cex.main=3, col=meanCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  meanTicksPixel = getLogScaleTicks(meanRangePixel)
  meanTickLabelsPixel = as.character(meanTicksPixel)
  meanTicks = pretty(meanRange, n=5)
  meanTickLabels = as.character(meanTicks)
  image.plot(zlim=range(log(meanRangePixel)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
             col=meanCols, add = TRUE, axis.args=list(at=log(meanTicksPixel), labels=meanTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.width=3, smallplot=c(.88,.91,.1,.9))
  
  # constituency level
  meanTicksConstituency = getLogScaleTicks(meanRangeConstituency)
  meanTickLabelsConstituency = as.character(meanTicksConstituency)
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeConstituency), ticks=meanTicksConstituency, tickLabels=meanTickLabelsConstituency, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar=TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # county level
  meanTicksCounty = getLogScaleTicks(meanRangeCounty)
  meanTickLabelsCounty = as.character(meanTicksCounty)
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeCounty), ticks=meanTicksCounty, tickLabels=meanTickLabelsCounty, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # pixel level
  
  quilt.plot(popGrid$lon, popGrid$lat, countCIWidthPixellcpb, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             nx=160, ny=160, main="", cex.main=3, col=sdCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  widthTicksPixel = getLogScaleTicks(c(.01, rangeCountCIWidthPixellcpb[2]))
  widthTickLabelsPixel = as.character(widthTicksPixel)
  widthRangePixel = rangeCountCIWidthPixel
  image.plot(zlim=range(log(widthRangePixel)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
             col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksPixel), labels=widthTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.width=3, smallplot=c(.88,.91,.1,.9))
  
  # constituency level
  widthTicksSubarea = getLogScaleTicks(rangeCountCIWidthConstituency)
  widthTickLabelsSubarea = as.character(widthTicksSubarea)
  widthRangeSubarea = rangeCountCIWidthConstituency
  plotMapDat(plotVar=countCIWidthConstituency, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, ticks=widthTicksSubarea, tickLabels=widthTickLabelsSubarea, 
             zlim=log(widthRangeSubarea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # county level
  widthTicksArea = getLogScaleTicks(rangeCountCIWidthCounty)
  widthTickLabelsArea = as.character(widthTicksArea)
  widthRangeArea = rangeCountCIWidthCounty
  plotMapDat(plotVar=countCIWidthCounty, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, ticks=widthTicksArea, tickLabels=widthTickLabelsArea, 
             zlim=log(widthRangeArea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # add title in the top margin
  mtext(side = 2, paste0("Mean Burden"), line = 2.5, cex=2.5, outer=TRUE, at=.76)
  mtext(side = 2, paste0(100*signif, "% CIs"), line = 2.5, cex=2.5, outer=TRUE, at=.26)
  
  mtext(side = 1, "Longitude", line = 2.5, cex=2, outer=TRUE, at=.167)
  mtext(side = 1, "Longitude", line = 2.5, cex=2, outer=TRUE, at=.5)
  mtext(side = 1, "Longitude", line = 2.5, cex=2, outer=TRUE, at=.833)
  
  mtext(side = 3, "5km Pixel Level", line = 2, cex=2.5, outer=TRUE, at=.167)
  mtext(side = 3, "Constituency Level", line = 2, cex=2.5, outer=TRUE, at=.5)
  mtext(side = 3, "County Level", line = 2, cex=2.5, outer=TRUE, at=.833)
  
  dev.off()
  
  # same plot but for relative prevalence. Don't include pixel level because 
  # pixel is either urban or rural, not both
  browser()
  constituencyMean = rowMeans(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence)
  countyMean = rowMeans(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPb$regionMatrices$p)
  meanRange = range(c(constituencyMean, countyMean)[is.finite(c(constituencyMean, countyMean))], na.rm=TRUE)
  meanRangeConstituency = range(constituencyMean[is.finite(constituencyMean)], na.rm=TRUE)
  meanRangeCounty = range(countyMean, na.rm=TRUE)
  meanTicks = pretty(meanRange, n=5)
  meanTickLabels = as.character(meanTicks)
  
  urbDivergingCols = makeGreenBlueDivergingColors(64, center=1, valRange=meanRange)
  urbDivergingCols = centerColorScale(64, valRange=meanRange, center=1, colScale=makeGreenBlueDivergingColors, scaleFun=log)
  
  png(paste0(figDirectory, "application/relativePrevMeanCIWidth", logisticText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c(3,5,3,3), mar=c(3.1, 5.1, 1.1, 5.1))
  
  # constituency level
  meanTicksConstituency = pretty(meanRangeConstituency, n=10)
  meanTickLabelsConstituency = as.character(meanTicksConstituency)
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=urbDivergingCols, zlim=log(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar=FALSE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2.5), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="Latitude")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # county level
  meanTicksCounty = pretty(meanRangeCounty, n=3)
  meanTickLabelsCounty = as.character(meanTicksCounty)
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=urbDivergingCols, zlim=log(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # constituency level
  widthTicksSubarea = c(getLogScaleTicks(rangeRelativePrevalenceCIWidthConstituency), 4)
  widthTickLabelsSubarea = as.character(widthTicksSubarea)
  widthRangeSubarea = rangeRelativePrevalenceCIWidthConstituency
  plotMapDat(plotVar=relativePrevalenceCIWidthConstituency, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, ticks=widthTicksSubarea, tickLabels=widthTickLabelsSubarea, 
             zlim=log(widthRangeSubarea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2.5), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="Latitude")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # county level
  widthTicksArea = getLogScaleTicks(rangeRelativePrevalenceCIWidthCounty)
  widthTickLabelsArea = as.character(widthTicksArea)
  widthRangeArea = rangeRelativePrevalenceCIWidthCounty
  plotMapDat(plotVar=relativePrevalenceCIWidthCounty, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, ticks=widthTicksArea, tickLabels=widthTickLabelsArea, 
             zlim=log(widthRangeArea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot=c(.88,.91,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # add title in the top margin
  mtext(side = 2, paste0("Mean Relative Prevalence"), line = 2.5, cex=2.5, outer=TRUE, at=.76)
  mtext(side = 2, paste0(100*signif, "% CIs"), line = 2.5, cex=2.5, outer=TRUE, at=.26)
  
  mtext(side = 1, "Longitude", line = 1, cex=2, outer=TRUE, at=.25)
  mtext(side = 1, "Longitude", line = 1, cex=2, outer=TRUE, at=.75)
  
  mtext(side = 3, "Constituency Level", line = 1, cex=2.5, outer=TRUE, at=.25)
  mtext(side = 3, "County Level", line = 1, cex=2.5, outer=TRUE, at=.75)
  
  dev.off()
  
  # table ----
  riskOut$interceptSummary
  riskOut$fixedEffectSummary
  riskOut$parameterSummaryTable
  
  tab = riskOut$fixedEffectSummary[,c(1, 2, 4, 3, 5)]
  riskOut$parameterSummaryTable = data.frame(riskOut$parameterSummaryTable)
  names(tab) = names(riskOut$parameterSummaryTable)
  tab = rbind(tab, riskOut$parameterSummaryTable)
  row.names(tab)[1:2] = c("Intercept", "Urban Effect")
  
  xtable(tab, 
         caption=paste0("Summary statistics for risk model parameters when ", 
                        "fit to NMR data from the 2014 KDHS for 2010-2014."), 
         label="tab:riskModTable", digits=3)
  
  # old plots ----
  
  # plot credible interval widths
  pixelWidth = prevalenceCIWidthPixel
  constituencyWidth = apply(subareaPop$pFineScalePrevalence, 1, getWidth)
  countyWidth = apply(areaPop$pFineScalePrevalence, 1, getWidth)
  # provinceWidth = apply(agg$aggregatedResultsLCPB$regionMatrices$p, 1, getWidth)
  widthRange = range(c(rangePrevalenceCIWidthConstituency, 
                       rangePrevalenceCIWidthCounty))
  widthRangePixel = rangePrevalenceCIWidthPixel
  
  png(paste0(figDirectory, "application/prevalenceCIWidth", logisticText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelWidth, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(widthRangePixel), nx=160, ny=160, main="", cex.main=3, col=sdCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
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
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyWidth, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  # plotMapDat(plotVar=provinceWidth, new = TRUE, mapDat=regionMap, 
  #            main="", scaleFun=log, scaleFunInverse=exp, 
  #            cols=sdCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
  #            xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
  #            legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
  #            plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
  #            ylab="")
  # # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Posterior Prevalence 95% CI Width", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  ## 2 x 2 plot of predictions (counts)
  
  # plot means
  pixelMean = rowMeans(pixelPop$ZFineScaleRisk)
  constituencyMean = rowMeans(subareaPop$ZFineScalePrevalence)
  countyMean = rowMeans(areaPop$ZFineScalePrevalence)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPB$regionMatrices$Z)
  meanRangePixel = range(pixelMean[pixelMean >= .05])
  meanRangeConstituency = rangeCountPredConstituency
  meanRangeCounty = rangeCountPredCounty
  # meanRangeProvince = rangeCountPredProvince
  
  png(paste0(figDirectory, "application/countMean", logisticText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(meanRangePixel), nx=160, ny=160, main="", cex.main=3, col=meanCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  meanTicksPixel = c(1, 10, 50, 100, pretty(meanRangePixel, n=5)[-c(1, 4, 6)])
  meanTickLabelsPixel = as.character(meanTicksPixel)
  meanTicksConstituency = c(100, 500, pretty(meanRangeConstituency, n=8)[-1])
  meanTickLabelsConstituency = as.character(meanTicksConstituency)
  meanTicksCounty = c(500, 1000, pretty(meanRangeCounty, n=8)[-c(1, 6)])
  meanTickLabelsCounty = as.character(meanTicksCounty)
  # meanTicksProvince = pretty(meanRangeProvince, n=8)[-1]
  # meanTickLabelsProvince = as.character(meanTicksProvince)
  image.plot(zlim=range(log(meanRangePixel)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
             col=meanCols, add = TRUE, axis.args=list(at=log(meanTicksPixel), labels=meanTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeConstituency), ticks=meanTicksConstituency, tickLabels=meanTickLabelsConstituency, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=meanCols, zlim=log(meanRangeCounty), ticks=meanTicksCounty, tickLabels=meanTickLabelsCounty, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # # province level
  # plotMapDat(plotVar=provinceMean, new = TRUE, mapDat=regionMap, 
  #            main="", scaleFun=log, scaleFunInverse=exp, 
  #            cols=meanCols, zlim=log(meanRangeProvince), ticks=meanTicksProvince, tickLabels=meanTickLabelsProvince, 
  #            xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
  #            legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), 
  #            plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
  #            ylab="")
  # # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Posterior Mean Burden", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  # plot credible interval means
  pixelWidth = countCIWidthPixel
  constituencyWidth = apply(subareaPop$ZFineScalePrevalence, 1, getWidth)
  countyWidth = apply(areaPop$ZFineScalePrevalence, 1, getWidth)
  # provinceWidth = apply(agg$aggregatedResultsLCPB$regionMatrices$Z, 1, getWidth)
  widthRangePixel = rangeCountCIWidthPixel # range(pixelWidth)
  widthRangeConstituency = rangeCountCIWidthConstituency
  widthRangeCounty = rangeCountCIWidthCounty
  # widthRangeProvince = rangeCountCIWidthProvince
  
  png(paste0(figDirectory, "application/countWidth", logisticText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popGrid$lon, popGrid$lat, pixelWidth, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
             zlim=log(widthRangePixel), nx=160, ny=160, main="", cex.main=3, col=sdCols, 
             add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
             xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  widthTicksPixel = c(1, 10, 50, 100, pretty(widthRangePixel, n=5)[-c(1)])
  widthTickLabelsPixel = as.character(widthTicksPixel)
  widthTicksConstituency = c(100, pretty(widthRangeConstituency, n=8)[-1])
  widthTickLabelsConstituency = as.character(widthTicksConstituency)
  widthTicksCounty = c(500, pretty(widthRangeCounty, n=5)[-c(1)])
  widthTickLabelsCounty = as.character(widthTicksCounty)
  # widthTicksProvince = pretty(widthRangeProvince, n=8)
  # widthTickLabelsProvince = as.character(widthTicksProvince)
  image.plot(zlim=range(log(widthRangePixel)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
             col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksPixel), labels=widthTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
             legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  # constituency level
  plotMapDat(plotVar=constituencyWidth, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeConstituency), ticks=widthTicksConstituency, tickLabels=widthTickLabelsConstituency, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
             xlab="", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # county level
  plotMapDat(plotVar=countyWidth, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeCounty), ticks=widthTicksCounty, tickLabels=widthTickLabelsCounty, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  # plotMapDat(plotVar=provinceWidth, new = TRUE, mapDat=regionMap, 
  #            main="", scaleFun=log, scaleFunInverse=exp, 
  #            cols=sdCols, zlim=log(widthRangeProvince), ticks=widthTicksProvince, tickLabels=widthTickLabelsProvince, 
  #            xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
  #            legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
  #            plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
  #            ylab="")
  # # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the top margin
  mtext(side = 3, "Total deaths 95% CI width", line = 0.5, cex=2.5, outer=TRUE)
  
  dev.off()
  
  ## 2 x 3 plot of predictions (relative prevalence)
  
  # plot mean
  constituencyMean = relativePrevalencePredConstituency
  countyMean = relativePrevalencePredCounty
  # provinceMean = relativePrevalencePredProvince
  meanRange = range(c(rangeRelativePrevalencePredConstituency, 
                      rangeRelativePrevalencePredCounty))
  urbDivergingCols = makeGreenBlueDivergingColors(64, center=1, valRange=meanRange)
  
  png(paste0(figDirectory, "application/relativePrevalence", logisticText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,4,4,7), mar=c(6.1, 6.5, 1.1, 2.5))
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=urbDivergingCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.2,1)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2.5), legend.mar=0, lwd=.25, border=rgb(.4,.4,.4), 
             xlab="", ylab="Latitude")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the left margin
  mtext(side = 2, "Rel. Prevalence Posterior Mean", line = 6, cex=2.5)
  
  # county level
  plotMapDat(plotVar=countyMean, new = TRUE, 
             main="", #scaleFun=logit, scaleFunInverse=expit, 
             cols=urbDivergingCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, ylab="", xlab="", 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2.5), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  # plotMapDat(plotVar=provinceMean, new = TRUE, mapDat=regionMap, 
  #            main="", #scaleFun=logit, scaleFunInverse=expit, 
  #            cols=urbDivergingCols, zlim=meanRange, # ticks=meanTicks, tickLabels=meanTickLabels, 
  #            xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
  #            legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
  #            plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
  #            ylab="", xlab="")
  # # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # now plot credible interval widths
  constituencyWidth = relativePrevalenceCIWidthConstituency
  countyWidth = relativePrevalenceCIWidthCounty
  # provinceWidth = relativePrevalenceCIWidthProvince
  widthRangeConstituency = range(rangeRelativePrevalenceCIWidthConstituency, na.rm=TRUE)
  widthRangeCounty = range(rangeRelativePrevalenceCIWidthCounty, na.rm=TRUE)
  # widthRangeProvince = range(rangeRelativePrevalenceCIWidthProvince, na.rm=TRUE)
  
  widthTicksConstituency = pretty(widthRangeConstituency, n=5)[-1]
  widthTickLabelsConstituency = as.character(widthTicksConstituency)
  widthTicksCounty = pretty(widthRangeCounty, n=5)[-c(1)]
  widthTickLabelsCounty = as.character(widthTicksCounty)
  # widthTicksProvince = pretty(widthRangeProvince, n=5)
  # widthTickLabelsProvince = as.character(widthTicksProvince)
  
  # constituency level
  plotMapDat(plotVar=constituencyWidth, mapDat=constituencyMap, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeConstituency), ticks=widthTicksConstituency, tickLabels=widthTickLabelsConstituency, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2.5), legend.mar=0, lwd=.25, border=rgb(.4,.4,.4), 
             xlab="Longitude", ylab="Latitude")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # add title in the left margin
  mtext(side = 2, "Rel. Prevalence 95% CI Width", line = 6, cex=2.5)
  
  # county level
  plotMapDat(plotVar=countyWidth, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=sdCols, zlim=log(widthRangeCounty), ticks=widthTicksCounty, tickLabels=widthTickLabelsCounty, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, ylab="", xlab="Longitude", 
             legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
             plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2.5), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  # province level
  # plotMapDat(plotVar=provinceWidth, new = TRUE, mapDat=regionMap, 
  #            main="", scaleFun=log, scaleFunInverse=exp, 
  #            cols=sdCols, zlim=log(widthRangeProvince), ticks=widthTicksProvince, tickLabels=widthTickLabelsProvince, 
  #            xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
  #            legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
  #            plotArgs=list(cex.main=3, cex.axis=2, cex.lab=3), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4), 
  #            ylab="")
  # # plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  # plotMapDat(mapDat=regionMap, lwd=2.5)
  
  dev.off()
  
  ##### Box plots of credible interval widths and relative widths to the stratified model
  # modelNames: S, SC, SCP (stratified, cluster, population)
  
  ## prevalence uncertainty
  
  # absolute prevalence uncertainty
  pdf(paste0(figDirectory, "application/prevalenceWidthBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  pixels = length(prevalenceCIWidthPixel)
  pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
                       width=c(prevalenceCIWidthPixellcpb, prevalenceCIWidthPixelLCPb, prevalenceCIWidthPixel))
  boxplot(width~modelName, data=pixelDat, ylab="", xlab="Model", main="Pixel", col="skyblue", log="y")
  mtext(side = 2, "95% CI width", line = 3, cex=1)
  
  constituencies = length(prevalenceCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                        width=c(prevalenceCIWidthConstituencylcpb, prevalenceCIWidthConstituencyLCPb, prevalenceCIWidthConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="95% CI Width", xlab="Model", main="Constituency", col="skyblue", log="y")
  
  counties = length(prevalenceCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                           width=c(prevalenceCIWidthCountylcpb, prevalenceCIWidthCountyLCPb, prevalenceCIWidthCounty))
  boxplot(width~modelName, data=countyDat, ylab="95% CI Width", xlab="Model", main="County", col="skyblue", log="y")
  mtext(side = 2, "95% CI width", line = 3, cex=1)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(prevalenceCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                              width=c(prevalenceCIWidthProvincelcpb, prevalenceCIWidthProvinceLCPb, prevalenceCIWidthProvince))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue", log="y")
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Prevalence 95% CI Width", line = 0.5, cex=1, outer=TRUE)
  dev.off()
  
  # relative prevalence uncertainty
  pdf(paste0(figDirectory, "application/prevalenceRelWidthBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  pixels = length(prevalenceCIWidthPixel)
  pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
                        width=c(100*(prevalenceCIWidthPixellcpb-prevalenceCIWidthPixellcpb)/prevalenceCIWidthPixellcpb, 
                                100*(prevalenceCIWidthPixelLCPb-prevalenceCIWidthPixellcpb)/prevalenceCIWidthPixellcpb, 
                                100*(prevalenceCIWidthPixel-prevalenceCIWidthPixellcpb)/prevalenceCIWidthPixellcpb))
  boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
          main="Pixel", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  constituencies = length(prevalenceCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                        width=c(100*(prevalenceCIWidthConstituencylcpb-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb, 
                                100*(prevalenceCIWidthConstituencyLCPb-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb, 
                                100*(prevalenceCIWidthConstituency-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  
  counties = length(prevalenceCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                               width=c(100*(prevalenceCIWidthCountylcpb-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb, 
                                       100*(prevalenceCIWidthCountyLCPb-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb, 
                                       100*(prevalenceCIWidthCounty-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  # provinces = length(prevalenceCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                        width=c(100*(prevalenceCIWidthProvincelcpb-prevalenceCIWidthProvincelcpb)/prevalenceCIWidthProvincelcpb, 
  #                                100*(prevalenceCIWidthProvinceLCPb-prevalenceCIWidthProvincelcpb)/prevalenceCIWidthProvincelcpb, 
  #                                100*(prevalenceCIWidthProvince-prevalenceCIWidthProvincelcpb)/prevalenceCIWidthProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Prevalence 95% CI Relative Width", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  browser()
  # relative prevalence uncertainty
  pdf(paste0(figDirectory, "application/prevalenceRelSDBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  pixels = length(prevalenceSDPixel)
  pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
                        width=c(100*(prevalenceSDPixellcpb-prevalenceSDPixellcpb)/prevalenceSDPixellcpb, 
                                100*(prevalenceSDPixelLCPb-prevalenceSDPixellcpb)/prevalenceSDPixellcpb, 
                                100*(prevalenceSDPixel-prevalenceSDPixellcpb)/prevalenceSDPixellcpb))
  boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
          main="Pixel", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  constituencies = length(prevalenceSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(prevalenceSDConstituencylcpb-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb, 
                                       100*(prevalenceSDConstituencyLCPb-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb, 
                                       100*(prevalenceSDConstituency-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  
  counties = length(prevalenceSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(prevalenceSDCountylcpb-prevalenceSDCountylcpb)/prevalenceSDCountylcpb, 
                                 100*(prevalenceSDCountyLCPb-prevalenceSDCountylcpb)/prevalenceSDCountylcpb, 
                                 100*(prevalenceSDCounty-prevalenceSDCountylcpb)/prevalenceSDCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  # provinces = length(prevalenceSDProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(prevalenceSDProvincelcpb-prevalenceSDProvincelcpb)/prevalenceSDProvincelcpb, 
  #                                  100*(prevalenceSDProvinceLCPb-prevalenceSDProvincelcpb)/prevalenceSDProvincelcpb, 
  #                                  100*(prevalenceSDProvince-prevalenceSDProvincelcpb)/prevalenceSDProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Prevalence Relative SD", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  
  ## count uncertainty
  
  # absolute count uncertainty
  pdf(paste0(figDirectory, "application/countWidthBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  pixels = length(countCIWidthPixel)
  pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
                        width=c(countCIWidthPixellcpb, countCIWidthPixelLCPb, countCIWidthPixel))
  boxplot(width~modelName, data=pixelDat, ylab="", xlab="Model", main="Pixel", col="skyblue", log="y")
  mtext(side = 2, "95% CI width", line = 3, cex=1)
  
  constituencies = length(countCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(countCIWidthConstituencylcpb, countCIWidthConstituencyLCPb, countCIWidthConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="95% CI Width", xlab="Model", main="Constituency", col="skyblue", log="y")
  
  counties = length(countCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(countCIWidthCountylcpb, countCIWidthCountyLCPb, countCIWidthCounty))
  boxplot(width~modelName, data=countyDat, ylab="95% CI Width", xlab="Model", main="County", col="skyblue", log="y")
  mtext(side = 2, "95% CI width", line = 3, cex=1)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(countCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(countCIWidthProvincelcpb, countCIWidthProvinceLCPb, countCIWidthProvince))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue", log="y")
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Total Deaths 95% CI Width", line = 0.5, cex=1, outer=TRUE)
  dev.off()
  
  # relative count uncertainty
  pdf(paste0(figDirectory, "application/countRelWidthBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  pixels = length(countCIWidthPixel)
  pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
                        width=c(100*(countCIWidthPixellcpb-countCIWidthPixellcpb)/countCIWidthPixellcpb, 
                                100*(countCIWidthPixelLCPb-countCIWidthPixellcpb)/countCIWidthPixellcpb, 
                                100*(countCIWidthPixel-countCIWidthPixellcpb)/countCIWidthPixellcpb))
  boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
          main="Pixel", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  constituencies = length(countCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(countCIWidthConstituencylcpb-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb, 
                                       100*(countCIWidthConstituencyLCPb-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb, 
                                       100*(countCIWidthConstituency-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  
  counties = length(countCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(countCIWidthCountylcpb-countCIWidthCountylcpb)/countCIWidthCountylcpb, 
                                 100*(countCIWidthCountyLCPb-countCIWidthCountylcpb)/countCIWidthCountylcpb, 
                                 100*(countCIWidthCounty-countCIWidthCountylcpb)/countCIWidthCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  # provinces = length(countCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(countCIWidthProvincelcpb-countCIWidthProvincelcpb)/countCIWidthProvincelcpb, 
  #                                  100*(countCIWidthProvinceLCPb-countCIWidthProvincelcpb)/countCIWidthProvincelcpb, 
  #                                  100*(countCIWidthProvince-countCIWidthProvincelcpb)/countCIWidthProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Total Deaths 95% CI Relative Width", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  
  pdf(paste0(figDirectory, "application/countRelSDBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  pixels = length(countSDPixel)
  pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
                        width=c(100*(countSDPixellcpb-countSDPixellcpb)/countSDPixellcpb, 
                                100*(countSDPixelLCPb-countSDPixellcpb)/countSDPixellcpb, 
                                100*(countSDPixel-countSDPixellcpb)/countSDPixellcpb))
  boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
          main="Pixel", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  constituencies = length(countSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(countSDConstituencylcpb-countSDConstituencylcpb)/countSDConstituencylcpb, 
                                       100*(countSDConstituencyLCPb-countSDConstituencylcpb)/countSDConstituencylcpb, 
                                       100*(countSDConstituency-countSDConstituencylcpb)/countSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  
  counties = length(countSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(countSDCountylcpb-countSDCountylcpb)/countSDCountylcpb, 
                                 100*(countSDCountyLCPb-countSDCountylcpb)/countSDCountylcpb, 
                                 100*(countSDCounty-countSDCountylcpb)/countSDCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  # provinces = length(countSDProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(countSDProvincelcpb-countSDProvincelcpb)/countSDProvincelcpb, 
  #                                  100*(countSDProvinceLCPb-countSDProvincelcpb)/countSDProvincelcpb, 
  #                                  100*(countSDProvince-countSDProvincelcpb)/countSDProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Total Deaths Relative SD", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  
  ## relative prevalence uncertainty
  
  # absolute relative prevalence uncertainty
  pdf(paste0(figDirectory, "application/relativePrevalenceWidthBoxplot", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  constituencies = length(relativePrevalenceCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(relativePrevalenceCIWidthConstituencylcpb, relativePrevalenceCIWidthConstituencyLCPb, relativePrevalenceCIWidthConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="95% CI Width", xlab="Model", main="Constituency", col="skyblue")
  mtext(side = 2, "95% CI width", line = 3, cex=1)
  
  counties = length(relativePrevalenceCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(relativePrevalenceCIWidthCountylcpb, relativePrevalenceCIWidthCountyLCPb, relativePrevalenceCIWidthCounty))
  boxplot(width~modelName, data=countyDat, ylab="95% CI Width", xlab="Model", main="County", col="skyblue")
  
  # provinces = length(relativePrevalenceCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
                           # width=c(relativePrevalenceCIWidthProvincelcpb, relativePrevalenceCIWidthProvinceLCPb, relativePrevalenceCIWidthProvince))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  
  # relative relativePrevalence uncertainty
  constituencies = length(relativePrevalenceCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(relativePrevalenceCIWidthConstituencylcpb-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb, 
                                       100*(relativePrevalenceCIWidthConstituencyLCPb-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb, 
                                       100*(relativePrevalenceCIWidthConstituency-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  counties = length(relativePrevalenceCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(relativePrevalenceCIWidthCountylcpb-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb, 
                                 100*(relativePrevalenceCIWidthCountyLCPb-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb, 
                                 100*(relativePrevalenceCIWidthCounty-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(relativePrevalenceCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(relativePrevalenceCIWidthProvincelcpb-relativePrevalenceCIWidthProvincelcpb)/relativePrevalenceCIWidthProvincelcpb, 
  #                                  100*(relativePrevalenceCIWidthProvinceLCPb-relativePrevalenceCIWidthProvincelcpb)/relativePrevalenceCIWidthProvincelcpb, 
  #                                  100*(relativePrevalenceCIWidthProvince-relativePrevalenceCIWidthProvincelcpb)/relativePrevalenceCIWidthProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the top margin
  mtext(side = 3, "Relative Prevalence", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  
  # absolute relative prevalence uncertainty
  pdf(paste0(figDirectory, "application/relativePrevalenceWidthBoxplotSD", logisticText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  onstituencies = length(relativePrevalenceSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(relativePrevalenceSDConstituencylcpb, relativePrevalenceSDConstituencyLCPb, relativePrevalenceSDConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="SD", xlab="Model", main="Constituency", col="skyblue")
  mtext(side = 2, "SD", line = 3, cex=1)
  
  counties = length(relativePrevalenceSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(relativePrevalenceSDCountylcpb, relativePrevalenceSDCountyLCPb, relativePrevalenceSDCounty))
  boxplot(width~modelName, data=countyDat, ylab="SD", xlab="Model", main="County", col="skyblue")
  
  # provinces = length(relativePrevalenceSDProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(relativePrevalenceSDProvincelcpb, relativePrevalenceSDProvinceLCPb, relativePrevalenceSDProvince))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  
  # relative relativePrevalence uncertainty
  constituencies = length(relativePrevalenceSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(relativePrevalenceSDConstituencylcpb-relativePrevalenceSDConstituencylcpb)/relativePrevalenceSDConstituencylcpb, 
                                       100*(relativePrevalenceSDConstituencyLCPb-relativePrevalenceSDConstituencylcpb)/relativePrevalenceSDConstituencylcpb, 
                                       100*(relativePrevalenceSDConstituency-relativePrevalenceSDConstituencylcpb)/relativePrevalenceSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  counties = length(relativePrevalenceSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(relativePrevalenceSDCountylcpb-relativePrevalenceSDCountylcpb)/relativePrevalenceSDCountylcpb, 
                                 100*(relativePrevalenceSDCountyLCPb-relativePrevalenceSDCountylcpb)/relativePrevalenceSDCountylcpb, 
                                 100*(relativePrevalenceSDCounty-relativePrevalenceSDCountylcpb)/relativePrevalenceSDCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(relativePrevalenceSDProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(relativePrevalenceSDProvincelcpb-relativePrevalenceSDProvincelcpb)/relativePrevalenceSDProvincelcpb, 
  #                                  100*(relativePrevalenceSDProvinceLCPb-relativePrevalenceSDProvincelcpb)/relativePrevalenceSDProvincelcpb, 
  #                                  100*(relativePrevalenceSDProvince-relativePrevalenceSDProvincelcpb)/relativePrevalenceSDProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the top margin
  mtext(side = 3, "Relative Prevalence", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
}

# for now, just gets the results for Nairobi
getMortGridResolutionResults = function() {
  seed=123
  set.seed(seed)
  
  # make sampling frame
  easpaN = makeDefaultEASPA()
  easpaN = easpaN[easpaN$area == "Nairobi", ]
  poppconN = poppcon
  poppconN = poppconN[poppconN$County == "Nairobi", ]
  
  # make population density grids
  resolutions = c(.1, .2, .5, 1, 2, 5, 10, 20)
  deltas = c(.005, .01, .025, .05, .1, .1, .2, .4)
  meanNeighbors = c(rep(500, 5), rep(50, 3))
  popGrids = list()
  popGridsAdjusted = list()
  for(i in 1:length(resolutions)) {
    print(paste0("Creating integration grid at ", resolutions[i], " km resolution"))
    thisPopGrid = makeInterpPopGrid(kmRes=resolutions[i], 
                                    mean.neighbor=meanNeighbors[i], 
                                    delta=deltas[i], conMap=adm2, poppcon=poppcon, 
                                    mapDat=adm1, polygonSubsetI=30) # Nairobi is the 30th one
    thisPopGrid$area = thisPopGrid$admin1
    thisPopGrid$constituency = thisPopGrid$admin2
    thisPopGridAdjusted = adjustPopGrid(thisPopGrid, poppcon, "Constituency")
    
    # subset grids to be over only Nairobi
    thisPopGrid = thisPopGrid[thisPopGrid$admin1 == "Nairobi",]
    thisPopGridAdjusted = thisPopGridAdjusted[thisPopGridAdjusted$admin1 == "Nairobi",]
    
    popGrids = c(popGrids, list(thisPopGrid))
    popGridsAdjusted = c(popGridsAdjusted, list(thisPopGridAdjusted))
  }
  sortI = sort(resolutions, index.return=TRUE)$ix
  allResolutions = resolutions[sortI]
  popGrids = popGrids[sortI]
  popGridsAdjusted = popGridsAdjusted[sortI]
  names(popGrids) = paste0("popGrid", allResolutions)
  names(popGridsAdjusted) = paste0("popGridAdjusted", allResolutions)
  ns = sapply(popGrids, nrow)
  sum(sapply(popGrids, nrow))
  endIs = cumsum(ns)
  startIs = c(1, endIs[-length(endIs)]+1)
  
  # combine the grids
  popMatCombined = do.call("rbind", popGrids)
  
  # fit SPDE cluster level risk model over all integration points
  spdeFitN = fitSPDEKenyaDat(mort, nPostSamples=10000, popMat=popMatCombined)
  
  # apply aggregation models at each resolution
  nSamples = c(c(500, 1000, 5000), rep(10000, length(resolutions)-3))
  aggResultsN = list()
  for(i in 1:length(popGrids)) {
    thisNSamples = nSamples[i]
    
    # obtain the grids at this resolution
    thisPopMat = popGrids[[i]]
    thisPopMatAdjusted = popGridsAdjusted[[i]]
    
    # obtain model output at this resolution
    thisResolutionI = startIs[i]:endIs[i]
    thisUDraws = spdeFitN$uDraws[thisResolutionI,1:thisNSamples]
    sigmaEpsilonDraws = spdeFitN$sigmaEpsilonDraws[1:thisNSamples]
    print(paste0("Generating aggregation results for grid resolution ", resolutions[i]))
    
    thisAggResultsN = modLCPB(thisUDraws, sigmaEpsilonDraws, easpaN, thisPopMat, 
                              thisPopMatAdjusted, doLCPb=TRUE, doIHMEModel=TRUE, 
                              constituencyPop=poppconN, ensureAtLeast1PerConstituency=TRUE, 
                              logisticApproximation=FALSE, verbose=TRUE, 
                              stopOnFrameMismatch=FALSE)
    
    aggResultsN = c(aggResultsN, list(thisAggResultsN))
  }
  names(aggResultsN) = paste0("aggResultsN", resolutions)
  save(popGrids, popGridsAdjusted, aggResultsN, file="savedOutput/application/gridResolutionTestNairobi.RData")
  
  invisible(NULL)
}

makeMortGridResolutionPlots = function() {
  out = load("savedOutput/application/gridResolutionTestNairobi.RData")
  
  # Calculate predictions, variance, 95% Coverage interval width
  predsConstituency = matrix(nrow=nrow(aggResultsN[[1]]$aggregatedResultsLCPB$constituencyMatrices$p), ncol=length(resolutions))
  residsConstituencySmoothRisk = list()
  residsConstituencyRisk = list()
  residsConstituencyPrevalence = list()
  residsConstituencyGriddedRisk = list()
  for(i in 1:length(resolutions)) {
    thesePreds = rowMeans(aggResultsN[[i]]$aggregatedResultslcpb$constituencyMatrices$p)
    theseResidsSmoothRisk = sweep(aggResultsN[[i]]$aggregatedResultslcpb$constituencyMatrices$p, 1, thesePreds, "-")
    theseResidsRisk = sweep(aggResultsN[[i]]$aggregatedResultsLCPb$constituencyMatrices$p, 1, thesePreds, "-")
    theseResidsPrevalence = sweep(aggResultsN[[i]]$aggregatedResultsLCPB$constituencyMatrices$p, 1, thesePreds, "-")
    theseResidsGriddedRisk = sweep(aggResultsN[[i]]$aggregatedResultsIHME$constituencyMatrices$p, 1, thesePreds, "-")
    predsConstituency[,i] = thesePreds
    residsConstituencySmoothRisk = c(residsConstituencySmoothRisk, list(theseResidsSmoothRisk))
    residsConstituencyRisk = c(residsConstituencyRisk, list(theseResidsRisk))
    residsConstituencyPrevalence = c(residsConstituencyPrevalence, list(theseResidsPrevalence))
    residsConstituencyGriddedRisk = c(residsConstituencyGriddedRisk, list(theseResidsGriddedRisk))
  }
  
  lowConstituencySmoothRisk = sapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
  lowConstituencyRisk = sapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
  lowConstituencyPrevalence = sapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
  lowConstituencyGriddedRisk = sapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
  highConstituencySmoothRisk = sapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
  highConstituencyRisk = sapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
  highConstituencyPrevalence = sapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
  highConstituencyGriddedRisk = sapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
  
  CIWidthSmoothRisk = highConstituencySmoothRisk - lowConstituencySmoothRisk
  CIWidthRisk = highConstituencyRisk - lowConstituencyRisk
  CIWidthPrevalence = highConstituencyPrevalence - lowConstituencyPrevalence
  CIWidthGriddedRisk = highConstituencyGriddedRisk - lowConstituencyGriddedRisk
  
  meanCIWidthSmoothRisk = colMeans(CIWidthSmoothRisk)
  meanCIWidthRisk = colMeans(CIWidthRisk)
  meanCIWidthPrevalence = colMeans(CIWidthPrevalence)
  meanCIWidthGriddedRisk = colMeans(CIWidthGriddedRisk)
  
  # Plot central predictions versus resolution
  ylim = range(c(predsConstituency))
  pdf("figures/application/gridResolutionTestPredictionVRes.pdf", width=5, height=5)
  cols = rainbow(4)
  thisFrame = data.frame(predsConstituency)
  boxplot(predsConstituency, names=resolutions, col="skyblue", 
          main="", xlab="Resolution", ylab="Prediction")
  dev.off()
  
  # Plot CI Widths versus resolution and model
  constituenciesN = poppcon$Constituency[poppcon$County=="Nairobi"]
  CIWidth = c(c(CIWidthSmoothRisk), c(CIWidthRisk), c(CIWidthPrevalence), c(CIWidthGriddedRisk))
  tempRes = resolutions[col(CIWidthSmoothRisk)]
  tempCon = factor(as.character(poppcon$Constituency[constituenciesN][col(CIWidthSmoothRisk)]))
  N=length(tempCon)
  nIntegrationPoints = sapply(popGrids, function(x) {table(x$admin2)})
  percentIncreaseWidthSmoothRisk = sweep(CIWidthSmoothRisk, 1, CIWidthSmoothRisk[,1], function(x,y){100*x/y - 100})
  CIWidthFrame = data.frame(Constituency=rep(tempCon, 4), 
                            Resolution=rep(tempRes, 4), 
                            Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                           rep("Prevalence", N), rep("Gridded risk", N)), 
                                         levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                            CIWidth=CIWidth, 
                            nIntegrationPoints=rep(c(nIntegrationPoints), 4))
  
  pdf("figures/application/gridResolutionTestCIWidthVRes.pdf", width=7, height=5)
  ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
    geom_boxplot(position="dodge2") + ylim(0, max(CIWidth)) + 
    labs(x="Grid resolution (km)", y="95% credible interval width", fill="Model") + 
    theme_classic()
  dev.off()
  
  pdf("figures/application/gridResolutionTestCIWidthVResLog.pdf", width=7, height=5)
  ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
    geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
    labs(x="Grid resolution (km)", y="95% credible interval width", fill="Model") + 
    theme_classic()
  dev.off()
  
  # set color scales
  blueGreenCols = rev(makeGreenBlueSequentialColors(64))
  yellowBlueCols = makeBlueGreenYellowSequentialColors(64)
  
  # plot uDraws
  ns = sapply(popGrids, nrow)
  sum(sapply(popGrids, nrow))
  endIs = cumsum(ns)
  startIs = c(1, endIs[-length(endIs)]+1)
  nSamples = c(c(500, 1000), rep(10000, length(resolutions)-1))
  cexs = c(.005, .04, .1, .2, .7, 1.2, 1.5, 2)
  zlim=range(colMeans(spdeFitN$uDraws))
  xlim=range(popGrids[[1]]$lon)
  ylim=range(popGrids[[1]]$lat)
  for(i in 1:length(resolutions)) {
    thisNSamples = nSamples[i]
    
    # obtain the grids at this resolution
    thisPopMat = popGrids[[i]]
    
    # obtain model output at this resolution
    thisResolutionI = startIs[i]:endIs[i]
    thisUDraws = spdeFitN$uDraws[thisResolutionI,1:thisNSamples]
    
    # make the uDraws plot
    pdf(paste0("figures/application/gridResolutionTestUDraws", i, ".pdf"), width=5, height=5)
    plotWithColor(thisPopMat$lon, thisPopMat$lat, rowMeans(thisUDraws), 
                  colScale=blueGreenCols, xlab="Longitude", 
                  ylab="Latitude", main=paste0("uDraws resolution ", resolutions[i], " km"), 
                  cex=cexs[i], pch=19, zlim=zlim, xlim=xlim, ylim=ylim, 
                  ordering="decreasing")
    plotMapDat(mapDat=adm2compressed, new=FALSE)
    points(mort$lon, mort$lat, cex=.5)
    dev.off()
  }
  
  invisible(NULL)
}





