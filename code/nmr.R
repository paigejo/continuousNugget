# script for getting all NMR results for the application
# 11GB for 10k posterior samples, 1.1Gb for 10k posterior samples under coarse resolution
getMortResults = function(seed=123, useCoarseGrid=FALSE, logisticApproximation=FALSE, nPostSamples=1000, 
                          resultType=c("std", "FBpop", "census2019", "censusJittered"), doGridLevel=resultType=="std") {
  resultType = match.arg(resultType)
  set.seed(seed)
  
  time1 = proc.time()[3]
  
  if(resultType == "std") {
    easpa = makeDefaultEASPA()
    poppsub = poppsubKenyaThresh
    
    if(useCoarseGrid) {
      popMat = popGridCoarseThresh
      popMatAdjusted = popGridCoarseAdjustedThresh
    } else {
      popMat = popGridThresh
      popMatAdjusted = popGridAdjustedThresh
    }
  } else if(resultType == "FBpop") {
    out = load("savedOutput/global/kenyaFacePopulationMats.RData")
    out = load("savedOutput/global/poppsubFace.RData")
    
    easpa = makeEASPAfacebook()
    poppsub = poppsubKenyaFaceThresh
    
    if(useCoarseGrid) {
      popMat = popMatCoarseFaceThresh
      popMatAdjusted = popMatCoarseAdjustedFaceThresh
    } else {
      popMat = popMatKenyaFaceThresh
      popMatAdjusted = popMatKenyaFaceNeonatalThresh
    }
  } else if(resultType == "census2019") {
    load("savedOutput/global/kenya2019PopulationMats.RData")
    load("savedOutput/global/poppsub2019.RData")
    
    easpa = makeEASPA2019()
    poppsub = poppsubKenya2019Thresh
    
    if(useCoarseGrid) {
      popMat = popMatCoarse2019Thresh
      popMatAdjusted = popMatCoarseAdjusted2019Thresh
    } else {
      popMat = popMatKenya2019Thresh
      popMatAdjusted = popMatKenya2019NeonatalThresh
    }
  } else if(resultType == "censusJittered") {
    load("savedOutput/global/kenyaJitteredPopulationMats.RData")
    load("savedOutput/global/poppsubJittered.RData")
    
    easpa = makeEASPAJittered()
    poppsub = poppsubKenyaJitteredThresh
    
    if(useCoarseGrid) {
      popMat = popMatCoarseJitteredThresh
      popMatAdjusted = popMatCoarseAdjustedJitteredThresh
    } else {
      popMat = popMatKenyaJitteredThresh
      popMatAdjusted = popMatKenyaJitteredNeonatalThresh
    }
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
                     stratifyByUrban=TRUE, subareaLevel=TRUE, gridLevel=doGridLevel, 
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
  fileName = paste0("savedOutput/application/finalMort", ifelse(useCoarseGrid, "Coarse", ""), "_", resultType, ".RData")
  save(riskOut, pixelPop, subareaPop, areaPop, logitDraws, aggregationTimings, rawTimes, totalTimes, resultType, file=fileName)
  
  list(riskOut=riskOut, pixelPop=pixelPop, subareaPop=subareaPop, 
       areaPop=areaPop, logitDraws=logitDraws, 
       aggregationTimings=aggregationTimings, 
       rawTimes=rawTimes, totalTimes=totalTimes, resultType=resultType)
}

# Make plots for the neonatal mortality application
makeMortPlots = function(logisticApproximation=FALSE, coarse=TRUE, signif=.95, 
                         resultType=c("std", "FBpop", "census2019", "censusJittered")) {
  alpha = 1 - signif
  
  # first load the model predictions
  if(logisticApproximation) {
    stop("logisticApproximation not currently used")
  }
  logisticText = ifelse(!logisticApproximation, "", "logisticApprox")
  coarseText = ifelse(!coarse, "", "Coarse")
  resultTypeText = paste0("_", resultType)
  out = load(paste0("savedOutput/application/finalMort", coarseText, logisticText, resultTypeText, ".RData"))
  
  if(resultType == "std") {
    easpa = makeDefaultEASPA()
    poppsub = poppsubKenyaThresh
    
    if(coarse) {
      popMat = popGridCoarseThresh
      popMatAdjusted = popGridCoarseAdjustedThresh
    } else {
      popMat = popGridThresh
      popMatAdjusted = popGridAdjustedThresh
    }
  } else if(resultType == "FBpop") {
    out = load("savedOutput/global/kenyaFacePopulationMats.RData")
    out = load("savedOutput/global/poppsubFace.RData")
    
    easpa = makeEASPAfacebook()
    poppsub = poppsubKenyaFaceThresh
    
    if(coarse) {
      popMat = popMatCoarseFaceThresh
      popMatAdjusted = popMatCoarseAdjustedFaceThresh
    } else {
      popMat = popMatKenyaFaceThresh
      popMatAdjusted = popMatKenyaFaceNeonatalThresh
    }
  } else if(resultType == "census2019") {
    load("savedOutput/global/kenya2019PopulationMats.RData")
    load("savedOutput/global/poppsub2019.RData")
    
    easpa = makeEASPA2019()
    poppsub = poppsubKenya2019Thresh
    
    if(coarse) {
      popMat = popMatCoarse2019Thresh
      popMatAdjusted = popMatCoarseAdjusted2019Thresh
    } else {
      popMat = popMatKenya2019Thresh
      popMatAdjusted = popMatKenya2019NeonatalThresh
    }
  } else if(resultType == "censusJittered") {
    load("savedOutput/global/kenyaJitteredPopulationMats.RData")
    load("savedOutput/global/poppsubJittered.RData")
    
    easpa = makeEASPAJittered()
    poppsub = poppsubKenyaJitteredThresh
    
    if(coarse) {
      popMat = popMatCoarseJitteredThresh
      popMatAdjusted = popMatCoarseAdjustedJitteredThresh
    } else {
      popMat = popMatKenyaJitteredThresh
      popMatAdjusted = popMatKenyaJitteredNeonatalThresh
    }
  }
  
  # calculate the range of predictions and CI widths ----
  rangePrevalencePredPixel = c()
  rangePrevalencePredConStrat = c()
  rangePrevalencePredConstituency = c()
  rangePrevalencePredCounty = c()
  rangePrevalencePredProvince = c()
  rangePrevalenceCIWidthPixel = c()
  rangePrevalenceCIWidthConStrat = c()
  rangePrevalenceCIWidthConstituency = c()
  rangePrevalenceCIWidthCounty = c()
  rangePrevalenceCIWidthProvince = c()
  rangePrevalenceSDPixel = c()
  rangePrevalenceSDConStrat = c()
  rangePrevalenceSDConstituency = c()
  rangePrevalenceSDCounty = c()
  rangePrevalenceSDProvince = c()
  
  rangeCountPredPixel = c()
  rangeCountPredConStrat = c()
  rangeCountPredConstituency = c()
  rangeCountPredCounty = c()
  rangeCountPredProvince = c()
  rangeCountCIWidthPixel = c()
  rangeCountCIWidthConStrat = c()
  rangeCountCIWidthConstituency = c()
  rangeCountCIWidthCounty = c()
  rangeCountCIWidthProvince = c()
  rangeCountSDPixel = c()
  rangeCountSDConStrat = c()
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
  
  areaLevels = c("pixel", "constrat", "constituency", "county") # province results no longer supported
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
    }  else if(thisLevel == "constrat") {
      urbanConstituencies = poppsub$County == "Nairobi" | poppsub$County == "Mombasa"
      undefinedPrevalenceConStrats = c(poppsub$popUrb == 0, poppsub$popRur == 0)
      rangePrevalencePredConStrat = range(c(rangePrevalencePredConStrat, 
                                            rowMeans(subareaPop$pUrbanFineScaleRisk, na.rm=TRUE), 
                                            rowMeans(subareaPop$pRuralFineScaleRisk, na.rm=TRUE)), na.rm=TRUE)
      rangeCountPredConStrat = range(c(rangeCountPredConStrat, 
                                       rowMeans(subareaPop$ZUrbanFineScaleRisk), 
                                       rowMeans(subareaPop$ZRuralFineScaleRisk)), na.rm=TRUE)
      
      prevalenceCIWidthConStrat = c(apply(subareaPop$pUrbanFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}), 
                                    apply(subareaPop$pRuralFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}))
      prevalenceCIWidthConStrat = prevalenceCIWidthConStrat[!undefinedPrevalenceConStrats]
      countCIWidthConStrat = c(apply(subareaPop$ZUrbanFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}), 
                               apply(subareaPop$ZRuralFineScalePrevalence, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}))
      countCIWidthConStrat = countCIWidthConStrat[!undefinedPrevalenceConStrats]
      rangePrevalenceCIWidthConStrat = range(prevalenceCIWidthConStrat)
      rangeCountCIWidthConStrat = range(countCIWidthConStrat)
      
      prevalenceSDConStrat = c(apply(subareaPop$pUrbanFineScalePrevalence, 1, sd, na.rm=TRUE), 
                               apply(subareaPop$pRuralFineScalePrevalence, 1, sd, na.rm=TRUE))
      prevalenceSDConStrat = prevalenceSDConStrat[!undefinedPrevalenceConStrats]
      countSDConStrat = c(apply(subareaPop$ZUrbanFineScalePrevalence, 1, sd, na.rm=TRUE), 
                          apply(subareaPop$ZRuralFineScalePrevalence, 1, sd, na.rm=TRUE))
      countSDConStrat = countSDConStrat[!undefinedPrevalenceConStrats]
      rangePrevalenceSDConStrat = range(prevalenceSDConStrat)
      rangeCountSDConStrat = range(countSDConStrat)
      
      # get credible interval widths for the lcpb model
      prevalenceCIWidthConStratlcpb = c(apply(subareaPop$pUrbanSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}), 
                                        apply(subareaPop$pRuralSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}))
      prevalenceCIWidthConStratlcpb = prevalenceCIWidthConStratlcpb[!undefinedPrevalenceConStrats]
      countCIWidthConStratlcpb = c(apply(subareaPop$ZUrbanSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}), 
                                   apply(subareaPop$ZRuralSmoothRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}))
      countCIWidthConStratlcpb = countCIWidthConStratlcpb[!undefinedPrevalenceConStrats]
      
      prevalenceSDConStratlcpb = c(apply(subareaPop$pUrbanSmoothRisk, 1, sd, na.rm=TRUE), 
                                   apply(subareaPop$pRuralSmoothRisk, 1, sd, na.rm=TRUE))
      prevalenceSDConStratlcpb = prevalenceSDConStratlcpb[!undefinedPrevalenceConStrats]
      countSDConStratlcpb = c(apply(subareaPop$ZUrbanSmoothRisk, 1, sd, na.rm=TRUE), 
                              apply(subareaPop$ZRuralSmoothRisk, 1, sd, na.rm=TRUE))
      countSDConStratlcpb = countSDConStratlcpb[!undefinedPrevalenceConStrats]
      
      # do the same for the LCPb model
      prevalenceCIWidthConStratLCPb = c(apply(subareaPop$pUrbanFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}), 
                                        apply(subareaPop$pRuralFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}))
      prevalenceCIWidthConStratLCPb = prevalenceCIWidthConStratLCPb[!undefinedPrevalenceConStrats]
      countCIWidthConStratLCPb = c(apply(subareaPop$ZUrbanFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}), 
                                   apply(subareaPop$ZRuralFineScaleRisk, 1, function(x) {diff(quantile(x, probs=c(alpha/2, 1-alpha/2), na.rm=TRUE))}))
      countCIWidthConStratLCPb = countCIWidthConStratLCPb[!undefinedPrevalenceConStrats]
      
      prevalenceSDConStratLCPb = c(apply(subareaPop$pUrbanFineScaleRisk, 1, sd, na.rm=TRUE), 
                                   apply(subareaPop$pRuralFineScaleRisk, 1, sd, na.rm=TRUE))
      prevalenceSDConStratLCPb = prevalenceSDConStratLCPb[!undefinedPrevalenceConStrats]
      countSDConStratLCPb = c(apply(subareaPop$ZUrbanFineScaleRisk, 1, sd, na.rm=TRUE), 
                              apply(subareaPop$ZRuralFineScaleRisk, 1, sd, na.rm=TRUE))
      countSDConStratLCPb = countSDConStratLCPb[!undefinedPrevalenceConStrats]
    } else if(thisLevel == "constituency") {
      urbanConstituencies = poppsub$County == "Nairobi" | poppsub$County == "Mombasa"
      undefinedRelativePrevalenceConstituencies = (poppsub$popUrb == 0) | (poppsub$popRur == 0)
      rangePrevalencePredConstituency = range(c(rangePrevalencePredConstituency, 
                                                rowMeans(subareaPop$pFineScaleRisk, na.rm=TRUE)))
      rangeCountPredConstituency = range(c(rangeCountPredConstituency, 
                                           rowMeans(subareaPop$ZFineScaleRisk)), na.rm=TRUE)
      relativePrevalencePredConstituency = rowMeans(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, na.rm=TRUE)
      undefinedRelativePrevalenceConstituencies = undefinedRelativePrevalenceConstituencies | !is.finite(relativePrevalencePredConstituency)
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
  
  # make plots ----
  
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
  
  ## 3 x 2 plot of predictions (prevalence mean and CI width) ----
  
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
  png(paste0(figDirectory, "application/prevalenceMeanCIWidth", logisticText, coarseText, resultTypeText, ".png"), width=1500, height=1000)
  par(mfrow=c(2,3), oma=c(5,5,5,5), mar=c(3.1, 5.1, 1.1, 7.1))
  
  # pixel level
  quilt.plot(popMat$lon, popMat$lat, pixelMean, FUN=function(x){logit(mean(x, na.rm=TRUE))}, 
             zlim=logit(meanRange), nx=155, ny=160, main="", cex.main=3, col=meanCols, 
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
  
  quilt.plot(popMat$lon, popMat$lat, prevalenceCIWidthPixellcpb, FUN=function(x){logit(mean(x, na.rm=TRUE))}, 
             nx=155, ny=160, main="", cex.main=3, col=sdCols, 
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
  ## same plot but for burden ----
  pixelMean = rowMeans(pixelPop$ZSmoothRisk)
  constituencyMean = rowMeans(subareaPop$ZFineScaleRisk)
  countyMean = rowMeans(areaPop$ZFineScaleRisk)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPb$regionMatrices$p)
  meanRange = c(.01, max(c(pixelMean, constituencyMean, countyMean)))
  meanRangePixel = c(.01, max(pixelMean))
  meanRangeConstituency = range(constituencyMean)
  meanRangeCounty = range(countyMean)
  
  png(paste0(figDirectory, "application/burdenMeanCIWidth", logisticText, coarseText, resultTypeText, ".png"), width=1500, height=1000)
  par(mfrow=c(2,3), oma=c(5,5,5,5), mar=c(3.1, 5.1, 1.1, 7.1))
  
  # pixel level
  quilt.plot(popMat$lon, popMat$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
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
  
  quilt.plot(popMat$lon, popMat$lat, countCIWidthPixellcpb, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
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
  
  ## same plot but for relative prevalence ----
  # Don't include pixel level because pixel is either urban or rural, not both
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
  
  png(paste0(figDirectory, "application/relativePrevMeanCIWidth", logisticText, coarseText, resultTypeText, ".png"), width=1000, height=1000)
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
  
  # old plots ----
  
  ## plot prevalence credible interval widths ----
  pixelWidth = prevalenceCIWidthPixel
  constituencyWidth = apply(subareaPop$pFineScalePrevalence, 1, getWidth)
  countyWidth = apply(areaPop$pFineScalePrevalence, 1, getWidth)
  # provinceWidth = apply(agg$aggregatedResultsLCPB$regionMatrices$p, 1, getWidth)
  widthRange = range(c(rangePrevalenceCIWidthConstituency, 
                       rangePrevalenceCIWidthCounty))
  widthRangePixel = rangePrevalenceCIWidthPixel
  
  png(paste0(figDirectory, "application/prevalenceCIWidth", logisticText, coarseText, resultTypeText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popMat$lon, popMat$lat, pixelWidth, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
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
  
  ## 2 x 2 plot of burden predictions ----
  
  # plot means
  pixelMean = rowMeans(pixelPop$ZFineScaleRisk)
  constituencyMean = rowMeans(subareaPop$ZFineScalePrevalence)
  countyMean = rowMeans(areaPop$ZFineScalePrevalence)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPB$regionMatrices$Z)
  meanRangePixel = range(pixelMean[pixelMean >= .05])
  meanRangeConstituency = rangeCountPredConstituency
  meanRangeCounty = rangeCountPredCounty
  # meanRangeProvince = rangeCountPredProvince
  
  png(paste0(figDirectory, "application/countMean", logisticText, coarseText, resultTypeText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popMat$lon, popMat$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
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
  
  ## plot credible interval widths burden ----
  pixelWidth = countCIWidthPixel
  constituencyWidth = apply(subareaPop$ZFineScalePrevalence, 1, getWidth)
  countyWidth = apply(areaPop$ZFineScalePrevalence, 1, getWidth)
  # provinceWidth = apply(agg$aggregatedResultsLCPB$regionMatrices$Z, 1, getWidth)
  widthRangePixel = rangeCountCIWidthPixel # range(pixelWidth)
  widthRangeConstituency = rangeCountCIWidthConstituency
  widthRangeCounty = rangeCountCIWidthCounty
  # widthRangeProvince = rangeCountCIWidthProvince
  
  png(paste0(figDirectory, "application/countWidth", logisticText, coarseText, resultTypeText, ".png"), width=1000, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,4,7), mar=c(6.1, 8.5, 1.1, 3.5))
  
  # pixel level
  quilt.plot(popMat$lon, popMat$lat, pixelWidth, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
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
  
  ## 2 x 3 plot of predictions (relative prevalence) ----
  
  # plot mean
  constituencyMean = relativePrevalencePredConstituency
  countyMean = relativePrevalencePredCounty
  # provinceMean = relativePrevalencePredProvince
  meanRange = range(c(rangeRelativePrevalencePredConstituency, 
                      rangeRelativePrevalencePredCounty))
  urbDivergingCols = makeGreenBlueDivergingColors(64, center=1, valRange=meanRange)
  
  png(paste0(figDirectory, "application/relativePrevalence", logisticText, coarseText, resultTypeText, ".png"), width=1000, height=1000)
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
  
  ## Boxplots of relative CI widths to the stratified model ----
  # modelNames: S, SC, SCP (stratified, cluster, population)
  
  ### prevalence (with abs width) ----
  
  # absolute prevalence uncertainty
  pdf(paste0(figDirectory, "application/prevalenceWidthBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=6, height=6)
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
  
  ### prevalence ----
  pdf(paste0(figDirectory, "application/prevalenceRelWidthBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=8, height=5)
  par(mfrow=c(1,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  # pixels = length(prevalenceCIWidthPixel)
  # pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
  #                       width=c(100*(prevalenceCIWidthPixellcpb-prevalenceCIWidthPixellcpb)/prevalenceCIWidthPixellcpb, 
  #                               100*(prevalenceCIWidthPixelLCPb-prevalenceCIWidthPixellcpb)/prevalenceCIWidthPixellcpb, 
  #                               100*(prevalenceCIWidthPixel-prevalenceCIWidthPixellcpb)/prevalenceCIWidthPixellcpb))
  # boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
  #         main="Pixel", col="skyblue")
  # abline(h=0, lty=2)
  
  constituencies = length(prevalenceCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(prevalenceCIWidthConstituencylcpb-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb, 
                                       100*(prevalenceCIWidthConstituencyLCPb-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb, 
                                       100*(prevalenceCIWidthConstituency-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  counties = length(prevalenceCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(prevalenceCIWidthCountylcpb-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb, 
                                 100*(prevalenceCIWidthCountyLCPb-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb, 
                                 100*(prevalenceCIWidthCounty-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
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
  ### prevalence (SD) ----
  pdf(paste0(figDirectory, "application/prevalenceRelSDBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=8, height=5)
  par(mfrow=c(1,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  # pixels = length(prevalenceSDPixel)
  # pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
  #                       width=c(100*(prevalenceSDPixellcpb-prevalenceSDPixellcpb)/prevalenceSDPixellcpb, 
  #                               100*(prevalenceSDPixelLCPb-prevalenceSDPixellcpb)/prevalenceSDPixellcpb, 
  #                               100*(prevalenceSDPixel-prevalenceSDPixellcpb)/prevalenceSDPixellcpb))
  # boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
  #         main="Pixel", col="skyblue")
  # abline(h=0, lty=2)
  
  constituencies = length(prevalenceSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(prevalenceSDConstituencylcpb-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb, 
                                       100*(prevalenceSDConstituencyLCPb-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb, 
                                       100*(prevalenceSDConstituency-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  counties = length(prevalenceSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(prevalenceSDCountylcpb-prevalenceSDCountylcpb)/prevalenceSDCountylcpb, 
                                 100*(prevalenceSDCountyLCPb-prevalenceSDCountylcpb)/prevalenceSDCountylcpb, 
                                 100*(prevalenceSDCounty-prevalenceSDCountylcpb)/prevalenceSDCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
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
  
  ### prevalence (SD, con strat) ----
  pdf(paste0(figDirectory, "application/prevalenceRelSDBoxplotConStrat", logisticText, coarseText, resultTypeText, ".pdf"), width=9, height=6)
  par(mfrow=c(2,3), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  conStrats = length(prevalenceSDConStratlcpb)
  conStratDat = data.frame(modelName=c(rep("Smooth latent", conStrats), rep("Latent", conStrats), rep("Empirical", conStrats)), 
                           width=c(prevalenceSDConStratlcpb, prevalenceSDConStratLCPb, prevalenceSDConStrat))
  boxplot(width~modelName, data=conStratDat, ylab="", xlab="", main="Constituency x stratum", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Prevalence SD", line = 3, cex=1)
  
  constituencies = length(prevalenceSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth latent", constituencies), rep("Latent", constituencies), rep("Empirical", constituencies)), 
                               width=c(prevalenceSDConstituencylcpb, prevalenceSDConstituencyLCPb, prevalenceSDConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  
  counties = length(prevalenceSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth latent", counties), rep("Latent", counties), rep("Empirical", counties)), 
                         width=c(prevalenceSDCountylcpb, prevalenceSDCountyLCPb, prevalenceSDCounty))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  
  conStrats = length(prevalenceSDConStratlcpb)
  conStratDatRel = data.frame(modelName=c(rep("Smooth latent", conStrats), rep("latent", conStrats), rep("Empirical", conStrats)), 
                              width=c(100*(prevalenceSDConStratlcpb-prevalenceSDConStratlcpb)/prevalenceSDConStratlcpb, 
                                      100*(prevalenceSDConStratLCPb-prevalenceSDConStratlcpb)/prevalenceSDConStratlcpb, 
                                      100*(prevalenceSDConStrat-prevalenceSDConStratlcpb)/prevalenceSDConStratlcpb))
  boxplot(width~modelName, data=conStratDatRel, ylab="", xlab="", main="Constituency x stratum", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase SD", line = 3, cex=1)
  
  constituencies = length(prevalenceSDConstituency)
  constituencyDatRel = data.frame(modelName=c(rep("Smooth latent", constituencies), rep("Latent", constituencies), rep("Empirical", constituencies)), 
                                  width=c(100*(prevalenceSDConstituencylcpb-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb, 
                                          100*(prevalenceSDConstituencyLCPb-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb, 
                                          100*(prevalenceSDConstituency-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDatRel, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  counties = length(prevalenceSDCounty)
  countyDatRel = data.frame(modelName=c(rep("Smooth latent", counties), rep("Latent", counties), rep("Empirical", counties)), 
                            width=c(100*(prevalenceSDCountylcpb-prevalenceSDCountylcpb)/prevalenceSDCountylcpb, 
                                    100*(prevalenceSDCountyLCPb-prevalenceSDCountylcpb)/prevalenceSDCountylcpb, 
                                    100*(prevalenceSDCounty-prevalenceSDCountylcpb)/prevalenceSDCountylcpb))
  boxplot(width~modelName, data=countyDatRel, ylab="", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(countSDProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(countSDProvincelcpb-countSDProvincelcpb)/countSDProvincelcpb, 
  #                                  100*(countSDProvinceLCPb-countSDProvincelcpb)/countSDProvincelcpb, 
  #                                  100*(countSDProvince-countSDProvincelcpb)/countSDProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  # mtext(side = 3, "Total Deaths Relative SD", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  
  ### burden (abs scale) ----
  pdf(paste0(figDirectory, "application/countWidthBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=8, height=5)
  par(mfrow=c(1,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  # pixels = length(countCIWidthPixel)
  # pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
  #                       width=c(countCIWidthPixellcpb, countCIWidthPixelLCPb, countCIWidthPixel))
  # boxplot(width~modelName, data=pixelDat, ylab="", xlab="Model", main="Pixel", col="skyblue", log="y")
  
  constituencies = length(countCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(countCIWidthConstituencylcpb, countCIWidthConstituencyLCPb, countCIWidthConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="95% CI Width", xlab="Model", main="Constituency", col="skyblue", log="y")
  mtext(side = 2, "95% CI width", line = 3, cex=1)
  
  counties = length(countCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(countCIWidthCountylcpb, countCIWidthCountyLCPb, countCIWidthCounty))
  boxplot(width~modelName, data=countyDat, ylab="95% CI Width", xlab="Model", main="County", col="skyblue", log="y")
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(countCIWidthProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(countCIWidthProvincelcpb, countCIWidthProvinceLCPb, countCIWidthProvince))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue", log="y")
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  mtext(side = 3, "Total Deaths 95% CI Width", line = 0.5, cex=1, outer=TRUE)
  dev.off()
  
  ### burden ----
  pdf(paste0(figDirectory, "application/countRelWidthBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=8, height=5)
  par(mfrow=c(1,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  # pixels = length(countCIWidthPixel)
  # pixelDat = data.frame(modelName=c(rep("Smooth Risk", pixels), rep("Risk", pixels), rep("Prevalence", pixels)), 
  #                       width=c(100*(countCIWidthPixellcpb-countCIWidthPixellcpb)/countCIWidthPixellcpb, 
  #                               100*(countCIWidthPixelLCPb-countCIWidthPixellcpb)/countCIWidthPixellcpb, 
  #                               100*(countCIWidthPixel-countCIWidthPixellcpb)/countCIWidthPixellcpb))
  # boxplot(width~modelName, data=pixelDat, ylab="", xlab="", 
  #         main="Pixel", col="skyblue")
  # abline(h=0, lty=2)
  
  constituencies = length(countCIWidthConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                               width=c(100*(countCIWidthConstituencylcpb-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb, 
                                       100*(countCIWidthConstituencyLCPb-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb, 
                                       100*(countCIWidthConstituency-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  counties = length(countCIWidthCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                         width=c(100*(countCIWidthCountylcpb-countCIWidthCountylcpb)/countCIWidthCountylcpb, 
                                 100*(countCIWidthCountyLCPb-countCIWidthCountylcpb)/countCIWidthCountylcpb, 
                                 100*(countCIWidthCounty-countCIWidthCountylcpb)/countCIWidthCountylcpb))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
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
  
  ### Burden (SD) ----
  pdf(paste0(figDirectory, "application/countRelSDBoxplotWithPixel", logisticText, coarseText, resultTypeText, ".pdf"), width=6, height=6)
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
  
  ### Burden (SD, no pixel) ----
  pdf(paste0(figDirectory, "application/countRelSDBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=6, height=6)
  par(mfrow=c(2,2), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  constituencies = length(countSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth Latent", constituencies), rep("Latent", constituencies), rep("Empirical", constituencies)), 
                               width=c(countSDConstituencylcpb, countSDConstituencyLCPb, countSDConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "SD", line = 3, cex=1)
  
  counties = length(countSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth Latent", counties), rep("Latent", counties), rep("Empirical", counties)), 
                         width=c(countSDCountylcpb, countSDCountyLCPb, countSDCounty))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  
  constituencies = length(countSDConstituency)
  constituencyDatRel = data.frame(modelName=c(rep("Smooth Risk", constituencies), rep("Risk", constituencies), rep("Prevalence", constituencies)), 
                                  width=c(100*(countSDConstituencylcpb-countSDConstituencylcpb)/countSDConstituencylcpb, 
                                          100*(countSDConstituencyLCPb-countSDConstituencylcpb)/countSDConstituencylcpb, 
                                          100*(countSDConstituency-countSDConstituencylcpb)/countSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDatRel, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase", line = 3, cex=1)
  
  counties = length(countSDCounty)
  countyDatRel = data.frame(modelName=c(rep("Smooth Risk", counties), rep("Risk", counties), rep("Prevalence", counties)), 
                            width=c(100*(countSDCountylcpb-countSDCountylcpb)/countSDCountylcpb, 
                                    100*(countSDCountyLCPb-countSDCountylcpb)/countSDCountylcpb, 
                                    100*(countSDCounty-countSDCountylcpb)/countSDCountylcpb))
  boxplot(width~modelName, data=countyDatRel, ylab="", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  ### Burden (SD, no pixel, con strat) ----
  pdf(paste0(figDirectory, "application/countRelSDBoxplotConStrat", logisticText, coarseText, resultTypeText, ".pdf"), width=9, height=6)
  par(mfrow=c(2,3), oma=c(3,3,2,0), mar=c(2, 2, 2, 1))
  
  conStrats = length(countSDConStratlcpb)
  conStratDat = data.frame(modelName=c(rep("Smooth latent", conStrats), rep("Latent", conStrats), rep("Empirical", conStrats)), 
                           width=c(countSDConStratlcpb, countSDConStratLCPb, countSDConStrat))
  boxplot(width~modelName, data=conStratDat, ylab="", xlab="", main="Constituency x stratum", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 2, "Burden SD", line = 3, cex=1)
  
  constituencies = length(countSDConstituency)
  constituencyDat = data.frame(modelName=c(rep("Smooth latent", constituencies), rep("Latent", constituencies), rep("Empirical", constituencies)), 
                               width=c(countSDConstituencylcpb, countSDConstituencyLCPb, countSDConstituency))
  boxplot(width~modelName, data=constituencyDat, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  
  counties = length(countSDCounty)
  countyDat = data.frame(modelName=c(rep("Smooth latent", counties), rep("Latent", counties), rep("Empirical", counties)), 
                         width=c(countSDCountylcpb, countSDCountyLCPb, countSDCounty))
  boxplot(width~modelName, data=countyDat, ylab="Percent Increase", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  
  conStrats = length(countSDConStratlcpb)
  conStratDatRel = data.frame(modelName=c(rep("Smooth latent", conStrats), rep("latent", conStrats), rep("Empirical", conStrats)), 
                              width=c(100*(countSDConStratlcpb-countSDConStratlcpb)/countSDConStratlcpb, 
                                      100*(countSDConStratLCPb-countSDConStratlcpb)/countSDConStratlcpb, 
                                      100*(countSDConStrat-countSDConStratlcpb)/countSDConStratlcpb))
  boxplot(width~modelName, data=conStratDatRel, ylab="", xlab="", main="Constituency x stratum", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  mtext(side = 2, "Percent increase SD", line = 3, cex=1)
  
  constituencies = length(countSDConstituency)
  constituencyDatRel = data.frame(modelName=c(rep("Smooth latent", constituencies), rep("Latent", constituencies), rep("Empirical", constituencies)), 
                                  width=c(100*(countSDConstituencylcpb-countSDConstituencylcpb)/countSDConstituencylcpb, 
                                          100*(countSDConstituencyLCPb-countSDConstituencylcpb)/countSDConstituencylcpb, 
                                          100*(countSDConstituency-countSDConstituencylcpb)/countSDConstituencylcpb))
  boxplot(width~modelName, data=constituencyDatRel, ylab="", xlab="", main="Constituency", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  counties = length(countSDCounty)
  countyDatRel = data.frame(modelName=c(rep("Smooth latent", counties), rep("Latent", counties), rep("Empirical", counties)), 
                            width=c(100*(countSDCountylcpb-countSDCountylcpb)/countSDCountylcpb, 
                                    100*(countSDCountyLCPb-countSDCountylcpb)/countSDCountylcpb, 
                                    100*(countSDCounty-countSDCountylcpb)/countSDCountylcpb))
  boxplot(width~modelName, data=countyDatRel, ylab="", xlab="", main="County", col="skyblue")
  abline(h=0, lty=2)
  mtext(side = 1, "Model", line = 3, cex=1)
  
  # provinces = length(countSDProvince)
  # provinceDat = data.frame(modelName=c(rep("Smooth Risk", provinces), rep("Risk", provinces), rep("Prevalence", provinces)), 
  #                          width=c(100*(countSDProvincelcpb-countSDProvincelcpb)/countSDProvincelcpb, 
  #                                  100*(countSDProvinceLCPb-countSDProvincelcpb)/countSDProvincelcpb, 
  #                                  100*(countSDProvince-countSDProvincelcpb)/countSDProvincelcpb))
  # boxplot(width~modelName, data=provinceDat, ylab="", xlab="", main="Province", col="skyblue")
  # abline(h=0, lty=2)
  # mtext(side = 1, "Model", line = 3, cex=1)
  
  # add title in the left margin
  # mtext(side = 3, "Total Deaths Relative SD", line = 0.5, cex=1, outer=TRUE)
  
  dev.off()
  
  ### Relative prevalence CI width (abs scale) ----
  pdf(paste0(figDirectory, "application/relativePrevalenceWidthBoxplot", logisticText, coarseText, resultTypeText, ".pdf"), width=6, height=6)
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
  
  ### relative prevalence (SD) ----
  pdf(paste0(figDirectory, "application/relativePrevalenceWidthBoxplotSD", logisticText, coarseText, resultTypeText, ".pdf"), width=6, height=6)
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
  
  ## N EAs vs % Increase SD ----
  
  # Calculate number of EAs per area
  easpa = makeDefaultEASPA()
  meanEAspSub = meanEAsPerCon()
  meanEAsConStrat = c(meanEAspSub$meanUrbanEAs, meanEAspSub$meanRuralEAs)[!undefinedPrevalenceConStrats]
  meanEAsConstituency = c(meanEAspSub$meanTotalEAs)
  meanEAsCounty = easpa$EATotal
  meanEAs = c(meanEAsConStrat, meanEAsConstituency, meanEAsCounty)
  
  # precompute plotting values
  cols = c("blue", "purple", "red") # for conStrat, constituency, county
  pchs = c(16, 15, 17)
  cexs=c(.35, .5, .9)
  allCols = c(rep(cols[1], length(meanEAsConStrat)), rep(cols[2], length(meanEAsConstituency)), rep(cols[3], length(meanEAsCounty)))
  allPchs = c(rep(pchs[1], length(meanEAsConStrat)), rep(pchs[2], length(meanEAsConstituency)), rep(pchs[3], length(meanEAsCounty)))
  
  ### prevalence ----
  pdf(paste0(figDirectory, "application/prevalencePctIncreaseSDvsEAs", logisticText, coarseText, resultTypeText, ".pdf"), width=5, height=5)
  
  conStratDat = 100*(prevalenceCIWidthConStrat-prevalenceCIWidthConStratlcpb)/prevalenceCIWidthConStratlcpb
  constituencyDat = 100*(prevalenceCIWidthConstituency-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb
  countyDat = 100*(prevalenceCIWidthCounty-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb
  allDat = c(conStratDat, constituencyDat, countyDat)
  minVal = min(allDat[allDat > 0])*.95
  # conStratDat[conStratDat < 0] = minVal
  # constituencyDat[constituencyDat < 0] = minVal
  # countyDat[countyDat < 0] = minVal
  xRange = range(meanEAs)
  yRange = c(minVal, max(allDat))
  plot(meanEAsConStrat[conStratDat>0], conStratDat[conStratDat>0], col=cols[1], pch=pchs[1], main="Prevalence SD vs. number EAs", 
       xlab="Mean number EAs", ylab="Percent increase SD", log="xy", xlim=xRange, ylim=yRange, cex=cexs[1])
  # points(meanEAsConStrat[conStratDat<=0], rep(minVal, sum(conStratDat<=0)), col=cols[1], pch=pchs[1], cex=cexs[1])
  points(meanEAsConstituency[constituencyDat>0], constituencyDat[constituencyDat>0], col=cols[2], pch=pchs[2], cex=cexs[2])
  # points(meanEAsConstituency[constituencyDat<=0], rep(minVal, sum(constituencyDat<=0)), col=cols[2], pch=pchs[2], cex=cexs[2])
  points(meanEAsCounty[countyDat>0], countyDat[countyDat>0], col=cols[3], pch=pchs[3], cex=cexs[3])
  # points(meanEAsCounty[countyDat<=0], rep(minVal, sum(countyDat<=0)), col=cols[3], pch=pchs[3], cex=cexs[3])
  
  legend("bottomleft", c("Constituency x stratum", "Constituency", "County"), 
         pch=pchs, col=cols, cex=.7)
  
  dev.off()
  
  ### burden ----
  pdf(paste0(figDirectory, "application/countPctIncreaseSDvsEAs", logisticText, coarseText, resultTypeText, ".pdf"), width=5, height=5)
  
  conStratDat = 100*(countCIWidthConStrat-countCIWidthConStratlcpb)/countCIWidthConStratlcpb
  constituencyDat = 100*(countCIWidthConstituency-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb
  countyDat = 100*(countCIWidthCounty-countCIWidthCountylcpb)/countCIWidthCountylcpb
  allDat = c(conStratDat, constituencyDat, countyDat)
  minVal = min(allDat[allDat > 0])*.95
  # conStratDat[conStratDat < 0] = minVal
  # constituencyDat[constituencyDat < 0] = minVal
  # countyDat[countyDat < 0] = minVal
  xRange = range(meanEAs)
  yRange = c(minVal, max(allDat))
  plot(meanEAsConStrat[conStratDat>0], conStratDat[conStratDat>0], col=cols[1], pch=pchs[1], main="Burden SD vs. number EAs", 
       xlab="Mean number EAs", ylab="Percent increase SD", log="xy", xlim=xRange, ylim=yRange, cex=cexs[1])
  # points(meanEAsConStrat[conStratDat<=0], rep(minVal, sum(conStratDat<=0)), col=cols[1], pch=pchs[1], cex=cexs[1])
  points(meanEAsConstituency[constituencyDat>0], constituencyDat[constituencyDat>0], col=cols[2], pch=pchs[2], cex=cexs[2])
  # points(meanEAsConstituency[constituencyDat<=0], rep(minVal, sum(constituencyDat<=0)), col=cols[2], pch=pchs[2], cex=cexs[2])
  points(meanEAsCounty[countyDat>0], countyDat[countyDat>0], col=cols[3], pch=pchs[3], cex=cexs[3])
  # points(meanEAsCounty[countyDat<=0], rep(minVal, sum(countyDat<=0)), col=cols[3], pch=pchs[3], cex=cexs[3])
  
  legend("bottomleft", c("Constituency x stratum", "Constituency", "County"), 
         pch=pchs, col=cols, cex=.7)
  
  dev.off()
  
  ### relative prevalence ----
  pdf(paste0(figDirectory, "application/relativePrevalencePctIncreaseSDvsEAs", logisticText, coarseText, resultTypeText, ".pdf"), width=5, height=5)
  
  constituencyDat = 100*(relativePrevalenceCIWidthConstituency-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb
  countyDat = 100*(relativePrevalenceCIWidthCounty-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb
  allDat = c(constituencyDat, countyDat)
  minVal = min(allDat)
  # conStratDat[conStratDat < 0] = minVal
  # constituencyDat[constituencyDat < 0] = minVal
  # countyDat[countyDat < 0] = minVal
  xRange = range(meanEAs)
  yRange = range(allDat, na.rm=TRUE)
  plot(meanEAsConstituency, constituencyDat, col=cols[2], pch=pchs[2], main="Relative prevalence SD vs. number EAs", 
       xlab="Mean number EAs", ylab="Percent increase SD", log="xy", xlim=xRange, ylim=yRange, cex=cexs[2])
  points(meanEAsCounty, countyDat, col=cols[3], pch=pchs[3], cex=cexs[3])
  
  legend("bottomleft", c("Constituency", "County"), 
         pch=pchs[-1], col=cols[-1], cex=.7, bg="white")
  
  dev.off()
  
  ### All ----
  pdf(paste0(figDirectory, "application/allPctIncreaseSDvsEAs", logisticText, coarseText, resultTypeText, ".pdf"), width=9, height=4)
  options(scipen=5)
  par(mfrow=c(1, 3), mar=c(4, 2.5, 3, 1), oma=c(0, 2, 0, 0))
  
  # set common y range
  conStratDat = 100*(prevalenceCIWidthConStrat-prevalenceCIWidthConStratlcpb)/prevalenceCIWidthConStratlcpb
  constituencyDat = 100*(prevalenceCIWidthConstituency-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb
  countyDat = 100*(prevalenceCIWidthCounty-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb
  allDat = c(conStratDat, constituencyDat, countyDat)
  minVal = min(allDat[allDat > 0])*.95
  yRange = c(minVal, max(allDat))
  
  conStratDat = 100*(countCIWidthConStrat-countCIWidthConStratlcpb)/countCIWidthConStratlcpb
  constituencyDat = 100*(countCIWidthConstituency-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb
  countyDat = 100*(countCIWidthCounty-countCIWidthCountylcpb)/countCIWidthCountylcpb
  allDat = c(conStratDat, constituencyDat, countyDat)
  minVal = min(allDat[allDat > 0])*.95
  yRange = range(c(yRange, minVal, max(allDat)))
  
  constituencyDat = 100*(relativePrevalenceCIWidthConstituency-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb
  countyDat = 100*(relativePrevalenceCIWidthCounty-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb
  allDat = c(constituencyDat, countyDat)
  yRange = range(c(yRange, range(allDat, na.rm=TRUE)))
  
  # prevalence
  conStratDat = 100*(prevalenceCIWidthConStrat-prevalenceCIWidthConStratlcpb)/prevalenceCIWidthConStratlcpb
  constituencyDat = 100*(prevalenceCIWidthConstituency-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb
  countyDat = 100*(prevalenceCIWidthCounty-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb
  allDat = c(conStratDat, constituencyDat, countyDat)
  minVal = min(allDat[allDat > 0])*.95
  # conStratDat[conStratDat < 0] = minVal
  # constituencyDat[constituencyDat < 0] = minVal
  # countyDat[countyDat < 0] = minVal
  xRange = range(meanEAs)
  # yRange = c(minVal, max(allDat))
  plot(meanEAsConStrat[conStratDat>0], conStratDat[conStratDat>0], col=cols[1], pch=pchs[1], main="Prevalence", 
       xlab="", ylab="", log="xy", xlim=xRange, ylim=yRange, cex=cexs[1])
  abline(v=300, lty=2)
  # points(meanEAsConStrat[conStratDat<=0], rep(minVal, sum(conStratDat<=0)), col=cols[1], pch=pchs[1], cex=cexs[1])
  points(meanEAsConstituency[constituencyDat>0], constituencyDat[constituencyDat>0], col=cols[2], pch=pchs[2], cex=cexs[2])
  # points(meanEAsConstituency[constituencyDat<=0], rep(minVal, sum(constituencyDat<=0)), col=cols[2], pch=pchs[2], cex=cexs[2])
  points(meanEAsCounty[countyDat>0], countyDat[countyDat>0], col=cols[3], pch=pchs[3], cex=cexs[3])
  # points(meanEAsCounty[countyDat<=0], rep(minVal, sum(countyDat<=0)), col=cols[3], pch=pchs[3], cex=cexs[3])
  mtext("Percent increase SD", side=2, line=3, cex=cexs[3])
  # legend("bottomleft", c("Constituency x stratum", "Constituency", "County"), 
  #        pch=pchs, col=cols, cex=.7)
  
  # burden
  conStratDat = 100*(countCIWidthConStrat-countCIWidthConStratlcpb)/countCIWidthConStratlcpb
  constituencyDat = 100*(countCIWidthConstituency-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb
  countyDat = 100*(countCIWidthCounty-countCIWidthCountylcpb)/countCIWidthCountylcpb
  allDat = c(conStratDat, constituencyDat, countyDat)
  minVal = min(allDat[allDat > 0])*.95
  # conStratDat[conStratDat < 0] = minVal
  # constituencyDat[constituencyDat < 0] = minVal
  # countyDat[countyDat < 0] = minVal
  xRange = range(meanEAs)
  # yRange = c(minVal, max(allDat))
  plot(meanEAsConStrat[conStratDat>0], conStratDat[conStratDat>0], col=cols[1], pch=pchs[1], main="Burden", 
       xlab="", ylab="", log="xy", xlim=xRange, ylim=yRange, cex=cexs[1])
  abline(v=300, lty=2)
  # points(meanEAsConStrat[conStratDat<=0], rep(minVal, sum(conStratDat<=0)), col=cols[1], pch=pchs[1], cex=cexs[1])
  points(meanEAsConstituency[constituencyDat>0], constituencyDat[constituencyDat>0], col=cols[2], pch=pchs[2], cex=cexs[2])
  # points(meanEAsConstituency[constituencyDat<=0], rep(minVal, sum(constituencyDat<=0)), col=cols[2], pch=pchs[2], cex=cexs[2])
  points(meanEAsCounty[countyDat>0], countyDat[countyDat>0], col=cols[3], pch=pchs[3], cex=cexs[3])
  # points(meanEAsCounty[countyDat<=0], rep(minVal, sum(countyDat<=0)), col=cols[3], pch=pchs[3], cex=cexs[3])
  
  legend("top", c("Constituency x stratum", "Constituency", "County"), 
         pch=pchs, col=cols, cex=.85, bg="white")
  mtext("Mean number EAs", side=1, line=3, cex=cexs[3])
  
  # relative prevalence
  constituencyDat = 100*(relativePrevalenceCIWidthConstituency-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb
  countyDat = 100*(relativePrevalenceCIWidthCounty-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb
  allDat = c(constituencyDat, countyDat)
  minVal = min(allDat)
  # conStratDat[conStratDat < 0] = minVal
  # constituencyDat[constituencyDat < 0] = minVal
  # countyDat[countyDat < 0] = minVal
  xRange = range(meanEAs)
  # yRange = range(allDat, na.rm=TRUE)
  plot(meanEAsConstituency, constituencyDat, col=cols[2], pch=pchs[2], main="Relative prevalence", 
       xlab="", ylab="", log="xy", xlim=xRange, ylim=yRange, cex=cexs[2])
  # abline(v=3000, lty=2)
  points(meanEAsCounty, countyDat, col=cols[3], pch=pchs[3], cex=cexs[3])
  
  # legend("bottomleft", c("Constituency", "County"), 
  #        pch=pchs[-1], col=cols[-1], cex=.7)
  
  dev.off()
  
  # tables ----
  
  # parameter table
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
  
  # table of mean pct increase SD for each type of area
  # 3 sets columns: prev, burd, relPrev
  # 3 rows: constituency x stratum, constituency, county
  conStratDat = 100*(prevalenceSDConStrat-prevalenceSDConStratlcpb)/prevalenceSDConStratlcpb
  constituencyDat = 100*(prevalenceSDConstituency-prevalenceSDConstituencylcpb)/prevalenceSDConstituencylcpb
  countyDat = 100*(prevalenceSDCounty-prevalenceSDCountylcpb)/prevalenceSDCountylcpb
  prevConStrat = mean(conStratDat)
  prevCon = mean(constituencyDat)
  prevCounty = mean(countyDat)
  prevMeans = c("conStrat"=prevConStrat, "constituency"=prevCon, "county"=prevCounty)
  
  conStratDat = 100*(countSDConStrat-countSDConStratlcpb)/countSDConStratlcpb
  constituencyDat = 100*(countSDConstituency-countSDConstituencylcpb)/countSDConstituencylcpb
  countyDat = 100*(countSDCounty-countSDCountylcpb)/countSDCountylcpb
  burdConStrat = mean(conStratDat)
  burdCon = mean(constituencyDat)
  burdCounty = mean(countyDat)
  burdMeans = c("conStrat"=burdConStrat, "constituency"=burdCon, "county"=burdCounty)
  
  constituencyDat = 100*(relativePrevalenceSDConstituency-relativePrevalenceSDConstituencylcpb)/relativePrevalenceSDConstituencylcpb
  countyDat = 100*(relativePrevalenceSDCounty-relativePrevalenceSDCountylcpb)/relativePrevalenceSDCountylcpb
  relPrevCon = mean(constituencyDat, na.rm=TRUE)
  relPrevCounty = mean(countyDat, na.rm=TRUE)
  relPrevMeans = c("conStrat"=NA, "constituency"=relPrevCon, "county"=relPrevCounty)
  
  tab = cbind(prevMeans, burdMeans, relPrevMeans)
  colnames(tab) = c("Prevalence", "Burden", "Relative prevalence")
  row.names(tab) = c("Constituency x stratum", "Constituency", "County")
  tab = round(tab, 0)
  
  print(xtable(tab, digits=0, 
               caption=paste0("Summary statistics for risk model parameters when ", 
                              "fit to NMR data from the 2014 KDHS for 2010-2014 ", 
                              tableText, "."), 
               label=paste0("tab:riskModTable_", resultType)))
  
  # same for CI width
  conStratDat = 100*(prevalenceCIWidthConStrat-prevalenceCIWidthConStratlcpb)/prevalenceCIWidthConStratlcpb
  constituencyDat = 100*(prevalenceCIWidthConstituency-prevalenceCIWidthConstituencylcpb)/prevalenceCIWidthConstituencylcpb
  countyDat = 100*(prevalenceCIWidthCounty-prevalenceCIWidthCountylcpb)/prevalenceCIWidthCountylcpb
  prevConStrat = mean(conStratDat)
  prevCon = mean(constituencyDat)
  prevCounty = mean(countyDat)
  prevMeans = c("conStrat"=prevConStrat, "constituency"=prevCon, "county"=prevCounty)
  
  conStratDat = 100*(countCIWidthConStrat-countCIWidthConStratlcpb)/countCIWidthConStratlcpb
  constituencyDat = 100*(countCIWidthConstituency-countCIWidthConstituencylcpb)/countCIWidthConstituencylcpb
  countyDat = 100*(countCIWidthCounty-countCIWidthCountylcpb)/countCIWidthCountylcpb
  burdConStrat = mean(conStratDat)
  burdCon = mean(constituencyDat)
  burdCounty = mean(countyDat)
  burdMeans = c("conStrat"=burdConStrat, "constituency"=burdCon, "county"=burdCounty)
  
  constituencyDat = 100*(relativePrevalenceCIWidthConstituency-relativePrevalenceCIWidthConstituencylcpb)/relativePrevalenceCIWidthConstituencylcpb
  countyDat = 100*(relativePrevalenceCIWidthCounty-relativePrevalenceCIWidthCountylcpb)/relativePrevalenceCIWidthCountylcpb
  relPrevCon = mean(constituencyDat, na.rm=TRUE)
  relPrevCounty = mean(countyDat, na.rm=TRUE)
  relPrevMeans = c("conStrat"=NA, "constituency"=relPrevCon, "county"=relPrevCounty)
  
  tab = cbind(prevMeans, burdMeans, relPrevMeans)
  colnames(tab) = c("Prevalence", "Burden", "Relative prevalence")
  row.names(tab) = c("Constituency x stratum", "Constituency", "County")
  tab = round(tab, 0)
  
  print(xtable(tab, digits=0, 
               caption=paste0("Summary statistics for risk model parameters when ", 
                              "fit to NMR data from the 2014 KDHS for 2010-2014 ", 
                              tableText, "."), 
               label=paste0("tab:riskModTable_", resultType)))
  
  # \begin{table}[ht]
  # \centering
  # \begin{tabular}{rrrr}
  # \hline
  # & Prevalence & Burden & Relative prevalence \\ 
  # \hline
  # Constituency x stratum & 31 & 45 &  \\ 
  # Constituency & 9 & 13 & 133 \\ 
  # County & 2 & 2 & 20 \\ 
  # \hline
  # \end{tabular}
  # \end{table}
  
  print("finished with plots")
  browser()
}

makeSensitivityPlots = function(logisticApproximation=FALSE, coarse=FALSE, signif=.95) {
  alpha = 1 - signif
  
  # first load the model predictions
  if(logisticApproximation) {
    stop("logisticApproximation not currently used")
  }
  logisticText = ifelse(!logisticApproximation, "", "logisticApprox")
  coarseText = ifelse(!coarse, "", "Coarse")
  
  resultTypes = c("std", "FBpop", "census2019", "censusJittered")
  for(resultType in resultTypes) {
    resultTypeText = paste0("_", resultType)
    out = load(paste0("savedOutput/application/finalMort", coarseText, logisticText, resultTypeText, ".RData"))
    
    if(resultType == "std") {
      easpa = makeDefaultEASPA()
      poppsub = poppsubKenyaThresh
      
      if(coarse) {
        popMat = popGridCoarseThresh
        popMatAdjusted = popGridCoarseAdjustedThresh
      } else {
        popMat = popGridThresh
        popMatAdjusted = popGridAdjustedThresh
      }
    } else if(resultType == "FBpop") {
      out = load("savedOutput/global/kenyaFacePopulationMats.RData")
      out = load("savedOutput/global/poppsubFace.RData")
      
      easpa = makeEASPAfacebook()
      poppsub = poppsubKenyaFaceThresh
      
      if(coarse) {
        popMat = popMatCoarseFaceThresh
        popMatAdjusted = popMatCoarseAdjustedFaceThresh
      } else {
        popMat = popMatKenyaFaceThresh
        popMatAdjusted = popMatKenyaFaceNeonatalThresh
      }
    } else if(resultType == "census2019") {
      load("savedOutput/global/kenya2019PopulationMats.RData")
      load("savedOutput/global/poppsub2019.RData")
      
      easpa = makeEASPA2019()
      poppsub = poppsubKenya2019Thresh
      
      if(coarse) {
        popMat = popMatCoarse2019Thresh
        popMatAdjusted = popMatCoarseAdjusted2019Thresh
      } else {
        popMat = popMatKenya2019Thresh
        popMatAdjusted = popMatKenya2019NeonatalThresh
      }
    } else if(resultType == "censusJittered") {
      load("savedOutput/global/kenyaJitteredPopulationMats.RData")
      load("savedOutput/global/poppsubJittered.RData")
      
      easpa = makeEASPAJittered()
      poppsub = poppsubKenyaJitteredThresh
      
      if(coarse) {
        popMat = popMatCoarseJitteredThresh
        popMatAdjusted = popMatCoarseAdjustedJitteredThresh
      } else {
        popMat = popMatKenyaJitteredThresh
        popMatAdjusted = popMatKenyaJitteredNeonatalThresh
      }
    }
    if(!is.null(poppsub$area)) {
      poppsub$County = poppsub$area
    }
    
    # calculate predictions and SDs ----
    
    areaLevels = c("pixel", "constrat", "constituency", "county") # province results no longer supported
    for(i in 1:length(areaLevels)) {
      thisLevel = areaLevels[i]
      
      if(thisLevel == "constrat") {
        undefinedPrevalenceConStrats = c(poppsub$popUrb == 0, poppsub$popRur == 0)
        
        prevalenceSDConStrat = c(apply(subareaPop$pUrbanFineScalePrevalence, 1, sd, na.rm=TRUE), 
                                 apply(subareaPop$pRuralFineScalePrevalence, 1, sd, na.rm=TRUE))
        # prevalenceSDConStrat = prevalenceSDConStrat[!undefinedPrevalenceConStrats]
        countSDConStrat = c(apply(subareaPop$ZUrbanFineScalePrevalence, 1, sd, na.rm=TRUE), 
                            apply(subareaPop$ZRuralFineScalePrevalence, 1, sd, na.rm=TRUE))
        # countSDConStrat = countSDConStrat[!undefinedPrevalenceConStrats]
        
        # get same for the lcpb model
        prevalenceSDConStratlcpb = c(apply(subareaPop$pUrbanSmoothRisk, 1, sd, na.rm=TRUE), 
                                     apply(subareaPop$pRuralSmoothRisk, 1, sd, na.rm=TRUE))
        # prevalenceSDConStratlcpb = prevalenceSDConStratlcpb[!undefinedPrevalenceConStrats]
        countSDConStratlcpb = c(apply(subareaPop$ZUrbanSmoothRisk, 1, sd, na.rm=TRUE), 
                                apply(subareaPop$ZRuralSmoothRisk, 1, sd, na.rm=TRUE))
        # countSDConStratlcpb = countSDConStratlcpb[!undefinedPrevalenceConStrats]
        
        # do the same for the LCPb model
        prevalenceSDConStratLCPb = c(apply(subareaPop$pUrbanFineScaleRisk, 1, sd, na.rm=TRUE), 
                                     apply(subareaPop$pRuralFineScaleRisk, 1, sd, na.rm=TRUE))
        # prevalenceSDConStratLCPb = prevalenceSDConStratLCPb[!undefinedPrevalenceConStrats]
        countSDConStratLCPb = c(apply(subareaPop$ZUrbanFineScaleRisk, 1, sd, na.rm=TRUE), 
                                apply(subareaPop$ZRuralFineScaleRisk, 1, sd, na.rm=TRUE))
        # countSDConStratLCPb = countSDConStratLCPb[!undefinedPrevalenceConStrats]
        
        # predictions: prevalence, burden, relPrev
        prevalencePredsConStrat = c(rowMeans(subareaPop$pUrbanFineScalePrevalence, na.rm=TRUE), 
                                    rowMeans(subareaPop$pRuralFineScalePrevalence, na.rm=TRUE))
        prevalencePredsConStratlcpb = c(rowMeans(subareaPop$pUrbanSmoothRisk, na.rm=TRUE), 
                                        rowMeans(subareaPop$pRuralSmoothRisk, na.rm=TRUE))
        prevalencePredsConStratLCPb = c(rowMeans(subareaPop$pUrbanFineScaleRisk, na.rm=TRUE), 
                                        rowMeans(subareaPop$pRuralFineScaleRisk, na.rm=TRUE))
        
        countPredsConStrat = c(rowMeans(subareaPop$ZUrbanFineScalePrevalence, na.rm=TRUE), 
                               rowMeans(subareaPop$ZRuralFineScalePrevalence, na.rm=TRUE))
        countPredsConStratlcpb = c(rowMeans(subareaPop$ZUrbanSmoothRisk, na.rm=TRUE), 
                                   rowMeans(subareaPop$ZRuralSmoothRisk, na.rm=TRUE))
        countPredsConStratLCPb = c(rowMeans(subareaPop$ZUrbanFineScaleRisk, na.rm=TRUE), 
                                   rowMeans(subareaPop$ZRuralFineScaleRisk, na.rm=TRUE))
      } else if(thisLevel == "constituency") {
        urbanConstituencies = poppsub$County == "Nairobi" | poppsub$County == "Mombasa"
        undefinedRelativePrevalenceConstituencies = (poppsub$popUrb == 0) | (poppsub$popRur == 0)
        relativePrevalencePredConstituency = rowMeans(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, na.rm=TRUE)
        undefinedRelativePrevalenceConstituencies = undefinedRelativePrevalenceConstituencies | !is.finite(relativePrevalencePredConstituency)
        relativePrevalencePredConstituency[undefinedRelativePrevalenceConstituencies] = NA
        
        prevalenceSDConstituency = apply(subareaPop$pFineScalePrevalence, 1, sd, na.rm=TRUE)
        countSDConstituency = apply(subareaPop$ZFineScalePrevalence, 1, sd, na.rm=TRUE)
        relativePrevalenceSDConstituency = apply(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, 1, sd, na.rm=TRUE)
        relativePrevalenceSDConstituency[undefinedRelativePrevalenceConstituencies] = NA
        
        # get SDs for the lcpb model
        prevalenceSDConstituencylcpb = apply(subareaPop$pSmoothRisk, 1, sd, na.rm=TRUE)
        countSDConstituencylcpb = apply(subareaPop$ZSmoothRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDConstituencylcpb = apply(subareaPop$pUrbanSmoothRisk/subareaPop$pRuralSmoothRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDConstituencylcpb[undefinedRelativePrevalenceConstituencies] = NA
        
        # do the same for the LCPb model
        prevalenceSDConstituencyLCPb = apply(subareaPop$pFineScaleRisk, 1, sd, na.rm=TRUE)
        countSDConstituencyLCPb = apply(subareaPop$ZFineScaleRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDConstituencyLCPb = apply(subareaPop$pUrbanFineScaleRisk/subareaPop$pRuralFineScaleRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDConstituencyLCPb[undefinedRelativePrevalenceConstituencies] = NA
        
        # predictions: prevalence, burden, relPrev
        prevalencePredsConstituency = rowMeans(subareaPop$pFineScalePrevalence, na.rm=TRUE)
        prevalencePredsConstituencylcpb = rowMeans(subareaPop$pSmoothRisk, na.rm=TRUE)
        prevalencePredsConstituencyLCPb = rowMeans(subareaPop$pFineScaleRisk, na.rm=TRUE)
        
        countPredsConstituency = rowMeans(subareaPop$ZFineScalePrevalence, na.rm=TRUE)
        countPredsConstituencylcpb = rowMeans(subareaPop$ZSmoothRisk, na.rm=TRUE)
        countPredsConstituencyLCPb = rowMeans(subareaPop$ZFineScaleRisk, na.rm=TRUE)
        
        relativePrevalencePredsConstituency = rowMeans(subareaPop$pUrbanFineScalePrevalence/subareaPop$pRuralFineScalePrevalence, na.rm=TRUE)
        relativePrevalencePredsConstituencylcpb = rowMeans(subareaPop$pUrbanSmoothRisk/subareaPop$pRuralSmoothRisk, na.rm=TRUE)
        relativePrevalencePredsConstituencyLCPb = rowMeans(subareaPop$pUrbanFineScaleRisk/subareaPop$pRuralFineScaleRisk, na.rm=TRUE)
        relativePrevalencePredsConstituency[undefinedRelativePrevalenceConstituencies] = NA
        relativePrevalencePredsConstituencylcpb[undefinedRelativePrevalenceConstituencies] = NA
        relativePrevalencePredsConstituencyLCPb[undefinedRelativePrevalenceConstituencies] = NA
      } else if(thisLevel == "county") {
        urbanCounties = sort(poppc$County) == "Nairobi" | sort(poppc$County) == "Mombasa"
        relativePrevalencePredCounty = rowMeans(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence, na.rm=TRUE)
        relativePrevalencePredCounty[urbanCounties] = NA
        
        prevalenceSDCounty = apply(areaPop$pFineScalePrevalence, 1, sd, na.rm=TRUE)
        countSDCounty = apply(areaPop$ZFineScalePrevalence, 1, sd, na.rm=TRUE)
        relativePrevalenceSDCounty = apply(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence, 1, sd, na.rm=TRUE)
        relativePrevalenceSDCounty[urbanCounties] = NA
        
        # get SDs for the lcpb model
        prevalenceSDCountylcpb = apply(areaPop$pSmoothRisk, 1, sd, na.rm=TRUE)
        countSDCountylcpb = apply(areaPop$ZSmoothRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDCountylcpb = apply(areaPop$pUrbanSmoothRisk/areaPop$pRuralSmoothRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDCountylcpb[urbanCounties] = NA
        
        # get SDs for the LCPb model
        prevalenceSDCountyLCPb = apply(areaPop$pFineScaleRisk, 1, sd, na.rm=TRUE)
        countSDCountyLCPb = apply(areaPop$ZFineScaleRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDCountyLCPb = apply(areaPop$pUrbanFineScaleRisk/areaPop$pRuralFineScaleRisk, 1, sd, na.rm=TRUE)
        relativePrevalenceSDCountyLCPb[urbanCounties] = NA
        
        # predictions: prevalence, burden, relPrev
        prevalencePredsCounty = rowMeans(areaPop$pFineScalePrevalence, na.rm=TRUE)
        prevalencePredsCountylcpb = rowMeans(areaPop$pSmoothRisk, na.rm=TRUE)
        prevalencePredsCountyLCPb = rowMeans(areaPop$pFineScaleRisk, na.rm=TRUE)
        
        countPredsCounty = rowMeans(areaPop$ZFineScalePrevalence, na.rm=TRUE)
        countPredsCountylcpb = rowMeans(areaPop$ZSmoothRisk, na.rm=TRUE)
        countPredsCountyLCPb = rowMeans(areaPop$ZFineScaleRisk, na.rm=TRUE)
        
        relativePrevalencePredsCounty = rowMeans(areaPop$pUrbanFineScalePrevalence/areaPop$pRuralFineScalePrevalence, na.rm=TRUE)
        relativePrevalencePredsCountylcpb = rowMeans(areaPop$pUrbanSmoothRisk/areaPop$pRuralSmoothRisk, na.rm=TRUE)
        relativePrevalencePredsCountyLCPb = rowMeans(areaPop$pUrbanFineScaleRisk/areaPop$pRuralFineScaleRisk, na.rm=TRUE)
        relativePrevalencePredsCounty[urbanCounties] = NA
        relativePrevalencePredsCountylcpb[urbanCounties] = NA
        relativePrevalencePredsCountyLCPb[urbanCounties] = NA
      }
    }
    
    if(resultType == "std") {
      urbanConstituencies_std = urbanConstituencies
      undefinedPrevalenceConStrats_std = undefinedPrevalenceConStrats
      prevalencePredsConStrat_std = prevalencePredsConStrat
      prevalenceSDConStrat_std = prevalenceSDConStrat
      countSDConStrat_std = countSDConStrat
      
      # get same for the lcpb model
      prevalencePredsConStratlcpb_std = prevalencePredsConStratlcpb
      prevalenceSDConStratlcpb_std = prevalenceSDConStratlcpb
      countSDConStratlcpb_std = countSDConStratlcpb
      
      # do the same for the LCPb model
      prevalencePredsConStratLCPb_std = prevalencePredsConStrat
      prevalenceSDConStratLCPb_std = prevalenceSDConStratLCPb
      countSDConStratLCPb_std = countSDConStratLCPb
      
      urbanConstituencies_std = urbanConstituencies
      undefinedRelativePrevalenceConstituencies_std = undefinedRelativePrevalenceConstituencies
      relativePrevalencePredConstituency_std = relativePrevalencePredConstituency
      
      prevalenceSDConstituency_std = prevalenceSDConstituency
      countSDConstituency_std = countSDConstituency
      relativePrevalenceSDConstituency_std = relativePrevalenceSDConstituency
      
      # get SDs for the lcpb model
      prevalenceSDConstituencylcpb_std = prevalenceSDConstituencylcpb
      countSDConstituencylcpb_std = countSDConstituencylcpb
      relativePrevalenceSDConstituencylcpb_std = relativePrevalenceSDConstituencylcpb
      
      # do the same for the LCPb model
      prevalenceSDConstituencyLCPb_std = prevalenceSDConstituencyLCPb
      countSDConstituencyLCPb_std = countSDConstituencyLCPb
      relativePrevalenceSDConstituencyLCPb_std = relativePrevalenceSDConstituencyLCPb
      
      urbanCounties_std = urbanCounties
      relativePrevalencePredCounty_std = relativePrevalencePredCounty
      
      prevalenceSDCounty_std = prevalenceSDCounty
      countSDCounty_std = countSDCounty
      relativePrevalenceSDCounty_std = relativePrevalenceSDCounty
      
      # get SDs for the lcpb model
      prevalenceSDCountylcpb_std = prevalenceSDCountylcpb
      countSDCountylcpb_std = countSDCountylcpb
      relativePrevalenceSDCountylcpb_std = relativePrevalenceSDCountylcpb
      
      # get SDs for the LCPb model
      prevalenceSDCountyLCPb_std = prevalenceSDCountyLCPb
      countSDCountyLCPb_std = countSDCountyLCPb
      relativePrevalenceSDCountyLCPb_std = relativePrevalenceSDCountyLCPb
      
      # predictions
      prevalencePredsConStrat_std = prevalencePredsConStrat
      prevalencePredsConStratlcpb_std = prevalencePredsConStratlcpb
      prevalencePredsConStratLCPb_std = prevalencePredsConStratLCPb
      
      countPredsConStrat_std = countPredsConStrat
      countPredsConStratlcpb_std = countPredsConStratlcpb
      countPredsConStratLCPb_std = countPredsConStratLCPb
      
      
      prevalencePredsConstituency_std = prevalencePredsConstituency
      prevalencePredsConstituencylcpb_std = prevalencePredsConstituencylcpb
      prevalencePredsConstituencyLCPb_std = prevalencePredsConstituencyLCPb
      
      countPredsConstituency_std = countPredsConstituency
      countPredsConstituencylcpb_std = countPredsConstituencylcpb
      countPredsConstituencyLCPb_std = countPredsConstituencyLCPb
      
      relativePrevalencePredsConstituency_std = relativePrevalencePredsConstituency
      relativePrevalencePredsConstituencylcpb_std = relativePrevalencePredsConstituencylcpb
      relativePrevalencePredsConstituencyLCPb_std = relativePrevalencePredsConstituencyLCPb
      
      
      prevalencePredsCounty_std = prevalencePredsCounty
      prevalencePredsCountylcpb_std = prevalencePredsCountylcpb
      prevalencePredsCountyLCPb_std = prevalencePredsCountyLCPb
      
      countPredsCounty_std = countPredsCounty
      countPredsCountylcpb_std = countPredsCountylcpb
      countPredsCountyLCPb_std = countPredsCountyLCPb
      
      relativePrevalencePredsCounty_std = relativePrevalencePredsCounty
      relativePrevalencePredsCountylcpb_std = relativePrevalencePredsCountylcpb
      relativePrevalencePredsCountyLCPb_std = relativePrevalencePredsCountyLCPb
    } else if(resultType == "FBpop") {
      urbanConstituencies_FBpop = urbanConstituencies
      undefinedPrevalenceConStrats_FBpop = undefinedPrevalenceConStrats
      
      prevalenceSDConStrat_FBpop = prevalenceSDConStrat
      countSDConStrat_FBpop = countSDConStrat
      
      # get SDss for the lcpb model
      prevalenceSDConStratlcpb_FBpop = prevalenceSDConStratlcpb
      countSDConStratlcpb_FBpop = countSDConStratlcpb
      
      # do the same for the LCPb model
      prevalenceSDConStratLCPb_FBpop = prevalenceSDConStratLCPb
      countSDConStratLCPb_FBpop = countSDConStratLCPb
      
      urbanConstituencies_FBpop = urbanConstituencies
      undefinedRelativePrevalenceConstituencies_FBpop = undefinedRelativePrevalenceConstituencies
      relativePrevalencePredConstituency_FBpop = relativePrevalencePredConstituency
      
      prevalenceSDConstituency_FBpop = prevalenceSDConstituency
      countSDConstituency_FBpop = countSDConstituency
      relativePrevalenceSDConstituency_FBpop = relativePrevalenceSDConstituency
      
      # get SDs for the lcpb model
      prevalenceSDConstituencylcpb_FBpop = prevalenceSDConstituencylcpb
      countSDConstituencylcpb_FBpop = countSDConstituencylcpb
      relativePrevalenceSDConstituencylcpb_FBpop = relativePrevalenceSDConstituencylcpb
      
      # do the same for the LCPb model
      prevalenceSDConstituencyLCPb_FBpop = prevalenceSDConstituencyLCPb
      countSDConstituencyLCPb_FBpop = countSDConstituencyLCPb
      relativePrevalenceSDConstituencyLCPb_FBpop = relativePrevalenceSDConstituencyLCPb
      
      urbanCounties_FBpop = urbanCounties
      relativePrevalencePredCounty_FBpop = relativePrevalencePredCounty
      
      prevalenceSDCounty_FBpop = prevalenceSDCounty
      countSDCounty_FBpop = countSDCounty
      relativePrevalenceSDCounty_FBpop = relativePrevalenceSDCounty
      
      # get SDs for the lcpb model
      prevalenceSDCountylcpb_FBpop = prevalenceSDCountylcpb
      countSDCountylcpb_FBpop = countSDCountylcpb
      relativePrevalenceSDCountylcpb_FBpop = relativePrevalenceSDCountylcpb
      
      # get SDs for the LCPb model
      prevalenceSDCountyLCPb_FBpop = prevalenceSDCountyLCPb
      countSDCountyLCPb_FBpop = countSDCountyLCPb
      relativePrevalenceSDCountyLCPb_FBpop = relativePrevalenceSDCountyLCPb
      
      # predictions
      prevalencePredsConStrat_FBpop = prevalencePredsConStrat
      prevalencePredsConStratlcpb_FBpop = prevalencePredsConStratlcpb
      prevalencePredsConStratLCPb_FBpop = prevalencePredsConStratLCPb
      
      countPredsConStrat_FBpop = countPredsConStrat
      countPredsConStratlcpb_FBpop = countPredsConStratlcpb
      countPredsConStratLCPb_FBpop = countPredsConStratLCPb
      
      
      prevalencePredsConstituency_FBpop = prevalencePredsConstituency
      prevalencePredsConstituencylcpb_FBpop = prevalencePredsConstituencylcpb
      prevalencePredsConstituencyLCPb_FBpop = prevalencePredsConstituencyLCPb
      
      countPredsConstituency_FBpop = countPredsConstituency
      countPredsConstituencylcpb_FBpop = countPredsConstituencylcpb
      countPredsConstituencyLCPb_FBpop = countPredsConstituencyLCPb
      
      relativePrevalencePredsConstituency_FBpop = relativePrevalencePredsConstituency
      relativePrevalencePredsConstituencylcpb_FBpop = relativePrevalencePredsConstituencylcpb
      relativePrevalencePredsConstituencyLCPb_FBpop = relativePrevalencePredsConstituencyLCPb
      
      
      prevalencePredsCounty_FBpop = prevalencePredsCounty
      prevalencePredsCountylcpb_FBpop = prevalencePredsCountylcpb
      prevalencePredsCountyLCPb_FBpop = prevalencePredsCountyLCPb
      
      countPredsCounty_FBpop = countPredsCounty
      countPredsCountylcpb_FBpop = countPredsCountylcpb
      countPredsCountyLCPb_FBpop = countPredsCountyLCPb
      
      relativePrevalencePredsCounty_FBpop = relativePrevalencePredsCounty
      relativePrevalencePredsCountylcpb_FBpop = relativePrevalencePredsCountylcpb
      relativePrevalencePredsCountyLCPb_FBpop = relativePrevalencePredsCountyLCPb
    }  else if(resultType == "census2019") {
      urbanConstituencies_2019 = urbanConstituencies
      undefinedPrevalenceConStrats_2019 = undefinedPrevalenceConStrats
      
      prevalenceSDConStrat_2019 = prevalenceSDConStrat
      countSDConStrat_2019 = countSDConStrat
      
      # get SDss for the lcpb model
      prevalenceSDConStratlcpb_2019 = prevalenceSDConStratlcpb
      countSDConStratlcpb_2019 = countSDConStratlcpb
      
      # do the same for the LCPb model
      prevalenceSDConStratLCPb_2019 = prevalenceSDConStratLCPb
      countSDConStratLCPb_2019 = countSDConStratLCPb
      
      urbanConstituencies_2019 = urbanConstituencies
      undefinedRelativePrevalenceConstituencies_2019 = undefinedRelativePrevalenceConstituencies
      relativePrevalencePredConstituency_2019 = relativePrevalencePredConstituency
      
      prevalenceSDConstituency_2019 = prevalenceSDConstituency
      countSDConstituency_2019 = countSDConstituency
      relativePrevalenceSDConstituency_2019 = relativePrevalenceSDConstituency
      
      # get SDs for the lcpb model
      prevalenceSDConstituencylcpb_2019 = prevalenceSDConstituencylcpb
      countSDConstituencylcpb_2019 = countSDConstituencylcpb
      relativePrevalenceSDConstituencylcpb_2019 = relativePrevalenceSDConstituencylcpb
      
      # do the same for the LCPb model
      prevalenceSDConstituencyLCPb_2019 = prevalenceSDConstituencyLCPb
      countSDConstituencyLCPb_2019 = countSDConstituencyLCPb
      relativePrevalenceSDConstituencyLCPb_2019 = relativePrevalenceSDConstituencyLCPb
      
      urbanCounties_2019 = urbanCounties
      relativePrevalencePredCounty_2019 = relativePrevalencePredCounty
      
      prevalenceSDCounty_2019 = prevalenceSDCounty
      countSDCounty_2019 = countSDCounty
      relativePrevalenceSDCounty_2019 = relativePrevalenceSDCounty
      
      # get SDs for the lcpb model
      prevalenceSDCountylcpb_2019 = prevalenceSDCountylcpb
      countSDCountylcpb_2019 = countSDCountylcpb
      relativePrevalenceSDCountylcpb_2019 = relativePrevalenceSDCountylcpb
      
      # get SDs for the LCPb model
      prevalenceSDCountyLCPb_2019 = prevalenceSDCountyLCPb
      countSDCountyLCPb_2019 = countSDCountyLCPb
      relativePrevalenceSDCountyLCPb_2019 = relativePrevalenceSDCountyLCPb
      
      # predictions
      prevalencePredsConStrat_2019 = prevalencePredsConStrat
      prevalencePredsConStratlcpb_2019 = prevalencePredsConStratlcpb
      prevalencePredsConStratLCPb_2019 = prevalencePredsConStratLCPb
      
      countPredsConStrat_2019 = countPredsConStrat
      countPredsConStratlcpb_2019 = countPredsConStratlcpb
      countPredsConStratLCPb_2019 = countPredsConStratLCPb
      
      
      prevalencePredsConstituency_2019 = prevalencePredsConstituency
      prevalencePredsConstituencylcpb_2019 = prevalencePredsConstituencylcpb
      prevalencePredsConstituencyLCPb_2019 = prevalencePredsConstituencyLCPb
      
      countPredsConstituency_2019 = countPredsConstituency
      countPredsConstituencylcpb_2019 = countPredsConstituencylcpb
      countPredsConstituencyLCPb_2019 = countPredsConstituencyLCPb
      
      relativePrevalencePredsConstituency_2019 = relativePrevalencePredsConstituency
      relativePrevalencePredsConstituencylcpb_2019 = relativePrevalencePredsConstituencylcpb
      relativePrevalencePredsConstituencyLCPb_2019 = relativePrevalencePredsConstituencyLCPb
      
      
      prevalencePredsCounty_2019 = prevalencePredsCounty
      prevalencePredsCountylcpb_2019 = prevalencePredsCountylcpb
      prevalencePredsCountyLCPb_2019 = prevalencePredsCountyLCPb
      
      countPredsCounty_2019 = countPredsCounty
      countPredsCountylcpb_2019 = countPredsCountylcpb
      countPredsCountyLCPb_2019 = countPredsCountyLCPb
      
      relativePrevalencePredsCounty_2019 = relativePrevalencePredsCounty
      relativePrevalencePredsCountylcpb_2019 = relativePrevalencePredsCountylcpb
      relativePrevalencePredsCountyLCPb_2019 = relativePrevalencePredsCountyLCPb
    } else if(resultType == "censusJittered") {
      urbanConstituencies_jittered = urbanConstituencies
      undefinedPrevalenceConStrats_jittered = undefinedPrevalenceConStrats
      
      prevalenceSDConStrat_jittered = prevalenceSDConStrat
      countSDConStrat_jittered = countSDConStrat
      
      # get SDss for the lcpb model
      prevalenceSDConStratlcpb_jittered = prevalenceSDConStratlcpb
      countSDConStratlcpb_jittered = countSDConStratlcpb
      
      # do the same for the LCPb model
      prevalenceSDConStratLCPb_jittered = prevalenceSDConStratLCPb
      countSDConStratLCPb_jittered = countSDConStratLCPb
      
      urbanConstituencies_jittered = urbanConstituencies
      undefinedRelativePrevalenceConstituencies_jittered = undefinedRelativePrevalenceConstituencies
      relativePrevalencePredConstituency_jittered = relativePrevalencePredConstituency
      
      prevalenceSDConstituency_jittered = prevalenceSDConstituency
      countSDConstituency_jittered = countSDConstituency
      relativePrevalenceSDConstituency_jittered = relativePrevalenceSDConstituency
      
      # get SDs for the lcpb model
      prevalenceSDConstituencylcpb_jittered = prevalenceSDConstituencylcpb
      countSDConstituencylcpb_jittered = countSDConstituencylcpb
      relativePrevalenceSDConstituencylcpb_jittered = relativePrevalenceSDConstituencylcpb
      
      # do the same for the LCPb model
      prevalenceSDConstituencyLCPb_jittered = prevalenceSDConstituencyLCPb
      countSDConstituencyLCPb_jittered = countSDConstituencyLCPb
      relativePrevalenceSDConstituencyLCPb_jittered = relativePrevalenceSDConstituencyLCPb
      
      urbanCounties_jittered = urbanCounties
      relativePrevalencePredCounty_jittered = relativePrevalencePredCounty
      
      prevalenceSDCounty_jittered = prevalenceSDCounty
      countSDCounty_jittered = countSDCounty
      relativePrevalenceSDCounty_jittered = relativePrevalenceSDCounty
      
      # get SDs for the lcpb model
      prevalenceSDCountylcpb_jittered = prevalenceSDCountylcpb
      countSDCountylcpb_jittered = countSDCountylcpb
      relativePrevalenceSDCountylcpb_jittered = relativePrevalenceSDCountylcpb
      
      # get SDs for the LCPb model
      prevalenceSDCountyLCPb_jittered = prevalenceSDCountyLCPb
      countSDCountyLCPb_jittered = countSDCountyLCPb
      relativePrevalenceSDCountyLCPb_jittered = relativePrevalenceSDCountyLCPb
      
      # predictions
      prevalencePredsConStrat_jittered = prevalencePredsConStrat
      prevalencePredsConStratlcpb_jittered = prevalencePredsConStratlcpb
      prevalencePredsConStratLCPb_jittered = prevalencePredsConStratLCPb
      
      countPredsConStrat_jittered = countPredsConStrat
      countPredsConStratlcpb_jittered = countPredsConStratlcpb
      countPredsConStratLCPb_jittered = countPredsConStratLCPb
      
      
      prevalencePredsConstituency_jittered = prevalencePredsConstituency
      prevalencePredsConstituencylcpb_jittered = prevalencePredsConstituencylcpb
      prevalencePredsConstituencyLCPb_jittered = prevalencePredsConstituencyLCPb
      
      countPredsConstituency_jittered = countPredsConstituency
      countPredsConstituencylcpb_jittered = countPredsConstituencylcpb
      countPredsConstituencyLCPb_jittered = countPredsConstituencyLCPb
      
      relativePrevalencePredsConstituency_jittered = relativePrevalencePredsConstituency
      relativePrevalencePredsConstituencylcpb_jittered = relativePrevalencePredsConstituencylcpb
      relativePrevalencePredsConstituencyLCPb_jittered = relativePrevalencePredsConstituencyLCPb
      
      
      prevalencePredsCounty_jittered = prevalencePredsCounty
      prevalencePredsCountylcpb_jittered = prevalencePredsCountylcpb
      prevalencePredsCountyLCPb_jittered = prevalencePredsCountyLCPb
      
      countPredsCounty_jittered = countPredsCounty
      countPredsCountylcpb_jittered = countPredsCountylcpb
      countPredsCountyLCPb_jittered = countPredsCountyLCPb
      
      relativePrevalencePredsCounty_jittered = relativePrevalencePredsCounty
      relativePrevalencePredsCountylcpb_jittered = relativePrevalencePredsCountylcpb
      relativePrevalencePredsCountyLCPb_jittered = relativePrevalencePredsCountyLCPb
    }
  }
  
  # Combine data into tables ----
  allPrevalenceSD_std = c(prevalenceSDConStrat_std, prevalenceSDConstituency_std, prevalenceSDCounty_std)
  allPrevalenceSDlcpb_std = c(prevalenceSDConStratlcpb_std, prevalenceSDConstituencylcpb_std, prevalenceSDCountylcpb_std)
  allPrevalenceSDLCPb_std = c(prevalenceSDConStratLCPb_std, prevalenceSDConstituencyLCPb_std, prevalenceSDCountyLCPb_std)
  
  allPrevalenceSD_FBpop = c(prevalenceSDConStrat_FBpop, prevalenceSDConstituency_FBpop, prevalenceSDCounty_FBpop)
  allPrevalenceSDlcpb_FBpop = c(prevalenceSDConStratlcpb_FBpop, prevalenceSDConstituencylcpb_FBpop, prevalenceSDCountylcpb_FBpop)
  allPrevalenceSDLCPb_FBpop = c(prevalenceSDConStratLCPb_FBpop, prevalenceSDConstituencyLCPb_FBpop, prevalenceSDCountyLCPb_FBpop)
  
  allPrevalenceSD_2019 = c(prevalenceSDConStrat_2019, prevalenceSDConstituency_2019, prevalenceSDCounty_2019)
  allPrevalenceSDlcpb_2019 = c(prevalenceSDConStratlcpb_2019, prevalenceSDConstituencylcpb_2019, prevalenceSDCountylcpb_2019)
  allPrevalenceSDLCPb_2019 = c(prevalenceSDConStratLCPb_2019, prevalenceSDConstituencyLCPb_2019, prevalenceSDCountyLCPb_2019)
  
  allPrevalenceSD_jittered = c(prevalenceSDConStrat_jittered, prevalenceSDConstituency_jittered, prevalenceSDCounty_jittered)
  allPrevalenceSDlcpb_jittered = c(prevalenceSDConStratlcpb_jittered, prevalenceSDConstituencylcpb_jittered, prevalenceSDCountylcpb_jittered)
  allPrevalenceSDLCPb_jittered = c(prevalenceSDConStratLCPb_jittered, prevalenceSDConstituencyLCPb_jittered, prevalenceSDCountyLCPb_jittered)
  
  
  allPrevalencePreds_std = c(prevalencePredsConStrat_std, prevalencePredsConstituency_std, prevalencePredsCounty_std)
  allPrevalencePredslcpb_std = c(prevalencePredsConStratlcpb_std, prevalencePredsConstituencylcpb_std, prevalencePredsCountylcpb_std)
  allPrevalencePredsLCPb_std = c(prevalencePredsConStratLCPb_std, prevalencePredsConstituencyLCPb_std, prevalencePredsCountyLCPb_std)
  
  allPrevalencePreds_FBpop = c(prevalencePredsConStrat_FBpop, prevalencePredsConstituency_FBpop, prevalencePredsCounty_FBpop)
  allPrevalencePredslcpb_FBpop = c(prevalencePredsConStratlcpb_FBpop, prevalencePredsConstituencylcpb_FBpop, prevalencePredsCountylcpb_FBpop)
  allPrevalencePredsLCPb_FBpop = c(prevalencePredsConStratLCPb_FBpop, prevalencePredsConstituencyLCPb_FBpop, prevalencePredsCountyLCPb_FBpop)
  
  allPrevalencePreds_2019 = c(prevalencePredsConStrat_2019, prevalencePredsConstituency_2019, prevalencePredsCounty_2019)
  allPrevalencePredslcpb_2019 = c(prevalencePredsConStratlcpb_2019, prevalencePredsConstituencylcpb_2019, prevalencePredsCountylcpb_2019)
  allPrevalencePredsLCPb_2019 = c(prevalencePredsConStratLCPb_2019, prevalencePredsConstituencyLCPb_2019, prevalencePredsCountyLCPb_2019)
  
  allPrevalencePreds_jittered = c(prevalencePredsConStrat_jittered, prevalencePredsConstituency_jittered, prevalencePredsCounty_jittered)
  allPrevalencePredslcpb_jittered = c(prevalencePredsConStratlcpb_jittered, prevalencePredsConstituencylcpb_jittered, prevalencePredsCountylcpb_jittered)
  allPrevalencePredsLCPb_jittered = c(prevalencePredsConStratLCPb_jittered, prevalencePredsConstituencyLCPb_jittered, prevalencePredsCountyLCPb_jittered)
  
  areaLevelPrevalence = c(rep("Constituency x stratum", length(prevalencePredsConStrat_std)), 
                          rep("Constituency", length(prevalencePredsConstituency_std)), 
                          rep("County", length(prevalencePredsCounty_std)))
  
  fullPrevalenceTab = data.frame(areaLevel=areaLevelPrevalence, 
                                 allPrevalencePreds_std, allPrevalencePredslcpb_std, allPrevalencePredsLCPb_std, 
                                 allPrevalencePreds_FBpop, allPrevalencePredslcpb_FBpop, allPrevalencePredsLCPb_FBpop, 
                                 allPrevalencePreds_2019, allPrevalencePredslcpb_2019, allPrevalencePredsLCPb_2019, 
                                 allPrevalencePreds_jittered, allPrevalencePredslcpb_jittered, allPrevalencePredsLCPb_jittered)
  
  fullPrevalenceTabSD = data.frame(areaLevel=areaLevelPrevalence, 
                                 allPrevalenceSD_std, allPrevalenceSDlcpb_std, allPrevalenceSDLCPb_std, 
                                 allPrevalenceSD_FBpop, allPrevalenceSDlcpb_FBpop, allPrevalenceSDLCPb_FBpop, 
                                 allPrevalenceSD_2019, allPrevalenceSDlcpb_2019, allPrevalenceSDLCPb_2019, 
                                 allPrevalenceSD_jittered, allPrevalenceSDlcpb_jittered, allPrevalenceSDLCPb_jittered)
  
  fullPrevalenceTabPct = data.frame(areaLevel=areaLevelPrevalence, 
                                 empiricial_FBpop=100*(allPrevalencePreds_FBpop-allPrevalencePreds_std)/allPrevalencePreds_std, 
                                   latent_FBpop=100*(allPrevalencePredslcpb_FBpop-allPrevalencePredslcpb_std)/allPrevalencePredslcpb_std, 
                                   smoothLatent_FBpop=100*(allPrevalencePredsLCPb_FBpop-allPrevalencePredsLCPb_std)/allPrevalencePredsLCPb_std, 
                                 empiricial_2019=100*(allPrevalencePreds_2019-allPrevalencePreds_std)/allPrevalencePreds_std, 
                                   latent_2019=100*(allPrevalencePredslcpb_2019-allPrevalencePredslcpb_std)/allPrevalencePredslcpb_std, 
                                   smoothLatent_2019=100*(allPrevalencePredsLCPb_2019-allPrevalencePredsLCPb_std)/allPrevalencePredsLCPb_std, 
                                 empiricial_jittered=100*(allPrevalencePreds_jittered-allPrevalencePreds_std)/allPrevalencePreds_std, 
                                   latent_jittered=100*(allPrevalencePredslcpb_jittered-allPrevalencePredslcpb_std)/allPrevalencePredslcpb_std, 
                                   smoothLatent_jittered=100*(allPrevalencePredsLCPb_jittered-allPrevalencePredsLCPb_std)/allPrevalencePredsLCPb_std)
  
  fullPrevalenceTabSDPct = data.frame(areaLevel=areaLevelPrevalence, 
                                    empiricial_FBpop=100*(allPrevalenceSD_FBpop-allPrevalenceSD_std)/allPrevalenceSD_std, 
                                      latent_FBpop=100*(allPrevalenceSDlcpb_FBpop-allPrevalenceSDlcpb_std)/allPrevalenceSDlcpb_std, 
                                      smoothLatent_FBpop=100*(allPrevalenceSDLCPb_FBpop-allPrevalenceSDLCPb_std)/allPrevalenceSDLCPb_std, 
                                    empiricial_2019=100*(allPrevalenceSD_2019-allPrevalenceSD_std)/allPrevalenceSD_std, 
                                      latent_2019=100*(allPrevalenceSDlcpb_2019-allPrevalenceSDlcpb_std)/allPrevalenceSDlcpb_std, 
                                      smoothLatent_2019=100*(allPrevalenceSDLCPb_2019-allPrevalenceSDLCPb_std)/allPrevalenceSDLCPb_std, 
                                    empiricial_jittered=100*(allPrevalenceSD_jittered-allPrevalenceSD_std)/allPrevalenceSD_std, 
                                      latent_jittered=100*(allPrevalenceSDlcpb_jittered-allPrevalenceSDlcpb_std)/allPrevalenceSDlcpb_std, 
                                      smoothLatent_jittered=100*(allPrevalenceSDLCPb_jittered-allPrevalenceSDLCPb_std)/allPrevalenceSDLCPb_std)
  
  allCountSD_std = c(countSDConStrat_std, countSDConstituency_std, countSDCounty_std)
  allCountSDlcpb_std = c(countSDConStratlcpb_std, countSDConstituencylcpb_std, countSDCountylcpb_std)
  allCountSDLCPb_std = c(countSDConStratLCPb_std, countSDConstituencyLCPb_std, countSDCountyLCPb_std)
  
  allCountSD_FBpop = c(countSDConStrat_FBpop, countSDConstituency_FBpop, countSDCounty_FBpop)
  allCountSDlcpb_FBpop = c(countSDConStratlcpb_FBpop, countSDConstituencylcpb_FBpop, countSDCountylcpb_FBpop)
  allCountSDLCPb_FBpop = c(countSDConStratLCPb_FBpop, countSDConstituencyLCPb_FBpop, countSDCountyLCPb_FBpop)
  
  allCountSD_2019 = c(countSDConStrat_2019, countSDConstituency_2019, countSDCounty_2019)
  allCountSDlcpb_2019 = c(countSDConStratlcpb_2019, countSDConstituencylcpb_2019, countSDCountylcpb_2019)
  allCountSDLCPb_2019 = c(countSDConStratLCPb_2019, countSDConstituencyLCPb_2019, countSDCountyLCPb_2019)
  
  allCountSD_jittered = c(countSDConStrat_jittered, countSDConstituency_jittered, countSDCounty_jittered)
  allCountSDlcpb_jittered = c(countSDConStratlcpb_jittered, countSDConstituencylcpb_jittered, countSDCountylcpb_jittered)
  allCountSDLCPb_jittered = c(countSDConStratLCPb_jittered, countSDConstituencyLCPb_jittered, countSDCountyLCPb_jittered)
  
  
  allCountPreds_std = c(countPredsConStrat_std, countPredsConstituency_std, countPredsCounty_std)
  allCountPredslcpb_std = c(countPredsConStratlcpb_std, countPredsConstituencylcpb_std, countPredsCountylcpb_std)
  allCountPredsLCPb_std = c(countPredsConStratLCPb_std, countPredsConstituencyLCPb_std, countPredsCountyLCPb_std)
  
  allCountPreds_FBpop = c(countPredsConStrat_FBpop, countPredsConstituency_FBpop, countPredsCounty_FBpop)
  allCountPredslcpb_FBpop = c(countPredsConStratlcpb_FBpop, countPredsConstituencylcpb_FBpop, countPredsCountylcpb_FBpop)
  allCountPredsLCPb_FBpop = c(countPredsConStratLCPb_FBpop, countPredsConstituencyLCPb_FBpop, countPredsCountyLCPb_FBpop)
  
  allCountPreds_2019 = c(countPredsConStrat_2019, countPredsConstituency_2019, countPredsCounty_2019)
  allCountPredslcpb_2019 = c(countPredsConStratlcpb_2019, countPredsConstituencylcpb_2019, countPredsCountylcpb_2019)
  allCountPredsLCPb_2019 = c(countPredsConStratLCPb_2019, countPredsConstituencyLCPb_2019, countPredsCountyLCPb_2019)
  
  allCountPreds_jittered = c(countPredsConStrat_jittered, countPredsConstituency_jittered, countPredsCounty_jittered)
  allCountPredslcpb_jittered = c(countPredsConStratlcpb_jittered, countPredsConstituencylcpb_jittered, countPredsCountylcpb_jittered)
  allCountPredsLCPb_jittered = c(countPredsConStratLCPb_jittered, countPredsConstituencyLCPb_jittered, countPredsCountyLCPb_jittered)
  
  areaLevelCount = c(rep("Constituency x stratum", length(countPredsConStrat_std)), 
                          rep("Constituency", length(countPredsConstituency_std)), 
                          rep("County", length(countPredsCounty_std)))
  
  fullCountTab = data.frame(areaLevel=areaLevelCount, 
                                 allCountPreds_std, allCountPredslcpb_std, allCountPredsLCPb_std, 
                                 allCountPreds_FBpop, allCountPredslcpb_FBpop, allCountPredsLCPb_FBpop, 
                                 allCountPreds_2019, allCountPredslcpb_2019, allCountPredsLCPb_2019, 
                                 allCountPreds_jittered, allCountPredslcpb_jittered, allCountPredsLCPb_jittered)
  
  fullCountTabSD = data.frame(areaLevel=areaLevelCount, 
                                   allCountSD_std, allCountSDlcpb_std, allCountSDLCPb_std, 
                                   allCountSD_FBpop, allCountSDlcpb_FBpop, allCountSDLCPb_FBpop, 
                                   allCountSD_2019, allCountSDlcpb_2019, allCountSDLCPb_2019, 
                                   allCountSD_jittered, allCountSDlcpb_jittered, allCountSDLCPb_jittered)
  
  fullCountTabPct = data.frame(areaLevel=areaLevelCount, 
                                    empiricial_FBpop=100*(allCountPreds_FBpop-allCountPreds_std)/allCountPreds_std, 
                                      latent_FBpop=100*(allCountPredslcpb_FBpop-allCountPredslcpb_std)/allCountPredslcpb_std, 
                                      smoothLatent_FBpop=100*(allCountPredsLCPb_FBpop-allCountPredsLCPb_std)/allCountPredsLCPb_std, 
                                    empiricial_2019=100*(allCountPreds_2019-allCountPreds_std)/allCountPreds_std, 
                                      latent_2019=100*(allCountPredslcpb_2019-allCountPredslcpb_std)/allCountPredslcpb_std, 
                                      smoothLatent_2019=100*(allCountPredsLCPb_2019-allCountPredsLCPb_std)/allCountPredsLCPb_std, 
                                    empiricial_jittered=100*(allCountPreds_jittered-allCountPreds_std)/allCountPreds_std, 
                                      latent_jittered=100*(allCountPredslcpb_jittered-allCountPredslcpb_std)/allCountPredslcpb_std, 
                                      smoothLatent_jittered=100*(allCountPredsLCPb_jittered-allCountPredsLCPb_std)/allCountPredsLCPb_std)
  
  fullCountTabSDPct = data.frame(areaLevel=areaLevelCount, 
                                      empiricial_FBpop=100*(allCountSD_FBpop-allCountSD_std)/allCountSD_std, 
                                        latent_FBpop=100*(allCountSDlcpb_FBpop-allCountSDlcpb_std)/allCountSDlcpb_std, 
                                        smoothLatent_FBpop=100*(allCountSDLCPb_FBpop-allCountSDLCPb_std)/allCountSDLCPb_std, 
                                      empiricial_2019=100*(allCountSD_2019-allCountSD_std)/allCountSD_std, 
                                        latent_2019=100*(allCountSDlcpb_2019-allCountSDlcpb_std)/allCountSDlcpb_std, 
                                        smoothLatent_2019=100*(allCountSDLCPb_2019-allCountSDLCPb_std)/allCountSDLCPb_std, 
                                      empiricial_jittered=100*(allCountSD_jittered-allCountSD_std)/allCountSD_std, 
                                        latent_jittered=100*(allCountSDlcpb_jittered-allCountSDlcpb_std)/allCountSDlcpb_std, 
                                        smoothLatent_jittered=100*(allCountSDLCPb_jittered-allCountSDLCPb_std)/allCountSDLCPb_std)
  
  allRelPrevSD_std = c(relativePrevalenceSDConstituency_std, relativePrevalenceSDCounty_std)
  allRelPrevSDlcpb_std = c(relativePrevalenceSDConstituencylcpb_std, relativePrevalenceSDCountylcpb_std)
  allRelPrevSDLCPb_std = c(relativePrevalenceSDConstituencyLCPb_std, relativePrevalenceSDCountyLCPb_std)
  
  allRelPrevSD_FBpop = c(relativePrevalenceSDConstituency_FBpop, relativePrevalenceSDCounty_FBpop)
  allRelPrevSDlcpb_FBpop = c(relativePrevalenceSDConstituencylcpb_FBpop, relativePrevalenceSDCountylcpb_FBpop)
  allRelPrevSDLCPb_FBpop = c(relativePrevalenceSDConstituencyLCPb_FBpop, relativePrevalenceSDCountyLCPb_FBpop)
  
  allRelPrevSD_2019 = c(relativePrevalenceSDConstituency_2019, relativePrevalenceSDCounty_2019)
  allRelPrevSDlcpb_2019 = c(relativePrevalenceSDConstituencylcpb_2019, relativePrevalenceSDCountylcpb_2019)
  allRelPrevSDLCPb_2019 = c(relativePrevalenceSDConstituencyLCPb_2019, relativePrevalenceSDCountyLCPb_2019)
  
  allRelPrevSD_jittered = c(relativePrevalenceSDConstituency_jittered, relativePrevalenceSDCounty_jittered)
  allRelPrevSDlcpb_jittered = c(relativePrevalenceSDConstituencylcpb_jittered, relativePrevalenceSDCountylcpb_jittered)
  allRelPrevSDLCPb_jittered = c(relativePrevalenceSDConstituencyLCPb_jittered, relativePrevalenceSDCountyLCPb_jittered)
  
  
  allRelPrevPreds_std = c(relativePrevalencePredsConstituency_std, relativePrevalencePredsCounty_std)
  allRelPrevPredslcpb_std = c(relativePrevalencePredsConstituencylcpb_std, relativePrevalencePredsCountylcpb_std)
  allRelPrevPredsLCPb_std = c(relativePrevalencePredsConstituencyLCPb_std, relativePrevalencePredsCountyLCPb_std)
  
  allRelPrevPreds_FBpop = c(relativePrevalencePredsConstituency_FBpop, relativePrevalencePredsCounty_FBpop)
  allRelPrevPredslcpb_FBpop = c(relativePrevalencePredsConstituencylcpb_FBpop, relativePrevalencePredsCountylcpb_FBpop)
  allRelPrevPredsLCPb_FBpop = c(relativePrevalencePredsConstituencyLCPb_FBpop, relativePrevalencePredsCountyLCPb_FBpop)
  
  allRelPrevPreds_2019 = c(relativePrevalencePredsConstituency_2019, relativePrevalencePredsCounty_2019)
  allRelPrevPredslcpb_2019 = c(relativePrevalencePredsConstituencylcpb_2019, relativePrevalencePredsCountylcpb_2019)
  allRelPrevPredsLCPb_2019 = c(relativePrevalencePredsConstituencyLCPb_2019, relativePrevalencePredsCountyLCPb_2019)
  
  allRelPrevPreds_jittered = c(relativePrevalencePredsConstituency_jittered, relativePrevalencePredsCounty_jittered)
  allRelPrevPredslcpb_jittered = c(relativePrevalencePredsConstituencylcpb_jittered, relativePrevalencePredsCountylcpb_jittered)
  allRelPrevPredsLCPb_jittered = c(relativePrevalencePredsConstituencyLCPb_jittered, relativePrevalencePredsCountyLCPb_jittered)
  
  areaLevelRelPrev = c(rep("Constituency", length(relativePrevalencePredsConstituency_std)), 
                       rep("County", length(relativePrevalencePredsCounty_std)))
  
  fullRelPrevTab = data.frame(areaLevel=areaLevelRelPrev, 
                              allRelPrevPreds_std, allRelPrevPredslcpb_std, allRelPrevPredsLCPb_std, 
                              allRelPrevPreds_FBpop, allRelPrevPredslcpb_FBpop, allRelPrevPredsLCPb_FBpop, 
                              allRelPrevPreds_2019, allRelPrevPredslcpb_2019, allRelPrevPredsLCPb_2019, 
                              allRelPrevPreds_jittered, allRelPrevPredslcpb_jittered, allRelPrevPredsLCPb_jittered)
  
  fullRelPrevTabSD = data.frame(areaLevel=areaLevelRelPrev, 
                                   allRelPrevSD_std, allRelPrevSDlcpb_std, allRelPrevSDLCPb_std, 
                                   allRelPrevSD_FBpop, allRelPrevSDlcpb_FBpop, allRelPrevSDLCPb_FBpop, 
                                   allRelPrevSD_2019, allRelPrevSDlcpb_2019, allRelPrevSDLCPb_2019, 
                                   allRelPrevSD_jittered, allRelPrevSDlcpb_jittered, allRelPrevSDLCPb_jittered)
  
  fullRelPrevTabPct = data.frame(areaLevel=areaLevelRelPrev, 
                                    empiricial_FBpop=100*(allRelPrevPreds_FBpop-allRelPrevPreds_std)/allRelPrevPreds_std, 
                                      latent_FBpop=100*(allRelPrevPredslcpb_FBpop-allRelPrevPredslcpb_std)/allRelPrevPredslcpb_std, 
                                      smoothLatent_FBpop=100*(allRelPrevPredsLCPb_FBpop-allRelPrevPredsLCPb_std)/allRelPrevPredsLCPb_std, 
                                    empiricial_2019=100*(allRelPrevPreds_2019-allRelPrevPreds_std)/allRelPrevPreds_std, 
                                      latent_2019=100*(allRelPrevPredslcpb_2019-allRelPrevPredslcpb_std)/allRelPrevPredslcpb_std, 
                                      smoothLatent_2019=100*(allRelPrevPredsLCPb_2019-allRelPrevPredsLCPb_std)/allRelPrevPredsLCPb_std, 
                                    empiricial_jittered=100*(allRelPrevPreds_jittered-allRelPrevPreds_std)/allRelPrevPreds_std, 
                                      latent_jittered=100*(allRelPrevPredslcpb_jittered-allRelPrevPredslcpb_std)/allRelPrevPredslcpb_std, 
                                      smoothLatent_jittered=100*(allRelPrevPredsLCPb_jittered-allRelPrevPredsLCPb_std)/allRelPrevPredsLCPb_std)
  
  fullRelPrevTabSDPct = data.frame(areaLevel=areaLevelRelPrev, 
                                      empiricial_FBpop=100*(allRelPrevSD_FBpop-allRelPrevSD_std)/allRelPrevSD_std, 
                                        latent_FBpop=100*(allRelPrevSDlcpb_FBpop-allRelPrevSDlcpb_std)/allRelPrevSDlcpb_std, 
                                        smoothLatent_FBpop=100*(allRelPrevSDLCPb_FBpop-allRelPrevSDLCPb_std)/allRelPrevSDLCPb_std, 
                                      empiricial_2019=100*(allRelPrevSD_2019-allRelPrevSD_std)/allRelPrevSD_std, 
                                        latent_2019=100*(allRelPrevSDlcpb_2019-allRelPrevSDlcpb_std)/allRelPrevSDlcpb_std, 
                                        smoothLatent_2019=100*(allRelPrevSDLCPb_2019-allRelPrevSDLCPb_std)/allRelPrevSDLCPb_std, 
                                      empiricial_jittered=100*(allRelPrevSD_jittered-allRelPrevSD_std)/allRelPrevSD_std, 
                                        latent_jittered=100*(allRelPrevSDlcpb_jittered-allRelPrevSDlcpb_std)/allRelPrevSDlcpb_std, 
                                        smoothLatent_jittered=100*(allRelPrevSDLCPb_jittered-allRelPrevSDLCPb_std)/allRelPrevSDLCPb_std)
  
  # Calculate plotting range ----
  fullRangePredsPrevalence = range(as.matrix(fullPrevalenceTab[,-1]), na.rm=TRUE)
  fullRangePredsPctPrevalence = range(as.matrix(fullPrevalenceTabPct[,-1]), na.rm=TRUE)
  fullRangeSDsPrevalence = range(as.matrix(fullPrevalenceTabSD[,-1]), na.rm=TRUE)
  fullRangeSDsPctPrevalance = range(as.matrix(fullPrevalenceTabSDPct[,-1]), na.rm=TRUE)
  
  fullRangePredsCount = range(as.matrix(fullCountTab[,-1]), na.rm=TRUE)
  fullRangePredsPctCount = range(as.matrix(fullCountTabPct[,-1]), na.rm=TRUE)
  fullRangeSDsCount = range(as.matrix(fullCountTabSD[,-1]), na.rm=TRUE)
  fullRangeSDsPctCount = range(as.matrix(fullCountTabSDPct[,-1]), na.rm=TRUE)
  
  fullRangePredsRelPrev = range(as.matrix(fullRelPrevTab[,-1]), na.rm=TRUE)
  fullRangePredsPctRelPrev = range(as.matrix(fullRelPrevTabPct[,-1]), na.rm=TRUE)
  fullRangeSDsRelPrev = range(as.matrix(fullRelPrevTabSD[,-1]), na.rm=TRUE)
  fullRangeSDsPctRelPrev = range(as.matrix(fullRelPrevTabSDPct[,-1]), na.rm=TRUE)
  
  # Make tables ----
  
  # table of mean pct increase SD for each type of area and scenario
  # 3 sets columns: prev, burd, relPrev
  # 3 rows: constituency x stratum, constituency, county
  
  browser()
  
  # Empirical model ----
  # FBpop
  prevConStrat = mean(abs(fullPrevalenceTabPct$empiricial_FBpop[fullPrevalenceTabPct$areaLevel == "Constituency x stratum"]), na.rm=TRUE)
  prevCon = mean(abs(fullPrevalenceTabPct$empiricial_FBpop[fullPrevalenceTabPct$areaLevel == "Constituency"]))
  prevCounty = mean(abs(fullPrevalenceTabPct$empiricial_FBpop[fullPrevalenceTabPct$areaLevel == "County"]))
  prevConStratSD = mean(abs(fullPrevalenceTabSDPct$empiricial_FBpop[fullPrevalenceTabSDPct$areaLevel == "Constituency x stratum"]), na.rm=TRUE)
  prevConSD = mean(abs(fullPrevalenceTabSDPct$empiricial_FBpop[fullPrevalenceTabSDPct$areaLevel == "Constituency"]))
  prevCountySD = mean(abs(fullPrevalenceTabSDPct$empiricial_FBpop[fullPrevalenceTabSDPct$areaLevel == "County"]))
  prevMeans = cbind(c("conStrat"=prevConStrat, "constituency"=prevCon, "county"=prevCounty), 
                c("conStratSD"=prevConStratSD, "constituencySD"=prevConSD, "countySD"=prevCountySD))
  
  countConStrat = niceMean(abs(fullCountTabPct$empiricial_FBpop[fullCountTabPct$areaLevel == "Constituency x stratum"]))
  countCon = mean(abs(fullCountTabPct$empiricial_FBpop[fullCountTabPct$areaLevel == "Constituency"]))
  countCounty = mean(abs(fullCountTabPct$empiricial_FBpop[fullCountTabPct$areaLevel == "County"]))
  countConStratSD = niceMean(abs(fullCountTabSDPct$empiricial_FBpop[fullCountTabSDPct$areaLevel == "Constituency x stratum"]))
  countConSD = mean(abs(fullCountTabSDPct$empiricial_FBpop[fullCountTabSDPct$areaLevel == "Constituency"]))
  countCountySD = mean(abs(fullCountTabSDPct$empiricial_FBpop[fullCountTabSDPct$areaLevel == "County"]))
  countMeans = cbind(c("conStrat"=countConStrat, "constituency"=countCon, "county"=countCounty), 
                c("conStratSD"=countConStratSD, "constituencySD"=countConSD, "countySD"=countCountySD))
  
  relPrevConStrat = NA
  relPrevCon = niceMean(abs(fullRelPrevTabPct$empiricial_FBpop[fullRelPrevTabPct$areaLevel == "Constituency"]))
  relPrevCounty = niceMean(abs(fullRelPrevTabPct$empiricial_FBpop[fullRelPrevTabPct$areaLevel == "County"]))
  relPrevConStratSD = NA
  relPrevConSD = niceMean(abs(fullRelPrevTabSDPct$empiricial_FBpop[fullRelPrevTabSDPct$areaLevel == "Constituency"]))
  relPrevCountySD = niceMean(abs(fullRelPrevTabSDPct$empiricial_FBpop[fullRelPrevTabSDPct$areaLevel == "County"]))
  relPrevMeans = cbind(c("conStrat"=relPrevConStrat, "constituency"=relPrevCon, "county"=relPrevCounty), 
                c("conStratSD"=relPrevConStratSD, "constituencySD"=relPrevConSD, "countySD"=relPrevCountySD))
  
  # 2019
  prevConStrat = mean(abs(fullPrevalenceTabPct$empiricial_2019[fullPrevalenceTabPct$areaLevel == "Constituency x stratum"]), na.rm=TRUE)
  prevCon = mean(abs(fullPrevalenceTabPct$empiricial_2019[fullPrevalenceTabPct$areaLevel == "Constituency"]))
  prevCounty = mean(abs(fullPrevalenceTabPct$empiricial_2019[fullPrevalenceTabPct$areaLevel == "County"]))
  prevConStratSD = mean(abs(fullPrevalenceTabSDPct$empiricial_2019[fullPrevalenceTabSDPct$areaLevel == "Constituency x stratum"]), na.rm=TRUE)
  prevConSD = mean(abs(fullPrevalenceTabSDPct$empiricial_2019[fullPrevalenceTabSDPct$areaLevel == "Constituency"]))
  prevCountySD = mean(abs(fullPrevalenceTabSDPct$empiricial_2019[fullPrevalenceTabSDPct$areaLevel == "County"]))
  prevMeans = rbind(prevMeans, 
                    cbind(c("conStrat"=prevConStrat, "constituency"=prevCon, "county"=prevCounty), 
                      c("conStratSD"=prevConStratSD, "constituencySD"=prevConSD, "countySD"=prevCountySD)))
  
  countConStrat = niceMean(abs(fullCountTabPct$empiricial_2019[fullCountTabPct$areaLevel == "Constituency x stratum"]))
  countCon = mean(abs(fullCountTabPct$empiricial_2019[fullCountTabPct$areaLevel == "Constituency"]))
  countCounty = mean(abs(fullCountTabPct$empiricial_2019[fullCountTabPct$areaLevel == "County"]))
  countConStratSD = niceMean(abs(fullCountTabSDPct$empiricial_2019[fullCountTabSDPct$areaLevel == "Constituency x stratum"]))
  countConSD = mean(abs(fullCountTabSDPct$empiricial_2019[fullCountTabSDPct$areaLevel == "Constituency"]))
  countCountySD = mean(abs(fullCountTabSDPct$empiricial_2019[fullCountTabSDPct$areaLevel == "County"]))
  countMeans = rbind(countMeans, 
                     cbind(c("conStrat"=countConStrat, "constituency"=countCon, "county"=countCounty), 
                       c("conStratSD"=countConStratSD, "constituencySD"=countConSD, "countySD"=countCountySD)))
    
  relPrevConStrat = NA
  relPrevCon = niceMean(abs(fullRelPrevTabPct$empiricial_2019[fullRelPrevTabPct$areaLevel == "Constituency"]))
  relPrevCounty = niceMean(abs(fullRelPrevTabPct$empiricial_2019[fullRelPrevTabPct$areaLevel == "County"]))
  relPrevConStratSD = NA
  relPrevConSD = niceMean(abs(fullRelPrevTabSDPct$empiricial_2019[fullRelPrevTabSDPct$areaLevel == "Constituency"]))
  relPrevCountySD = niceMean(abs(fullRelPrevTabSDPct$empiricial_2019[fullRelPrevTabSDPct$areaLevel == "County"]))
  relPrevMeans = rbind(relPrevMeans, 
                       cbind(c("conStrat"=relPrevConStrat, "constituency"=relPrevCon, "county"=relPrevCounty), 
                         c("conStratSD"=relPrevConStratSD, "constituencySD"=relPrevConSD, "countySD"=relPrevCountySD)))
  
  # jittered
  prevConStrat = niceMean(abs(fullPrevalenceTabPct$empiricial_jittered[fullPrevalenceTabPct$areaLevel == "Constituency x stratum"]))
  prevCon = mean(abs(fullPrevalenceTabPct$empiricial_jittered[fullPrevalenceTabPct$areaLevel == "Constituency"]))
  prevCounty = mean(abs(fullPrevalenceTabPct$empiricial_jittered[fullPrevalenceTabPct$areaLevel == "County"]))
  prevConStratSD = niceMean(abs(fullPrevalenceTabSDPct$empiricial_jittered[fullPrevalenceTabSDPct$areaLevel == "Constituency x stratum"]))
  prevConSD = mean(abs(fullPrevalenceTabSDPct$empiricial_jittered[fullPrevalenceTabSDPct$areaLevel == "Constituency"]))
  prevCountySD = mean(abs(fullPrevalenceTabSDPct$empiricial_jittered[fullPrevalenceTabSDPct$areaLevel == "County"]))
  prevMeans = rbind(prevMeans, 
                cbind(c("conStrat"=prevConStrat, "constituency"=prevCon, "county"=prevCounty), 
                  c("conStratSD"=prevConStratSD, "constituencySD"=prevConSD, "countySD"=prevCountySD)))
  
  countConStrat = niceMean(abs(fullCountTabPct$empiricial_jittered[fullCountTabPct$areaLevel == "Constituency x stratum"]))
  countCon = mean(abs(fullCountTabPct$empiricial_jittered[fullCountTabPct$areaLevel == "Constituency"]))
  countCounty = mean(abs(fullCountTabPct$empiricial_jittered[fullCountTabPct$areaLevel == "County"]))
  countConStratSD = niceMean(abs(fullCountTabSDPct$empiricial_jittered[fullCountTabSDPct$areaLevel == "Constituency x stratum"]))
  countConSD = mean(abs(fullCountTabSDPct$empiricial_jittered[fullCountTabSDPct$areaLevel == "Constituency"]))
  countCountySD = mean(abs(fullCountTabSDPct$empiricial_jittered[fullCountTabSDPct$areaLevel == "County"]))
  countMeans = rbind(countMeans, 
                 cbind(c("conStrat"=countConStrat, "constituency"=countCon, "county"=countCounty), 
                   c("conStratSD"=countConStratSD, "constituencySD"=countConSD, "countySD"=countCountySD)))
  
  relPrevConStrat = NA
  relPrevCon = niceMean(abs(fullRelPrevTabPct$empiricial_jittered[fullRelPrevTabPct$areaLevel == "Constituency"]))
  relPrevCounty = niceMean(abs(fullRelPrevTabPct$empiricial_jittered[fullRelPrevTabPct$areaLevel == "County"]))
  relPrevConStratSD = NA
  relPrevConSD = niceMean(abs(fullRelPrevTabSDPct$empiricial_jittered[fullRelPrevTabSDPct$areaLevel == "Constituency"]))
  relPrevCountySD = niceMean(abs(fullRelPrevTabSDPct$empiricial_jittered[fullRelPrevTabSDPct$areaLevel == "County"]))
  relPrevMeans = rbind(relPrevMeans, 
                   cbind(c("conStrat"=relPrevConStrat, "constituency"=relPrevCon, "county"=relPrevCounty), 
                     c("conStratSD"=relPrevConStratSD, "constituencySD"=relPrevConSD, "countySD"=relPrevCountySD)))
  
  # combine into single table
  tab = cbind(prevMeans[,1], countMeans[,1], relPrevMeans[,1], 
              prevMeans[,2], countMeans[,2], relPrevMeans[,2])
  colnames(tab) = c("Prevalence", "Burden", "Relative prevalence", 
                    "PrevalenceS", "BurdenS", "Relative prevalenceS")
  row.names(tab) = rep(c("Constituency x stratum", "Constituency", "County"), 3)
  tab = round(tab, 0)
  
  print(xtable(tab, digits=0, 
               caption=paste0("Mean percent difference of empirical model posterior ", 
                              "mean and SD under sensitivity analysis scenarios ", 
                              "compared to under the original data."), 
               label=paste0("tab:pctIncreaseSensitivity_empirical")))
  
  
  # Make plots ----
  
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

logScaleAverage = function(x, y, maxPctDiff=20) {
  
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

# for each area (Admin1 and Admin2), check if estimated variance of pemp(A) 
# matches empirical estimate
testVarRatios = function(logisticApproximation=FALSE, coarse=TRUE) {
  # first load the model predictions
  if(logisticApproximation) {
    stop("logisticApproximation not currently used")
  }
  logisticText = ifelse(!logisticApproximation, "", "logisticApprox")
  coarseText = ifelse(!coarse, "", "Coarse")
  out = load(paste0("savedOutput/application/finalMort", coarseText, logisticText, ".RData"))
  
  # calculate r_sl^(2) since we'll need it
  varNames = row.names(riskOut$parameterSummaryTable)
  sigmaEps = riskOut$parameterSummaryTable[varNames == "errorSD",1]
  smoothRiskSqDraws = matrix(logitNormSqMean(cbind(c(as.matrix(logitDraws)), rep(sigmaEps, length(logitDraws)))), nrow=nrow(logitDraws))
  
  # get the population integration grid and check that it is consistent with pop totals
  if(coarse) {
    popMat = popGridCoarse
    targetPopMat = popGridCoarseAdjusted
  } else {
    popMat = popGrid
    targetPopMat = popGridAdjusted
  }
  # make sure the population integration grids match the population frame
  easpa = makeDefaultEASPA()
  poppsub = poppsubKenya
  out = checkPopFrameAndIntWeights(popMat=popMat, targetPopMat=targetPopMat, 
                                   easpa=easpa, poppsub=poppsub, 
                                   stopOnFrameMismatch=FALSE, 
                                   tol=1e-3)
  popMat = out$popMat
  targetPopMat = out$targetPopMat
  
  # Start with Admin1 areas ----
  Murb = easpa$EAUrb
  Mrur = easpa$EARur
  Nurb = easpa$popUrb
  Nrur = easpa$popRur
  
  estimatedVarRSLAdmin1 = numeric(nrow(poppaKenya))
  estimatedVarPempAdmin1 = numeric(nrow(poppaKenya))
  sampleVarRSLAdmin1 = numeric(nrow(poppaKenya))
  sampleVarPempAdmin1 = numeric(nrow(poppaKenya))
  
  estimatedVarBSLAdmin1 = numeric(nrow(poppaKenya))
  estimatedVarBempAdmin1 = numeric(nrow(poppaKenya))
  sampleVarBSLAdmin1 = numeric(nrow(poppaKenya))
  sampleVarBempAdmin1 = numeric(nrow(poppaKenya))
  for(i in 1:nrow(poppaKenya)) {
    thisArea = easpa$area[i]
    thisAreaInds = targetPopMat$area == thisArea
    
    q = targetPopMat[thisAreaInds,]$pop
    urbVec = targetPopMat[thisAreaInds,]$urban
    
    out = varPrevEmpStrat(Murb=Murb[i], Mrur=Mrur[i], Nurb=Nurb[i], Nrur=Nrur[i], 
                          q=q, urbVec=urbVec, smoothRiskSqDraws=smoothRiskSqDraws[thisAreaInds,], 
                          smoothRiskDraws=pixelPop$pSmoothRisk[thisAreaInds,], 
                          returnVarRSL=TRUE)
    estimatedVarRSLAdmin1[i] = out[2]
    estimatedVarPempAdmin1[i] = out[1]
    
    out = varBurdEmpStrat(Murb=Murb[i], Mrur=Mrur[i], Nurb=Nurb[i], Nrur=Nrur[i], 
                          q=q, urbVec=urbVec, smoothRiskSqDraws=smoothRiskSqDraws[thisAreaInds,], 
                          smoothRiskDraws=pixelPop$pSmoothRisk[thisAreaInds,], 
                          returnVarBSL=TRUE)
    estimatedVarBSLAdmin1[i] = out[2]
    estimatedVarBempAdmin1[i] = out[1]
    
    # TODO include relative prevalence
    
    # calculate sample variances
    sampleVarRSLAdmin1[i] = var(areaPop$pSmoothRisk[i,])
    sampleVarPempAdmin1[i] = var(areaPop$pFineScalePrevalence[i,])
    
    sampleVarBSLAdmin1[i] = var(areaPop$ZSmoothRisk[i,])
    sampleVarBempAdmin1[i] = var(areaPop$ZFineScalePrevalence[i,])
    
    # TODO include relative prevalence
  }
  
  # Do the same with Admin2 areas ----
  meanEAspsub = meanEAsPerCon()
  Murb = meanEAspsub$meanUrbanEAs
  Mrur = meanEAspsub$meanRuralEAs
  Nurb = meanEAspsub$popUrb
  Nrur = meanEAspsub$popRur
  
  estimatedVarRSLAdmin2 = numeric(nrow(poppsubKenya))
  estimatedVarPempAdmin2 = numeric(nrow(poppsubKenya))
  sampleVarRSLAdmin2 = numeric(nrow(poppsubKenya))
  sampleVarPempAdmin2 = numeric(nrow(poppsubKenya))
  
  estimatedVarBSLAdmin2 = numeric(nrow(poppsubKenya))
  estimatedVarBempAdmin2 = numeric(nrow(poppsubKenya))
  sampleVarBSLAdmin2 = numeric(nrow(poppsubKenya))
  sampleVarBempAdmin2 = numeric(nrow(poppsubKenya))
  for(i in 1:nrow(poppsubKenya)) {
    thisSubrea = poppsubKenya$subarea[i]
    thisSubareaInds = targetPopMat$subarea == thisSubrea
    
    q = targetPopMat[thisSubareaInds,]$pop
    urbVec = targetPopMat[thisSubareaInds,]$urban
    
    out = varPrevEmpStrat(Murb=Murb[i], Mrur=Mrur[i], Nurb=Nurb[i], Nrur=Nrur[i], 
                          q=q, urbVec=urbVec, smoothRiskSqDraws=smoothRiskSqDraws[thisSubareaInds,], 
                          smoothRiskDraws=pixelPop$pSmoothRisk[thisSubareaInds,], 
                          returnVarRSL=TRUE)
    estimatedVarRSLAdmin2[i] = out[2]
    estimatedVarPempAdmin2[i] = out[1]
    
    out = varBurdEmpStrat(Murb=Murb[i], Mrur=Mrur[i], Nurb=Nurb[i], Nrur=Nrur[i], 
                          q=q, urbVec=urbVec, smoothRiskSqDraws=smoothRiskSqDraws[thisSubareaInds,], 
                          smoothRiskDraws=pixelPop$pSmoothRisk[thisSubareaInds,], 
                          returnVarBSL=TRUE)
    estimatedVarBSLAdmin2[i] = out[2]
    estimatedVarBempAdmin2[i] = out[1]
    
    # TODO include relative prevalence
    
    # calculate sample variances
    sampleVarRSLAdmin2[i] = var(subareaPop$pSmoothRisk[i,])
    sampleVarPempAdmin2[i] = var(subareaPop$pFineScalePrevalence[i,])
    
    sampleVarBSLAdmin2[i] = var(subareaPop$ZSmoothRisk[i,])
    sampleVarBempAdmin2[i] = var(subareaPop$ZFineScalePrevalence[i,])
    
    # TODO include relative prevalence
  }
  
  # print comparisons ----
  print("Admin1 comparisons: ")
  print(paste0("mean((sampleVarRSLAdmin1 - estimatedVarRSLAdmin1)^2): ", mean((sampleVarRSLAdmin1 - estimatedVarRSLAdmin1)^2)))
  print(paste0("mean((sampleVarPempAdmin1 - estimatedVarPempAdmin1)^2): ", mean((sampleVarPempAdmin1 - estimatedVarPempAdmin1)^2)))
  print(paste0("mean((sampleVarPRatioAdmin1 - estimatedVarPRatioEmpAdmin1)^2): ", mean((sampleVarPempAdmin1/sampleVarRSLAdmin1 - 
                                                                                          estimatedVarPempAdmin1/estimatedVarRSLAdmin1)^2)))
  print(paste0("mean((sampleVarBSLAdmin1 - estimatedVarBSLAdmin1)^2): ", mean((sampleVarBSLAdmin1 - estimatedVarBSLAdmin1)^2)))
  print(paste0("mean((sampleVarBempAdmin1 - estimatedVarBempAdmin1)^2): ", mean((sampleVarBempAdmin1 - estimatedVarBempAdmin1)^2)))
  print(paste0("mean((sampleVarBRatioAdmin1 - estimatedVarBRatioEmpAdmin1)^2): ", mean((sampleVarBempAdmin1/sampleVarBSLAdmin1 - 
                                                                                          estimatedVarBempAdmin1/estimatedVarBSLAdmin1)^2)))
  
  print("Admin2 comparisons: ")
  print(paste0("mean((sampleVarRSLAdmin2 - estimatedVarRSLAdmin2)^2): ", mean((sampleVarRSLAdmin2 - estimatedVarRSLAdmin2)^2)))
  print(paste0("mean((sampleVarPempAdmin2 - estimatedVarPempAdmin2)^2): ", mean((sampleVarPempAdmin2 - estimatedVarPempAdmin2)^2)))
  print(paste0("mean((sampleVarPRatioAdmin2 - estimatedVarPRatioEmpAdmin2)^2): ", mean((sampleVarPempAdmin2/sampleVarRSLAdmin2 - 
                                                                                          estimatedVarPempAdmin2/estimatedVarRSLAdmin2)^2)))
  print(paste0("mean((sampleVarBSLAdmin2 - estimatedVarBSLAdmin2)^2): ", mean((sampleVarBSLAdmin2 - estimatedVarBSLAdmin2)^2)))
  print(paste0("mean((sampleVarBempAdmin2 - estimatedVarBempAdmin2)^2): ", mean((sampleVarBempAdmin2 - estimatedVarBempAdmin2)^2)))
  print(paste0("mean((sampleVarBRatioAdmin2 - estimatedVarBRatioEmpAdmin2)^2): ", mean((sampleVarBempAdmin2/sampleVarBSLAdmin2 - 
                                                                                          estimatedVarBempAdmin2/estimatedVarBSLAdmin2)^2)))
  
  # plot comparisons
  browser()
}




