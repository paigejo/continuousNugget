# IWSM plots

# Make plots for the neonatal mortality application
makeIWSMPlots = function(logisticApproximation=FALSE, signif=.95) {
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
  countyMap = adm1
  constituencyMap = adm2
  
  # make color scales
  meanCols=makeRedBlueDivergingColors(64, rev = TRUE)
  sdCols=makeBlueGreenYellowSequentialColors(64)
  popCols=makePurpleYellowSequentialColors(64, rev=TRUE)
  urbCols=makeGreenBlueSequentialColors(64)
  
  getWidth = function(x) {
    diff(quantile(x, prob=c(alpha/2, 1-alpha/2), na.rm=TRUE))
  }
  
  ## 1 x 2 plot of predictions (prevalence mean and CI width)
  
  # plot mean
  pixelMean = rowMeans(pixelPop$pSmoothRisk)
  constituencyMean = rowMeans(subareaPop$pFineScaleRisk)
  countyMean = rowMeans(areaPop$pFineScaleRisk)
  # provinceMean = rowMeans(agg$aggregatedResultsLCPb$regionMatrices$p)
  meanRange = range(pixelMean, constituencyMean, countyMean)
  # widthRange = range(c(rangePrevalenceCIWidthPixel, 
  #                     rangePrevalenceCIWidthConstituency, 
  #                     rangePrevalenceCIWidthCounty))
  
  browser()
  pdf(paste0(figDirectory, "IWSM/prevalenceMeanCIWidth", logisticText, ".pdf"), width=10, height=5)
  par(mfrow=c(1,2), mar=c(4.2, 4.0, 1.2, 3.1), oma=c(0, 0, 0, 0.7))
  
  meanTicks = pretty(meanRange, n=5)
  meanTickLabels = as.character(meanTicks)
  
  # constituency level
  plotMapDat(plotVar=constituencyMean, mapDat=constituencyMap, new = TRUE, 
             main="Mean Prevalence", scaleFun=logit, scaleFunInverse=expit, 
             cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
             xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.88,.91,.2,.9)), legend.width=3, 
             plotArgs=list(cex.main=1, cex.axis=1, cex.lab=1), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="Longitude", ylab="Latitude")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  # constituency level
  widthTicksSubarea = pretty(rangePrevalenceCIWidthConstituency, n=8)
  widthTickLabelsSubarea = as.character(widthTicksSubarea)
  widthRangeSubarea = rangePrevalenceCIWidthConstituency
  plotMapDat(plotVar=prevalenceCIWidthConstituency, mapDat=constituencyMap, new = TRUE, 
             main=paste0("Prevalence ", 100*signif, "% CI Width"), scaleFun=logit, scaleFunInverse=expit, 
             cols=sdCols, ticks=widthTicksSubarea, tickLabels=widthTickLabelsSubarea, 
             zlim=logit(widthRangeSubarea), xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.88,.91,.2,.9)), legend.width=3, 
             plotArgs=list(cex.main=1, cex.axis=1, cex.lab=1), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="Longitude", ylab="")
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  points(mort$lon, mort$lat, pch=19, cex=.1)
  
  dev.off()
}

plotGridResolutionFig = function() {
  out = load("savedOutput/simpleExample/gridResolutionTestNairobi_1_100.RData")
  
  ##### Plot results ----
  
  # truePrevalenceConstituencyKenya
  # truePrevalenceCountyKenya
  
  # for each truth, calculate coverage and other metrics.
  # In aggregate, generate:
  #   boxplot of individual truth coverages versus res.
  #   scatterplot of average (over all truths) coverage vs res.
  #   pair plot (over all truths and resolutions) of central predictions
  #   boxplot (over all truths) of 80% and 95% CI widths
  
  maxSamples = max(nSamples)
  maxSamples = min(nSamples)
  maxSamples = 1000
  allPredsPrevalence = list()
  allPredsGriddedRisk = list()
  allCIWidthsPrevalence = list()
  allCIWidthsGriddedRisk = list()
  allCoveragesPrevalence = list()
  allCoveragesGriddedRisk = list()
  for(i in 1:length(truths)) {
    # get truth
    truePrevalenceConstituencyKenya = truths[[i]]$truePrevalenceConstituencyKenya
    
    # Calculate RMSE, 95% Coverage
    predsPrevalence = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
    predsGriddedRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
    residsConstituencyPrevalence = list()
    residsConstituencyGriddedRisk = list()
    for(j in 1:length(resolutions)) {
      thisMaxSamples = min(c(nSamples[j], maxSamples))
      predsPrevalence[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScalePrevalence[,1:thisMaxSamples])
      predsGriddedRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pGriddedRisk[,1:thisMaxSamples])
      theseResidsPrevalence = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScalePrevalence[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
      theseResidsGriddedRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pGriddedRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
      residsConstituencyPrevalence = c(residsConstituencyPrevalence, list(theseResidsPrevalence))
      residsConstituencyGriddedRisk = c(residsConstituencyGriddedRisk, list(theseResidsGriddedRisk))
    }
    allPredsPrevalence = c(allPredsPrevalence, list(predsPrevalence))
    allPredsGriddedRisk = c(allPredsGriddedRisk, list(predsGriddedRisk))
    
    lowConstituencyPrevalence = lapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
    lowConstituencyGriddedRisk = lapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
    highConstituencyPrevalence = lapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
    highConstituencyGriddedRisk = lapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
    
    # CIWidthSmoothRisk = highConstituencySmoothRisk - lowConstituencySmoothRisk
    # CIWidthRisk = highConstituencyRisk - lowConstituencyRisk
    # CIWidthPrevalence = highConstituencyPrevalence - lowConstituencyPrevalence
    # CIWidthGriddedRisk = highConstituencyGriddedRisk - lowConstituencyGriddedRisk
    CIWidthPrevalence = lapply(residsConstituencyPrevalence, 
                               function(mat) {
                                 apply(mat, 1, function(x) {
                                   quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                               })
    CIWidthGriddedRisk = lapply(residsConstituencyGriddedRisk, 
                                function(mat) {
                                  apply(mat, 1, function(x) {
                                    quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                                })
    allCIWidthsPrevalence = c(allCIWidthsPrevalence, list(CIWidthPrevalence))
    allCIWidthsGriddedRisk = c(allCIWidthsGriddedRisk, list(CIWidthGriddedRisk))
    
    # this is never used anyway, so commented out:
    # meanCIWidthSmoothRisk = colMeans(CIWidthSmoothRisk)
    # meanCIWidthRisk = colMeans(CIWidthRisk)
    # meanCIWidthPrevalence = colMeans(CIWidthPrevalence)
    # meanCIWidthGriddedRisk = colMeans(CIWidthGriddedRisk)
    
    # inCISmoothRisk = (0 <= highConstituencySmoothRisk) & (0 >= lowConstituencySmoothRisk)
    # inCIRisk = (0 <= highConstituencyRisk) & (0 >= lowConstituencyRisk)
    # inCIPrevalence = (0 <= highConstituencyPrevalence) & (0 >= lowConstituencyPrevalence)
    # inCIGriddedRisk = (0 <= highConstituencyGriddedRisk) & (0 >= lowConstituencyGriddedRisk)
    
    inCIPrevalence = lapply(residsConstituencyPrevalence, 
                            function(mat) {
                              apply(mat, 1, function(x) {
                                (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                            })
    inCIGriddedRisk = lapply(residsConstituencyGriddedRisk, 
                             function(mat) {
                               apply(mat, 1, function(x) {
                                 (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                             })
    
    # coverageSmoothRisk = colMeans(inCISmoothRisk)
    # coverageRisk = colMeans(inCIRisk)
    # coveragePrevalence = colMeans(inCIPrevalence)
    # coverageGriddedRisk = colMeans(inCIGriddedRisk)
    coveragePrevalence = sapply(inCIPrevalence, rowMeans)
    coverageGriddedRisk = sapply(inCIGriddedRisk, rowMeans)
    allCoveragesPrevalence = c(allCoveragesPrevalence, list(coveragePrevalence))
    allCoveragesGriddedRisk = c(allCoveragesGriddedRisk, list(coverageGriddedRisk))
  }
  
  # compile relevant results
  
  ## central predictions
  allPredsPrevalenceMat = do.call("rbind", allPredsPrevalence)
  allPredsGriddedRiskMat = do.call("rbind", allPredsGriddedRisk)
  
  ## CI widths
  allCIWidthsPrevalence95 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[1,]})
  })
  allCIWidthsPrevalence95 = do.call("rbind", allCIWidthsPrevalence95)
  
  allCIWidthsGriddedRisk95 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[1,]})
  })
  allCIWidthsGriddedRisk95 = do.call("rbind", allCIWidthsGriddedRisk95)
  
  ## coverages
  allCoveragesPrevalence95 = sapply(allCoveragesPrevalence, function(x) {x[1,]})
  allCoveragesGriddedRisk95 = sapply(allCoveragesGriddedRisk, function(x) {x[1,]})
  
  # Plot CI Widths versus resolution and model
  
  ## 95% CIs
  CIWidth = c(c(allCIWidthsPrevalence95), c(allCIWidthsGriddedRisk95))
  tempRes = resolutions[col(allCIWidthsPrevalence95)]
  tempCon = factor(as.character(poppsubSimple$subarea[col(allCIWidthsPrevalence95)]))
  N=length(tempCon)
  CIWidthFrame = data.frame(Constituency=rep(tempCon, 2), 
                            Resolution=rep(tempRes, 2), 
                            Model=factor(c(rep("Prevalence", N), rep("Gridded risk", N)), 
                                         levels=c("Prevalence", "Gridded risk")), 
                            CIWidth=CIWidth)
  
  pdf(paste0("figures/gridResolutionTest/CIWidthVRes95.pdf"), width=7, height=5)
  ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
    geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
    labs(x="Grid resolution (km)", y="95% credible interval width", fill="Model") + 
    theme_classic()
  dev.off()
}





