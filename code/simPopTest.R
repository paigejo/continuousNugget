
pixelPopToAreaTest = function(pixelLevelPop, eaSamples, areas, stratifyByUrban=TRUE, targetPopMat=NULL, 
                          doFineScaleRisk=!is.null(pixelLevelPop$fineScaleRisk$p), 
                          doSmoothRisk=!is.null(pixelLevelPop$smoothRisk$p)) {
  
  # fine scale prevalence aggregation model
  nSamples = pixelLevelPop$NFineScalePrevalence
  zSamples = pixelLevelPop$ZFineScalePrevalence
  zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works
  out = aggPixelPredsTest(Zg=zSamples, Ng=nSamples, areas=areas, targetPopMat=targetPopMat, 
                      useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
  aggregationResults = out$aggregationResults
  aggregationMatrices = out$aggregationMatrices
  names(aggregationResults)[-1] = paste(names(aggregationResults)[-1], "FineScalePrevalence", sep="")
  names(aggregationMatrices) = paste(names(aggregationMatrices), "FineScalePrevalence", sep="")
  
  # fine scale risk aggregation model
  if(doFineScaleRisk) {
    nSamplesFineScaleRisk = pixelLevelPop$NFineScaleRisk
    zSamplesFineScaleRisk = pixelLevelPop$ZFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPredsTest(Zg=zSamplesFineScaleRisk, Ng=nSamplesFineScaleRisk, areas=areas, targetPopMat=targetPopMat, 
                        useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
    aggregationResultsFineScaleRisk = out$aggregationResults
    aggregationMatricesFineScaleRisk = out$aggregationMatrices
    names(aggregationResultsFineScaleRisk)[-1] = paste(names(aggregationResultsFineScaleRisk)[-1], "FineScaleRisk", sep="")
    names(aggregationMatricesFineScaleRisk) = paste(names(aggregationMatricesFineScaleRisk), "FineScaleRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsFineScaleRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsFineScaleRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesFineScaleRisk)
  }
  
  if(doSmoothRisk) {
    # NOTE: although useDensity is set to FALSE, that's only because the density is already 
    # directly incorporated into nSamplesSmoothRisk
    nSamplesSmoothRisk = pixelLevelPop$NSmoothRisk
    zSamplesSmoothRisk = pixelLevelPop$ZSmoothRisk
    zSamplesSmoothRisk[is.na(zSamplesSmoothRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPredsTest(Zg=zSamplesSmoothRisk, Ng=nSamplesSmoothRisk, areas=areas, targetPopMat=targetPopMat, 
                        useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
    aggregationResultsSmoothRisk = out$aggregationResults
    aggregationMatricesSmoothRisk = out$aggregationMatrices
    names(aggregationResultsSmoothRisk)[-1] = paste(names(aggregationResultsSmoothRisk)[-1], "SmoothRisk", sep="")
    names(aggregationMatricesSmoothRisk) = paste(names(aggregationMatricesSmoothRisk), "SmoothRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsSmoothRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsSmoothRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesSmoothRisk)
  }
  
  list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
}

#' Helper function of \code{\link{pixelPopToArea}}
#' 
#' Aggregates population from the 
#' pixel level to the level of the area of interest.
#' 
#' @param Zg nIntegrationPoint x nsim matrix of simulated response (population numerators) for each pixel and sample
#' @param Ng nIntegrationPoint x nsim matrix of simulated counts (population denominators) for each pixel and sample
#' @param areas nIntegrationPoint length character vector of areas (or subareas) 
#' @param urban nIntegrationPoint length vector of indicators specifying whether or not pixels are urban or rural
#' @param targetPopMat same as in \code{\link{simPopCustom}}
#' @param useDensity whether to use population density as aggregation weights. 
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param normalize if TRUE, pixel level aggregation weights within specified area are normalized to sum to 1. This produces an 
#' average of the values in Zg rather than a sum. In general, should only be set to TRUE for smooth integrals of risk.
aggPixelPredsTest = function(Zg, Ng, areas, urban=targetPopMat$urban, targetPopMat=NULL, useDensity=FALSE, 
                         stratifyByUrban=TRUE, normalize=useDensity) {
  
  if(useDensity && !normalize) {
    stop("if useDensity is set to TRUE, normalize must be set to TRUE as well")
  }
  predsUrban = urban
  predsArea = areas
  
  # set NAs and pixels without any sample size to 0
  Ng[is.na(Ng)] = 0
  if(!useDensity) {
    # is useDensity is true, then Zg is really a set of probabilities, so no need to set to 0
    Zg[Ng == 0] = 0
  }
  
  # function to aggregate predictions over the 
  # population density grid. Use the following function to get numerical 
  # integration matrix for a given level of areal aggregation. returned 
  # matrices have dimension length(unique(areaNames)) x length(areaNames)
  # areaNames: 
  # urbanProportions: DEPRACATED vector giving proportion of population urban for each unique area in areaNames. 
  #                   If specified, ensure that urban and rural parts of the full integration 
  #                   matrix have the appropriate relative weights for each area. Used for population 
  #                   density based integration
  # normalize: whether or not to normalize the rows of the matrices to sum to 1 or to instead 
  #            contain only binary values (or non-binary values based on the binary values if 
  #            urbanProportions is not NULL)
  getIntegrationMatrix = function(areaNames, urbanProportions=NULL, normalize=FALSE) {
    
    if(useDensity) {
      popDensities = targetPopMat$pop
      densities = popDensities
    } else {
      equalDensities = rep(1, nrow(Zg))
      densities = equalDensities
    }
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL, thisUseDensity=useDensity, thisNormalize=normalize) {
      areaI = areaNames == uniqueNames[i]
      
      if(thisUseDensity) {
        theseDensities = popDensities
      } else {
        theseDensities = equalDensities
      }
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban != thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0 && thisNormalize)
        theseDensities * (1/thisSum)
      else if(thisSum == 0)
        rep(0, length(theseDensities))
      else
        theseDensities
    }
    
    if(!stratifyByUrban) {
      integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      
      integrationMatrix
    } else {
      integrationMatrixUrban = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE), ncol=length(uniqueNames)))
      integrationMatrixRural = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE), ncol=length(uniqueNames)))
      if(!is.null(urbanProportions)) {
        integrationMatrix = sweep(integrationMatrixUrban, 1, urbanProportions, "*") + sweep(integrationMatrixRural, 1, 1-urbanProportions, "*")
      } else {
        integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      }
      
      rownames(integrationMatrix) = uniqueNames
      rownames(integrationMatrixUrban) = uniqueNames
      rownames(integrationMatrixRural) = uniqueNames
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames, normalize=normalize)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!stratifyByUrban) {
      ZAggregated = A %*% Zg
      NAggregated = A %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      aggregationResults = list(p=pAggregated, Z=ZAggregated, N=NAggregated)
      aggregationMatrices = list(A=A, AUrban=NULL, ARural=NULL)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      # first aggregate the numerator. The denominator will depend on the aggregation method
      ZAggregated = A %*% Zg
      ZAggregatedUrban = AUrban %*% Zg
      ZAggregatedRural = ARural %*% Zg
      
      if(useDensity) {
        # for population density aggregation, we integrate probabilities rather than aggregate 
        # empirical proportions
        NAggregated = NULL
        pAggregated = ZAggregated
        ZAggregated = NULL
        
        NAggregatedUrban = NULL
        pAggregatedUrban = ZAggregatedUrban
        ZAggregatedUrban = NULL
        
        NAggregatedRural = NULL
        pAggregatedRural = ZAggregatedRural
        ZAggregatedRural = NULL
      } else {
        # if we do not use density, we must also aggregate the denominator to calculate 
        # the aggregated empirical proportions
        NAggregated = A %*% Ng
        pAggregated = ZAggregated / NAggregated
        pAggregated[NAggregated == 0] = NA
        
        NAggregatedUrban = AUrban %*% Ng
        pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
        pAggregatedUrban[NAggregatedUrban == 0] = NA
        
        NAggregatedRural = ARural %*% Ng
        pAggregatedRural = ZAggregatedRural / NAggregatedRural
        pAggregatedRural[NAggregatedRural == 0] = NA
      }
      
      aggregationResults = list(region=sort(unique(areas)), p=pAggregated, Z=ZAggregated, N=NAggregated, 
                                pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                                pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
      aggregationMatrices = list(A=A, AUrban=AUrban, ARural=ARural)
    }
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  areaResults = getIntegratedPredictions(predsArea)
  
  # return results
  areaResults
}

#' @describeIn aggPop Aggregate areal populations to another areal level
#' @export
areaPopToAreaTest = function(areaLevelPop, areasFrom, areasTo, 
                         stratifyByUrban=TRUE, 
                         doFineScaleRisk=!is.null(areaLevelPop$aggregationResults$pFineScaleRisk), 
                         doSmoothRisk=!is.null(areaLevelPop$aggregationResults$pSmoothRisk)) {
  
  if(length(areasFrom) != length(unique(areasFrom))) {
    stop("areasFrom must contain only unique names of areas to which we want to aggregate")
  }
  
  uniqueNames = sort(unique(areasTo))
  
  # construct row of the aggregation matrix given toArea index
  getMatrixHelper = function(i) {
    # get areasTo associated with this fromArea
    thisToArea = uniqueNames[i]
    thisFromAreas = unique(areasFrom[areasTo == thisToArea])
    areaI = areasFrom %in% thisFromAreas
    
    areaI
  }
  
  # construct the aggregation matrix from areasFrom to areasTo
  A = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
  rownames(A) = uniqueNames
  colnames(A) = areasFrom
  
  ##### aggregate populations
  # fine scale prevalence aggregation model
  getaggregationResults = function(resultNameRoot="FineScalePrevalence") {
    if(! (paste("N", resultNameRoot, sep="") %in% names(areaLevelPop$aggregationResults))) {
      stop(paste0(resultNameRoot, " was not computed in input areaLevelPop"))
    }
    
    nSamples = areaLevelPop$aggregationResults[[paste("N", resultNameRoot, sep="")]]
    zSamples = areaLevelPop$aggregationResults[[paste("Z", resultNameRoot, sep="")]]
    zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    ZAggregated =  A %*% zSamples
    NAggregated =  A %*% nSamples
    pAggregated = ZAggregated / NAggregated
    pAggregated[NAggregated == 0] = NA
    thisA=A %*% areaLevelPop$aggregationMatrices[[paste("A", resultNameRoot, sep="")]]
    rownames(thisA) = uniqueNames
    
    if(stratifyByUrban) {
      nSamplesUrban = areaLevelPop$aggregationResults[[paste("NUrban", resultNameRoot, sep="")]]
      zSamplesUrban = areaLevelPop$aggregationResults[[paste("ZUrban", resultNameRoot, sep="")]]
      zSamplesUrban[is.na(zSamplesUrban)] = 0 # must set to zero temporarily so matrix multiplication works out
      
      nSamplesRural = areaLevelPop$aggregationResults[[paste("NRural", resultNameRoot, sep="")]]
      zSamplesRural = areaLevelPop$aggregationResults[[paste("ZRural", resultNameRoot, sep="")]]
      zSamplesRural[is.na(zSamplesRural)] = 0 # must set to zero temporarily so matrix multiplication works out
      
      ZAggregatedUrban =  A %*% zSamplesUrban
      NAggregatedUrban =  A %*% nSamplesUrban
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = NA
      thisAUrban=A %*% areaLevelPop$aggregationMatrices[[paste("AUrban", resultNameRoot, sep="")]]
      rownames(thisAUrban) = uniqueNames
      
      ZAggregatedRural =  A %*% zSamplesRural
      NAggregatedRural =  A %*% nSamplesRural
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = NA
      thisARural=A %*% areaLevelPop$aggregationMatrices[[paste("ARural", resultNameRoot, sep="")]]
      rownames(thisARural) = uniqueNames
    } else {
      ZAggregatedUrban = NULL
      NAggregatedUrban = NULL
      pAggregatedUrban = NULL
      thisAUrban=NULL
      
      ZAggregatedRural = NULL
      NAggregatedRural = NULL
      pAggregatedRural = NULL
      thisARural=NULL
    }
    
    aggregationResults = list(region=uniqueNames, p=pAggregated, Z=ZAggregated, N=NAggregated, 
                              pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                              pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
    aggregationMatrices = list(A=thisA, AUrban=thisAUrban, ARural=thisARural)
    
    capitalResultNameRoot = resultNameRoot
    # capitalResultNameRoot = paste(toupper(substr(capitalResultNameRoot, 1, 1)), substr(capitalResultNameRoot, 2, nchar(capitalResultNameRoot)), sep="")
    names(aggregationResults)[-1] = paste(names(aggregationResults)[-1], capitalResultNameRoot, sep="")
    names(aggregationMatrices) = paste(names(aggregationMatrices)[-1], capitalResultNameRoot, sep="")
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  # fine scale prevalence model 
  out = getaggregationResults("FineScalePrevalence")
  aggregationResults = out$aggregationResults
  aggregationMatrices = out$aggregationMatrices
  
  # fine scale risk aggregation model
  if(doFineScaleRisk) {
    out = getaggregationResults("FineScaleRisk")
    resFineScaleRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resFineScaleRisk, by="region")
    aggregationResults = c(aggregationResults, resFineScaleRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  if(doSmoothRisk) {
    out = getaggregationResults("SmoothRisk")
    resSmoothRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resSmoothRisk, by="region")
    aggregationResults = c(aggregationResults, resSmoothRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
}

#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses custom spatial logit risk function and can include iid cluster 
#' level effect.
#' 
#' @export
simPopCustomTest = function(logitRiskDraws, sigmaEpsilonDraws, easpa, popMat, targetPopMat, 
                        stratifyByUrban=TRUE, validationPixelI=NULL, validationClusterI=NULL, 
                        clustersPerPixel=NULL, 
                        doFineScaleRisk=FALSE, doSmoothRisk=FALSE, doGriddedRisk=FALSE, 
                        doSmoothRiskLogisticApprox=FALSE, 
                        poppsub=NULL, subareaLevel=FALSE, gridLevel=FALSE, 
                        min1PerSubarea=TRUE, 
                        fixPopPerEA=NULL, fixHHPerEA=NULL, fixPopPerHH=NULL, 
                        returnEAinfo=FALSE, epsc=NULL, stopOnFrameMismatch=FALSE, tol=1e-3, 
                        verbose=TRUE, doGC=FALSE) {
  
  time1 = proc.time()[3]
  
  if(!is.null(validationPixelI) || !is.null(validationClusterI) || !is.null(clustersPerPixel)) {
    stop("validationPixelI, validationClusterI, and clustersPerPixel not yet fully implemented")
  }
  
  if(any(c(!is.null(fixPopPerEA), !is.null(fixHHPerEA), !is.null(fixPopPerHH)))) {
    if(verbose) {
      print("adjusting easpa population and households based on fixPopPerEA, fixHHPerEA, and fixPopPerHH")
    }
  }
  if(!is.null(fixPopPerEA)) {
    easpa$popUrb = easpa$EAUrb * fixPopPerEA
    easpa$popRur = easpa$EARur * fixPopPerEA
    easpa$popTotal = easpa$EATotal * fixPopPerEA
  }
  if(!is.null(fixHHPerEA)) {
    easpa$HHUrb = easpa$EAUrb * fixHHPerEA
    easpa$HHRur = easpa$EARur * fixHHPerEA
    easpa$HHTotal = easpa$EATotal * fixHHPerEA
  }
  if(!is.null(fixPopPerHH)) {
    if(fixPopPerHH * fixHHPerEA != fixPopPerEA) {
      stop("fixPopPerHH * fixHHPerEA != fixPopPerEA")
    }
  }
  
  # make sure the population integration grids match the population frame
  out = checkPopFrameAndIntWeights(popMat=popMat, targetPopMat=targetPopMat, 
                                   easpa=easpa, poppsub=poppsub, 
                                   stopOnFrameMismatch=stopOnFrameMismatch, 
                                   tol=tol)
  popMat = out$popMat
  targetPopMat = out$targetPopMat
  
  nDraws = ncol(logitRiskDraws)
  
  # set default inputs
  totalEAs = sum(easpa$EATotal)
  if(!is.null(clustersPerPixel)) {
    emptyPixels = clustersPerPixel == 0
    if(totalEAs != sum(clustersPerPixel))
      stop("sum(easpa$EATotal) != sum(clustersPerPixel)")
  }
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # determine if we care about subareas (smallest areas we care about. No info of EAs per subarea)
  sampleBySubarea = !is.null(popMat$subarea) && !is.null(poppsub)
  
  ##### Line 1 (of the algorithm): take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  time2 = proc.time()[3]
  
  print("drawing EAs")
  pixelProbs = popMat$pop
  
  # take draws from the stratified binomial process for each posterior sample
  if(is.null(clustersPerPixel)) {
    if(verbose) {
      print("Sampling EAs per pixel...")
    }
    if(sampleBySubarea) {
      eaSamples = rStratifiedMultnomialBySubareaTest(nDraws, popMat, easpa, stratifyByUrban, poppsub=poppsub, 
                                                          min1PerSubarea=min1PerSubarea)
    } else {
      eaSamples = rStratifiedMultnomialTest(nDraws, popMat, easpa, stratifyByUrban)
    }
  }
  
  if(!is.null(clustersPerPixel) && !exists("eaSamples")) {
    eaSamples = matrix(rep(clustersPerPixel, nDraws), ncol=nDraws)
  }
  time3 = proc.time()[3]
  
  # make matrix (or list) of pixel indices mapping matrices of EA values to matrices of pixel values
  if(!is.null(clustersPerPixel)) {
    pixelIndices = rep(1:nrow(popMat), times=clustersPerPixel) # this contains repetitions and has length == nEAs
    uniquePixelIndices = sort(unique(pixelIndices))
  } else {
    pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  }
  
  # determine which EAs are urban if necessary
  if(stratifyByUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    if(!is.null(clustersPerPixel)) {
      urbanVals = popMat$urban[pixelIndices]
      uniqueUrbanVals = popMat$urban[uniquePixelIndices]
    }else {
      urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
    }
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  if(!is.null(clustersPerPixel)) {
    areaVals = popMat$area[pixelIndices]
    uniqueAreaVals = popMat$area[uniquePixelIndices]
  } else {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  }
  time4 = proc.time()[3]
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  if(verbose) {
    print("simulating EA level risks, numerators, and denominators")
  }
  if(is.null(epsc)) {
    epsc = matrix(stats::rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
  }
  time5 = proc.time()[3]
  
  ##### Line 3: draw EA population denominators, N
  
  if(!is.null(clustersPerPixel)) {
    if(is.null(validationPixelI))
      stop("clustersPerPixel must only be set for validation, but validationPixelI is NULL")
    
    # in this case, every left out cluster has exactly 25 households. Simply sample target population 
    # with equal probability from each cluster/faux EA
    if(is.null(fixPopPerEA)) {
      Ncs = sampleNMultilevelMultinomialFixedTest(clustersPerPixel, nDraws=nDraws, pixelIndices=pixelIndices, 
                                                       urbanVals=urbanVals, areaVals=areaVals, easpa=easpa, popMat=popMat, stratifyByUrban=stratifyByUrban, 
                                                       verbose=verbose)
    } else {
      if(verbose) {
        print(paste0("drawing Ns for each EA..."))
      }
      Ncs = matrix(rep(fixPopPerEA, totalEAs*nDraws), ncol=nDraws)
    }
  } else {
    out = sampleNMultilevelMultinomialTest(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                                popMat=popMat, stratifyByUrban=stratifyByUrban, verbose=verbose, returnEAinfo=returnEAinfo, 
                                                fixPopPerHH=fixPopPerHH)
    if(returnEAinfo) {
      householdDraws = out$householdDraws
      Ncs = out$targetPopDraws
    } else {
      Ncs = out
      householdDraws = NULL
    }
  }
  time6 = proc.time()[3]
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  if(verbose) {
    print("Calculating EA level risk and number of deaths...")
  }
  if(!is.null(clustersPerPixel)) {
    uc = logitRiskDraws[pixelIndices,]
    muc = expit(uc + epsc)
  } else {
    uc = matrix(logitRiskDraws[cbind(rep(rep(1:nrow(logitRiskDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
    muc = expit(uc + epsc)
  }
  time7 = proc.time()[3]
  
  # calculate Z_{ic} for each EA in each pixel
  Zcs = matrix(stats::rbinom(n=totalEAs * nDraws, size=Ncs, prob=as.matrix(muc)), ncol=nDraws)
  time8 = proc.time()[3]
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, vals, applyFun=sum, popWeightMatrix=NULL) {
    # calculate levels over which to aggregate
    if(!is.null(clustersPerPixel)) {
      indices = pixelIndices
    } else {
      indices = factor(as.character(pixelIndexMat[,i]))
    }
    
    # in this case (the LCPb model), we calculate weighted means within factor levels using popWeightMatrix
    if(!is.null(popWeightMatrix)) {
      stop("using popWeightMatrix is no longer support, since this is much slower than calculating normalized 
           weights separately and multiplying values by them outside this function")
      # applyFun = function(x) {stats::weighted.mean(x, popWeightMatrix[,i], na.rm=TRUE)}
      
      Data = data.frame(v=vals[,i], w=popWeightMatrix[,i])
      out = sapply(split(Data, indices), function(x) stats::weighted.mean(x$v,x$w))
    } else {
      if(!is.null(clustersPerPixel)) {
        # out = tapply(vals[,i], factor(as.character(pixelIndices)), FUN=applyFun)
        out = tapply(vals[,i], indices, FUN=applyFun)
      } else {
        out = tapply(vals[,i], indices, FUN=applyFun)
      }
    }
    
    if(!is.null(clustersPerPixel)) {
      returnValues = out
    } else {
      indices = as.numeric(names(out))
      
      returnValues = rep(NA, nrow(logitRiskDraws))
      returnValues[indices] = out
    }
    
    returnValues
  }
  
  ##### Line 5: We already did this, resulting in logitRiskDraws input
  
  ##### Line 6: aggregate population denominators for each grid cell to get N_{ig}
  if(verbose) {
    print("Aggregating from EA level to the pixel level...")
  }
  
  # Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, vals=Ncs)
  # Ng[is.na(Ng)] = 0
  # 
  # ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  # Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=Zc)
  # 
  # ##### Line 8: Calculate empirical mortality proportions for each grid cell, p_{ig}. 
  # #####         Whenever N_{ig} is 0, set p_{ig} to NA as well
  # pg = Zg / Ng
  # pg[Ng == 0] = NA
  
  # even if we aren't producing grid level results in general, gridded and smooth risk 
  # require it
  griddedRisk = NULL
  zSamplesGriddedRisk = NULL
  nSamplesGriddedRisk = NULL
  smoothRisk = NULL
  zSamplesSmoothRisk = NULL
  nSamplesSmoothRisk = NULL
  if(doGriddedRisk) {
    # NOTE: these cluster effects aren't necessarily in agreement with the 
    #       cluster effects of the fine scale models, since they're drawn 
    #       independently
    pixelEps = matrix(rnorm(length(logitRiskDraws), mean=0, sd=rep(sigmaEpsilonDraws, each=nrow(logitRiskDraws))), ncol=nDraws)
    griddedRisk = expit(logitRiskDraws + pixelEps)
    
    # get expected numerator and denominators
    nSamplesGriddedRisk = matrix(rep(targetPopMat$pop, nDraws), ncol=nDraws)
    zSamplesGriddedRisk = sweep(griddedRisk, 1, targetPopMat$pop, "*")
  }
  time9 = proc.time()[3]
  
  if(doSmoothRisk) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(logitNormMean(cbind(c(as.matrix(logitRiskDraws)), rep(sigmaEpsilonDraws, each=nrow(logitRiskDraws))), logisticApprox=doSmoothRiskLogisticApprox), nrow=nrow(logitRiskDraws))
    
    # get expected numerator and denominators
    nSamplesSmoothRisk = matrix(rep(targetPopMat$pop, nDraws), ncol=nDraws)
    zSamplesSmoothRisk = sweep(smoothRisk, 1, targetPopMat$pop, "*")
  }
  time10 = proc.time()[3]
  
  ZcsFineScaleRisk = NULL
  NcsSamplesFineScaleRisk = NULL
  if(gridLevel) {
    # calculate pixel/grid level results
    # time1=system.time(out <- aggPredsVariablePerArea(popNumerators=Zcs, popDenominators=Ncs, 
    #                               areaMat=pixelIndexMat, areaLevels=1:nrow(popMat)))
    # time2=system.time(out <- aggPredsVariablePerAreaDT(popNumerators=Zcs, popDenominators=Ncs, 
    #                               areaMat=pixelIndexMat, areaLevels=1:nrow(popMat)))
    out <- aggPredsVariablePerArea(popNumerators=Zcs, popDenominators=Ncs, 
                                   areaMat=pixelIndexMat, areaLevels=1:nrow(popMat))
    Ng = out$N
    Zg = out$Z
    pg = out$p
    
    time11 = proc.time()[3]
    
    ##### calculate results for the other models if necessary
    
    if(doFineScaleRisk) {
      ## the following code is for a previous model that is no longer used, and so is commented out:
      # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
      # fineScaleRisk = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
      # fineScaleRisk[!is.finite(fineScaleRisk)] = NA
      # nPerEA = SUMMER:::getExpectedNperEA(easpa, targetPopMat)
      # nSamplesFineScaleRisk = sweep(eaSamples, 1, nPerEA, "*")
      # zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0
      
      # We aggregate just like for the fine scale prevalence model, but using the expected population 
      # numerators conditional on the true denominators.
      NcsFineScaleRisk = Ncs
      ZcsFineScaleRisk = muc * NcsFineScaleRisk
      out = aggPredsVariablePerArea(popNumerators=ZcsFineScaleRisk, popDenominators=Ncs, 
                                    areaMat=pixelIndexMat, areaLevels=1:nrow(popMat))
      nSamplesFineScaleRisk = out$N
      zSamplesFineScaleRisk = out$Z
      fineScaleRisk = out$p
    } else {
      nSamplesFineScaleRisk = NULL
      zSamplesFineScaleRisk = NULL
    }
    time12 = proc.time()[3]
    
    ##### Extra steps: collect draws at each level and generate:
    ##### areas, preds, 
    pixelLevelPop = list(pFineScalePrevalence=pg, ZFineScalePrevalence=Zg, NFineScalePrevalence=Ng, 
                         pFineScaleRisk=fineScaleRisk, ZFineScaleRisk=zSamplesFineScaleRisk, NFineScaleRisk=nSamplesFineScaleRisk, 
                         pSmoothRisk=smoothRisk, ZSmoothRisk=zSamplesSmoothRisk, NSmoothRisk=nSamplesSmoothRisk, 
                         pGriddedRisk=griddedRisk, ZGriddedRisk=zSamplesGriddedRisk, NGriddedRisk=nSamplesGriddedRisk)
  } else if(doGriddedRisk || doSmoothRisk) {
    pixelLevelPop = list(pFineScalePrevalence=NULL, ZFineScalePrevalence=NULL, NFineScalePrevalence=NULL, 
                         pFineScaleRisk=NULL, ZFineScaleRisk=NULL, NFineScaleRisk=NULL, 
                         pSmoothRisk=smoothRisk, ZSmoothRisk=zSamplesSmoothRisk, NSmoothRisk=nSamplesSmoothRisk, 
                         pGriddedRisk=griddedRisk, ZGriddedRisk=zSamplesGriddedRisk, NGriddedRisk=nSamplesGriddedRisk)
    
    time11 = proc.time()[3]
    time12 = proc.time()[3]
  } else {
    pixelLevelPop = NULL
    
    time11 = proc.time()[3]
    time12 = proc.time()[3]
  }
  
  if(subareaLevel) {
    # aggregate up from the most aggregated set of results we have: pixel/grid level if we have them, 
    # EA level results otherwise.
    if(verbose) {
      print("Aggregating to the subarea level...")
    }
    
    subareas = sort(unique(popMat$subarea))
    if(gridLevel) {
      # first get results for fine scale prevalence
      subareaLevelPop = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$subarea, 
                                       stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                                       doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                                       doGriddedRisk=doGriddedRisk)
      time13 = subareaLevelPop$rawTimes[2]
      time14 = subareaLevelPop$rawTimes[3]
      time15 = subareaLevelPop$rawTimes[4]
      time16 = subareaLevelPop$rawTimes[5]
      time17 = subareaLevelPop$rawTimes[6]
    } else {
      # we aggregate EA level results to subarea level. This takes a bit more leg work
      subareaLevelPop = list()
      
      ## create matrices of relevant variable values at the EA level
      subareaMat = matrix(popMat$subarea[pixelIndexMat], ncol=ncol(pixelIndexMat))
      if(stratifyByUrban) {
        urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=ncol(pixelIndexMat))
      } else {
        urbanMat = NULL
      }
      
      time13 = proc.time()[3]
      
      # fine scale prevalence
      fineScalePrevalenceSubarea <- aggPredsVariablePerArea(popNumerators=Zcs, popDenominators=Ncs, 
                                                            areaMat=subareaMat, urbanMat=urbanMat, 
                                                            areaLevels=sort(as.character(poppsub$subarea)))
      # time2 = system.time(fineScalePrevalenceSubarea <- aggPredsVariablePerAreaDT(popNumerators=Zcs, popDenominators=Ncs,
      # areaMat=subareaMat, urbanMat=urbanMat,
      # areaLevels=sort(as.character(poppsub$subarea))))
      names(fineScalePrevalenceSubarea)[-1] = paste(names(fineScalePrevalenceSubarea)[-1], "FineScalePrevalence", sep="")
      subareaLevelPop = c(subareaLevelPop, fineScalePrevalenceSubarea)
      time14 = proc.time()[3]
      
      # fine scale risk
      if(doFineScaleRisk) {
        NcsFineScaleRisk = Ncs
        ZcsFineScaleRisk = NcsFineScaleRisk * muc
        fineScaleRiskSubarea = aggPredsVariablePerArea(popNumerators=ZcsFineScaleRisk, 
                                                       popDenominators=NcsFineScaleRisk, 
                                                       areaMat=subareaMat, urbanMat=urbanMat, 
                                                       areaLevels=sort(as.character(poppsub$subarea)))
        names(fineScaleRiskSubarea)[-1] = paste(names(fineScaleRiskSubarea)[-1], "FineScaleRisk", sep="")
        subareaLevelPop = c(subareaLevelPop, fineScaleRiskSubarea)
      }
      time15 = proc.time()[3]
      
      # smooth risk
      if(doSmoothRisk || doGriddedRisk) {
        # we have pixel level results, so use them instead of EA level results since it's faster
        
        smoothGriddedRiskSubarea = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$subarea, 
                                                  stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                                                  doFineScalePrevalence=FALSE, doFineScaleRisk=FALSE, 
                                                  doSmoothRisk=doSmoothRisk, doGriddedRisk=doGriddedRisk)
        smoothGriddedRiskSubarea$aggregationResults$region = NULL
        subareaLevelPop = c(subareaLevelPop, smoothGriddedRiskSubarea$aggregationResults)
        
        time16 = smoothGriddedRiskSubarea$rawTimes[5]
        time17 = smoothGriddedRiskSubarea$rawTimes[6]
      }
      
      # aggregation matrices are not applicable for aggregating from EA level to subarea level
      subareaLevelPop = list(aggregationResults=subareaLevelPop, aggregationMatrices=NULL)
    }
  } else {
    subareaLevelPop = NULL
  }
  
  # aggregate to the area level. Aggregate up from largest possible level to save time
  if(verbose) {
    print("Aggregating to the area level...")
  }
  if(subareaLevel) {
    areasFrom = subareaLevelPop$aggregationResults$region
    areasTo = poppsub$area[match(subareaLevelPop$aggregationResults$region, 
                                 poppsub$subarea)]
    areaLevelPop = areaPopToArea(subareaLevelPop, areasFrom=areasFrom, areasTo=areasTo, 
                                 stratifyByUrban=stratifyByUrban, doFineScaleRisk=doFineScaleRisk, 
                                 doSmoothRisk=doSmoothRisk, doGriddedRisk=doGriddedRisk)
    
    time18 = areaLevelPop$rawTimes[2]
    time19 = areaLevelPop$rawTimes[3]
    time20 = areaLevelPop$rawTimes[4]
    time21 = areaLevelPop$rawTimes[5]
    time22 = areaLevelPop$rawTimes[6]
  } else if(gridLevel) {
    areaLevelPop = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$areas, 
                                  stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                                  doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                                  doGriddedRisk=doGriddedRisk)
    
    time18 = areaLevelPop$rawTimes[2]
    time19 = areaLevelPop$rawTimes[3]
    time20 = areaLevelPop$rawTimes[4]
    time21 = areaLevelPop$rawTimes[5]
    time22 = areaLevelPop$rawTimes[6]
  } else {
    # In this case we must aggregate up from EA level for all risk/prevalence types except 
    # those only defined at the pixel level, and als aggregate those from the pixel level
    
    # First aggregate from EA level
    EALevelPop = list(pixelIndexMat=pixelIndexMat, 
                      ZFineScalePrevalence=Zcs, NFineScalePrevalence=Ncs, 
                      ZFineScaleRisk=ZcsFineScaleRisk, NFineScaleRisk=Ncs)
    
    areaLevelPop = EAPopToArea(EALevelPop, areaMat=areaMat, areaLevels=sort(as.character(poppsub$subarea)), 
                               urbanMat=urbanMat, doFineScalePrevalence=TRUE, 
                               doFineScaleRisk=doFineScaleRisk)
    time18 = areaLevelPop$rawTimes[2]
    time19 = areaLevelPop$rawTimes[3]
    time20 = areaLevelPop$rawTimes[4]
    
    # now aggregate grid level models
    if(doSmoothRisk || doGriddedRisk) {
      out = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$area, 
                           stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                           doFineScalePrevalence=FALSE, doFineScaleRisk=FALSE, 
                           doSmoothRisk=doSmoothRisk, doGriddedRisk=doGriddedRisk)
      out = out$aggregationResults
      out$region = NULL
      
      time21 = out$rawTimes[5]
      time22 = out$rawTimes[6]
    }
    
    areaLevelPop = c(areaLevelPop, out)
  }
  
  # setup computation time table
  env = environment()
  rawTimes = sapply(paste("time", 1:22, sep=""), get, envir=env)
  allTimings = diff(rawTimes)
  names(allTimings) = c("setup", "samplingEAs", "EAtoPixelIndexMapAndVariables", 
                        "drawingClusterEffects", "drawingNc", "drawing_uc_muc", 
                        "drawingZc", "getPixelGriddedRisk", "getPixelSmoothRisk", 
                        "getPixelPrevalence", "getPixelRisk", "subareaAggSetup", 
                        "getSubareaPrevalence", "getSubareaRisk", 
                        "getSubareaSmoothRisk", "getSubareaGriddedRisk", 
                        "areaAggSetup", "getAreaPrevalence", "getAreaRisk", 
                        "getAreaSmoothRisk", "getAreaGriddedRisk")
  temp = as.list(allTimings)
  processedTimings = c(prevalence = temp$setup + temp$samplingEAs + 
                         temp$EAtoPixelIndexMapAndVariables + 
                         temp$drawingClusterEffects + temp$drawingNc + 
                         temp$drawing_uc_muc + temp$drawingZc + 
                         temp$getPixelPrevalence + temp$subareaAggSetup +
                         temp$getSubareaPrevalence + temp$areaAggSetup + 
                         temp$getAreaPrevalence, 
                       risk = temp$setup + temp$samplingEAs + 
                         temp$EAtoPixelIndexMapAndVariables + 
                         temp$drawingClusterEffects + temp$drawingNc + 
                         temp$drawing_uc_muc + 
                         temp$getPixelRisk + temp$subareaAggSetup +
                         temp$getSubareaRisk + temp$areaAggSetup + 
                         temp$getAreaRisk, 
                       smoothRisk = temp$setup + 
                         temp$drawing_uc_muc + 
                         temp$getPixelSmoothRisk + temp$subareaAggSetup +
                         temp$getSubareaSmoothRisk + temp$areaAggSetup + 
                         temp$getAreaSmoothRisk, 
                       griddedRisk = temp$setup + 
                         temp$drawing_uc_muc + 
                         temp$getPixelGriddedRisk + temp$subareaAggSetup +
                         temp$getSubareaGriddedRisk + temp$areaAggSetup + 
                         temp$getAreaGriddedRisk, 
                       total = sum(allTimings))
  
  if(!returnEAinfo) {
    totalTime = sum(allTimings)
    allTimings = c(allTimings, totalTime=totalTime)
    allTimings = cbind(allTimings, allTimings/totalTime)
    list(pixelPop=pixelLevelPop, subareaPop=subareaLevelPop, areaPop=areaLevelPop, 
         allTimings=allTimings, processedTimings=processedTimings)
  } else {
    
    # return list of eaDat objects
    getEADat = function(i) {
      theseI = pixelIndexMat[,i]
      
      eaDat = data.frame(lon=popMat$lon[theseI], lat=popMat$lat[theseI], 
                         area=popMat$area[theseI], subarea=rep("temp", length(theseI)), 
                         urban=popMat$urban[theseI], east=popMat$east[theseI], north=popMat$north[theseI], 
                         popDensity=popMat$pop[theseI], popDensityTarget=targetPopMat$pop[theseI], pixelIs=theseI, 
                         nHH=householdDraws[,i], N=Ncs[,i], Z=Zcs[,i], 
                         pFineScalePrevalence=Zcs[,i]/Ncs[,i])
      if(doFineScaleRisk) {
        eaDat$pFineScaleRisk = muc[,i]
        eaDat$ZFineScaleRisk = muc[,i] * Ncs[,i]
        eaDat$NFineScaleRisk = Ncs[,i]
      }
      if(doSmoothRisk) {
        eaDat$pSmoothRisk=smoothRisk[theseI,i]
        # N and Z not defined at the EA level
      }
      if(doGriddedRisk) {
        eaDat$pGriddedRisk=griddedRisk[theseI,i]
        # N and Z not defined at the EA level
      }
      
      if("subarea" %in% names(popMat)) {
        eaDat$subarea = popMat$subarea[theseI]
      } else {
        eaDat$subarea = NULL
      }
      
      eaDat$pFineScalePrevalence[eaDat$N == 0] = NA
      
      eaDat
    }
    
    if(verbose) {
      print("Constructing list of simulated EA data.frames...")
    }
    
    eaDatList = lapply(1:nDraws, getEADat)
    time23 = proc.time()[3]
    
    env = environment()
    rawTimes = sapply(paste("time", 1:23, sep=""), get, envir=env)
    allTimings = diff(rawTimes)
    names(allTimings) = c("setup", "samplingEAs", "EAtoPixelIndexMapAndVariables", 
                          "drawingClusterEffects", "drawingNc", "drawing_uc_muc", 
                          "drawingZc", "getPixelGriddedRisk", "getPixelSmoothRisk", 
                          "getPixelPrevalence", "getPixelRisk", "subareaAggSetup", 
                          "getSubareaPrevalence", "getSubareaRisk", 
                          "getSubareaSmoothRisk", "getSubareaGriddedRisk", 
                          "areaAggSetup", "getAreaPrevalence", "getAreaRisk", 
                          "getAreaSmoothRisk", "getAreaGriddedRisk", "getEAdat")
    temp = as.list(allTimings)
    processedTimings = c(prevalence = temp$setup + temp$samplingEAs + 
                           temp$EAtoPixelIndexMapAndVariables + 
                           temp$drawingClusterEffects + temp$drawingNc + 
                           temp$drawing_uc_muc + temp$drawingZc + 
                           temp$getPixelPrevalence + temp$subareaAggSetup +
                           temp$getSubareaPrevalence + temp$areaAggSetup + 
                           temp$getAreaPrevalence, 
                         risk = temp$setup + temp$samplingEAs + 
                           temp$EAtoPixelIndexMapAndVariables + 
                           temp$drawingClusterEffects + temp$drawingNc + 
                           temp$drawing_uc_muc + 
                           temp$getPixelRisk + temp$subareaAggSetup +
                           temp$getSubareaRisk + temp$areaAggSetup + 
                           temp$getAreaRisk, 
                         smoothRisk = temp$setup + 
                           temp$drawing_uc_muc + 
                           temp$getPixelSmoothRisk + temp$subareaAggSetup +
                           temp$getSubareaSmoothRisk + temp$areaAggSetup + 
                           temp$getAreaSmoothRisk, 
                         griddedRisk = temp$setup + 
                           temp$drawing_uc_muc + 
                           temp$getPixelGriddedRisk + temp$subareaAggSetup +
                           temp$getSubareaGriddedRisk + temp$areaAggSetup + 
                           temp$getAreaGriddedRisk, 
                         total = sum(allTimings))
    
    totalTime = sum(allTimings)
    allTimings = c(allTimings, totalTime=totalTime)
    allTimings = cbind(allTimings, allTimings/totalTime)
    
    list(pixelPop=pixelLevelPop, subareaPop=subareaLevelPop, areaPop=areaLevelPop, 
         eaDatList=eaDatList, eaSamples=eaSamples, 
         allTimings=allTimings, processedTimings=processedTimings)
  }
}


#' Internal functions for population simulation
#' 
#' Functions for calculating valuable quantities and for drawing from important 
#' distributions for population simulation.
#' 
#' @param easpa Census frame. See \code{\link{simPopCustomTest}} for details
#' @param popMat data.frame of pixellated grid of population densities. See \code{\link{simPopCustomTest}} for details
#' @param i Index
#' @param urban If TRUE, calculate only for urban part of the area. If FALSE, for only rural part
#' @param stratifyByUrban whether or not to stratify calculations by urban/rural classification
#' @param validationPixelI CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for which we want to simulate populations (used for pixel level validation)
#' @param n Number of samples
#' @param poppsub Population per subarea. See \code{\link{simPopCustomTest}} for details
#' @param min1PerSubarea Whether or not to ensure there is at least 1 EA per subarea. See \code{\link{simPopCustom}} for details
#' @param method If min1PerSubarea is TRUE, the sampling method for the truncated multinomial to use with rmulitnom1. rmultinom1Test automatically 
#'         switches between them depending on the number of expected samples. The methods are:
#' \describe{
#'   \item{mult1}{rejection sampling from multinomial plus 1 in each category}
#'   \item{mult}{rejection sampling from multinomial if any category has zero count}
#'   \item{indepMH}{independent Metropolis-Hastings using multinomial plus 1 distribution as proposal}
#' }
#' @param minSample The minimum number of samples per `chunk` of samples for truncated multinomial sampling. Defaults to 1
#' @param easpsub This could either be total EAs per subarea, or subarea crossed with urban or 
#'          rural if stratifyByUrban is TRUE
#' @param size Multinomial size parameter. See \code{\link[stats]{rmultinom}}
#' @param prob Multinomial probability vector parameter. See \code{\link[stats]{rmultinom}}
#' @param maxSize The maximum number of elements in a matrix drawn from the proposal distribution per sample chunk. 
#' @param maxExpectedSizeBeforeSwitch Max expected number of samples / k, the number of categories, before switching method
#' @param init Initial sample if method is `indepMH`
#' @param burnIn Number of initial samples before samples are collected if method is `indepMH`
#' @param filterEvery Store only every filterEvery samples if method is i`indepMH`
#' @param zeroProbZeroSamples If TRUE, set samples for parts of prob vector that are zero to zero. Otherwise they are set to one.
#' @param allowSizeLessThanK If TRUE, then if size < the number of categories (k), returns matrix where each 
#'                     column is vector of size ones and k - size zeros. If FALSE, throws an error if size < k
#' @param clustersPerPixel CURRENTLY FOR TESTING PURPOSES ONLY a vector of length nIntegrationPoints specifying the number of clusters per pixel if they are fixed
#' @param pixelIndices A nEA x n matrix of pixel indices associated with each EA per simulation/draw
#' @param urbanVals A nEA x n matrix of urbanicities associated with each EA per simulation/draw
#' @param areaVals A nEA x n matrix of area names associated with each EA per simulation/draw
#' @param easpaList A list of length n with each element being of the format of easpa 
#'            giving the number of households and EAs 
#'            per stratum. It is assumed that the number of EAs per stratum is 
#'            the same in each list element. If easpaList is a data frame, 
#'            number of households per stratum is assumed constant
#' @param nDraws Number of draws
#' @param pixelIndexMat Matrix of pixel indices associated with each EA and draw. Not 
#' required by getExpectedNperEA unless level == "EA"
#' @param urbanMat Matrix of urbanicities associated with each EA and draw
#' @param areaMat Matrix of areas associated with each EA and draw
#' @param verbose Whether to print progress as the function proceeds
#' @param returnEAinfo Whether a data frame at the EA level is desired
#' @param minHHPerEA The minimum number of households per EA (defaults to 25, since 
#' that is the number of households sampled per DHS cluster)
#' @param fixHHPerEA If not NULL, the fixed number of households per EA
#' @param fixPopPerHH If not NULL, the fixed target population per household
#' @param level Whether to calculate results at the integration grid or EA level
#' @name simPopInternal
NULL

#' @describeIn simPopInternal Calculates expected denominator per enumeration area.
getExpectedNperEATest = function(easpa, popMat, level=c("grid", "EA"), pixelIndexMat=NULL) {
  level = match.arg(level)
  
  # calculate the expected denominator per enumeration area in each stratum. 
  nPerEAUrban = easpa$popUrb / easpa$EAUrb
  nPerEARural = easpa$popRur / easpa$EARur
  
  # expanded the expected denominator values vector to be of length equal 
  # to the number of grid cells
  uniqueAreas = sort(unique(popMat$area))
  outUrban = numeric(nrow(popMat))
  outRural = numeric(nrow(popMat))
  for(i in 1:length(uniqueAreas)) {
    urbanI = getSortIndicesTest(i, urban=TRUE, popMat=popMat, stratifyByUrban=TRUE)
    ruralI = getSortIndicesTest(i, urban=FALSE, popMat=popMat, stratifyByUrban=TRUE)
    outUrban[urbanI] = nPerEAUrban[i]
    outRural[ruralI] = nPerEARural[i]
  }
  
  out = outUrban + outRural
  
  if(level == "EA") {
    if(is.null(pixelIndexMat)) {
      stop("if level == EA, must specify pixelIndexMat")
    }
    out = matrix(outPixel[pixelIndexMat], ncol=ncol(pixelIndexMat))
  }
  
  out
}

#' @describeIn simPopInternal For recombining separate multinomials into the draws over all grid points
getSortIndicesTest = function(i, urban=TRUE, popMat, stratifyByUrban=TRUE, validationPixelI=NULL) {
  
  # get area names
  areas = sort(unique(popMat$area))
  
  # determine which pixels and how many EAs are in this stratum
  if(stratifyByUrban) {
    includeI = popMat$area == areas[i] & popMat$urban == urban
  }
  else {
    includeI = popMat$area == areas[i]
  }
  
  # include only indices included within validation if necessary
  if(!is.null(validationPixelI)) {
    
    # convert validationPixelI into a logical
    temp = rep(FALSE, length(includeI))
    temp[validationPixelI] = TRUE
    
    # include only indices we are interested in for the validation
    includeI = includeI & temp
  }
  
  which(includeI)
}

#' @describeIn simPopInternal Gives nIntegrationPoints x n matrix of draws from the stratified multinomial with values 
#' corresponding to the value of |C^g| for each pixel, g (the number of EAs/pixel)
rStratifiedMultnomialTest = function(n, popMat, easpa, stratifyByUrban=TRUE) {
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum. Start by 
  # creating matrix of all draws of |C^g|
  eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
  
  # now draw multinomials
  if(stratifyByUrban) {
    # draw for each area crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomialTest, n=n, urban=TRUE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpa=easpa))
    ruralSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomialTest, n=n, urban=FALSE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    urbanIndices = unlist(sapply(1:length(areas), getSortIndicesTest, urban=TRUE, popMat=popMat, stratifyByUrban=stratifyByUrban))
    ruralIndices = unlist(sapply(1:length(areas), getSortIndicesTest, urban=FALSE, popMat=popMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples
    eaSamples[urbanIndices,] = urbanSamples
    eaSamples[ruralIndices,] = ruralSamples
  } else {
    # draw for each area
    stratumSamples = rbind(sapply(1:length(areas), n=n, rMyMultinomialTest, 
                                  stratifyByUrban=stratifyByUrban, popMat=popMat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    stratumIndices = c(sapply(1:length(areas), getSortIndicesTest, popMat=popMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples
    eaSamples[stratumIndices,] = stratumSamples
  }
  
  # return results
  eaSamples
}

#' @describeIn simPopInternal Gives nIntegrationPoints x n matrix of draws from the stratified multinomial with values 
# corresponding to the number of EAs in each pixel
rStratifiedMultnomialBySubareaTest = function(n, popMat, easpa, stratifyByUrban=TRUE, poppsub=NULL, 
                                          min1PerSubarea=TRUE, minSample=1) {
  if(is.null(poppsub)) {
    # use popMat to calculate poppsub
    stop("Calculating poppsub with popMat not yet implemented")
    
    # poppsub: a table with the following variables:
    # subarea
    # area
    # popUrb
    # popRur
    # popTotal
  }
  
  # get area names
  areas = sort(unique(popMat$area))
  subareas = sort(unique(popMat$subarea))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum
  
  # create temporary popMat, except with one row for each constituency
  popSubareaMat = popMat[1:length(subareas),]
  popSubareaMat$area = poppsub$area
  popSubareaMat$subarea = poppsub$subarea
  if(stratifyByUrban) {
    popSubareaMat$urban = FALSE
    popSubareaMat = rbind(popSubareaMat, popSubareaMat)
    popSubareaMat$urban[1:length(subareas)] = TRUE
    popSubareaMat$pop[1:length(subareas)] = poppsub$popUrb
    popSubareaMat$pop[(length(subareas) + 1):(2 * length(subareas))] = poppsub$popRur
  } else {
    popSubareaMat$pop = poppsub$popTotal
  }
  # browser()
  # now draw multinomials
  if(stratifyByUrban) {
    # draw for each constituency in each area crossed with urban/rural
    urbanSamplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomialTest, n=n, urban=TRUE, 
                                              stratifyByUrban=stratifyByUrban, popMat=popSubareaMat, easpa=easpa, 
                                              min1PerSubarea=min1PerSubarea, method="mult", 
                                              minSample=minSample))
    ruralSamplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomialTest, n=n, urban=FALSE, 
                                              stratifyByUrban=stratifyByUrban, popMat=popSubareaMat, easpa=easpa, 
                                              min1PerSubarea=min1PerSubarea, method="mult", 
                                              minSample=minSample))
    
    # get the indices used to recombine into the full set of draws for the subareas
    urbanIndicesCon = unlist(sapply(1:length(areas), getSortIndicesTest, urban=TRUE, popMat=popSubareaMat, stratifyByUrban=stratifyByUrban))
    ruralIndicesCon = unlist(sapply(1:length(areas), getSortIndicesTest, urban=FALSE, popMat=popSubareaMat, stratifyByUrban=stratifyByUrban)) - length(urbanIndicesCon)
    
    # recombine into eaSamples for the subareas
    urbanSamplesCon[urbanIndicesCon,] = urbanSamplesCon
    ruralSamplesCon[ruralIndicesCon,] = ruralSamplesCon
    
    # draw for each pixel crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(subareas), rMyMultinomialSubareaTest, n=n, urban=TRUE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpsub=urbanSamplesCon))
    ruralSamples = do.call("rbind", lapply(1:length(subareas), rMyMultinomialSubareaTest, n=n, urban=FALSE, 
                                           stratifyByUrban=stratifyByUrban, popMat=popMat, easpsub=ruralSamplesCon))
    
    # get the indices used to recombine into the full set of draws
    tempPopMat = popMat
    tempPopMat$area = tempPopMat$subarea
    urbanIndices = unlist(sapply(1:length(subareas), getSortIndicesTest, urban=TRUE, popMat=tempPopMat, stratifyByUrban=stratifyByUrban))
    ruralIndices = unlist(sapply(1:length(subareas), getSortIndicesTest, urban=FALSE, popMat=tempPopMat, stratifyByUrban=stratifyByUrban))
    
    # create matrix of all draws of |C^g| and recombine urban/rural draws into eaSamples
    eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
    eaSamples[urbanIndices,] = urbanSamples
    eaSamples[ruralIndices,] = ruralSamples
  } else {
    # draw for each constituency in each area crossed with urban/rural
    samplesCon = do.call("rbind", lapply(1:length(areas), rMyMultinomialTest, n=n, urban=TRUE, 
                                         stratifyByUrban=stratifyByUrban, popMat=popSubareaMat, easpa=easpa, 
                                         min1PerSubarea=min1PerSubarea, method="mult", 
                                         minSample=minSample))
    
    # get the indices used to recombine into the full set of draws for the subareas
    indicesCon = unlist(sapply(1:length(areas), getSortIndicesTest, popMat=popSubareaMat, stratifyByUrban=stratifyByUrban))
    
    # recombine into eaSamples for the subareas
    samplesCon[indicesCon,] = samplesCon
    
    # draw for each pixel in each constituency
    stratumSamples = rbind(sapply(1:length(subareas), n=n, rMyMultinomialSubareaTest, 
                                  stratifyByUrban=stratifyByUrban, popMat=popMat, easpsub=samplesCon))
    
    # get the indices used to recombine into the full set of draws
    tempPopMat = popMat
    tempPopMat$area = tempPopMat$subarea
    stratumIndices = c(sapply(1:length(subareas), getSortIndicesTest, popMat=tempPopMat, stratifyByUrban=stratifyByUrban))
    
    # create matrix of all draws of |C^g|
    eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
    eaSamples[stratumIndices,] = stratumSamples
  }
  
  # return results
  eaSamples
}

#' @describeIn simPopInternal 
rMyMultinomialTest = function(n, i, stratifyByUrban=TRUE, urban=TRUE, popMat=NULL, easpa=NULL, min1PerSubarea=FALSE, 
                          method=c("mult1", "mult", "indepMH"), minSample=1) {
  method = match.arg(method)
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # determine which pixels and how many EAs are in this stratum
  if(stratifyByUrban) {
    includeI = popMat$area == areas[i] & popMat$urban == urban
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EARur[i])
  }
  else {
    includeI = popMat$area == areas[i]
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EATotal[i])
  }
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0)
    return(matrix(nrow=0, ncol=n))
  thesePixelProbs = popMat$pop[includeI]
  if(any(thesePixelProbs > 0)) {
    if(!min1PerSubarea) {
      stats::rmultinom(n, nEA, prob=thesePixelProbs)
    } else {
      rmultinom1Test(n, nEA, prob=thesePixelProbs, method=method, allowSizeLessThanK=TRUE, minSample=minSample)
    }
  } else {
    matrix(0, nrow=length(thesePixelProbs), ncol=n)
  }
}

#' @describeIn simPopInternal 
rMyMultinomialSubareaTest = function(n, i, easpsub, stratifyByUrban=TRUE, urban=TRUE, popMat=NULL) {
  
  # get constituency names
  subareas = sort(unique(popMat$subarea))
  
  # determine which pixels and how many EAs are in this stratum
  if(stratifyByUrban) {
    includeI = popMat$subarea == subareas[i] & popMat$urban == urban
  }
  else {
    includeI = popMat$subarea == subareas[i]
  }
  nEA = easpsub[i,]
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0){
    if(any(nEA != 0))
      stop(paste0("no valid pixels to put EAs in for constituency ", as.character(subareas[i]), " and urban level ", urban))
    return(matrix(nrow=0, ncol=n))
  }
  thesePixelProbs = popMat$pop[includeI]
  sapply(nEA, stats::rmultinom, n=1, prob=thesePixelProbs)
}

#' @describeIn simPopInternal Random (truncated) multinomial draws conditional on the number of each type being at least one
rmultinom1Test = function(n=1, size, prob, maxSize=5000*5000, method=c("mult1", "mult", "indepMH"), verbose=FALSE, minSample=100, 
                      maxExpectedSizeBeforeSwitch=1000*1e7, init=NULL, burnIn=floor(n/4), filterEvery=10, zeroProbZeroSamples=TRUE, 
                      allowSizeLessThanK=FALSE) {
  method = match.arg(method)
  prob = prob*(1/sum(prob))
  
  if(zeroProbZeroSamples && any(prob == 0)) {
    zero = prob == 0
    out = matrix(0, nrow=length(prob), ncol=n)
    
    if(sum(!zero) > 0) {
      out[!zero,] = rmultinom1Test(n, size, prob[!zero], maxSize, method, verbose, minSample, maxExpectedSizeBeforeSwitch, init, burnIn, 
                               filterEvery, zeroProbZeroSamples, allowSizeLessThanK)
    }
    
    return(out)
  }
  
  k = length(prob)
  if(allowSizeLessThanK && (size <= k)) {
    return(replicate(n, as.numeric(1:k %in% sample(1:k, size, replace=FALSE))))
  } else if(size < k) {
    stop("size < k but rmultinom1Test requires at least 1 sample per multinomial type")
  }
  
  maxSamples = floor(maxSize / k)
  averageProbMult = prod((size/k)*prob)
  
  if(method != "indepMH")
    samples = matrix(NA, nrow=k, ncol=n)
  else
    samples = matrix(NA, nrow=k, ncol=round(n*filterEvery))
  if(method == "mult1") {
    averagex = 1 + (size-k)*prob
    averageProb = (size-k) / (prod(averagex))
    
    while(any(is.na(samples))) {
      # calculate the number of remaining samples
      samplesLeft = sum(apply(samples, 2, function(x) {any(is.na(x))}))
      
      # approximate expected number of samples so that, after some are rejected, we will 
      # have the right number of samples
      expectedSamples = ceiling(samplesLeft/averageProb)
      
      if(expectedSamples*k > maxExpectedSizeBeforeSwitch) {
        warning("too many samples expected with method=='mult1'. Switching to method=='indepMH'")
        return(rmultinom1Test(n, size, prob, maxSize, method="indepMH", verbose, minSample, 
                          maxExpectedSizeBeforeSwitch, init, burnIn, filterEvery, 
                          zeroProbZeroSamples, allowSizeLessThanK))
      }
      
      # sample expectedSamples times a fudge factor, but make sure we don't get past memory limit
      thisNumberOfSamples = max(minSample, min(maxSamples, expectedSamples * 1.1))
      if(verbose)
        print(paste0("Sampling ", thisNumberOfSamples, ". Sampled ", n-samplesLeft, "/", n, ". Expected remaining samples: ", expectedSamples))
      thisSamples = 1 + stats::rmultinom(thisNumberOfSamples, size-k, prob=prob)
      
      # calculate accept probabilities
      thisProbs = (size-k) / apply(thisSamples, 2, prod)
      if(verbose) {
        print(paste0("Max sampled accept prob: ", max(thisProbs), ". Mean sampled accept prob: ", mean(thisProbs)))
        print(paste0("Max theoretical accept prob: ", 1, ". Mean 'theoretical' accept prob: ", averageProb))
      }
      
      # reject relevant samples
      u = stats::runif(thisNumberOfSamples)
      thisSamples = thisSamples[,u<thisProbs]
      
      # remove excess samples if necessary
      totalSamples = ncol(thisSamples) + n - samplesLeft
      if(totalSamples > n) {
        thisSamples = thisSamples[,1:samplesLeft]
      }
      
      # add in accepted samples, if any
      if(ncol(thisSamples) > 0) {
        samples[,(n-samplesLeft+1):(n-samplesLeft+ncol(thisSamples))] = thisSamples
      } else {
        warning(paste0("no samples accepted this round out of ", thisNumberOfSamples, " total. Redrawing..."))
      }
    }
  } else if(method == "mult") {
    
    while(any(is.na(samples))) {
      # calculate the number of remaining samples
      samplesLeft = sum(apply(samples, 2, function(x) {any(is.na(x))}))
      
      # approximate expected number of samples so that, after some are rejected, we will 
      # have the right number of samples
      expectedSamples = ceiling(samplesLeft/averageProbMult)
      
      if(expectedSamples*k > maxExpectedSizeBeforeSwitch) {
        warning("too many samples expected with method=='mult'. Switching to method=='mult1'")
        return(rmultinom1Test(n, size, prob, maxSize, method="mult1", verbose, minSample, maxExpectedSizeBeforeSwitch, 
                          init, burnIn, filterEvery, 
                          zeroProbZeroSamples, allowSizeLessThanK))
      }
      
      # sample expectedSamples times a fudge factor, but make sure we don't get past memory limit
      thisNumberOfSamples = max(minSample, min(maxSamples, expectedSamples * 1.1))
      if(verbose)
        print(paste0("Sampling ", thisNumberOfSamples, ". Sampled ", n-samplesLeft, "/", n, ". Expected remaining samples: ", expectedSamples))
      thisSamples = matrix(stats::rmultinom(thisNumberOfSamples, size, prob=prob), ncol=thisNumberOfSamples)
      
      # reject relevant samples
      accept = apply(thisSamples, 2, function(x) {all(x>0)})
      thisSamples = matrix(thisSamples[,accept], nrow=length(prob))
      
      # remove excess samples if necessary
      totalSamples = ncol(thisSamples) + n - samplesLeft
      if(totalSamples > n) {
        thisSamples = matrix(thisSamples[,1:samplesLeft], nrow=length(prob))
      }
      
      # add in accepted samples, if any
      if(ncol(thisSamples) > 0) {
        samples[,(n-samplesLeft+1):(n-samplesLeft+ncol(thisSamples))] = thisSamples
      } else {
        warning(paste0("no samples accepted this round out of ", thisNumberOfSamples, " total. Redrawing..."))
      }
    }
  } else if(method == "indepMH") {
    # we use the mult1 method for independent proposals with independent Metropolis-Hastings
    # (https://www.statlect.com/fundamentals-of-statistics/Metropolis-Hastings-algorithm#hid9)
    
    # initialize at something reasonable, if not set by user
    if(is.null(init)) {
      init = 1 + floor(size*prob)
      while(sum(init) > size) {
        tooMuchBy = sum(init) - size
        numberReducible = sum(init > 1)
        reduceNumber = min(tooMuchBy, numberReducible)
        init[init > 1][1:reduceNumber] = init[init > 1][1:reduceNumber] - 1
      }
    }
    if(sum(init) != size)
      stop("sum(init) != size")
    
    # approximate target log-density
    lp <- log(prob)
    lf <- function(x) {
      if(any(x < 1) || sum(x) != size)
        return(-Inf)
      sum(lp*x - lfactorial(x))
    }
    
    # true proposal log-density
    lq = function(x) {
      if(sum(x) != size)
        return(-Inf)
      sum(lp*(x-1) - lfactorial(x-1)) + lfactorial(size-k)
    }
    
    # proposal function
    q <- function(x) {
      1 + stats::rmultinom(1, size-k, prob)
    }
    
    # do the sampling
    tmp <- init
    ar <- 0
    for (i in 1:burnIn) {
      proposal <- q(tmp)
      p <- exp((lf(proposal) - lq(proposal)) - (lf(tmp) - lq(tmp)))
      if (stats::runif(1) < p) {
        tmp <- proposal
      }
    }
    for (i in 1:ncol(samples)) {
      proposal <- q(tmp)
      p <- exp((lf(proposal) - lq(proposal)) - (lf(tmp) - lq(tmp)))
      if (stats::runif(1) < p) {
        tmp <- proposal
        ar <- ar + 1
      }
      samples[,i] <- tmp
    }
    
    # calculated acceptance percentage
    if(verbose) {
      print(paste0("acceptance percentage: ", ar/ncol(samples)))
    }
    
    # filter out samples to reduce autocorrelation
    samples = samples[,seq(from=1, to=ncol(samples), by=filterEvery)]
  }
  
  samples
}

#' @describeIn simPopInternal Take multilevel multinomial draws first from joint distribution of 
#' number of households per EA given the total per stratum, and then from the joint 
#' distribution of the total target population per household given 
#' the total per stratum
sampleNMultilevelMultinomialTest = function(nDraws = ncol(pixelIndexMat), pixelIndexMat=NULL, urbanMat=NULL, areaMat=NULL, easpaList, 
                                        popMat, stratifyByUrban=TRUE, verbose=TRUE, returnEAinfo=FALSE, 
                                        minHHPerEA=25, fixHHPerEA=NULL, fixPopPerHH=NULL) {
  
  if(length(easpaList) == 1) {
    easpaList = replicate(nDraws, easpaList[[1]], simplify=FALSE)
  }
  areas = easpaList[[1]]$area
  
  if((is.null(areaMat) || is.null(urbanMat)) && is.null(pixelIndexMat)) {
    stop("user must either supply pixelIndexMat or both areaMat and urbanMat")
  }
  if(is.null(areaMat)) {
    areaIs = match(popMat$area, sort(unique(areas))) # convert from character to indices to save memory
    areaMat = matrix(areaIs[pixelIndexMat], ncol=nDraws)
  }
  if(is.null(urbanMat)) {
    urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
  }
  
  # start by drawing the totals, then divide households amongst EAs, then divide target population amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the totals
  
  # get the total number of enumeration areas per stratum (this does not change between draws)
  areasIs = match(areas, sort(unique(areas))) # convert from character to indices
  totalEAsUrban = easpaList[[1]]$EAUrb
  totalEAsRural = easpaList[[1]]$EARur
  totalEAs = easpaList[[1]]$EATotal
  nEAs = sum(totalEAs)
  
  ##### draw the total target population per enumeration area
  
  targetPopDraws = matrix(nrow=nEAs, ncol=nDraws)
  if(returnEAinfo) {
    householdDraws = matrix(nrow=nEAs, ncol=nDraws)
  }
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum minHHPerEA)
  if(stratifyByUrban) {
    totalHouseholdsUrban = sweep(matrix(sapply(easpaList, function(x) {x$HHUrb}), ncol=length(easpaList)), 1, -minHHPerEA*totalEAsUrban, "+")
    totalHouseholdsRural = sweep(matrix(sapply(easpaList, function(x) {x$HHRur}), ncol=length(easpaList)), 1, -minHHPerEA*totalEAsRural, "+")
    totalChildrenUrban = matrix(sapply(easpaList, function(x) {x$popUrb}), ncol=length(easpaList))
    totalChildrenRural = matrix(sapply(easpaList, function(x) {x$popRur}), ncol=length(easpaList))
  } else {
    totalHouseholds = sweep(sapply(matrix(easpaList, function(x) {x$HHTotal}), ncol=length(easpaList)), 1, -minHHPerEA*totalEAs, "+")
    totalChildren = matrix(sapply(easpaList, function(x) {x$popHHTotal}), ncol=length(easpaList))
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the target population amongst the households, also with a multinomial distribution
  for(i in 1:length(areasIs)) {
    thisArea = areasIs[i]
    thisAreaL = areaMat==thisArea
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(stratifyByUrban) {
      if(totalEAsUrban[i] != 0) {
        householdDrawsUrban <- matrix(sapply(totalHouseholdsUrban[i,], stats::rmultinom, n=1, prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])), nrow=totalEAsUrban[i], ncol=nDraws) + minHHPerEA
      }
      if(totalEAsRural[i] != 0) {
        householdDrawsRural <- matrix(sapply(totalHouseholdsRural[i,], stats::rmultinom, n=1, prob=rep(1/totalEAsRural[i], totalEAsRural[i])), nrow=totalEAsRural[i], ncol=nDraws) + minHHPerEA
      }
      
      # if we must return EA info, we must return the household draws for each EA:
      if(returnEAinfo) {
        if(totalEAsUrban[i] != 0) {
          householdDraws[thisAreaL & urbanMat] = householdDrawsUrban
        }
        if(totalEAsRural[i] != 0) {
          householdDraws[thisAreaL & !urbanMat] = householdDrawsRural
        }
      }
    } else {
      householdDraws = matrix(sapply(totalHouseholds[i,], stats::rmultinom, n=1, prob=rep(1/totalEAs[i], totalEAs[i])), nrow=totalEAs[i], ncol=nDraws) + minHHPerEA
    }
    
    # drawing target population per EA
    if(is.null(fixPopPerHH)) {
      # draw target pop per EA with probability proportional to the number of households
      if(stratifyByUrban) {
        if(totalEAsUrban[i] != 0) {
          probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
          targetPopDraws[thisAreaL & urbanMat] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
        }
        
        if(totalEAsRural[i] != 0) {
          probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
          targetPopDraws[thisAreaL & !urbanMat] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
        }
      } else {
        probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
        targetPopDraws[thisAreaL] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildren[i,j], probs[,j])})
      }
    } else {
      # set target pop per EA based on fixed number per household
      if(stratifyByUrban) {
        if(totalEAsUrban[i] != 0) {
          targetPopDraws[thisAreaL & urbanMat] = fixPopPerHH*householdDrawsUrban
        }
        if(totalEAsRural[i] != 0) {
          targetPopDraws[thisAreaL & !urbanMat] = fixPopPerHH*householdDrawsRural
        }
      } else {
        targetPopDraws[thisAreaL] = fixPopPerHH*householdDraws
      }
    }
  }
  
  ##### Return results
  if(!returnEAinfo) {
    targetPopDraws
  } else {
    list(householdDraws=householdDraws, targetPopDraws=targetPopDraws)
  }
}

#' @describeIn simPopInternal Same as sampleNMultilevelMultinomialTest, except the number of EAs per pixel is fixed
sampleNMultilevelMultinomialFixedTest = function(clustersPerPixel, nDraws=ncol(pixelIndices), pixelIndices=NULL, 
                                             urbanVals=NULL, areaVals=NULL, easpa, popMat, stratifyByUrban=TRUE, 
                                             verbose=TRUE) {
  
  # set default inputs
  if((is.null(areaVals) || is.null(urbanVals)) && is.null(pixelIndices)) {
    stop("user must either supply pixelIndices or both areaVals and urbanVals")
  }
  if(is.null(areaVals)) {
    areaVals = matrix(popMat$area[pixelIndices], ncol=nDraws)
  }
  if(is.null(urbanVals)) {
    urbanVals = matrix(popMat$urban[pixelIndices], ncol=nDraws)
  }
  
  # start by drawing the totals, then divide households amongst EAs, then divide target population amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the totals
  
  # get the total number of enumeration areas per stratum (this does not change between draws)
  areas = easpa$area
  totalEAsUrban = easpa$EAUrb
  totalEAsRural = easpa$EARur
  totalEAs = easpa$EATotal
  nEAs = sum(totalEAs)
  
  if(nEAs != sum(clustersPerPixel)) {
    stop("sum(easpa$EATotal) != sum(clustersPerPixel)")
  }
  
  ##### draw the total target population per EA
  targetPopDraws = matrix(nrow=nEAs, ncol=nDraws)
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(stratifyByUrban) {
    totalHouseholdsUrban = easpa$HHUrb -25*totalEAsUrban
    totalHouseholdsRural = easpa$HHRur -25*totalEAsRural
    totalChildrenUrban = easpa$popUrb
    totalChildrenRural = easpa$popRur
  } else {
    totalHouseholds = easpa$HHTotal - 25*totalEAs
    totalChildren = easpa$popHHTotal
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the target population amongst the households, also with a multinomial distribution
  for(i in 1:length(areas)) {
    thisArea = areas[i]
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(stratifyByUrban) {
      if(totalEAsUrban[i] != 0) {
        if(any(totalHouseholdsUrban != 0)) {
          householdDrawsUrban = stats::rmultinom(n=nDraws, size=totalHouseholdsUrban[i], prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])) + 25
        } else {
          householdDrawsUrban = matrix(rep(25, totalEAsUrban[i]*nDraws), ncol=nDraws)
        }
      }
      
      if(totalEAsRural[i] != 0) {
        if(any(totalHouseholdsRural != 0)) {
          householdDrawsRural = stats::rmultinom(n=nDraws, size=totalHouseholdsRural[i], prob=rep(1/totalEAsRural[i], totalEAsRural[i])) + 25
        } else {
          householdDrawsRural = matrix(rep(25, totalEAsRural[i]*nDraws), ncol=nDraws)
        }
      }
    } else {
      if(any(totalHouseholdsRural != 0)) {
        householdDraws = stats::rmultinom(n=nDraws, size=totalHouseholds[i], prob=rep(1/totalEAs[i], totalEAs[i])) + 25
      } else {
        householdDraws = matrix(rep(25, totalEAs[i]*nDraws), ncol=nDraws)
      }
    }
    
    # drawing target population per EA, with probability proportional to the number of households
    if(stratifyByUrban) {
      
      if(totalEAsUrban[i] != 0) {
        probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
        targetPopDraws[areaVals==thisArea & urbanVals] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenUrban[i], probsUrban[,j])})
      }
      
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        targetPopDraws[areaVals==thisArea & !urbanVals] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildrenRural[i], probsRural[,j])})
      }
      
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      targetPopDraws[areaVals==thisArea] = sapply(1:nDraws, function(j) {stats::rmultinom(1, totalChildren[i], probs[,j])})
    }
  }
  
  ##### Return results
  targetPopDraws
}