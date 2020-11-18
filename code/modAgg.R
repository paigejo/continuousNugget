### Aggregation Models
## Each aggregation model must at least use the following inputs: 
# uDraws: posterior draws from the spatial field on logit scale. nLocs x nsim dimension matrix.
# sigmaEpsilonDraws: posterior draws of cluster effect SD (joint draws with uDraws). Of length nsim
# easpa: table of number of EAs per area similar to easpc
# popMat: nLoc row matrix of population densities/EA densities (not necessarily normalized) and strata 
#         for pixel grid such as from the popGrid variable. Eventually, this may include 'subarea' ID 
#         for aggregation to
# empiricalDistributions: A list of empirical distributions for the EA denominators for each urban/rural 
#                         value 
# includeUrban: if TRUE, uses EAsUrb and EAsRur in easpa.  Otherwise uses EAsTotal
## In addition, here are a few other possible inputs:
# maxEmptyFraction: the maximum fraction of posterior samples for which a given subarea can have 
#                   no EAs before separate marginal approximations are performed. If 1, never 
#                   perform separate marginal approximations
# doModifiedPixelLevel: include draws at the pixel level conditional on there being at least one EA per pixel
## NOTE1: For aggregation consistent models: 'areas' should be the finest spatial areas for which number 
##        of urban and rural EAs is known. For coarser areas, further aggregate the draws of each 
##        aggregation function as necessary. The same goes for the pixel grid. For larger 'subareas' 
##        containing an unknown number of EAs, in which case aggregate pixel predictions as need be. 
## NOTE2: Aggregation consistent models:
#           pixels/subareas (unknown EAs): LCPB, lcpb (any others?)
#           areas/regions (known EAs): all of them?

# Use the popGrid variable to construct a matrix suitable for input into the aggregation models. 
# population matrix must have the following values: 
# lon: longitude
# lat: latitude
# east: easting (km)
# north: northing (km)
# pop: proportional to population density for each grid cell
# area: an id or area name in which the grid cell corresponding to each row resides
# urban: whether the grid cell is urban or rural
# constituency: the sub-area
# province: the super-area
makeDefaultPopMat = function(getAdjusted=FALSE) {
  if(getAdjusted) {
    out = popGridAdjusted
  } else {
    out = popGrid
  }
  out$pop = out$popOrig
  out$area = out$admin1
  out$popOrig = NULL
  out$admin1 = NULL
  out$constituency = out$admin2
  out$province = out$region
  
  out
}

# Use the easpc variable to construct a table suitable for input into the aggregation models. 
# output must have the following values: 
# area: the name or id of the area
# EAUrb: the number of EAs in the urban part of the area
# EARur: the number of EAs in the rural part of the area
# EATotal: the number of EAs in the the area
# HHUrb: the number of households in the urban part of the area
# HHRur: the number of households in the rural part of the area
# HHTotal: the number of households in the the area
# popUrb: the number of people in the urban part of the area
# popRur: the number of people in the rural part of the area
# popTotal: the number of people in the the area
# Also, the rows of this table will be ordered by Area
makeDefaultEASPA = function(dataType=c("children", "women"), validationClusterI=NULL, 
                            useClustersAsEAs=!is.null(validationClusterI), usePixelUrban=TRUE) {
  dataType = match.arg(dataType)
  
  if(!useClustersAsEAs) {
    out = merge(easpc, adjustPopulationPerCountyTable(dataType))
  } else {
    if(is.null(validationClusterI)) {
      # This case will likely not be used, and corresponds with leave one county out
      out = merge(easpcMort, poppcMort)
    } else {
      ## construct faux population based on left out observations
      thisMort = mort[validationClusterI,]
      
      # change urbanicity to take the pixel level values rather than the cluster level values if necessary:
      test = getPixelIndex(cbind(thisMort$east, thisMort$north), popGrid, thisMort$admin1)
      temp = thisMort
      if(usePixelUrban) {
        temp$urban = popGrid$urban[test]
      }
      
      # faux poppc
      out = aggregate(temp$n, by=list(admin1=temp$admin1, urban=temp$urban), FUN=sum, drop=FALSE)
      out[is.na(out[,3]), 3] = 0
      # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
      out2 = cbind(out, rural=0)[48:94,]
      out2[, 4] = out$x[1:47]
      unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
      poppcMort = poppc
      poppcMort[unsortedToSorted, 2:3] = out2[,3:4]
      poppcMort$popTotal = poppcMort$popUrb + poppcMort$popRur
      poppcMort$pctTotal = 100 * poppcMort$popTotal * (1 / sum(poppcMort$popTotal))
      poppcMort$pctUrb = 100 * poppcMort$popUrb / poppcMort$popTotal
      
      # faux easpc
      out = aggregate(temp$n, by=list(admin1=temp$admin1, urban=temp$urban), FUN=length, drop=FALSE)
      out[is.na(out[,3]), 3] = 0
      # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
      out2 = cbind(out, rural=0)[48:94,]
      out2[, 4] = out$x[1:47]
      unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
      easpcMort = clustpc
      names(easpcMort) = names(easpc)
      easpcMort[unsortedToSorted, 2:3] = out2[,3:4]
      easpcMort$EATotal = easpcMort$EAUrb + easpcMort$EARur
      
      easpcMort$HHUrb = 25 * easpcMort$EAUrb
      easpcMort$HHRur = 25 * easpcMort$EARur
      easpcMort$HHTotal = easpcMort$HHUrb + easpcMort$HHRur
      
      # combine results
      out = merge(easpcMort, poppcMort)
    }
  }
  
  names(out)[1] = "area"
  
  # sort by area
  sortI = sort(out$area, index.return=TRUE)$ix
  out = out[sortI,]
  
  out
}

# Use the eaDat to construct easpa in same format as makeDefaultEASPA()
# output must have the following values: 
# area: the name or id of the area
# EAUrb: the number of EAs in the urban part of the area
# EARur: the number of EAs in the rural part of the area
# EATotal: the number of EAs in the the area
# HHUrb: the number of households in the urban part of the area
# HHRur: the number of households in the rural part of the area
# HHTotal: the number of households in the the area
# popUrb: the number of people in the urban part of the area
# popRur: the number of people in the rural part of the area
# popTotal: the number of people in the the area
# Also, the rows of this table will be ordered by Area name alphabetically
makeEASPAFromEADat = function(eaDat) {
  
  # first get the number of enumeration areas per stratum and county
  out = aggregate(eaDat$n, by=list(admin1=eaDat$admin1, urban=eaDat$urban), FUN=length, drop=FALSE)
  out[is.na(out[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  out2 = cbind(out, rural=0)[48:94,]
  out2[, 4] = out$x[1:47]
  easpa = makeDefaultEASPA()
  easpa[, 2:3] = out2[,3:4]
  easpa$EATotal = easpa$EAUrb + easpa$EARur
  
  # now get the target population count per stratum and county (and percent population total and urban)
  out = aggregate(eaDat$n, by=list(admin1=eaDat$admin1, urban=eaDat$urban), FUN=sum, drop=FALSE)
  out[is.na(out[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  out2 = cbind(out, rural=0)[48:94,]
  out2[, 4] = out$x[1:47]
  easpa[, 8:9] = out2[,3:4]
  easpa$popTotal = easpa$popUrb + easpa$popRur
  easpa$pctTotal = easpa$popTotal / sum(easpa$popTotal) * 100
  easpa$pctUrb = easpa$popUrb / easpa$popTotal * 100
  
  # now get the number of households per stratum and county
  out = aggregate(eaDat$nHH, by=list(admin1=eaDat$admin1, urban=eaDat$urban), FUN=sum, drop=FALSE)
  out[is.na(out[,3]), 3] = 0
  # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
  out2 = cbind(out, rural=0)[48:94,]
  out2[, 4] = out$x[1:47]
  easpa[, 5:6] = out2[,3:4]
  easpa$HHTotal = easpa$HHUrb + easpa$HHRur
  
  easpa
}

# in this model, we assume a binomial process for the EA locations. Follows algorithm 2 from the 
# outline
# validationPixelI: a set of indices of pixels for which we want predictions in validation.
# urbanEffect: only used if constituency level results are desired. Helps produce urban/rural predictions in 
#              constituencies without urban/rural pixels
modLCPB = function(uDraws, sigmaEpsilonDraws=NULL, easpa=NULL, popMat=NULL, adjustedPopMat=NULL, empiricalDistributions=NULL, 
                   includeUrban=TRUE, maxEmptyFraction=1, clusterLevel=TRUE, pixelLevel=TRUE, constituencyLevel=TRUE, countyLevel=TRUE, 
                   regionLevel=TRUE, nationalLevel=TRUE, doModifiedPixelLevel=TRUE, validationPixelI=NULL, validationClusterI=NULL, 
                   onlyDoModifiedPixelLevel=FALSE, clustersPerPixel=NULL, 
                   doLcpb=FALSE, doLCpb=FALSE, doLCPb=FALSE, constituencyPop=poppcon, 
                   ensureAtLeast1PerConstituency=FALSE, urbanEffect=NULL, 
                   returnEAinfo=FALSE, epsc=NULL, urbanEffectDraws=NULL) {
  nDraws = ncol(uDraws)
  
  startTime0 = proc.time()[3]
  
  # set default inputs
  if(is.null(easpa)) {
    if(!is.null(validationPixelI) || !is.null(validationClusterI))
      stop("Must include both validationClusterI and validationPixelI or easpa")
    
    easpa = makeDefaultEASPA(validationClusterI=validationClusterI, useClustersAsEAs=!is.null(validationClusterI))
    # area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
    # HHUrb: the number of households in the urban part of the area
    # HHRur: the number of households in the rural part of the area
    # HHTotal: the number of households in the the area
    # popUrb: the number of people in the urban part of the area
    # popRur: the number of people in the rural part of the area
    # popTotal: the number of people in the the area
  }
  
  totalEAs = sum(easpa$EATotal)
  if(!is.null(clustersPerPixel)) {
    emptyPixels = clustersPerPixel == 0
    if(totalEAs != sum(clustersPerPixel))
      stop("sum(easpa$EATotal) != sum(clustersPerPixel)")
    if(doModifiedPixelLevel)
      stop("doModifiedPixelLevel cannot be set to TRUE when clustersPerPixel is specified")
  }
  
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  if(is.null(empiricalDistributions)) {
    out = load(paste0(globalDirectory, "empiricalDistributions.RData"))
    # list(households=householdDistribution, mothers=motherDistribution, children=childrenDistribution,
    #      householdsUrban=householdDistributionUrban, mothersUrban=motherDistributionUrban, childrenUrban=childrenDistributionUrban,
    #      householdsRural=householdDistributionRural, mothersRural=motherDistributionRural, childrenRural=childrenDistributionRural)
    # May also have popUrban, popRural, and popTotal, in which case they are used
    
    # make edfun versions of the ecdfs
    fastDistributions = list()
    for(i in 1:length(empiricalDistributions)) {
      fastDistributions = c(fastDistributions, list(ecdf2edfun(empiricalDistributions[[i]])))
    }
    names(fastDistributions) = names(empiricalDistributions)
  }
  
  if(!is.null(urbanEffect)) {
    if(!is.null(urbanEffectDraws) && any(urbanEffectDraws != urbanEffect)) {
      stop("noncompatible urbanEffect and urbanEffectDraws given as input. Use either one or the other")
    }
    urbanEffectDraws = rep(urbanEffect, ncol(uDraws))
  }
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  finishedSetupTime1 = proc.time()[3]
  
  ##### Line 1: take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  pixelProbs = popMat$pop
  
  # take draws from the stratified binomial process for each posterior sample
  if(is.null(clustersPerPixel)) {
    if(!onlyDoModifiedPixelLevel) {
      if(constituencyLevel) {
        eaSamples = rStratifiedMultnomialByConstituency(nDraws, popMat, easpa, includeUrban, constituencyPop=constituencyPop, 
                                                        ensureAtLeast1PerConstituency=ensureAtLeast1PerConstituency)
      } else {
        eaSamples = rStratifiedMultnomial(nDraws, popMat, easpa, includeUrban)
      }
    }
    if(doModifiedPixelLevel) {
      if(constituencyLevel) 
        stop("constituencyLevel and doModifiedPixelLevel both being set to TRUE is not yet supported")
      eaSamplesMod = rStratifiedBinomial1(nDraws, popMat, easpa, includeUrban, validationPixelI=validationPixelI)
    }
  }
  
  if(!is.null(clustersPerPixel) && !exists("eaSamples")) {
    eaSamples = matrix(rep(clustersPerPixel, nDraws), ncol=nDraws)
  }
  
  # make matrix (or list) of pixel indices mapping matrices of EA values to matrices of pixel values
  if(!is.null(clustersPerPixel)) {
    pixelIndices = rep(1:nrow(popMat), times=clustersPerPixel) # this contains repetitions and has length == nEAs
    uniquePixelIndices = sort(unique(pixelIndices))
  } else if(!onlyDoModifiedPixelLevel) {
    pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  } else if(doModifiedPixelLevel) {
    if(is.null(validationPixelI)) {
      pixelIndexListMod = lapply(1:nDraws, function(j) {rep(1:nrow(popMat), times=eaSamplesMod[,j])})
    } else {
      pixelIndexListMod = lapply(1:nDraws, function(j) {rep(sort(validationPixelI), times=eaSamplesMod[,j])})
    }
  }
  
  # determine which EAs are urban if necessary
  if(includeUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    if(!is.null(clustersPerPixel)) {
      urbanVals = popMat$urban[pixelIndices]
      uniqueUrbanVals = popMat$urban[uniquePixelIndices]
    }else if(!onlyDoModifiedPixelLevel) {
      urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
    } else if(doModifiedPixelLevel) {
      urbanListMod = lapply(1:nDraws, function(j) {popMat$urban[pixelIndexListMod[[j]]]})
    }
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  if(!is.null(clustersPerPixel)) {
    areaVals = popMat$area[pixelIndices]
    uniqueAreaVals = popMat$area[uniquePixelIndices]
  } else if(!onlyDoModifiedPixelLevel) {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  } else if(doModifiedPixelLevel) {
    areaListMod = lapply(1:nDraws, function(j) {popMat$area[pixelIndexListMod[[j]]]})
  }
  
  finishedDrawingEAsTime2 = proc.time()[3]
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  if(!onlyDoModifiedPixelLevel || !is.null(clustersPerPixel)) {
    if(is.null(epsc)) {
      epsc = matrix(rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
    }
  } else if(doModifiedPixelLevel) {
    epscMod = lapply(1:nDraws, function(j) {rnorm(length(areaListMod[[j]]), sd=sigmaEpsilonDraws[j])})
  }
  
  finishedDrawingClusterEffectsTime3 = proc.time()[3]
  
  ##### Line 3: draw EA population denominators, N
  # NOTE1: we assume the values for each EA within urban/rural boundaries are iid
  # if(includeUrban) {
  #   load(paste0(globalDirectory, "NcsUrbanRural.RData"))
  #   
  #   Ncs = matrix(0, nrow=nrow(epsc), ncol=ncol(epsc))
  #   Ncs[urbanMat] = NcsUrban
  #   Ncs[!urbanMat] = NcsRural
  # } else {
  #   load(paste0(globalDirectory, "Ncs.RData"))
  # }
  if(!is.null(clustersPerPixel)) {
    if(is.null(validationPixelI))
      stop("clustersPerPixel must only be set for validation, but validationPixelI is NULL")
    
    # in this case, every left out cluster has exactly 25 households. Simply sample children 
    # with equal probability from each cluster/faux EA
    Ncs = sampleNPoissonMultinomialFixed(clustersPerPixel, nDraws=nDraws, pixelIndices=pixelIndices, 
                                         urbanVals=urbanVals, areaVals=areaVals, easpa=easpa, popMat=popMat, includeUrban=includeUrban, 
                                         verbose=TRUE)
  } else if(!onlyDoModifiedPixelLevel) {
    if(returnEAinfo) {
      out = sampleNPoissonMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                popMat=popMat, includeUrban=includeUrban, verbose=TRUE, returnEAinfo=returnEAinfo)
      householdDraws = out$householdDraws
      Ncs = out$childrenDraws
    } else {
      Ncs <- sampleNPoissonMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                       popMat=popMat, includeUrban=includeUrban, verbose=TRUE, returnEAinfo=returnEAinfo)
      householdDraws = NULL
    }
  } else if(doModifiedPixelLevel) {
    NcsMod = sampleNPoissonBinomial(eaSamplesMod=eaSamplesMod, pixelIndexListMod=pixelIndexListMod, areaListMod=areaListMod, urbanListMod=urbanListMod, 
                                    includeUrban=includeUrban, easpa=easpa, popMat=popMat)
  }
  
  finishedDrawingNsTime4 = proc.time()[3]
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  if(!is.null(clustersPerPixel)) {
    uc = uDraws[pixelIndices,]
    muc = expit(uc + epsc)
  } else if(!onlyDoModifiedPixelLevel) {
    uc = matrix(uDraws[cbind(rep(rep(1:nrow(uDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
    muc = expit(uc + epsc)
  }
  if(doModifiedPixelLevel) {
    if(is.null(validationPixelI)) {
      ucMod = lapply(1:nDraws, function(j) {uDraws[rep(1:nrow(uDraws), times=eaSamplesMod[,j]), j]})
    } else {
      ucMod = lapply(1:nDraws, function(j) {uDraws[rep(sort(validationPixelI), times=eaSamplesMod[,j]), j]})
    }
    
    mucMod = lapply(1:nDraws, function(j) {expit(ucMod[[j]] + epscMod[[j]])})
  }
  
  finishedCalculatingMusTime5 = proc.time()[3]
  
  # calculate Z_{ic} for each EA in each pixel
  if(!onlyDoModifiedPixelLevel || !is.null(clustersPerPixel)) {
    Zc = matrix(rbinom(n=totalEAs * nDraws, size=Ncs, prob=as.matrix(muc)), ncol=nDraws)
  }
  if(doModifiedPixelLevel) {
    ZcMod = lapply(1:nDraws, function(j) {rbinom(n=length(mucMod[[j]]), size=NcsMod[[j]], prob=mucMod[[j]])})
  }
  
  finishedDrawingZsTime6 = proc.time()[3]
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, vals, applyFun=sum, doMod=FALSE, popWeightMatrix=NULL) {
    # calculate levels over which to aggregate
    if(!doMod) {
      if(!is.null(clustersPerPixel)) {
        indices = pixelIndices
      } else {
        indices = factor(as.character(pixelIndexMat[,i]))
      }
    } else {
      indices = factor(as.character(pixelIndexListMod[[i]]))
    }
    
    # in this case (the LCPb model), we calculate weighted means within factor levels using popWeightMatrix
    if(!is.null(popWeightMatrix)) {
      stop("using popWeightMatrix is no longer support, since this is much slower than calculating normalized 
           weights separately and multiplying values by them outside this function")
      # applyFun = function(x) {weighted.mean(x, popWeightMatrix[,i], na.rm=TRUE)}
      
      Data = data.frame(v=vals[,i], w=popWeightMatrix[,i])
      out = sapply(split(Data, indices), function(x) weighted.mean(x$v,x$w))
    } else {
      if(!doMod) {
        if(!is.null(clustersPerPixel)) {
          # out = tapply(vals[,i], factor(as.character(pixelIndices)), FUN=applyFun)
          out = tapply(vals[,i], indices, FUN=applyFun)
        } else {
          out = tapply(vals[,i], indices, FUN=applyFun)
        }
      } else {
        # in this case
        out = tapply(vals[[i]], indices, FUN=applyFun)
      }
    }
    
    if(!is.null(clustersPerPixel)) {
      returnValues = out
    } else {
      indices = as.numeric(names(out))
      
      returnValues = rep(NA, nrow(uDraws))
      returnValues[indices] = out
    }
    
    returnValues
  }
  
  ##### Line 5: We already did this, resulting in uDraws input
  
  ##### Line 6: aggregate population denominators for each grid cell to get N_{ig}
  if(!onlyDoModifiedPixelLevel || !is.null(clustersPerPixel)) {
    Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, vals=Ncs)
    Ng[is.na(Ng)] = 0
  }
  if(doModifiedPixelLevel) {
    NgMod <- sapply(1:nDraws, getPixelColumnFromEAs, vals=NcsMod, doMod=TRUE)
    NgMod[is.na(NgMod)] = 0
  }
  
  finishedPixelAggregationNsTime7 = proc.time()[3]
  
  ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  if(!onlyDoModifiedPixelLevel || !is.null(clustersPerPixel)) {
    Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=Zc)
  }
  if(doModifiedPixelLevel) {
    ZgMod <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=ZcMod, doMod=TRUE)
  }
  
  finishedPixelAggregationZsTime8 = proc.time()[3]
  
  ##### Line 8: Calculate empirical mortality proportions for each grid cell, p_{ig}. 
  #####         Whenever Z_{ig} is 0, set p_{ig} to 0 as well
  if(!onlyDoModifiedPixelLevel) {
    pg = Zg / Ng
    pg[Zg == 0] = 0
  }
  
  if(doModifiedPixelLevel) {
    pgMod = ZgMod / NgMod
    pgMod[is.na(ZgMod)] = 0
  }
  
  finishedPixelAggregationPsTime9 = proc.time()[3]
  
  ##### calculate results for the other models if necessary
  
  # uses the logistic approximation for speedup. Even if dolcpb is FALSE we still use this for calculating the mean (central predictions)
  lcpbSwitchedUrban = NULL
  if(is.null(clustersPerPixel)) {
    lcpb = matrix(logitNormMean(cbind(c(as.matrix(uDraws)), rep(sigmaEpsilonDraws, each=nrow(uDraws)))), nrow=nrow(uDraws))
    if(constituencyLevel && !is.null(urbanEffectDraws)) {
      lcpbSwitchedUrban = lcpb
      zeroConUrban = constituencyPop$popUrb == 0
      zeroConRural = constituencyPop$popRur == 0
      allZeroCon = constituencyPop$Constituency[zeroConUrban | zeroConRural]
      uDrawsSwitchedUrban = uDraws
      uDrawsSwitchedUrban[(popMat$admin2 %in% allZeroCon) & popMat$urban,] = sweep(uDrawsSwitchedUrban[(popMat$admin2 %in% allZeroCon) & popMat$urban,], 2, urbanEffectDraws, "-")
      uDrawsSwitchedUrban[(popMat$admin2 %in% allZeroCon) & !popMat$urban,] = sweep(uDrawsSwitchedUrban[(popMat$admin2 %in% allZeroCon) & !popMat$urban,], 2, urbanEffectDraws, "+")
      lcpbSwitchedUrban[popMat$admin2 %in% allZeroCon,] = matrix(logitNormMean(cbind(c(as.matrix(uDrawsSwitchedUrban[popMat$admin2 %in% allZeroCon,])), rep(sigmaEpsilonDraws, each=sum(popMat$admin2 %in% allZeroCon)))), nrow=sum(popMat$admin2 %in% allZeroCon))
    }
  } else {
    lcpb = matrix(logitNormMean(cbind(c(as.matrix(uDraws)[uniquePixelIndices,]), rep(sigmaEpsilonDraws, each=length(uniquePixelIndices)))), nrow=length(uniquePixelIndices))
  }
  
  finishedClusterIntegrationTime10 = proc.time()[3]
  
  if(doLcpb) {
    # the only difference between this model and lcpb is in the areal predictions assuming pixels are sufficiently small. 
    # For those, we must aggregate lcpb predictions with weights proportional to the number of sampled EAs per pixel. 
    # No steps necessary for pixel level predictions though
    Lcpb = lcpb
  } else {
    Lcpb = NULL
  }
  
  if(doLCpb) {
    LCpb = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
    LCpb[!is.finite(LCpb)] = NA
    
    if(doModifiedPixelLevel) {
      LCpbMod = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=mucMod, applyFun=function(x) {mean(x, na.rm=TRUE)}, doMod=TRUE)
      LCpbMod[!is.finite(mugMod)] = NA
    }
  } else {
    LCpb = NULL
  }
  
  finishedPixelAggregationLCpbTime11 = proc.time()[3]
  
  if(doLCPb) {
    # first make population weights matrix
    popWeightMat = Ncs
    
    if(!is.null(clustersPerPixel)) {
      pixelTotals = apply(Ng, 2, function(x) {rep(x, clustersPerPixel[clustersPerPixel != 0])})
    } else {
      pixelTotals = sapply(1:nDraws, function(i) {rep(Ng[,i], eaSamples[,i])})
    }
    if(doModifiedPixelLevel) {
      pixelTotalsMod = lapply(1:nDraws, function(i) {rep(Ng[,i], eaSamplesMod[[i]])})
    } 
    
    popWeightMat = popWeightMat / pixelTotals
    
    # now calculate weighted average
    LCPb = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc*popWeightMat)
    LCPb[!is.finite(LCPb)] = NA
    
    if(doModifiedPixelLevel) {
      temp = lapply(1:nDraws, function(i) {mucMod[[i]] * pixelTotalsMod[[i]]})
      LCPbMod = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=temp, doMod=TRUE)
      LCPbMod[!is.finite(mugMod)] = NA
    }
  } else {
    LCPb = NULL
  }
  
  finishedPixelAggregationLCPbTime12 = proc.time()[3]
  
  # just for testing purposes:
  if(FALSE) {
    mug = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
    mug[!is.finite(mug)] = NA
    
    if(doModifiedPixelLevel) {
      mugMod = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=mucMod, applyFun=function(x) {mean(x, na.rm=TRUE)}, doMod=TRUE)
      mugMod[!is.finite(mugMod)] = NA
    }
    
    # this takes a veeeeery long time. Just use the logistic approximation instead
    lcpb = matrix(logitNormMean(cbind(c(as.matrix(uDraws)), rep(sigmaEpsilonDraws, each=nrow(uDraws)))), nrow=nrow(uDraws))
    
    # load shape files for plotting
    require(maptools)
    # regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
    out = load("../LK-INLA/regionMap.RData")
    out = load("../U5MR/adminMapData.RData")
    kenyaMap = adm0
    countyMap = adm1
    
    meanCols=makeRedBlueDivergingColors(64, rev = TRUE)
    sdCols=makeBlueYellowSequentialColors(64)
    popCols=makePurpleYellowSequentialColors(64, rev=TRUE)
    urbCols=makeGreenBlueSequentialColors(64)
    
    # plot means
    png(paste0(figDirectory, "exploratoryAnalysis/pixelMean.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    meanRange = quantile(probs=c(.005, .995), c(rowMeans(mug, na.rm=TRUE), rowMeans(pg, na.rm=TRUE), rowMeans(lcpb, na.rm=TRUE)), na.rm=TRUE)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(lcpb, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (lcpb model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(mug, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCpb model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(pg, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCPB model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    meanTicks = pretty(meanRange, n=5)[-1]
    meanTickLabels = as.character(meanTicks)
    image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
    
    # compare means between pixel level sampling methods
    png(paste0(figDirectory, "exploratoryAnalysis/pixelMeanMod.png"), width=1500, height=1200)
    par(mfrow=c(2,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    meanRange = quantile(probs=c(.005, .995), c(rowMeans(mug, na.rm=TRUE), rowMeans(mugMod, na.rm=TRUE), 
                                                rowMeans(pg, na.rm=TRUE), rowMeans(pgMod, na.rm=TRUE), 
                                                rowMeans(lcpb, na.rm=TRUE)), na.rm=TRUE)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(lcpb, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (lcpb model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(mug, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCpb model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(mugMod, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCpb model, mod.)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    meanTicks = pretty(meanRange, n=5)[-1]
    meanTickLabels = as.character(meanTicks)
    image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    
    plot.new() # empty plot on bottom left
    
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(pg, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCPB model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(pgMod, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCPB model, mod.)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    meanTicks = pretty(meanRange, n=5)[-1]
    meanTickLabels = as.character(meanTicks)
    image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
    
    # plot SDs
    SDsUg = apply(lcpb, 1, sd, na.rm=TRUE)
    SDsMug = apply(mug, 1, sd, na.rm=TRUE)
    SDsPg = apply(pg, 1, sd, na.rm=TRUE)
    allSDs = c(SDsUg, SDsMug, SDsPg)
    sdRange = quantile(probs=c(.001, .995), allSDs, na.rm=TRUE)
    png(paste0(figDirectory, "exploratoryAnalysis/pixelSD.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    quilt.plot(popMat$lon, popMat$lat, log(SDsUg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (lcpb model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsMug), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (LCpb model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsPg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (LCPB model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    sdTicks = pretty(sdRange, n=5)[-1]
    sdTickLabels = as.character(sdTicks)
    image.plot(zlim=range(log(sdRange)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
    
    # compare SDs between pixel level sampling methods
    SDsUg = apply(lcpb, 1, sd, na.rm=TRUE)
    SDsMug = apply(mug, 1, sd, na.rm=TRUE)
    SDsMugMod = sqrt(rowMeans(sweep(mugMod, 1, rowMeans(lcpb), "-")^2))
    SDsPg = apply(pg, 1, sd, na.rm=TRUE)
    SDsPgMod = sqrt(rowMeans(sweep(pgMod, 1, rowMeans(lcpb), "-")^2))
    allSDs = c(SDsUg, SDsMug, SDsPg, SDsMugMod, SDsPgMod)
    sdRangeMod = quantile(probs=c(.001, 1), allSDs, na.rm=TRUE)
    png(paste0(figDirectory, "exploratoryAnalysis/pixelSDMod.png"), width=1500, height=1200)
    par(mfrow=c(2,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    quilt.plot(popMat$lon, popMat$lat, log(SDsUg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRangeMod), nx=160, ny=160, main=TeX("Posterior SD (lcpb model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, log(SDsMug), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRangeMod), nx=160, ny=160, main=TeX("Posterior SD (LCpb model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, log(SDsMugMod), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRangeMod), nx=160, ny=160, main=TeX("Posterior SD (LCpb model, mod.)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    sdTicks = pretty(sdRangeMod, n=5)[-1]
    sdTickLabels = as.character(sdTicks)
    image.plot(zlim=range(log(sdRangeMod)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    
    plot.new() # empty plot on bottom left
    
    quilt.plot(popMat$lon, popMat$lat, log(SDsPg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRangeMod), nx=160, ny=160, main=TeX("Posterior SD (LCPB model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, log(SDsPgMod), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRangeMod), nx=160, ny=160, main=TeX("Posterior SD (LCPB model, mod.)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    sdTicks = pretty(sdRangeMod, n=5)[-1]
    sdTickLabels = as.character(sdTicks)
    image.plot(zlim=range(log(sdRangeMod)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
    
    # pair plot of SDs
    SDsUg = apply(lcpb, 1, sd, na.rm=TRUE)
    SDsMug = apply(mug, 1, sd, na.rm=TRUE)
    SDsMugMod = sqrt(rowMeans(sweep(mugMod, 1, rowMeans(lcpb), "-")^2))
    SDsPg = apply(pg, 1, sd, na.rm=TRUE)
    SDsPgMod = sqrt(rowMeans(sweep(pgMod, 1, rowMeans(lcpb), "-")^2))
    allSDs = cbind(SDsUg, SDsMug, SDsPg, SDsMugMod, SDsPgMod)
    sdRangeMod = quantile(probs=c(.001, 1), allSDs, na.rm=TRUE)
    
    my_line <- function(x,y,...){
      abline(a = 0,b = 1,...)
      points(x,y,..., col=theseCols)
    }
    
    numberModels = 3
    width = 3 * numberModels
    valMat = cbind(SDsUg, SDsMugMod, SDsPgMod)
    theseCols = rep(urbCols[1], nrow(valMat))
    theseCols[popMat$urban] = urbCols[64]
    zlim = quantile(probs=c(.001, 1), valMat, na.rm=TRUE)
    lims = rep(list(zlim), numberModels)
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelSDPairs.pdf"), width=width, height=width)
    
    lims = rep(list(zlim), numberModels)
    modelNames = c("lcpb", "LCpb", "LCPB")
    myPairs(valMat, 
            modelNames, 
            pch=19, cex=.2, lower.panel=my_line, upper.panel = my_line, 
            main=paste0("Pixel level posterior SD comparisons"), 
            lims=lims, oma=c(3,3,6,7))
    legend("topleft", c("Urban", "Rural"), col=c(urbCols[64], urbCols[1]), pch=19)
    dev.off()
    
    mean(SDsMugMod / SDsUg) # 1.417688
    mean(SDsPgMod / SDsUg) # 3.464197
    mean(SDsPgMod / SDsMugMod) # 2.403102
    
    mean(SDsMugMod[popMat$urban] / SDsUg[popMat$urban]) # 1.030838
    mean(SDsPgMod[popMat$urban] / SDsUg[popMat$urban]) # 1.378257
    mean(SDsPgMod[popMat$urban] / SDsMugMod[popMat$urban]) # 1.32478
    
    mean(SDsMugMod[!popMat$urban] / SDsUg[!popMat$urban]) # 1.425775
    mean(SDsPgMod[!popMat$urban] / SDsUg[!popMat$urban]) # 3.507805
    mean(SDsPgMod[!popMat$urban] / SDsMugMod[!popMat$urban]) # 2.425645
    
    # plot Ns, pop density
    png(paste0(figDirectory, "exploratoryAnalysis/pixelN.png"), width=1500, height=600)
    par(mfrow=c(1,3))
    quilt.plot(popMat$lon, popMat$lat, popMat$pop, FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("log_{10} Population surface"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, rowMeans(eaSamples), FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("log_{10} Avg. Num. EAs (min 20/1000)"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3, zlim=c(log10(.1+20/1000), 3.1))
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, rowMeans(Ng), FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("Posterior Mean of $log_{10}(N)$"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    dev.off()
    
    # compare Ns and pop density between pixel level sampling methods
    densityRange = log10(.1+range(c(popMat$pop, rowMeans(Ng), rowMeans(NgMod))))
    png(paste0(figDirectory, "exploratoryAnalysis/pixelNMod.png"), width=1500, height=1200)
    par(mfrow=c(2,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    quilt.plot(popMat$lon, popMat$lat, popMat$pop, FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("log_{10} Population surface"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3, zlim=densityRange)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, rowMeans(eaSamples), FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("log_{10} Avg. Num. EAs (min 20/1000)"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3, zlim=c(log10(.1+20/1000), 3.3))
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, rowMeans(eaSamplesMod), FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("log_{10} Avg. Num. EAs (min 20/1000, mod.)"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3, zlim=c(log10(.1+20/1000), 3.3))
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    plot.new() # empty plot on bottom left
    
    quilt.plot(popMat$lon, popMat$lat, rowMeans(Ng), FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("Posterior Mean of $log_{10}(N)$"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3, zlim=densityRange)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    quilt.plot(popMat$lon, popMat$lat, rowMeans(NgMod), FUN=function(x){log10(mean(x, na.rm=TRUE)+.1)}, 
               nx=160, ny=160, main=TeX("Posterior Mean of $log_{10}(N)$ (mod.)"), cex.main=3, col=popCols, 
               legend.cex=4, legend.width=3, zlim=densityRange)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    dev.off()
    
    
    
    ## new plots: testing all models
    cols = rep("green", length(uniquePixelIndices))
    cols[urbanPixels] = "blue"
    my_line <- function(x,y,...){
      abline(a = 0,b = 1,...)
      points(x,y,..., col=cols)
    }
    testSDs = data.frame(apply(lcpb, 1, sd), apply(Lcpb, 1, sd), apply(LCpb, 1, sd), apply(LCPb, 1, sd), apply(pg, 1, sd))
    
    numberModels = 5
    width = 3 * numberModels
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelSDPairsAllExtended.pdf"), width=width, height=width)
    myPairs(testSDs, pch=19, cex=.2, labels=c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB"), 
            lower.panel=my_line, upper.panel = my_line, ylim=range(c(testSDs), na.rm=TRUE), 
            xlim=range(c(testSDs), na.rm=TRUE), main="Prediction SD (Pixel Level)")
    dev.off()
    
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelSDPairsAll.pdf"), width=width, height=width)
    myPairs(testSDs, pch=19, cex=.2, labels=c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB"), 
            lower.panel=my_line, upper.panel = my_line, main="Prediction SD (Pixel Level)")
    dev.off()
    
    names(testSDs) = c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB")
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelSDBoxplotAll.pdf"), width=5, height=5)
    boxplot(testSDs, col="skyblue", pch=19, cex=.3, main="Predictive SDs (Pixel Level)", ylim=c(0, max(testSDs, na.rm=TRUE)))
    dev.off()
    
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelVarPairsAllExtended.pdf"), width=width, height=width)
    testVars = data.frame(apply(lcpb, 1, var), apply(Lcpb, 1, var), apply(LCpb, 1, var), apply(LCPb, 1, var), apply(pg, 1, var))
    myPairs(testVars, pch=19, cex=.2, labels=c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB"), 
            lower.panel=my_line, upper.panel = my_line, ylim=range(c(testVars), na.rm=TRUE), 
            xlim=range(c(testVars), na.rm=TRUE), main="Prediction Variance (Pixel Level)")
    dev.off()
    
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelVarPairsAll.pdf"), width=5, height=5)
    myPairs(testVars, pch=19, cex=.2, labels=c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB"), 
            lower.panel=my_line, upper.panel = my_line, main="Prediction Variance (Pixel Level)")
    dev.off()
    
    names(testVars) = c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB")
    pdf(file=paste0(figDirectory, "exploratoryAnalysis/pixelVarBoxplotAll.pdf"), width=5, height=5)
    boxplot(testVars, col="skyblue", pch=19, cex=.3, main="Predictive Variances (Pixel Level)", ylim=c(0, max(testVars, na.rm=TRUE)))
    dev.off()
  }
  
  ##### Line 9: done aggregating to the pixel level
  
  ##### Line 10: aggregate from pixel level to relevant areal levels if necessary
  if(!onlyDoModifiedPixelLevel) {
    # subset eaSamples matrix to pixels of interest if need be
    if(!is.null(clustersPerPixel)) {
      eaSamples = eaSamples[uniquePixelIndices,]
    }
    
    # lcpb model (always do this, since all models agree in expectation, and this has lowest variance)
    if(onlyDoModifiedPixelLevel && (constituencyLevel || countyLevel || regionLevel || nationalLevel)) {
      stop("onlyDoModifiedPixelLevel must be set to FALSE if dolcpb is set to TRUE and aggregations coarser than pixel level is desired")
    }
    
    # construct population density surface, adjusted based on known or estimated stratum populations (or "faux" stratum populations)
    if(is.null(adjustedPopMat)) {
      adjustedPopMat = adjustPopGrid(popMat, poppcAdjusted=easpa)
    }
    
    finishedArealAggregationSetupTime13 = proc.time()[3]
    
    # LCPB model
    aggregatedResults = aggregatePixelPredictions(Zg, Ng, popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
                                                  constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
                                                  separateUrbanRural=TRUE)
    
    finishedArealAggregationLCPBTime14 = proc.time()[3]
    
    # aggregatedResultslcpb = aggregatePixelPredictions(lcpb*eaSamples, eaSamples, popGrid=popMat, useDensity=FALSE, 
    #                                                   countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
    #                                                   separateUrbanRural=TRUE, normalize=FALSE)
    # aggregatedResultslcpb = aggregatePixelPredictions(lcpb, eaSamples, popMatAdjusted=adjustedPopMat, useDensity=TRUE, 
    #                                                   constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
    #                                                   separateUrbanRural=TRUE, normalize=TRUE, lcpbSwitchedUrban=lcpbSwitchedUrban)
    aggregatedResultslcpb = aggregatePixelPredictions(sweep(lcpb, 1, adjustedPopMat$pop, "*"), matrix(rep(adjustedPopMat$pop, nDraws), ncol=nDraws), popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
                                                      constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
                                                      separateUrbanRural=TRUE, normalize=FALSE, lcpbSwitchedUrban=lcpbSwitchedUrban)
    
    finishedArealAggregationlcpbTime15 = proc.time()[3]
    
    # Lcpb model
    if(doLcpb) {
      # unlike for the lcpb model, we weight by the number of EAs per pixel rather than population density. 
      # 
      
      # aggregatedResultsLcpb = aggregatePixelPredictions(Lcpb*eaSamples, eaSamples, popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
      #                                                   constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
      #                                                   separateUrbanRural=TRUE, normalize=FALSE)
      
      # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
      nPerEA = getExpectedNperEA(easpa, adjustedPopMat)
      thisNSamples = sweep(eaSamples, 1, nPerEA, "*")
      aggregatedResultsLcpb = aggregatePixelPredictions(Lcpb*thisNSamples, thisNSamples, popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
                                                        constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
                                                        separateUrbanRural=TRUE, normalize=FALSE)
    } else {
      aggregatedResultsLcpb = NULL
    }
    
    finishedArealAggregationLcpbTime16 = proc.time()[3]
    
    # LCpb model
    if(doLCpb) {
      # aggregatedResultsLCpb = aggregatePixelPredictions(LCpb*eaSamples, eaSamples, popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
      #                                                   constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
      #                                                   separateUrbanRural=TRUE, normalize=FALSE, lcpbSwitchedUrban=lcpbSwitchedUrban)
      
      # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
      nPerEA = getExpectedNperEA(easpa, adjustedPopMat)
      thisNSamples = sweep(eaSamples, 1, nPerEA, "*")
      aggregatedResultsLCpb = aggregatePixelPredictions(LCpb*thisNSamples, thisNSamples, popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
                                                        constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
                                                        separateUrbanRural=TRUE, normalize=FALSE, lcpbSwitchedUrban=lcpbSwitchedUrban)
    } else {
      aggregatedResultsLCpb = NULL
    }
    
    finishedArealAggregationLCpbTime17 = proc.time()[3]
    
    # LCPb model
    if(doLCPb) {
      aggregatedResultsLCPb = aggregatePixelPredictions(LCPb*Ng, Ng, popMatAdjusted=adjustedPopMat, useDensity=FALSE, 
                                                        constituencyLevel=constituencyLevel, countyLevel=countyLevel, regionLevel=regionLevel, nationalLevel=nationalLevel, 
                                                        separateUrbanRural=TRUE, normalize=FALSE, lcpbSwitchedUrban=lcpbSwitchedUrban)
    } else {
      aggregatedResultsLCPb = NULL
    }
    
    finishedArealAggregationLCPbTime18 = proc.time()[3]
    
    # determine which constituencies don't have urban or rural populations
    if(constituencyLevel) {
      zeroConstituenciesUrban = apply(aggregatedResults$constituencyMatrices$NUrban, 1, function(x) {all(x == 0)})
      zeroConstituenciesRural = apply(aggregatedResults$constituencyMatrices$NRural, 1, function(x) {all(x == 0)})
    }
    
    # use risk predictions rather than prevalence for constituency strata if there is no population denominator
    if(constituencyLevel && (any(zeroConstituenciesUrban) || any(zeroConstituenciesRural))) {
      # LCPB
      aggregatedResults$constituencyMatrices$pUrban[zeroConstituenciesUrban,] = aggregatedResultslcpb$constituencyMatrices$pUrban[zeroConstituenciesUrban,]
      aggregatedResults$constituencyMatrices$pRural[zeroConstituenciesRural,] = aggregatedResultslcpb$constituencyMatrices$pRural[zeroConstituenciesRural,]
      
      # Lcpb
      if(doLcpb) {
        aggregatedResultsLcpb$constituencyMatrices$pUrban[zeroConstituenciesUrban,] = aggregatedResultslcpb$constituencyMatrices$pUrban[zeroConstituenciesUrban,]
        aggregatedResultsLcpb$constituencyMatrices$pRural[zeroConstituenciesRural,] = aggregatedResultslcpb$constituencyMatrices$pRural[zeroConstituenciesRural,]
      }
      
      # LCpb
      if(doLCpb) {
        aggregatedResultsLCpb$constituencyMatrices$pUrban[zeroConstituenciesUrban,] = aggregatedResultslcpb$constituencyMatrices$pUrban[zeroConstituenciesUrban,]
        aggregatedResultsLCpb$constituencyMatrices$pRural[zeroConstituenciesRural,] = aggregatedResultslcpb$constituencyMatrices$pRural[zeroConstituenciesRural,]
      }
      
      # LCPb
      if(doLCPb) {
        aggregatedResultsLCPb$constituencyMatrices$pUrban[zeroConstituenciesUrban,] = aggregatedResultslcpb$constituencyMatrices$pUrban[zeroConstituenciesUrban,]
        aggregatedResultsLCPb$constituencyMatrices$pRural[zeroConstituenciesRural,] = aggregatedResultslcpb$constituencyMatrices$pRural[zeroConstituenciesRural,]
      }
    }
    
    finishedArealAggregationRiskSubstitutionTime19 = proc.time()[3]
    
    # browser()
    
    # constituency level
    Zconstituency = aggregatedResults$constituencyMatrices$Z
    Nconstituency = aggregatedResults$constituencyMatrices$N
    Pconstituency = aggregatedResults$constituencyMatrices$p
    ZconstituencyUrban = aggregatedResults$constituencyMatrices$ZUrban
    NconstituencyUrban = aggregatedResults$constituencyMatrices$NUrban
    PconstituencyUrban = aggregatedResults$constituencyMatrices$pUrban
    ZconstituencyRural = aggregatedResults$constituencyMatrices$ZRural
    NconstituencyRural = aggregatedResults$constituencyMatrices$NRural
    PconstituencyRural = aggregatedResults$constituencyMatrices$pRural
    
    # county level
    Zcounty = aggregatedResults$countyMatrices$Z
    Ncounty = aggregatedResults$countyMatrices$N
    Pcounty = aggregatedResults$countyMatrices$p
    ZcountyUrban = aggregatedResults$countyMatrices$ZUrban
    NcountyUrban = aggregatedResults$countyMatrices$NUrban
    PcountyUrban = aggregatedResults$countyMatrices$pUrban
    ZcountyRural = aggregatedResults$countyMatrices$ZRural
    NcountyRural = aggregatedResults$countyMatrices$NRural
    PcountyRural = aggregatedResults$countyMatrices$pRural
    
    # province level
    Zregion = aggregatedResults$regionMatrices$Z
    Nregion = aggregatedResults$regionMatrices$N
    Pregion = aggregatedResults$regionMatrices$p
    ZregionUrban = aggregatedResults$regionMatrices$ZUrban
    NregionUrban = aggregatedResults$regionMatrices$NUrban
    PregionUrban = aggregatedResults$regionMatrices$pUrban
    ZregionRural = aggregatedResults$regionMatrices$ZRural
    NregionRural = aggregatedResults$regionMatrices$NRural
    PregionRural = aggregatedResults$regionMatrices$pRural
    
    # national level
    Znation = aggregatedResults$nationalMatrices$Z
    Nnation = aggregatedResults$nationalMatrices$N
    Pnation = aggregatedResults$nationalMatrices$p
    ZnationUrban = aggregatedResults$nationalMatrices$ZUrban
    NnationUrban = aggregatedResults$nationalMatrices$NUrban
    PnationUrban = aggregatedResults$nationalMatrices$pUrban
    ZnationRural = aggregatedResults$nationalMatrices$ZRural
    NnationRural = aggregatedResults$nationalMatrices$NRural
    PnationRural = aggregatedResults$nationalMatrices$pRural
    
    # plots for testing purposes
    if(FALSE) {
      ## lcpb
      # aggregatedResultslcpb = aggregatePixelPredictions(lcpb, Ng, popGrid=popMat, useDensity=TRUE, countyLevel=countyLevel, 
      #                                               regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
      lcpbc = matrix(logitNormMean(cbind(c(logit(as.matrix(uc))), rep(sigmaEpsilonDraws, each=nrow(uc)))), nrow=nrow(uc))
      # this takes quite a while (15 mins?)
      aggregatedResultslcpb = aggregateEAPredictions(lcpbc, Ncs, areaMat, urbanMat, easpa=easpa, countyLevel=countyLevel, 
                                                     regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
      
      Ng2 = Ng
      Ng2[Ng2 != 0] = 1
      
      # NOTE: the following MUST be INCORRECT. set useDensity to TRUE, and don't use eaSamples at all!
      aggregatedResultslcpb2 = aggregatePixelPredictions(lcpb*eaSamples, eaSamples, popMatAdjusted=adjustedPopMat, useDensity=FALSE, countyLevel=countyLevel, 
                                                         regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=FALSE)
      
      # county level
      # Pcountylcpb = aggregatedResultslcpb$countyResults$Z
      # PcountyUrbanlcpb = aggregatedResultslcpb$countyResults$ZUrban
      # PcountyRurallcpb = aggregatedResultslcpb$countyResults$ZRural
      
      Pcountylcpb = aggregatedResultslcpb2$countyMatrices$p
      PcountyUrbanlcpb2 = aggregatedResultslcpb$countyResults$pUrban
      PcountyRurallcpb2 = aggregatedResultslcpb$countyResults$pRural
      
      # province level
      Pregionlcpb = aggregatedResultslcpb$regionResults$Z
      PregionUrbanlcpb = aggregatedResultslcpb$regionResults$ZUrban
      PregionRurallcpb = aggregatedResultslcpb$regionResults$ZRural
      
      # national level
      Pnationlcpb = aggregatedResultslcpb$nationalResults$Z
      PnationUrbanlcpb = aggregatedResultslcpb$nationalResults$ZUrban
      PnationRurallcpb = aggregatedResultslcpb$nationalResults$ZRural
      
      ## LCpb
      # this takes quite a while (15 mins?)
      aggregatedResultsLCpb = aggregateEAPredictions(lcpbc, Ncs, areaMat, urbanMat, easpa=easpa, countyLevel=countyLevel, 
                                                     regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
      mug2 = mug
      mug2[!is.finite(mug2)] = 0
      aggregatedResultsLCpb = aggregatePixelPredictions(mug*eaSamples, eaSamples, popMatAdjusted=adjustedPopMat, useDensity=FALSE, countyLevel=countyLevel, 
                                                        regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=FALSE)
      
      # county level
      PcountyLCpb = aggregatedResultsLCpb$countyMatrices$p
      PcountyUrbanLCpb = aggregatedResultsLCpb$countyResults$pUrban
      PcountyRuralLCpb = aggregatedResultsLCpb$countyResults$pRural
      
      # province level
      PregionLCpb = aggregatedResultsLCpb$regionResults$p
      PregionUrbanLCpb = aggregatedResultsLCpb$regionResults$pUrban
      PregionRuralLCpb = aggregatedResultsLCpb$regionResults$pRural
      
      # national level
      PnationLCpb = aggregatedResultsLCpb$nationalResults$Z
      PnationUrbanLCpb = aggregatedResultsLCpb$nationalResults$ZUrban
      PnationRuralLCpb = aggregatedResultsLCpb$nationalResults$ZRural
      
      save(Pcountylcpb, PcountyLCpb, Pcounty, file=paste0(outputDirectory, "test/clusterTest.RData"))
      out = load(paste0(outputDirectory, "test/clusterTest.RData"))
      
      # load shape files for plotting
      require(maptools)
      # regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
      out = load("../LK-INLA/regionMap.RData")
      out = load("../U5MR/adminMapData.RData")
      kenyaMap = adm0
      countyMap = adm1
      
      meanCols=makeRedBlueDivergingColors(64, rev = TRUE)
      sdCols=makeBlueYellowSequentialColors(64)
      popCols=makePurpleYellowSequentialColors(64, rev=TRUE)
      
      ## plot county results
      # plot mean
      # meanRange = range(c(rowMeans(Pcountylcpb, na.rm=TRUE), rowMeans(PcountyLCpb, na.rm=TRUE), rowMeans(Pcounty, na.rm=TRUE)), na.rm=TRUE)
      # meanTicks = pretty(meanRange, n=5)
      # meanTickLabels = as.character(meanTicks)
      meanRange = quantile(probs=c(.005, .995), c(rowMeans(mug, na.rm=TRUE), rowMeans(pg, na.rm=TRUE), rowMeans(lcpb, na.rm=TRUE)), na.rm=TRUE)
      meanTicks = pretty(meanRange, n=5)[-1]
      meanTickLabels = as.character(meanTicks)
      
      png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugMeanCounty.png"), width=1500, height=600)
      par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
      
      plotMapDat(plotVar=rowMeans(Pcountylcpb, na.rm=TRUE), new = TRUE, 
                 main="Posterior mean (lcpb model)", scaleFun=logit, scaleFunInverse=expit, 
                 cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
                 legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                 plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
      plotMapDat(mapDat=regionMap, lwd=2.5)
      
      plotMapDat(plotVar=rowMeans(PcountyLCpb, na.rm=TRUE), new = TRUE, 
                 main="Posterior mean (LCpb model)", scaleFun=logit, scaleFunInverse=expit, 
                 cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
                 legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                 plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
      plotMapDat(mapDat=regionMap, lwd=2.5)
      
      plotMapDat(plotVar=rowMeans(Pcounty, na.rm=TRUE), new = TRUE, 
                 main="Posterior mean (LCPB model)", scaleFun=logit, scaleFunInverse=expit, 
                 cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
                 legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                 plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
      plotMapDat(mapDat=regionMap, lwd=2.5)
      
      # image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
      #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
      #            legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
      dev.off()
      
      # plot SDs
      SDslcpb = apply(Pcountylcpb, 1, sd, na.rm=TRUE)
      SDsLCpb = apply(PcountyLCpb, 1, sd, na.rm=TRUE)
      SDsPcounty = apply(Pcounty, 1, sd, na.rm=TRUE)
      allSDs = c(SDslcpb, SDsLCpb, SDsPcounty)
      sdRange = range(allSDs, na.rm=TRUE)
      sdTicks = pretty(sdRange, n=5)
      sdTickLabels = as.character(sdTicks)
      
      png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugSDCounty.png"), width=1500, height=600)
      par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
      plotMapDat(plotVar=SDslcpb, new = TRUE, 
                 main="Posterior SD (lcpb model)", scaleFun=log, scaleFunInverse=exp, 
                 cols=sdCols, zlim=log(sdRange), ticks=sdTicks, tickLabels=sdTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
                 legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                 plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
      plotMapDat(mapDat=regionMap, lwd=2.5)
      
      plotMapDat(plotVar=SDsLCpb, new = TRUE, 
                 main="Posterior SD (LCpb model)", scaleFun=log, scaleFunInverse=exp, 
                 cols=sdCols, zlim=log(sdRange), ticks=sdTicks, tickLabels=sdTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
                 legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                 plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
      plotMapDat(mapDat=regionMap, lwd=2.5)
      
      plotMapDat(plotVar=SDsPcounty, new = TRUE, 
                 main="Posterior SD (LCPB model)", scaleFun=log, scaleFunInverse=exp, 
                 cols=sdCols, zlim=log(sdRange), ticks=sdTicks, tickLabels=sdTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = TRUE, 
                 legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                 plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
      plotMapDat(mapDat=regionMap, lwd=2.5)
      
      dev.off()
    }
    
    # aggregate timings to the time required for each model:
    # startTime0
    # finishedSetupTime1
    # finishedDrawingEAsTime2
    # finishedDrawingClusterEffectsTime3
    # finishedDrawingNsTime4
    # finishedCalculatingMusTime5
    # finishedDrawingZsTime6
    # finishedPixelAggregationNsTime7
    # finishedPixelAggregationZsTime8
    # finishedPixelAggregationPsTime9
    # finishedClusterIntegrationTime10
    # finishedPixelAggregationLCpbTime11
    # finishedPixelAggregationLCPbTime12
    # finishedArealAggregationSetupTime13
    # finishedArealAggregationLCPBTime14
    # finishedArealAggregationlcpbTime15
    # finishedArealAggregationLcpbTime16
    # finishedArealAggregationLCpbTime17
    # finishedArealAggregationLCPbTime18
    # finishedArealAggregationRiskSubstitutionTime19
    setupTime = finishedSetupTime1 - startTime0
    drawEAsTime = finishedDrawingEAsTime2 - finishedSetupTime1
    drawClustersTime = finishedDrawingClusterEffectsTime3 - finishedDrawingEAsTime2
    drawNsTime = finishedDrawingNsTime4 - finishedDrawingClusterEffectsTime3
    getMusTime = finishedCalculatingMusTime5 - finishedDrawingNsTime4
    drawZsTime = finishedDrawingZsTime6 - finishedCalculatingMusTime5
    pixelAggregationNsTime = finishedPixelAggregationNsTime7 - finishedDrawingZsTime6
    pixelAggregationZsTime = finishedPixelAggregationZsTime8 - finishedPixelAggregationNsTime7
    pixelAggregationPsTime = finishedPixelAggregationPsTime9 - finishedPixelAggregationZsTime8
    pixelAggregationTime = finishedPixelAggregationPsTime9 - finishedDrawingZsTime6
    clusterIntegrationTime = finishedClusterIntegrationTime10 - finishedPixelAggregationPsTime9
    pixelAggregationLCpbTime = finishedPixelAggregationLCpbTime11 - finishedClusterIntegrationTime10
    pixelAggregationLCPbTime = finishedPixelAggregationLCPbTime12 - finishedPixelAggregationLCpbTime11
    arealAggregationSetupTime = finishedArealAggregationSetupTime13 - finishedPixelAggregationLCPbTime12
    LCPBArealAggregationTime = finishedArealAggregationLCPBTime14 - finishedArealAggregationSetupTime13
    lcpbArealAggregationTime = finishedArealAggregationlcpbTime15 - finishedArealAggregationLCPBTime14
    LcpbArealAggregationTime = finishedArealAggregationLcpbTime16 - finishedArealAggregationlcpbTime15
    LCpbArealAggregationTime = finishedArealAggregationLCpbTime17 - finishedArealAggregationLcpbTime16
    LCPbArealAggregationTime = finishedArealAggregationLCPbTime18 - finishedArealAggregationLCpbTime17
    riskSubstitutionTime = finishedArealAggregationRiskSubstitutionTime19 - finishedArealAggregationLCPbTime18
    
    lcpbTime = setupTime + clusterIntegrationTime + arealAggregationSetupTime + lcpbArealAggregationTime
    LcpbTime = setupTime + drawEAsTime + clusterIntegrationTime + arealAggregationSetupTime + LcpbArealAggregationTime + riskSubstitutionTime
    LCpbTime = setupTime + drawEAsTime + drawClustersTime + getMusTime + pixelAggregationLCpbTime + arealAggregationSetupTime + LCpbArealAggregationTime + riskSubstitutionTime
    LCPbTime = setupTime + drawEAsTime + drawClustersTime + getMusTime + pixelAggregationNsTime + drawNsTime + pixelAggregationLCPbTime + arealAggregationSetupTime + LCPbArealAggregationTime + riskSubstitutionTime
    LCPBTime = setupTime + drawEAsTime + drawClustersTime + getMusTime + pixelAggregationTime + drawNsTime + drawZsTime + arealAggregationSetupTime + LCPBArealAggregationTime + riskSubstitutionTime
    
    totalTime = finishedArealAggregationRiskSubstitutionTime19 - startTime0
    timings = matrix(c(lcpbTime=lcpbTime, LcpbTime=LcpbTime, LCpbTime=LCpbTime, LCPbTime=LCPbTime, LCPBTime=LCPBTime, 
                       setupTime=setupTime, drawEAsTime=drawEAsTime, clusterIntegrationTime=clusterIntegrationTime, 
                       drawClustersTime=drawClustersTime, drawNsTime=drawNsTime, getMusTime=getMusTime, drawZsTime=drawZsTime, 
                       pixelAggregationNsTime=pixelAggregationNsTime, pixelAggregationZsTime=pixelAggregationZsTime, pixelAggregationPsTime=pixelAggregationPsTime, 
                       pixelAggregationTime=pixelAggregationTime, pixelAggregationLCpbTime=pixelAggregationLCpbTime, 
                       pixelAggregationLCPbTime=pixelAggregationLCPbTime, arealAggregationSetupTime=arealAggregationSetupTime, 
                       lcpbArealAggregationTime=lcpbArealAggregationTime, 
                       LcpbArealAggregationTime=LcpbArealAggregationTime, 
                       LCpbArealAggregationTime=LCpbArealAggregationTime, 
                       LCPbArealAggregationTime=LCPbArealAggregationTime, 
                       LCPBArealAggregationTime=LCPBArealAggregationTime, 
                       riskSubstitutionTime=riskSubstitutionTime, 
                       totalTime = totalTime), nrow=1)
    timings = rbind(timings, timings / totalTime)
    row.names(timings) = c("Time (s)", "Prop Total")
    colnames(timings) = c("lcpbTime", "LcpbTime", "LCpbTime", "LCPbTime", "LCPBTime", 
                          "setupTime", "drawEAsTime", "clusterIntegrationTime", 
                          "drawClustersTime", "drawNsTime", "getMusTime", "drawZsTime", 
                          "pixelAggregationNsTime", "pixelAggregationZsTime", "pixelAggregationPsTime", 
                          "pixelAggregationTime", "pixelAggregationLCpbTime", "pixelAggregationLCPbTime", 
                          "arealAggregationSetupTime", 
                          "lcpbArealAggregationTime", 
                          "LcpbArealAggregationTime", 
                          "LCpbArealAggregationTime", 
                          "LCPbArealAggregationTime", 
                          "LCPBArealAggregationTime", 
                          "riskSubstitutionTime", 
                          "totalTime")
    
    timingslcpb = matrix(c(lcpbTime=lcpbTime, 
                           setupTime=setupTime, clusterIntegrationTime=clusterIntegrationTime, 
                           arealAggregationSetupTime=arealAggregationSetupTime, 
                           lcpbArealAggregationTime=lcpbArealAggregationTime), nrow=1)
    timingslcpb = rbind(timingslcpb, timingslcpb / lcpbTime)
    row.names(timingslcpb) = c("Time (s)", "Prop Total")
    colnames(timingslcpb) = c("lcpbTime", 
                               "setupTime", "clusterIntegrationTime", 
                               "arealAggregationSetupTime", 
                               "lcpbArealAggregationTime")
    
    timingsLcpb = matrix(c(LcpbTime=LcpbTime, 
                           setupTime=setupTime, drawEAsTime=drawEAsTime, clusterIntegrationTime=clusterIntegrationTime, 
                           arealAggregationSetupTime=arealAggregationSetupTime, LcpbArealAggregationTime=LcpbArealAggregationTime, 
                           riskSubstitutionTime=riskSubstitutionTime), nrow=1)
    timingsLcpb = rbind(timingsLcpb, timingsLcpb / LcpbTime)
    row.names(timingsLcpb) = c("Time (s)", "Prop Total")
    colnames(timingsLcpb) = c("LcpbTime", 
                               "setupTime", "drawEAsTime", "clusterIntegrationTime", 
                               "arealAggregationSetupTime", "LcpbArealAggregationTime", 
                               "riskSubstitutionTime")
    
    timingsLCpb = matrix(c(LCpbTime=LCpbTime, 
                           setupTime=setupTime, drawEAsTime=drawEAsTime, drawClustersTime=drawClustersTime, getMusTime=getMusTime, 
                           pixelAggregationLCpbTime=pixelAggregationLCpbTime, arealAggregationSetupTime=arealAggregationSetupTime, 
                           LCpbArealAggregationTime=LCpbArealAggregationTime, 
                           riskSubstitutionTime=riskSubstitutionTime), nrow=1)
    timingsLCpb = rbind(timingsLCpb, timingsLCpb / LCpbTime)
    row.names(timingsLCpb) = c("Time (s)", "Prop Total")
    colnames(timingsLCpb) = c("LCpbTime", 
                               "setupTime", "drawEAsTime", "drawClustersTime", "getMusTime", 
                               "arealAggregationSetupTime", "pixelAggregationLCpbTime", "LCpbArealAggregationTime", 
                               "riskSubstitutionTime")
    
    timingsLCPb = matrix(c(LCPbTime=LCPbTime, 
                           setupTime=setupTime, drawEAsTime=drawEAsTime, drawClustersTime=drawClustersTime, getMusTime=getMusTime, 
                           pixelAggregationNsTime=pixelAggregationNsTime, drawNsTime=drawNsTime, 
                           pixelAggregationLCPbTime=pixelAggregationLCPbTime, arealAggregationSetupTime=arealAggregationSetupTime, LCPbArealAggregationTime=LCPbArealAggregationTime, 
                           riskSubstitutionTime=riskSubstitutionTime), nrow=1)
    timingsLCPb = rbind(timingsLCPb, timingsLCPb / LCPbTime)
    row.names(timingsLCPb) = c("Time (s)", "Prop Total")
    colnames(timingsLCPb) = c("LCPbTime", 
                               "setupTime", "drawEAsTime", "drawClustersTime", "getMusTime", 
                               "pixelAggregationNsTime", "drawNsTime", "pixelAggregationLCPbTime", 
                               "arealAggregationSetupTime", "LCPbArealAggregationTime", 
                               "riskSubstitutionTime")
    
    timingsLCPB = matrix(c(LCPBTime=LCPBTime, 
                           setupTime=setupTime, drawEAsTime=drawEAsTime, drawClustersTime=drawClustersTime, getMusTime=getMusTime, 
                           pixelAggregationTime=pixelAggregationTime, drawNsTime=drawNsTime, drawZsTime=drawZsTime, 
                           arealAggregationSetupTime=arealAggregationSetupTime, LCPBArealAggregationTime=LCPBArealAggregationTime, 
                           riskSubstitutionTime=riskSubstitutionTime), nrow=1)
    timingsLCPB = rbind(timingsLCPB, timingsLCPB / LCPBTime)
    row.names(timingsLCPB) = c("Time (s)", "Prop Total")
    colnames(timingsLCPB) = c("LCPBTime", 
                           "setupTime", "drawEAsTime", "drawClustersTime", "getMusTime", 
                           "pixelAggregationTime", "drawNsTime", "drawZsTime", 
                           "arealAggregationSetupTime", "LCPBArealAggregationTime", 
                           "riskSubstitutionTime")
    
    allTimings = list(timings=timings, 
                      timingslcpb=timingslcpb, 
                      timingsLcpb=timingsLcpb, 
                      timingsLCpb=timingsLCpb, 
                      timingsLCPb=timingsLCPb, 
                      timingsLCPB=timingsLCPB)
    
    ##### Extra steps: collect draws at each level and generate:
    ##### areas, preds, 
    if(pixelLevel) {
      allMatrices = list(pixelMatricesLCPB=list(p=pg, Z=Zg, N=Ng), aggregatedResultsLCPB=aggregatedResults, 
                         pixelMatriceslcpb=list(p=lcpb), aggregatedResultslcpb=aggregatedResultslcpb, 
                         pixelMatricesLcpb=list(p=Lcpb), aggregatedResultsLcpb=aggregatedResultsLcpb, 
                         pixelMatricesLCpb=list(p=LCpb), aggregatedResultsLCpb=aggregatedResultsLCpb, 
                         pixelMatricesLCPb=list(p=LCPb), aggregatedResultsLCPb=aggregatedResultsLCPb, 
                         allTimings=allTimings)
    } else {
      allMatrices = list(aggregatedResultsLCPB=aggregatedResults, 
                         aggregatedResultslcpb=aggregatedResultslcpb, 
                         aggregatedResultsLcpb=aggregatedResultsLcpb, 
                         aggregatedResultsLCpb=aggregatedResultsLCpb, 
                         aggregatedResultsLCPb=aggregatedResultsLCPb, 
                         allTimings=allTimings)
    }
  } else {
    allTimings=NULL
    
    # for the if we only care about the modified pixel level predictions, just include those
    if(is.null(clustersPerPixel)) {
      allMatrices = list(pixelMatricesLCPB=list(p=pgMod, Z=ZgMod, N=NgMod), 
                         pixelMatriceslcpb=list(p=lcpb), allTimings=allTimings)
    } else {
      allMatrices = list(pixelMatricesLCPB=list(p=pg, Z=Zg, N=Ng), 
                         pixelMatriceslcpb=list(p=lcpb), allTimings=allTimings)
    }
  }
  
  if(!returnEAinfo) {
    allMatrices
  } else {
    theseI = pixelIndexMat[,1]
    eaDat = data.frame(lon=popGrid$lon[theseI], lat=popGrid$lat[theseI], popOverall=popGrid$popOrig[theseI], 
                       region=popGrid$region[theseI], admin1=popGrid$admin1[theseI], admin2=popGrid$admin2[theseI], 
                       urban=popGrid$urban[theseI], east=popGrid$east[theseI], north=popGrid$north[theseI], 
                       popTarget=popGridAdjusted$popOrig[theseI], pixelIs=theseI, 
                       nHH=householdDraws[,1], n=Ncs[,1], y=Zc[,1], 
                       plcpb=lcpb[theseI,1], pLcpb=lcpb[theseI,1], pLCpb=muc[,1], pLCPb=muc[,1], pLCPB=Zc[,1]/Ncs[,1])
    eaDat$pLCPB[eaDat$n == 0] = muc[eaDat$n == 0,1]
    c(allMatrices, list(eaDat=eaDat))
  }
}

# use the fitSPDEKenyaDat function to validate SPDE model to binomial data within Kenya with leave one 
# region out validation, and prediction at cluster level
validateAggregationModelsKenyaDat = function(dat=NULL, dataType=c("mort", "ed"), 
                                             mesh=getSPDEMeshKenya(), prior=getSPDEPrior(mesh), 
                                             significanceCI=.8, int.strategy="ccd", strategy="gaussian", 
                                             nPostSamples=1000, verbose=FALSE, link=1, seed=123, 
                                             urbanEffect=TRUE, clusterEffect=TRUE, kmres=5, 
                                             loadPreviousFit=TRUE, saveResults=TRUE, 
                                             sampleTable=NULL, stratifiedValidation=TRUE, 
                                             doPixelLevelValidation=TRUE, doCountyLevelValidation=TRUE, 
                                             loadPreviousResults=FALSE, 
                                             doLCPb=TRUE, doLCpb=TRUE, doLcpb=TRUE) {
  
  if(!is.null(seed))
    set.seed(seed)
  
  # load observations
  dataType = match.arg(dataType)
  if(is.null(dat)) {
    if(dataType == "mort") {
      dat = mort
      dataType2 = "children"
    }
    else {
      dat = ed
      dataType2 = "women"
    }
  }
  if(dataType == "mort")
    fileNameRoot = "Mort"
  else
    fileNameRoot = "Ed"
  
  # first fit the full model (we will use this to initialize the model during the validation fits for each left out county)
  
  fileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationFull", 
                    "_urb", urbanEffect, ".RData")
  if(!loadPreviousFit || !file.exists(fileName)) {
    
    tempFileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationFullPointwise", 
                          "_urb", urbanEffect, ".RData")
    if(!file.exists(tempFileName)) {
      print("Fitting full point level model")
      timePoint = system.time(fit <- fitSPDEKenyaDat(dat, dataType, mesh, prior, significanceCI, int.strategy, strategy, nPostSamples,
                                                     verbose, link, NULL, urbanEffect, clusterEffect=TRUE,
                                                     kmres=kmres, doValidation=TRUE, family="binomial"))[3]
      
      if(saveResults)
        save(fit, timePoint, file=tempFileName)
    } else {
      print("Loading full point level model")
      load(tempFileName)
    }
    
    if(urbanEffect) {
      urbanEffectValue = fit$mod$summary.fixed[2,1]
    } else {
      urbanEffectValue = NULL
    }
    
    print("Fitting full population aggregation model")
    popMat = makeDefaultPopMat()
    truthTableFull = getPixelLevelTruth(dat=dat, popMat=popMat, targetPop="children")
    clustersPerPixel = rep(0, nrow(popMat))
    clustersPerPixel[truthTableFull$pixelI] = truthTableFull$nClusters
    easpa = makeDefaultEASPA(dataType2)
    easpaMort = makeDefaultEASPA(dataType2, NULL, useClustersAsEAs=TRUE)
    easpaMortClusterUrban = makeDefaultEASPA(dataType2, NULL, useClustersAsEAs=TRUE, usePixelUrban=FALSE)
    
    # do the pixel level validation
    browser()
    timeAggregationPixel = system.time(fit2 <- modLCPB(uDraws=fit$uDraws, fit$sigmaEpsilonDraws, easpa=easpaMort, popMat=popMat, empiricalDistributions=NULL, 
                                                  includeUrban=urbanEffect, clusterLevel=FALSE, pixelLevel=TRUE, constituencyLevel=FALSE, countyLevel=FALSE, 
                                                  regionLevel=FALSE, nationalLevel=FALSE, doModifiedPixelLevel=FALSE, validationPixelI=truthTableFull$pixelI, 
                                                  onlyDoModifiedPixelLevel=FALSE, clustersPerPixel=clustersPerPixel, 
                                                  doLCPb=doLCPb, doLCpb=doLCpb, doLcpb=doLcpb))[3]
    
    # do the county and province level validation
    adjustedPopMat = makeDefaultPopMat(getAdjusted=TRUE)
    
    # do the same but making population level estimates rather than estimates of the data
    timeAggregationCountyProvincePopulation = system.time(fit4 <- modLCPB(uDraws=fit$uDraws, fit$sigmaEpsilonDraws, easpa=easpa, popMat=popMat, 
                                                                          adjustedPopMat=adjustedPopMat, empiricalDistributions=NULL, 
                                                                          includeUrban=urbanEffect, clusterLevel=FALSE, pixelLevel=FALSE, constituencyLevel=TRUE, countyLevel=TRUE, 
                                                                          regionLevel=TRUE, nationalLevel=FALSE, doModifiedPixelLevel=FALSE, 
                                                                          onlyDoModifiedPixelLevel=FALSE, 
                                                                          doLCPb=doLCPb, doLCpb=doLCpb, doLcpb=doLcpb, urbanEffect=urbanEffectValue))[3]
    browser()
    timeAggregationCountyProvince = system.time(fit3 <- modLCPB(uDraws=fit$uDraws, fit$sigmaEpsilonDraws, easpa=easpaMortClusterUrban, popMat=popMat, 
                                                                adjustedPopMat=adjustedPopMat, empiricalDistributions=NULL, 
                                                                includeUrban=urbanEffect, clusterLevel=FALSE, pixelLevel=FALSE, constituencyLevel=TRUE, countyLevel=TRUE, 
                                                                regionLevel=TRUE, nationalLevel=FALSE, doModifiedPixelLevel=FALSE, 
                                                                onlyDoModifiedPixelLevel=FALSE, 
                                                                doLCPb=doLCPb, doLCpb=doLCpb, doLcpb=doLcpb))[3]
    
    ## pixel level
    # get observations and prediction summary statistics
    # easpaFull = makeDefaultEASPA(validationClusterI=NULL, useClustersAsEAs=TRUE)
    truthFull = truthTableFull$p
    obsUrban = truthTableFull$urban
    
    cpo = fit$mod$cpo$cpo
    cpoFailure = fit$mod$cpo$failure
    dic = fit$mod$dic$dic
    waic = fit$mod$waic$waic
    modelFit = fit$mod
    
    est = rowMeans(fit2$pixelMatriceslcpb$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit2$pixelMatricesLCPB$p
    # vars = apply(fit2$pixelMatricesLCPB$p, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    # estMatBinomial = addBinomialVar(estMat, dat$n)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPB, pixel):")
    fullPooledScoresPixelLCPB = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel))
    print(fullPooledScoresPixelLCPB)
    print("Rural scores (LCPB, pixel):")
    fullRuralScoresPixelLCPB = data.frame(c(getScores(truthFull[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresPixelLCPB)
    print("Urban scores (LCPB, pixel):")
    fullUrbanScoresPixelLCPB = data.frame(c(getScores(truthFull[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresPixelLCPB)
    
    # est = rowMeans(fit2$pixelMatriceslcpb$p)
    # vars = apply(fit2$pixelMatricesLCPb$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit2$pixelMatricesLCPb$p
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    # estMatBinomial = addBinomialVar(estMat, dat$n)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPb, pixel):")
    fullPooledScoresPixelLCPb = data.frame(c(getScores(truthFull, est, vars, lower, upper, as.matrix(estMat), doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    # fullPooledScoresPixelLCPb2 = data.frame(c(getScores(truthFull[truthFull != 0], est[truthFull != 0], vars[truthFull != 0], lower[truthFull != 0], upper[truthFull != 0], as.matrix(estMat)[truthFull != 0,], doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresPixelLCPb)
    print("Rural scores (LCPb, pixel):")
    fullRuralScoresPixelLCPb = data.frame(c(getScores(truthFull[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresPixelLCPb)
    print("Urban scores (LCPb, pixel):")
    fullUrbanScoresPixelLCPb = data.frame(c(getScores(truthFull[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresPixelLCPb)
    
    # est = rowMeans(fit2$pixelMatriceslcpb$p)
    # vars = apply(fit2$pixelMatricesLCpb$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit2$pixelMatricesLCpb$p
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    # estMatBinomial = addBinomialVar(estMat, dat$n)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCpb, pixel):")
    fullPooledScoresPixelLCpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresPixelLCpb)
    print("Rural scores (LCpb, pixel):")
    fullRuralScoresPixelLCpb = data.frame(c(getScores(truthFull[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresPixelLCpb)
    print("Urban scores (LCpb, pixel):")
    fullUrbanScoresPixelLCpb = data.frame(c(getScores(truthFull[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresPixelLCpb)
    
    # est = rowMeans(fit2$pixelMatriceslcpb$p)
    # vars = apply(fit2$pixelMatricesLcpb$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit2$pixelMatricesLcpb$p
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    # estMatBinomial = addBinomialVar(estMat, dat$n)
    
    # # calculate validation scoring rules
    print("Pooled scores (Lcpb, pixel):")
    fullPooledScoresPixelLcpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresPixelLcpb)
    print("Rural scores (Lcpb, pixel):")
    fullRuralScoresPixelLcpb = data.frame(c(getScores(truthFull[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresPixelLcpb)
    print("Urban scores (Lcpb, pixel):")
    fullUrbanScoresPixelLcpb = data.frame(c(getScores(truthFull[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresPixelLcpb)
    
    # est = rowMeans(fit2$pixelMatriceslcpb$p)
    # vars = apply(fit2$pixelMatriceslcpb$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit2$pixelMatriceslcpb$p
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    # estMatBinomial = addBinomialVar(estMat, dat$n)
    
    # # calculate validation scoring rules
    print("Pooled scores (lcpb, pixel):")
    fullPooledScoresPixellcpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresPixellcpb)
    print("Rural scores (lcpb, pixel):")
    fullRuralScoresPixellcpb = data.frame(c(getScores(truthFull[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresPixellcpb)
    print("Urban scores (lcpb, pixel):")
    fullUrbanScoresPixellcpb = data.frame(c(getScores(truthFull[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresPixellcpb)
    
    fullScoresPixel = list(fullPooledScoreslcpb=fullPooledScoresPixellcpb, fullRuralScoreslcpb=fullRuralScoresPixellcpb, fullUrbanScoreslcpb=fullUrbanScoresPixellcpb, 
                           fullPooledScoresLcpb=fullPooledScoresPixelLcpb, fullRuralScoresLcpb=fullRuralScoresPixelLcpb, fullUrbanScoresLcpb=fullUrbanScoresPixelLcpb, 
                           fullPooledScoresLCpb=fullPooledScoresPixelLCpb, fullRuralScoresLCpb=fullRuralScoresPixelLCpb, fullUrbanScoresLCpb=fullUrbanScoresPixelLCpb, 
                           fullPooledScoresLCPb=fullPooledScoresPixelLCPb, fullRuralScoresLCPb=fullRuralScoresPixelLCPb, fullUrbanScoresLCPb=fullUrbanScoresPixelLCPb, 
                           fullPooledScoresLCPB=fullPooledScoresPixelLCPB, fullRuralScoresLCPB=fullRuralScoresPixelLCPB, fullUrbanScoresLCPB=fullUrbanScoresPixelLCPB)
    
    ## Constituency level
    browser()
    truthTableFullConstituency = getConstituencyLevelTruth(dat=dat, targetPop="children")
    
    # get observations and prediction summary statistics
    # easpaFull = makeDefaultEASPA(validationClusterI=NULL, useClustersAsEAs=TRUE)
    truthFull = truthTableFullConstituency$pStratified
    truthFullUrban = truthTableFullConstituency$pUrban
    truthFullRural = truthTableFullConstituency$pRural
    
    truthDirectFull = truthTableFullConstituency$pDirect
    truthDirectFullUrban = truthTableFullConstituency$pDirectUrban
    truthDirectFullRural = truthTableFullConstituency$pDirectRural
    truthDirectFullLogitVar = truthTableFullConstituency$varDirectLogit
    truthDirectFullUrbanLogitVar = truthTableFullConstituency$varDirectUrbanLogit
    truthDirectFullRuralLogitVar = truthTableFullConstituency$varDirectRuralLogit
    keepFullI = is.finite(truthDirectFullLogitVar) # & (poppcon$popTotal != 0)
    keepFullUrbanI = (is.finite(truthDirectFullUrbanLogitVar) & (truthDirectFullUrbanLogitVar > 1e-15)) # & (poppcon$popUrb != 0)
    keepFullRuralI = is.finite(truthDirectFullRuralLogitVar) & (truthDirectFullRuralLogitVar > 1e-15) # & (poppcon$popRur != 0)
    # truthDirectFull = truthTableFullConstituency$pDirect[keepFullI]
    # truthDirectFullUrban = truthTableFullConstituency$pDirectUrban[keepFullUrbanI]
    # truthDirectFullRural = truthTableFullConstituency$pDirectRural[keepFullRuralI]
    # truthDirectFullLogitVar = truthTableFullConstituency$varDirectLogit[keepFullI]
    # truthDirectFullUrbanLogitVar = truthTableFullConstituency$varDirectUrbanLogit[keepFullUrbanI]
    # truthDirectFullRuralLogitVar = truthTableFullConstituency$varDirectRuralLogit[keepFullRuralI]
    
    if(FALSE) {
      makeRankPlots(fit4$aggregatedResultslcpb$constituencyMatrices$p, admin="admin2", 
                    plotNameSuffix = "Risk", savePlots = TRUE)
      makeRankPlots(fit4$aggregatedResultsLCPB$constituencyMatrices$p, admin="admin2", 
                    plotNameSuffix = "CPAM", savePlots = TRUE)
      
      widths80lcpb = apply(fit4$aggregatedResultslcpb$constituencyMatrices$p, 1, function(x) {diff(quantile(x, prob=c(.1, .9)))})
      widths80LCpb = apply(fit4$aggregatedResultsLCpb$constituencyMatrices$p, 1, function(x) {diff(quantile(x, prob=c(.1, .9)))})
      widths80LCPB = apply(fit4$aggregatedResultsLCPB$constituencyMatrices$p, 1, function(x) {diff(quantile(x, prob=c(.1, .9)))})
      mean((widths80LCPB - widths80lcpb) / widths80lcpb) # 0.03904018 (these comparisons are all using the logistic approximation for the lcpb model)
      mean((widths80LCpb - widths80lcpb) / widths80lcpb) # -0.003937166
      
      widths80lcpbZ = apply(fit4$aggregatedResultslcpb$constituencyMatrices$Z, 1, function(x) {diff(quantile(x, prob=c(.1, .9)))})
      widths80LCpbZ = apply(fit4$aggregatedResultsLCpb$constituencyMatrices$Z, 1, function(x) {diff(quantile(x, prob=c(.1, .9)))})
      widths80LCPBZ = apply(fit4$aggregatedResultsLCPB$constituencyMatrices$Z, 1, function(x) {diff(quantile(x, prob=c(.1, .9)))})
      mean((widths80LCPBZ - widths80lcpbZ) / widths80lcpbZ) # 0.1084088
      mean((widths80LCpbZ - widths80lcpbZ) / widths80lcpbZ) # 0.06681593
      
      lcpbRR = fit4$aggregatedResultslcpb$constituencyMatrices$pUrban / fit4$aggregatedResultslcpb$constituencyMatrices$pRural
      LCPBRR = fit4$aggregatedResultsLCPB$constituencyMatrices$pUrban / fit4$aggregatedResultsLCPB$constituencyMatrices$pRural
      range(lcpbRR) # 0.4498809 2.4913617
      range(LCPBRR[is.finite(LCPBRR)]) # 0.000000 9.040873
      hist(lcpbRR)
      hist(LCPBRR)
      out = lcpbRR - LCPBRR
      mean(out[is.finite(out)], na.rm=TRUE) # -0.01962926
      
      widths80lcpbRR = apply(fit4$aggregatedResultslcpb$constituencyMatrices$pUrban / fit4$aggregatedResultslcpb$constituencyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      widths80LCpbRR = apply(fit4$aggregatedResultsLCpb$constituencyMatrices$pUrban / fit4$aggregatedResultsLCpb$constituencyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      widths80LCPbRR = apply(fit4$aggregatedResultsLCPb$constituencyMatrices$pUrban / fit4$aggregatedResultsLCPb$constituencyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      widths80LCPBRR = apply(fit4$aggregatedResultsLCPB$constituencyMatrices$pUrban / fit4$aggregatedResultsLCPB$constituencyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      out = (widths80LCPBRR - widths80lcpbRR) / widths80lcpbRR
      mean(out[is.finite(out)], na.rm=TRUE) # 4.908996
      median(out[is.finite(out)], na.rm=TRUE) # 1.012442
      
      out = (widths80LCPbRR - widths80lcpbRR) / widths80lcpbRR
      mean(out[is.finite(out)], na.rm=TRUE) # 2.219203
      median(out[is.finite(out)], na.rm=TRUE) # 0.3092367
      
      out = (widths80LCpbRR - widths80lcpbRR) / widths80lcpbRR
      mean(out[is.finite(out)], na.rm=TRUE) # 2.194338]
      median(out[is.finite(out)], na.rm=TRUE) # 0.3021015
      
      widths80lcpbRR = apply(fit4$aggregatedResultslcpb$countyMatrices$pUrban / fit4$aggregatedResultslcpb$countyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      widths80LCpbRR = apply(fit4$aggregatedResultsLCpb$countyMatrices$pUrban / fit4$aggregatedResultsLCpb$countyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      widths80LCPbRR = apply(fit4$aggregatedResultsLCPb$countyMatrices$pUrban / fit4$aggregatedResultsLCPb$countyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      widths80LCPBRR = apply(fit4$aggregatedResultsLCPB$countyMatrices$pUrban / fit4$aggregatedResultsLCPB$countyMatrices$pRural, 1, function(x) {diff(quantile(x, prob=c(.1, .9), na.rm=TRUE))})
      out = (widths80LCPBRR - widths80lcpbRR) / widths80lcpbRR
      mean(out[is.finite(out)], na.rm=TRUE) # 0.2112706
      median(out[is.finite(out)], na.rm=TRUE) # 0.1803033
      
      out = (widths80LCPbRR - widths80lcpbRR) / widths80lcpbRR
      mean(out[is.finite(out)], na.rm=TRUE) # 0.08914678
      median(out[is.finite(out)], na.rm=TRUE) # 0.07005627
      
      out = (widths80LCpbRR - widths80lcpbRR) / widths80lcpbRR
      mean(out[is.finite(out)], na.rm=TRUE) # 0.09089862
      median(out[is.finite(out)], na.rm=TRUE) # 0.06804808
    }
    
    # obsUrban = truthTableFullCounty$urban
    
    cpo = fit$mod$cpo$cpo
    cpoFailure = fit$mod$cpo$failure
    dic = fit$mod$dic$dic
    waic = fit$mod$waic$waic
    modelFit = fit$mod
    
    est = rowMeans(fit3$aggregatedResultslcpb$constituencyMatrices$pUrban) * truthTableFullConstituency$urbanProp + rowMeans(fit3$aggregatedResultslcpb$constituencyMatrices$pRural) * (1 - truthTableFullConstituency$urbanProp)
    estUrban = rowMeans(fit3$aggregatedResultslcpb$constituencyMatrices$pUrban)
    estRural = rowMeans(fit3$aggregatedResultslcpb$constituencyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCPB$constituencyMatrices$pUrban, 1, truthTableFullConstituency$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCPB$constituencyMatrices$pRural, 1, 1-truthTableFullConstituency$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCPB$constituencyMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCPB$constituencyMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPB, constituency, stratified):")
    fullPooledScoresConstituencyLCPB = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationConstituencyProvince=timeAggregationConstituencyProvince))
    print(fullPooledScoresConstituencyLCPB)
    print("Rural scores (LCPB, constituency, stratified):")
    fullRuralScoresConstituencyLCPB = data.frame(c(getScores(truthFullRural[truthTableFullConstituency$urbanProp!=1], estRural[truthTableFullConstituency$urbanProp!=1], varsRural[truthTableFullConstituency$urbanProp!=1], lower, upper, estMatRural[truthTableFullConstituency$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationConstituencyProvince=timeAggregationConstituencyProvince))
    print(fullRuralScoresConstituencyLCPB)
    print("Urban scores (LCPB, constituency, stratified):")
    fullUrbanScoresConstituencyLCPB = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationConstituencyProvince=timeAggregationConstituencyProvince))
    print(fullUrbanScoresConstituencyLCPB)
    
    # do the same for the direct estimate validation
    est = rowMeans(fit4$aggregatedResultsLCpb$constituencyMatrices$p) #[keepFullI]
    estUrban = rowMeans(fit4$aggregatedResultsLCpb$constituencyMatrices$pUrban) #[keepFullUrbanI]
    estRural = rowMeans(fit4$aggregatedResultsLCpb$constituencyMatrices$pRural) #[keepFullRuralI]
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit4$aggregatedResultsLCPB$constituencyMatrices$p #[keepFullI,]
    # vars = apply(estMat, 1, var)
    vars = rep(NA, length(est))
    varsUrban = rep(NA, length(estUrban))
    varsRural = rep(NA, length(estRural))
    vars[keepFullI] = rowMeans(expit(sweep(logit(estMat[keepFullI,]) + matrix(rnorm(n=length(estMat[keepFullI,]), sd=sqrt(truthDirectFullLogitVar[keepFullI])), ncol=ncol(estMat)), 1, logit(est[keepFullI]), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultsLCPB$constituencyMatrices$pUrban #[keepFullUrbanI,]
    varsUrban[keepFullUrbanI] = rowMeans(expit(sweep(logit(estMatUrban[keepFullUrbanI,]) + matrix(rnorm(n=length(estMatUrban[keepFullUrbanI,]), sd=sqrt(truthDirectFullUrbanLogitVar[keepFullUrbanI])), ncol=ncol(estMatUrban)), 1, logit(estUrban[keepFullUrbanI]), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultsLCPB$constituencyMatrices$pRural #[keepFullRuralI,]
    varsRural[keepFullRuralI] = rowMeans(expit(sweep(logit(estMatRural[keepFullRuralI,]) + matrix(rnorm(n=length(estMatRural[keepFullRuralI,]), sd=sqrt(truthDirectFullRuralLogitVar[keepFullRuralI])), ncol=ncol(estMatRural)), 1, logit(estRural[keepFullRuralI]), FUN="-"))^2)
    # vars[is.nan(vars)] = 0
    # varsUrban[is.nan(varsUrban)] = 0
    # varsRural[is.nan(varsRural)] = 0
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPB, constituency, direct):")
    fullPooledScoresConstituencyLCPBDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvincePopulation=timeAggregationCountyProvincePopulation))
    print(fullPooledScoresConstituencyLCPBDirect)
    print("Rural scores (LCPB, constituency, direct):")
    fullRuralScoresConstituencyLCPBDirect = data.frame(c(getScores(truthDirectFullRural, estRural, varsRural, lower, upper, estMatRural, doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvincePopulation=timeAggregationCountyProvincePopulation))
    print(fullRuralScoresConstituencyLCPBDirect)
    print("Urban scores (LCPB, constituency, direct):")
    fullUrbanScoresConstituencyLCPBDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvincePopulation=timeAggregationCountyProvincePopulation))
    print(fullUrbanScoresConstituencyLCPBDirect)
    
    ## lcpb
    est = rowMeans(fit4$aggregatedResultslcpb$constituencyMatrices$p) #[keepFullI]
    estUrban = rowMeans(fit4$aggregatedResultslcpb$constituencyMatrices$pUrban) #[keepFullUrbanI]
    estRural = rowMeans(fit4$aggregatedResultslcpb$constituencyMatrices$pRural) #[keepFullRuralI]
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit4$aggregatedResultslcpb$constituencyMatrices$p #[keepFullI,]
    # vars = apply(estMat, 1, var)
    vars = rep(NA, length(est))
    varsUrban = rep(NA, length(estUrban))
    varsRural = rep(NA, length(estRural))
    vars[keepFullI] = rowMeans(expit(sweep(logit(estMat[keepFullI,]) + matrix(rnorm(n=length(estMat[keepFullI,]), sd=sqrt(truthDirectFullLogitVar[keepFullI])), ncol=ncol(estMat)), 1, logit(est[keepFullI]), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultslcpb$constituencyMatrices$pUrban #[keepFullUrbanI,]
    varsUrban[keepFullUrbanI] = rowMeans(expit(sweep(logit(estMatUrban[keepFullUrbanI,]) + matrix(rnorm(n=length(estMatUrban[keepFullUrbanI,]), sd=sqrt(truthDirectFullUrbanLogitVar[keepFullUrbanI])), ncol=ncol(estMatUrban)), 1, logit(estUrban[keepFullUrbanI]), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultslcpb$constituencyMatrices$pRural #[keepFullRuralI,]
    varsRural[keepFullRuralI] = rowMeans(expit(sweep(logit(estMatRural[keepFullRuralI,]) + matrix(rnorm(n=length(estMatRural[keepFullRuralI,]), sd=sqrt(truthDirectFullRuralLogitVar[keepFullRuralI])), ncol=ncol(estMatRural)), 1, logit(estRural[keepFullRuralI]), FUN="-"))^2)
    # vars[is.nan(vars)] = 0
    # varsUrban[is.nan(varsUrban)] = 0
    # varsRural[is.nan(varsRural)] = 0
    
    # # calculate validation scoring rules
    print("Pooled scores (lcpb, constituency, direct):")
    fullPooledScoresConstituencylcpbDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvincePopulation=timeAggregationCountyProvincePopulation))
    print(fullPooledScoresConstituencylcpbDirect)
    print("Rural scores (lcpb, constituency, direct):")
    fullRuralScoresConstituencylcpbDirect = data.frame(c(getScores(truthDirectFullRural, estRural, varsRural, lower, upper, estMatRural, doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvincePopulation=timeAggregationCountyProvincePopulation))
    print(fullRuralScoresConstituencylcpbDirect)
    print("Urban scores (lcpb, constituency, direct):")
    fullUrbanScoresConstituencylcpbDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvincePopulation=timeAggregationCountyProvincePopulation))
    print(fullUrbanScoresConstituencylcpbDirect)
    
    ##### County level
    browser()
    truthTableFullCounty = getCountyLevelTruth(dat=dat, easpa=easpa, targetPop="children")
    
    # get observations and prediction summary statistics
    # easpaFull = makeDefaultEASPA(validationClusterI=NULL, useClustersAsEAs=TRUE)
    truthFull = truthTableFullCounty$pStratified
    truthFullUrban = truthTableFullCounty$pUrban
    truthFullRural = truthTableFullCounty$pRural
    
    truthDirectFull = truthTableFullCounty$pDirect
    truthDirectFullUrban = truthTableFullCounty$pDirectUrban
    truthDirectFullRural = truthTableFullCounty$pDirectRural
    truthDirectFullLogitVar = truthTableFullCounty$varDirectLogit
    truthDirectFullUrbanLogitVar = truthTableFullCounty$varDirectUrbanLogit
    truthDirectFullRuralLogitVar = truthTableFullCounty$varDirectRuralLogit
    # obsUrban = truthTableFullCounty$urban
    
    cpo = fit$mod$cpo$cpo
    cpoFailure = fit$mod$cpo$failure
    dic = fit$mod$dic$dic
    waic = fit$mod$waic$waic
    modelFit = fit$mod
    
    browser()
    est = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pUrban) * truthTableFullCounty$urbanProp + rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pRural) * (1 - truthTableFullCounty$urbanProp)
    estUrban = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pUrban)
    estRural = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCPB$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCPB$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCPB$countyMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCPB$countyMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPB, county, stratified):")
    fullPooledScoresCountyLCPB = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLCPB)
    print("Rural scores (LCPB, county, stratified):")
    fullRuralScoresCountyLCPB = data.frame(c(getScores(truthFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLCPB)
    print("Urban scores (LCPB, county, stratified):")
    fullUrbanScoresCountyLCPB = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLCPB)
    
    # do the same for the direct estimate validation
    est = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$p)
    estUrban = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban)
    estRural = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit4$aggregatedResultsLCPB$countyMatrices$p
    # vars = apply(estMat, 1, var)
    vars = rowMeans(expit(sweep(logit(estMat) + matrix(rnorm(n=length(estMat), sd=sqrt(truthDirectFullLogitVar)), ncol=ncol(estMat)), 1, logit(est), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultsLCPB$countyMatrices$pUrban
    varsUrban = rowMeans(expit(sweep(logit(estMatUrban) + matrix(rnorm(n=length(estMatUrban), sd=sqrt(truthDirectFullUrbanLogitVar)), ncol=ncol(estMatUrban)), 1, logit(estUrban), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultsLCPB$countyMatrices$pRural
    varsRural = rowMeans(expit(sweep(logit(estMatRural) + matrix(rnorm(n=length(estMatRural), sd=sqrt(truthDirectFullRuralLogitVar)), ncol=ncol(estMatRural)), 1, logit(estRural), FUN="-"))^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPB, county, direct):")
    fullPooledScoresCountyLCPBDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLCPBDirect)
    print("Rural scores (LCPB, county, direct):")
    fullRuralScoresCountyLCPBDirect = data.frame(c(getScores(truthDirectFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar[truthTableFullCounty$urbanProp!=1]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLCPBDirect)
    print("Urban scores (LCPB, county, direct):")
    fullUrbanScoresCountyLCPBDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLCPBDirect)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCPb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCPb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCPb$countyMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCPb$countyMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPb, county):")
    fullPooledScoresCountyLCPb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLCPb)
    print("Rural scores (LCPb, county):")
    fullRuralScoresCountyLCPb = data.frame(c(getScores(truthFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLCPb)
    print("Urban scores (LCPb, county):")
    fullUrbanScoresCountyLCPb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLCPb)
    
    # do the same for the direct estimate validation
    est = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban) * truthTableFullCounty$urbanProp + rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pRural) * (1 - truthTableFullCounty$urbanProp)
    estUrban = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban)
    estRural = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit4$aggregatedResultsLCPb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit4$aggregatedResultsLCPb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(expit(sweep(logit(estMat) + matrix(rnorm(n=length(estMat), sd=sqrt(truthDirectFullLogitVar)), ncol=ncol(estMat)), 1, logit(est), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultsLCPb$countyMatrices$pUrban
    varsUrban = rowMeans(expit(sweep(logit(estMatUrban) + matrix(rnorm(n=length(estMatUrban), sd=sqrt(truthDirectFullUrbanLogitVar)), ncol=ncol(estMatUrban)), 1, logit(estUrban), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultsLCPb$countyMatrices$pRural
    varsRural = rowMeans(expit(sweep(logit(estMatRural) + matrix(rnorm(n=length(estMatRural), sd=sqrt(truthDirectFullRuralLogitVar)), ncol=ncol(estMatRural)), 1, logit(estRural), FUN="-"))^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPb, county, direct):")
    fullPooledScoresCountyLCPbDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLCPbDirect)
    print("Rural scores (LCPb, county, direct):")
    fullRuralScoresCountyLCPbDirect = data.frame(c(getScores(truthDirectFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar[truthTableFullCounty$urbanProp!=1]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLCPbDirect)
    print("Urban scores (LCPb, county, direct):")
    fullUrbanScoresCountyLCPbDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLCPbDirect)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCpb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCpb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCpb$countyMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCpb$countyMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCpb, county):")
    fullPooledScoresCountyLCpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLCpb)
    print("Rural scores (LCpb, county):")
    fullRuralScoresCountyLCpb = data.frame(c(getScores(truthFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLCpb)
    print("Urban scores (LCpb, county):")
    fullUrbanScoresCountyLCpb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLCpb)
    
    # do the same for the direct estimate validation
    est = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban) * truthTableFullCounty$urbanProp + rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pRural) * (1 - truthTableFullCounty$urbanProp)
    estUrban = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban)
    estRural = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit4$aggregatedResultsLCpb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit4$aggregatedResultsLCpb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(expit(sweep(logit(estMat) + matrix(rnorm(n=length(estMat), sd=sqrt(truthDirectFullLogitVar)), ncol=ncol(estMat)), 1, logit(est), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultsLCpb$countyMatrices$pUrban
    varsUrban = rowMeans(expit(sweep(logit(estMatUrban) + matrix(rnorm(n=length(estMatUrban), sd=sqrt(truthDirectFullUrbanLogitVar)), ncol=ncol(estMatUrban)), 1, logit(estUrban), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultsLCpb$countyMatrices$pRural
    varsRural = rowMeans(expit(sweep(logit(estMatRural) + matrix(rnorm(n=length(estMatRural), sd=sqrt(truthDirectFullRuralLogitVar)), ncol=ncol(estMatRural)), 1, logit(estRural), FUN="-"))^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCpb, county, direct):")
    fullPooledScoresCountyLCpbDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLCpbDirect)
    print("Rural scores (LCpb, county, direct):")
    fullRuralScoresCountyLCpbDirect = data.frame(c(getScores(truthDirectFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar[truthTableFullCounty$urbanProp!=1]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLCpbDirect)
    print("Urban scores (LCpb, county, direct):")
    fullUrbanScoresCountyLCpbDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLCpbDirect)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLcpb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLcpb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLcpb$countyMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLcpb$countyMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (Lcpb, county):")
    fullPooledScoresCountyLcpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLcpb)
    print("Rural scores (Lcpb, county):")
    fullRuralScoresCountyLcpb = data.frame(c(getScores(truthFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLcpb)
    print("Urban scores (Lcpb, county):")
    fullUrbanScoresCountyLcpb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLcpb)
    
    # do the same for the direct estimate validation
    est = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban) * truthTableFullCounty$urbanProp + rowMeans(fit3$aggregatedResultslcpb$countyMatrices$pRural) * (1 - truthTableFullCounty$urbanProp)
    estUrban = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban)
    estRural = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit4$aggregatedResultsLcpb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit4$aggregatedResultsLcpb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(expit(sweep(logit(estMat) + matrix(rnorm(n=length(estMat), sd=sqrt(truthDirectFullLogitVar)), ncol=ncol(estMat)), 1, logit(est), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultsLcpb$countyMatrices$pUrban
    varsUrban = rowMeans(expit(sweep(logit(estMatUrban) + matrix(rnorm(n=length(estMatUrban), sd=sqrt(truthDirectFullUrbanLogitVar)), ncol=ncol(estMatUrban)), 1, logit(estUrban), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultsLcpb$countyMatrices$pRural
    varsRural = rowMeans(expit(sweep(logit(estMatRural) + matrix(rnorm(n=length(estMatRural), sd=sqrt(truthDirectFullRuralLogitVar)), ncol=ncol(estMatRural)), 1, logit(estRural), FUN="-"))^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (Lcpb, county, direct):")
    fullPooledScoresCountyLcpbDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountyLcpbDirect)
    print("Rural scores (Lcpb, county, direct):")
    fullRuralScoresCountyLcpbDirect = data.frame(c(getScores(truthDirectFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar[truthTableFullCounty$urbanProp!=1]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountyLcpbDirect)
    print("Urban scores (Lcpb, county, direct):")
    fullUrbanScoresCountyLcpbDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountyLcpbDirect)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$countyMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultslcpb$countyMatrices$pUrban, 1, truthTableFullCounty$urbanProp, FUN="*") + sweep(fit3$aggregatedResultslcpb$countyMatrices$pRural, 1, 1-truthTableFullCounty$urbanProp, FUN="*")
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultslcpb$countyMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultslcpb$countyMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (lcpb, county):")
    fullPooledScoresCountylcpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountylcpb)
    print("Rural scores (lcpb, county):")
    fullRuralScoresCountylcpb = data.frame(c(getScores(truthFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountylcpb)
    print("Urban scores (lcpb, county):")
    fullUrbanScoresCountylcpb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountylcpb)
    
    # do the same for the direct estimate validation
    est = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$p)
    estUrban = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pUrban)
    estRural = rowMeans(fit4$aggregatedResultslcpb$countyMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit4$aggregatedResultslcpb$countyMatrices$p
    # vars = apply(estMat, 1, var)
    vars = rowMeans(expit(sweep(logit(estMat) + matrix(rnorm(n=length(estMat), sd=sqrt(truthDirectFullLogitVar)), ncol=ncol(estMat)), 1, logit(est), FUN="-"))^2)
    estMatUrban = fit4$aggregatedResultslcpb$countyMatrices$pUrban
    varsUrban = rowMeans(expit(sweep(logit(estMatUrban) + matrix(rnorm(n=length(estMatUrban), sd=sqrt(truthDirectFullUrbanLogitVar)), ncol=ncol(estMatUrban)), 1, logit(estUrban), FUN="-"))^2)
    estMatRural = fit4$aggregatedResultslcpb$countyMatrices$pRural
    varsRural = rowMeans(expit(sweep(logit(estMatRural) + matrix(rnorm(n=length(estMatRural), sd=sqrt(truthDirectFullRuralLogitVar)), ncol=ncol(estMatRural)), 1, logit(estRural), FUN="-"))^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (lcpb, county, direct):")
    fullPooledScoresCountylcpbDirect = data.frame(c(getScores(truthDirectFull, est, vars, lower, upper, estMat, doRandomReject=TRUE, logitTruthVar=truthDirectFullLogitVar), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresCountylcpbDirect)
    print("Rural scores (lcpb, county, direct):")
    fullRuralScoresCountylcpbDirect = data.frame(c(getScores(truthDirectFullRural[truthTableFullCounty$urbanProp!=1], estRural[truthTableFullCounty$urbanProp!=1], varsRural[truthTableFullCounty$urbanProp!=1], lower, upper, estMatRural[truthTableFullCounty$urbanProp!=1,], doRandomReject=TRUE, logitTruthVar=truthDirectFullRuralLogitVar[truthTableFullCounty$urbanProp!=1]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresCountylcpbDirect)
    print("Urban scores (lcpb, county, direct):")
    fullUrbanScoresCountylcpbDirect = data.frame(c(getScores(truthDirectFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE, logitTruthVar=truthDirectFullUrbanLogitVar), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresCountylcpbDirect)
    
    fullScoresCounty = list(fullPooledScoreslcpb=fullPooledScoresCountylcpb, fullRuralScoreslcpb=fullRuralScoresCountylcpb, fullUrbanScoreslcpb=fullUrbanScoresCountylcpb, 
                            fullPooledScoresLcpb=fullPooledScoresCountyLcpb, fullRuralScoresLcpb=fullRuralScoresCountyLcpb, fullUrbanScoresLcpb=fullUrbanScoresCountyLcpb, 
                            fullPooledScoresLCpb=fullPooledScoresCountyLCpb, fullRuralScoresLCpb=fullRuralScoresCountyLCpb, fullUrbanScoresLCpb=fullUrbanScoresCountyLCpb, 
                            fullPooledScoresLCPb=fullPooledScoresCountyLCPb, fullRuralScoresLCPb=fullRuralScoresCountyLCPb, fullUrbanScoresLCPb=fullUrbanScoresCountyLCPb, 
                            fullPooledScoresLCPB=fullPooledScoresCountyLCPB, fullRuralScoresLCPB=fullRuralScoresCountyLCPB, fullUrbanScoresLCPB=fullUrbanScoresCountyLCPB, 
                            fullPooledScoreslcpbDirect=fullPooledScoresCountylcpbDirect, fullRuralScoreslcpbDirect=fullRuralScoresCountylcpbDirect, fullUrbanScoreslcpbDirect=fullUrbanScoresCountylcpbDirect, 
                            fullPooledScoresLcpbDirect=fullPooledScoresCountyLcpbDirect, fullRuralScoresLcpbDirect=fullRuralScoresCountyLcpbDirect, fullUrbanScoresLcpbDirect=fullUrbanScoresCountyLcpbDirect, 
                            fullPooledScoresLCpbDirect=fullPooledScoresCountyLCpbDirect, fullRuralScoresLCpbDirect=fullRuralScoresCountyLCpbDirect, fullUrbanScoresLCpbDirect=fullUrbanScoresCountyLCpbDirect, 
                            fullPooledScoresLCPbDirect=fullPooledScoresCountyLCPbDirect, fullRuralScoresLCPbDirect=fullRuralScoresCountyLCPbDirect, fullUrbanScoresLCPbDirect=fullUrbanScoresCountyLCPbDirect, 
                            fullPooledScoresLCPBDirect=fullPooledScoresCountyLCPBDirect, fullRuralScoresLCPBDirect=fullRuralScoresCountyLCPBDirect, fullUrbanScoresLCPBDirect=fullUrbanScoresCountyLCPBDirect)
    
    ## Province level
    truthTableFullProvince = getProvinceLevelTruth(dat=dat, easpa=easpa, targetPop="children")
    
    # get observations and prediction summary statistics
    # easpaFull = makeDefaultEASPA(validationClusterI=NULL, useClustersAsEAs=TRUE)
    truthFull = truthTableFullProvince$pStratified
    truthFullUrban = truthTableFullProvince$pUrban
    truthFullRural = truthTableFullProvince$pRural
    # obsUrban = truthTableFullProvince$urban
    
    cpo = fit$mod$cpo$cpo
    cpoFailure = fit$mod$cpo$failure
    dic = fit$mod$dic$dic
    waic = fit$mod$waic$waic
    modelFit = fit$mod
    
    est = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$pUrban) * truthTableFullProvince$urbanProp + rowMeans(fit3$aggregatedResultslcpb$regionMatrices$pRural) * (1 - truthTableFullProvince$urbanProp)
    estUrban = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$pUrban)
    estRural = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$pRural)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCPB$regionMatrices$pUrban, 1, truthTableFullProvince$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCPB$regionMatrices$pRural, 1, 1-truthTableFullProvince$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCPB$regionMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCPB$regionMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPB, province):")
    fullPooledScoresProvinceLCPB = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresProvinceLCPB)
    print("Rural scores (LCPB, province):")
    fullRuralScoresProvinceLCPB = data.frame(c(getScores(truthFullRural[truthTableFullProvince$urbanProp!=1], estRural[truthTableFullProvince$urbanProp!=1], varsRural[truthTableFullProvince$urbanProp!=1], lower, upper, estMatRural[truthTableFullProvince$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresProvinceLCPB)
    print("Urban scores (LCPB, province):")
    fullUrbanScoresProvinceLCPB = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresProvinceLCPB)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCPb$regionMatrices$pUrban, 1, truthTableFullProvince$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCPb$regionMatrices$pRural, 1, 1-truthTableFullProvince$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCPb$regionMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCPb$regionMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCPb, province):")
    fullPooledScoresProvinceLCPb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresProvinceLCPb)
    print("Rural scores (LCPb, province):")
    fullRuralScoresProvinceLCPb = data.frame(c(getScores(truthFullRural[truthTableFullProvince$urbanProp!=1], estRural[truthTableFullProvince$urbanProp!=1], varsRural[truthTableFullProvince$urbanProp!=1], lower, upper, estMatRural[truthTableFullProvince$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresProvinceLCPb)
    print("Urban scores (LCPb, province):")
    fullUrbanScoresProvinceLCPb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresProvinceLCPb)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLCpb$regionMatrices$pUrban, 1, truthTableFullProvince$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLCpb$regionMatrices$pRural, 1, 1-truthTableFullProvince$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLCpb$regionMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLCpb$regionMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (LCpb, province):")
    fullPooledScoresProvinceLCpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresProvinceLCpb)
    print("Rural scores (LCpb, province):")
    fullRuralScoresProvinceLCpb = data.frame(c(getScores(truthFullRural[truthTableFullProvince$urbanProp!=1], estRural[truthTableFullProvince$urbanProp!=1], varsRural[truthTableFullProvince$urbanProp!=1], lower, upper, estMatRural[truthTableFullProvince$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresProvinceLCpb)
    print("Urban scores (LCpb, province):")
    fullUrbanScoresProvinceLCpb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresProvinceLCpb)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultsLcpb$regionMatrices$pUrban, 1, truthTableFullProvince$urbanProp, FUN="*") + sweep(fit3$aggregatedResultsLcpb$regionMatrices$pRural, 1, 1-truthTableFullProvince$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultsLcpb$regionMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultsLcpb$regionMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (Lcpb, province):")
    fullPooledScoresProvinceLcpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresProvinceLcpb)
    print("Rural scores (Lcpb, province):")
    fullRuralScoresProvinceLcpb = data.frame(c(getScores(truthFullRural[truthTableFullProvince$urbanProp!=1], estRural[truthTableFullProvince$urbanProp!=1], varsRural[truthTableFullProvince$urbanProp!=1], lower, upper, estMatRural[truthTableFullProvince$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresProvinceLcpb)
    print("Urban scores (Lcpb, province):")
    fullUrbanScoresProvinceLcpb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresProvinceLcpb)
    
    # est = rowMeans(fit3$aggregatedResultslcpb$regionMatrices$p)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = sweep(fit3$aggregatedResultslcpb$regionMatrices$pUrban, 1, truthTableFullProvince$urbanProp, FUN="*") + sweep(fit3$aggregatedResultslcpb$regionMatrices$pRural, 1, 1-truthTableFullProvince$urbanProp, FUN="*")
    # vars = apply(estMat, 1, var)
    vars = rowMeans(sweep(estMat, 1, est, FUN="-")^2)
    estMatUrban = fit3$aggregatedResultslcpb$regionMatrices$pUrban
    varsUrban = rowMeans(sweep(estMatUrban, 1, estUrban, FUN="-")^2)
    estMatRural = fit3$aggregatedResultslcpb$regionMatrices$pRural
    varsRural = rowMeans(sweep(estMatRural, 1, estRural, FUN="-")^2)
    
    # # calculate validation scoring rules
    print("Pooled scores (lcpb, province):")
    fullPooledScoresProvincelcpb = data.frame(c(getScores(truthFull, est, vars, lower, upper, estMat, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullPooledScoresProvincelcpb)
    print("Rural scores (lcpb, province):")
    fullRuralScoresProvincelcpb = data.frame(c(getScores(truthFullRural[truthTableFullProvince$urbanProp!=1], estRural[truthTableFullProvince$urbanProp!=1], varsRural[truthTableFullProvince$urbanProp!=1], lower, upper, estMatRural[truthTableFullProvince$urbanProp!=1,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullRuralScoresProvincelcpb)
    print("Urban scores (lcpb, province):")
    fullUrbanScoresProvincelcpb = data.frame(c(getScores(truthFullUrban, estUrban, varsUrban, lower, upper, estMatUrban, doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince))
    print(fullUrbanScoresProvincelcpb)
    
    fullScoresProvince = list(fullPooledScoreslcpb=fullPooledScoresProvincelcpb, fullRuralScoreslcpb=fullRuralScoresProvincelcpb, fullUrbanScoreslcpb=fullUrbanScoresProvincelcpb, 
                              fullPooledScoresLcpb=fullPooledScoresProvinceLcpb, fullRuralScoresLcpb=fullRuralScoresProvinceLcpb, fullUrbanScoresLcpb=fullUrbanScoresProvinceLcpb, 
                              fullPooledScoresLCpb=fullPooledScoresProvinceLCpb, fullRuralScoresLCpb=fullRuralScoresProvinceLCpb, fullUrbanScoresLCpb=fullUrbanScoresProvinceLCpb, 
                              fullPooledScoresLCPb=fullPooledScoresProvinceLCPb, fullRuralScoresLCPb=fullRuralScoresProvinceLCPb, fullUrbanScoresLCPb=fullUrbanScoresProvinceLCPb, 
                              fullPooledScoresLCPB=fullPooledScoresProvinceLCPB, fullRuralScoresLCPB=fullRuralScoresProvinceLCPB, fullUrbanScoresLCPB=fullUrbanScoresProvinceLCPB)
    
    if(saveResults) {
      print(paste0("Saving full model results to: ", fileName))
      save(truthTableFull=truthTableFull, truthTableFullCounty=truthTableFullCounty, truthTableFullProvince=truthTableFullProvince, 
           timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince, 
           fit=fit, fit2=fit2, fit3=fit3, 
           fullScoresPixel, fullScoresCounty, fullScoresProvince, 
           file=fileName)
    }
  }
  else {
    print("Loading previous full model fit")
    load(fileName)
  }
  previousFit = fit
  browser()
  # set up sample table of indices if using stratified validation
  if(stratifiedValidation && is.null(sampleTable)) {
    out = getValidationI(dat=dat, dataType=dataType, pixelLevel=TRUE)
    sampleTable = out$sampleMatrix
    sampleListPixel = out$pixelIndexList
    # sampleMatrix
    # pixelIndexList
  }
  
  # get region names
  allRegions = countyToRegion(dat$admin1)
  regions = sort(unique(allRegions))
  if(!stratifiedValidation)
    nFold = length(regions)
  else
    nFold = 10
  
  # calculate bins for nearest neighbor distances
  distanceMax = 0
  for(i in 1:nFold) {
    if(!stratifiedValidation) {
      thisRegion = regions[i]
      thisSampleI = allRegions == thisRegion
      leaveOutI = NULL
    } else {
      thisRegion = NULL
      leaveOutI = sampleTable[,i]
      thisSampleI = leaveOutI
    }
    
    ##### Break scores down by distance
    outOfSamplePixelIndices = sampleListPixel[[i]]
    inSamplePixelIndices = truthTableFull$pixelI[-match(outOfSamplePixelIndices, truthTableFull$pixelI)]
    predPts = cbind(popMat$east[outOfSamplePixelIndices], popMat$north[outOfSamplePixelIndices])
    obsCoords = cbind(popMat$east[inSamplePixelIndices], popMat$north[inSamplePixelIndices])
    predUrban = popMat$urban[outOfSamplePixelIndices]
    obsUrban = popMat$urban[inSamplePixelIndices]
    
    # first calculate all distances, broken down by urban, rural, and all aggregated observations
    distMatuu = rdist(obsCoords[!obsUrban,], predPts[!predUrban,])
    distMatuU = rdist(obsCoords[!obsUrban,], predPts[predUrban,])
    distMatUu = rdist(obsCoords[obsUrban,], predPts[!predUrban,])
    distMatUU = rdist(obsCoords[obsUrban,], predPts[predUrban,])
    distMatAu = rdist(obsCoords, predPts[!predUrban,])
    distMatAU = rdist(obsCoords, predPts[predUrban,])
    distMatuA = rdist(obsCoords[!obsUrban,], predPts)
    distMatUA = rdist(obsCoords[obsUrban,], predPts)
    distMatAA = rdist(obsCoords, predPts)
    
    # now calculate nearest distances
    nndistsuu = apply(distMatuu, 2, function(x) {min(x[x != 0])})
    nndistsuU = apply(distMatuU, 2, function(x) {min(x[x != 0])})
    nndistsUu = apply(distMatUu, 2, function(x) {min(x[x != 0])})
    nndistsUU = apply(distMatUU, 2, function(x) {min(x[x != 0])})
    nndistsAu = apply(distMatAu, 2, function(x) {min(x[x != 0])})
    nndistsAU = apply(distMatAU, 2, function(x) {min(x[x != 0])})
    nndistsuA = apply(distMatuA, 2, function(x) {min(x[x != 0])})
    nndistsUA = apply(distMatUA, 2, function(x) {min(x[x != 0])})
    nndistsAA = apply(distMatAA, 2, function(x) {min(x[x != 0])})
    tempMax = c(nndistsuu, nndistsuU, nndistsUu, nndistsUU, nndistsAu, nndistsAU, nndistsuA, nndistsuA, nndistsUA, nndistsUA, nndistsAA, nndistsAA)
    distanceMax = max(distanceMax, tempMax)
  }
  distanceBreaks = seq(0, distanceMax+1, l=20)
  
  
  
  appendScores = function(thisEst, thisPredMatrices, thisTruthTab, thisScoreTabList, aggregationLevel=c("cluster", "pixel", "county", "province"), printRoot="") {
    
    aggregationLevel = match.arg(aggregationLevel)
    if(aggregationLevel %in% c("cluster", "pixel")) {
      stratifiedByUrban = TRUE
    } else {
      stratifiedByUrban = FALSE
    }
    
    # get observations and prediction summary statistics
    thisTruth = thisTruthTab$p
    obsUrban = thisTruthTab$urban
    est = thisEst
    vars = apply(thisPredMatrices$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = thisPredMatrices$p
    
    # calculate validation scoring rules
    print(paste0(printRoot, " pooled scores:"))
    if(!stratifiedValidation)
      thisPooledScores = data.frame(c(list(Region=thisRegion), getScores(thisTruth, est, vars, lower, upper, estMat, doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
    else
      thisPooledScores = data.frame(c(list(Fold=i), getScores(thisTruth, est, vars, lower, upper, estMat, doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
    print(thisPooledScores)
    
    # calculate stratified validation scoring rules by urban/rural if necessary
    if(stratifiedByUrban) {
      if(stratifiedValidation || thisRegion != "Nairobi") {
        print(paste0(printRoot, " rural scores:"))
        if(!stratifiedValidation)
          thisRuralScores = data.frame(c(list(Region=thisRegion), getScores(thisTruth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
        else
          thisRuralScores = data.frame(c(list(Fold=i), getScores(thisTruth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
        print(thisRuralScores)
      } else {
        thisRuralScores = thisPooledScores
        thisRuralScores[,2:(ncol(thisRuralScores)-1)] = NA
        thisRuralScores = thisRuralScores
      }
      
      print(paste0(printRoot, " urban scores:"))
      if(!stratifiedValidation)
        thisUrbanScores = data.frame(c(list(Region=thisRegion), getScores(thisTruth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
      else
        thisUrbanScores = data.frame(c(list(Fold=i), getScores(thisTruth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
      print(thisUrbanScores)
    }
    
    
    # append scoring rule tables
    thisScoreTabList$completeScoreTable = rbind(thisScoreTabList$completeScoreTable, thisPooledScores)
    thisScoreTabList$pooledScoreTable = rbind(thisScoreTabList$pooledScoreTable, thisPooledScores)
    
    if(stratifiedByUrban) {
      thisScoreTabList$completeScoreTable = rbind(thisScoreTabList$completeScoreTable, thisRuralScores)
      thisScoreTabList$completeScoreTable = rbind(thisScoreTabList$completeScoreTable, thisUrbanScores)
      
      thisScoreTabList$ruralScoreTable = rbind(thisScoreTabList$ruralScoreTable, thisRuralScores)
      thisScoreTabList$urbanScoreTable = rbind(thisScoreTabList$urbanScoreTable, thisUrbanScores)
    }
    
    ##### Break scores down by distance if necessary
    if(stratifiedByUrban) {
      outOfSamplePixelIndices = sampleListPixel[[i]]
      inSamplePixelIndices = truthTableFull$pixelI[-match(outOfSamplePixelIndices, truthTableFull$pixelI)]
      predPts = cbind(popMat$east[outOfSamplePixelIndices], popMat$north[outOfSamplePixelIndices])
      obsCoords = cbind(popMat$east[inSamplePixelIndices], popMat$north[inSamplePixelIndices])
      predUrban = popMat$urban[outOfSamplePixelIndices]
      obsUrban = popMat$urban[inSamplePixelIndices]
      
      # first calculate all distances, broken down by urban, rural, and all aggregated observations
      distMatuU = rdist(obsCoords[!obsUrban,], predPts[predUrban,])
      distMatUU = rdist(obsCoords[obsUrban,], predPts[predUrban,])
      distMatAU = rdist(obsCoords, predPts[predUrban,])
      distMatuA = rdist(obsCoords[!obsUrban,], predPts)
      distMatUA = rdist(obsCoords[obsUrban,], predPts)
      distMatAA = rdist(obsCoords, predPts)
      
      # now calculate nearest distances
      nndistsuU = apply(distMatuU, 2, function(x) {min(x[x != 0])})
      nndistsUU = apply(distMatUU, 2, function(x) {min(x[x != 0])})
      nndistsAU = apply(distMatAU, 2, function(x) {min(x[x != 0])})
      nndistsuA = apply(distMatuA, 2, function(x) {min(x[x != 0])})
      nndistsUA = apply(distMatUA, 2, function(x) {min(x[x != 0])})
      nndistsAA = apply(distMatAA, 2, function(x) {min(x[x != 0])})
      if(stratifiedValidation || thisRegion != "Nairobi") {
        distMatuu = rdist(obsCoords[!obsUrban,], predPts[!predUrban,])
        distMatUu = rdist(obsCoords[obsUrban,], predPts[!predUrban,])
        distMatAu = rdist(obsCoords, predPts[!predUrban,])
        
        nndistsuu = apply(distMatuu, 2, function(x) {min(x[x != 0])})
        nndistsUu = apply(distMatUu, 2, function(x) {min(x[x != 0])})
        nndistsAu = apply(distMatAu, 2, function(x) {min(x[x != 0])})
        
        binnedScoringRulesuu = getScores(thisTruth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsuu, breaks=distanceBreaks)$binnedResults
        binnedScoringRulesUu = getScores(thisTruth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsUu, breaks=distanceBreaks)$binnedResults
        binnedScoringRulesAu = getScores(thisTruth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsAu, breaks=distanceBreaks)$binnedResults
      } else {
        binnedScoringRulesuu = NULL
        binnedScoringRulesUu = NULL
        binnedScoringRulesAu = NULL
      }
      
      # calculate scores accounting for binomial variation using fuzzy reject intervals
      binnedScoringRulesuU = getScores(thisTruth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsuU, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesUU = getScores(thisTruth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsUU, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesAU = getScores(thisTruth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsAU, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesuA = getScores(thisTruth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsuA, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesUA = getScores(thisTruth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsUA, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesAA = getScores(thisTruth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsAA, breaks=distanceBreaks)$binnedResults
      
      # concatenate binned scoring rule results
      thisScoreTabList$binnedScoringRulesuuAll = c(thisScoreTabList$binnedScoringRulesuuAll, list(binnedScoringRulesuu))
      thisScoreTabList$binnedScoringRulesuUAll = c(thisScoreTabList$binnedScoringRulesuUAll, list(binnedScoringRulesuU))
      thisScoreTabList$binnedScoringRulesUuAll = c(thisScoreTabList$binnedScoringRulesUuAll, list(binnedScoringRulesUu))
      thisScoreTabList$binnedScoringRulesUUAll = c(thisScoreTabList$binnedScoringRulesUUAll, list(binnedScoringRulesUU))
      thisScoreTabList$binnedScoringRulesAuAll = c(thisScoreTabList$binnedScoringRulesAuAll, list(binnedScoringRulesAu))
      thisScoreTabList$binnedScoringRulesAUAll = c(thisScoreTabList$binnedScoringRulesAUAll, list(binnedScoringRulesAU))
      thisScoreTabList$binnedScoringRulesuAAll = c(thisScoreTabList$binnedScoringRulesuAAll, list(binnedScoringRulesuA))
      thisScoreTabList$binnedScoringRulesUAAll = c(thisScoreTabList$binnedScoringRulesUAAll, list(binnedScoringRulesUA))
      thisScoreTabList$binnedScoringRulesAAAll = c(thisScoreTabList$binnedScoringRulesAAAll, list(binnedScoringRulesAA))
    }
    
    ##### Calculate individual scoring rules
    # calculate the scoring rules, and add nearest neighbor distances for each stratum
    if(!stratifiedValidation) {
      if(stratifiedByUrban) {
        thisSingleScores = data.frame(c(list(Region=thisRegion, dataI=thisTruthTab[,1], NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(thisTruth, est, vars, lower, upper, estMat, getAverage=FALSE), timePoint=timePoint, timeAggregation=timeAggregation))
      } else {
        thisSingleScores = data.frame(c(list(Region=thisRegion, dataI=thisTruthTab[,1]), getScores(thisTruth, est, vars, lower, upper, estMat, getAverage=FALSE), timePoint=timePoint, timeAggregation=timeAggregation))
      }
    }
    else {
      if(stratifiedByUrban) {
        thisSingleScores = data.frame(c(list(Fold=i, dataI=thisTruthTab[,1], NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(thisTruth, est, vars, lower, upper, estMat, getAverage=FALSE), timePoint=timePoint, timeAggregation=timeAggregation))
      } else {
        thisSingleScores = data.frame(c(list(Fold=i, dataI=thisTruthTab[,1]), getScores(thisTruth, est, vars, lower, upper, estMat, getAverage=FALSE), timePoint=timePoint, timeAggregation=timeAggregation))
      }
    }
    
    # concatenate the results
    thisScoreTabList$singleScores = rbind(thisScoreTabList$singleScores, thisSingleScores)
    
    thisScoreTabList
  }
  
  scoresPixelLCPB = list(completeScoreTable = c(), pooledScoreTable = c(),urbanScoreTable = c(), ruralScoreTable = c(), 
                         binnedScoringRulesuuAll = list(), binnedScoringRulesuUAll = list(), binnedScoringRulesUuAll = list(), 
                         binnedScoringRulesUUAll = list(), binnedScoringRulesAuAll = list(), binnedScoringRulesAUAll = list(), 
                         binnedScoringRulesuAAll = list(), binnedScoringRulesUAAll = list(), binnedScoringRulesAAAll = list(), 
                         singleScores = c())
  scoresPixelLCPb = list(completeScoreTable = c(), pooledScoreTable = c(),urbanScoreTable = c(), ruralScoreTable = c(), 
                         binnedScoringRulesuuAll = list(), binnedScoringRulesuUAll = list(), binnedScoringRulesUuAll = list(), 
                         binnedScoringRulesUUAll = list(), binnedScoringRulesAuAll = list(), binnedScoringRulesAUAll = list(), 
                         binnedScoringRulesuAAll = list(), binnedScoringRulesUAAll = list(), binnedScoringRulesAAAll = list(), 
                         singleScores = c())
  scoresPixelLCpb = list(completeScoreTable = c(), pooledScoreTable = c(),urbanScoreTable = c(), ruralScoreTable = c(), 
                         binnedScoringRulesuuAll = list(), binnedScoringRulesuUAll = list(), binnedScoringRulesUuAll = list(), 
                         binnedScoringRulesUUAll = list(), binnedScoringRulesAuAll = list(), binnedScoringRulesAUAll = list(), 
                         binnedScoringRulesuAAll = list(), binnedScoringRulesUAAll = list(), binnedScoringRulesAAAll = list(), 
                         singleScores = c())
  scoresPixelLcpb = list(completeScoreTable = c(), pooledScoreTable = c(),urbanScoreTable = c(), ruralScoreTable = c(), 
                         binnedScoringRulesuuAll = list(), binnedScoringRulesuUAll = list(), binnedScoringRulesUuAll = list(), 
                         binnedScoringRulesUUAll = list(), binnedScoringRulesAuAll = list(), binnedScoringRulesAUAll = list(), 
                         binnedScoringRulesuAAll = list(), binnedScoringRulesUAAll = list(), binnedScoringRulesAAAll = list(), 
                         singleScores = c())
  scoresPixellcpb = list(completeScoreTable = c(), pooledScoreTable = c(),urbanScoreTable = c(), ruralScoreTable = c(), 
                         binnedScoringRulesuuAll = list(), binnedScoringRulesuUAll = list(), binnedScoringRulesUuAll = list(), 
                         binnedScoringRulesUUAll = list(), binnedScoringRulesAuAll = list(), binnedScoringRulesAUAll = list(), 
                         binnedScoringRulesuAAll = list(), binnedScoringRulesUAAll = list(), binnedScoringRulesAAAll = list(), 
                         singleScores = c())
  
  startFrom = 1
  
  # load previous results if necessary
  fileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationAllTemp", 
                    "_urb", urbanEffect, ".RData")
  if(loadPreviousResults && file.exists(fileName)) {
    load(fileName)
    startFrom = i+1
  }
  
  for(i in startFrom:nFold) {
    if(!stratifiedValidation) {
      thisRegion = regions[i]
      thisSampleI = allRegions == thisRegion
      print(paste0("Validating SPDE model with urban=", urbanEffect, 
                   ", region=", thisRegion, " (", i, "/", length(regions), " regions)"))
      leaveOutI = NULL
    } else {
      thisRegion = NULL
      print(paste0("Validating SPDE model with urban=", urbanEffect, 
                   ", (", i, "/", nFold, " folds)"))
      leaveOutI = sampleTable[,i]
      thisSampleI = leaveOutI
    }
    thisPixelI = sampleListPixel[[i]]
    
    # subset the data
    thisDat = dat[!leaveOutI,]
    
    # fit the point level model
    timePoint = system.time(fit <- fitSPDEKenyaDat(thisDat, dataType, mesh, prior, significanceCI, int.strategy, strategy, nPostSamples, 
                                                   verbose, link, NULL, urbanEffect, clusterEffect, thisRegion,
                                                   kmres, TRUE, previousFit, family="binomial"))[3]
    
    # construction the faux population data frame for the left out data
    thiseaspa = makeDefaultEASPA(dataType2, validationClusterI=leaveOutI, useClustersAsEAs=TRUE)
    
    # get the aggregation model predictions
    thisTruthTable = truthTableFull[match(thisPixelI, truthTableFull$pixelI),]
    clustersPerPixel = rep(0, nrow(popMat))
    clustersPerPixel[thisTruthTable$pixelI] = thisTruthTable$nClusters
    
    timeAggregation = system.time(aggregatedPreds <- modLCPB(fit$uDraws, fit$sigmaEpsilonDraws, easpa=thiseaspa, popMat=popMat, empiricalDistributions=NULL, 
                                                             includeUrban=urbanEffect, maxEmptyFraction=1, clusterLevel=FALSE, pixelLevel=TRUE, countyLevel=FALSE, 
                                                             regionLevel=FALSE, nationalLevel=FALSE, doModifiedPixelLevel=FALSE, 
                                                             validationPixelI=thisPixelI, clustersPerPixel=clustersPerPixel, 
                                                             doLCPb=doLCPb, doLCpb=doLCpb, doLcpb=doLcpb))[3]
    
    # get observations and prediction summary statistics
    thisTruth = thisTruthTable$p
    obsUrban = thisTruthTable$urban
    est = rowMeans(aggregatedPreds$pixelMatriceslcpb$p)
    vars = apply(aggregatedPreds$pixelMatricesLCPB$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = aggregatedPreds$pixelMatricesLCPB$p
    
    scoresPixelLCPB = appendScores(thisEst, thisPredMatrices, thisTruthTab=thisTruthTable, scoresPixelLCPB, aggregationLevel="pixel", printRoot="LCPB")
    
    # calculate validation scoring rules
    print("Pooled scores:")
    if(!stratifiedValidation)
      thisPooledScores = data.frame(c(list(Region=thisRegion), getScores(thisTruth, est, vars, lower, upper, estMat, doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
    else
      thisPooledScores = data.frame(c(list(Fold=i), getScores(thisTruth, est, vars, lower, upper, estMat, doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
    print(thisPooledScores)
    
    if(stratifiedValidation || thisRegion != "Nairobi") {
      print("Rural scores:")
      if(!stratifiedValidation)
        thisRuralScores = data.frame(c(list(Region=thisRegion), getScores(thisTruth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
      else
        thisRuralScores = data.frame(c(list(Fold=i), getScores(thisTruth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
      print(thisRuralScores)
    } else {
      thisRuralScores = thisPooledScores
      thisRuralScores[,2:(ncol(thisRuralScores)-1)] = NA
      thisRuralScores = thisRuralScores
    }
    
    print("Urban scores:")
    if(!stratifiedValidation)
      thisUrbanScores = data.frame(c(list(Region=thisRegion), getScores(thisTruth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
    else
      thisUrbanScores = data.frame(c(list(Fold=i), getScores(thisTruth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), timePoint=timePoint, timeAggregation=timeAggregation))
    print(thisUrbanScores)
    
    # append scoring rule tables
    completeScoreTable = rbind(completeScoreTable, thisPooledScores)
    completeScoreTable = rbind(completeScoreTable, thisRuralScores)
    completeScoreTable = rbind(completeScoreTable, thisUrbanScores)
    
    pooledScoreTable = rbind(pooledScoreTable, thisPooledScores)
    ruralScoreTable = rbind(ruralScoreTable, thisRuralScores)
    urbanScoreTable = rbind(urbanScoreTable, thisUrbanScores)
    
    ##### Break scores down by distance
    outOfSamplePixelIndices = sampleListPixel[[i]]
    inSamplePixelIndices = truthTableFull$pixelI[-match(outOfSamplePixelIndices, truthTableFull$pixelI)]
    predPts = cbind(popMat$east[outOfSamplePixelIndices], popMat$north[outOfSamplePixelIndices])
    obsCoords = cbind(popMat$east[inSamplePixelIndices], popMat$north[inSamplePixelIndices])
    predUrban = popMat$urban[outOfSamplePixelIndices]
    obsUrban = popMat$urban[inSamplePixelIndices]
    
    # first calculate all distances, broken down by urban, rural, and all aggregated observations
    distMatuU = rdist(obsCoords[!obsUrban,], predPts[predUrban,])
    distMatUU = rdist(obsCoords[obsUrban,], predPts[predUrban,])
    distMatAU = rdist(obsCoords, predPts[predUrban,])
    distMatuA = rdist(obsCoords[!obsUrban,], predPts)
    distMatUA = rdist(obsCoords[obsUrban,], predPts)
    distMatAA = rdist(obsCoords, predPts)
    
    # now calculate nearest distances
    nndistsuU = apply(distMatuU, 2, function(x) {min(x[x != 0])})
    nndistsUU = apply(distMatUU, 2, function(x) {min(x[x != 0])})
    nndistsAU = apply(distMatAU, 2, function(x) {min(x[x != 0])})
    nndistsuA = apply(distMatuA, 2, function(x) {min(x[x != 0])})
    nndistsUA = apply(distMatUA, 2, function(x) {min(x[x != 0])})
    nndistsAA = apply(distMatAA, 2, function(x) {min(x[x != 0])})
    if(stratifiedValidation || thisRegion != "Nairobi") {
      distMatuu = rdist(obsCoords[!obsUrban,], predPts[!predUrban,])
      distMatUu = rdist(obsCoords[obsUrban,], predPts[!predUrban,])
      distMatAu = rdist(obsCoords, predPts[!predUrban,])
      
      nndistsuu = apply(distMatuu, 2, function(x) {min(x[x != 0])})
      nndistsUu = apply(distMatUu, 2, function(x) {min(x[x != 0])})
      nndistsAu = apply(distMatAu, 2, function(x) {min(x[x != 0])})
      
      binnedScoringRulesuu = getScores(thisTruth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsuu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesUu = getScores(thisTruth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsUu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesAu = getScores(thisTruth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsAu, breaks=distanceBreaks)$binnedResults
    } else {
      binnedScoringRulesuu = NULL
      binnedScoringRulesUu = NULL
      binnedScoringRulesAu = NULL
    }
    
    # calculate scores accounting for binomial variation using fuzzy reject intervals
    binnedScoringRulesuU = getScores(thisTruth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsuU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUU = getScores(thisTruth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsUU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAU = getScores(thisTruth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsAU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesuA = getScores(thisTruth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsuA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUA = getScores(thisTruth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsUA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAA = getScores(thisTruth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsAA, breaks=distanceBreaks)$binnedResults
    
    # concatenate binned scoring rule results
    binnedScoringRulesuuAll = c(binnedScoringRulesuuAll, list(binnedScoringRulesuu))
    binnedScoringRulesuUAll = c(binnedScoringRulesuUAll, list(binnedScoringRulesuU))
    binnedScoringRulesUuAll = c(binnedScoringRulesUuAll, list(binnedScoringRulesUu))
    binnedScoringRulesUUAll = c(binnedScoringRulesUUAll, list(binnedScoringRulesUU))
    binnedScoringRulesAuAll = c(binnedScoringRulesAuAll, list(binnedScoringRulesAu))
    binnedScoringRulesAUAll = c(binnedScoringRulesAUAll, list(binnedScoringRulesAU))
    binnedScoringRulesuAAll = c(binnedScoringRulesuAAll, list(binnedScoringRulesuA))
    binnedScoringRulesUAAll = c(binnedScoringRulesUAAll, list(binnedScoringRulesUA))
    binnedScoringRulesAAAll = c(binnedScoringRulesAAAll, list(binnedScoringRulesAA))
    
    ##### Calculate individual scoring rules
    # calculate the scoring rules, and add nearest neighbor distances for each stratum
    if(!stratifiedValidation) {
      thisSingleScores = data.frame(c(list(Region=thisRegion, dataI=thisPixelI, NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(thisTruth, est, vars, lower, upper, estMat, getAverage=FALSE), timePoint=timePoint, timeAggregation=timeAggregation))
    }
    else {
      thisSingleScores = data.frame(c(list(Fold=i, dataI=thisPixelI, NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(thisTruth, est, vars, lower, upper, estMat, getAverage=FALSE), timePoint=timePoint, timeAggregation=timeAggregation))
    }
    
    # concatenate the results
    singleScores = rbind(singleScores, thisSingleScores)
    
    # save results so far
    print(paste0("Saving validation fold ", i, "/", nFold, " results to: ", fileName))
    save(completeScoreTable, pooledScoreTable, ruralScoreTable, urbanScoreTable, 
         binnedScoringRulesuuAll, binnedScoringRulesuUAll, binnedScoringRulesUuAll, binnedScoringRulesUUAll, 
         binnedScoringRulesAuAll, binnedScoringRulesAUAll, binnedScoringRulesuAAll, binnedScoringRulesUAAll, 
         binnedScoringRulesAAAll, 
         singleScores, 
         i, file=fileName)
  }
  
  list(completeScoreTable=completeScoreTable, 
       pooledScoreTable=pooledScoreTable, 
       ruralScoreTable=ruralScoreTable, 
       urbanScoreTable=urbanScoreTable, 
       inSamplePooledScores=fullPooledScores, 
       inSampleUrbanScores=fullUrbanScores, 
       inSampleRuralScores=fullRuralScores, 
       
       binnedScoringRulesuuAll=averageBinnedScores(binnedScoringRulesuuAll), binnedScoringRulesuUAll=averageBinnedScores(binnedScoringRulesuUAll), 
       binnedScoringRulesUuAll=averageBinnedScores(binnedScoringRulesUuAll), binnedScoringRulesUUAll=averageBinnedScores(binnedScoringRulesUUAll), 
       binnedScoringRulesAuAll=averageBinnedScores(binnedScoringRulesAuAll), binnedScoringRulesAUAll=averageBinnedScores(binnedScoringRulesAUAll), 
       binnedScoringRulesuAAll=averageBinnedScores(binnedScoringRulesuAAll), binnedScoringRulesUAAll=averageBinnedScores(binnedScoringRulesUAAll), 
       binnedScoringRulesAAAll=averageBinnedScores(binnedScoringRulesAAAll), 
       
       singleScores=singleScores, 
       
       fullModelFit=previousFit)
}

# Use Poisson-Multinomial method to draw first from joint distribution of 
# number of households per EA given the total per stratum, and then from the joint 
# distribution of the number of people (children) per household given 
# the total per stratum
# n: the number of draws from the joint distribution of N
# eaPixelSamples: this is eaSamples from modLCPB (nPixel x n matrix of number of EAs per pixel)
# easpaList: a list of length n with each element being of the format given 
#            in makeDefaultEASPA() giving the number of households and EAs 
#            per stratum. It is assumed that the number of EAs per stratum is 
#            the same in each list element. If easpaList is a data frame, 
#            number of households per stratum is assumed constant
sampleNPoissonMultinomial = function(nDraws = ncol(pixelIndexMat), pixelIndexMat=NULL, urbanMat=NULL, areaMat=NULL, easpaList=NULL, 
                                     popMat=NULL, includeUrban=TRUE, verbose=TRUE, returnEAinfo=FALSE, useRcpp=FALSE) {
  
  # set default inputs
  if(is.null(easpaList)) {
    easpaList = lapply(1:nDraws, function(x){makeDefaultEASPA()})
    # area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
    # HHUrb: the number of households in the urban part of the area
    # HHRur: the number of households in the rural part of the area
    # HHTotal: the number of households in the the area
    # popUrb: the number of people in the urban part of the area
    # popRur: the number of people in the rural part of the area
    # popTotal: the number of people in the the area
    # pctTotal: the percentage of people in this area out of the total national population
    # pctUrb: the percentage of people in this area that live in urban areas
  } else if(length(easpaList) == 1) {
    easpaList = replicate(nDraws, easpaList[[1]], simplify=FALSE)
  }
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  if((is.null(areaMat) || is.null(urbanMat)) && is.null(pixelIndexMat)) {
    stop("user must either supply pixelIndexMat or both areaMat and urbanMat")
  }
  if(is.null(areaMat)) {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  }
  if(is.null(urbanMat)) {
    urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
  }
  
  # subset areaMat by urban/rural if necessary
  # if(includeUrban) {
  #   areaMatUrban = matrix(areaMat[urbanMat], ncol=nDraws)
  #   areaMatRural = matrix(areaMat[!urbanMat], ncol=nDraws)
  #   urbanStringMat = matrix(sapply(urbanMat, function(x) {ifelse(x, "u", "r")}), ncol=nDraws)
  #   areaUrbanicityMat = matrix(paste(areaMat, urbanStringMat, sep=","), ncol=nDraws)
  # }
  
  # start by drawing the totals, then divide households amongst EAs, then divide children amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the totals
  
  # get the total number of enumeration areas per stratum (this does not change between draws)
  areas = easpaList[[1]]$area
  totalEAsUrban = easpaList[[1]]$EAUrb
  totalEAsRural = easpaList[[1]]$EARur
  totalEAs = easpaList[[1]]$EATotal
  nEAs = sum(totalEAs)
  
  # ##### draw number of households per EA
  # 
  # # first preallocate the draws
  # if(includeUrban) {
  #   householdDrawsUrban = matrix(nrow=totalEAsUrban, ncol=nDraws)
  #   householdDrawsRural = matrix(nrow=totalEAsRural, ncol=nDraws)
  # } else {
  #   householdDraws = matrix(nrow=totalEAs, ncol=nDraws)
  # }
  # 
  # # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  # if(includeUrban) {
  #   totalHouseholdsUrban = sweep(sapply(easpaList, function(x) {x$HHUrb}), 1, -25*totalEAsUrban, "+")
  #   totalHouseholdsRural = sweep(sapply(easpaList, function(x) {x$HHRur}), 1, -25*totalEAsRural, "+")
  # } else {
  #   totalHouseholds = sweep(sapply(easpaList, function(x) {x$HHTotal}), 1, -25*totalEAs, "+")
  # }
  # 
  # # distribute the households throughout the enumeration areas with multinomial distribution
  # for(i in 1:length(areas)) {
  #   thisArea = areas[i]
  #   
  #   if(includeUrban) {
  #     householdDrawsUrban[areaMatUrban==thisArea] = rmultinom(nDraws, totalHouseholdsUrban[i], rep(1/totalEAsUrban[i], totalEAsUrban[i]))
  #     householdDrawsRural[areaMatRural==thisArea] = rmultinom(nDraws, totalHouseholdsRural[i], rep(1/totalEAsRural[i], totalEAsUrban[i]))
  #   } else {
  #     householdDraws[areaMat==thisArea] = rmultinom(nDraws, totalHouseholds[i], rep(1/totalEAs[i], totalEAs[i]))
  #   }
  # }
  # 
  # # add in the minimum 25 number of households
  # if(includeUrban) {
  #   householdDrawsUrban = householdDrawsUrban + 25
  #   householdDrawsRural = householdDrawsRural + 25
  # } else {
  #   householdDraws = householdDraws + 25
  # }
  
  ##### draw the number of children per enumeration area
  
  # first preallocate the draws
  # if(includeUrban) {
  #   childrenDrawsUrban = matrix(nrow=totalEAsUrban, ncol=nDraws)
  #   childrenDrawsRural = matrix(nrow=totalEAsRural, ncol=nDraws)
  # } else {
  #   childrenDraws = matrix(nrow=totalEAs, ncol=nDraws)
  # }
  childrenDraws = matrix(nrow=nEAs, ncol=nDraws)
  if(returnEAinfo) {
    householdDraws = matrix(nrow=nEAs, ncol=nDraws)
  }
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(includeUrban) {
    totalHouseholdsUrban = sweep(sapply(easpaList, function(x) {x$HHUrb}), 1, -25*totalEAsUrban, "+")
    totalHouseholdsRural = sweep(sapply(easpaList, function(x) {x$HHRur}), 1, -25*totalEAsRural, "+")
    totalChildrenUrban = sapply(easpaList, function(x) {x$popUrb})
    totalChildrenRural = sapply(easpaList, function(x) {x$popRur})
  } else {
    totalHouseholds = sweep(sapply(easpaList, function(x) {x$HHTotal}), 1, -25*totalEAs, "+")
    totalChildren = sapply(easpaList, function(x) {x$popHHTotal})
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the children amongst the households, also with a multinomial distribution
  for(i in 1:length(areas)) {
    thisArea = areas[i]
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(includeUrban) {
      householdDrawsUrban <- matrix(sapply(totalHouseholdsUrban[i,], rmultinom, n=1, prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])), nrow=totalEAsUrban[i], ncol=nDraws) + 25
      # householdDrawsUrban <- rmultinomSizeVec_rcpp(totalHouseholdsUrban[i,], rep(1/totalEAsUrban[i], totalEAsUrban[i])) + 25
      if(totalEAsRural[i] != 0) {
        householdDrawsRural <- matrix(sapply(totalHouseholdsRural[i,], rmultinom, n=1, prob=rep(1/totalEAsRural[i], totalEAsRural[i])), nrow=totalEAsRural[i], ncol=nDraws) + 25
        # householdDrawsRural <- rmultinomSizeVec_rcpp(totalHouseholdsRural[i,], rep(1/totalEAsRural[i], totalEAsRural[i])) + 25
      }
      
      # if we must return EA info, we must return the household draws for each EA:
      if(returnEAinfo) {
        householdDraws[areaMat==thisArea & urbanMat] = householdDrawsUrban
        if(totalEAsRural[i] != 0) {
          householdDraws[areaMat==thisArea & !urbanMat] = householdDrawsRural
        }
      }
    } else {
      householdDraws = matrix(sapply(totalHouseholds[i,], rmultinom, n=1, prob=rep(1/totalEAs[i], totalEAs[i])), nrow=totalEAs[i], ncol=nDraws) + 25
    }
    
    # drawing children per EA, with probability proportional to the number of households
    if(includeUrban) {
      probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
      # childrenDrawsUrban[areaMatUrban==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
      if(useRcpp) {
        childrenDraws[areaMat==thisArea & urbanMat] <- rmultinomSizeVecProbMat_rcpp(totalChildrenUrban[i,], probsUrban)
      } else {
        childrenDraws[areaMat==thisArea & urbanMat] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
      }
      
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        # childrenDrawsRural[areaMatRural==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
        if(useRcpp) {
          childrenDraws[areaMat==thisArea & !urbanMat] <- rmultinomSizeVecProbMat_rcpp(totalChildrenRural[i,], probsRural)
        } else {
          childrenDraws[areaMat==thisArea & !urbanMat] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
        }
        
      }
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      # childrenDraws[areaMat==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildren[i,j], probs[,j])})
      childrenDraws[areaMat==thisArea] = rmultinomSizeVecProbMat_rcpp(totalChildren[i,], probs)
    }
  }
  
  ##### Return results
  if(!returnEAinfo) {
    childrenDraws
  } else {
    list(householdDraws=householdDraws, childrenDraws=childrenDraws)
  }
}

# Use Poisson-Binomial method to draw first from marginal distributions of 
# number of households per EA given the total per stratum, and then from the joint 
# distribution of the number of people (children) per household given 
# the total per stratum
# n: the number of draws from the joint distribution of N
# eaPixelSamples: this is eaSamples from modLCPB (nPixel x n matrix of number of EAs per pixel)
# easpaList: a list of length n with each element being of the format given 
#            in makeDefaultEASPA() giving the number of households and EAs 
#            per stratum. It is assumed that the number of EAs per stratum is 
#            the same in each list element. If easpaList is a data frame, 
#            number of households per stratum is assumed constant
sampleNPoissonBinomial = function(eaSamplesMod, pixelIndexListMod, areaListMod, urbanListMod, nDraws=length(pixelIndexListMod), 
                                  easpa=NULL, easpaList=NULL, popMat=NULL, includeUrban=TRUE, verbose=TRUE, validationPixelI=NULL) {
  
  # set default inputs
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA(validationPixelI=validationPixelI)
    # area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
    # HHUrb: the number of households in the urban part of the area
    # HHRur: the number of households in the rural part of the area
    # HHTotal: the number of households in the the area
    # popUrb: the number of people in the urban part of the area
    # popRur: the number of people in the rural part of the area
    # popTotal: the number of people in the the area
    # pctTotal: the percentage of people in this area out of the total national population
    # pctUrb: the percentage of people in this area that live in urban areas
  }
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  if(is.null(easpaList)) {
    # nEAs = nEAsByStratum(areaListMod, urbanListMod)
    # easpaList = lapply(1:nDraws, function(j) {
    #   thiseaspa = easpa
    #   thiseaspa[,2:4] = nEAs[[j]][,1:3]
    #   thiseaspa
    # })
    
    nEAs = nEAsByStratum2(eaSamplesMod, popMat)
    easpaList = lapply(1:nDraws, function(j) {
      thiseaspa = easpa
      thiseaspa[,2:4] = nEAs[[j]][,2:4]
      thiseaspa
    })
  }
  
  # get popMat info
  popMatAreas = as.character(popMat$area)
  popMatUrban = popMat$urban
  areas = easpa$area
  
  ##### draw the number of children per enumeration area
  
  # initialize the children draws object
  childrenDraws = lapply(1:nDraws, function(j) {rep(NA, sum(nEAs[[j]]$EATotal))})
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the children amongst the households, also with a multinomial distribution
  for(i in 1:length(areas)) {
    thisArea = areas[i]
    
    thisPopMatI = popMatAreas == thisArea
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(includeUrban) {
      # get pixels for this area
      thisPopMatUrban = thisPopMatI & popMatUrban
      nUrbanPixels = sum(thisPopMatUrban)
      
      # calculate the probabilities of drawing a household in each urban and rural pixel
      urbanPixelHouseholdProbs = matrix(eaSamplesMod[thisPopMatUrban,], ncol=nDraws)
      urbanPixelHouseholdProbs = urbanPixelHouseholdProbs * (1/easpa$EAUrb[i])
      
      # draw households in each pixel
      householdDrawsUrbanPixel = 25*eaSamplesMod[thisPopMatUrban,] + matrix(rbinom(n=nUrbanPixels*nDraws, size=easpa$HHUrb[i] - 25*easpa$EAUrb[i], prob=urbanPixelHouseholdProbs), ncol=nDraws)
      
      # draw children in each pixel
      childrenDrawsUrbanPixel = matrix(rbinom(n=nUrbanPixels*nDraws, size=round(easpa$popUrb[i]), prob=householdDrawsUrbanPixel * (1 / easpa$HHUrb[i])), ncol=nDraws)
      
      # split children up amongst the EAs in the pixel for each pixel
      childrenDrawsUrban = lapply(1:nDraws, 
                                  function(j) {
                                    unlist(lapply(1:nUrbanPixels, function(i) {
                                      thisnEAs = eaSamplesMod[thisPopMatUrban,j][i]
                                      rmultinom(1, size=childrenDrawsUrbanPixel[i,j], prob=rep(1/thisnEAs, thisnEAs))}))
                                  })
      
      if(easpa$EARur[i] != 0) {
        ## do the same for the rural pixels in this area if need be
        
        thisPopMatRural = thisPopMatI & !popMatUrban
        nRuralPixels = sum(thisPopMatRural)
        ruralPixelHouseholdProbs = matrix(eaSamplesMod[thisPopMatRural,], ncol=nDraws)
        ruralPixelHouseholdProbs = ruralPixelHouseholdProbs * (1/easpa$EARur[i])
        
        householdDrawsRuralPixel = 25*eaSamplesMod[thisPopMatRural,] + matrix(rbinom(n=nRuralPixels*nDraws, size=easpa$HHRur[i] - 25*easpa$EARur[i], prob=ruralPixelHouseholdProbs), ncol=nDraws)
        childrenDrawsRuralPixel = matrix(rbinom(n=nRuralPixels*nDraws, size=round(easpa$popRur[i]), prob=householdDrawsRuralPixel * (1 / easpa$HHRur[i])), ncol=nDraws)
        childrenDrawsRural = lapply(1:nDraws, 
                                    function(j) {
                                      unlist(lapply(1:nRuralPixels, function(i) {
                                        thisnEAs = eaSamplesMod[thisPopMatRural,j][i]
                                        rmultinom(1, size=childrenDrawsRuralPixel[i,j], prob=rep(1/thisnEAs, thisnEAs))}))
                                    })
        
        # update childrenDraws object with both urban and rural children draws
        childrenDraws = lapply(1:nDraws, 
                               function(j) {
                                 childrenDraws[[j]][areaListMod[[j]] == thisArea & urbanListMod[[j]]] = childrenDrawsUrban[[j]]
                                 childrenDraws[[j]][areaListMod[[j]] == thisArea & !urbanListMod[[j]]] = childrenDrawsRural[[j]]
                                 childrenDraws[[j]]
                               })
      } else {
        # update childrenDraws object with only urban children draws
        childrenDraws = lapply(1:nDraws, 
                               function(j) {
                                 childrenDraws[[j]][areaListMod[[j]] == thisArea] = childrenDrawsUrban[[j]]
                                 childrenDraws[[j]]
                               })
      }
    } else {
      # get urban and rural pixels for this area
      thisPopMat = thisPopMatI
      nPixels = sum(thisPopMat)
      
      # calculate the probabilities of drawing a household in each urban and rural pixel
      pixelHouseholdProbs = matrix(eaSamplesMod[thisPopMat,], ncol=nDraws)
      pixelHouseholdProbs = pixelHouseholdProbs * (1/easpa$EATotal[i])
      
      # draw households in each pixel
      householdDrawsPixel = 25*eaSamplesMod[thisPopMat,] + matrix(rbinom(n=nPixels*nDraws, size=easpa$HHTotal[i] - 25*easpa$EATotal[i], prob=PixelHouseholdProbs), ncol=nDraws)
      
      # draw children in each pixel
      childrenDrawsPixel = matrix(rbinom(n=nPixels*nDraws, size=round(easpa$popTotal[i]), prob=householdDrawsPixel * (1 / easpa$HHTotal[i])), ncol=nDraws)
      
      # split children up amongst the EAs in the pixel for each pixel
      childrenDraws = lapply(1:nDraws, 
                                  function(j) {
                                    unlist(lapply(1:nPixels, function(i) {
                                      thisnEAs = eaSamplesMod[thisPopMat,j][i]
                                      rmultinom(1, size=childrenDrawsPixel[i,j], prob=rep(1/thisnEAs, thisnEAs))}))
                                  })
    }
  }
  
  ##### Return results
  childrenDraws
}

# Same as sampleNPoissonMultinomial, except the number of EAs per pixel is fixed
# n: the number of draws from the joint distribution of N
# eaPixelSamples: this is eaSamples from modLCPB (nPixel x n matrix of number of EAs per pixel)
# easpaList: a list of length n with each element being of the format given 
#            in makeDefaultEASPA() giving the number of households and EAs 
#            per stratum. It is assumed that the number of EAs per stratum is 
#            the same in each list element. If easpaList is a data frame, 
#            number of households per stratum is assumed constant
sampleNPoissonMultinomialFixed = function(clustersPerPixel, nDraws=ncol(pixelIndices), pixelIndices=NULL, 
                                          urbanVals=NULL, areaVals=NULL, easpa=NULL, popMat=NULL, includeUrban=TRUE, 
                                          verbose=TRUE) {
  
  # set default inputs
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA()
    # area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
    # HHUrb: the number of households in the urban part of the area
    # HHRur: the number of households in the rural part of the area
    # HHTotal: the number of households in the the area
    # popUrb: the number of people in the urban part of the area
    # popRur: the number of people in the rural part of the area
    # popTotal: the number of people in the the area
    # pctTotal: the percentage of people in this area out of the total national population
    # pctUrb: the percentage of people in this area that live in urban areas
  }
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  if((is.null(areaVals) || is.null(urbanVals)) && is.null(pixelIndices)) {
    stop("user must either supply pixelIndices or both areaVals and urbanVals")
  }
  if(is.null(areaVals)) {
    areaVals = matrix(popMat$area[pixelIndices], ncol=nDraws)
  }
  if(is.null(urbanVals)) {
    urbanVals = matrix(popMat$urban[pixelIndices], ncol=nDraws)
  }
  
  # # subset areaVals by urban/rural if necessary
  # if(includeUrban) {
  #   areaValsUrban = areaVals[urbanVals]
  #   areaValsRural = areaVals[!urbanVals]
  #   urbanStringVals = sapply(urbanVals, function(x) {ifelse(x, "u", "r")})
  #   areaUrbanicityVals = paste(areaVals, urbanStringVals, sep=",")
  # }
  
  # start by drawing the totals, then divide households amongst EAs, then divide children amongst households. 
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
  
  # ##### draw number of households per EA
  # 
  # # first preallocate the draws
  # if(includeUrban) {
  #   householdDrawsUrban = matrix(nrow=totalEAsUrban, ncol=nDraws)
  #   householdDrawsRural = matrix(nrow=totalEAsRural, ncol=nDraws)
  # } else {
  #   householdDraws = matrix(nrow=totalEAs, ncol=nDraws)
  # }
  # 
  # # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  # if(includeUrban) {
  #   totalHouseholdsUrban = sweep(sapply(easpaList, function(x) {x$HHUrb}), 1, -25*totalEAsUrban, "+")
  #   totalHouseholdsRural = sweep(sapply(easpaList, function(x) {x$HHRur}), 1, -25*totalEAsRural, "+")
  # } else {
  #   totalHouseholds = sweep(sapply(easpaList, function(x) {x$HHTotal}), 1, -25*totalEAs, "+")
  # }
  # 
  # # distribute the households throughout the enumeration areas with multinomial distribution
  # for(i in 1:length(areas)) {
  #   thisArea = areas[i]
  #   
  #   if(includeUrban) {
  #     householdDrawsUrban[areaMatUrban==thisArea] = rmultinom(nDraws, totalHouseholdsUrban[i], rep(1/totalEAsUrban[i], totalEAsUrban[i]))
  #     householdDrawsRural[areaMatRural==thisArea] = rmultinom(nDraws, totalHouseholdsRural[i], rep(1/totalEAsRural[i], totalEAsUrban[i]))
  #   } else {
  #     householdDraws[areaMat==thisArea] = rmultinom(nDraws, totalHouseholds[i], rep(1/totalEAs[i], totalEAs[i]))
  #   }
  # }
  # 
  # # add in the minimum 25 number of households
  # if(includeUrban) {
  #   householdDrawsUrban = householdDrawsUrban + 25
  #   householdDrawsRural = householdDrawsRural + 25
  # } else {
  #   householdDraws = householdDraws + 25
  # }
  
  ##### draw the number of children per enumeration area
  
  # first preallocate the draws
  # if(includeUrban) {
  #   childrenDrawsUrban = matrix(nrow=totalEAsUrban, ncol=nDraws)
  #   childrenDrawsRural = matrix(nrow=totalEAsRural, ncol=nDraws)
  # } else {
  #   childrenDraws = matrix(nrow=totalEAs, ncol=nDraws)
  # }
  childrenDraws = matrix(nrow=nEAs, ncol=nDraws)
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(includeUrban) {
    # totalHouseholdsUrban = sweep(sapply(easpaList, function(x) {x$HHUrb}), 1, -25*totalEAsUrban, "+")
    # totalHouseholdsRural = sweep(sapply(easpaList, function(x) {x$HHRur}), 1, -25*totalEAsRural, "+")
    # totalChildrenUrban = sapply(easpaList, function(x) {x$popUrb})
    # totalChildrenRural = sapply(easpaList, function(x) {x$popRur})
    totalHouseholdsUrban = easpa$HHUrb -25*totalEAsUrban
    totalHouseholdsRural = easpa$HHRur -25*totalEAsRural
    totalChildrenUrban = easpa$popUrb
    totalChildrenRural = easpa$popRur
  } else {
    # totalHouseholds = sweep(sapply(easpaList, function(x) {x$HHTotal}), 1, -25*totalEAs, "+")
    # totalChildren = sapply(easpaList, function(x) {x$popHHTotal})
    totalHouseholds = easpa$HHTotal - 25*totalEAs
    totalChildren = easpa$popHHTotal
  }
  
  # distribute the households throughout the enumeration areas with multinomial distribution, then 
  # distribute the children amongst the households, also with a multinomial distribution
  for(i in 1:length(areas)) {
    thisArea = areas[i]
    
    # print progress if in verbose mode
    if(verbose) {
      print(paste0("drawing Ns for each EA for area ", thisArea, " (", i, "/", length(areas), ")"))
    }
    
    # draw households per EA (make sure there are any rural EAs)
    if(includeUrban) {
      if(totalEAsUrban[i] != 0) {
        if(any(totalHouseholdsUrban != 0)) {
          householdDrawsUrban = rmultinom(n=nDraws, size=totalHouseholdsUrban[i], prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])) + 25
        } else {
          householdDrawsUrban = matrix(rep(25, totalEAsUrban[i]*nDraws), ncol=nDraws)
        }
      }
      
      if(totalEAsRural[i] != 0) {
        if(any(totalHouseholdsRural != 0)) {
          householdDrawsRural = rmultinom(n=nDraws, size=totalHouseholdsRural[i], prob=rep(1/totalEAsRural[i], totalEAsRural[i])) + 25
        } else {
          householdDrawsRural = matrix(rep(25, totalEAsRural[i]*nDraws), ncol=nDraws)
        }
      }
    } else {
      if(any(totalHouseholdsRural != 0)) {
        householdDraws = rmultinom(n=nDraws, size=totalHouseholds[i], prob=rep(1/totalEAs[i], totalEAs[i])) + 25
      } else {
        householdDraws = matrix(rep(25, totalEAs[i]*nDraws), ncol=nDraws)
      }
    }
    
    # drawing children per EA, with probability proportional to the number of households
    if(includeUrban) {
      
      if(totalEAsUrban[i] != 0) {
        probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
        # childrenDrawsUrban[areaMatUrban==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
        childrenDraws[areaVals==thisArea & urbanVals] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenUrban[i], probsUrban[,j])})
      }
      
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        # childrenDrawsRural[areaMatRural==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
        childrenDraws[areaVals==thisArea & !urbanVals] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenRural[i], probsRural[,j])})
      }
      
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      childrenDraws[areaVals==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildren[i], probs[,j])})
    }
  }
  
  ##### Return results
  childrenDraws
}


##### Simulation study functions

# this function generates results for the simulation study for the SPDE model
# input arguments:
#   argument specifying the dataset type
resultsSPDE_LCPB = function(randomSeeds=NULL, gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), 
                            effRange=400, beta0=-2.9, surveyI=1, 
                            maxDataSets=NULL, seed=surveyI, representativeSampling=TRUE) {
  # make strings representing the simulation parameters
  dataID = paste0("Beta", round(beta0, 4), "rho", round(rho, 4), "sigmaEps", 
                  round(sigmaEpsilon, 4), "gamma", round(gamma, 4))
  out = load(paste0("savedOutput/simDataSets/simDataMulti", dataID, ".RData"))
  
  if(representativeSampling) {
    clustDat = SRSDat$clustDat
    eaDat = SRSDat$eaDat
  } else {
    clustDat = overSampDat$clustDat
    eaDat = overSampDat$eaDat
  }
  
  if(is.null(maxDataSets)) {
    maxDataSets = length(clustDat)
  }
  
  if(surveyI>maxDataSets) {
    return(NULL)
  }
  
  # get this data set and population frame
  thisData = clustDat[[surveyI]]
  thiseaspa = makeEASPAFromEADat(eaDat)
  
  # generate results for the specified data sets and return results (TODO: otherVariableNames)
  timeSPDE = system.time(resultsSPDE <- fitSPDEKenyaDat(dat=thisData, dataType=c("mort", "ed"), 
                                                       significanceCI=.8, 
                                                       nPostSamples=1000, verbose=TRUE, seed=seed, 
                                                       urbanEffect=TRUE, clusterEffect=TRUE, 
                                                       leaveOutRegionName=NULL, doValidation=FALSE))[3]
  # resultsSPDE = fitModelToDataSets(fitSPDE, dataSets, randomSeeds=randomSeeds, otherArgs=list(mesh=mesh), 
  #                                  maxDataSets=maxDataSets)
  
  # aggregate predictions of the SPDE model
  # browser()
  timeAllAgg = system.time(agg <- modLCPB(uDraws=resultsSPDE$uDraws, resultsSPDE$sigmaEpsilonDraws, easpa=thiseaspa, 
                                       includeUrban=TRUE, clusterLevel=FALSE, pixelLevel=TRUE, constituencyLevel=TRUE, countyLevel=TRUE, 
                                       regionLevel=TRUE, nationalLevel=TRUE, doModifiedPixelLevel=FALSE, 
                                       onlyDoModifiedPixelLevel=FALSE, 
                                       doLCPb=TRUE, doLCpb=TRUE, doLcpb=TRUE, urbanEffectDraws=resultsSPDE$fixedEffectDraws[2,]))[3]
  
  # get scores?
  
  # save results
  fileName = paste0("savedOutput/simStudyResults/resLCPB_", dataID, "repSamp", representativeSampling, "surveyI", surveyI, "Of", maxDataSets, ".RData")
  save(agg, timeSPDE, timeAllAgg, file=fileName)
  
  agg
}




