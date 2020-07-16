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
makeDefaultPopMat = function() {
  out = popGrid
  out$pop = out$popOrig
  out$area = out$admin1
  out$popOrig = NULL
  out$admin1 = NULL
  
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
makeDefaultEASPA = function(dataType=c("children", "women"), validationClusterI=NULL, useClustersAsEAs=!is.null(validationClusterI)) {
  dataType = match.arg(dataType)
  
  if(!useClustersAsEAs) {
    out = merge(easpc, adjustPopulationPerCountyTable(dataType))
  } else {
    if(is.null(validationClusterI)) {
      # This case will likely not be used, and corresponds with leave one county out
      out = merge(easpcMort, poppcMort)
    } else {
      ## construct faux population based on left out observations
      thisMort = mort[validationClusterI]
      
      # faux poppc
      out = aggregate(thisMort$n, by=list(admin1=thisMort$admin1, urban=thisMort$urban), FUN=sum, drop=FALSE)
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
      thisMort = mort[validationClusterI]
      out = aggregate(thisMort$n, by=list(admin1=thisMort$admin1, urban=thisMort$urban), FUN=length, drop=FALSE)
      out[is.na(out[,3]), 3] = 0
      # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
      out2 = cbind(out, rural=0)[48:94,]
      out2[, 4] = out$x[1:47]
      unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
      easpcMort = clustpc
      names(easpcMort) = names(easpc)
      easpcMort[unsortedToSorted, 2:3] = out2[,3:4]
      easpcMort$EATotal = easpcMort$EAUrb + easpcMort$EARur
      
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

# in this model, we assume a binomial process for the EA locations. Follows algorithm 2 from the 
# outline
# validationPixelI: a set of indices of pixels for which we want predictions in validation.
modLCPB = function(uDraws, sigmaEpsilonDraws, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                   includeUrban=TRUE, maxEmptyFraction=1, clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, 
                   regionLevel=TRUE, doModifiedPixelLevel=TRUE, validationPixelI=NULL, 
                   onlyDoModifiedPixelLevel=!is.null(validationPixelI)) {
  nDraws = ncol(uDraws)
  
  # set default inputs
  if(is.null(easpa)) {
    if(!is.null(validationPixelI))
      stop("Must include validationClusterI or easpa in this case")
    
    easpa = makeDefaultEASPA(useClustersAsEAs=!is.null(validationPixelI))
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
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  ##### Line 1: take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  pixelProbs = popMat$pop
  
  # take draws from the stratified binomial process for each posterior sample
  if(!onlyDoModifiedPixelLevel) {
    eaSamples = rStratifiedMultnomial(nDraws, popMat, easpa, includeUrban)
  }
  if(doModifiedPixelLevel) {
    eaSamplesMod = rStratifiedBinomial1(nDraws, popMat, easpa, includeUrban, validationPixelI=validationPixelI)
  }
  
  # make matrix of pixel indices mapping matrices of EA values to matrices of pixel values
  if(!onlyDoModifiedPixelLevel) {
    pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  }
  if(doModifiedPixelLevel) {
    if(is.null(validationPixelI)) {
      pixelIndexListMod = lapply(1:nDraws, function(j) {rep(1:nrow(popMat), times=eaSamplesMod[,j])})
    } else {
      pixelIndexListMod = lapply(1:nDraws, function(j) {rep(sort(validationPixelI), times=eaSamplesMod[,j])})
    }
  }
  
  # determine which EAs are urban if necessary
  if(includeUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    if(!onlyDoModifiedPixelLevel) {
      urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
    }
    if(doModifiedPixelLevel) {
      urbanListMod = lapply(1:nDraws, function(j) {popMat$urban[pixelIndexListMod[[j]]]})
    }
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  if(!onlyDoModifiedPixelLevel) {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  }
  if(doModifiedPixelLevel) {
    areaListMod = lapply(1:nDraws, function(j) {popMat$area[pixelIndexListMod[[j]]]})
  }
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  if(!onlyDoModifiedPixelLevel) {
    epsc = matrix(rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
  }
  if(doModifiedPixelLevel) {
    epscMod = lapply(1:nDraws, function(j) {rnorm(length(areaListMod[[j]]), sd=sigmaEpsilonDraws[j])})
  }
  
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
  if(!onlyDoModifiedPixelLevel) {
    Ncs = sampleNPoissonMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, 
                                    popMat=popMat, includeUrban=includeUrban, verbose=TRUE)
  }
  if(doModifiedPixelLevel) {
    NcsMod = sampleNPoissonBinomial(eaSamplesMod=eaSamplesMod, pixelIndexListMod=pixelIndexListMod, areaListMod=areaListMod, urbanListMod=urbanListMod, includeUrban=includeUrban)
  }
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  if(!onlyDoModifiedPixelLevel) {
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
  
  # calculate Z_{ic} for each EA in each pixel
  if(!onlyDoModifiedPixelLevel) {
    Zc = matrix(rbinom(n=totalEAs * nDraws, size=Ncs, prob=muc), ncol=nDraws)
  }
  if(doModifiedPixelLevel) {
    ZcMod = lapply(1:nDraws, function(j) {rbinom(n=length(mucMod[[j]]), size=NcsMod[[j]], prob=mucMod[[j]])})
  }
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, vals, applyFun=sum, doMod=FALSE) {
    if(!doMod) {
      out = tapply(vals[,i], factor(as.character(pixelIndexMat[,i])), FUN=applyFun)
    } else {
      # in this case
      out = tapply(vals[[i]], factor(as.character(pixelIndexListMod[[i]])), FUN=applyFun)
    }
    
    indices = as.numeric(names(out))
    
    returnValues = rep(NA, nrow(uDraws))
    returnValues[indices] = out
    returnValues
  }
  
  ##### Line 5: We already did this, resulting in uDraws input
  
  ##### Line 6: aggregate population denominators for each grid cell to get N_{ig}
  if(!onlyDoModifiedPixelLevel) {
    Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, vals=Ncs)
    Ng[is.na(Ng)] = 0
  }
  
  if(doModifiedPixelLevel) {
    NgMod <- sapply(1:nDraws, getPixelColumnFromEAs, vals=NcsMod, doMod=TRUE)
    NgMod[is.na(NgMod)] = 0
  }
  
  ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  if(!onlyDoModifiedPixelLevel) {
    Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=Zc)
  }
  
  if(doModifiedPixelLevel) {
    ZgMod <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=ZcMod, doMod=TRUE)
  }
  
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
  }
  
  ##### Line 9: done aggregating to the pixel level
  
  ##### Line 10: aggregate from pixel level to relevant areal levels if necessary
  if(!onlyDoModifiedPixelLevel) {
    aggregatedResults = aggregatePixelPredictions(Zg, Ng, popGrid=popMat, useDensity=FALSE, countyLevel=countyLevel, 
                                                  regionLevel=regionLevel, separateUrbanRural=TRUE)
    
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
    
    browser()
    
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
      aggregatedResultslcpb2 = aggregatePixelPredictions(lcpb*eaSamples, eaSamples, popGrid=popMat, useDensity=FALSE, countyLevel=countyLevel, 
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
      aggregatedResultsLCpb = aggregatePixelPredictions(mug*eaSamples, eaSamples, popGrid=popMat, useDensity=FALSE, countyLevel=countyLevel, 
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
    
    ##### Extra steps: collect draws at each level and generate:
    ##### areas, preds, 
    allMatrices = c(pixelMatrices=list(p=pg, Z=Zg, N=Ng), aggregatedResults)
  } else {
    # for the pixel level validation, return only the pixel results
    allMatrices = c(pixelMatrices=list(p=pgMod, Z=ZgMod, N=NgMod))
  }
  
  allMatrices
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
                                             doPixelLevelValidation=TRUE, doCountyLevelValidation=FALSE, 
                                             loadPreviousResults=FALSE, 
                                             doLCPB=TRUE, doLCPb=FALSE, doLCpb=FALSE, doLcpb=FALSE, dolcpb=FALSE) {
  
  if(doLCPb || doLCpb || doLcpb || dolcpb) {
    stop("Only LCPB aggregation model is currently supported in validation")
  }
  
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
    print("Fitting full point level model")
    timePoint = system.time(fit <- fitSPDEKenyaDat(dat, dataType, mesh, prior, significanceCI, int.strategy, strategy, nPostSamples,
                                                   verbose, link, NULL, urbanEffect, clusterEffect=TRUE,
                                                   kmres=kmres, doValidation=TRUE, family="binomial"))[3]
    browser()
    print("Fitting full population aggregation model")
    truthTableFull = getPixelLevelTruth(dat=dat, popMat=NULL, targetPop="children")
    easpaMort = makeDefaultEASPA(dataType2, NULL, useClustersAsEAs=TRUE)
    timeAggregation = system.time(fit2 <- modLCPB(uDraws=fit$uDraws, fit$sigmaEpsilonDraws, easpa=easpaMort, popMat=NULL, empiricalDistributions=NULL, 
                                                  includeUrban=urbanEffect, clusterLevel=FALSE, pixelLevel=TRUE, countyLevel=FALSE, 
                                                  regionLevel=FALSE, doModifiedPixelLevel=TRUE, validationPixelI=truthTableFull$pixelI, 
                                                  onlyDoModifiedPixelLevel=TRUE))[3]
    
    # get observations and prediction summary statistics
    # easpaFull = makeDefaultEASPA(validationClusterI=NULL, useClustersAsEAs=TRUE)
    truthFull = truthTableFull$p
    truthUrban = truthTableFull$urban
    est = rowMeans(fit2$pixelMatrices$p)
    vars = apply(fit2$pixelMatrices$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = fit2$pixelMatrices$p
    # estMatBinomial = addBinomialVar(estMat, dat$n)
    
    cpo = fit$mod$cpo$cpo
    cpoFailure = fit$mod$cpo$failure
    dic = fit$mod$dic$dic
    waic = fit$mod$waic$waic
    modelFit = fit$mod
    
    # # calculate validation scoring rules
    # print("Pooled scores:")
    # fullPooledScoresBinomial = data.frame(c(getScores(truth, est, vars, lower, upper, estMatBinomial, doRandomReject=TRUE), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), Time=time[3]))
    # print(fullPooledScoresBinomial)
    # print("Rural scores:")
    # fullRuralScoresBinomial = data.frame(c(getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMatBinomial[!obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), Time=time[3]))
    # print(fullRuralScoresBinomial)
    # print("Urban scores:")
    # fullUrbanScoresBinomial = data.frame(c(getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMatBinomial[obsUrban,], doRandomReject=TRUE), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), Time=time[3]))
    # print(fullUrbanScoresBinomial)
    # 
    # fullPooledScores = data.frame(c(getScores(truth, est, vars, lower, upper, estMat), WAIC=waic, DIC=dic, CPO=mean(cpo, na.rm=TRUE), Time=time[3]))
    # fullRuralScores = data.frame(c(getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,]), WAIC=NA, DIC=NA, CPO=mean(cpo[!obsUrban], na.rm=TRUE), Time=time[3]))
    # fullUrbanScores = data.frame(c(getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,]), WAIC=NA, DIC=NA, CPO=mean(cpo[obsUrban], na.rm=TRUE), Time=time[3]))
    
    if(saveResults) {
      save(timePoint=timePoint, timeAggregation=timeAggregation, fit=fit, fit2=fit2, file=fileName)}
  }
  else {
    print("Loading previous full model fit")
    load(fileName)
  }
  previousFit = fit
  
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
    obsUrban = popMat$east[inSamplePixelIndices]
    
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
  
  completeScoreTableBinomial = c()
  pooledScoreTableBinomial = c()
  urbanScoreTableBinomial = c()
  ruralScoreTableBinomial = c()
  completeScoreTable = c()
  pooledScoreTable = c()
  urbanScoreTable = c()
  ruralScoreTable = c()
  binnedScoringRulesuuAll = list()
  binnedScoringRulesuUAll = list()
  binnedScoringRulesUuAll = list()
  binnedScoringRulesUUAll = list()
  binnedScoringRulesAuAll = list()
  binnedScoringRulesAUAll = list()
  binnedScoringRulesuAAll = list()
  binnedScoringRulesUAAll = list()
  binnedScoringRulesAAAll = list()
  singleScores = c()
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
    
    # fit the point level model
    timePoint = system.time(fit <- fitSPDEKenyaDat(dat, dataType, mesh, prior, significanceCI, int.strategy, strategy, nPostSamples, 
                                                   verbose, link, NULL, urbanEffect, clusterEffect, thisRegion,
                                                   kmres, TRUE, previousFit, family="binomial", leaveOutI=leaveOutI))[3]
    
    # get the aggregation model predictions
    aggregatedPreds = modLCPB(fit$uDraws, fit$sigmaEpsilonDraws, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                              includeUrban=urbanEffect, maxEmptyFraction=1, clusterLevel=FALSE, pixelLevel=TRUE, countyLevel=FALSE, 
                              regionLevel=FALSE, doModifiedPixelLevel=TRUE, validationPixelI=)
    
    # get observations and prediction summary statistics
    truth = (dat$y / dat$n)[thisSampleI]
    obsUrban = dat$urban[thisSampleI]
    est = fit$preds
    vars = fit$sigmas^2
    # lower = fit$lower
    # upper = fit$upper
    lower = NULL
    upper = NULL
    estMat = fit$predMat
    
    thisTruthTable = truthTable[match(thisPixelI, truthTable$pixelI),]
    thisTruth = thisTruthTable$p
    truthUrban = thisTruthTable$urban
    est = rowMeans(aggregatedPreds$pixelMatrices$p)
    vars = apply(aggregatedPreds$pixelMatrices$p, 1, var)
    # lower = fit$obsLower
    # upper = fit$obsUpper
    lower = NULL
    upper = NULL
    estMat = aggregatedPreds$pixelMatrices$p
    
    # calculate validation scoring rules
    print("Pooled scores:")
    if(!stratifiedValidation)
      thisPooledScores = data.frame(c(list(Region=thisRegion), getScores(truth, est, vars, lower, upper, estMat, doRandomReject=TRUE), Time=time[3]))
    else
      thisPooledScores = data.frame(c(list(Fold=i), getScores(truth, est, vars, lower, upper, estMat, doRandomReject=TRUE), Time=time[3]))
    print(thisPooledScores)
    
    if(stratifiedValidation || thisRegion != "Nairobi") {
      print("Rural scores:")
      if(!stratifiedValidation)
        thisRuralScores = data.frame(c(list(Region=thisRegion), getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), Time=time[3]))
      else
        thisRuralScores = data.frame(c(list(Fold=i), getScores(truth[!obsUrban], est[!obsUrban], vars[!obsUrban], lower[!obsUrban], upper[!obsUrban], estMat[!obsUrban,], doRandomReject=TRUE), Time=time[3]))
      print(thisRuralScores)
    } else {
      thisRuralScores = thisPooledScores
      thisRuralScores[,2:(ncol(thisRuralScores)-1)] = NA
      thisRuralScores = thisRuralScores
    }
    
    print("Urban scores:")
    if(!stratifiedValidation)
      thisUrbanScores = data.frame(c(list(Region=thisRegion), getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), Time=time[3]))
    else
      thisUrbanScores = data.frame(c(list(Fold=i), getScores(truth[obsUrban], est[obsUrban], vars[obsUrban], lower[obsUrban], upper[obsUrban], estMat[obsUrban,], doRandomReject=TRUE), Time=time[3]))
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
    obsUrban = popMat$east[inSamplePixelIndices]
    
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
      
      binnedScoringRulesuu = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsuu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesUu = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsUu, breaks=distanceBreaks)$binnedResults
      binnedScoringRulesAu = getScores(truth[!predUrban], est[!predUrban], vars[!predUrban], estMat=estMat[!predUrban,], doRandomReject=TRUE, distances=nndistsAu, breaks=distanceBreaks)$binnedResults
    } else {
      binnedScoringRulesuu = NULL
      binnedScoringRulesUu = NULL
      binnedScoringRulesAu = NULL
    }
    
    # calculate scores accounting for binomial variation using fuzzy reject intervals
    binnedScoringRulesuU = getScores(truth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsuU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUU = getScores(truth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsUU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAU = getScores(truth[predUrban], est[predUrban], vars[predUrban], estMat=estMat[predUrban,], doRandomReject=TRUE, distances=nndistsAU, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesuA = getScores(truth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsuA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesUA = getScores(truth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsUA, breaks=distanceBreaks)$binnedResults
    binnedScoringRulesAA = getScores(truth, est, vars, estMat=estMat, doRandomReject=TRUE, distances=nndistsAA, breaks=distanceBreaks)$binnedResults
    
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
      thisSingleScores = data.frame(c(list(Region=thisRegion, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMat, getAverage=FALSE), Time=time[3]))
    }
    else {
      thisSingleScores = data.frame(c(list(Fold=i, dataI=which(thisSampleI), NNDist=nndistsAA, NNDistU=nndistsUA, NNDistu=nndistsuA), getScores(truth, est, vars, lower, upper, estMat, getAverage=FALSE), Time=time[3]))
    }
    
    # concatenate the results
    singleScores = rbind(singleScores, thisSingleScores)
    
    # save results so far
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
                                     popMat=NULL, includeUrban=TRUE, verbose=TRUE) {
  
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
  if(includeUrban) {
    areaMatUrban = matrix(areaMat[urbanMat], ncol=nDraws)
    areaMatRural = matrix(areaMat[!urbanMat], ncol=nDraws)
    urbanStringMat = matrix(sapply(urbanMat, function(x) {ifelse(x, "u", "r")}), ncol=nDraws)
    areaUrbanicityMat = matrix(paste(areaMat, urbanStringMat, sep=","), ncol=nDraws)
  }
  
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
      householdDrawsUrban = sapply(totalHouseholdsUrban[i,], rmultinom, n=1, prob=rep(1/totalEAsUrban[i], totalEAsUrban[i])) + 25
      if(totalEAsRural[i] != 0) {
        householdDrawsRural = sapply(totalHouseholdsRural[i,], rmultinom, n=1, prob=rep(1/totalEAsRural[i], totalEAsRural[i])) + 25
      }
    } else {
      householdDraws = sapply(totalHouseholds[i,], rmultinom, n=1, prob=rep(1/totalEAs[i], totalEAs[i]))
    }
    
    # drawing children per EA, with probability proportional to the number of households
    if(includeUrban) {
      probsUrban = sweep(householdDrawsUrban, 2, 1 / colSums(householdDrawsUrban), "*")
      # childrenDrawsUrban[areaMatUrban==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
      childrenDraws[areaMat==thisArea & urbanMat] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenUrban[i,j], probsUrban[,j])})
      if(totalEAsRural[i] != 0) {
        probsRural = sweep(householdDrawsRural, 2, 1 / colSums(householdDrawsRural), "*")
        # childrenDrawsRural[areaMatRural==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
        childrenDraws[areaMat==thisArea & !urbanMat] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildrenRural[i,j], probsRural[,j])})
      }
    } else {
      probs = sweep(householdDraws, 2, 1 / colSums(householdDraws), "*")
      childrenDraws[areaMat==thisArea] = sapply(1:nDraws, function(j) {rmultinom(1, totalChildren[i,j], probs[,j])})
    }
  }
  
  ##### Return results
  childrenDraws
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
                                  easpa=NULL, easpaList=NULL, popMat=NULL, includeUrban=TRUE, verbose=TRUE) {
  
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







