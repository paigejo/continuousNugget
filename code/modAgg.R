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
makeDefaultEASPA = function(dataType=c("children", "women")) {
  dataType = match.arg(dataType)
  
  out = merge(easpc, adjustPopulationPerCountyTable(dataType))
  names(out)[1] = "area"
  
  # sort by area
  sortI = sort(out$area, index.return=TRUE)$ix
  out = out[sortI,]
  
  out
}

# in this model, we assume a binomial process for the EA locations. Follows algorithm 2 from the 
# outline
modLCPB = function(uDraws, sigmaEpsilonDraws, results, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                   includeUrban=TRUE, maxEmptyFraction=1, clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, 
                   regionLevel=TRUE, doModifiedPixelLevel=TRUE) {
  nDraws = ncol(uDraws)
  
  # set default inputs
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA()
    # Area: the name or id of the area
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
  eaSamples = rStratifiedMultnomial(nDraws, popMat, easpa, includeUrban)
  if(doModifiedPixelLevel) {
    eaSamplesMod = rStratifiedBinomial1(nDraws, popMat, easpa, includeUrban)
  }
  
  # make matrix of pixel indices mapping matrices of EA values to matrices of pixel values
  pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  if(doModifiedPixelLevel) {
    pixelIndexListMod = lapply(1:nDraws, function(j) {rep(1:nrow(popMat), times=eaSamplesMod[,j])})
  }
  
  # determine which EAs are urban if necessary
  if(includeUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
    if(doModifiedPixelLevel) {
      urbanListMod = lapply(1:nDraws, function(j) {popMat$urban[pixelIndexListMod[[j]]]})
    }
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  if(doModifiedPixelLevel) {
    areaListMod = lapply(1:nDraws, function(j) {popMat$area[pixelIndexListMod[[j]]]})
  }
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  epsc = matrix(rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
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
  Ncs = sampleNPoissonMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, 
                                  popMat=popMat, includeUrban=includeUrban, verbose=TRUE)
  if(doModifiedPixelLevel) {
    NcsMod = sampleNPoissonBinomial(pixelIndexListMod=pixelIndexListMod, areaListMod=areaListMod, urbanListMod=urbanListMod, includeUrban=includeUrban)
  }
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  uc = matrix(uDraws[cbind(rep(rep(1:nrow(uDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
  muc = expit(uc + epsc)
  if(doModifiedPixelLevel) {
    ucMod = lapply(1:nDraws, function(j) {uDraws[rep(1:nrow(uDraws), times=eaSamplesMod[,j]), j]})
    mucMod = lapply(1:nDraws, function(j) {expit(ucMod[[j]] + epscMod[[j]])})
  }
  
  # calculate Z_{ic} for each EA in each pixel
  Zc = matrix(rbinom(n=totalEAs * nDraws, size=Ncs, prob=muc), ncol=nDraws)
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
  Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, valMat=Ncs)
  Ng[is.na(Ng)] = 0
  
  if(doModifiedPixelLevel) {
    NgMod <- sapply(1:nDraws, getPixelColumnFromEAs, valMat=NcsMod, doMod=TRUE)
    NgMod[is.na(NgMod)] = 0
  }
  
  ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, valMat=Zc)
  
  if(doModifiedPixelLevel) {
    ZgMod <- sapply(1:ncol(Zc), getPixelColumnFromEAs, valMat=Zc, doMod=TRUE)
  }
  
  ##### Line 8: Calculate empirical mortality proportions for each grid cell, p_{ig}. 
  #####         Whenever Z_{ig} is 0, set p_{ig} to 0 as well
  pg = Zg / Ng
  pg[Zg == 0] = 0
  
  if(doModifiedPixelLevel) {
    pgMod = ZgMod / NgMod
    pgMod[is.na(ZgMod)] = 0
  }
  
  # just for testing purposes:
  if(FALSE) {
    mug = sapply(1:ncol(muc), getPixelColumnFromEAs, valMat=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
    mug[!is.finite(mug)] = NA
    
    if(doModifiedPixelLevel) {
      mugMod = sapply(1:ncol(muc), getPixelColumnFromEAs, valMat=mucMod, applyFun=function(x) {mean(x, na.rm=TRUE)}, doMod=TRUE)
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
    
    # plot means
    png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugMean.png"), width=1500, height=600)
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
    
    # plot SDs
    SDsUg = apply(lcpb, 1, sd, na.rm=TRUE)
    SDsMug = apply(mug, 1, sd, na.rm=TRUE)
    SDsPg = apply(pg, 1, sd, na.rm=TRUE)
    allSDs = c(SDsUg, SDsMug, SDsPg)
    sdRange = quantile(probs=c(.001, .995), allSDs, na.rm=TRUE)
    png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugSD.png"), width=1500, height=600)
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
    
    # plot Ns, pop density
    png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugN.png"), width=1500, height=600)
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
    
    ## now do the same but for the modified pixel sampling scheme with >= 1 EA per pixel
    # plot means
    png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugMeanMod.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    # meanRange = quantile(probs=c(.005, .995), c(rowMeans(mugMod, na.rm=TRUE), rowMeans(pgMod, na.rm=TRUE), rowMeans(lcpb, na.rm=TRUE)), na.rm=TRUE)
    meanRange = range(c(rowMeans(mugMod, na.rm=TRUE), rowMeans(pgMod, na.rm=TRUE), rowMeans(lcpb, na.rm=TRUE)), na.rm=TRUE)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(lcpb, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (lcpb model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(mugMod, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (LCpb model, mod.)"), cex.main=3, col=meanCols, add.legend=FALSE)
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
    SDsMugMod = apply(mugMod, 1, sd, na.rm=TRUE)
    SDsPgMod = apply(pgMod, 1, sd, na.rm=TRUE)
    allSDs = c(SDsUg, SDsMugMod, SDsPgMod)
    # sdRange = quantile(probs=c(.001, .995), allSDs, na.rm=TRUE)
    sdRange = range(allSDs, na.rm=TRUE)
    png(paste0(figDirectory, "exploratoryAnalysis/LCPBpgMugSDMod.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    quilt.plot(popMat$lon, popMat$lat, log(SDsUg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (lcpb model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsMugMod), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (LCpb model, mod.)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsPgMod), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (LCPB model, mod.)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    sdTicks = pretty(sdRange, n=5)[-1]
    sdTickLabels = as.character(sdTicks)
    image.plot(zlim=range(log(sdRange)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
  }
  
  ##### Line 9: done aggregating to the pixel level
  
  ##### Line 10: aggregate from pixel level to relevant areal levels
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
    # Area: the name or id of the area
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
sampleNPoissonBinomial = function(pixelIndexListMod, areaListMod, urbanListMod, nDraws=length(pixelIndexListMod), 
                                  easpa=NULL, easpaList=NULL, popMat=NULL, includeUrban=TRUE, verbose=TRUE) {
  
  # set default inputs
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA()
    # Area: the name or id of the area
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
  if(is.null(easpaList)) {
    nEAs = nEAsByStratum(areaListMod, urbanListMod)
    easpaList = lapply(1:nDraws, function(j) {
      thiseaspa = easpa
      thiseaspa[,2:4] = nEAs[[j]][,1:3]
      thiseaspa
    })
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
  
  # subset areaMat by urban/rural if necessary
  if(includeUrban) {
    areaListModUrban = lapply(urbanListMod, function(isUrban) {areaListMod[isUrban]})
    areaListModRural = lapply(urbanListMod, function(isUrban) {areaListMod[!isUrban]})
    urbanStringListMod = lapply(urbanListMod, function(isUrban) {sapply(isUrban, function(singleUrban) {ifelse(singleUrban, "u", "r")})})
    areaUrbanicityListMod = lapply(1:length(areaListMod), function(j) {paste(areaListMod[[j]], urbanStringListMod, sep=",")})
  }
  
  # start by drawing the totals, then divide households amongst EAs, then divide children amongst households. 
  # Make sure there are at least 25 households per EA (only divide the rest randomly)
  
  ##### Draw the total EAs
  
  # get the approximate total number of enumeration areas per stratum (actual total varies between draws)
  areas = easpa$Area
  approxTotalEAsUrban = easpa$EAUrb
  approxTotalEAsRural = easpa$EARur
  approxTotalEAs = easpa$EATotal
  approxnEAs = sum(totalEAs)
  
  ##### draw the number of children per enumeration area
  
  childrenDraws = lapply(1:nDraws, function(j) {rep(NA, sum(nEAs[[j]]$EATotal))})
  
  # Draw the number of households per stratum area that will be randomly distributed (total minus the minimum 25)
  if(includeUrban) {
    # totalHouseholdsUrban = sweep(sapply(easpaList, function(x) {x$HHUrb}), 1, -25*approxTotalEAsUrban, "+")
    # totalHouseholdsRural = sweep(sapply(easpaList, function(x) {x$HHRur}), 1, -25*approxTotalEAsRural, "+")
    totalHouseholdsUrban = sapply(1:nDraws, function(j) {easpaList[[j]]$HHUrb - 25*easpaList[[j]]$EAUrb})
    totalHouseholdsRural = sapply(1:nDraws, function(j) {easpaList[[j]]$HHRur - 25*easpaList[[j]]$EARur})
    totalChildrenUrban = sapply(easpaList, function(x) {x$popUrb})
    totalChildrenRural = sapply(easpaList, function(x) {x$popRur})
    eaUrban = sapply(easpaList, function(x) {x$EAUrb})
    eaRural = sapply(easpaList, function(x) {x$EARur})
  } else {
    totalHouseholds = sapply(1:nDraws, function(j) {easpaList[[j]]$HHTotal - 25*easpaList[[j]]$EATotal})
    totalChildren = sapply(easpaList, function(x) {x$popHHTotal})
    eaTotal = sapply(easpaList, function(x) {x$EATotal})
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
      # householdDrawsUrban = sapply(approxTotalHouseholdsUrban[i,], rbinom, n=approxTotalEAsUrban[i], prob=1/approxTotalEAsUrban[i]) + 25
      householdDrawsUrban = lapply(1:nDraws, function(j) {25 + rbinom(n=eaUrban[i,j], size=totalHouseholdsUrban[i,j], prob=1 / eaUrban[i,j])})
      if(approxTotalEAsRural[i] != 0) {
        householdDrawsRural = lapply(1:nDraws, function(j) {25 + rbinom(n=eaRural[i,j], size=totalHouseholdsRural[i,j], prob=1 / eaRural[i,j])})
      }
    } else {
      # householdDraws = sapply(approxTotalHouseholds[i,], rbinom, n=approxTotalEAs[i], prob=1/approxTotalEAs[i]) + 25
      householdDraws = lapply(1:nDraws, function(j) {25 + rbinom(n=eaTotal[i,j], size=totalHouseholds[i,j], prob=1 / eaTotal[i,j])})
    }
    
    # drawing children per EA, with probability proportional to the number of households
    if(includeUrban) {
      probsUrban = lapply(1:nDraws, function(j) {householdDrawsUrban[[j]] * (1 / sum(householdDrawsUrban[[j]]))})
      thisUrbanDraws = lapply(1:nDraws, function(j) {rbinom(n=eaUrban[i,j], size=totalChildrenUrban[i,j], prob=probsUrban[[j]])})
      childrenDraws = lapply(1:nDraws, function(j) {childrenDraws[[j]][areaListMod[[j]] == thisArea & urbanListMod[[j]]] = thisUrbanDraws})
      
      if(approxTotalEAsRural[i] != 0) {
        probsRural = lapply(1:nDraws, function(j) {householdDrawsRural[[j]] * (1 / sum(householdDrawsRural[[j]]))})
        thisRuralDraws = lapply(1:nDraws, function(j) {rbinom(n=eaRural[i,j], size=totalChildrenRural[i,j], prob=probsRural[[j]])})
        childrenDraws = lapply(1:nDraws, function(j) {childrenDraws[[j]][areaListMod[[j]] == thisArea & !urbanListMod[[j]]] = thisRuralDraws})
      }
    } else {
      probs = lapply(1:nDraws, function(j) {householdDraws[[j]] * (1 / sum(householdDraws[[j]]))})
      thisDraws = lapply(1:nDraws, function(j) {rbinom(n=eaTotal[i,j], size=totalChildren[i,j], prob=probs[[j]])})
      childrenDraws = lapply(1:nDraws, function(j) {childrenDraws[[j]][areaListMod[[j]] == thisArea] = thisDraws})
    }
  }
  
  ##### Return results
  childrenDraws
}







