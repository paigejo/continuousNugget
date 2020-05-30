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
## NOTE1: For aggregation consistent models: 'areas' should be the finest spatial areas for which number 
##        of urban and rural EAs is known. For coarser areas, further aggregate the draws of each 
##        aggregation function as necessary. The same goes for the pixel grid. For larger 'subareas' 
##        containing an unknown number of EAs, in which case aggregate pixel predictions as need be. 
## NOTE2: Aggregation consistent models:
#           pixels/subareas (unknown EAs): CPBL, cpbl (any others?)
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
modCPBL = function(uDraws, sigmaEpsilonDraws, results, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                   includeUrban=TRUE, maxEmptyFraction=1, clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, 
                   regionLevel=TRUE) {
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
  
  # make matrix of pixel indices mapping matrices of EA values to matrices of pixel values
  pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  
  # determine which EAs are urban if necessary
  if(includeUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  epsc = matrix(rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
  
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
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  uc = matrix(uDraws[cbind(rep(rep(1:nrow(uDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
  muc = expit(uc + epsc)
  
  # calculate Z_{ic} for each EA in each pixel
  Zc = matrix(rbinom(n=totalEAs * nDraws, size=Ncs, prob=muc), ncol=nDraws)
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, valMat, applyFun=sum) {
    out = tapply(valMat[,i], factor(as.character(pixelIndexMat[,i])), FUN=applyFun)
    indices = as.numeric(names(out))
    
    vals = rep(NA, nrow(uDraws))
    vals[indices] = out
    vals
  }
  
  ##### Line 5: We already did this, resulting in uDraws input
  
  ##### Line 6: aggregate population denominators for each grid cell to get N_{ig}
  Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, valMat=Ncs)
  Ng[is.na(Ng)] = 0
  
  ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, valMat=Zc)
  
  ##### Line 8: Calculate empirical mortality proportions for each grid cell, p_{ig}. 
  #####         Whenever Z_{ig} is 0, set p_{ig} to 0 as well
  pg = Zg / Ng
  pg[Zg == 0] = 0
  
  # just for testing purposes:
  if(FALSE) {
    mug = sapply(1:ncol(muc), getPixelColumnFromEAs, valMat=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
    mug[!is.finite(mug)] = NA
    
    # this takes a veeeeery long time. Just use the logistic approximation instead
    cpbl = matrix(logitNormMean(cbind(logit(c(as.matrix(uDraws))), rep(sigmaEpsilonDraws, each=nrow(uDraws)))), nrow=nrow(uDraws))
    
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
    png(paste0(figDirectory, "exploratoryAnalysis/CPBLpgMugMean.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    meanRange = quantile(probs=c(.005, .995), c(rowMeans(mug, na.rm=TRUE), rowMeans(pg, na.rm=TRUE), rowMeans(cpbl, na.rm=TRUE)), na.rm=TRUE)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(cpbl, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (cpbl model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(mug, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (CpbL model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, logit(rowMeans(pg, na.rm=TRUE)), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=range(logit(meanRange)), nx=160, ny=160, main=TeX("Posterior mean (CPBL model)"), cex.main=3, col=meanCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    meanTicks = pretty(meanRange, n=5)[-1]
    meanTickLabels = as.character(meanTicks)
    image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
    
    # plot SDs
    SDsUg = apply(cpbl, 1, sd, na.rm=TRUE)
    SDsMug = apply(mug, 1, sd, na.rm=TRUE)
    SDsPg = apply(pg, 1, sd, na.rm=TRUE)
    allSDs = c(SDsUg, SDsMug, SDsPg)
    sdRange = quantile(probs=c(.001, .995), allSDs, na.rm=TRUE)
    png(paste0(figDirectory, "exploratoryAnalysis/CPBLpgMugSD.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    quilt.plot(popMat$lon, popMat$lat, log(SDsUg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (cpbl model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsMug), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (CpbL model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsPg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (CPBL model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    sdRange = zlim
    sdTicks = pretty(sdRange, n=5)[-1]
    sdTickLabels = as.character(sdTicks)
    image.plot(zlim=range(log(sdRange)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
    dev.off()
    
    # plot Ns, pop density
    png(paste0(figDirectory, "exploratoryAnalysis/CPBLpgMugN.png"), width=1500, height=600)
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
    ## cpbl
    # aggregatedResultscpbl = aggregatePixelPredictions(cpbl, Ng, popGrid=popMat, useDensity=TRUE, countyLevel=countyLevel, 
    #                                               regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
    cpblc = matrix(logitNormMean(cbind(c(logit(as.matrix(uc))), rep(sigmaEpsilonDraws, each=nrow(uc)))), nrow=nrow(uc))
    # this takes quite a while (15 mins?)
    aggregatedResultscpbl = aggregateEAPredictions(cpblc, Ncs, areaMat, urbanMat, easpa=easpa, countyLevel=countyLevel, 
                                                   regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
    aggregatedResultscpbl2 = aggregatePixelPredictions(cpbl, Ng2, popGrid=popMat, useDensity=TRUE, countyLevel=countyLevel, 
                                                       regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
    
    # county level
    Pcountycpbl = aggregatedResultscpbl$countyResults$Z
    PcountyUrbancpbl = aggregatedResultscpbl$countyResults$ZUrban
    PcountyRuralcpbl = aggregatedResultscpbl$countyResults$ZRural
    
    Pcountycpbl2 = aggregatedResultscpbl2$countyMatrices$p
    PcountyUrbancpbl2 = aggregatedResultscpbl$countyResults$pUrban
    PcountyRuralcpbl2 = aggregatedResultscpbl$countyResults$pRural
    
    # province level
    Pregioncpbl = aggregatedResultscpbl$regionResults$Z
    PregionUrbancpbl = aggregatedResultscpbl$regionResults$ZUrban
    PregionRuralcpbl = aggregatedResultscpbl$regionResults$ZRural
    
    # national level
    Pnationcpbl = aggregatedResultscpbl$nationalResults$Z
    PnationUrbancpbl = aggregatedResultscpbl$nationalResults$ZUrban
    PnationRuralcpbl = aggregatedResultscpbl$nationalResults$ZRural
    
    ## CpbL
    # this takes quite a while (15 mins?)
    aggregatedResultsCpbL = aggregateEAPredictions(cpblc, Ncs, areaMat, urbanMat, easpa=easpa, countyLevel=countyLevel, 
                                                   regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
    aggregatedResultsCpbL = aggregatePixelPredictions(mug, Ng2, popGrid=popMat, useDensity=FALSE, countyLevel=countyLevel, 
                                                      regionLevel=regionLevel, separateUrbanRural=TRUE, normalize=TRUE)
    
    # county level
    PcountyCpbL = aggregatedResultsCpbL$countyMatrices$p
    PcountyUrbanCpbL = aggregatedResultsCpbL$countyResults$pUrban
    PcountyRuralCpbL = aggregatedResultsCpbL$countyResults$pRural
    
    # province level
    PregionCpbL = aggregatedResultsCpbL$regionResults$p
    PregionUrbanCpbL = aggregatedResultsCpbL$regionResults$pUrban
    PregionRuralCpbL = aggregatedResultsCpbL$regionResults$pRural
    
    # national level
    PnationCpbL = aggregatedResultsCpbL$nationalResults$Z
    PnationUrbanCpbL = aggregatedResultsCpbL$nationalResults$ZUrban
    PnationRuralCpbL = aggregatedResultsCpbL$nationalResults$ZRural
    
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
    meanRange = quantile(probs=c(.005, .995), c(rowMeans(Pcountycpbl2, na.rm=TRUE), rowMeans(PcountyCpbL, na.rm=TRUE), rowMeans(Pcounty, na.rm=TRUE)), na.rm=TRUE)
    meanTicks = pretty(meanRange, n=5)[-1]
    meanTickLabels = as.character(meanTicks)
    
    png(paste0(figDirectory, "exploratoryAnalysis/CPBLpgMugMeanCounty.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    
    plotMapDat(plotVar=rowMeans(Pcountycpbl, na.rm=TRUE), new = TRUE, 
               main="Posterior mean (cpbl model)", scaleFun=logit, scaleFunInverse=expit, 
               cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
               xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
               legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
               plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    plotMapDat(plotVar=rowMeans(PcountyCpbL, na.rm=TRUE), new = TRUE, 
               main="Posterior mean (CpbL model)", scaleFun=logit, scaleFunInverse=expit, 
               cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
               xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = FALSE, 
               legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
               plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    plotMapDat(plotVar=rowMeans(Pcounty, na.rm=TRUE), new = TRUE, 
               main="Posterior mean (CPBL model)", scaleFun=logit, scaleFunInverse=expit, 
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
    SDsUg = apply(cpbl, 1, sd, na.rm=TRUE)
    SDsMug = apply(mug, 1, sd, na.rm=TRUE)
    SDsPg = apply(pg, 1, sd, na.rm=TRUE)
    allSDs = c(SDsUg, SDsMug, SDsPg)
    sdRange = quantile(probs=c(.001, .995), allSDs, na.rm=TRUE)
    png(paste0(figDirectory, "exploratoryAnalysis/CPBLpgMugSD.png"), width=1500, height=600)
    par(mfrow=c(1,3), oma=c( 0,0,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    quilt.plot(popMat$lon, popMat$lat, log(SDsUg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (cpbl model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsMug), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (CpbL model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    quilt.plot(popMat$lon, popMat$lat, log(SDsPg), FUN=function(x){mean(x, na.rm=TRUE)}, 
               zlim=log(sdRange), nx=160, ny=160, main=TeX("Posterior SD (CPBL model)"), cex.main=3, 
               col=sdCols, add.legend=FALSE)
    plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
    plotMapDat(mapDat=regionMap, lwd=2.5)
    
    sdRange = zlim
    sdTicks = pretty(sdRange, n=5)[-1]
    sdTickLabels = as.character(sdTicks)
    image.plot(zlim=range(log(sdRange)), nlevel=length(sdCols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
               legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
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
# eaPixelSamples: this is eaSamples from modCPBL (nPixel x n matrix of number of EAs per pixel)
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
      householdDraws[areaMat==thisArea] = sapply(totalHouseholds[i,], rmultinom, n=1, prob=rep(1/totalEAs[i], totalEAs[i]))
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







