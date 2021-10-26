# this script is for comparing performance between models
doSetup = TRUE

# download data ----
# download 2014 Kenya population density and associated TIF file
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/pop.rda?raw=true")
popFilename = paste0(tempDirectory, "/pop.rda")
if(!file.exists(popFilename)) {
  download.file(githubURL,popFilename)
}

githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/worldpop_total_1y_2014_00_00.tif?raw=true")
popTIFFilename = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")
if(!file.exists(popTIFFilename)) {
  download.file(githubURL,popTIFFilename)
}

# load it in
require(raster)
out = load(popFilename)
out

# make sure this is correct for re-projections
pop@file@name = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")

# get frame/pop info ----
popMatSimple = makePopIntegrationTab(kmRes=5, pop=pop, domainPoly=kenyaPoly, 
                                     eastLim=eastLim, northLim=northLim, 
                                     mapProjection=projKenya, poppa=poppaKenya, 
                                     poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                     areaMapDat=adm1, subareaMapDat=adm2, 
                                     areaPolygonSubsetI=30)

popMatSimpleNeonatal = adjustPopMat(popMatSimple, poppaTarget=poppsubKenyaNeonatal, adjustBy="subarea")
easpaSimple = makeDefaultEASPA()
easpaSimple = easpaSimple[easpaSimple$area == "Nairobi",]
poppsubSimple = poppsubKenya
poppsubSimple = poppsubSimple[poppsubSimple$area == "Nairobi",]

# simulate truths ----
constituenciesN = poppsubKenya$area == "Nairobi"
countyN = sort(unique(poppaKenya$area)) == "Nairobi"
seeds = sample(1:100000, 100, replace=FALSE)
inlaSeeds = sample(1:100000, 100, replace=FALSE)
if(doSetup) {
  truths = list()
  dats = list()
  for(i in 1:100) {
    # simulate population over all of Kenya and generate survey from the EAs
    thisTime = system.time(simDatKenya <- generateSimDataSetsLCPB2(nsim=1, targetPopMat=popMatKenyaNeonatal, 
                                                                   popMat=popMatKenya, 
                                                                   doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
                                                                   gridLevel=FALSE, subareaLevel=TRUE, 
                                                                   fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                                                   logisticApproximation=FALSE, 
                                                                   dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                                                   seed=seeds[i], inla.seed=inlaSeeds[i], simPopOnly=FALSE, returnEAinfo=TRUE, 
                                                                   easpa=easpaKenya, poppsub=poppsubKenya))
    
    simDatKenya$overSampDat = NULL
    simDatKenya$simulatedEAs$eaDat = NULL
    simDatKenya$simulatedEAs$eaSamples = NULL
    simDatKenya$simulatedEAs$clustDat = NULL
    simDatKenya$simulatedEAs$thisclustpc = NULL
    
    # get the data and the true population
    dat = simDatKenya$SRSDat$clustDat[[1]]
    dat$y = dat$Z
    dat$n = dat$N
    dat$admin1 = dat$area
    dat$admin2 = dat$subarea
    dats[[i]] = dat
    
    EAsN = simDatKenya$simulatedEAs$aggregatedPop$eaPop$eaDatList[[1]]$area == "Nariobi"
    truePrevalenceEAsKenya = simDatKenya$simulatedEAs$aggregatedPop$eaPop$eaDatList[[1]]$pFineScalePrevalence
    truePrevalenceEAsKenya = truePrevalenceEAsKenya[EAsN]
    truePrevalenceConstituencyKenya = simDatKenya$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence
    truePrevalenceConstituencyKenya = truePrevalenceConstituencyKenya[constituenciesN]
    truePrevalenceCountyKenya = truePrevalenceCountyKenya[countyN]
    truePrevalenceCountyKenya = simDatKenya$simulatedEAs$aggregatedPop$areaPop$aggregationResults$pFineScalePrevalence
    
    truths[[i]]  = list(EAsN=EAsN, truePrevalenceEAsKenya=truePrevalenceEAsKenya, 
                        truePrevalenceConstituencyKenya=truePrevalenceConstituencyKenya, 
                        truePrevalenceCountyKenya=truePrevalenceCountyKenya)
  }
  save(dats, truths, file="savedOutput/simpleExample/gridResTest_datsAndTruths.RData")
} else {
  out = load("savedOutput/simpleExample/gridResTest_datsAndTruths.RData")
}



# Generate grids ----
# construct integration grids at different resolutions
# resolutions = c(1, 10, 20, 40, 60, 80)
# deltas = c(.05, .2, .4, .8, 1.2, 1.4)
# meanNeighbors = c(500, rep(50, 5))
# resolutions = c(.1, .2, .5, 1, 2, 5, 10, 20)
resolutions = c(.2, 1, 5, 25)
# deltas = c(.005, .01, .025, .05, .1, .1, .2, .4)
deltas = c(.01, .05, .1, .5)
# meanNeighbors = c(rep(500, 5), rep(50, 3))
meanNeighbors = c(rep(500, 2), rep(50, 2))
if(doSetup) {
  popGrids = list()
  popGridsAdjusted = list()
  for(i in 1:length(resolutions)) {
    print(paste0("Creating integration grid at ", resolutions[i], " km resolution"))
    # thisPopGrid = makePopIntegrationTab(kmRes=resolutions[i], mapProjection=SUMMER::projKenya, 
    #                                     domainPoly=kenyaPoly, 
    #                                     subareaMapDat=adm2, poppsub=poppsub, 
    #                                     eastLim=eastLim, northLim=northLim, 
    #                                     mean.neighbor=meanNeighbors[i], 
    #                                     delta=deltas[i], 
    #                                     areaMapDat=adm1, areaPolygonSubsetI=30) # Nairobi is the 30th one
    # thisPopGrid$area = thisPopGrid$area
    # thisPopGrid$constituency = thisPopGrid$subarea
    # thisPopGridAdjusted = adjustPopGrid(thisPopGrid, poppconAdjusted, "Constituency")
    
    thisPopGrid = makePopIntegrationTab(kmRes=resolutions[i], pop=pop, domainPoly=kenyaPoly, 
                                        eastLim=eastLim, northLim=northLim, 
                                        mapProjection=SUMMER::projKenya, poppa=poppaKenya, 
                                        poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                        areaMapDat=adm1, subareaMapDat=adm2, 
                                        mean.neighbor=meanNeighbors[i], 
                                        delta=deltas[i], 
                                        areaPolygonSubsetI=30)
    
    thisPopGridAdjusted = adjustPopMat(thisPopGrid, poppaTarget=poppsubKenyaNeonatal, adjustBy="subarea")
    
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
  
  save(popGrids, popGridsAdjusted, ns, endIs, startIs, file="savedOutput/simpleExample/popGridResTest.RData")
} else {
  out = load("savedOutput/simpleExample/popGridResTest.RData")
}

# combine the grids into 1 matrix for evaluating the SPDE model all at once
popMatCombined = do.call("rbind", popGrids)

# Run main results ----
# fit SPDE cluster level risk model over all integration points, for each truth, 
# and then generate predictions from the aggregation model conditional on the 
# risk model predictions
# spde$effRange, spde$margVar, familyPrec, clusterPrec, beta
fixedParameters = list(spde=list(effRange=400, margVar=(1/3)^2), clusterPrec=2.5, beta=c(-2.9, -1))
nSamples = rep(1000, length(resolutions))
allAggResultsN = list()
startI = 1
endI = length(truths)
for(i in startI:endI) {
  print(paste0("Running analysis for truth ", i, "/", length(truths), ", endI=", endI))
  
  dat = dats[[i]]
  spdeFitN = fitSPDEKenyaDat(dat, nPostSamples=max(nSamples), popMat=popMatCombined, 
                             fixedParameters=fixedParameters, prior=NULL)
  
  # remove unnecessary parts of the spdeFitN object to save space
  spdeFitN$mod$.args = NULL
  spdeFitN$mod$all.hyper = NULL
  
  # apply aggregation models at each resolution
  # nSamples = c(c(500, 1000), rep(10000, length(resolutions)-1))
  separateUDraws = list()
  for(j in 1:length(popGrids)) {
    # obtain the grids at this resolution
    thisNSamples = nSamples[j]
    
    # obtain model output at this resolution
    thisResolutionI = startIs[j]:endIs[j]
    separateUDraws = c(separateUDraws, spdeFitN$uDraws[thisResolutionI,1:thisNSamples])
  }
  spdeFitN$uDraws = NULL
  
  # save spdeFitN and its relevant uDraws
  # save(spdeFitN, separateUDraws, file="savedOutput/simpleExample/spdeFitNandUDraws.RData")
  # out = load("savedOutput/simpleExample/spdeFitNandUDraws.RData")
  
  # run the actual analysis
  aggResultsN = list()
  for(j in 1:length(popGrids)) {
    # obtain the grids at this resolution
    thisNSamples = nSamples[j]
    thisPopMat = popGrids[[j]]
    thisPopMatAdjusted = popGridsAdjusted[[j]]
    print(paste0("Running analysis for grid ", j, "/", length(popGrids), " with ", nrow(popGridsAdjusted[[j]]), " points"))
    startTime = proc.time()[3]
    
    # obtain model output at this resolution
    thisResolutionI = startIs[j]:endIs[j]
    # thisUDraws = spdeFitN$uDraws[thisResolutionI,1:thisNSamples]
    thisUDraws = separateUDraws[[j]]
    sigmaEpsilonDraws = spdeFitN$sigmaEpsilonDraws[1:thisNSamples]
    
    thisAggResultsN = simPopCustom(thisUDraws, sigmaEpsilonDraws, easpaSimple, thisPopMat, 
                                   thisPopMatAdjusted, doFineScaleRisk=TRUE, doGriddedRisk=TRUE, 
                                   doSmoothRisk=TRUE, subareaLevel=TRUE, gridLevel=FALSE, 
                                   poppsub=poppsubSimple, min1PerSubarea=TRUE, 
                                   doSmoothRiskLogisticApprox=FALSE, returnEAinfo=TRUE, 
                                   fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                   verbose=FALSE)
    
    # remove unnecessary components to save space
    thisAggResultsN = list(subareaPop = thisAggResultsN$subareaPop)
    thisAggResultsN$subareaPop$aggregationMatrices = NULL
    
    # compile results
    aggResultsN = c(aggResultsN, list(thisAggResultsN))
    endTime = proc.time()[3]
    print(paste0("Took ", (endTime-startTime)/60, " minutes"))
  }
  names(aggResultsN) = paste0("aggResultsN", resolutions)
  
  # add results for this truth
  allAggResultsN = c(allAggResultsN, list(aggResultsN))
}

save(allAggResultsN, 
     file=paste0("savedOutput/simpleExample/gridResolutionTestNairobi_", 
                 startI, "_", endI, ".RData"))

# concatenate all results if need be:
tempAllAggResultsN = list()
if(FALSE) {
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_1_10.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_11_20.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_21_30.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_31_40.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_41_50.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_51_60.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_61_70.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_71_80.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_81_90.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  out = load(paste0("savedOutput/simpleExample/gridResolutionTestNairobi_91_100.RData"))
  tempAllAggResultsN = c(tempAllAggResultsN, allAggResultsN)
  
  allAggResultsN = tempAllAggResultsN
}

out = load("savedOutput/simpleExample/gridResolutionTestNairobi.RData")

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
allPredsSmoothRisk = list()
allPredsRisk = list()
allPredsPrevalence = list()
allPredsGriddedRisk = list()
allCIWidthsSmoothRisk = list()
allCIWidthsRisk = list()
allCIWidthsPrevalence = list()
allCIWidthsGriddedRisk = list()
allCoveragesSmoothRisk = list()
allCoveragesRisk = list()
allCoveragesPrevalence = list()
allCoveragesGriddedRisk = list()
for(i in 1:length(truths)) {
  # get truth
  truePrevalenceConstituencyKenya = truths[[i]]$truePrevalenceConstituencyKenya
  
  # Calculate RMSE, 80% and 95% Coverage
  predsSmoothRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
  predsRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
  predsPrevalence = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
  predsGriddedRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
  residsConstituencySmoothRisk = list()
  residsConstituencyRisk = list()
  residsConstituencyPrevalence = list()
  residsConstituencyGriddedRisk = list()
  for(j in 1:length(resolutions)) {
    thisMaxSamples = min(c(nSamples[j], maxSamples))
    predsSmoothRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pSmoothRisk[,1:thisMaxSamples])
    predsRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScaleRisk[,1:thisMaxSamples])
    predsPrevalence[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScalePrevalence[,1:thisMaxSamples])
    predsGriddedRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pGriddedRisk[,1:thisMaxSamples])
    theseResidsSmoothRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pSmoothRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
    theseResidsRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScaleRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
    theseResidsPrevalence = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScalePrevalence[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
    theseResidsGriddedRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pGriddedRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
    residsConstituencySmoothRisk = c(residsConstituencySmoothRisk, list(theseResidsSmoothRisk))
    residsConstituencyRisk = c(residsConstituencyRisk, list(theseResidsRisk))
    residsConstituencyPrevalence = c(residsConstituencyPrevalence, list(theseResidsPrevalence))
    residsConstituencyGriddedRisk = c(residsConstituencyGriddedRisk, list(theseResidsGriddedRisk))
  }
  allPredsSmoothRisk = c(allPredsSmoothRisk, list(predsSmoothRisk))
  allPredsRisk = c(allPredsRisk, list(predsRisk))
  allPredsPrevalence = c(allPredsPrevalence, list(predsPrevalence))
  allPredsGriddedRisk = c(allPredsGriddedRisk, list(predsGriddedRisk))
  
  lowConstituencySmoothRisk = lapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
  lowConstituencyRisk = lapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
  lowConstituencyPrevalence = lapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
  lowConstituencyGriddedRisk = lapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
  highConstituencySmoothRisk = lapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
  highConstituencyRisk = lapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
  highConstituencyPrevalence = lapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
  highConstituencyGriddedRisk = lapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
  
  # CIWidthSmoothRisk = highConstituencySmoothRisk - lowConstituencySmoothRisk
  # CIWidthRisk = highConstituencyRisk - lowConstituencyRisk
  # CIWidthPrevalence = highConstituencyPrevalence - lowConstituencyPrevalence
  # CIWidthGriddedRisk = highConstituencyGriddedRisk - lowConstituencyGriddedRisk
  CIWidthSmoothRisk = lapply(residsConstituencySmoothRisk, 
                             function(mat) {
                               apply(mat, 1, function(x) {
                                 quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                               })
  CIWidthRisk = lapply(residsConstituencyRisk, 
                       function(mat) {
                         apply(mat, 1, function(x) {
                           quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                       })
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
  allCIWidthsSmoothRisk = c(allCIWidthsSmoothRisk, list(CIWidthSmoothRisk))
  allCIWidthsRisk = c(allCIWidthsRisk, list(CIWidthRisk))
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
  inCISmoothRisk = lapply(residsConstituencySmoothRisk, 
                          function(mat) {
                            apply(mat, 1, function(x) {
                              (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                          })
  inCIRisk = lapply(residsConstituencyRisk, 
                    function(mat) {
                      apply(mat, 1, function(x) {
                        (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                    })
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
  coverageSmoothRisk = sapply(inCISmoothRisk, rowMeans)
  coverageRisk = sapply(inCIRisk, rowMeans)
  coveragePrevalence = sapply(inCIPrevalence, rowMeans)
  coverageGriddedRisk = sapply(inCIGriddedRisk, rowMeans)
  allCoveragesSmoothRisk = c(allCoveragesSmoothRisk, list(coverageSmoothRisk))
  allCoveragesRisk = c(allCoveragesRisk, list(coverageRisk))
  allCoveragesPrevalence = c(allCoveragesPrevalence, list(coveragePrevalence))
  allCoveragesGriddedRisk = c(allCoveragesGriddedRisk, list(coverageGriddedRisk))
}

# compile relevant results

## central predictions
allPredsSmoothRiskMat = do.call("rbind", allPredsSmoothRisk)
allPredsRiskMat = do.call("rbind", allPredsRisk)
allPredsPrevalenceMat = do.call("rbind", allPredsPrevalence)
allPredsGriddedRiskMat = do.call("rbind", allPredsGriddedRisk)

## CI widths
allCIWidthsSmoothRisk80 = lapply(allCIWidthsSmoothRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[3,]})
})
allCIWidthsSmoothRisk80 = do.call("rbind", allCIWidthsSmoothRisk80)
allCIWidthsSmoothRisk90 = lapply(allCIWidthsSmoothRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[2,]})
})
allCIWidthsSmoothRisk90 = do.call("rbind", allCIWidthsSmoothRisk90)
allCIWidthsSmoothRisk95 = lapply(allCIWidthsSmoothRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[1,]})
})
allCIWidthsSmoothRisk95 = do.call("rbind", allCIWidthsSmoothRisk95)

allCIWidthsRisk80 = lapply(allCIWidthsRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[3,]})
})
allCIWidthsRisk80 = do.call("rbind", allCIWidthsRisk80)
allCIWidthsRisk90 = lapply(allCIWidthsRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[2,]})
})
allCIWidthsRisk90 = do.call("rbind", allCIWidthsRisk90)
allCIWidthsRisk95 = lapply(allCIWidthsRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[1,]})
})
allCIWidthsRisk95 = do.call("rbind", allCIWidthsRisk95)

allCIWidthsPrevalence80 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[3,]})
})
allCIWidthsPrevalence80 = do.call("rbind", allCIWidthsPrevalence80)
allCIWidthsPrevalence90 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[2,]})
})
allCIWidthsPrevalence90 = do.call("rbind", allCIWidthsPrevalence90)
allCIWidthsPrevalence95 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[1,]})
})
allCIWidthsPrevalence95 = do.call("rbind", allCIWidthsPrevalence95)

allCIWidthsGriddedRisk80 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[3,]})
  })
allCIWidthsGriddedRisk80 = do.call("rbind", allCIWidthsGriddedRisk80)
allCIWidthsGriddedRisk90 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[2,]})
})
allCIWidthsGriddedRisk90 = do.call("rbind", allCIWidthsGriddedRisk90)
allCIWidthsGriddedRisk95 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
  sapply(listOfResolutions, function(x) {x[1,]})
})
allCIWidthsGriddedRisk95 = do.call("rbind", allCIWidthsGriddedRisk95)

## coverages
allCoveragesSmoothRisk80 = sapply(allCoveragesSmoothRisk, function(x) {x[3,]})
allCoveragesSmoothRisk90 = sapply(allCoveragesSmoothRisk, function(x) {x[2,]})
allCoveragesSmoothRisk95 = sapply(allCoveragesSmoothRisk, function(x) {x[1,]})
meanCoveragesSmoothRisk80 = rowMeans(allCoveragesSmoothRisk80)
meanCoveragesSmoothRisk90 = rowMeans(allCoveragesSmoothRisk90)
meanCoveragesSmoothRisk95 = rowMeans(allCoveragesSmoothRisk95)

allCoveragesRisk80 = sapply(allCoveragesRisk, function(x) {x[3,]})
allCoveragesRisk90 = sapply(allCoveragesRisk, function(x) {x[2,]})
allCoveragesRisk95 = sapply(allCoveragesRisk, function(x) {x[1,]})
meanCoveragesRisk80 = rowMeans(allCoveragesRisk80)
meanCoveragesRisk90 = rowMeans(allCoveragesRisk90)
meanCoveragesRisk95 = rowMeans(allCoveragesRisk95)

allCoveragesPrevalence80 = sapply(allCoveragesPrevalence, function(x) {x[3,]})
allCoveragesPrevalence90 = sapply(allCoveragesPrevalence, function(x) {x[2,]})
allCoveragesPrevalence95 = sapply(allCoveragesPrevalence, function(x) {x[1,]})
meanCoveragesPrevalence80 = rowMeans(allCoveragesPrevalence80)
meanCoveragesPrevalence90 = rowMeans(allCoveragesPrevalence90)
meanCoveragesPrevalence95 = rowMeans(allCoveragesPrevalence95)

allCoveragesGriddedRisk80 = sapply(allCoveragesGriddedRisk, function(x) {x[3,]})
allCoveragesGriddedRisk90 = sapply(allCoveragesGriddedRisk, function(x) {x[2,]})
allCoveragesGriddedRisk95 = sapply(allCoveragesGriddedRisk, function(x) {x[1,]})
meanCoveragesGriddedRisk80 = rowMeans(allCoveragesGriddedRisk80)
meanCoveragesGriddedRisk90 = rowMeans(allCoveragesGriddedRisk90)
meanCoveragesGriddedRisk95 = rowMeans(allCoveragesGriddedRisk95)

# Plot central predictions versus resolution
tempRes = resolutions[col(allPredsSmoothRiskMat)]
tempCon = factor(as.character(poppsubSimple$subarea[col(allPredsSmoothRiskMat)]))
N=length(tempCon)
preds = c(allPredsSmoothRiskMat, allPredsRiskMat, 
          allPredsPrevalenceMat, allPredsGriddedRiskMat)
predFrame = data.frame(Constituency=rep(tempCon, 4), 
                          Resolution=rep(tempRes, 4), 
                          Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                         rep("Prevalence", N), rep("Gridded risk", N)), 
                                       levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                          preds=preds)
pdf(paste0("figures/gridResolutionTest/predictionVRes.pdf"), width=5, height=5)
ggplot(predFrame, aes(factor(Resolution), preds, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
  labs(x="Grid resolution (km)", y="Central Predictions", fill="Model") + 
  theme_classic()
dev.off()

# Plot CI Widths versus resolution and model
## 80% CIs
CIWidth = c(c(allCIWidthsSmoothRisk80), c(allCIWidthsRisk80), 
            c(allCIWidthsPrevalence80), c(allCIWidthsGriddedRisk80))
tempRes = resolutions[col(allCIWidthsSmoothRisk80)]
tempCon = factor(as.character(poppsubSimple$subarea[col(allCIWidthsSmoothRisk80)]))
N=length(tempCon)
CIWidthFrame = data.frame(Constituency=rep(tempCon, 4), 
                          Resolution=rep(tempRes, 4), 
                          Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                         rep("Prevalence", N), rep("Gridded risk", N)), 
                                       levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                          CIWidth=CIWidth)

pdf(paste0("figures/gridResolutionTest/CIWidthVRes80.pdf"), width=7, height=5)
ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
  labs(x="Grid resolution (km)", y="80% credible interval width", fill="Model") + 
  theme_classic()
dev.off()

## 90% CIs
CIWidth = c(c(allCIWidthsSmoothRisk90), c(allCIWidthsRisk90), 
            c(allCIWidthsPrevalence90), c(allCIWidthsGriddedRisk90))
tempRes = resolutions[col(allCIWidthsSmoothRisk90)]
tempCon = factor(as.character(poppsubSimple$subarea[col(allCIWidthsSmoothRisk90)]))
N=length(tempCon)
CIWidthFrame = data.frame(Constituency=rep(tempCon, 4), 
                          Resolution=rep(tempRes, 4), 
                          Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                         rep("Prevalence", N), rep("Gridded risk", N)), 
                                       levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                          CIWidth=CIWidth)

pdf(paste0("figures/gridResolutionTest/CIWidthVRes90.pdf"), width=7, height=5)
ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
  labs(x="Grid resolution (km)", y="90% credible interval width", fill="Model") + 
  theme_classic()
dev.off()

## 95% CIs
CIWidth = c(c(allCIWidthsSmoothRisk95), c(allCIWidthsRisk95), 
            c(allCIWidthsPrevalence95), c(allCIWidthsGriddedRisk95))
tempRes = resolutions[col(allCIWidthsSmoothRisk95)]
tempCon = factor(as.character(poppsubSimple$subarea[col(allCIWidthsSmoothRisk95)]))
N=length(tempCon)
CIWidthFrame = data.frame(Constituency=rep(tempCon, 4), 
                          Resolution=rep(tempRes, 4), 
                          Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                         rep("Prevalence", N), rep("Gridded risk", N)), 
                                       levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                          CIWidth=CIWidth)

pdf(paste0("figures/gridResolutionTest/CIWidthVRes95.pdf"), width=7, height=5)
ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
  labs(x="Grid resolution (km)", y="95% credible interval width", fill="Model") + 
  theme_classic()
dev.off()

### mean coverages

pdf(paste0("figures/gridResolutionTest/meanCoverage80VRes.pdf"), width=5, height=5)
pchs = 15:18
cols = rainbow(4)
ylim = range(c(meanCoveragesSmoothRisk80), c(meanCoveragesRisk80), 
             c(meanCoveragesPrevalence80), c(meanCoveragesGriddedRisk80))
ylim = c(0, 1)
plot(resolutions*.97, meanCoveragesSmoothRisk80, pch=pchs[1], col=cols[1], 
     ylim=ylim, ylab="80% coverage", xlab="Grid resolution (km)", 
     log="x")
abline(a=.8, b=0, lty=2)
points(resolutions*.99, meanCoveragesRisk80, pch=pchs[2], col=cols[2])
points(resolutions*1.01, meanCoveragesPrevalence80, pch=pchs[3], col=cols[3])
points(resolutions*1.03, meanCoveragesGriddedRisk80, pch=pchs[4], col=cols[4])
legend("right", c("Smooth risk", "Risk", "Prevalence", "Gridded risk"), 
       pch=pchs, col=cols)
dev.off()

pdf(paste0("figures/gridResolutionTest/meanCoverage90VRes.pdf"), width=5, height=5)
pchs = 15:18
cols = rainbow(4)
ylim = range(c(meanCoveragesSmoothRisk90), c(meanCoveragesRisk90), 
             c(meanCoveragesPrevalence90), c(meanCoveragesGriddedRisk90))
ylim = c(0, 1)
plot(resolutions*.97, meanCoveragesSmoothRisk90, pch=pchs[1], col=cols[1], 
     ylim=ylim, ylab="90% coverage", xlab="Grid resolution (km)", 
     log="x")
abline(a=.8, b=0, lty=2)
points(resolutions*.99, meanCoveragesRisk90, pch=pchs[2], col=cols[2])
points(resolutions*1.01, meanCoveragesPrevalence90, pch=pchs[3], col=cols[3])
points(resolutions*1.03, meanCoveragesGriddedRisk90, pch=pchs[4], col=cols[4])
legend("right", c("Smooth risk", "Risk", "Prevalence", "Gridded risk"), 
       pch=pchs, col=cols)
dev.off()

pdf(paste0("figures/gridResolutionTest/meanCoverage95VRes.pdf"), width=5, height=5)
pchs = 15:18
cols = rainbow(4)
ylim = range(c(meanCoveragesSmoothRisk95), c(meanCoveragesRisk95), 
             c(meanCoveragesPrevalence95), c(meanCoveragesGriddedRisk95))
ylim = c(0, 1)
plot(resolutions*.97, meanCoveragesSmoothRisk95, pch=pchs[1], col=cols[1], 
     ylim=ylim, ylab="95% coverage", xlab="Grid resolution (km)", 
     log="x")
abline(a=.8, b=0, lty=2)
points(resolutions*.99, meanCoveragesRisk95, pch=pchs[2], col=cols[2])
points(resolutions*1.01, meanCoveragesPrevalence95, pch=pchs[3], col=cols[3])
points(resolutions*1.03, meanCoveragesGriddedRisk95, pch=pchs[4], col=cols[4])
legend("right", c("Smooth risk", "Risk", "Prevalence", "Gridded risk"), 
       pch=pchs, col=cols)
dev.off()

### boxplots of all coverages

coverages = 100*c(c(allCoveragesSmoothRisk80), c(allCoveragesRisk80), 
            c(allCoveragesPrevalence80), c(allCoveragesGriddedRisk80))
tempRes = resolutions[row(allCoveragesSmoothRisk80)]
tempTruth = factor(as.character(col(allCoveragesSmoothRisk80)))
N=length(tempRes)
coverageFrame = data.frame(Truth=rep(tempTruth, 4), 
                          Resolution=rep(tempRes, 4), 
                          Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                         rep("Prevalence", N), rep("Gridded risk", N)), 
                                       levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                          coverages=coverages)

pdf(paste0("figures/gridResolutionTest/allCoveragesVRes80.pdf"), width=7, height=5)
ggplot(coverageFrame, aes(factor(Resolution), coverages, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + geom_hline(yintercept=80, linetype = 'dotted') + 
  labs(x="Grid resolution (km)", y="80% credible interval coverage", fill="Model") + 
  theme_classic()
dev.off()

coverages = 100*c(c(allCoveragesSmoothRisk90), c(allCoveragesRisk90), 
                  c(allCoveragesPrevalence90), c(allCoveragesGriddedRisk90))
tempRes = resolutions[row(allCoveragesSmoothRisk90)]
tempTruth = factor(as.character(col(allCoveragesSmoothRisk90)))
N=length(tempRes)
coverageFrame = data.frame(Truth=rep(tempTruth, 4), 
                           Resolution=rep(tempRes, 4), 
                           Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                          rep("Prevalence", N), rep("Gridded risk", N)), 
                                        levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                           coverages=coverages)

pdf(paste0("figures/gridResolutionTest/allCoveragesVRes90.pdf"), width=7, height=5)
ggplot(coverageFrame, aes(factor(Resolution), coverages, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + geom_hline(yintercept=90, linetype = 'dotted') + 
  labs(x="Grid resolution (km)", y="90% credible interval coverage", fill="Model") + 
  theme_classic()
dev.off()

coverages = 100*c(c(allCoveragesSmoothRisk95), c(allCoveragesRisk95), 
                  c(allCoveragesPrevalence95), c(allCoveragesGriddedRisk95))
tempRes = resolutions[row(allCoveragesSmoothRisk95)]
tempTruth = factor(as.character(col(allCoveragesSmoothRisk95)))
N=length(tempRes)
coverageFrame = data.frame(Truth=rep(tempTruth, 4), 
                           Resolution=rep(tempRes, 4), 
                           Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                          rep("Prevalence", N), rep("Gridded risk", N)), 
                                        levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                           coverages=coverages)

pdf(paste0("figures/gridResolutionTest/allCoveragesVRes95.pdf"), width=7, height=5)
ggplot(coverageFrame, aes(factor(Resolution), coverages, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + geom_hline(yintercept=95, linetype = 'dotted') + 
  labs(x="Grid resolution (km)", y="95% credible interval coverage", fill="Model") + 
  theme_classic()
dev.off()




