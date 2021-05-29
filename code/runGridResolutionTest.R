# this script is for comparing performance between models

simDatKenya = generateSimDataSetsLCPB(nsim=1, adjustedPopMat=popMatSimpleAdjusted, 
                                      fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                      logisticApproximation=FALSE, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      seed=1, inla.seed=1L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                      easpa=NULL, popMat=NULL, constituencyPop=poppcon)

# get the data and the true population
dat = simDatKenya$SRSDat$clustDat[[1]]
truePrevalenceConstituencyKenya = simDatKenya$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[,1]
truePrevalenceCountyKenya = simDatKenya$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$countyMatrices$p[,1]

# construct integration grids at different resolutions
# resolutions = c(1, 10, 20, 40, 60, 80)
# deltas = c(.05, .2, .4, .8, 1.2, 1.4)
# meanNeighbors = c(500, rep(50, 5))
resolutions = c(.5, 1, 2, 5, 10, 20)
deltas = c(.025, .05, .1, .1, .2, .4)
meanNeighbors = c(rep(500, 3), rep(50, 3))
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
nSamples = 10000
spdeFitN = fitSPDEKenyaDat(dat, nPostSamples=nSamples, popMat=popMatCombined)
sigmaEpsilonDraws = spdeFitN$sigmaEpsilonDraws

# apply aggregation models at each resolution
aggResultsN = list()
for(i in 1:length(popGrids)) {
  # obtain the grids at this resolution
  thisPopMat = popGrids[[i]]
  thisPopMatAdjusted = popGridsAdjusted[[i]]
  
  # obtain model output at this resolution
  thisResolutionI = startIs[i]:endIs[i]
  thisUDraws = spdeFitN$uDraws[thisResolutionI,]
  
  thisAggResultsN = modLCPB(thisUDraws, sigmaEpsilonDraws, easpaN, thisPopMat, 
                            thisPopMatAdjusted, doLCPb=TRUE, doIHMEModel=TRUE, 
                            constituencyPop=poppconN, ensureAtLeast1PerConstituency=TRUE, 
                            logisticApproximation=FALSE, verbose=TRUE, 
                            fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                            stopOnFrameMismatch=FALSE)
  
  aggResultsN = c(aggResultsN, list(thisAggResultsN))
}
names(aggResultsN) = paste0("aggResultsN", c(5, resolutions))
save(popGrids, popGridsAdjusted, aggResultsN, file="savedOutput/simpleExample/gridResolutionTestNairobi.RData")





