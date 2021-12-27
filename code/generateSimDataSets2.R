# simulate and save datasets used for the simulation study with the given model parameters
# NOTE: paired with the dataset using the passed parameters will be another dataset from the 
#       same model without a nugget/cluster effect
# nsim: number of surveys taken from the true latent population in the standard size survey collections
# nsimBig: number of surveys taken from the true latent population in the large size survey collections
# seeds: random number seeds used for making the latent population and generating surveys respectively
# beta0: latent gaussian model intercept
# rho: marginal variance of the spatial field
# tausq: the nugget/cluster effect variance
# gamma: latent gaussian model urban effect
# HHoldVar: household effect variance
# effRange: spatial range
# urbanOverSamplefrac: the proportion with which to inflate the amount of urban samples in the surveys
generateSimDataSets = function(nsim=10, rho=0.243, sigmaEpsilon=sqrt(0.463), 
                               gamma = 0.009, effRange = 406.51, beta0=-3.922, 
                               figureSaveDirectory="~/git/continuousNugget/figures/simDataSets/", 
                               dataSaveDirectory="~/git/continuousNugget/savedOutput/simDataSets/", seed=123) {
  tausq = sigmaEpsilon^2
  set.seed(seed)
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Beta", round(beta0, 4), "rho", round(rho, 4), "sigmaEps", 
                  round(sigmaEpsilon, 4), "gamma", round(gamma, 4))
  
  # there should be 1 true population, but many simulated cluster surveys
  load(paste0(globalDirectory, "empiricalDistributions.RData"))
  simulatedEAs = simDatEmpirical(empiricalDistributions, kenyaEAs, clustDat=NULL, nsim=1, 
                                 beta0=beta0, margVar=rho, urbanOverSamplefrac=0, 
                                 tausq=tausq, gamma=gamma, HHoldVar=0, effRange=effRange)
  kenyaEAs = simulatedEAs$eaDat
  kenyaEAs$eaIs = 1:nrow(kenyaEAs)
  kenyaEAsLong = kenyaEAs[rep(1:nrow(kenyaEAs), kenyaEAs$nHH),]
  
  # simulate the cluster sampling and add to the data sets
  overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsim, NULL, urbanOverSamplefrac, verbose=FALSE)
  clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, overSampClustDat$eaIs, overSampClustDat$sampleWeights)
  overSampDat = list(eaDat=kenyaEAs, clustDat=clustList)
  
  SRSClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsim, NULL, SRS=TRUE, verbose=FALSE)
  clustList = genAndreaFormatFromEAIs(kenyaEAs, SRSClustDat$eaIs, SRSClustDat$sampleWeights)
  SRSDat = list(eaDat=kenyaEAs, clustDat=clustList) # the only thing different is the sampling of the clusters
  
  # plot the first simulation of the over sampled and simple random sample data sets
  clustDat = SRSDat$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  eaDat = overSampDat$eaDat
  pdf(paste0(figureSaveDirectory, "/unstratifiedSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # plot the first simulation of the over sampled and simple random sample data sets
  clustDat = overSampDat$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  eaDat = overSampDat$eaDat
  pdf(paste0(figureSaveDirectory, "/stratifiedSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  save(overSampDat, SRSDat, file=paste0(dataSaveDirectory, "simDataMulti", dataID, ".RData"))
  
  invisible(NULL)
}

generateSimPopulationsLCPB = function(nsim=10, rho=0.243, sigmaEpsilon=sqrt(0.463), 
                                      gamma = 0.009, effRange = 406.51, beta0=-3.922, 
                                      figureSaveDirectory="~/git/continuousNugget/figures/simDataSets/", 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simDataSets/", 
                                      seed=1) {
  set.seed(seed)
  
  
}

# simulate and save datasets used for the simulation study with the given model parameters from LCPB model
# NOTE: paired with the dataset using the passed parameters will be another dataset from the 
#       same model without a nugget/cluster effect
# nsim: number of surveys taken from the true latent population in the standard size survey collections
# nsimBig: number of surveys taken from the true latent population in the large size survey collections
# seeds: random number seeds used for making the latent population and generating surveys respectively
# beta0: latent gaussian model intercept
# rho: marginal variance of the spatial field
# tausq: the nugget/cluster effect variance
# gamma: latent gaussian model urban effect
# HHoldVar: household effect variance
# effRange: spatial range
# urbanOverSamplefrac: the proportion with which to inflate the amount of urban samples in the surveys
# fixPopPerEA: if not NULL, fix the target population in each EA to be this number
# spreadEAsInPixels: whether or not to draw exact EA locations uniformly from the pixel from which they are drawn, 
#                    or just to take the centroid of the pixel. NOTE: setting this to TRUE may jitter the EAs out 
#                    of the administrative area assigned to the pixel, but the area assigned to the EA will stay 
#                    the same
generateSimDataSetsLCPB2 = function(nsim=100, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), 
                                   gamma=-1, effRange = 400, beta0=-2.9, 
                                   gridLevel=FALSE, subareaLevel=TRUE, 
                                   doFineScaleRisk=TRUE, doSmoothRisk=FALSE, 
                                   fixPopPerEA=NULL, fixHHPerEA=NULL, fixPopPerHH=NULL, 
                                   logisticApproximation=FALSE, spreadEAsInPixels=TRUE, targetPopMat=NULL, 
                                   figureSaveDirectory="~/git/continuousNugget/figures/simDataSets/", 
                                   dataSaveDirectory="~/git/continuousNugget/savedOutput/simDataSets/", 
                                   seed=NULL, inla.seed=0L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                   easpa=NULL, popMat=NULL, poppsub=NULL, 
                                   verbose=TRUE, stopOnFrameMismatch=TRUE, thisclustpc=NULL, 
                                   nEAsFac=1, nClustFac=1, representativeSampling=FALSE) {
  tausq = sigmaEpsilon^2
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  if(!simPopOnly && !returnEAinfo) {
    warning("simPopOnly is FALSE, but returnEAinfo was set to FALSE. Setting returnEAinfo to TRUE instead")
    returnEAinfo = TRUE
  }
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Bet", signif(beta0, 3), "rho", signif(rho, 3), "sigEps", 
                  signif(sigmaEpsilon, 3), "gam", signif(gamma, 3), 
                  "nEAsFac", signif(nEAsFac, 3), "nClustFac", signif(nClustFac, 3))
  
  # there should be 1 true population, but many simulated cluster surveys
  # load(paste0(globalDirectory, "empiricalDistributions.RData"))
  # simulatedEAs = simDatEmpirical(empiricalDistributions, kenyaEAs, clustDat=NULL, nsim=1, 
  #                                beta0=beta0, margVar=rho, urbanOverSamplefrac=0, 
  #                                tausq=tausq, gamma=gamma, HHoldVar=0, effRange=effRange)
  simOut = simDatLCPB2(nsim=nsim, margVar=rho, sigmaEpsilon=sigmaEpsilon, 
                      gamma=gamma, effRange=effRange, beta0=beta0, 
                      stratifyByUrban=TRUE, targetPopMat=targetPopMat, 
                      gridLevel=gridLevel, subareaLevel=subareaLevel, 
                      doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                      min1PerSubarea=TRUE, inla.seed=inla.seed, 
                      spreadEAsInPixels=spreadEAsInPixels, 
                      fixPopPerEA=fixPopPerEA, fixHHPerEA=fixHHPerEA, fixPopPerHH=fixPopPerHH, 
                      logisticApproximation=logisticApproximation, 
                      simPopOnly=simPopOnly, returnEAinfo=returnEAinfo, 
                      verbose=verbose, stopOnFrameMismatch=stopOnFrameMismatch, 
                      easpa=easpa, popMat=popMat, poppsub=poppsub, 
                      nEAsFac=nEAsFac, nClustFac=nClustFac, 
                      representativeSampling=representativeSampling)
  
  if(!simPopOnly) {
    kenyaEAs = simOut$eaDat
    kenyaEAs$eaIs = 1:nrow(kenyaEAs)
    kenyaEAsLong = kenyaEAs[rep(1:nrow(kenyaEAs), kenyaEAs$nHH),]
    thisclustpc = simOut$thisclustpc 
    
    # simulate the cluster sampling and add to the data sets
    clustList = simOut$clustList
    overSampDat = list(eaDat=kenyaEAs, clustDat=clustList, aggregatedPop=simOut$aggregatedPop)
    
    # SRSClustDat = simClustersSRS(nsim, kenyaEAs, kenyaEAsLong, nEASampled=sum(thisclustpc$clustTotal))
    # clustList = genAndreaFormatFromEAIs(kenyaEAs, SRSClustDat$eaIs, SRSClustDat$sampleWeights)
    # SRSDat = list(eaDat=kenyaEAs, clustDat=clustList, aggregatedPop=simOut$aggregatedPop) # the only thing different is the sampling of the clusters
    
    # plot the first simulation of the over sampled and simple random sample data sets
    # clustDat = SRSDat$clustDat[[1]]
    # # clustDat = overSampDat$clustDat[[1]]
    # # eaDat = overSampDat$eaDat
    # pdf(paste0(figureSaveDirectory, "/unstratifiedSimulation", dataID, ".pdf"), width=8, height=8)
    # par(mfrow =c(2, 2), mar=c(5, 4, 4, 8))
    # obsCoords = cbind(clustDat$east, clustDat$north)
    # obsNs = clustDat$N
    # obsCounts = clustDat$Z
    # zlim = c(0, quantile(c(eaDat$pFineScalePrevalence, clustDat$pFineScalePrevalence, 
    #                        eaDat$pFineScaleRisk), probs=.975))
    # quilt.plot(eaDat$east, eaDat$north, eaDat$pFineScalePrevalence, main="Population Prevalences (LCPB)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Prevalences (LCPB)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # quilt.plot(eaDat$east, eaDat$north, eaDat$pFineScaleRisk, main="Population Risks (LCPb)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # quilt.plot(obsCoords, clustDat$pFineScaleRisk, main="Sample Risks (LCPb)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # # quilt.plot(eaDat$east, eaDat$north, eaDat$pLCpb, main="Population Risks (LCpb)", 
    # #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # # plotMapDat(project=TRUE)
    # # quilt.plot(obsCoords, clustDat$pLCpb, main="Sample Risks (LCpb)", 
    # #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # # plotMapDat(project=TRUE)
    # # quilt.plot(eaDat$east, eaDat$north, eaDat$pLcpb, main="Population Risks (Lcpb)", 
    # #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # # plotMapDat(project=TRUE)
    # # quilt.plot(obsCoords, clustDat$pLcpb, main="Sample Risks (Lcpb)", 
    # #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # # plotMapDat(project=TRUE)
    # dev.off()
    
    # plot the first simulation of the over sampled and simple random sample data sets
    clustDat = overSampDat$clustDat[[1]]
    # clustDat = overSampDat$clustDat[[1]]
    eaDat = overSampDat$eaDat
    pdf(paste0(figureSaveDirectory, "/stratifiedSimulation", dataID, ".pdf"), width=8, height=8)
    par(mfrow =c(2, 2), mar=c(5, 4, 4, 8))
    obsCoords = cbind(clustDat$east, clustDat$north)
    obsNs = clustDat$N
    obsCounts = clustDat$Z
    zlim = c(0, quantile(c(eaDat$pFineScalePrevalence, clustDat$pFineScalePrevalence, 
                           eaDat$pFineScaleRisk), probs=.975))
    quilt.plot(eaDat$east, eaDat$north, eaDat$pFineScalePrevalence, main="Population Prevalences (LCPB)", 
               xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    plotMapDat(project=TRUE)
    quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Prevalences (LCPB)", 
               xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    plotMapDat(project=TRUE)
    quilt.plot(eaDat$east, eaDat$north, eaDat$pFineScaleRisk, main="Population Risks (Fine Scale Risk)", 
               xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    plotMapDat(project=TRUE)
    quilt.plot(obsCoords, clustDat$pFineScaleRisk, main="Sample Risks (Fine Scale Risk)", 
               xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    plotMapDat(project=TRUE)
    # quilt.plot(eaDat$east, eaDat$north, eaDat$pLCpb, main="Population Risks (LCpb)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # quilt.plot(obsCoords, clustDat$pLCpb, main="Sample Risks (LCpb)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # quilt.plot(eaDat$east, eaDat$north, eaDat$pLcpb, main="Population Risks (Lcpb)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    # quilt.plot(obsCoords, clustDat$pLcpb, main="Sample Risks (Lcpb)", 
    #            xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
    # plotMapDat(project=TRUE)
    dev.off()
    
    save(overSampDat, SRSDat, file=paste0(dataSaveDirectory, "simDataMulti", dataID, ".RData"))
    
    return(invisible(list(simulatedEAs=simOut, overSampDat=overSampDat, SRSDat=SRSDat)))
  } else {
    invisible(simOut)
  }
}

## TODO: modify function to take in vector of seeds, make different population replication 
##       from same parameters for each seed
generateAllDataSets = function() {
  # generateSimDataSets(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-3.9)
  # generateSimDataSets(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-3.9)
  # generateSimDataSets(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  # generateSimDataSets(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  # generateSimDataSetsLCPB(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-2.9)
  # generateSimDataSetsLCPB(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-3.9)
  # generateSimDataSetsLCPB(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  # generateSimDataSetsLCPB(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  set.seed(1)
  generateSimDataSetsLCPB(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=150, beta0=-2.9)
  generateSimDataSetsLCPB(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=150, beta0=-3.9)
  generateSimDataSetsLCPB(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=150, beta0=0)
  generateSimDataSetsLCPB(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=150, beta0=0)
}

## TODO: modify function to take in vector of seeds, make different population replication 
##       from same parameters for each seed
generateAllDataSets2 = function() {
  # generateSimDataSets(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-3.9)
  # generateSimDataSets(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-3.9)
  # generateSimDataSets(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  # generateSimDataSets(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  # generateSimDataSetsLCPB(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-2.9)
  # generateSimDataSetsLCPB(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=-3.9)
  # generateSimDataSetsLCPB(gamma=-1, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  # generateSimDataSetsLCPB(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), effRange=400, beta0=0)
  
  out = load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
  
  nPar = length(spde_prevRiskSimStudyCommandArgs)
  set.seed(1)
  allSeeds = sample(1:1000000, nPar, replace=FALSE)
  for(i in nPar) {
    set.seed(allSeeds[i])
    
    theseArgs = spde_prevRiskSimStudyCommandArgs[[i]]
    do.call("generateSimDataSetsLCPB2", theseArgs)
  }
}

# simulate and save datasets used for the simulation study with the given model parameters
# NOTE: paired with the dataset using the passed parameters will be another dataset from the 
#       same model without a nugget/cluster effect
# nsim: number of surveys taken from the true latent population in the standard size survey collections
# nsimBig: number of surveys taken from the true latent population in the large size survey collections
# seeds: random number seeds used for making the latent population and generating surveys respectively
# beta0: latent gaussian model intercept
# margVar: marginal variance of the spatial field
# tausq: the nugget/cluster effect variance
# gamma: latent gaussian model urban effect
# HHoldVar: household effect variance
# effRange: spatial range
# urbanOverSamplefrac: the proportion with which to inflate the amount of urban samples in the surveys
plotSimDataSets = function(nsim=100, nsimBig = 250, seeds=c(580252, 1234), beta0 = -1.75, margVar = .15^2, 
                           tausq = .1^2, gamma = -1, HHoldVar = 0, effRange = 150, 
                           urbanOverSamplefrac = 0, colorScale=makeRedBlueDivergingColors(64, rev = TRUE), 
                           kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  set.seed(seeds[1])
  wd = getwd()
  setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
  
  rangeText = ""
  if(effRange != 150)
    rangeText = paste0("Range", effRange)
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                  round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                  "urbanOverSamplefrac", round(urbanOverSamplefrac, 4), rangeText)
  dataID0 = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                   round(0, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                   "urbanOverSamplefrac", round(urbanOverSamplefrac, 4), rangeText)
  
  # there should be 1 true data set, but many simulated cluster samples
  load("empiricalDistributions.RData")
  
  # save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "Big.RData"))
  load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # plot the first simulation of the over sampled and simple random sample data sets
  clustDat = SRSDat$clustDat[[1]]
  clustDat = overSampDat$clustDat[[1]]
  eaDat = overSampDat$eaDat
  pdf(paste0("figures/exampleSRSSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  png(file=paste0("figures/exampleOverSampSimulation", dataID, ".png"), width=900, height=500)
  par(oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6), mfrow =c(1, 2))
  
  plot(cbind(eaDat$lon, eaDat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
       xlab="Longitude", ylab="Latitude", main=paste0("Simulated NMRs"), asp=1)
  quilt.plot(eaDat$lon, eaDat$lat, eaDat$died/eaDat$numChildren, nx=150, ny=150, col=colorScale, zlim=zlim, add=TRUE)
  # world(add=TRUE)
  plotMapDat(adm1)
  
  plot(cbind(clustDat$lon, clustDat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
       xlab="Longitude", ylab="Latitude", main=paste0("Example Survey NMRs"), asp=1)
  quilt.plot(clustDat$lon, clustDat$lat, clustDat$died/clustDat$numChildren, nx=150, ny=150, col=colorScale, zlim=zlim, add=TRUE)
  # world(add=TRUE)
  plotMapDat(adm1)
  dev.off()
  
  clustDat = overSampDat$clustDat[[1]]
  pdf(paste0("figures/exampleUrbanOverSampSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  
  setwd(wd)
  
  invisible(NULL)
}

# function for calculating the average number of urban and rural clusters per simulated dataset
getAverageUrbanRuralClusters = function() {
  out = system("ls ~/Google\\ Drive/UW/Wakefield/WakefieldShared/U5MR/simDataMultiBeta*0Big.RData", intern=TRUE)
  numUrban = c()
  numRural = c()
  for(i in 1:length(out)) {
    print(paste0("Loading ", out[i], "... (", i, "/", length(out), ")"))
    load(out[i])
    numUrban = c(numUrban, sapply(SRSDat$clustDat, function(x){sum(x$urban)}))
    numRural = c(numRural, sapply(SRSDat$clustDat, function(x){sum(!x$urban)}))
  }
  
  print(paste0("Average number of urban clusters: ", mean(numUrban)))
  print(paste0("Average number of rural clusters: ", mean(numRural)))
  
  invisible(list(numUrban=numUrban, numRural=numRural))
}


