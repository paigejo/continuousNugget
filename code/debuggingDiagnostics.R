meanTruths = 1:100
meanDatsN = 1:100
meanDats = 1:100
meanDatsUrban = 1:100
meanDatsRural = 1:100
meanEsts = 1:100
for(i in 1:length(truths)) {
  thisTruth = truths[[i]]
  thisDat = dats[[i]]
  
  truePrevalenceConstituencyKenya = thisTruth$truePrevalenceConstituencyKenya
  truePrevalenceDatKenyaN = thisDat$pFineScalePrevalence[thisDat$area == "Nairobi"]
  truePrevalenceDatKenyaUrban = thisDat$pFineScalePrevalence[thisDat$urban]
  truePrevalenceDatKenyaRural = thisDat$pFineScalePrevalence[!thisDat$urban]
  
  predsSmoothRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
  for(j in 1:length(resolutions)) {
    thisMaxSamples = 1000
    predsSmoothRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pSmoothRisk[,1:thisMaxSamples])
  }
  predsSmoothRisk = rowMeans(predsSmoothRisk)
  
  meanTruths[i] = mean(truePrevalenceConstituencyKenya)
  meanDatsN[i] = mean(truePrevalenceDatKenyaN)
  meanDats[i] = mean(thisDat$pFineScalePrevalence)
  meanDatsUrban[i] = mean(truePrevalenceDatKenyaUrban)
  meanDatsRural[i] = mean(truePrevalenceDatKenyaRural)
  meanEsts[i] = mean(predsSmoothRisk)
}
sigma = sqrt((1/3)^2 + 1/2.5)
mean(meanDatsRural) # 0.06303736
mean(meanTruths) # 0.02428097
mean(meanDatsN) # 0.02432309
mean(meanDatsUrban) # 0.02446275
mean(meanEsts) # 0.02775667
logitNormMean(rbind(c(-2.9, sigma),  # 0.06391782 0.02507139
                    c(-3.9, sigma)), logisticApproximation=FALSE)
quantile(meanEsts - meanTruths, probs=c(.025, .05, .1, .5, .9, .95, .975))
mean(meanEsts > meanTruths)
mean(meanEsts > meanDatsN)

temp = do.call("rbind", separateUDraws)
mean(temp)

logitNormMean(rbind(c(-2.9, sqrt(1/2.5)), 
                    c(-3.9, sqrt(1/2.5))), logisticApproximation=FALSE)
logitNormMean(matrix(c(-4.221521, sqrt(1/2.5)), nrow=1), logisticApproximation=FALSE)

# test variance of predictions
for(i in 1:length(truths)) {
  
}


# test SPDE predictions alone:
fixedParameters = list(spde=list(effRange=400, margVar=(1/3)^2), clusterPrec=2.5, beta=c(-2.9, -1))
fixedParameters = list(spde=list(effRange=400, margVar=(1/3)^2), clusterPrec=250, beta=c(logit(.8), 0))
fixedParameters = list(spde=list(effRange=400, margVar=(1/3)^2), clusterPrec=250, beta=c(logit(.2), 0))
fixedParameters = list(spde=list(effRange=400, margVar=(1/3)^2), clusterPrec=2.5, beta=c(logit(.8), logit(.2)-logit(.8)))
fixedParameters = list(spde=list(effRange=400, margVar=(1/3)^2), clusterPrec=2.5, beta=c(logit(.8), logit(.2)-logit(.8)))
nSamples = rep(1000, length(resolutions))
allSPDEResultsN = list()
startI = 1
endI = length(truths)
endI = 100
allTempDatMeansUrban = startI:endI
allTempDatMeansRural = startI:endI
for(i in startI:endI) {
  print(paste0("Running analysis for truth ", i, "/", length(truths), ", endI=", endI))
  
  dat = dats[[i]]
  tempDat = dat
  tempDat$Z[tempDat$urban] = rbinom(sum(tempDat$urban), 25, expit(-3.9+rnorm(sum(tempDat$urban), sd=sqrt(1/2.5))))
  tempDat$y[tempDat$urban] = tempDat$Z[tempDat$urban]
  tempDat$Z[!tempDat$urban] = rbinom(sum(!tempDat$urban), 25, expit(-2.9+rnorm(sum(!tempDat$urban), sd=sqrt(1/2.5))))
  tempDat$y[!tempDat$urban] = tempDat$Z[!tempDat$urban]
  allTempDatMeansUrban[i] = mean(tempDat$Z[tempDat$urban]/25)
  allTempDatMeansRural[i] = mean(tempDat$Z[!tempDat$urban]/25)
  spdeFitN = fitSPDEKenyaDat(tempDat, nPostSamples=max(nSamples), popMat=popGrids[[2]],
                             fixedParameters=fixedParameters, prior=NULL, strategy="simplified.laplace")
  # spdeFitN = fitSPDEKenyaDat(tempDat, nPostSamples=max(nSamples), popMat=popGrids[[2]])
  
  # apply SPDE model
  j = length(popGrids)
  # obtain the grids at this resolution
  thisNSamples = nSamples[j]
  
  # obtain SPDE model output at this resolution
  thisResolutionI = startIs[j]:endIs[j]
  allSPDEResultsN = c(allSPDEResultsN, spdeFitN$uDraws)
}

out = logitNormMoments(c(-3.9, sqrt((1/3)^2/100)))
c(out[1], sqrt(out[2]))

temp = c(sapply(allSPDEResultsN, as.numeric))
tempMeans = sapply(allSPDEResultsN, mean)
tempMeansProb = logitNormMean(cbind(tempMeans, sqrt(1/2.5)), logisticApproximation = FALSE)
tempMean = mean(temp)
vars = c(sapply(allSPDEResultsN, function(x) {var(as.numeric(x))}))
mean(temp)
# [1] -3.805226
withinVar = mean(vars)
logitNormMean(cbind(tempMean, sqrt(1/2.5 + withinVar)), logisticApproximation = FALSE)
# 0.02623688
logitNormMean(cbind(logit(.2), sqrt(1/2.5 + withinVar)), logisticApproximation = FALSE)
logitNormMean(cbind(-3.9, sqrt(1/2.5 + withinVar)), logisticApproximation = FALSE)
# 0.02393966
logitNormMean(cbind(-3.9, sqrt(1/2.5)), logisticApproximation = FALSE)
mean(allTempDatMeansUrban) # 0.02306388
mean(allTempDatMeansUrban <= logitNormMean(cbind(-3.9, sqrt(1/2.5)), logisticApproximation = FALSE))
# 0.87
mean(allTempDatMeansUrban <= logitNormMean(cbind(-3.9, sqrt(1/2.5+withinVar)), logisticApproximation = FALSE))
# 0.9

# using simplified Laplace:
# > mean(temp)
# [1] -3.926438
# > withinVar = mean(vars)
# > logitNormMean(cbind(tempMean, sqrt(1/2.5 + withinVar)), logisticApproximation = FALSE)
# [1] 0.02333609
# > logitNormMean(cbind(-3.9, sqrt(1/2.5 + withinVar)), logisticApproximation = FALSE)
# [1] 0.02393951

tempMeansProb2 = logitNormMean(cbind(tempMeans, sqrt(1/2.5 + withinVar)), logisticApproximation = FALSE)

# average logit scale predictions for Nairobi, fit to data, have .0061 
# probability to be as or more extreme under null distribution:
mean(sapply(allSPDEResultsN, mean)) # -3.805226
2*(1-pnorm(mean(temp), -3.9, sqrt(((1/3)^2 + withinVar)/100)))
# [1] 0.00608827
qnorm(c(.025, .05, .1, .9, .95, .975), -3.9, sqrt((1/3)^2/100))
# [1] -3.965332 -3.954828 -3.942718 -3.857282 -3.845172 -3.834668
qnorm(c(.025, .05, .1, .9, .95, .975), -3.9, sqrt(((1/3)^2 + withinVar)/100))
# [1] -3.967720 -3.956832 -3.944280 -3.855720 -3.843168 -3.832280



# check data urbanicity
sum(dat$urban)
sum(clustpc$clustUrb)

out = aggregate(dat$urban, by=list(area=dat$area), FUN=sum)
out2 = aggregate(!dat$urban, by=list(area=dat$area), FUN=sum)
cbind(out, out2, easpa25$popUrb, easpa25$popRur, easpa25$popTotal)



##### TEST 2
simDatKenya <- generateSimDataSetsLCPB2(nsim=1, targetPopMat=popMatKenyaNeonatal,
                                        popMat=popMatKenya,
                                        doFineScaleRisk=TRUE, doSmoothRisk=TRUE,
                                        gridLevel=FALSE, subareaLevel=TRUE,
                                        fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1,
                                        logisticApproximation=FALSE,
                                        dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/",
                                        seed=seeds[i], inla.seed=inlaSeeds[i],
                                        simPopOnly=FALSE, returnEAinfo=TRUE,
                                        easpa=easpaKenya, poppsub=poppsubKenya)

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
                     easpa=easpa, popMat=popMat, poppsub=poppsub)

outLCPB = simPopSPDE(nsim=nsim, easpa=easpa, popMat=popMat, targetPopMat=targetPopMat, 
                     poppsub=poppsub, spdeMesh=spdeMesh, 
                     margVar=margVar, sigmaEpsilon=sigmaEpsilon, 
                     gamma=gamma, effRange=effRange, beta0=beta0, 
                     seed=seed, inla.seed=inla.seed, nHHSampled=nHHSampled, 
                     stratifyByUrban=stratifyByUrban, 
                     subareaLevel=subareaLevel, gridLevel=gridLevel, 
                     doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                     doSmoothRiskLogisticApprox=logisticApproximation, 
                     min1PerSubarea=min1PerSubarea, 
                     fixPopPerEA=fixPopPerEA, fixHHPerEA=fixHHPerEA, fixPopPerHH=fixPopPerHH)








test = generateSimDataSetsLCPB2(nsim=1000, targetPopMat=popMatKenyaNeonatal, 
                                popMat=popMatKenya, 
                                doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
                                gridLevel=FALSE, subareaLevel=TRUE, 
                                fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                logisticApproximation=FALSE, 
                                dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                seed=seeds[i], inla.seed=inlaSeeds[i], 
                                simPopOnly=TRUE, returnEAinfo=TRUE, 
                                easpa=easpaKenya, poppsub=poppsubKenya)
mean(sapply(test$eaDatList, function(x) {mean(x$pFineScaleRisk[x$urban], na.rm=TRUE)}))
sigma = sqrt(1/2.5 + 1/9)
logitNormMean(rbind(c(-3.9, sigma), 
                    c(-2.9, sigma)), logisticApproximation=FALSE) # 0.02507139 0.06391782
logitNormMean(rbind(c(-3.9, sqrt(1/2.5)), 
                    c(-2.9, sqrt(1/2.5))), logisticApproximation=FALSE) # 0.02385813 0.06130256
logitNormMean(rbind(c(-3.9, sqrt(1/9)), 
                    c(-2.9, sqrt(1/9))), logisticApproximation=FALSE) # 0.02089987 0.05463981
mean(sapply(out$eaDatList, function(x) {mean(x$pSmoothRisk[x$urban], na.rm=TRUE)})) # 0.02480162
mean(sapply(out$eaDatList, function(x) {mean(x$pSmoothRisk[!x$urban], na.rm=TRUE)})) # 0.06332789

test = simVals + matrix(stats::rnorm(length(simVals), sd=sigmaEpsilon), ncol=nsim)


