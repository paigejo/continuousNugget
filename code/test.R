# functions for testing purposes

# simulate and save a population and set of surveys for testing purposes with parameters similar to 
# Kenya NMR model fits
# Simulate with: 
#     SPDE
#     Same hyperparameters as in NMR paper (sigmaEpsilon = 0.463, sigmaS^2 = 0.243)
#     Same fixed effects as in NMR paper (intercept = -4, Urban = 0.079)


# calculate expectation and variance of the logit-normal distribution
# code based on logitnorm package
logitNormMoments = function(muSigmaMat, parClust=NULL, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      t(apply(muSigmaMat, 1, logitNormMoments, ...))
    else
      t(parApply(parClust, muSigmaMat, 1, logitNormMoments, ...))
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
    expectation = integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
    
    fVar = function(x) (plogis(x) - expectation)^2 * dnorm(x, mean = mu, sd = sigma)
    variance = integrate(fVar, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
    c(mean = expectation, var = variance)
  }
}

# 
getLogitNormalApproximation = function(u=-1, sigmaEpsilon=0.1, n=10, qs = expit(rnorm(n))) {
  qs = qs/sum(qs)
  qs2 = qs^2
  varEpsilon = sigmaEpsilon^2
  
  multivariateExpitJacobian = function(z) {
    # based on results from WolframAlpha
    normalizer = (1 + sum(exp(z)))^2
    mat = outer(z, z, FUN = function(zi, zj) {-exp(zi + zj)})
    diag(mat) = diag(mat) + exp(z)*(1+sum(exp(z)))
    mat * (1 / normalizer)
  }
  
  ## calculate mean and variance on the probability scale
  # C = 16 * sqrt(3) / (15 * pi)
  # C2 = C^2
  # P = expit(u / sqrt(1 + C2 * varEpsilon))
  
  out = logitNormMoments(cbind(u, varEpsilon))
  P = out[1]
  V = out[2] * sum(qs2)
  V = sum(qs2) * c(multivariateExpitJacobian(u))^2 * varEpsilon
  ## determine logit normal distribution corresponding to this mean and variance
  cost = function(pars) {
    m = pars[1]
    v = exp(pars[2])
    out = logitNormMoments(cbind(m, v))
    thisP = out[1]
    thisV = out[2]
    (P - thisP)^2 + (V - thisV)^2
  }
  opt = optim(c(0, 1), cost)
  browser()
}

# test the CPBL model using the SPDE fit
testLCPB = function(fit=NULL) {
  # load SPDE UC model (fit it if necessary)
  if(is.null(fit)) {
    fitFile = paste0(outputDirectory, "test/spdeMortFitUC.RData")
    fitFileCompact = paste0(outputDirectory, "test/spdeMortFitUCCompact.RData")
    if(file.exists(fitFileCompact)) {
      out = load(fitFileCompact)
    } else if(file.exists(fitFile)) {
      out = load(fitFile)
      fit$mod = NULL
      save(fit, file=fitFileCompact)
    } else {
      fit = fitSPDEKenyaDat()
      save(fit, file=fitFile)
    }
  }
  
  uDraws = logit(fit$uDraws)
  sigmaEpsilonDraws = fit$sigmaEpsilonDraws
  
  # aggregate using CPBL model
  out = modLCPB(uDraws, sigmaEpsilonDraws, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                includeUrban=TRUE, maxEmptyFraction=1)
  
  browser()
}

testrbinom1 = function(n=1000000, size=100, prob=.02, maxPlot=size) {
  simulations = rbinom1(n, size, prob)
  
  xs = 1:size
  ps = dbinom1(xs, size, prob)
  
  ylim=range(ps)
  hist(simulations, breaks=seq(from=0.5, to=size+0.5, by=1), freq=FALSE, col="skyblue", xlim=c(0, maxPlot))
  lines(xs, ps, type="l", ylim=ylim)
  
  browser()
}
# testrbinom1(prob=20 / (2011 * 1000), size=2011, maxPlot=5)
# head(ps)
# [1] 9.900382e-01 9.895557e-03 6.590543e-05 3.290388e-07 1.313550e-09 4.367659e-12

# testrbinom1(prob=100 / (211 * 1000), size=211, maxPlot=5)
# head(ps)
# [1] 9.510471e-01 4.734943e-02 1.564095e-03 3.856470e-05 7.570311e-07 1.232404e-08

# testrbinom1(prob=1000 / (2011 * 1000), size=2011, maxPlot=10)
# head(ps)
# [1] 0.5820372099 0.2910186050 0.0969579399 0.0242153661 0.0048358447 0.0008043702

# testrbinom1(prob=1000 / (211 * 1000), size=211, maxPlot=10)
# head(ps)
# [1] 0.5825546364 0.2912773182 0.0966300944 0.0239274520 0.0047171262 0.0007712127

testrbinomTrunc = function(n=1000000, size=100, prob=.02, maxPlot=size, low=1, high=size) {
  simulations = rbinomTrunc(n, size, prob, low, high)
  
  xs = 1:size
  ps = dbinomTrunc(xs, size, prob, low, high)
  
  ylim=range(ps)
  hist(simulations, breaks=seq(from=0.5, to=size+0.5, by=1), freq=FALSE, col="skyblue", xlim=c(0, maxPlot))
  lines(xs, ps, type="l", ylim=ylim)
  
  browser()
}

testrmultinom1 = function(n=100000, size=400, k=30, probs=1.1^(1:k), ylim=NULL) {
  probs = probs / sum(probs)
  
  timeMult = system.time(simulationsMult <- rmultinom1(n, size, probs, method="mult", verbose=TRUE))[3]
  timeMult1 = system.time(simulationsMult1 <- rmultinom1(n, size, probs, method="mult1", verbose=TRUE))[3]
  timeIndepMH = system.time(simulationsIndepMH <- rmultinom1(n, size, probs, method="indepMH", verbose=TRUE))[3]
  
  allTimes = matrix(c(timeMult, timeMult1, timeIndepMH), nrow=1)
  colnames(allTimes) = c("timeMult", "timeMult1", "timeIndepMH")
  rownames(allTimes) = paste("n =", n, "; size =", size, "; k =", k)
  print(allTimes)
  
  par(mfrow=c(1, 3))
  maxX1 = max(simulationsMult[1,], simulationsMult1[1,], simulationsIndepMH[1,])
  maxXMid = max(simulationsMult[round(k/2),], simulationsMult1[round(k/2),], simulationsIndepMH[round(k/2),])
  maxXk = max(simulationsMult[k,], simulationsMult1[k,], simulationsIndepMH[k,])
  
  meanTestFun = function(prob) {
    sum((1:(size-k+1))*dbinom1(1:(size-k+1), prob=prob, size=size-k+1))
  }
  meansTest = sapply(probs, meanTestFun)[c(1, round(k/2), k)]
  
  pdf(paste0("figures/illustrations/truncatedMultinomialTestProb1_n", n, "_k", k, ".pdf"), 
      width=8, height=5)
  hist(simulationsMult[1,], main=paste0("Standard Multinomial, prob = ", round(probs[1], digits=4)), breaks=seq(from=0.5, to=maxX1+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxX1+0.5))
  abline(v=mean(simulationsMult[1,]), col="purple")
  abline(v=meansTest[1], col="green")
  legend("topright", c("Sample mean", "Truncated binomial mean"), col=c("purple", "green"), lty=1)
  hist(simulationsMult1[1,], main=paste0("Multinomial Plus 1, prob = ", round(probs[1], digits=4)), breaks=seq(from=0.5, to=maxX1+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxX1+0.5))
  abline(v=mean(simulationsMult1[1,]), col="purple")
  abline(v=meansTest[1], col="green")
  hist(simulationsIndepMH[1,], main=paste0("Independent MH, prob = ", round(probs[1], digits=4)), breaks=seq(from=0.5, to=maxX1+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxX1+0.5))
  abline(v=mean(simulationsIndepMH[1,]), col="purple")
  abline(v=meansTest[1], col="green")
  dev.off()
  
  pdf(paste0("figures/illustrations/truncatedMultinomialTestProbMid_n", n, "_k", k, ".pdf"), 
      width=8, height=5)
  hist(simulationsMult[round(k/2),], main=paste0("Standard Multinomial, prob = ", round(probs[round(k/2)], digits=4)), breaks=seq(from=0.5, to=maxXMid+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxXMid+0.5))
  abline(v=mean(simulationsMult[round(k/2),]), col="purple")
  abline(v=meansTest[2], col="green")
  legend("topright", c("Sample mean", "Truncated binomial mean"), col=c("purple", "green"), lty=1)
  hist(simulationsMult1[round(k/2),], main=paste0("Multinomial Plus 1, prob = ", round(probs[round(k/2)], digits=4)), breaks=seq(from=0.5, to=maxXMid+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxXMid+0.5))
  abline(v=mean(simulationsMult1[round(k/2),]), col="purple")
  abline(v=meansTest[2], col="green")
  hist(simulationsIndepMH[round(k/2),], main=paste0("Independent MH, prob = ", round(probs[round(k/2)], digits=4)), breaks=seq(from=0.5, to=maxXMid+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxXMid+0.5))
  abline(v=mean(simulationsIndepMH[round(k/2),]), col="purple")
  abline(v=meansTest[2], col="green")
  dev.off()
  
  pdf(paste0("figures/illustrations/truncatedMultinomialTestProbk_n", n, "_k", k, ".pdf"), 
      width=8, height=5)
  hist(simulationsMult[k,], main=paste0("Standard Multinomial, prob = ", round(probs[k], digits=4)), breaks=seq(from=0.5, to=maxXk+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxXk+0.5))
  abline(v=mean(simulationsMult[k,]), col="purple")
  abline(v=meansTest[3], col="green")
  legend("topright", c("Sample mean", "Truncated binomial mean"), col=c("purple", "green"), lty=1)
  hist(simulationsMult1[k,], main=paste0("Multinomial Plus 1, prob = ", round(probs[k], digits=4)), breaks=seq(from=0.5, to=maxXk+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxXk+0.5))
  abline(v=mean(simulationsMult1[k,]), col="purple")
  abline(v=meansTest[3], col="green")
  hist(simulationsIndepMH[k,], main=paste0("Independent MH, prob = ", round(probs[k], digits=4)), breaks=seq(from=0.5, to=maxXk+0.5, by=1), 
       freq=FALSE, col="skyblue", xlim=c(0, maxXk+0.5))
  abline(v=mean(simulationsIndepMH[k,]), col="purple")
  abline(v=meansTest[3], col="green")
  dev.off()
  
  browser()
}
# probs=1.1^(1:k): 
#                              timeMult timeMult1 timeIndepMH
# n = 1000 ; size = 40 ; k = 5    0.028     1.759       0.259

# [1] "acceptance percentage: 0.880217"
#                               timeMult timeMult1 timeIndepMH
# n = 1e+05 ; size = 40 ; k = 5   10.265   147.408      25.035

# [1] "acceptance percentage: 0.6484"
#                               timeMult timeMult1 timeIndepMH
# n = 1000 ; size = 40 ; k = 20   14.167    12.948       0.379

# [1] "acceptance percentage: 0.7805"
#                               timeMult  timeMult1 timeIndepMH
# n = 1000 ; size = 40 ; k = 30    1.248*     1.141       0.385
# *: method switched to mult1 since too many samples were expected to be required

# [1] "acceptance percentage: 0.3123"
#                                 timeMult  timeMult1 timeIndepMH
# n = 1000 ; size = 400 ; k = 30    0.593*     0.366*       0.363
# *: method switched to timeIndepMH since too many samples were expected to be required

# [1] "acceptance percentage: 0.319578"
#                                  timeMult  timeMult1 timeIndepMH
# n = 1e+05 ; size = 400 ; k = 30   42.137*    42.732*      40.646
# *: method switched to timeIndepMH since too many samples were expected to be required

# test rmultinom1 on the actual parameters relevant for validation
testrmultinom1Empirical = function(n=1000, sampleTable=NULL, stratifiedValidation=TRUE, ylim=NULL) {
  
  # set up sample table of indices if using stratified validation
  if(stratifiedValidation && is.null(sampleTable)) {
    out = getValidationI(dat=mort, pixelLevel=TRUE)
    sampleTable = out$sampleMatrix
    sampleListPixel = out$pixelIndexList
  }
  
  # initialize relevant vectors, matrices
  tmp = aggregate(mort$y, by=list(county=mort$admin1, urban=mort$urban), FUN=length)
  strata = tmp[,1:2]
  ks = matrix(NA, nrow=nrow(strata), ncol=ncol(sampleTable))
  sizes = matrix(NA, nrow=nrow(strata), ncol=ncol(sampleTable))
  times = matrix(NA, nrow=nrow(strata), ncol=ncol(sampleTable))
  acceptanceProbs = matrix(NA, nrow=nrow(strata), ncol=ncol(sampleTable))
  easpa = makeDefaultEASPA(validationClusterI=validationClusterI, useClustersAsEAs=!is.null(validationClusterI))
  
  timeMult = system.time(simulationsMult <- rmultinom1(n, size, probs, method="mult", verbose=TRUE))[3]
  timeMult1 = system.time(simulationsMult1 <- rmultinom1(n, size, probs, method="mult1", verbose=TRUE))[3]
  timeIndepMH = system.time(simulationsIndepMH <- rmultinom1(n, size, probs, method="indepMH", verbose=TRUE))[3]
}

testrMultinomTrunc = function(n=1000000, maxPlot=size, low=1, high=size) {
  size=20
  prob = rep(1/7, 7)
  x1 = rbinomTrunc(n, size, prob[1], low, min(high, size-7+1))
  x2 = rbinomTrunc(n, size-x1, prob[2]/(1 - prob[1]), low, sapply(x1, function(x) {min(high, size-6+1-x)}))
  x3 = rbinomTrunc(n, size-x1-x2, prob[3]/(1 - sum(prob[1:2])), low, sapply(x1+x2, function(x) {min(high, size-5+1-x)}))
  x4 = rbinomTrunc(n, size-x1-x2-x3, prob[4]/(1 - sum(prob[1:3])), low, sapply(x1+x2+x3, function(x) {min(high, size-4+1-x)}))
  x5 = rbinomTrunc(n, size-x1-x2-x3-x4, prob[5]/(1 - sum(prob[1:4])), low, sapply(x1+x2+x3+x4, function(x) {min(high, size-3+1-x)}))
  x6 = rbinomTrunc(n, size-x1-x2-x3-x4-x5, prob[6]/(1 - sum(prob[1:5])), low, sapply(x1+x2+x3+x4+x5, function(x) {min(high, size-2+1-x)}))
  x7 = rbinomTrunc(n, size-x1-x2-x3-x4-x5-x6, prob[7]/(1 - sum(prob[1:6])), low, sapply(x1+x2+x3+x4+x5+x6, function(x) {min(high, size-1+1-x)}))
  xTest = 1 + rbinom(n, size-7, prob=prob[1])
  
  binTruncMean1 = sum(low:high*dbinomTrunc(low:high, size, prob[1], low, min(high, size-7+1)))
  binTruncMean2 = sum(low:high*dbinomTrunc(low:high, size, prob[2], low, min(high, size-7+1)))
  binTruncMean3 = sum(low:high*dbinomTrunc(low:high, size, prob[3], low, min(high, size-7+1)))
  binTruncMean4 = sum(low:high*dbinomTrunc(low:high, size, prob[4], low, min(high, size-7+1)))
  binTruncMean5 = sum(low:high*dbinomTrunc(low:high, size, prob[5], low, min(high, size-7+1)))
  binTruncMean6 = sum(low:high*dbinomTrunc(low:high, size, prob[6], low, min(high, size-7+1)))
  binTruncMean7 = sum(low:high*dbinomTrunc(low:high, size, prob[7], low, min(high, size-7+1)))
  testMean = sum(low:(size-7+1)*dbinom(0:(size-7), size-7, prob[1]))
  
  means = c(mean(x1), mean(x2), mean(x3), mean(x4), mean(x5), mean(x6), mean(x7), mean(xTest))
  
  xlim = c(1,max(c(x1, x2, x3, x4, x5, x6, x7)))
  ylim=c(0, .33)
  par(mfrow=c(2,4))
  hist(x1, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[1], col="blue")
  abline(v=binTruncMean1, col="black")
  hist(x2, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[2], col="blue")
  abline(v=binTruncMean2, col="black")
  hist(x3, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[3], col="blue")
  abline(v=binTruncMean3, col="black")
  
  hist(x4, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[4], col="blue")
  abline(v=binTruncMean4, col="black")
  hist(x5, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[5], col="blue")
  abline(v=binTruncMean5, col="black")
  hist(x6, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[6], col="blue")
  abline(v=binTruncMean6, col="black")
  
  hist(x7, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[7], col="blue")
  abline(v=binTruncMean7, col="black")
  
  hist(xTest, breaks=seq(from=low-0.5, to=high+0.5, by=1), freq=FALSE, col="skyblue", ylim=ylim, xlim=xlim)
  abline(v=means[8], col="blue")
  abline(v=testMean, col="black")
  
  browser()
}
# testrbinomTrunc(prob=20 / (2011 * 1000), size=2011, maxPlot=5)
# head(ps)
# 9.900382e-01 9.895557e-03 6.590543e-05 3.290388e-07 1.313550e-09 4.367659e-12

# testrbinomTrunc(prob=100 / (211 * 1000), size=211, maxPlot=5)
# head(ps)
# 9.510471e-01 4.734943e-02 1.564095e-03 3.856470e-05 7.570311e-07 1.232404e-08

# testrbinomTrunc(prob=.5, size=10, maxPlot=10, low=2, high=5)
# head(ps)

# testrbinomTrunc(prob=1000 / (211 * 1000), size=211, maxPlot=10)
# head(ps)


testrpois1 = function(n=1000000, prob=20 / (211 * 1000), maxPlot=qpois1(.99, prob)) {
  simulations = rpois1(n, prob)
  
  xs = 1:max(simulations)
  ps = dpois1(xs, prob)
  
  ylim=range(ps)
  hist(simulations, breaks=seq(from=0.5, to=max(simulations)+0.5, by=1), freq=FALSE, col="skyblue", xlim=c(0, maxPlot))
  lines(xs, ps, type="l", ylim=ylim)
  
  browser()
}

testValidationIndices = function(validationOut=NULL) {
  popMat = makeDefaultPopMat()
  
  if(is.null(validationOut)) {
    validationOut = getValidationI(pixelLevel=TRUE, popMat=popMat)
  }
  
  pixelIndices = getPixelIndex(cbind(mort$east, mort$north), popMat=popMat)
  
  indexMatrix = validationOut$sampleMatrix
  pixelIndexMatrix = validationOut$pixelIndexList
  browser()
  # load shape files for plotting
  require(maptools)
  regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  
  # do the plotting
  par(mfrow=c(1,2))
  cols = rainbow(10)
  for(i in 1:10) {
    theseIndices = indexMatrix[,i]
    
    if(i == 1) {
      plot(mort$lon, mort$lat, type="n", xlim=kenyaLonRange, ylim=kenyaLatRange)
    }
    
    points(mort$lon[theseIndices], mort$lat[theseIndices], pch=19, col=cols[i], cex=.2)
  }
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  for(i in 1:10) {
    thesePixelIndices = pixelIndexMatrix[[i]]
    
    if(i == 1) {
      plot(mort$lon, mort$lat, type="n", xlim=kenyaLonRange, ylim=kenyaLatRange)
    }
    
    points(popMat$lon[thesePixelIndices], popMat$lat[thesePixelIndices], pch=19, col=cols[i], cex=.2)
  }
  plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
  plotMapDat(mapDat=regionMap, lwd=2.5)
  
  for(i in 1:10) {
    thesePixelIndices = pixelIndexMatrix[[i]]
    theseClusterIndices = indexMatrix[,i]
    
    thesePixelAreas = popMat$area[sort(thesePixelIndices)]
    theseClusterAreas = popMat$area[sort(unique(pixelIndices[theseClusterIndices]))]
    
    print(mean(thesePixelAreas == theseClusterAreas))
    if(any(thesePixelAreas != theseClusterAreas)) {
      stop()
    }
  }
}

testValidationWeights = function() {
  nus = runif(47)
  sigmas = runif(47)
  lambda = -2/47 * sum(sigmas^2 / (sigmas^2 + nus^2)) / sum(1/(nus^2 + sigmas^2))
  ais = (nus^2/47 -lambda/2)/(nus^2 + sigmas^2)
  
  # this should be 1:
  print(sum(ais))
  
  plot((nus^2/(nus^2 + sigmas^2))/sum(nus^2/(nus^2 + sigmas^2)), ais, pch=19, col="blue", 
       main="ais versus proportion of signal variance", xlim=range(ais))
  abline(a=0, b=1)
  
}







