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
}









