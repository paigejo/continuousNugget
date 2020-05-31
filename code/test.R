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
testCPBL = function(fit=NULL) {
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
  out = modCPBL(uDraws, sigmaEpsilonDraws, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                includeUrban=TRUE, maxEmptyFraction=1)
  
  browser()
}

testrbinom1 = function(n=1000000, size=100, prob=.02) {
  simulations = rbinom1(n, size, prob)
  
  xs = 1:size
  ps = dbinom1(xs, size, prob)
  
  ylim=range(ps)
  hist(simulations, breaks=seq(from=0.5, to=size+0.5, by=1), freq=FALSE, col="skyblue")
  lines(xs, ps, type="l", ylim=ylim)
  
  browser()
}









