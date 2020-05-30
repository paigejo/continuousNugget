# this script contains some miscellaneous, but useful functions

logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  res = exp(x)/(1+exp(x))
  res[x > 100] = 1
  res[x < -100] = 0
  res
}

# Do precomputations for computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
# latticeInfo: an object returned by makeLatGrids
precomputationsQ = function(latticeInfo, thisLayer=1) {
  
  # save layer info quantities
  nLayer = length(latticeInfo)
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  xRange=latticeInfo[[thisLayer]]$xRangeKnots
  yRange=latticeInfo[[thisLayer]]$yRangeKnots
  knotPts = latticeInfo[[thisLayer]]$latCoords
  
  # make (Laplacian) differential operators
  Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
  Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
  
  # generate x and y (Laplacian) differential operators
  Inx = Diagonal(n=nx)
  Iny = Diagonal(n=ny)
  Bx = kronecker(Iny, Dnx)
  By = kronecker(Dny, Inx)
  Bxy = Bx + By
  
  ## now construct relevant row of A for value at midpoint at this layer, and Q matrix
  
  # make mid lattice point
  # xi = ceiling(nx/2)
  # yi = ceiling(ny/2)
  # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
  midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
  
  Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
  
  # return results
  # If multiple layers, return block diagonal sparse matrix
  if(thisLayer == nLayer) {
    return(list(list(Bxy=Bxy, Ai=Ai)))
  }
  else {
    return(c(list(list(Bxy=Bxy, Ai=Ai)), precomputationsQ(latticeInfo, thisLayer+1)))
  }
}

# Do precomputations for computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers. Return pre- constructed block diagonal matrices instead of 
# individual matrices going into the block diagonal
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
# latticeInfo: an object returned by makeLatGrids
precomputationsQ2 = function(latticeInfo) {
  
  out = precomputationsQ(latticeInfo, 1)
  Bxys = lapply(out, function(x) {x$Bxy})
  Ais = lapply(out, function(x) {x$Ai})
  Bxy = bdiag(Bxys)
  mBxymBxyT = -Bxy - t(Bxy) # - Bxy - t(Bxy)
  BxyTBxy = t(Bxy) %*% Bxy # t(Bxy) %*% Bxy
  A = bdiag(Ais)
  At = t(A)
  ms = sapply(out, function(x) {length(x$Ai)})
  
  list(mBxymBxyT=mBxymBxyT, BxyTBxy=BxyTBxy, A=A, At=At, ms=ms)
}

# precompute normalization constants using natural smoothing spline on log log scale
precomputeNormalization = function(xRangeDat=c(-1,1), yRangeDat=c(-1,1), effRangeRange=NULL, nLayer=3, NC=13, 
                                   nBuffer=5, nKnots=NULL, saveResults=FALSE, doFinalTest=!saveResults, 
                                   latticeInfo=NULL, plotNormalizationSplines=TRUE) {
  # construct lattice info if necessary, precompute relevant matrices
  if(is.null(latticeInfo))
    latticeInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  else {
    xRangeDat = latticeInfo[[1]]$xRangeDat
    yRangeDat = latticeInfo[[1]]$yRangeDat
    nLayer = length(latticeInfo)
    NC = latticeInfo[[1]]$NC
    nBuffer = latticeInfo[[1]]$nBuffer
  }
  precomputationsQ = precomputationsQ2(latticeInfo)
  
  # if there are any missing layers, then we know that we are allowing for multiple kappas
  latticeWidths = sapply(latticeInfo, function(x) {x$latWidth})
  if(any(latticeWidths[2:length(latticeWidths)] != latticeWidths[1:(length(latticeWidths)-1)] / 2))
    singleKappa = FALSE
  else
    singleKappa = TRUE
  
  # set the range of effRange if necessary to be between half of the finest layer's lattice width and the data domain diameter
  if(is.null(effRangeRange))
    effRangeRange = c(latticeInfo[[nLayer]]$latWidth / 5, max(c(diff(xRangeDat), diff(yRangeDat))))
  
  # set the number of knot points so that the second to largest point is roughly 95% of the maximum knot point if necessary
  if(is.null(nKnots)) {
    width = abs(log(.95))
    nKnots = ceiling(diff(log(effRangeRange))/width) + 1
  }
  
  # set the values of effRange between which we want to interpolate
  effRangeKnots = exp(seq(log(effRangeRange[1]), log(effRangeRange[2]), l=nKnots))
  if(singleKappa) {
    kappaKnots = sqrt(8)/effRangeKnots * latticeInfo[[1]]$latWidth
  } else {
    effRangeKnots = matrix(effRangeKnots, nrow=1)
    kappaKnots = outer(sapply(latticeInfo, function(x) {x$latWidth}), c(sqrt(8)/effRangeKnots), "*")
  }
  
  # compute ctilde vector for each value of effective range. The ctilde value for one layer is independent of the 
  # kappa value in another layer. Hence, only univariate splines are necessary to precompute
  # NOTE: multiply these by alphas = c(1/nLayer, ..., 1/nLayer) now, then divide by alpha later depending on alpha
  if(singleKappa) {
    ctildes = sapply(kappaKnots, makeQPrecomputed, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                     alphas=rep(1/nLayer, nLayer), normalized=TRUE, fastNormalize=TRUE, returnctildes=TRUE) / nLayer
  } else {
    ctildes = apply(kappaKnots, 2, makeQPrecomputed, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                    alphas=rep(1/nLayer, nLayer), normalized=TRUE, fastNormalize=TRUE, returnctildes=TRUE) / nLayer
  }
  
  # estimate splines:
  # a^T %*% Q^(-1) %*% a /(rho * alpha) = ctildes
  getSplineFun = function(cts) {
    logFun = splinefun(log(effRangeKnots), log(cts), method = "hyman")
    function(x, alpha) {exp(logFun(log(x)))/alpha}
  }
  funs = apply(ctildes, 1, getSplineFun)
  fullFun = function(effRange, alphas) {
    # sapply(alphas, funs, x=effRange)
    res = numeric(nLayer)
    for(i in 1:nLayer) {
      if(length(effRange) == 1)
        res[i] = funs[[i]](effRange, alphas[i])
      else
        res[i] = funs[[i]](effRange[i], alphas[i])
    }
    res
  }
  
  # plot functions
  for(i in 1:nLayer) {
    effRanges = seq(effRangeRange[1], effRangeRange[2], l=500)
    thisctildes = funs[[i]](effRanges, alpha=1/nLayer)
    if(plotNormalizationSplines) {
      par(mfrow=c(1,1))
      plot(effRanges, thisctildes, type="l", col="blue", main=paste0("Layer ", i), xlab="Effective Range", ylab="ctilde")
      points(effRangeKnots, ctildes[i,]*nLayer, pch=19, cex=.3)
    }
  }
  
  for(i in 1:nLayer) {
    effRanges = seq(effRangeRange[1], effRangeRange[2], l=500)
    thisctildes = funs[[i]](effRanges, alpha=1/nLayer)
    if(plotNormalizationSplines) {
      par(mfrow=c(1,1))
      plot(log(effRanges), log(thisctildes), type="l", col="blue", main=paste0("Layer ", i), xlab="Log Effective Range", ylab="Log ctilde")
      points(log(effRangeKnots), log(ctildes[i,]*nLayer), pch=19, cex=.3)
    }
  }
  
  if(doFinalTest) {
    # the true value of alphas doesn't matter as long as it is different than what was used in the precomputation
    alphas = getAlphas(nLayer)
    
    if(singleKappa) {
      i = round(nKnots/2)
      thisKappa = kappaKnots[i]
      thisEffRange = effRangeKnots[i]
      thisctildes = fullFun(thisEffRange, alphas)
    } else {
      i = round(nKnots/2)
      thisKappa = kappaKnots[,i]
      thisEffRange = effRangeKnots[i]
      thisctildes = fullFun(thisEffRange, alphas)
    }
    Q = makeQPrecomputed(kappa=thisKappa, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                         alphas=alphas, normalized=TRUE, fastNormalize=TRUE)
    Q2 = makeQPrecomputed(kappa=thisKappa, precomputedMatrices=precomputationsQ, latticeInfo=latticeInfo, 
                          alphas=alphas, normalized=TRUE, fastNormalize=TRUE, ctildes=thisctildes)
    if(plotNormalizationSplines)
      print(mean(abs(Q - Q2)))
  }
  
  # save functions
  if(saveResults) {
    allNCs = sapply(latticeInfo, function(x) {x$NC})
    if(length(allNCs) == 1)
      ncText = paste0("_NC", allNCs)
    else {
      tempText = do.call("paste0", as.list(c(allNCs[1], paste0("_", allNCs[-1]))))
      ncText = paste0("_NC", tempText)
    }
    
    save(list(funs=funs, fullFun=fullFun, latInfo=latticeInfo), 
         file=paste0("ctildeSplines_nLayer", nLayer, ncText, 
                     "_xmin", round(xRangeDat[1], 1), "_xmax", round(xRangeDat[2], 1), 
                     "_ymin", round(yRangeDat[1], 1), "_ymax", round(yRangeDat[2], 1), 
                     ".RData"))
  }
  
  list(funs=funs, fullFun=fullFun)
}

# get the marginal variance for multi-resolution process
# tod: theta/delta, or theta/latticeWidth
# either nu or alphas must be non-null
getMultiMargVar = function(kappa=1, rho=1, tod=2.5, nLayer=3, nu=NULL, alphas=NULL, nx=NULL, ny=NULL, 
                           xRange=c(0,1), yRange=xRange, xRangeDat=c(-2,1), yRangeDat=xRangeDat, nBuffer=5) {
  # set alphas if nu has been set
  if(!is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  
  # set nx and ny if necessary and add buffer to avoid edge effects
  if(is.null(nx) || is.null(ny)) {
    maxPt = ceiling(tod)*4 + 1
    nx = maxPt
    ny = maxPt
  }
  
  # generate knot lattice locations and filter out locations 
  # too far outside of the data domain
  origNX = nx
  origNY = ny
  knotXs = seq(xRange[1], xRange[2], l=nx)
  knotYs = seq(yRange[1], yRange[2], l=ny)
  if(sum(knotXs > xRangeDat[2]) > nBuffer)
    knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  if(sum(knotXs < xRangeDat[1]) > nBuffer)
    knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  if(sum(knotYs > yRangeDat[2]) > nBuffer)
    knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  if(sum(knotYs < yRangeDat[1]) > nBuffer)
    knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  nx = length(knotXs)
  ny = length(knotYs)
  
  # sum the variances of each layer weighted by alphas
  totalMargVar = c()
  for(l in 1:nLayer) {
    # get the layer marginal variances
    layerMargVar = as.numeric(getMargVar(kappa, rho, tod, origNX*2^(l-1), origNY*2^(l-1), xRange=xRange, yRange=yRange, 
                                         xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))
    layerMargVar[1:2] = layerMargVar[1:2]*alphas[l]
    if(l == 1)
      totalMargVar = layerMargVar[1:2]
    else
      totalMargVar = totalMargVar + layerMargVar[1:2]
  }
  
  # add in a variance ratio column
  totalMargVar = c(totalMargVar, totalMargVar[1]/totalMargVar[2])
  names(totalMargVar) = c("actualVar", "theorVar", "inflation")
  
  totalMargVar
}

# compute the marginal variance for a given resolution layer
# tod: theta/delta, or theta/latticeWidth
getMargVar = function(kappa=1, rho=1, tod=2.5, nx=NULL, ny=NULL, xRange=c(-1,2), yRange=xRange, 
                      xRangeDat=c(0,1), yRangeDat=xRangeDat, nBuffer=5) {
  # set nx and ny if necessary and add buffer to avoid edge effects
  if(is.null(nx) || is.null(ny)) {
    maxPt = ceiling(tod)*4 + 1
    nx = maxPt
    ny = maxPt
  }
  
  # generate knot lattice locations and filter out locations 
  # too far outside of the data domain
  knotXs = seq(xRange[1], xRange[2], l=nx)
  knotYs = seq(yRange[1], yRange[2], l=ny)
  delta = knotXs[2]-knotXs[1]
  if(sum(knotXs > xRangeDat[2]) > nBuffer)
    knotXs = knotXs[1:(length(knotXs) - (sum(knotXs > xRangeDat[2]) - nBuffer))]
  if(sum(knotXs < xRangeDat[1]) > nBuffer)
    knotXs = knotXs[(1 + sum(knotXs < xRangeDat[1]) - nBuffer):length(knotXs)]
  if(sum(knotYs > yRangeDat[2]) > nBuffer)
    knotYs = knotYs[1:(length(knotYs) - (sum(knotYs > yRangeDat[2]) - nBuffer))]
  if(sum(knotYs < yRangeDat[1]) > nBuffer)
    knotYs = knotYs[(1 + sum(knotYs < yRangeDat[1]) - nBuffer):length(knotYs)]
  knotPts = make.surface.grid(list(x=knotXs, y=knotYs))
  
  # take the middle knot location
  midPt = matrix(c(knotXs[ceiling(length(knotXs)/2)], 
                   knotYs[ceiling(length(knotYs)/2)]), nrow=1)
  
  # compute variance of process at the middle knot
  A = as.matrix(makeA(midPt, xRange, nx, yRange, ny, tod*delta, 
                      xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer))
  Q = makeQ(kappa, rho, xRange, yRange, nx, ny, 
            xRangeDat=xRangeDat, yRangeDat=yRangeDat, nBuffer=nBuffer)
  # VarC = as.matrix(solve(Q))
  varMidPt = A %*% inla.qsolve(Q, t(A))
  
  # compare with theoretical marginal variance
  sigma2 = rho/(4*pi * kappa^2)
  inflation = varMidPt/sigma2
  
  # return results
  list(actualVar=varMidPt, theorVar=sigma2, inflation=inflation)
}


# estimates effective range for a given LatticeKrig model
getEffRange = function(predPts=NULL, xRangeKnot=c(0,1), xNKnot=10, yRangeKnot=c(0,1), yNKnot=10, theta=NULL, nLayer=1, thisLayer=1, 
                       xRangeDat=xRangeKnot, yRangeDat=yRangeKnot, nBuffer=5, mx=20, my=20) {
  
}

##### simulate from a latticeKrig model
# return a function with argument nsim for simulating a number of realizations from a latticeKrig model with no nugget.  
# coords: coordinates at which to simulate
# all other arguments: same as general latticeKrig arguments
LKSimulator = function(coords, NC=5, kappa=1, rho=1, nu=1.5, nBuffer=5, nLayer=3, normalize=TRUE) {
  # first make the grid on which to set the basis functions
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  knotGrid = makeLatGrid(xRange=xRangeDat, yRange=yRangeDat, NC=NC, nBuffer=nBuffer)
  xRangeKnots=knotGrid$xRangeKnots
  nx=knotGrid$nx
  yRangeKnots=knotGrid$yRangeKnots
  ny=knotGrid$ny
  
  # generate layer variance weights
  alphas = getAlphas(nLayer, nu)
  
  # generate basis function and precision matrices
  AObs = makeA(coords, xRangeKnots, nx, yRangeKnots, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeBasis, yRange=yRangeBasis, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize)
  L = as.matrix(t(chol(solve(Q))))
  zsim = matrix(rnorm(nrow(Q)), ncol=1)
  fieldSim = L %*% zsim
  
  fieldSims
}

# return a function with argument nsim for simulating a number of realizations from a latticeKrig model with no nugget.  
# same as LKSimulator, but uses marginal variance, effective range parameterization instead of rho, kappa.
# coords: coordinates at which to simulate
# all other arguments: same as general latticeKrig arguments
LKSimulator2 = function(coords, nsim=1, NC=5, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1, 
                        nu=1.5, nBuffer=5, nLayer=3, normalize=TRUE) {
  
  # first make the grid on which to set the basis functions
  xRangeDat = range(coords[,1])
  yRangeDat = range(coords[,2])
  knotGrid = makeLatGrid(xRange=xRangeDat, yRange=yRangeDat, NC=NC, nBuffer=nBuffer)
  xRangeKnots=knotGrid$xRangeKnots
  nx=knotGrid$nx
  yRangeKnots=knotGrid$yRangeKnots
  ny=knotGrid$ny
  
  # convert from effRange, margVar to rho, kappa
  latticeWidth = (xRangeKnots[2] - xRangeKnots[1])/(nx-1)
  kappa = sqrt(8)/effRange * latticeWidth
  
  # since we are normalizing the process, rho is just sigmaSq
  rho = margVar
  
  # generate layer variance weights
  alphas = getAlphas(nLayer, nu)
  
  # generate basis function and precision matrices
  AObs = makeA(coords, xRangeKnots, nx, yRangeKnots, ny, nLayer=nLayer, xRangeDat=xRangeDat, 
               yRangeDat=yRangeDat, nBuffer=nBuffer)
  Q = makeQ(kappa=kappa, rho=rho, xRange=xRangeKnots, yRange=yRangeKnots, nx=nx, ny=ny, 
            nLayer=nLayer, alphas=alphas, xRangeDat=xRangeDat, yRangeDat=yRangeDat, 
            nBuffer=nBuffer, normalized = normalize)
  L = as.matrix(t(chol(solve(Q))))
  zsims = matrix(rnorm(nrow(Q)*nsim), ncol=nsim)
  fieldSims = matrix(as.numeric(AObs %*% L %*% zsims), ncol=nsim)
  
  fieldSims
}

simSPDE = function(coords, nsim=1, mesh=NULL, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1, kenya=FALSE) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    if(kenya)
      mesh = getSPDEMeshKenya(coords, doPlot = FALSE)
    else
      mesh = getSPDEMesh(doPlot = FALSE)
  }
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate A and Q precision matrix
  Q = inla.spde2.precision(spde, theta = theta)
  A = inla.spde.make.A(mesh, coords)
  
  # generate simulations
  simField = inla.qsample(nsim, Q)
  simDat = as.matrix(A %*% simField)
  
  simDat
}

makeQSPDE = function(mesh, effRange, margVar=1) {
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate Q precision matrix
  inla.spde2.precision(spde, theta = theta)
}

# computing precision matrix for a single layer or a block diagonal sparse 
# precision matrix for multiple layers
# kappa: scale of Matern covariance with smoothness 1
# xRange, yRange: x and y intervals in space over which basis elements are placed
# nx, ny: number of basis elements when counting along the lattive in 
#         x and y directions respectively (for the first layer)
# rho: sill (marginal variance) for first layer
# nLayer: number of lattice layers
# thisLayer: user should always set to 1 to get full block diagonal precision matrix 
#            for all layers
# alphas: weights on the variances of each layer.  Scales rho for each layer.
# fastNormalize: simple normalization to make marginal variance = rho in center. Basis coefficients may have different variances
# assume there is a single kappa, rho is adjusted with alpha when there are multiple layers
# latticeInfo: an object returned by makeLatGrids
makeQ = function(kappa=1, rho=1, latticeInfo, thisLayer=1, alphas=NULL, nu=NULL, 
                 normalized=FALSE, fastNormalize=FALSE, precomputedMatrices=NULL) {
  require(Matrix)
  require(spam)
  require(fields)
  
  # save base layer input quantities
  origRho = rho
  nLayer = length(latticeInfo)
  
  # make alphas according to nu relation if nu is set
  if(is.null(alphas) && !is.null(nu)) {
    alphas = getAlphas(nLayer, nu)
  }
  else if(is.null(nu) && nLayer != 1 && is.null(alphas)) {
    warning("Both alphas and nu are NULL. Defaulting to exponential covariance.")
    nu = 0.5
    alphas = getAlphas(nLayer, nu)
  }
  
  nx = latticeInfo[[thisLayer]]$nx
  ny = latticeInfo[[thisLayer]]$ny
  if(is.null(precomputedMatrices)) {
    # make (Laplacian) differential operators
    Dnx = bandSparse(nx, k=0:1, diag=list(rep(-2, nx), rep(1, nx-1)), symmetric=TRUE)
    Dny = bandSparse(ny, k=0:1, diag=list(rep(-2, ny), rep(1, ny-1)), symmetric=TRUE)
    
    # generate x and y (Laplacian) differential operators
    Inx = Diagonal(n=nx)
    Iny = Diagonal(n=ny)
    Bx = kronecker(Iny, Dnx)
    By = kronecker(Dny, Inx)
    Bxy = Bx + By
  }
  else
    Bxy = precomputedMatrices[[thisLayer]]$Bxy
  
  # make B, SAR regression matrix for Bc = e
  B = Diagonal(n=nx*ny, x=kappa^2) - Bxy
  
  # compute precision matrix
  Q = t(B) %*% B
  
  if(normalized) {
    
    # now construct relevant row of A for value at midpoint at this layer, and Q matrix
    if(is.null(precomputedMatrices)) {
      xRange=latticeInfo[[thisLayer]]$xRangeKnots
      yRange=latticeInfo[[thisLayer]]$yRangeKnots
      
      # make mid lattice point
      # xi = ceiling(nx/2)
      # yi = ceiling(ny/2)
      # midPt = matrix(c(seq(xRange[1], xRange[2], l=nx)[xi], seq(yRange[1], yRange[2], l=ny)[yi]), nrow=1)
      midPt = matrix(c((xRange[1] + xRange[2])/2, (yRange[1] + yRange[2])/2), nrow=1)
      
      Ai = makeA(midPt, latticeInfo, thisLayer=thisLayer, maxLayer=thisLayer)
    }
    else
      Ai = precomputedMatrices[[thisLayer]]$Ai
    
    # # test
    # sds2 = 1/diag(Q)
    # sdMat = Diagonal(x=sds2)
    # # Qnorm2 = sweep(sweep(Q, 1, sds2, "*"), 2, sds2, "*")
    # Qnorm2 = sdMat %*% Q %*% sdMat
    # 
    # # QnormInv = sweep(sweep(Qinv, 1, 1/sds), 2, 1/sds)
    # procVar2 = as.numeric(Ai %*% inla.qsolve(Qnorm2, t(Ai)))
    # # procVar = as.numeric(Ai %*% QnormInv %*% t(Ai))
    # Q2 = Qnorm2 * (procVar2 / rho)
    # # Q = Q2 # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): ~5.8s
    
    #  test 2
    if(fastNormalize) {
      ctilde = as.numeric(Ai %*% inla.qsolve(Q, t(Ai)))
      Qtilde = ctilde * Q
      Q = Qtilde # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE, newnormalize=TRUE)): 2.72
    }
    else {
      # renormalize basis coefficients to have constant variance, and the process to have unit variance
      Qinv = inla.qsolve(Q, diag(nrow(Q)))
      sds = sqrt(diag(Qinv))
      sdMat = Diagonal(x=sds)
      # Qnorm = sweep(sweep(Q, 1, sds, "*"), 2, sds, "*")
      Qnorm = sdMat %*% Q %*% sdMat
      # QnormInv = sweep(sweep(Qinv, 1, 1/sds), 2, 1/sds)
      procVar = as.numeric(Ai %*% inla.qsolve(Qnorm, t(Ai)))
      # procVar = as.numeric(Ai %*% QnormInv %*% t(Ai))
      Q = Qnorm * procVar # system.time(out <- makeQ(nLayer=3, nx=15, ny=15, nu=1, normalized=TRUE) ~5.87
    }
    # # compute how good an approximation it was
    # hist(diag(solve(Q)))
    # hist(diag(solve(Q2)))
    # hist(diag(solve(Qtilde)))
    # image(Q)
    # image(Q2)
    # image(Qtilde)
  }
  if(nLayer == 1 && is.null(alphas))
    alphas = 1
  Q = Q * (1/alphas[thisLayer])
  
  # return results
  # If multiple layers, return block diagonal sparse matrix
  if(thisLayer == nLayer) {
    if(thisLayer == 1)
      Q = Q * (1/rho)
    return(Q)
  }
  else if((thisLayer == 1) && (nLayer != 1)) {
    Q = bdiag(c(list(Q), makeQ(kappa, origRho, latticeInfo, thisLayer+1, alphas, normalized=normalized, 
                               fastNormalize=fastNormalize)))
    Q = Q * (1/rho)
    return(Q)
  }
  else {
    return(c(list(Q), makeQ(kappa, origRho, latticeInfo, thisLayer+1, alphas, normalized=normalized, 
                            fastNormalize=fastNormalize)))
  }
}

meanSegmentLength = function(mesh, filterLargerThan=NULL) {
  t.sub = 1:nrow(mesh$graph$tv)
  idx = cbind(mesh$graph$tv[t.sub, c(1:3, 1), drop = FALSE], NA)
  x = mesh$loc[t(idx), 1]
  y = mesh$loc[t(idx), 2]
  indices = 1:4 + rep(seq(from=0, to=length(x)-5, by=5), each=4)
  segx = x[indices]
  segy = y[indices]
  coords = cbind(segx, segy)
  dists = rdist.vec(coords[1:(length(segx) - 1),], coords[2:length(segx),])
  dists = dists[-seq(from=4, to=length(dists), by=4)]
  
  if(!is.null(filterLargerThan))
    dists = dists[dists <= filterLargerThan]
  
  mean(dists)
}
# meanSegmentLength(getSPDEMeshKenya(), 50)
# meanSegmentLength(getSPDEMesh(), 0.01+1e-6)

LKINLA.cov = function(x1, x2, latticeInfo, kappa, alphas, rho=1, normalize=TRUE, 
                      fastNormalize=TRUE, 
                      precomputedMatrices=NULL, precomputedA1=NULL, precomputedA2=NULL, 
                      precomputationsFileNameRoot="") {
  ctildes = NULL
  
  if(precomputationsFileNameRoot != "") {
    # load in precomputations
    load(paste0("savedOutput/precomputations/", precomputationsFileNameRoot, ".RData"))
    
    if(!is.null(precomputedNormalizationFun)) {
      if(length(kappa) == 1)
        effectiveCor = sqrt(8) / kappa * latticeInfo[[1]]$latWidth
      else
        effectiveCor = sqrt(8) / kappa * sapply(latticeInfo, function(x) {x$latWidth})
      ctildes = precomputedNormalizationFun$fullFun(effectiveCor, alphas)
    }
    
  }
  if(!is.null(precomputedMatrices)) {
    Q = makeQPrecomputed(precomputedMatrices, kappa, rho, latticeInfo, alphas, normalized=normalize, 
                         fastNormalize=fastNormalize, ctildes=ctildes)
  } else {
    Q = makeQ(kappa, rho, latticeInfo, 1, alphas, normalized=normalize, fastNormalize=fastNormalize)
  }
  
  if(is.null(precomputedA1))
    A1 = makeA(x1, latticeInfo, 1)
  else
    A1 = precomputedA1
  if(is.null(precomputedA2))
    A2 = makeA(x2, latticeInfo, 1)
  else
    A2 = precomputedA2
  
  if(any(alphas == 0)) {
    # in this case, we must be careful to only consider the layers with nonzero weights
    includeI = rep(TRUE, nrow(Q))
    thisI = 1
    for(i in 1:length(latticeInfo)) {
      startI = thisI
      endI = startI - 1 + latticeInfo[[i]]$nx * latticeInfo[[i]]$ny
      if(alphas[i] == 0) {
        includeI[startI:endI] = FALSE
      }
      thisI = endI + 1
    }
    A1 = matrix(as.matrix(A1)[,includeI], nrow=nrow(A1))
    A2 = matrix(as.matrix(A2)[,includeI], nrow=nrow(A2))
    Q = Q[includeI,includeI]
  }
  # use A1 %*% inla.qsolve(Q, t(A2)) instead?
  A1 %*% inla.qsolve(Q, t(A2))
}

LK.cov = function(x1, x2, LKinfo, a.wght, alphas, lambda, sigma, rho) {
  # domainCoords = LKinfo$latticeInfo$grid.info$range
  # nBuffer = LKinfo$latticeInfo$NC.buffer
  # NC = max(LKinfo$latticeInfo$mxDomain[1,])
  # LKrigSetup(domainCoords, nlevel=length(alphas), a.wght=a.wght, normalize=normalize, 
  #            lambda=lambda, sigma=sigma, rho=rho)
  if(length(a.wght) > 1)
    a.wght = as.list(a.wght)
  LKinfo = LKinfoUpdate(LKinfo, a.wght=a.wght, alpha=alphas, lambda=lambda, sigma=sigma, rho=rho)
  LKrig.cov(x1, x2, LKinfo)
}

# make method for calculating individual covariance function
getLKInlaCovarianceFun = function(kappa, rho, nuggetVar, alphas, NP=200, latticeInfo, normalize=TRUE, fastNormalize=TRUE, 
                                  precomputedMatrices=NULL, precomputedAcenter=NULL, precomputedAx=NULL, precomputedAy=NULL) {
  # generate test locations based on code from LKrig.cov.plot
  xlim <- latticeInfo[[1]]$xRangeDat
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latticeInfo[[1]]$yRangeDat
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  # calculate covariances
  x1 <- cbind(ux, rep(center[2], NP))
  x2 <- rbind(center)
  d <- c(rdist(x1, x2))
  y <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                             precomputedMatrices, precomputedAx, precomputedAcenter))
  y[NP/2] = y[NP/2] + nuggetVar
  x1 <- cbind(rep(center[1], NP), uy)
  d2 <- c(rdist(x1, x2))
  y2 <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                              precomputedMatrices, precomputedAy, precomputedAcenter))
  y2[NP/2] = y2[NP/2] + nuggetVar
  
  # average x and y covariances
  sortXI = sort(d, index.return=TRUE)$ix
  d = d[sortXI]
  y = y[sortXI]
  sortYI = sort(d2, index.return=TRUE)$ix
  d2 = d2[sortYI]
  y2 = y2[sortYI]
  d = rowMeans(cbind(d, d2))
  y = rowMeans(cbind(y, y2))
  return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
}

covarianceDistributionLKINLA = function(latticeInfo, kappaVals, rhoVals=rep(1, length(kappaVals)), nuggetVarVals=rep(0, length(kappaVals)), 
                                        alphaMat, maxSamples=100, significanceCI=.8, normalize=TRUE, fastNormalize=TRUE, seed=NULL, NP = 200, 
                                        precomputationsFileNameRoot="", maxRadius=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  nLayer = length(latticeInfo)
  
  # get hyperparameter samples
  sampleI = sample(1:length(rhoVals), min(maxSamples, length(rhoVals)))
  if(!is.null(dim(kappaVals))) {
    kappaVals = kappaVals[,sampleI]
    separateRanges = TRUE
    effectiveRanges = sweep(sqrt(8)/kappaVals, 1, sapply(latticeInfo, function(x){x$latWidth}), "*")
    minRange = min(apply(effectiveRanges, 1, min))
    maxRange = max(apply(cbind(5 * sapply(latticeInfo, function(x){x$latWidth}), effectiveRanges), 1, max))
  } else {
    kappaVals = kappaVals[sampleI]
    separateRanges = FALSE
    effectiveRanges = sqrt(8) * latticeInfo[[1]]$latWidth / kappaVals
    minRange = min(effectiveRanges) / 2^(length(latticeInfo) - 1)
    maxRange = max(c(effectiveRanges, latticeInfo[[1]]$latWidth * 5))
  }
  rhoVals = rhoVals[sampleI]
  nuggetVarVals = nuggetVarVals[sampleI]
  alphaMat = matrix(alphaMat[,sampleI], ncol=length(sampleI))
  
  # generate test locations based on code from LKrig.cov.plot (modify sampling points to be the correct resolution)
  # xlim <- latticeInfo[[1]]$xRangeDat
  # ux <- seq(xlim[1], xlim[2], , NP)
  # ylim <- latticeInfo[[1]]$yRangeDat
  # uy <- seq(ylim[1], ylim[2], , NP)
  # center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  if(is.null(maxRadius))
    maxRadius = maxRange * 2
  minStep = minRange / 10
  ThisNP = 2 * maxRadius / minStep
  xlim <- latticeInfo[[1]]$xRangeDat
  ylim <- latticeInfo[[1]]$yRangeDat
  # centerX = mean(xlim)
  # widthX = xlim[2] - centerX
  # deltaX = min(widthX, maxRadius)
  # xlim = c(centerX - deltaX, centerX + deltaX)
  # ux <- seq(xlim[1], xlim[2], , NP)
  # ylim <- latticeInfo[[1]]$yRangeDat
  # centerY = mean(ylim)
  # widthY = ylim[2] - centerY
  # deltaY = min(widthY, maxRadius)
  # ylim = c(centerY - deltaY, centerY + deltaY)
  # uy <- seq(ylim[1], ylim[2], , NP)
  # center <- rbind(c(ux[NP/2], uy[NP/2]))
  centerX = mean(xlim)
  widthX = xlim[2] - centerX
  deltaX = min(widthX, maxRadius)
  centerY = mean(ylim)
  widthY = ylim[2] - centerY
  deltaY = min(widthY, maxRadius)
  delta = min(deltaX, deltaY)
  xlim = c(centerX - delta, centerX + delta)
  # ux <- seq(xlim[1], xlim[2], , NP)
  ux <- c(seq(xlim[1], centerX-minRange/100, l=NP/2), centerX, seq(centerX+minRange/100, xlim[2], l=NP/2))
  ylim = c(centerY - delta, centerY + delta)
  # uy <- seq(ylim[1], ylim[2], , NP)
  uy <- c(seq(ylim[1], centerY-minRange/100, l=NP/2), centerY, seq(centerY+minRange/100, ylim[2], l=NP/2))
  center <- rbind(c(centerX, centerY))
  
  # precompute relevant matrices
  Qprecomputations = precomputationsQ2(latticeInfo)
  Acenter = makeA(center, latticeInfo)
  Ax = makeA(cbind(ux, center[2]), latticeInfo)
  Ay = makeA(cbind(center[1], uy), latticeInfo)
  
  # make method for calculating individual covariance function
  getOneCovariance = function(parameters) {
    # get relevant parameters
    if(!separateRanges) {
      kappa = parameters[1]
      rho = parameters[2]
      nuggetVar = parameters[3]
      alphas = parameters[-(1:3)]
    } else {
      kappa = parameters[1:nLayer]
      rho = parameters[1 + nLayer]
      nuggetVar = parameters[2 + nLayer]
      alphas = parameters[-(1:(2 + nLayer))]
    }
    
    # calculate covariances
    x1 <- cbind(ux, rep(center[2], NP+1))
    x2 <- rbind(center)
    d <- c(rdist(x1, x2))
    y <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                               Qprecomputations, Ax, Acenter, precomputationsFileNameRoot=precomputationsFileNameRoot))
    y[NP/2+1] = y[NP/2+1] + nuggetVar
    x1 <- cbind(rep(center[1], NP+1), uy)
    d2 <- c(rdist(x1, x2))
    y2 <- as.numeric(LKINLA.cov(x1, x2, latticeInfo, kappa, alphas, rho, normalize, fastNormalize, 
                                Qprecomputations, Ay, Acenter, precomputationsFileNameRoot=precomputationsFileNameRoot))
    y2[NP/2+1] = y2[NP/2+1] + nuggetVar
    
    # average x and y covariances
    # sortXI = sort(d, index.return=TRUE)$ix
    # d = d[sortXI]
    # y = y[sortXI]
    # sortYI = sort(d2, index.return=TRUE)$ix
    # d2 = d2[sortYI]
    # y2 = y2[sortYI]
    # d = rowMeans(cbind(d, d2))
    # y = rowMeans(cbind(y, y2))
    d = c(0, rowMeans(cbind(rev(d[1:(NP/2)]),  d[(NP/2+2):length(d)],  rev(d2[1:(NP/2)]),  d2[(NP/2+2):length(d2)])))
    y = c(mean(y[NP/2+1], y2[NP/2+1]), rowMeans(cbind(rev(y[1:(NP/2)]),  y[(NP/2+2):length(y)],  rev(y2[1:(NP/2)]),  y2[(NP/2+2):length(y2)])))
    sortXI = sort(d, index.return=TRUE)$ix
    d = d[sortXI]
    y = y[sortXI]
    return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
  }
  
  # calculate covariances for each sample from the posterior
  if(!separateRanges)
    parameterMat = cbind(kappaVals, rhoVals, nuggetVarVals, t(alphaMat))
  else
    parameterMat = cbind(t(kappaVals), rhoVals, nuggetVarVals, t(alphaMat))
  # browser()
  out = apply(parameterMat, 1, getOneCovariance)
  d = out[1:(NP/2+1),1]
  covMat = out[(NP/2+2):(2*(NP/2+1)),]
  corMat = out[(2*(NP/2+1)+1):(3*(NP/2+1)),]
  
  # calculate summary statistics
  meanCov = rowMeans(covMat)
  lowerCov = apply(covMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCov = apply(covMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCor = rowMeans(corMat)
  lowerCor = apply(corMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCor = apply(corMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  
  # return results
  list(d=d, 
       cov=meanCov, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
       cor=meanCor, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
}

covarianceDistributionLK = function(latticeInfo, alphaVals, lambdaVals, a.wghtVals, rhoVals, 
                                    maxSamples=100, significanceCI=.8, normalize=TRUE, seed=NULL) {
  NP = 200
  if(!is.null(seed))
    set.seed(seed)
  
  # get number of layers
  nLayer = latticeInfo$nlevel
  a.wghtVals = matrix(a.wghtVals, ncol=length(lambdaVals))
  
  # get hyperparameter samples
  nuggetVarVals = rhoVals * lambdaVals
  sampleI = sample(1:length(rhoVals), maxSamples)
  alphaMat = alphaVals[,sampleI]
  lambdaVals = lambdaVals[sampleI]
  a.wghtVals = matrix(a.wghtVals[,sampleI], ncol=maxSamples)
  nuggetVarVals = nuggetVarVals[sampleI]
  rhoVals = rhoVals[sampleI]
  
  # generate test locations based on code from LKrig.cov.plot
  xlim <- latticeInfo$latticeInfo$rangeLocations[,1]
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latticeInfo$latticeInfo$rangeLocations[,2]
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  # make method for calculating individual covariance function
  getOneCovariance = function(parameters) {
    # get relevant parameters
    a.wght = parameters[1:nrow(a.wghtVals)]
    rho = parameters[nrow(a.wghtVals) + 1]
    nuggetVar = parameters[nrow(a.wghtVals) + 2]
    lambda = parameters[nrow(a.wghtVals) + 3]
    alphas = parameters[-(1:(nrow(a.wghtVals) + 3))]
    
    # calculate covariances
    x1 <- cbind(ux, rep(center[2], NP))
    x2 <- rbind(center)
    d <- c(rdist(x1, x2))
    y <- as.numeric(LK.cov(x1, x2, latticeInfo, a.wght, alphas, lambda, sqrt(nuggetVar), rho))
    y[NP/2] = y[NP/2] + nuggetVar
    x1 <- cbind(rep(center[1], NP), uy)
    d2 <- c(rdist(x1, x2))
    y2 <- as.numeric(LK.cov(x1, x2, latticeInfo, a.wght, alphas, lambda, sqrt(nuggetVar), rho))
    y2[NP/2] = y2[NP/2] + nuggetVar
    
    # calculate covariances excluding nugget
    y3 = y
    y3[NP/2] = y[NP/2] - nuggetVar
    y4 = y2
    y4[NP/2] = y2[NP/2] - nuggetVar
    
    # average x and y covariances
    sortXI = sort(d, index.return=TRUE)$ix
    d = d[sortXI]
    y = y[sortXI]
    y3 = y3[sortXI]
    sortYI = sort(d2, index.return=TRUE)$ix
    d2 = d2[sortYI]
    y2 = y2[sortYI]
    y4 = y4[sortYI]
    d = rowMeans(cbind(d, d2))
    y = rowMeans(cbind(y, y2))
    yNoNugget = rowMeans(cbind(y3, y4))
    return(cbind(d=d, cov=y, cor=y * (1 / max(y)), covNoNugget=yNoNugget, norNoNugget=yNoNugget * (1 / max(yNoNugget))))
  }
  
  # calculate covariances for each sample from the posterior
  parameterMat = cbind(t(a.wghtVals), rhoVals, nuggetVarVals, lambdaVals, t(alphaMat))
  # browser()
  out = apply(parameterMat, 1, getOneCovariance)
  d = out[1:200,1]
  covMat = out[201:400,]
  corMat = out[401:600,]
  covMatNoNugget = out[601:800,]
  corMatNoNugget = out[801:1000,]
  
  # calculate summary statistics
  meanCov = rowMeans(covMat)
  lowerCov = apply(covMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCov = apply(covMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCor = rowMeans(corMat)
  lowerCor = apply(corMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCor = apply(corMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCovNoNugget = rowMeans(covMatNoNugget)
  lowerCovNoNugget = apply(covMatNoNugget, 1, quantile, probs=(1-significanceCI)/2)
  upperCovNoNugget = apply(covMatNoNugget, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCorNoNugget = rowMeans(corMatNoNugget)
  lowerCorNoNugget = apply(corMatNoNugget, 1, quantile, probs=(1-significanceCI)/2)
  upperCorNoNugget = apply(corMatNoNugget, 1, quantile, probs=1 - (1-significanceCI)/2)
  
  # return results
  list(d=d, 
       cov=meanCov, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
       cor=meanCor, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat, 
       covNoNugget=meanCovNoNugget, upperCovNoNugget=upperCovNoNugget, lowerCovNoNugget=lowerCovNoNugget, covMatNoNugget=covMatNoNugget, 
       corNoNugget=meanCorNoNugget, upperCorNoNugget=upperCorNoNugget, lowerCorNoNugget=lowerCorNoNugget, corMatNoNugget=corMatNoNugget)
}

covarianceDistributionSPDE = function(effectiveRangeVals, rhoVals=rep(1, length(effectiveRangeVals)), nuggetVarVals=rep(0, length(rhoVals)), 
                                      mesh, maxSamples=100, significanceCI=c(.8, .95), seed=NULL, xRangeDat=NULL, yRangeDat=NULL, NP = 200, 
                                      maxRadius=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  # get hyperparameter samples
  sampleI = sample(1:length(effectiveRangeVals), maxSamples)
  effectiveRangeVals = effectiveRangeVals[sampleI]
  rhoVals = rhoVals[sampleI]
  nuggetVarVals = nuggetVarVals[sampleI]
  
  # generate test locations based on code from LKrig.cov.plot
  if(is.null(xRangeDat) || is.null(yRangeDat)) {
    idx = unique(c(mesh$segm$int$idx[,1], mesh$segm$int$idx[,2]))
    locs = mesh$loc[idx,]
    xlim = range(locs[,1])
    ylim = range(locs[,2])
  } else {
    xlim = xRangeDat
    ylim = yRangeDat
  }
  # ux <- seq(xlim[1], xlim[2], , NP)
  # uy <- seq(ylim[1], ylim[2], , NP)
  # center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  maxRange = max(effectiveRangeVals)
  minRange = min(effectiveRangeVals)
  if(is.null(maxRadius))
    maxRadius = maxRange * 2
  minStep = minRange / 10
  ThisNP = 2 * maxRadius / minStep
  centerX = mean(xlim)
  widthX = xlim[2] - centerX
  deltaX = min(widthX, maxRadius)
  centerY = mean(ylim)
  widthY = ylim[2] - centerY
  deltaY = min(widthY, maxRadius)
  delta = min(deltaX, deltaY)
  xlim = c(centerX - delta, centerX + delta)
  # ux <- seq(xlim[1], xlim[2], , NP)
  ux <- c(seq(xlim[1], centerX-minRange/100, l=NP/2), centerX, seq(centerX+minRange/100, xlim[2], l=NP/2))
  ylim = c(centerY - delta, centerY + delta)
  # uy <- seq(ylim[1], ylim[2], , NP)
  uy <- c(seq(ylim[1], centerY-minRange/100, l=NP/2), centerY, seq(centerY+minRange/100, ylim[2], l=NP/2))
  center <- rbind(c(centerX, centerY))
  
  # precompute relevant matrices
  Acenter = inla.spde.make.A(mesh, center)
  Ax = inla.spde.make.A(mesh, cbind(ux, center[2]))
  Ay = inla.spde.make.A(mesh, cbind(center[1], uy))
  
  # make method for calculating individual covariance function
  getOneCovariance = function(parameters) {
    # get relevant parameters
    effectiveRange = parameters[1]
    rho = parameters[2]
    nuggetVar = parameters[3]
    i = parameters[4]
    print(paste0("Calculating covariances for iteration ", i, "/", length(rhoVals)))
    
    # calculate covariances
    x1 <- cbind(ux, rep(center[2], NP+1))
    x2 <- rbind(center)
    d <- c(rdist(x1, x2))
    Q = makeQSPDE(mesh, effectiveRange, rho)
    y = as.numeric(Acenter %*% inla.qsolve(Q, t(Ax)))
    y[NP/2+1] = y[NP/2+1] + nuggetVar
    x1 <- cbind(rep(center[1], NP+1), uy)
    d2 <- c(rdist(x1, x2))
    y2 = as.numeric(Acenter %*% inla.qsolve(Q, t(Ay)))
    y2[NP/2+1] = y2[NP/2+1] + nuggetVar
    
    # average x and y covariances
    # sortXI = sort(d, index.return=TRUE)$ix
    # d = d[sortXI]
    # y = y[sortXI]
    # sortYI = sort(d2, index.return=TRUE)$ix
    # d2 = d2[sortYI]
    # y2 = y2[sortYI]
    # d = rowMeans(cbind(d, d2))
    # y = rowMeans(cbind(y, y2))
    # spatialCov = y[NP/2+1]
    # return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
    d = c(0, rowMeans(cbind(rev(d[1:(NP/2)]),  d[(NP/2+2):length(d)],  rev(d2[1:(NP/2)]),  d2[(NP/2+2):length(d2)])))
    y = c(mean(y[NP/2+1], y2[NP/2+1]), rowMeans(cbind(rev(y[1:(NP/2)]),  y[(NP/2+2):length(y)],  rev(y2[1:(NP/2)]),  y2[(NP/2+2):length(y2)])))
    sortXI = sort(d, index.return=TRUE)$ix
    d = d[sortXI]
    y = y[sortXI]
    return(cbind(d=d, cov=y, cor=y * (1 / max(y))))
  }
  
  # calculate covariances for each sample from the posterior
  parameterMat = cbind(effectiveRangeVals, rhoVals, nuggetVarVals, 1:length(rhoVals))
  # browser()
  out = apply(parameterMat, 1, getOneCovariance)
  d = out[1:(NP/2+1),1]
  covMat = out[(NP/2+2):(2*(NP/2+1)),]
  corMat = out[(2*(NP/2+1)+1):(3*(NP/2+1)),]
  
  # calculate summary statistics
  meanCov = rowMeans(covMat)
  lowerCov = apply(covMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCov = apply(covMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  meanCor = rowMeans(corMat)
  lowerCor = apply(corMat, 1, quantile, probs=(1-significanceCI)/2)
  upperCor = apply(corMat, 1, quantile, probs=1 - (1-significanceCI)/2)
  
  # return results
  list(d=d, 
       cov=meanCov, upperCov=upperCov, lowerCov=lowerCov, covMat=covMat, 
       cor=meanCor, upperCor=upperCor, lowerCor=lowerCor, corMat=corMat)
}

# test how close we can get to the spatial correlation function:
getTrueLKEffectiveRange = function(nLayer=3, NP=200, sigma2 = 0, rho=1, 
                                   nBuffer=5, normalize=TRUE, fastNormalize=TRUE, NC=13, 
                                   latInfo=NULL, effectiveRange=1, alphas=rep(1/nLayer, nLayer-1)) {
  
  
  # construct the lattice
  if(is.null(latInfo)) {
    xRangeDat = c(-1, 1)
    yRangeDat = c(-1, 1)
    latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer)
  }
  # set true parameter values
  kappa = (sqrt(8) * latInfo[[1]]$latWidth /effectiveRange)
  
  xlim <- latInfo[[1]]$xRangeDat
  ux <- seq(xlim[1], xlim[2], , NP)
  ylim <- latInfo[[1]]$yRangeDat
  uy <- seq(ylim[1], ylim[2], , NP)
  center <- rbind(c(ux[NP/2], uy[NP/2]))
  
  test = getLKInlaCovarianceFun(kappa, rho, sigma2, alphas, latticeInfo = latInfo, 
                                normalize=normalize, fastNormalize=fastNormalize)
  ds = test[,1]
  firstI = match(TRUE, test[,3] <= .1)
  if(is.na(firstI)) {
    warning(paste0("effective range larger than domain diameter, returned value is too small.  Correlation is ", 
                   test[nrow(test),3], " at distance ", ds[length(ds)]))
    firstI = length(ds)
  }
  ds[firstI]
}

# calculate the minimum value of NC in order to have at least a certain number of basis functions total 
# for a fixed value of nLayer and nBuffer
getMinimalLattice = function(nMin=5000, locs=cbind(c(-1, 1), c(-1, 1)), nLayer=3, nBuffer=5) {
  # set starting NC so that the minimal dimension has at least two basis functions
  xRange = range(locs[,1])
  yRange = range(locs[,2])
  minRange = min(c(diff(xRange), diff(yRange)))
  maxRange = max(c(diff(xRange), diff(yRange)))
  NC = ceiling(2 * maxRange / minRange)
  
  # keep increasing NC until there are enough basis functions
  N = 0
  while(N < nMin) {
    latticeInfo = makeLatGrids(xRange, yRange, NC, nBuffer, nLayer)
    N = sum(sapply(latticeInfo, function(x) {x$ny * x$nx}))
    print(paste0(N, " basis functions for NC=", NC))
    NC = NC + 1
  }
  NC = NC - 1
  
  # return results
  c(NC=NC, nbasis=N)
}

rpcvar = function(n, alpha=.01, u=1) {
  1 / inla.pc.rprec(n, alpha=alpha, u=u)
}

dpcvar = function(x, alpha=.01, u=1) {
  inla.pc.dprec(1 / x, alpha=alpha, u=u) / x^2
}

qpcvar = function(p, alpha=.01, u=1) {
  1 / inla.pc.qprec(1-p, alpha=alpha, u=u)
}

ppcvar = function(q, alpha=.01, u=1, tol = 1e-10) {
  fun = function(x) {dpcvar(x, alpha=alpha, u=u)}
  integrate(fun, lower = tol, upper=q)$value
}

# generate simulations on rectangular domain (e.g. unit square) from a 
# Matern Gaussian Process with zero mean
genSimsMatern = function(xRange=c(0,1), yRange=c(0,1), n=100, nsim=100, 
                         beta=(xRange[2]-xRange[1])/10, nu=1.5, sigmaSq=1, tauSq=sigmaSq/10) {
  # generate locations of data
  xs = matrix(runif(n*nsim)*(xRange[2] - xRange[1]) + xRange[1], ncol=nsim)
  ys = matrix(runif(n*nsim)*(yRange[2] - yRange[1]) + yRange[1], ncol=nsim)
  
  # simulate from standard normal
  zSims = matrix(rnorm(n*nsim), ncol=nsim)
  errs = matrix(rnorm(n*nsim, sd=sqrt(tauSq)), ncol=nsim)
  
  # generate a simulation
  genSim = function(i) {
    # print progress
    if(i %% 10 == 0)
      print(paste0("iteration ", i, "/", nsim))
    
    # use Cholesky decomposition of covariance matrix to simulate
    Sigma = stationary.cov(cbind(xs[,i], ys[,i]), Covariance="MaternLR", theta=beta, phi=sigmaSq, nu=nu)
    L = t(chol(Sigma))
    L %*% zSims[,i]
  }
  trueMat = sapply(1:nsim, genSim)
  obsMat = trueMat + errs
  
  list(xs=xs, ys=ys, trueMat=trueMat, obsMat=obsMat)
}

# modified version of fields packages Matern function to code LR2011
# parameterization of the Matern covariance
# range = beta (this is the effective range for 99% correlation, roughly)
# all other parameters are the same as in ?Matern from fields package
MaternLR = function (d, range = 1, beta=range, alpha = 1/beta, smoothness = 0.5, nu = smoothness, 
                     phi = 1) {
  Matern(sqrt(8*nu)*d*alpha, range, alpha, smoothness, nu, phi)
}

# construct numerical integration matrix, constructing mx*my aggregation regions by 
# diving xRange and yRange into mx x my grid of regions
# predPts: Ideally, predPoints should be a regular grid of locations, since they are 
# all weighed equally when being aggregated in each region
makeNumericalIntegralMat = function(predPts, xRange=c(-1, 1), yRange=c(-1, 1), mx=3, my=3) {
  # construct aggregation matrix for predictions by testing which prediction locations 
  # are in which aggregation regions
  xRegionGrid = seq(xRange[1], xRange[2], l=mx + 1)[-1]
  yRegionGrid = seq(yRange[1], yRange[2], l=my + 1)[-1]
  xRegion = function(x) {
    match(TRUE, x <= xRegionGrid)
  }
  yRegion = function(y) {
    match(TRUE, y <= yRegionGrid)
  }
  xRegionI = sapply(predPts[,1], xRegion)
  yRegionI = sapply(predPts[,2], yRegion)
  regionI = (yRegionI-1)*mx + xRegionI
  getARow = function(ai) {regionI == ai}
  
  A = t(sapply(1:(mx*my), getARow))
  A = sweep(A, 1, rowSums(A), "/")
  
  A
}

# given model output, aggregates predictions to the requested levels
# NOTE: for validation, all "obs____" variables must be modified to include the left out cluster 
#       predictions in the correct order, 
aggregateModelResultsKenya = function(dat, results, popGrid, clusterLevel=TRUE, pixelLevel=TRUE, countyLevel=TRUE, 
                                      regionLevel=TRUE) {
  stop("aggregateModelResultsKenya function is incomplete")
  # obsValues = dat$y/dat$n
  # obsCoords = cbind(dat$east, dat$north)
  # obsNs = dat$n
  # xObs = matrix(rep(1, length(obsValues)), ncol=1)
  obsUrban = dat$urban
  
  predPts = cbind(popGrid$east, popGrid$north)
  predsUrban = popGrid$urban
  predsCounty = popGrid$admin1
  predsRegion = countyToRegion(predsCounty)
  
  # Cluster level predictions (no aggregation required)
  if(clusterLevel) {
    clusterPredictions = data.frame(areaName=1:length(results$obsPreds), 
                                    urban=results$obsUrban, 
                                    preds=results$obsPreds, 
                                    SDs=results$obsSDs, 
                                    Q10=results$obsLower, 
                                    Q50=results$obsMedian, 
                                    Q90=results$obsUpper)
  } else {
    clusterPredictions = NULL
  }
  
  # Pixel level predictions (no aggregation required)
  if(pixelLevel) {
    pixelPredictions = data.frame(areaName=1:length(results$preds), 
                                  urban=predsUrban, 
                                  preds=results$preds, 
                                  SDs=results$sigmas, 
                                  Q10=results$lower, 
                                  Q50=results$median, 
                                  Q90=results$upper)
  } else {
    pixelPredictions = NULL
  }
  
  # From here onwards, we will need to aggregate predictions over the 
  # population density grid. Use the following function to get  numerical 
  # integration matrix for a given level of areal aggregation
  getIntegrationMatrix = function(areaNames) {
    densities = popGrid$popOrig
    uniqueNames = sort(unique(areaNames))
    
    integrationMatrix = sapply(1:length(uniqueNames), function(i) {
      areaI = areaNames == uniqueNames[i]
      theseDensities = densities
      theseDensities[!areaI] = 0
      theseDensities * (1/sum(theseDensities))
    })
    
    matrix(integrationMatrix, nrow=length(uniqueNames))
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames, urbanProportions) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames)
    
    # aggregate the prediction matrix
    newPredMat = t(A) %*% results$predMat
    
    # calculate relevant summary statistics
    data.frame(areaName=sort(unique(areaNames)), 
               urban=urbanProportions, 
               preds=rowMeans(newPredMat), 
               SDs=apply(newPredMat, 1, sd), 
               Q10=apply(newPredMat, 1, quantile, probs=.1), 
               Q50=apply(newPredMat, 1, quantile, probs=.5), 
               Q90=apply(newPredMat, 1, quantile, probs=.9))
  }
  
  # County level predictions
  if(countyLevel) {
    load(paste0(globalDirectory, "poppc.RData"))
    urbanProportions = poppc$popUrb / poppc$popTotal
    sortI = sort(poppc$County, index.return=TRUE)$ix
    urbanProportions = urbanProportions[sortI]
    
    countyPredictions = getIntegratedPredictions(predsCounty, urbanProportions)
  } else {
    countyPredictions = NULL
  }
  
  # Region level predictions
  if(regionLevel) {
    load(paste0(globalDirectory, "poppr.RData"))
    urbanProportions = poppr$popUrb / poppr$popTotal
    sortI = sort(as.character(poppr$Region), index.return=TRUE)$ix
    urbanProportions = urbanProportions[sortI]
    
    regionPredictions = getIntegratedPredictions(predsRegion, urbanProportions)
  } else {
    regionPredictions = NULL
  }
  
  # combine fixed effects and hyperparameter summary tables into one single table
  names(results$fixedEffectSummary) = colnames(results$parameterSummaryTable)
  parameterSummary = rbind(results$fixedEffectSummary, 
                           results$parameterSummaryTable)
  
  # return results
  list(predictions = list(clusterPredictions=clusterPredictions, 
                          pixelPredictions=pixelPredictions, 
                          countyPredictions=countyPredictions, 
                          regionPredictions=regionPredictions), 
       parameterSummary = parameterSummary)
}

# aggregates pixel predictions to the requested levels
# pg: matrix of empirical proportion draws at the pixel level
# Ng: matrix of population denominator draws at the pixel level
# popGrid: of similar format to that returned by makeDefaultPopMat() (the default)
# useDensity: whether to use population density weighting from popGrid 
#             or even weighting within the strata
# separateUrbanRural: whether or not to produce estimates in urban and 
#                     rural parts of areas separately
# popRegionNames: vector of region names with same length as nrow(popGrid)
# countyLevel, regionLevel, nationalLevel: whether or not to aggregate to the 
#                                          given levels
# countyPopTotal: the total population in the county
# normalize: whether or not to normalize the rowSums of the aggregation matrices. 
#            If set to TRUE, calculates average of Zg and Ng over the areas for any single draw, and 
#            pg is calculated as average of Zg over average of Ng
aggregatePixelPredictions = function(Zg, Ng, popGrid=NULL, useDensity=FALSE, 
                                     separateUrbanRural=TRUE, popRegionNames=NULL, 
                                     countyLevel=TRUE, regionLevel=TRUE, nationalLevel=TRUE, 
                                     countyPopTotal=NULL, regionPopTotal=NULL, nationalPopTotal=NULL, 
                                     countyPopUrban=NULL, regionPopUrban=NULL, nationalPopUrban=NULL, 
                                     country="Kenya", normalize=FALSE) {
  
  if(is.null(popGrid)) {
    popGrid = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
    
    assumeKenya = TRUE
  }
  
  predPts = cbind(popGrid$east, popGrid$north)
  predsUrban = popGrid$urban
  predsCounty = popGrid$area
  predsCountry = rep(country, length(predsCounty))
  if(is.null(popRegionNames) && regionLevel) {
    predsRegion = countyToRegion(predsCounty)
    assumeKenya = TRUE
  } else if(!is.null(popRegionNames) && regionLevel) {
    predsRegion = popRegionNames
  }
  
  if(assumeKenya) {
    warning("Since popGrid and predsRegion were not both specified, aggregation assumes country is Kenya")
  }
  
  # set NAs and pixels without any sample size to 0
  Ng[is.na(Ng)] = 0
  Zg[Ng == 0] = 0
  
  # From here onwards, we will need to aggregate predictions over the 
  # population density grid. Use the following function to get numerical 
  # integration matrix for a given level of areal aggregation. returned 
  # matrices have dimension length(unique(areaNames)) x length(areaNames)
  # areaNames: 
  # urbanProportions: vector giving proportion of population urban for each unique area in areaNames. 
  #                   If specified, ensure that urban and rural parts of the full integration 
  #                   matrix have the appropriate relative weights for each area
  # normalize: whether or not to normalize the rows of the matrices to sum to 1 or to instead 
  #            contain only binary values (or non-binary values based on the binary values if 
  #            urbanProportions is not NULL)
  getIntegrationMatrix = function(areaNames, urbanProportions=NULL, normalize=FALSE) {
    if(useDensity && normalize) {
      densities = popGrid$pop
    } else if(!useDensity) {
      densities = rep(1, nrow(popGrid))
    } else
      stop("case when useDensity == TRUE and normalize == FALSE is not currently supported")
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL) {
      areaI = areaNames == uniqueNames[i]
      theseDensities = densities
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban == thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0 && normalize)
        theseDensities * (1/thisSum)
      else if(thisSum == 0)
        rep(0, length(theseDensities))
      else
        theseDensities
    }
    
    if(!separateUrbanRural) {
      integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      
      integrationMatrix
    } else {
      integrationMatrixUrban = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE), ncol=length(uniqueNames)))
      integrationMatrixRural = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE), ncol=length(uniqueNames)))
      if(!is.null(urbanProportions)) {
        integrationMatrix = sweep(integrationMatrixUrban, 1, urbanProportions, "*") + sweep(integrationMatrixRural, 1, 1-urbanProportions, "*")
      } else {
        integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      }
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames, normalize=normalize)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!separateUrbanRural) {
      ZAggregated = A %*% Zg
      NAggregated = A %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[ZAggregated == 0] = 0
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated, A=A)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      ZAggregated = A %*% Zg
      NAggregated = A %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = 0
      
      ZAggregatedUrban = AUrban %*% Zg
      NAggregatedUrban = AUrban %*% Ng
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = 0
      
      ZAggregatedRural = ARural %*% Zg
      NAggregatedRural = ARural %*% Ng
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = 0
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated, 
           pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
           pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural, 
           A=A, AUrban=AUrban, ARural=ARural)
    }
  }
  
  # County level predictions
  if(countyLevel) {
    # load(paste0(globalDirectory, "poppc.RData"))
    # urbanProportions = poppc$popUrb / poppc$popTotal
    # sortI = sort(poppc$County, index.return=TRUE)$ix
    # urbanProportions = urbanProportions[sortI]
    
    countyMatrices = getIntegratedPredictions(predsCounty)
  } else {
    countyMatrices = NULL
  }
  
  # Region level predictions
  if(regionLevel) {
    # load(paste0(globalDirectory, "poppr.RData"))
    # urbanProportions = poppr$popUrb / poppr$popTotal
    # sortI = sort(as.character(poppr$Region), index.return=TRUE)$ix
    # urbanProportions = urbanProportions[sortI]
    
    regionMatrices = getIntegratedPredictions(predsRegion)
  } else {
    regionMatrices = NULL
  }
  
  # National level predictions
  if(nationalLevel) {
    nationalMatrices = getIntegratedPredictions(predsCountry)
  } else {
    nationalMatrices = NULL
  }
  
  # return results
  list(countyMatrices=countyMatrices, regionMatrices=regionMatrices, nationalMatrices=nationalMatrices)
}

# aggregates EA predictions to the requested levels
# pg: matrix of empirical proportion draws at the EA level
# Ng: matrix of population denominator draws at the EA level
# easpa: of similar format to that returned by makeDefaultEASPA() (the default)
# uniqueRegions: vector of length nrow(easpa) of region names associated with the areas easpa
# areaMat: a matrix of area names for each draw and each EA
# urbanMat: a matrix of urbanicity classifications for each draw and each EA
# separateUrbanRural: whether or not to produce estimates in urban and 
#                     rural parts of areas separately
# popRegionNames: vector of region names with same length as nrow(popGrid)
# countyLevel, regionLevel, nationalLevel: whether or not to aggregate to the 
#                                          given levels
# countyPopTotal: the total population in the county
# normalize: if TRUE, calculates areal means rather than areal totals
aggregateEAPredictions = function(Zc, Nc, areaMat, urbanMat, easpa=NULL, uniqueRegions=NULL, 
                                  separateUrbanRural=TRUE, 
                                  countyLevel=TRUE, regionLevel=TRUE, nationalLevel=TRUE, 
                                  countyPopTotal=NULL, regionPopTotal=NULL, nationalPopTotal=NULL, 
                                  countyPopUrban=NULL, regionPopUrban=NULL, nationalPopUrban=NULL, 
                                  country="Kenya", normalize=FALSE) {
  
  if(is.null(easpa) || (is.null(uniqueRegions) && regionLevel)) {
    warning("Since easpa or uniqueRegions were not both specified, regional aggregation assumes country is Kenya")
  }
  
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
  
  # set NAs to 0
  Zc[is.na(Zc)] = 0
  Nc[is.na(Nc)] = 0
  
  ##### 
  # We need to sort the nEA x nDraws matrices so that strata will always contain EAs in the 
  # same index range in the matrices. They are currently ordered by the pixel index instead.
  
  # this function sorts fromMat to be in the given stratum ordering
  sortMatrix = function(fromMat, fromAreaUrbanMat=NULL, toAreaUrbanMat=NULL) {
    toMat = fromMat
    uniqueStrata = sort(unique(fromAreaUrbanMat[,1]))
    totalPerStratum = table(fromAreaUrbanMat[,1])
    sortI = sort(names(totalPerStratum))
    totalPerStratum = totalPerStratum[sortI]
    
    # set toAreaUrbanMat if necessary
    if(is.null(toAreaUrbanMat) && is.null(fromAreaUrbanMat)) {
      stop("must either set toAreaUrbanMat or fromAreaUrbanMat")
    }
    if(is.null(toAreaUrbanMat)) {
      toAreaUrbanMat = matrix(rep(uniqueStrata, times=totalPerStratum), nc=ncol(fromMat), nr=nrow(fromMat))
    }
    if(is.null(fromAreaUrbanMat)) {
      fromAreaUrbanMat = matrix(rep(uniqueStrata, times=totalPerStratum), nc=ncol(fromMat), nr=nrow(fromMat))
    }
    
    # set relevant parts of the output matrix
    for(i in 1:length(uniqueStrata)) {
      thisStratum = uniqueStrata[i]
      fromI = which(fromAreaUrbanMat == thisStratum)
      toI = which(toAreaUrbanMat == thisStratum)
      toMat[toI] = fromMat[fromI]
    }
    
    toMat
  }
  
  # construct the stratum ordering matrices (nEA x nDraws)
  nDraws = ncol(Nc)
  urbanStringMat = matrix(sapply(urbanMat, function(x) {ifelse(x, "u", "r")}), ncol=nDraws)
  fromAreaUrbanMat = matrix(paste(areaMat, urbanStringMat, sep=","), ncol=nDraws)
  uniqueStrata = sort(unique(fromAreaUrbanMat[,1]))
  totalPerStratum = table(fromAreaUrbanMat[,1])
  sortI = sort(names(totalPerStratum))
  totalPerStratum = totalPerStratum[sortI]
  toAreaUrbanMat = matrix(rep(uniqueStrata, times=totalPerStratum), nc=ncol(fromAreaUrbanMat), nr=nrow(fromAreaUrbanMat))
  
  # get the (now consistent) ordering of area names and urbanicity classifications for integration
  predsArea = sapply(toAreaUrbanMat[,1], function(x) {strsplit(x, ",")[[1]][1]})
  predsUrban = sapply(toAreaUrbanMat[,1], function(x) {strsplit(x, ",")[[1]][2]}) == "u"
  predsCountry = rep(country, length(predsArea))
  
  # get the region names if necessary
  if(regionLevel && !is.null(uniqueRegions)) {
    predsRegion = predsArea
    for(i in 1:nrow(easpa)) {
      predsRegion[predsArea == easpa$area[i]] = uniqueRegions[i]
    }
  } else if(regionLevel) {
    predsRegion = countyToRegion(predsArea)
  }
  
  # reorder the input matrices
  Zc = sortMatrix(Zc, fromAreaUrbanMat, toAreaUrbanMat)
  Nc = sortMatrix(Nc, fromAreaUrbanMat, toAreaUrbanMat)
  
  # From here onwards, we will need to aggregate predictions over the 
  # areas. Use the following function to get numerical 
  # integration matrix for a given level of areal aggregation. returned 
  # matrices have dimension length(unique(areaNames)) x length(areaNames)
  # areaNames: a vector of area names with length equal to nEAs
  # urbanProportions: vector giving proportion of population urban for each unique area in areaNames. 
  #                   If specified, ensure that urban and rural parts of the full integration 
  #                   matrix have the appropriate relative weights for each area
  # normalize: whether or not to normalize the rows of the matrices to sum to 1 or to instead 
  #            contain only binary values (or non-binary values based on the binary values if 
  #            urbanProportions is not NULL)
  getIntegrationMatrix = function(areaNames, urbanProportions=NULL, normalize=FALSE) {
    densities = rep(1, nrow(Zc))
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL) {
      areaI = areaNames == uniqueNames[i]
      theseDensities = densities
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban == thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0 && normalize)
        theseDensities * (1/thisSum)
      else if(thisSum == 0)
        rep(0, length(theseDensities))
      else
        theseDensities
    }
    
    if(!separateUrbanRural) {
      integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      
      integrationMatrix
    } else {
      integrationMatrixUrban = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE), ncol=length(uniqueNames)))
      integrationMatrixRural = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE), ncol=length(uniqueNames)))
      if(!is.null(urbanProportions)) {
        integrationMatrix = sweep(integrationMatrixUrban, 1, urbanProportions, "*") + sweep(integrationMatrixRural, 1, 1-urbanProportions, "*")
      } else {
        integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      }
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames, normalize=normalize)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!separateUrbanRural) {
      ZAggregated = A %*% Zc
      NAggregated = A %*% Nc
      pAggregated = ZAggregated / NAggregated
      pAggregated[ZAggregated == 0] = 0
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      ZAggregated = A %*% Zc
      NAggregated = A %*% Nc
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = 0
      
      ZAggregatedUrban = AUrban %*% Zc
      NAggregatedUrban = AUrban %*% Nc
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = 0
      
      ZAggregatedRural = ARural %*% Zc
      NAggregatedRural = ARural %*% Nc
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = 0
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated, 
           pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
           pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
    }
  }
  
  # County level predictions
  if(countyLevel) {
    # load(paste0(globalDirectory, "poppc.RData"))
    # urbanProportions = poppc$popUrb / poppc$popTotal
    # sortI = sort(poppc$County, index.return=TRUE)$ix
    # urbanProportions = urbanProportions[sortI]
    
    countyResults = getIntegratedPredictions(predsArea)
  } else {
    countyResults = NULL
  }
  
  # Region level predictions
  if(regionLevel) {
    # load(paste0(globalDirectory, "poppr.RData"))
    # urbanProportions = poppr$popUrb / poppr$popTotal
    # sortI = sort(as.character(poppr$Region), index.return=TRUE)$ix
    # urbanProportions = urbanProportions[sortI]
    
    regionResults = getIntegratedPredictions(predsRegion)
  } else {
    regionResults = NULL
  }
  
  # National level predictions
  if(nationalLevel) {
    nationalResults = getIntegratedPredictions(predsCountry)
  } else {
    nationalResults = NULL
  }
  
  # return results
  list(countyResults=countyResults, regionResults=regionResults, nationalResults=nationalResults)
}

# aggregates values at EA level to the requested levels
# pg: matrix of empirical proportion draws at the pixel level
# Ng: matrix of population denominator draws at the pixel level
aggregateEAsKenya = function(Zg, Ng, popGrid=NULL, pixelLevel=TRUE, countyLevel=TRUE, 
                                          regionLevel=TRUE, separateUrbanRural=TRUE) {
  
  if(is.null(popGrid)) {
    popGrid = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  
  predPts = cbind(popGrid$east, popGrid$north)
  predsUrban = popGrid$urban
  predsCounty = popGrid$area
  predsRegion = countyToRegion(predsCounty)
  
  # From here onwards, we will need to aggregate predictions over the 
  # population density grid. Use the following function to get  numerical 
  # integration matrix for a given level of areal aggregation
  
  getIntegrationMatrix = function(areaNames) {
    if(useDensity) {
      densities = popGrid$popOrig
    } else {
      densities = rep(1, nrow(popGrid))
    }
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL) {
      areaI = areaNames == uniqueNames[i]
      theseDensities = densities
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban == thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0)
        theseDensities * (1/thisSum)
      else
        rep(0, length(theseDensities))
    }
    
    if(!separateUrbanRural) {
      integrationMatrix = sapply(1:length(uniqueNames), getMatrixHelper)
      
      integrationMatrix
    } else {
      integrationMatrix = sapply(1:length(uniqueNames), getMatrixHelper)
      integrationMatrixUrban = sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE)
      integrationMatrixRural = sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE)
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!separateUrbanRural) {
      ZAggregated = t(A) %*% Zg
      NAggregated = t(A) %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[ZAggregated == 0] = 0
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated)
    } else {
      A = A$integrationMatrix
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      
      ZAggregated = t(A) %*% Zg
      NAggregated = t(A) %*% Ng
      pAggregated = ZAggregated / NAggregated
      pAggregated[ZAggregated == 0] = 0
      
      ZAggregatedUrban = t(AUrban) %*% Zg
      NAggregatedUrban = t(AUrban) %*% Ng
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[ZAggregatedUrban == 0] = 0
      
      ZAggregatedRural = t(ARural) %*% Zg
      NAggregatedRural = t(ARural) %*% Ng
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[ZAggregatedRural == 0] = 0
      
      list(p=pAggregated, Z=ZAggregated, N=NAggregated, 
           pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
           pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
    }
  }
  
  # County level predictions
  if(countyLevel) {
    # load(paste0(globalDirectory, "poppc.RData"))
    # urbanProportions = poppc$popUrb / poppc$popTotal
    # sortI = sort(poppc$County, index.return=TRUE)$ix
    # urbanProportions = urbanProportions[sortI]
    
    countyMatrices = getIntegratedPredictions(predsCounty)
  } else {
    countyMatrices = NULL
  }
  
  # Region level predictions
  if(regionLevel) {
    # load(paste0(globalDirectory, "poppr.RData"))
    # urbanProportions = poppr$popUrb / poppr$popTotal
    # sortI = sort(as.character(poppr$Region), index.return=TRUE)$ix
    # urbanProportions = urbanProportions[sortI]
    
    regionMatrices = getIntegratedPredictions(predsRegion, urbanProportions)
  } else {
    regionMatrices = NULL
  }
  
  # return results
  list(countyMatrices=countyMatrices, regionMatrices=regionMatrices)
}

# Divide the domain into four parts: everything to the left of what's missing, everything to the right, and 
# the two rectangles above and below what's missing. Randomly draw how many points come from each, and sample 
# from each independently
runifsqMissingRectangle = function(n, xRange=c(-1, 1), yRange=c(-1, 1), xRangeMissing=c(-1/3, 1/3), yRangeMissing=c(-1/3, 1/3), randomPointOrdering=TRUE) {
  # Calculate the areas of the four rectangles. Probability of point being drawn from each is proportional to these areas
  areaLeft = diff(yRange) * (xRangeMissing[1] - xRange[1])
  areaRight = diff(yRange) * (xRange[2] - xRangeMissing[2])
  areaAbove = diff(xRangeMissing) * (yRange[2] - yRangeMissing[2])
  areaBelow = diff(xRangeMissing) * (yRangeMissing[1] - yRange[1])
  totalArea = areaLeft + areaRight + areaAbove + areaBelow
  probabilityLeft = areaLeft / totalArea
  probabilityRight = areaRight / totalArea
  probabilityAbove = areaAbove / totalArea
  probabilityBelow = areaBelow / totalArea
  
  # determine how many points are in each of the rectangles
  out = sample(1:4, n, replace=TRUE, prob=c(probabilityLeft, probabilityRight, probabilityAbove, probabilityBelow))
  nLeft = sum(out == 1)
  nRight = sum(out == 2)
  nAbove = sum(out == 3)
  nBelow = sum(out == 4)
  
  # sample the points
  out = rbind(runifsq(nLeft, c(xRange[1], xRangeMissing[1]), yRange), 
              runifsq(nRight, c(xRangeMissing[2], xRange[2]), yRange), 
              runifsq(nAbove, xRangeMissing, c(yRangeMissing[2], yRange[2])),
              runifsq(nBelow, xRangeMissing, c(yRange[1], yRangeMissing[1]))
  )
  
  # scramble the points if necessary
  if(randomPointOrdering) {
    reordering = sample(1:n, n, replace=FALSE)
    out = out[reordering,]
  }
  
  out
}

# generate uniform observations on a rectangle
runifsq = function(n, xRange=c(-1, 1), yRange=c(-1, 1)) {
  xs = runif(n, xRange[1], xRange[2])
  ys = runif(n, yRange[1], yRange[2])
  cbind(xs, ys)
}

# function for loading Kenya data, specifically the EA and cluster locations
# by default, use datasets with urban areas oversampled. This generates nTest testing locations 
# per simulation by sampling non-training EAs with SRS. Testing points are sampled in turn from 
# all non-training EAs, rural non-training EAs, and urban non-training EAs
# NOTE: this does not do leave one region out for training data. Instead, they are uniformly 
#       sampled from EAs not included in the respective datasets
# NOTE2: observations are generated elsewhere, since existing Kenya datasets simulate probabilities 
#        from a Matern using a SPDE model
# urbanOverSamplefrac: either 0 or 0.5
loadKenyaData = function(getOverSampledData=TRUE, nTestPerCounty=10, saveResults=TRUE, seed=123, urbanOverSamplefrac=0) {
  set.seed(seed)
  
  # load the data (the specific population doesn't matter since observations will be resimulated)
  fileName = paste0("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac", round(urbanOverSamplefrac, 4), ".RData")
  out = load(paste0("../U5MR/", fileName))
  if(getOverSampledData) {
    thisData = overSampDat
  } else {
    thisData = SRSDat
  }
  
  # get relevant data components
  clustX = sapply(thisData$clustDat, function(x) {x$east})
  clustY = sapply(thisData$clustDat, function(x) {x$north})
  clustLon = sapply(thisData$clustDat, function(x) {x$lon})
  clustLat = sapply(thisData$clustDat, function(x) {x$lat})
  clustN = sapply(thisData$clustDat, function(x) {x$numChildren})
  clustEAI = sapply(thisData$clustDat, function(x) {x$eaIs})
  clustCounties = sapply(thisData$clustDat, function(x) {x$admin1})
  clustUrban = sapply(thisData$clustDat, function(x) {x$urban})
  # clustWeight = sapply(thisData$clustDat, function(x) {x$samplingWeight}) # don't need this
  
  # get data range
  xTrain = clustX
  yTrain = clustY
  kenyaLatRange=c(-4.6, 5)
  kenyaLonRange=c(33.5, 42.0)
  eastNorthRange = projKenya(kenyaLonRange, kenyaLatRange)
  xRange = eastNorthRange[,1]
  yRange = eastNorthRange[,2]
  
  # create stratified random sampler function, sampling equal number of EAs per county
  allEAIs = thisData$eaDat$eaIs
  allCounties = thisData$eaDat$admin1
  allUrban = thisData$eaDat$urban
  
  getSamples = function(i, thisEAIs, thisCounties) {
    # this function is the same as getSamples, except makes sure to sample different EAs than the overall test samples
    thisUniqueCounties = unique(thisCounties)
    alreadyChosenEAIs = match(clustEAI[,i], thisEAIs)
    alreadyChosenEAIs = alreadyChosenEAIs[!is.na(alreadyChosenEAIs)]
    thisEAIs = thisEAIs[-alreadyChosenEAIs]
    thisCounties = thisCounties[-alreadyChosenEAIs]
    c(sapply(thisUniqueCounties, function(x) {sample(thisEAIs[thisCounties == x], nTestPerCounty, replace=FALSE)}))
  }
  
  # get testing points (sampling over all non-training EAs with stratified random sampling within counties)
  testEAI = sapply(1:100, getSamples, thisEAIs=allEAIs, thisCounties=allCounties)
  testX = apply(testEAI, 2, function(x) {thisData$eaDat$east[x]})
  testY = apply(testEAI, 2, function(x) {thisData$eaDat$north[x]})
  testLon = apply(testEAI, 2, function(x) {thisData$eaDat$lon[x]})
  testLat = apply(testEAI, 2, function(x) {thisData$eaDat$lat[x]})
  testN = apply(testEAI, 2, function(x) {thisData$eaDat$numChildren[x]})
  testCounties = apply(testEAI, 2, function(x) {thisData$eaDat$admin1[x]})
  testUrban = apply(testEAI, 2, function(x) {thisData$eaDat$urban[x]})
  
  getRestSamples = function(i, thisEAIs, thisCounties) {
    # this function is the same as getSamples, except makes sure to sample different EAs than the overall test samples
    thisUniqueCounties = unique(thisCounties)
    alreadyChosenEAIs = match(c(testEAI[,i], clustEAI[,i]), thisEAIs)
    alreadyChosenEAIs = alreadyChosenEAIs[!is.na(alreadyChosenEAIs)]
    thisEAIs = thisEAIs[-alreadyChosenEAIs]
    thisCounties = thisCounties[-alreadyChosenEAIs]
    c(sapply(thisUniqueCounties, function(x) {sample(thisEAIs[thisCounties == x], nTestPerCounty, replace=FALSE)}))
  }
  
  # get rural testing points (sampling over all rural non-training EAs with SRS)
  allEAIsRural = allEAIs[!allUrban]
  allCountiesRural = allCounties[!allUrban]
  # testEAIRural = replicate(100, getSamples(allEAIsRural, allCountiesRural))
  testEAIRural = sapply(1:100, getRestSamples, thisEAIs=allEAIsRural, thisCounties=allCountiesRural)
  testXRural = apply(testEAIRural, 2, function(x) {thisData$eaDat$east[x]})
  testYRural = apply(testEAIRural, 2, function(x) {thisData$eaDat$north[x]})
  testLonRural = apply(testEAIRural, 2, function(x) {thisData$eaDat$lon[x]})
  testLatRural = apply(testEAIRural, 2, function(x) {thisData$eaDat$lat[x]})
  testNRural = apply(testEAIRural, 2, function(x) {thisData$eaDat$numChildren[x]})
  testCountiesRural = apply(testEAIRural, 2, function(x) {thisData$eaDat$admin1[x]})
  testUrbanRural = apply(testEAIRural, 2, function(x) {allUrban[x]})
  
  # get urban testing points (sampling over all urban non-training EAs with SRS)
  allEAIsUrban = allEAIs[allUrban]
  allCountiesUrban = allCounties[allUrban]
  # testEAIUrban = replicate(100, getSamples(allEAIsUrban, allCountiesUrban))
  testEAIUrban = sapply(1:100, getRestSamples, thisEAIs=allEAIsUrban, thisCounties=allCountiesUrban)
  testXUrban = apply(testEAIUrban, 2, function(x) {thisData$eaDat$east[x]})
  testYUrban = apply(testEAIUrban, 2, function(x) {thisData$eaDat$north[x]})
  testLonUrban = apply(testEAIUrban, 2, function(x) {thisData$eaDat$lon[x]})
  testLatUrban = apply(testEAIUrban, 2, function(x) {thisData$eaDat$lat[x]})
  testNUrban = apply(testEAIUrban, 2, function(x) {thisData$eaDat$numChildren[x]})
  testCountiesUrban = apply(testEAIUrban, 2, function(x) {thisData$eaDat$admin1[x]})
  testUrbanUrban = apply(testEAIUrban, 2, function(x) {allUrban[x]})
  
  # get uniform grid of points over kenya
  out = getKenyaGrid()
  testXGrid = out$east
  testYGrid = out$north
  testLonGrid = projKenya(cbind(testXGrid, testYGrid), inverse=TRUE)[,1]
  testLatGrid = projKenya(cbind(testXGrid, testYGrid), inverse=TRUE)[,2]
  testCountiesGrid = out$admin1
  testGridUrban = out$urban
  
  # datasets eventually ideally contain the following:
  # list(xTrain=xTrain, yTrain=yTrain, zTrain=zTrain, xTest=xTest, yTest=yTest, zTest=zTest, 
  #      xGrid=xGrid, yGrid=yGrid, zGrid=zGrid, xValuesGrid=xValuesGrid, yValuesGrid=yValuesGrid, nx=nx, ny=ny, 
  #      corFun=corFun, marginalVar=marginalVar, errorVar=errorVar, xRange=xRange, yRange=yRange)
  # this dataset has:
  dataPointsKenya = list(xTrain=xTrain, yTrain=yTrain, xTest=testX, yTest=testY, 
                         xTestRural=testXRural, yTestRural=testYRural, xTestUrban=testXUrban, yTestUrban=testYUrban, 
                         xGrid=testXGrid, yGrid=testYGrid, 
                         urbanTrain=clustUrban, urbanTest=testUrban, urbanTestGrid=testGridUrban, 
                         xRange=xRange, yRange=yRange)
  
  # save results if necessary
  if(saveResults) {
    if(urbanOverSamplefrac != 0)
      oversampleText = as.character(round(urbanOverSamplefrac, 4))
    else
      oversampleText = ""
    save(dataPointsKenya, file=paste0("dataPointsKenya", oversampleText, ".RData"))
  }
  
  # plot first set of sample points for non-stratified test points
  plot(dataPointsKenya$xTrain[,1], dataPointsKenya$yTrain[,1], xlim=xRange, ylim=yRange, pch=19, cex=.1, 
       xlab="East (km)", ylab="North (km)", main="Train (Black) and Nonstratified Test (Red) Points")
  points(dataPointsKenya$xTest[,1], dataPointsKenya$yTest[,1], xlim=xRange, ylim=yRange, pch=19, cex=.2, col="red")
  
  # plot first set of sample points for rural test points
  plot(dataPointsKenya$xTrain[,1], dataPointsKenya$yTrain[,1], xlim=xRange, ylim=yRange, pch=19, cex=.1, 
       xlab="East (km)", ylab="North (km)", main="Train (Black) and Rural Test (Red) Points")
  points(dataPointsKenya$xTestRural[,1], dataPointsKenya$yTestRural[,1], xlim=xRange, ylim=yRange, pch=19, cex=.2, col="red")
  
  # plot first set of sample points for urban test points
  plot(dataPointsKenya$xTrain[,1], dataPointsKenya$yTrain[,1], xlim=xRange, ylim=yRange, pch=19, cex=.1, 
       xlab="East (km)", ylab="North (km)", main="Train (Black) and Urban Test (Red) Points")
  points(dataPointsKenya$xTestUrban[,1], dataPointsKenya$yTestUrban[,1], xlim=xRange, ylim=yRange, pch=19, cex=.2, col="red")
  
  # plot set of sample grid points
  plot(dataPointsKenya$xTrain[,1], dataPointsKenya$yTrain[,1], xlim=xRange, ylim=yRange, pch=19, cex=.1, 
       xlab="East (km)", ylab="North (km)", main="Train (Black) and Grid Test (Red) Points")
  points(dataPointsKenya$xGrid, dataPointsKenya$yGrid, xlim=xRange, ylim=yRange, pch=19, cex=.2, col="red")
  
  # return results
  dataPointsKenya
}

##### project from lat/lon to UTM northing/easting in kilometers.  Use epsg=21097
# either pass lon/east and lat/north, or a matrix with 2 columns: first being lon/east, second being lat/north
# inverse: if FALSE, projects from lon/lat to easting/northing.  Else from easting/northing to lon/lat
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    lonLatCoords = SpatialPoints(cbind(lon, lat), proj4string=CRS("+proj=longlat"))
    coordsUTM = spTransform(lonLatCoords, CRS("+init=epsg:21097 +units=km"))
    out = attr(coordsUTM, "coords")
  }
  else {
    # from easting/northing coords to lon/lat
    east = lon
    north = lat
    coordsUTM = SpatialPoints(cbind(east, north), proj4string=CRS("+init=epsg:21097 +units=km"))
    lonLatCoords = spTransform(coordsUTM, CRS("+proj=longlat"))
    out = attr(lonLatCoords, "coords")
  }
  
  out
}

# generate a grid with a fixed spatial resolution.
# either set km per cell (res), number of grid cells on longest axis (nc), 
# or number of grid cells on both axes (nx and ny).
getKenyaGrid = function(res=25, nc=NULL, nx=NULL, ny=NULL, getUrban=TRUE) {
  load("../U5MR/adminMapData.RData")
  
  kenyaLatRange=c(-4.6, 5)
  kenyaLonRange=c(33.5, 42.0)
  eastNorthRange = projKenya(kenyaLonRange, kenyaLatRange)
  xRange = eastNorthRange[,1]
  yRange = eastNorthRange[,2]
  eastLen = diff(xRange)
  northLen = diff(yRange)
  
  # set individual axis grids
  if(!is.null(nc))
    res = northLen/nc
  if(is.null(nx) && is.null(ny)) {
    # we must go off of res
    xs = seq(xRange[1], xRange[2], by=res)
    ys = seq(yRange[1], yRange[2], by=res)
  }
  else {
    xs = seq(xRange[1], xRange[2], l=nx)
    ys = seq(yRange[1], yRange[2], l=ny)
  }
  
  # get full grid
  fullGrid = make.surface.grid(list(x=xs, y=ys))
  
  # get subset of points inside Kenya polygon
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  kenyaPolyProj = projKenya(kenyaPoly)
  inKenya = in.poly(fullGrid, kenyaPolyProj)
  out = fullGrid[inKenya,]
  
  # calculate whether each point is urban or rural and which county it's in:
  urban = getUrbanRural(out)
  out = data.frame(east=out[,1], north=out[,2], urban=urban, admin1=getRegion(out)$regionNames)
  
  # plot the coordinates
  plot(out$east, out$north, pch=19, cex=.1, xlab="Easting (km)", ylab="Northing (km)", type="n")
  points(out$east[!urban], out$north[!urban], pch=19, cex=.1, col="green")
  points(out$east[urban], out$north[urban], pch=19, cex=.1, col="blue")
  
  out
}

# for computing what administrative regions the given points are in
# project: project to longitude/latitude coordinates
getRegion = function(points, project=FALSE) {
  load("../U5MR/adminMapData.RData")
  mapDat = adm1
  
  # project points to lon/lat coordinate system if user specifies
  if(project)
    points = projKenya(points, inverse=TRUE)
  
  regionNames = mapDat@data$NAME_1
  
  # make sure county names are consistent for mapDat == adm1
  regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
  regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  
  # get region map polygons and set helper function for testing if points are in the regions
  polys = mapDat@polygons
  inRegion = function(i) {
    countyPolys = polys[[i]]@Polygons
    inside = sapply(1:length(countyPolys), function(x) {in.poly(points, countyPolys[[x]]@coords, inflation=0)})
    insideAny = apply(inside, 1, any)
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inRegion)
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  regionID = apply(out, 1, function(vals) {match(1, vals != 0)})
  regionNameVec = regionNames[regionID]
  list(regionID=regionID, regionNames=regionNameVec, multipleRegs=multipleRegs)
}

# convert county to region
countyToRegion = function(countyNames) {
  load(paste0(globalDirectory, "kenya-prov-county-map.RData"))
  regionIs = match(as.character(countyNames), as.character(ctp[,1]))
  as.character(ctp[regionIs,2])
}

getUrbanRural = function(utmGrid) {
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("../U5MR/Kenya2014Pop/pop.RData")
  load("../U5MR/adminMapData.RData")
  load("../U5MR/Kenya2014Pop/pop.RData")
  
  # project coordinates into lat/lon
  lonLatGrid = projKenya(utmGrid, inverse=TRUE)
  
  # get population density at those coordinates
  interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  
  # compute counties associated with locations
  counties = getRegion(lonLatGrid, FALSE)$regionNames
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, admin1=counties))
  threshes = setThresholds()
  popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$counties == newPop$admin1[i]]})
  badPoints = which(sapply(popThreshes, length) > 1)
  if(length(badPoints) >= 1) {
    for(i in 1:length(badPoints)) {
      warnings(paste0("Point ", badPoints[i], " not in Kenya"))
      popThreshes[[badPoints[i]]] = NA
    }
  }
  popThreshes = unlist(popThreshes)
  urban = newPop$popOrig > popThreshes
  
  urban
}

# set thresholds within each county based on percent population urban
setThresholds = function() {
  require(raster)
  
  load("../U5MR/kenyaPopProj.RData")
  load("../U5MR/adminMapData.RData")
  load("../U5MR/poppc.RData")
  
  getCountyThresh = function(countyName) {
    # if Nairobi or Mombasa, always urban
    if((countyName == "Nairobi") || (countyName == "Mombasa"))
      return(-Inf)
    
    # do the setup
    thisCounty = as.character(kenyaPop$admin1) == countyName
    thisPop = kenyaPop$popOrig[thisCounty]
    thisTot = sum(thisPop)
    pctUrb = poppc$pctUrb[poppc$County == countyName]/100
    pctRural = 1 - pctUrb
    
    # objective function to minimize
    # objFun = function(thresh) {
    #   curPctUrb = sum(thisPop[thisPop > thresh])/thisTot
    #   (curPctUrb - pctUrb)^2
    # }
    
    # do optimization
    # out = optim(10, objFun)
    # thresh = out$par
    # out = optimize(objFun, c(.01, 50000))
    # thresh = out$par
    
    # calculate threshold by integrating ecdf via sorted value cumulative sum
    sortedPop = sort(thisPop)
    cumsumPop = cumsum(sortedPop)
    threshI = match(1, cumsumPop >= thisTot*pctRural)
    thresh = sortedPop[threshI]
    
    # print(paste0("pctUrb: ", pctUrb, "; resPctUrb: ", sum(thisPop[thisPop > thresh])/thisTot, "; thresh: ", thresh, "; obj: ", out$objective))
    thresh
  }
  
  # compute threshold for each county
  counties = poppc$County
  threshes = sapply(counties, getCountyThresh)
  
  list(counties=counties, threshes=threshes)
}

# for plotting administration data
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using projKenya function.  This can be used when plotting the projected `east' 
#          and `north' variables in kenyaEAs for instance.
# ...: arguments to polygon function
plotMapDat = function(plotVar=NULL, varCounties=NULL, zlim=NULL, project=FALSE, cols=tim.colors(), 
                      legend.mar=7, new=FALSE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                      ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, mapDat = NULL, addColorBar=TRUE, 
                      legendArgs=list(), leaveRoomForLegend=TRUE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), ...) {
  # load necessary data
  if(is.null(mapDat)) {
    if(length(plotVar) == 47) {
      out = load("../U5MR/adminMapData.RData")
      mapDat = adm1
    } else if(length(plotVar) == 8) {
      # shape file found at: https://jlinden.carto.com/tables/kenya_region_shapefile/public
      require(maptools)
      mapDat = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
    } else {
      out = load("../U5MR/adminMapData.RData")
      mapDat = adm0
    }
  }
  if(is.null(varCounties)) {
    if(length(mapDat) != 8) {
      out = load("../U5MR/kenyaData.RData")
      varCounties=sort(as.character(unique(mort$admin1)))
    } else {
      out = load("../U5MR/poppr.RData")
      varCounties = sort(as.character(poppr$Region))
    }
  }
  
  # do setup for ploting data by county if necessary
  if(!is.null(plotVar)) {
    if(is.null(zlim)) {
      zlim = range(plotVar)
    }
    
    # get region names from map data
    if(!is.null(mapDat@data$NAME_1)) {
      regionNames = mapDat@data$NAME_1
    } else if(!is.null(mapDat@data$name_1)) {
      regionNames = as.character(mapDat@data$name_1)
    } else {
      stop("mapDat has unrecognized region names")
    }
    
    # make sure county names are consistent for mapDat == adm1
    regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
    regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
    
    # make sure county names are consistent for plotting regions rather than counties
    regionNames[regionNames == "North-Eastern"] = "North Eastern"
  }
  
  # generate new plot if necessary
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(project) {
      if(is.null(xlab))
        xlab = "East (km)"
      if(is.null(xlim))
        xlim = eastLim
      if(is.null(ylab))
        ylab = "North (km)"
      if(is.null(ylim))
        ylim = northLim
    }
    else {
      if(is.null(xlab))
        xlab = "Longitude"
      if(is.null(xlim))
        xlim = kenyaLonRange
      if(is.null(ylab))
        ylab = "Latitude"
      if(is.null(ylim))
        ylim = kenyaLatRange
    }
    if(is.null(main))
      main = ""
    
    if(is.null(plotArgs)) {
      plotArgs = list(main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=asp)
    } else {
      plotArgs$main = main
      plotArgs$xlab = xlab
      plotArgs$ylab = ylab
      plotArgs$xlim = xlim
      plotArgs$ylim = ylim
      plotArgs$asp = asp
    }
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    do.call("plot", c(list(1, 2, type="n"), plotArgs))
  }
  
  # add polygons to plot
  polys = mapDat@polygons
  plotCounty = function(i) {
    countyPolys = polys[[i]]@Polygons
    
    if(is.null(plotVar)) {
      if(!project)
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords), list(...)))})
      else
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(projKenya(countyPolys[[x]]@coords)), list(...)))})
    }
    else {
      # get index of plotVar corresponding to this county
      thisI = which(varCounties == regionNames[i])
      
      # get color to plot
      vals = c(zlim, scaleFun(plotVar[thisI]))
      vals = vals-vals[1]
      vals = vals/(vals[2] - vals[1])
      col = cols[round(vals[3]*(length(cols)-1))+1]
      
      if(!project)
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords, col=col), list(...)))})
      else
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(projKenya(countyPolys[[x]]@coords), col=col), list(...)))})
    }
    
  }
  sapply(1:length(polys), plotCounty)
  
  if(!is.null(plotVar) && addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(scaleFunInverse(zlim), n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=zlim
    legendArgs$nlevel=length(cols)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=cols
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

# generate the population density surface along with urbanicity estimates
makeInterpPopGrid = function(kmRes=5, adjustPopSurface=FALSE, targetPop=c("children", "women")) {
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("../U5MR/Kenya2014Pop/pop.RData")
  load("../U5MR/lims.RData")
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1
  
  # get a rectangular grid
  eastGrid = seq(eastLim[1], eastLim[2], by=kmRes)
  northGrid = seq(northLim[1], northLim[2], by=kmRes)
  utmGrid = make.surface.grid(list(east=eastGrid, north=northGrid))
  
  # project coordinates into lat/lon
  lonLatGrid = projKenya(utmGrid, inverse=TRUE)
  
  # subset grid so it's in Kenya
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(lonLatGrid, kenyaPoly)
  utmGrid = utmGrid[inKenya,]
  lonLatGrid = lonLatGrid[inKenya,]
  
  # get population density at those coordinates
  interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  
  # compute counties associated with locations
  counties = getRegion(lonLatGrid)$regionNames
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, admin1=counties))
  threshes = setThresholds()
  popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$counties == newPop$admin1[i]]})
  urban = newPop$popOrig > popThreshes
  newPop$urban = urban
  
  newPop$east = utmGrid[,1]
  newPop$north = utmGrid[,2]
  
  # if necessary, adjust the population surface so that it better represents the the child population density 
  # rather than the total population density
  if(adjustPopSurface) {
    targetPop = match.arg(targetPop)
    
    # sort easpc by county name alphabetically
    counties=sort(unique(poppc$County))
    sortI = sort(easpc$County, index.return=TRUE)$ix
    temp = easpc[sortI,]
    
    # calculate the number of children per stratum using true total eas and empirical children per ea from census data
    load("../U5MR/empiricalDistributions.RData")
    if(targetPop == "children") {
      targetPopPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
        ecdfExpectation(empiricalDistributions$childrenUrban)
      targetPopPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
        ecdfExpectation(empiricalDistributions$childrenRural)
    }
    else {
      targetPopPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$womenUrban)
      targetPopPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$womenRural)
    }
    
    # generate 2 47 x nPixels matrices for urban and rural strata integrating pixels with respect to population density to get county estimates
    getCountyStratumIntegrationMatrix = function(getUrban=TRUE) {
      counties = as.character(counties)
      
      mat = t(sapply(counties, function(countyName) {newPop$admin1 == countyName}))
      mat = sweep(mat, 2, newPop$popOrig, "*")
      sweep(mat, 2, newPop$urban == getUrban, "*")
    }
    urbanIntegrationMat = getCountyStratumIntegrationMatrix()
    ruralIntegrationMat = getCountyStratumIntegrationMatrix(FALSE)
    
    # calculate number of people per stratum by integrating the population density surface
    urbanPopulations = rowSums(urbanIntegrationMat)
    ruralPopulations = rowSums(ruralIntegrationMat)
    
    # adjust each row of the integration matrices to get the correct expected number of children per stratum
    urbanIntegrationMat = sweep(urbanIntegrationMat, 1, targetPopPerStratumUrban / urbanPopulations, "*")
    ruralIntegrationMat = sweep(ruralIntegrationMat, 1, targetPopPerStratumRural / ruralPopulations, "*")
    ruralIntegrationMat[ruralPopulations == 0,] = 0
    
    # the column sums of the matrices give the correct modified population densities
    newPop$popOrig = colSums(urbanIntegrationMat) + colSums(ruralIntegrationMat)
  }
  
  newPop
}

# add in binomial variation to the probability sampling matrix
addBinomialVar = function(probMatrix, ns) {
  simulatedObservations = matrix(rbinom(n=length(probMatrix), size=rep(ns, ncol(probMatrix)), prob=c(as.matrix(probMatrix))), nrow=nrow(probMatrix))
  sweep(simulatedObservations, 1, 1/ns, "*")
}

# create list of inclusion logicals, partitioning the given data set into 8 (by default) pieces. These 
# pieces are of roughly equal size, and are chosen by a stratified sampling without replacement
getValidationI = function(dat=NULL, dataType=c("mort", "ed"), allCountyNames=NULL, urban=NULL, nFold=8, seed=123) {
  set.seed(seed)
  dataType = match.arg(dataType)
  
  if(is.null(allCountyNames) || is.null(urban)) {
    if(is.null(dat)) {
      if(dataType == "mort") {
        out = load("../U5MR/kenyaData.RData")
        dat = mort
      }
      else {
        out = load("../U5MR/kenyaDataEd.RData")
        dat = ed
      }
    }
    
    allCountyNames = dat$admin1
    urban = dat$urban
  }
  
  strata = sapply(1:length(urban), function(i) {paste0(allCountyNames[i], ifelse(urban[i], "U", "R"))})
  uniqueStrata = unique(strata)
  nPerStrata = table(strata)
  nPerStrata = nPerStrata[sort(match(names(nPerStrata), uniqueStrata), index.return=TRUE)$ix]
  nPerStrataFoldLow = floor(nPerStrata / nFold)
  nPerStrataFoldHigh = nPerStrataFoldLow + 1
  nHiFoldsPerStrata = nPerStrata %% nFold
  nLowFoldsPerStrata = nFold - nHiFoldsPerStrata
  
  # for each stratum, and each fold, determine the number of samples
  nSamplesTable = matrix(nPerStrataFoldLow, nrow=length(uniqueStrata), ncol=nFold)
  for(i in 1:nrow(nSamplesTable)) {
    if(nLowFoldsPerStrata[i] < nFold)
      nSamplesTable[i,(nLowFoldsPerStrata[i]+1):ncol(nSamplesTable)] = nSamplesTable[i,(nLowFoldsPerStrata[i]+1):ncol(nSamplesTable)] + 1
  }
  
  # construct a list of vectors, one for each fold, each vector containing the indices of the observations in the fold
  foldIndices = as.list(1:nFold)
  for(i in 1:length(uniqueStrata)) {
    # for each unique stratum, get the list of observation indices, and partition them into the folds
    stratumIndices = which(strata == uniqueStrata[i])
    
    for(j in 1:nFold) {
      
      # sample the fold observations from the remaining on sampled observations in the stratum
      if(length(stratumIndices) != 1)
        theseIndices = sample(stratumIndices, nSamplesTable[i,j], replace=FALSE)
      else
        theseIndices = stratumIndices # in this case, sample function does something different
      
      if(i == 1) {
        foldIndices[[j]] = theseIndices
      } else {
        foldIndices[[j]] = c(foldIndices[[j]], theseIndices)
      }
      
      # remove the sampled observations from the vector of observation indices in the stratum
      stratumIndices = stratumIndices[-match(theseIndices, stratumIndices)]
    }
  }
  
  # convert foldIndices to a matrix of logical
  sampleMatrix = matrix(FALSE, nrow=length(urban), ncol=nFold)
  for(i in 1:nFold) {
    sampleMatrix[foldIndices[[i]],i] = TRUE
  }
  
  # return results
  sampleMatrix
}

getAveragePredictionDistance = function(areaLevel=c("Region", "County", "Pixel", "Cluster")) {
  areaLevel = match.arg(areaLevel)
}

getArea = function(areaLevel=c("Region", "County")) {
  areaLevel = match.arg(areaLevel)
  require(shapefiles)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  getOneArea = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    # thisCentroid = centroid(poly@coords)
    # allPoints = rbind(thisCentroid, 
    #                   poly@coords)
    allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    # downsample spatial polygons as necessary
    if(nrow(allProjected) >= 6000) {
      # simplify the polygon
      newProjected = allProjected
      tolerance = .05
      while(nrow(newProjected) >= 6000) {
        newProjected = do.call("cbind", dp(list(x=allProjected[,1], y=allProjected[,2]), tolerance=tolerance))
        tolerance = tolerance * 2
      }
      allProjected = newProjected
    }
    
    maps:::area.polygon(allProjected)
  }
  # calculate areas in km^2
  areas = sapply(thisMap@polygons, getOneArea)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  areas = areas[sortI]
  names(areas) = areaNames[sortI]
  areas
}

getAreaPerObservation = function(areaLevel=c("Region", "County"), dataType=c("ed", "mort"), doLog = TRUE) {
  areaLevel = match.arg(areaLevel)
  areas = getArea(areaLevel)
  
  # load the dataset
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  
  # get area names
  allCounties = as.character(dat$admin1)
  counties = sort(unique(as.character(dat$admin1)))
  allRegions = countyToRegion(allCounties)
  regions = sort(unique(allRegions))
  if(areaLevel == "Region") {
    allAreas = allRegions
    uniqueAreas = regions
  } else if(areaLevel == "County") {
    allAreas = allCounties
    uniqueAreas = counties
  }
  
  # count the number of observations per area
  counts = table(allAreas)
  sortI = sort(names(counts), index.return=TRUE)$ix
  counts = counts[sortI]
  
  if(!doLog) {
    areas / counts
  } else {
    log(areas) - log(counts)
  }
}

getAreaPerObservationTicks = function(areaLevel=c("Region", "County"), dataType=c("ed", "mort")) {
  areaLevel = match.arg(areaLevel)
  dataType = match.arg(dataType)
  
  if(dataType != "ed")
    stop("Only education dataset currently supported")
  
  if(areaLevel == "Region") {
    c(50, 100, 200, 400, 800, 1200)
  } else {
    c(10, 25, 50, 100, 200, 400, 800, 1600, 2400)
  }
}

getRadius = function(areaLevel=c("Region", "County")) {
  require(geosphere)
  require(shapefiles)
  areaLevel = match.arg(areaLevel)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  # thisMap@polygons[[1]]@Polygons[[1]]@coords
  # thisMap@polygons[[1]]@Polygons[[1]]@area
  getOneRadius = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    # thisCentroid = centroid(poly@coords)
    # allPoints = rbind(thisCentroid, 
    #                   poly@coords)
    allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    # downsample spatial polygons as necessary
    if(nrow(allProjected) >= 6000) {
      # simplify the polygon
      newProjected = allProjected
      tolerance = .05
      while(nrow(newProjected) >= 6000) {
        newProjected = do.call("cbind", dp(list(x=allProjected[,1], y=allProjected[,2]), tolerance=tolerance))
        tolerance = tolerance * 2
      }
      allProjected = newProjected
    }
    max(rdist(allProjected)) / 2
  }
  # calculate radii in km
  radii = sapply(thisMap@polygons, getOneRadius)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  radii = radii[sortI]
  names(radii) = areaNames[sortI]
  radii
}

getMeanRadius = function(areaLevel=c("Region", "County")) {
  require(geosphere)
  require(shapefiles)
  areaLevel = match.arg(areaLevel)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  # thisMap@polygons[[1]]@Polygons[[1]]@coords
  # thisMap@polygons[[1]]@Polygons[[1]]@area
  getOneRadius = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    thisCentroid = centroid(poly@coords)
    allPoints = rbind(thisCentroid,
                      poly@coords)
    # allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    mean(rdist.vec(matrix(allProjected[1,], nrow=nrow(allProjected)-1, ncol=2, byrow=TRUE), allProjected[-1,]))
  }
  # calculate radii in km
  radii = sapply(thisMap@polygons, getOneRadius)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  radii = radii[sortI]
  names(radii) = areaNames[sortI]
  radii
}

getPredictionDistance = function(doLog=TRUE, dataType=c("ed", "mort")) {
  # load the pixel/prediction grid
  out = load("../U5MR/popGrid.RData")
  predictionPoints = cbind(popGrid$east, popGrid$north)
  
  # load the dataset
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  observationPoints = cbind(dat$east, dat$north)
  
  # calculate distance of nearest observation to each prediction point
  distances = rdist(observationPoints, predictionPoints)
  distances = apply(distances, 2, min)
  
  if(doLog)
    log(distances)
  else
    distances
}

getPredictionDistanceTicks = function(dataType=c("ed", "mort")) {
  dataType = match.arg(dataType)
  
  if(dataType != "ed")
    stop("Only education dataset currently supported")
  
  c(.1, 1, 5, 10, 20, 40, 80)
}

# simulate from the ELK model with the following parameters :
# nsim: the number of simulations
# x: matrix of spatial locations, each row being a single location
# alpha: proportion of spatial variance going into each layer
# effectiveRange: effective range of each layer, or just the first layer for the ELK-F model
# sigmaEpsilon: nugget standard deviation
# rho: spatial variance (or proportional to spatial variance if normalize==FALSE)
# gridInfo: latticeInfo object
# kappa: cap up parameter from the original LatticeKrig model
# normalize: whether or not to normalize the spatial process so that rho is the spatial variance
# fastNormalize: if FALSE, sets a variance at each basis knot to be rho. If TRUE, only sets 
#                variance in the center of the domain to be rho.
rELK = function(nsim=1, x, alpha, effectiveRange=NULL, sigmaEpsilon, rho=1, gridInfo, kappa=NULL, 
                normalize=TRUE, fastNormalize=TRUE) {
  # first do any reparameterization necessary
  deltas = sapply(gridInfo, function(x) {x$latWidth})
  kappa = effectiveRange / sqrt(8) * deltas
  
  # construct precision matrix, simulate basis coefficients from it
  Q = makeQ(kappa, rho, gridInfo, alphas=alpha, normalized=normalize, fastNormalize=fastNormalize)
  simCoefficients = inla.qsample(nsim, Q)
  
  # construct basis matrix, construct simulations at the spatial locations
  A = makeA(x, gridInfo)
  A %*% simCoefficients
}

# draw random numbers from an ecdf object
recdf = function(n, distribution) {
  probs = runif(n)
  # quantile(distribution, probs, type=1) # type=1 signifies inverse ecdf
  if(!("function" %in% class(distribution)) && !is.null(distribution$rfun))
    distribution$rfun(n)
  else
    ecdfQuantile(probs, distribution)
}

# draw random numbers from a composition of ecdf objects. I.e. we want children per EA and 
# we have households/EA, mothers/household, and children/mother
recdfComposed = function(n, distributions) {
  results = rep(1, n)
  for(i in 1:length(distributions)) {
    results = sapply(results, function(x) {sum(recdf(x, distributions[[i]]))})
  }
  
  results
}

# draw random numbers from an ecdf object (this could be implemented efficiently)
decdf = function(x, distribution) {
  distributionKnots = knots(distribution)
  masses = diff(c(0, evalq(y, environment(distribution))))
  i = match(x, distributionKnots)
  if(length(i) == 1 && is.na(i))
    return(0)
  else {
    out = masses[i]
    out[is.na(i)] = 0
    return(0)
  }
}

# get expectation of an ecdf object
ecdfExpectation = function(distribution) {
  distributionKnots = knots(distribution)
  distributionKnots = c(distributionKnots[1] - 1, distributionKnots)
  probs = distribution(distributionKnots[2:length(distributionKnots)]) - distribution(distributionKnots[1:(length(distributionKnots) - 1)])
  sum(distributionKnots[2:length(distributionKnots)] * probs)
}
# ecdfExpectation(empiricalDistributions$households) * ecdfExpectation(empiricalDistributions$mothers) * ecdfExpectation(empiricalDistributions$children)

# Becuase stats:::quantile.ecdf is terribly inefficient...
ecdfQuantile = function(p, distribution) {
  distributionKnots = knots(distribution)
  cumulativeMass = evalq(y, environment(distribution))
  indices= sapply(p, function(x) {match(TRUE, x <= cumulativeMass)})
  distributionKnots[indices]
}

ecdf2edfun = function(distribution) {
  samples = evalq(rep.int(x, diff(c(0, round(nobs * y)))), environment(distribution))
  edfun(samples)
}

# takes the poppc table, containing the proportion of population that is urban and rural in each stratum, and
# adjusts it to be representative of the children in urban and rural areas per stratum based on census data
adjustPopulationPerCountyTable = function(dataType=c("children", "women")) {
  dataType = match.arg(dataType)
  
  # calculate the number of childrenor women per stratum using true total eas and empirical children per ea from census data
  load(paste0(globalDirectory, "empiricalDistributions.RData"))
  if(dataType == "children") {
    targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
      ecdfExpectation(empiricalDistributions$childrenUrban)
    targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
      ecdfExpectation(empiricalDistributions$childrenRural)
  }
  else {
    targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$womenUrban)
    targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$womenRural)
  }
  
  
  # adjust poppc table to be representative of the number of children per stratum
  newPopTable = poppc
  targetPopPerCounty = targetPopPerStratumUrban + targetPopPerStratumRural
  newPopTable$popUrb = targetPopPerStratumUrban
  newPopTable$popRur = targetPopPerStratumRural
  newPopTable$popTotal = targetPopPerCounty
  newPopTable$pctUrb = newPopTable$popUrb / targetPopPerCounty  * 100
  newPopTable$pctTotal = newPopTable$popTotal/sum(newPopTable$popTotal) * 100
  
  # return results
  newPopTable
}

rMyMultinomial = function(n, i, includeUrban=TRUE, urban=TRUE, popMat=NULL, easpa=NULL) {
  
  # set default inputs
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
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # determine which pixels and how many EAs are in this stratum
  if(includeUrban) {
    includeI = popMat$area == areas[i] & popMat$urban == urban
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EARur[i])
  }
  else {
    includeI = popMat$area == areas[i]
    nEA = ifelse(urban, easpa$EAUrb[i], easpa$EATotal[i])
  }
  
  # sample from the pixels if this stratum exists
  if(sum(includeI) == 0)
    return(matrix(nrow=0, ncol=n))
  thesePixelProbs = popMat$pop[includeI]
  rmultinom(n, nEA, prob=thesePixelProbs)
}

# function for determining how to recombine separate multinomials into the draws over all pixels
getSortIndices = function(i, urban=TRUE, popMat=NULL, includeUrban=TRUE) {
  
  # set default inputs
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
  
  # get area names
  areas = sort(unique(popMat$area))
  
  # determine which pixels and how many EAs are in this stratum
  if(includeUrban) {
    includeI = popMat$area == areas[i] & popMat$urban == urban
  }
  else {
    includeI = popMat$area == areas[i]
  }
  
  which(includeI)
}

# gives nPixels x n matrix of draws from the stratified multinomial with values 
# corresponding to the value of |C^g| for each pixel, g (the number of EAs/pixel)
rStratifiedMultnomial = function(n, popMat=NULL, easpa=NULL, includeUrban=TRUE) {
  
  # set default inputs
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
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # we will need to draw separate multinomial for each stratum. Start by 
  # creating matrix of all draws of |C^g|
  eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
  
  # now draw multinomials
  if(includeUrban) {
    # draw for each area crossed with urban/rural
    urbanSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE, 
                                           includeUrban=includeUrban, popMat=popMat, easpa=easpa))
    ruralSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=FALSE, 
                                           includeUrban=includeUrban, popMat=popMat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    urbanIndices = unlist(sapply(1:length(areas), getSortIndices, urban=TRUE, popMat=popMat, includeUrban=includeUrban))
    ruralIndices = unlist(sapply(1:length(areas), getSortIndices, urban=FALSE, popMat=popMat, includeUrban=includeUrban))
    
    # recombine into eaSamples
    eaSamples[urbanIndices,] = urbanSamples
    eaSamples[ruralIndices,] = ruralSamples
  } else {
    # draw for each area
    stratumSamples = rbind(sapply(1:length(areas), n=n, rMyMultinomial, 
                                  includeUrban=includeUrban, popMat=popMat, easpa=easpa))
    
    # get the indices used to recombine into the full set of draws
    stratumIndices = c(sapply(1:length(areas), getSortIndices, popMat=popMat, includeUrban=includeUrban))
    
    # recombine into eaSamples
    eaSamples[stratumIndices,] = stratumSamples
  }
  
  # return results
  eaSamples
}

# taken from logitnorm package.  Calculates the mean of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat has a mean and standard deviation 
# on the logit scale
logitNormMean = function(muSigmaMat, parClust=NULL, logisticApproximation=TRUE, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      apply(muSigmaMat, 1, logitNormMean, ...)
    else
      parApply(parClust, muSigmaMat, 1, logitNormMean, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if(sigma == 0)
      expit(mu)
    else {
      if(any(is.na(c(mu, sigma))))
        NA
      else if(!logisticApproximation) {
        # numerically calculate the mean
        fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
        integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
      } else {
        # use logistic approximation
        k = 16 * sqrt(3) / (15 * pi)
        expit(mu / sqrt(1 + k^2 * sigma^2))
      }
    }
  }
}

# calculate the expected value of a summation under a Poisson distribution
expectSummationPoisson = function() {
  1
}

