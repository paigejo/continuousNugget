# same as runBYM2, except fits a single data set (the ed global data frame)
# doPredsAtPostMean: if TRUE, fix all model hyperparameters at the posterior mean
# getPosteriorDensity: EXPERIMENTAL: evaluate the posterior density of direct estimates using multivariate normal approximation
# directLogitEsts: must be specified, and used only if, getPosteriorDensity is TRUE. posterior density 
#                  is evaluated at the direct estimates.
modBYM2 = function(dat=mort, urbanEffect=TRUE, clusterEffect=TRUE, saveResults=FALSE, fileNameRoot="", 
                   previousResult=NULL, doValidation=FALSE, predCountyI=NULL, counties=sort(unique(poppc$County)), 
                   doPredsAtPostMean=FALSE, directLogitEsts=NULL, fixedParameters=FALSE, getPosteriorDensity=FALSE,
                   arealLevel=c("County", "Constituency"), poppc=poppc, poppcon=poppcon) {
  arealLevel = match.arg(arealLevel)
  
  # remove data from the given county for validation if necessary
  if(!is.null(predCountyI))
    dat$y[as.character(dat$admin1)==counties[predCountyI]] = NA
  
  includeUrban = urbanEffect
  
  # Get true ratios of urban/rural
  if(fileNameRoot == "Ed")
    thisPopTable = adjustPopulationPerCountyTable("women")
  else
    thisPopTable = adjustPopulationPerCountyTable("children")
  urbRatio = vector('numeric', length = 47)
  counties = sort(unique(as.character(dat$admin1)))
  urbRatio = thisPopTable$popUrb / thisPopTable$popTotal
  sortI = matchMultiple(counties, thisPopTable$County)
  urbRatio = urbRatio[sortI]
  
  # Define formula
  if(urbanEffect) {
    if(clusterEffect) {
      formula = y ~ urban +
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 2/3), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(1,0.01))))
    } else {
      formula = y ~ urban + 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 2/3), prior="pc")))
    }
  } else {
    if(clusterEffect) {
      formula = y ~ f(idx, model="bym2",
                      graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
                      hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), 
                                 phi=list(param=c(0.5, 2/3), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(1,0.01))))
    } else {
      formula = y ~ 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 2/3), prior="pc")))
    }
  }
  
  # Number of simulations for producing results
  Nsim = 1000
  
  # Help functions
  logit = function(x){
    return(log(x/(1-x)))
  }
  expit = function(x){
    return(1/(1+exp(-x)))
  }
  
  if(clusterEffect)
    parNames = c("Intercept", "Urban", "Cluster Var", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "Cluster SD", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  else
    parNames = c("Intercept", "Urban", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  includeI = c(1, rep(2, includeUrban), 3:length(parNames))
  parNames = parNames[includeI]
  
  # Go through education data set
  sampCountyDat = matrix(NA, nrow=47, ncol=Nsim)
  sampClusterDatUrban = matrix(NA, nrow=47, ncol=Nsim)
  sampClusterDatRural = matrix(NA, nrow=47, ncol=Nsim)
  sampPixelDatUrban = matrix(NA, nrow=47, ncol=Nsim)
  sampPixelDatRural = matrix(NA, nrow=47, ncol=Nsim)
  
  # for the final parameters to store, 2 fixed effects, 2 + clusterEffect estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountyDatPar = numeric(length(parNames))
  sampCountyDatSD = numeric(length(parNames))
  sampCountyDat10 = numeric(length(parNames))
  sampCountyDat50 = numeric(length(parNames))
  sampCountyDat90 = numeric(length(parNames))
  
  names(sampCountyDatPar) = parNames
  names(sampCountyDatSD) = parNames
  names(sampCountyDat10) = parNames
  names(sampCountyDat50) = parNames
  names(sampCountyDat90) = parNames
  
  # Extract data
  dat$admin1 = factor(dat$admin1)
  
  # INLA data
  dat = list(y = dat$y,
             Ntrials = dat$n,
             urban = dat$urban,
             idx = as.numeric(dat$admin1),
             idxEps = 1:length(dat$y))
  
  # Add unobserved data to make sampling easier
  dat$y = c(rep(NA, 47*2), dat$y)
  dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
  dat$urban = c(rep(c(1,0), each = 47), dat$urban)
  dat$idx = c(rep(1:47, 2), dat$idx)
  dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
  
  # set posterior approximation for calculating validation quantities
  
  if(!doPredsAtPostMean) {
    # based on http://www.r-inla.org/faq#TOC-How-can-I-compute-cross-validation-or-predictive-measures-of-fit-
    # when calculating log posterior density, we use fewer integration points to get better estimates of the 
    # covariance matrices at each integration point
    if(!getPosteriorDensity)
      control.inla = list(strategy="laplace", int.strategy="grid", diff.logdens=4, npoints=21)
    else
      control.inla = list(strategy="laplace", int.strategy="grid", diff.logdens=4, npoints=9)
  }
  else {
    # we need to generate predictions in this case based on a fixed set of hyperparameters
    control.inla = list(strategy="laplace", int.strategy="eb") 
  }
  
  # initialize the fitting process based on a previous optimum if necessary
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousResult)) {
    # initialize the fitting process based on a previous optimum
    # modeControl$result = previousResult
    modeControl$theta = previousResult$mode$theta
    modeControl$x = previousResult$mode$x
    modeControl$restart = !fixedParameters
    modeControl$fixed = fixedParameters
  }
  
  # Run model
  print("fitting BYM model...")
  result = inla(formula = formula, 
                family="binomial",
                Ntrials = Ntrials,
                data=dat, 
                control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                quantiles=c(0.1, 0.5, 0.9), 
                control.mode=modeControl)
  
  ## include parameter estimates in the table
  # fixed effects
  sampCountyDatPar[1:(1 + includeUrban)] = result$summary.fixed[,1]
  sampCountyDatSD[1:(1 + includeUrban)] = result$summary.fixed[,2]
  sampCountyDat10[1:(1 + includeUrban)] = result$summary.fixed[,3]
  sampCountyDat50[1:(1 + includeUrban)] = result$summary.fixed[,4]
  sampCountyDat90[1:(1 + includeUrban)] = result$summary.fixed[,5]
  
  # BYM2 hyperparameter phi
  sampCountyDatPar[(2 + clusterEffect + includeUrban)] = result$summary.hyperpar[2,1]
  sampCountyDatSD[(2 + clusterEffect + includeUrban)] = result$summary.hyperpar[2,2]
  sampCountyDat10[(2 + clusterEffect + includeUrban)] = result$summary.hyperpar[2,3]
  sampCountyDat50[(2 + clusterEffect + includeUrban)] = result$summary.hyperpar[2,4]
  sampCountyDat90[(2 + clusterEffect + includeUrban)] = result$summary.hyperpar[2,5]
  
  ## transformed hyperparameters
  # sample the hyperparameters, using the marginals to improve the sampling
  out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
  transformFunction = function(x) {c(1/x[3], x[2], 1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[3]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
  if(!clusterEffect)
    transformFunction = function(x) {c(1/x[1], x[2], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
  transformedOut = apply(out, 1, transformFunction)
  
  # now calculate the summary statistics of the transformed BYM2 hyperparameters
  sampCountyDatPar[(2 + includeUrban):length(parNames)] = rowMeans(transformedOut)
  sampCountyDatSD[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, sd)
  sampCountyDat10[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, quantile, probs=.1)
  sampCountyDat50[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, quantile, probs=.5)
  sampCountyDat90[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, quantile, probs=.9)
  
  if(clusterEffect)
    sampClusterSigmaSRS = 1:Nsim
  
  sortI = matchMultiple(counties, easpc$County)
  clustersPerUrban = easpc$EAUrb[sortI]
  clustersPerRural = easpc$EARur[sortI]
  
  # Simulate from posterior
  if(urbanEffect) {
    samp = inla.posterior.sample(n = Nsim, result = result)
    sampRural = matrix(NA, nrow = 47, ncol = Nsim)
    sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
    sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatUrban = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatRural = matrix(NA, nrow = 47, ncol = Nsim)
    
    for(j in 1:Nsim){
      sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
      sampUrban[, j] = samp[[j]]$latent[1:47]
      if(clusterEffect) {
        # if cluster effect is included, must debias predictions in each modeled strata
        clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
        # muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
        # muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
        # sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
        # sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        sampRuralMod[, j] = sapply(1:length(clustersPerRural), function(i) {mean(expit(sampRural[i, j] + rnorm(clustersPerRural[i], sd=clusterSigma)))})
        sampRuralMod[clustersPerRural == 0,] = 0 # we just need to set these values to something other than NaN
        sampUrbanMod[, j] = sapply(1:length(clustersPerUrban), function(i) {mean(expit(sampUrban[i, j] + rnorm(clustersPerUrban[i], sd=clusterSigma)))})
        sampClusterDatRural[, j] = sampRural[, j] + rnorm(47, sd=clusterSigma)
        sampClusterDatUrban[, j] = sampUrban[, j] + rnorm(47, sd=clusterSigma)
      }
    }
    sampCountyDat = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
    sampCountyDatMod = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
    if(!clusterEffect) {
      sampClusterDatRural = sampRural
      sampClusterDatUrban = sampUrban
    } else {
      sampPixelDatUrbanMod = sampUrbanMod
      sampPixelDatRuralMod = sampRuralMod
    }
    sampPixelDatUrban = sampUrban
    sampPixelDatRural = sampRural
  } else {
    samp = inla.posterior.sample(n = Nsim, result = result)
    sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
    sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
    clustersPerCounty = rowSums(cbind(clustersPerUrban, clustersPerRural))
    for(j in 1:Nsim){
      sampCounty[, j] = samp[[j]]$latent[1:47]
      if(clusterEffect) {
        # if cluster effect is included, must debias predictions in each modeled strata
        clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
        # muSigmaMat = cbind(sampCounty[, j], clusterSigma)
        # sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat
        sampCountyMod[, j] = sapply(1:length(clustersPerCounty), function(i) {mean(expit(sampCounty[i, j] + rnorm(clustersPerCounty[i], sd=clusterSigma)))})
        sampClusterDatRural[, j] = sampCounty[, j] + rnorm(47, sd=clusterSigma)
        sampClusterDatUrban[, j] = sampClusterDatRural[, j]
      }
    }
    sampCountyDat = sampCounty
    sampCountyDatMod = logit(sampCountyMod)
    sampPixelDatUrban = sampCounty
    sampPixelDatRural = sampCounty
    if(!clusterEffect) {
      sampClusterDatRural = sampCounty
      sampClusterDatUrban = sampCounty
    } else {
      sampPixelDatUrbanMod = sampCountyMod
      sampPixelDatRuralMod = sampCountyMod
    }
  }
  
  processSamples = function(samp){
    # 80% credible intervals
    CI = t(apply(X = samp, MARGIN = 1, FUN = quantile, probs = c(0.1, 0.5, 0.9)))
    mm = rowMeans(samp)
    ss = apply(X = samp, MARGIN = 1, FUN = sd)
    
    return(list(logit = list(CI = CI,
                             mean = mm,
                             stddev = ss),
                prob = list(CI = exp(CI))))
  }
  
  ## generate predictions
  Q10 = numeric(47)
  Q50 = numeric(47)
  Q90 = numeric(47)
  mm = numeric(47)
  ss = numeric(47)
  
  tmp = processSamples(sampCountyDat)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDat = list(Q10 = Q10,
                Q50 = Q50,
                Q90 = Q90,
                mean = mm,
                stddev = ss)
  
  if(!is.null(predCountyI))
    resDat = data.frame(resDat)[predCountyI,]
  
  # calculate summary statistics for cluster and pixel predictions
  tmp = processSamples(sampClusterDatRural)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatClusterRural = list(Q10 = Q10,
                            Q50 = Q50,
                            Q90 = Q90,
                            mean = mm,
                            stddev = ss)
  
  if(!is.null(predCountyI))
    resDatClusterRural = data.frame(resDatClusterRural)[predCountyI,]
  
  tmp = processSamples(sampClusterDatUrban)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatClusterUrban = list(Q10 = Q10,
                            Q50 = Q50,
                            Q90 = Q90,
                            mean = mm,
                            stddev = ss)
  
  if(!is.null(predCountyI))
    resDatClusterUrban = data.frame(resDatClusterUrban)[predCountyI,]
  
  tmp = processSamples(sampPixelDatUrban)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatPixelUrban = list(Q10 = Q10,
                          Q50 = Q50,
                          Q90 = Q90,
                          mean = mm,
                          stddev = ss)
  
  if(!is.null(predCountyI))
    resDatPixelUrban = data.frame(resDatPixelUrban)[predCountyI,]
  
  tmp = processSamples(sampPixelDatRural)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatPixelRural = list(Q10 = Q10,
                          Q50 = Q50,
                          Q90 = Q90,
                          mean = mm,
                          stddev = ss)
  
  if(!is.null(predCountyI))
    resDatPixelRural = data.frame(resDatPixelRural)[predCountyI,]
  
  if(clusterEffect) {
    tmp = processSamples(sampCountyDatMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatMod = list(Q10 = Q10,
                     Q50 = Q50,
                     Q90 = Q90,
                     mean = mm,
                     stddev = ss)
    
    if(!is.null(predCountyI))
      resDatMod = data.frame(resDatMod)[predCountyI,]
    
    tmp = processSamples(sampPixelDatUrbanMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatPixelUrbanMod = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
    
    if(!is.null(predCountyI))
      resDatPixelUrbanMod = data.frame(resDatPixelUrbanMod)[predCountyI,]
    
    tmp = processSamples(sampPixelDatRuralMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatPixelRuralMod = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
    
    if(!is.null(predCountyI))
      resDatPixelRuralMod = data.frame(resDatPixelRuralMod)[predCountyI,]
    
  } else {
    resDatMod = NULL
    resDatPixelUrbanMod = NULL
    resDatPixelRuralMod = NULL
  }
  
  
  
  ## now collect the parameters
  # make rural parameter the urban parameter
  # if(includeUrban) {
  #   sampCountyDatPar[2] = -sampCountyDatPar[2]
  #   sampCountyDat10[2] = -sampCountyDat10[2]
  #   sampCountyDat50[2] = -sampCountyDat50[2]
  #   sampCountyDat90[2] = -sampCountyDat90[2]
  # }
  mm = sampCountyDatPar
  ss = sampCountyDatSD
  Q10 = sampCountyDat10
  Q50 = sampCountyDat50
  Q90 = sampCountyDat90
  resDatPar = data.frame(list(mean = mm,
                              stddev = ss, 
                              Q10 = Q10,
                              Q50 = Q50,
                              Q90 = Q90))
  
  # Full result
  designRes = list(predictions = resDat,
                   parameters = resDatPar, 
                   predictionsPixelUrban = resDatPixelUrban, 
                   predictionsPixelRural = resDatPixelRural, 
                   predictionsClusterUrban = resDatClusterUrban, 
                   predictionsClusterRural = resDatClusterRural)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      urbanEffect, '.RData'), designRes = designRes)
  validationText = ""
  if(doValidation)
    validationText = "ValidationFull"
  else if(!is.null(predCountyI))
    validationText = paste0("Validation", predCountyI) # really we are setting saveResults to FALSE in this case
  if(saveResults) {
    save(file = paste0('bym2', fileNameRoot, validationText, 'UrbRur',urbanEffect, 'Cluster', clusterEffect, '.RData'), 
         designRes = designRes)
  }
  
  # include the debiased results if cluster effect is included
  if(clusterEffect) {
    temp = designRes
    designRes = list(predictions = resDatMod,
                     parameters = resDatPar, 
                     predictionsPixelUrban = resDatPixelUrbanMod, 
                     predictionsPixelRural = resDatPixelRuralMod, 
                     predictionsClusterUrban = resDatClusterUrban, 
                     predictionsClusterRural = resDatClusterRural)
    
    if(saveResults) {
      save(file = paste0('bym2', fileNameRoot, validationText, 'UrbRur',urbanEffect, 'Cluster', clusterEffect, 'debiased.RData'), 
           designRes = designRes)
    }
    
    designResMod = designRes
    designRes = temp
  }
  
  # in order to compute DIC or WAIC it is necessary to generate predictive distribution at the 
  # posterior mean for each county. We do it here on a logit scale:
  if(getPosteriorDensity) {
    if(is.null(directLogitEsts))
      stop("Must supply directLogitEsts if getPosteriorDensity == TRUE")
    
    # approximate the posterior county samples using a multivariate gaussian on the logit scale
    theseCountySamples = sampCountyDat
    mu = resDat$mean
    theseCountyResiduals = sweep(theseCountySamples, 1, mu, "-")
    Sigma = (1 / (Nsim - 1)) * theseCountyResiduals %*% t(theseCountyResiduals)
    if(clusterEffect) {
      theseCountySamplesMod = sampCountyDatMod
      muMod = resDatMod$mean
      theseCountyResidualsMod = sweep(theseCountySamplesMod, 1, muMod, "-")
      SigmaMod = (1 / (Nsim - 1)) * theseCountyResidualsMod %*% t(theseCountyResidualsMod)
    }
    
    # calculate posterior density of the direct estimates
    logLik = logLikGP(directLogitEsts - mu, chol(Sigma))
    if(clusterEffect)
      logLikMod = logLikGP(directLogitEsts - muMod, chol(SigmaMod))
    else
      logLikMod = NULL
  }
  
  # compute cluster predictive standard deviation for each strata (square root of the the strata variance plus the cluster variance)
  if(clusterEffect) {
    clustPredSD = sqrt(designRes$predictions$stddev^2 + designRes$parameters$stddev[2 + includeUrban]^2)
    designRes$predictions$clustPredSD = clustPredSD
    designResMod$predictions$clustPredSD = clustPredSD
  }
  
  # compile results
  if(!doValidation)
    out = designRes
  else if(!getPosteriorDensity)
    out = list(designRes=designRes, mod=result)
  else
    out = list(designRes=designRes, mod=result, directLogLik=logLik, directLogLikMod=logLikMod)
  if(clusterEffect)
    out$designResMod = designResMod
  
  out
}