
# estimate of var(p_{emp}(A))
# 
varPrevEmp = function(etaSDraws, MDraws=300, NDraws=50*MDraws, sigmaEps, varianceRiskSL=NULL, q=NULL) {
  if(!is.matrix(etaSDraws)) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  nDraws = ncol(etaSDraws)
  
  # reformat/convert from fixed parameters to draws if need be
  if(length(MDraws) == 1) {
    MDraws = rep(MDraws, nDraws)
  }
  if(length(NDraws) == 1) {
    NDraws = matrix(rep(1/MDraws, each=nrow(etaSDraws)), ncol=nDraws)
  } else if(is.matrix(NDraws) && (ncol(NDraws) == 1)) {
    NDraws = matrix(rep(NDraws, nDraws), ncol=nDraws)
  }
  
  if(is.null(varianceRiskSL)) {
    varianceRiskSL = varRiskSL(etaSDraws, sigmaEps, q)
  }
  
  ENiSq = mean(5)
}

# simplified estimate of var(p_{emp}(A)) assuming M and N are fixed
# 
varPrevEmpSimple = function(etaSDraws, N, M, sigmaEps, 
                            varianceRiskSL=NULL, 
                            expectRiskSL2A=NULL, expectRiskSLA2=NULL, 
                            expectRiskSLA=NULL, q=NULL) {
  if(!is.matrix(etaSDraws)) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  nDraws = ncol(etaSDraws)
  
  if(is.null(expectRiskSLA)) {
    expectRiskSLA = ERiskSLA(etaSDraws, sigmaEps, q)
  }
  if(is.null(expectRiskSLA2)) {
    expectRiskSLA2 = ERiskSLA2(etaSDraws, sigmaEps, q)
  }
  if(is.null(varianceRiskSL)) {
    varianceRiskSL = varRiskSL(etaSDraws, sigmaEps, q)
  }
  if(is.null(expectRiskSL2A)) {
    expectRiskSL2A = ERiskSL2A(etaSDraws, sigmaEps, q)
  }
  
  E1 = (1/M - 1/N) * expectRiskSL2A
  E2 = 1/M * expectRiskSLA2
  E3 = 1/N * expectRiskSLA
  E1 + E2 + E3 + varianceRiskSL
}

ERiskSLA = function(etaSDraws, sigmaEps, q, smoothRisk=NULL) {
  if(is.null(smoothRisk)) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws))), 
                                              logisticApprox=FALSE), nrow=nrow(etaSDraws))
  }
  
  # get expected numerator and denominators
  nSamplesSmoothRisk = matrix(rep(q, nDraws), ncol=nDraws)
  zSamplesSmoothRisk = sweep(smoothRisk, 1, q, "*")
  sum(zSamplesSmoothRisk)/sum(nSamplesSmoothRisk)
}

# estimate of var(r_{SL}(A))
# est_s is the "smooth" component of eta, the linear predictor. It is everything 
# in the linear predictor but the epsilon. sigmaEps is the SD of the nugget, and 
# q is the normalized target population density for each point for which etaS is 
# drawn
varRiskSL = function(etaSDraws, sigmaEps, q, smoothRisk=NULL) {
  if(!is.matrix(etaSDraws)) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  nDraws = ncol(etaSDraws)
  
  if(is.null(smoothRisk)) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws))), 
                                              logisticApprox=FALSE), nrow=nrow(etaSDraws))
  }
  
  # get expected numerator and denominators
  nSamplesSmoothRisk = matrix(rep(q, nDraws), ncol=nDraws)
  zSamplesSmoothRisk = sweep(smoothRisk, 1, q, "*")
  riskSL = colsums(zSamplesSmoothRisk)/colsums(nSamplesSmoothRisk)
  var(riskSL)
}

# estimate of E[r_{SL}(A)^2]
# est_s is the "smooth" component of eta, the linear predictor. It is everything 
# in the linear predictor but the epsilon. sigmaEps is the SD of the nugget, and 
# q is the normalized target population density for each point for which etaS is 
# drawn
ERiskSLA2 = function(etaSDraws, sigmaEps, q, smoothRisk=NULL) {
  if(!is.matrix(etaSDraws)) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  nDraws = ncol(etaSDraws)
  
  if(is.null(smoothRisk)) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws))), 
                                              logisticApprox=FALSE), nrow=nrow(etaSDraws))
  }
  
  # get expected numerator and denominators
  nSamplesSmoothRisk = matrix(rep(q, nDraws), ncol=nDraws)
  zSamplesSmoothRisk = sweep(smoothRisk, 1, q, "*")
  riskSL = colsums(zSamplesSmoothRisk)/colsums(nSamplesSmoothRisk)
  mean(riskSL^2)
}

ERiskSL2A = function(etaSDraws, sigmaEps, q, smoothRisk=NULL) {
  if(!is.matrix(etaSDraws)) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  nDraws = ncol(etaSDraws)
  
  if(is.null(smoothRisk)) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws))), 
                                              logisticApprox=FALSE), nrow=nrow(etaSDraws))
  }
  
  # get expected numerator and denominators
  nSamplesSmoothRisk = matrix(rep(q, nDraws), ncol=nDraws)
  zSamplesSmoothRisk = sweep(smoothRisk, 1, q, "*")
  riskSL = colsums(zSamplesSmoothRisk)/colsums(nSamplesSmoothRisk)
  mean(riskSL^2)
}

testVarPrevEmp = function(M=300) {
  
}

varRiskSLStrat = function(smoothRiskDraws, q=NULL) {
  if(is.null(q)) {
    q = rep(1, nrow(smoothRiskDraws))
  }
  if(sum(q) != 1) {
    q = q*(1/sum(q))
  }
  
  # get urban and rural areal smooth risks
  smoothRiskDrawsAreal = t(smoothRiskDraws) %*% q
  
  var(smoothRiskDrawsAreal)
}

varPrevEmpStrat = function(etaSDraws=NULL, sigmaEps=NULL, Murb=160, Mrur=160, Nurb=37*Murb, Nrur=37*Mrur, 
                           q=NULL, smoothRiskDraws=NULL, smoothRiskSqDraws=NULL, urbVec=c(TRUE, FALSE), 
                           Mprop=NULL, returnVarRSL=TRUE) {
  if(!is.null(etaSDraws) && (!is.matrix(etaSDraws))) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  if(!is.null(smoothRiskDraws) && (!is.matrix(smoothRiskDraws))) {
    smoothRiskDraws = matrix(smoothRiskDraws, nrow=1)
  }
  if(!is.null(smoothRiskSqDraws) && (!is.matrix(smoothRiskSqDraws))) {
    smoothRiskSqDraws = matrix(smoothRiskSqDraws, nrow=1)
  }
  
  N = Nurb + Nrur
  M = Murb + Mrur
  
  if(is.null(smoothRiskDraws)) {
    if(is.null(etaSDraws) || is.null(sigmaEps)) {
      stop("must either provide etaSDraws and sigmaEps or smoothRiskDraws")
    }
    
    # integrate out the cluster effect in the risk
    smoothRiskDraws = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws))), 
                                                   logisticApprox=FALSE), nrow=nrow(etaSDraws))
  }
  
  if(is.null(smoothRiskSqDraws)) {
    if(is.null(etaSDraws) || is.null(sigmaEps)) {
      stop("must either provide etaSDraws and sigmaEps or smoothRiskSqDraws")
    }
    
    # integrate out the cluster effect in the squared risk
    smoothRiskSqDraws = matrix(logitNormSqMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws)))), nrow=nrow(etaSDraws))
  }
  
  if(is.null(q)) {
    q = rep(1, length(urbVec))
    q[urbVec] = Nurb / sum(q[urbVec])
    q[!urbVec] = Nrur / sum(q[!urbVec])
  }
  if(sum(q) != 1) {
    q = q*(1/sum(q))
  }
  qUrb = q[urbVec]
  qRur = q[!urbVec]
  Qurb = sum(qUrb)
  Qrur = sum(qRur)
  qUrb = qUrb * (1/sum(qUrb))
  qRur = qRur * (1/sum(qRur))
  
  # get stratum and areal smooth risks
  nUrbanPixels = sum(urbVec)
  nRuralPixels = sum(!urbVec)
  smoothRiskDrawsUrb = matrix(smoothRiskDraws[urbVec,], nrow=nUrbanPixels)
  smoothRiskDrawsRur = matrix(smoothRiskDraws[!urbVec,], nrow=nRuralPixels)
  smoothRiskDrawsUrbAreal = t(smoothRiskDrawsUrb) %*% qUrb
  smoothRiskDrawsRurAreal = t(smoothRiskDrawsRur) %*% qRur
  if(nRuralPixels == 0) {
    smoothRiskDrawsAreal = smoothRiskDrawsUrbAreal
  } else if(nUrbanPixels == 0) {
    smoothRiskDrawsAreal = smoothRiskDrawsRurAreal
  } else {
    smoothRiskDrawsAreal = smoothRiskDrawsUrbAreal * Qurb + smoothRiskDrawsRurAreal * Qrur
  }
  
  # do the same for squared smooth risks
  smoothRiskSqDrawsUrb = matrix(smoothRiskSqDraws[urbVec,], nrow=nUrbanPixels)
  smoothRiskSqDrawsRur = matrix(smoothRiskSqDraws[!urbVec,], nrow=nRuralPixels)
  if(nUrbanPixels != 0) {
    smoothRisk2DrawsUrbAreal = t(smoothRiskSqDrawsUrb) %*% qUrb
  } else {
    smoothRisk2DrawsUrbAreal = numeric(0)
  }
  if(nRuralPixels != 0) {
    smoothRisk2DrawsRurAreal = t(smoothRiskSqDrawsRur) %*% qRur
  } else {
    smoothRisk2DrawsRurAreal = numeric(0)
  }
  
  if((nUrbanPixels != 0) && (nRuralPixels != 0)) {
    # estimate moments of urban and rural smooth risk
    varUrb = var(smoothRiskDrawsUrbAreal)
    varRur = var(smoothRiskDrawsRurAreal)
    varA = var(smoothRiskDrawsAreal)
    covUrbRur = cov(smoothRiskDrawsUrbAreal, smoothRiskDrawsRurAreal)
    ERslUrb = mean(smoothRiskDrawsUrbAreal)
    ERslRur = mean(smoothRiskDrawsRurAreal)
    ERslUrb2 = mean(smoothRiskDrawsUrbAreal^2)
    ERslRur2 = mean(smoothRiskDrawsRurAreal^2)
    ERslUrbRur = mean(smoothRiskDrawsUrbAreal * smoothRiskDrawsRurAreal)
    ERsl2Urb = mean(smoothRisk2DrawsUrbAreal)
    ERsl2Rur = mean(smoothRisk2DrawsRurAreal)
    
    # estimate variances and covariances of the "N" terms
    varBinProp = Qurb * (1 - Qurb) / N
    covBinPropUrbRur = -varBinProp
  } else if(nUrbanPixels != 0) {
    # estimate moments of urban and rural smooth risk
    varUrb = var(smoothRiskDrawsUrbAreal)
    varRur = 0
    varA = var(smoothRiskDrawsAreal)
    covUrbRur = 0
    ERslUrb = mean(smoothRiskDrawsUrbAreal)
    ERslRur = 0
    ERslUrb2 = mean(smoothRiskDrawsUrbAreal^2)
    ERslRur2 = 0
    ERslUrbRur = 0
    ERsl2Urb = mean(smoothRisk2DrawsUrbAreal)
    ERsl2Rur = 0
    
    # estimate variances and covariances of the "N" terms
    varBinProp = 0
    covBinPropUrbRur = 0
    
    # set to 1 to avoid numerical problems dividing by 0
    Mrur = 1
  } else if(nRuralPixels != 0) {
    # estimate moments of urban and rural smooth risk
    varUrb = 0
    varRur = var(smoothRiskDrawsRurAreal)
    varA = var(smoothRiskDrawsAreal)
    covUrbRur = 0
    ERslUrb = 0
    ERslRur = mean(smoothRiskDrawsRurAreal)
    ERslUrb2 = 0
    ERslRur2 = mean(smoothRiskDrawsRurAreal^2)
    ERslUrbRur = 0
    ERsl2Urb = 0
    ERsl2Rur = mean(smoothRisk2DrawsRurAreal)
    
    # estimate variances and covariances of the "N" terms
    varBinProp = 0
    covBinPropUrbRur = 0
    
    # set to 1 to avoid numerical problems dividing by 0
    Murb = 1
  }
  
  # calculate var(E[pemp(A) | Mvec, Nvec, u])
  varEterm = ERslUrb * varBinProp + 
    ERslRur * varBinProp + 
    ERslUrbRur * covBinPropUrbRur + varA
  
  # calculate E[var(pemp(A) | Mvec, Nvec, u)]
  EvarTerm = Qurb/N * (ERslUrb - ERsl2Urb) + Qrur/N * (ERslRur - ERsl2Rur) + 
    Qurb^2/Murb * (ERsl2Urb - ERslUrb2) + Qrur^2/Mrur * (ERsl2Rur - ERslRur2)
  
  if(!is.null(Mprop)) {
    # Mprop is the probability of an EA being in this area
    MurbProb = Mprop * Murb / M
    MrurProb = Mprop * Mrur / M
    expectInversePosBinom()
    EvarTerm2 = Qurb/N * (ERslUrb - ERsl2Urb) + Qrur/N * (ERslRur - ERsl2Rur) + 
      Qurb^2/Murb * (ERsl2Urb - ERslUrb2) + Qrur^2/Mrur * (ERsl2Rur - ERslRur2)
  }
  
  if(!returnVarRSL) {
    varEterm + EvarTerm
  } else {
    c(prevEmp=varEterm + EvarTerm, riskSL=varA)
  }
}

varBurdEmpStrat = function(etaSDraws=NULL, sigmaEps=NULL, Murb=160, Mrur=160, Nurb=37*Murb, Nrur=37*Mrur, 
                           q=NULL, smoothRiskDraws=NULL, smoothRiskSqDraws=NULL, urbVec=c(TRUE, FALSE), returnVarBSL=TRUE) {
  if(!is.null(etaSDraws) && (!is.matrix(etaSDraws))) {
    etaSDraws = matrix(etaSDraws, nrow=1)
  }
  if(!is.null(smoothRiskDraws) && (!is.matrix(smoothRiskDraws))) {
    smoothRiskDraws = matrix(smoothRiskDraws, nrow=1)
  }
  if(!is.null(smoothRiskSqDraws) && (!is.matrix(smoothRiskSqDraws))) {
    smoothRiskSqDraws = matrix(smoothRiskSqDraws, nrow=1)
  }
  
  N = Nurb + Nrur
  M = Murb + Mrur
  
  if(is.null(smoothRiskDraws)) {
    if(is.null(etaSDraws) || is.null(sigmaEps)) {
      stop("must either provide etaSDraws and sigmaEps or smoothRiskDraws")
    }
    
    # integrate out the cluster effect in the risk
    smoothRiskDraws = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws))), 
                                                   logisticApprox=FALSE), nrow=nrow(etaSDraws))
  }
  
  if(is.null(smoothRiskSqDraws)) {
    if(is.null(etaSDraws) || is.null(sigmaEps)) {
      stop("must either provide etaSDraws and sigmaEps or smoothRiskSqDraws")
    }
    
    # integrate out the cluster effect in the squared risk
    smoothRiskSqDraws = matrix(logitNormSqMean(cbind(c(as.matrix(etaSDraws)), rep(sigmaEps, length(etaSDraws)))), nrow=nrow(etaSDraws))
  }
  
  if(is.null(q)) {
    q = rep(1, length(urbVec))
    q[urbVec] = Nurb / sum(q[urbVec])
    q[!urbVec] = Nrur / sum(q[!urbVec])
  }
  if(sum(q) != 1) {
    q = q*(1/sum(q))
  }
  qUrb = q[urbVec]
  qRur = q[!urbVec]
  Qurb = sum(qUrb)
  Qrur = sum(qRur)
  qUrb = qUrb * (1/sum(qUrb))
  qRur = qRur * (1/sum(qRur))
  
  # get stratum and areal smooth risks
  nUrbanPixels = sum(urbVec)
  nRuralPixels = sum(!urbVec)
  smoothRiskDrawsUrb = matrix(smoothRiskDraws[urbVec,], nrow=nUrbanPixels)
  smoothRiskDrawsRur = matrix(smoothRiskDraws[!urbVec,], nrow=nRuralPixels)
  smoothRiskDrawsUrbAreal = t(smoothRiskDrawsUrb) %*% qUrb
  smoothRiskDrawsRurAreal = t(smoothRiskDrawsRur) %*% qRur
  if(nRuralPixels == 0) {
    smoothRiskDrawsAreal = smoothRiskDrawsUrbAreal
  } else if(nUrbanPixels == 0) {
    smoothRiskDrawsAreal = smoothRiskDrawsRurAreal
  } else {
    smoothRiskDrawsAreal = smoothRiskDrawsUrbAreal * Qurb + smoothRiskDrawsRurAreal * Qrur
  }
  
  # do the same for squared smooth risks
  smoothRiskSqDrawsUrb = matrix(smoothRiskSqDraws[urbVec,], nrow=nUrbanPixels)
  smoothRiskSqDrawsRur = matrix(smoothRiskSqDraws[!urbVec,], nrow=nRuralPixels)
  if(nUrbanPixels != 0) {
    smoothRisk2DrawsUrbAreal = t(smoothRiskSqDrawsUrb) %*% qUrb
  } else {
    smoothRisk2DrawsUrbAreal = numeric(0)
  }
  if(nRuralPixels != 0) {
    smoothRisk2DrawsRurAreal = t(smoothRiskSqDrawsRur) %*% qRur
  } else {
    smoothRisk2DrawsRurAreal = numeric(0)
  }
  
  if((nUrbanPixels != 0) && (nRuralPixels != 0)) {
    # estimate moments of urban and rural smooth risk
    varUrb = var(smoothRiskDrawsUrbAreal)
    varRur = var(smoothRiskDrawsRurAreal)
    varBurdenA = var(smoothRiskDrawsUrbAreal * Nurb + smoothRiskDrawsRurAreal * Nrur)
    covUrbRur = cov(smoothRiskDrawsUrbAreal, smoothRiskDrawsRurAreal)
    ERslUrb = mean(smoothRiskDrawsUrbAreal)
    ERslRur = mean(smoothRiskDrawsRurAreal)
    ERslUrb2 = mean(smoothRiskDrawsUrbAreal^2)
    ERslRur2 = mean(smoothRiskDrawsRurAreal^2)
    ERslUrbRur = mean(smoothRiskDrawsUrbAreal * smoothRiskDrawsRurAreal)
    ERsl2Urb = mean(smoothRisk2DrawsUrbAreal)
    ERsl2Rur = mean(smoothRisk2DrawsRurAreal)
    
    # estimate variances and covariances of the "N" terms
    varBin = Qurb * (1 - Qurb) * N
    covBinUrbRur = -varBin
  } else if(nUrbanPixels != 0) {
    varUrb = var(smoothRiskDrawsUrbAreal)
    varRur = 0
    varBurdenA = var(smoothRiskDrawsUrbAreal * Nurb)
    covUrbRur = 0
    ERslUrb = mean(smoothRiskDrawsUrbAreal)
    ERslRur = 0
    ERslUrb2 = mean(smoothRiskDrawsUrbAreal^2)
    ERslRur2 = 0
    ERslUrbRur = 0
    ERsl2Urb = mean(smoothRisk2DrawsUrbAreal)
    ERsl2Rur = 0
    
    # estimate variances and covariances of the "N" terms
    varBin = 0
    covBinUrbRur = 0
    
    # set to 1 to avoid numerical problems dividing by 0
    Mrur = 1
  } else if(nRuralPixels != 0) {
    varUrb = 0
    varRur = var(smoothRiskDrawsRurAreal)
    varBurdenA = var(smoothRiskDrawsRurAreal * Nrur)
    covUrbRur = 0
    ERslUrb = 0
    ERslRur = mean(smoothRiskDrawsRurAreal)
    ERslUrb2 = 0
    ERslRur2 = mean(smoothRiskDrawsRurAreal^2)
    ERslUrbRur = 0
    ERsl2Urb = 0
    ERsl2Rur = mean(smoothRisk2DrawsRurAreal)
    
    # estimate variances and covariances of the "N" terms
    varBin = 0
    covBinUrbRur = 0
    
    # set to 1 to avoid numerical problems dividing by 0
    Murb = 1
  }
  
  # calculate var(E[bemp(A) | Mvec, Nvec, u])
  varEterm = ERslUrb * varBin + ERslRur * varBin + 
    ERslUrbRur * covBinUrbRur + varBurdenA
  
  # calculate E[var(bemp(A) | Mvec, Nvec, u)]
  EvarTerm = Nurb * (ERslUrb - ERsl2Urb) + Nrur * (ERslRur - ERsl2Rur) + 
    Nurb^2/Murb * (ERsl2Urb - ERslUrb2) + Nrur^2/Mrur * (ERsl2Rur - ERslRur2)
  
  if(!returnVarBSL) {
    varEterm + EvarTerm
  } else {
    c(burdEmp=varEterm + EvarTerm, burdSL=varBurdenA)
  }
}

expectInversePosBinom = function(n, p) {
  probNonzero = 1 - dbinom(0, n, p)
  norm = 1/probNonzero
  norm * sum(dbinom(1:n, n, p) * (1/(1:n)))
}

testVarPrevEmpStrat = function(etaSDraws=NULL, sigmaEps=sqrt(0.4403118), Murb=100.71360, Mrur=163.35209, Nurb=4353.257, Nrur=10943.831, 
                               q=NULL, smoothRiskDraws=NULL, urbVec=c(TRUE, FALSE)) {
  # Use Ainabkoi, Uasin Gishu as an example for N, M:
  # head(poppconAdjusted)
  #        subarea        area   popUrb    popRur  popTotal
  # 1          805     Baringo    0.000  2430.858  2430.858
  # 2     Ainabkoi Uasin Gishu 4353.257 10943.831 15297.087
  # test = meanEAsPerCon2()
  # head(test)
  #        subarea        area   popUrb    popRur  popTotal   pctUrb   pctTotal meanUrbanEAs meanRuralEAs meanTotalEAs
  # 1          805     Baringo     0.00  20238.42  20238.42  0.00000 0.04706609      0.00000     64.87799     64.87799
  # 2     Ainabkoi Uasin Gishu 38372.38  90733.00 129105.37 29.72175 0.30024506    100.71360    163.35209    264.06570
  
  if(is.null(etaSDraws)) {
    covMat = matrix(c(1, .9, .9, 1), nrow=2) * 0.2012871
    L = t(chol(covMat))
    NDraws = 1000
    etaSDraws = matrix(-4.128279711 + c(0, -0.002455231) + L %*% matrix(rnorm(2*NDraws), ncol=NDraws), ncol=NDraws)
  }
  
}






