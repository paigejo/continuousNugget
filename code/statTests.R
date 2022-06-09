
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
                           q=NULL, smoothRiskDraws=NULL, urbVec=c(TRUE, FALSE), returnVarRSL=TRUE) {
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
  smoothRiskDrawsUrb = smoothRiskDraws[urbVec,]
  smoothRiskDrawsRur = smoothRiskDraws[!urbVec,]
  smoothRiskDrawsUrbAreal = t(smoothRiskDrawsUrb) %*% qUrb
  smoothRiskDrawsRurAreal = t(smoothRiskDrawsRur) %*% qRur
  smoothRiskDrawsAreal = smoothRiskDrawsUrbAreal * Qurb + smoothRiskDrawsRurAreal * Qrur
  
  # do the same for squared smooth risks
  smoothRisk2DrawsUrbAreal = t(smoothRiskDrawsUrb^2) %*% qUrb
  smoothRisk2DrawsRurAreal = t(smoothRiskDrawsRur^2) %*% qRur
  
  # estimate moments of urban and rural smooth risk
  varUrb = var(smoothRiskDrawsUrbAreal)
  varRur = var(smoothRiskDrawsRurAreal)
  varA = var(smoothRiskDrawsAreal)
  covUrbRur = cov(smoothRiskDrawsUrbAreal, smoothRiskDrawsRurAreal)
  ERslUrb = mean(smoothRiskDrawsUrbAreal)
  ERslRur = mean(smoothRiskDrawsRurAreal)
  ERslUrbRur = mean(smoothRiskDrawsUrbAreal * smoothRiskDrawsRurAreal)
  ERsl2Urb = mean(smoothRisk2DrawsUrbAreal)
  ERsl2Rur = mean(smoothRisk2DrawsRurAreal)
  
  # estimate variances and covariances of the "N" terms
  varBinProp = Qurb * (1 - Qurb) / N
  covBinPropUrbRur = -varBinProp
  
  # calculate var(E[pemp(A) | Mvec, Nvec, u])
  varEterm = ERslUrb * varBinProp + 
    ERslRur * varBinProp + 
    ERslUrbRur * covBinPropUrbRur + varA
  
  # calculate E[var(pemp(A) | Mvec, Nvec, u)]
  EvarTerm = Qurb * (ERslUrb - ERsl2Urb) + Qrur * (ERslRur - ERsl2Rur)
  
  if(!returnVarRSL) {
    varEterm + EvarTerm
  } else {
    c(prevEmp=varEterm + EvarTerm, riskSL=varA)
  }
}

varBurdEmpStrat = function(etaSDraws=NULL, sigmaEps=NULL, Murb=160, Mrur=160, Nurb=37*Murb, Nrur=37*Mrur, 
                           q=NULL, smoothRiskDraws=NULL, urbVec=c(TRUE, FALSE), returnVarBSL=TRUE) {
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
  
  # get stratum and areal smooth risks
  smoothRiskDrawsUrbAreal = t(smoothRiskDrawsUrb) %*% qUrb
  smoothRiskDrawsRurAreal = t(smoothRiskDrawsRur) %*% qRur
  smoothBurdenDrawsAreal = smoothRiskDrawsUrbAreal * Nurb + smoothRiskDrawsRurAreal * Nrur
  
  # do the same for squared smooth risks
  smoothRisk2DrawsUrbAreal = t(smoothRiskDrawsUrb^2) %*% qUrb
  smoothRisk2DrawsRurAreal = t(smoothRiskDrawsRur^2) %*% qRur
  
  # estimate moments of urban and rural smooth risk
  varUrb = var(smoothRiskDrawsUrbAreal)
  varRur = var(smoothRiskDrawsRurAreal)
  varBurdenA = var(smoothBurdenDrawsAreal)
  covUrbRur = cov(smoothRiskDrawsUrbAreal, smoothRiskDrawsRurAreal)
  ERslUrb = mean(smoothRiskDrawsUrbAreal)
  ERslRur = mean(smoothRiskDrawsRurAreal)
  ERslUrbRur = mean(smoothRiskDrawsUrbAreal * smoothRiskDrawsRurAreal)
  ERsl2Urb = mean(smoothRisk2DrawsUrbAreal)
  ERsl2Rur = mean(smoothRisk2DrawsRurAreal)
  
  # estimate variances and covariances of the "N" terms
  varBin = Qurb * (1 - Qurb) * N
  covBinUrbRur = -varBin
  
  # calculate var(E[pemp(A) | Mvec, Nvec, u])
  varEterm = ERslUrb * varBin + 
    ERslRur * varBin + 
    ERslUrbRur * covBinUrbRur + varBurdenA
  
  # calculate E[var(pemp(A) | Mvec, Nvec, u)]
  EvarTerm = Nurb * (ERslUrb - ERsl2Urb) + Nrur * (ERslRur - ERsl2Rur)
  
  if(!returnVarBSL) {
    varEterm + EvarTerm
  } else {
    c(burdEmp=varEterm + EvarTerm, burdSL=varBurdenA)
  }
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






