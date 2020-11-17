# this file contains functions for generating lists of command arguments for each model results function, 
# saving the files to be used later

# make the command arguments file for resultsLCPB
# TODO: add in seeds
getSPDE_LCPBCommandArgs = function(gamma=c(-1, 0), rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), 
                                   effRange=400, beta0=c(-3.9, 0), surveyI=1:10, representativeSampling=c(TRUE, FALSE)) {
  
  spde_lcpbCommandArgs = list()
  i = 1
  for(i1 in 1:length(gamma)) {
    thisGamma=gamma[i1]
    
    for(i2 in 1:length(rho)) {
      thisRho = rho[i2]
      
      for(i3 in 1:length(sigmaEpsilon)) {
        thisSigmaEpsilon = sigmaEpsilon[i3]
        
        for(i4 in 1:length(effRange)) {
          thisEffRange = effRange[i4]
          
          for(i5 in 1:length(beta0)) {
            thisBeta0 = beta0[i5]
            
            for(i6 in 1:length(surveyI)) {
              thisSurveyI = surveyI[i6]
              
              for(i7 in 1:length(representativeSampling)) {
                thisRepresentativeSampling = representativeSampling[i7]
                
                spde_lcpbCommandArgs[[i]] = list(gamma=thisGamma, rho=thisRho, sigmaEpsilon=thisSigmaEpsilon, 
                                                 effRange=thisEffRange, beta0=thisBeta0, surveyI=thisSurveyI, 
                                                 representativeSampling=thisRepresentativeSampling)
                i=i+1
              }
            }
          }
        }
      }
    }
  }
  
  save(spde_lcpbCommandArgs, file="savedOutput/simDataSets/spde_lcpbCommandArgs.RData")
}

# make the command arguments file for compareModelsSimulationStudy
getSPDE_LCPBSimStudyCommandArgs = function(gamma=c(-1, 0), rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), 
                                   effRange=400, beta0=c(-3.9, 0), representativeSampling=c(TRUE, FALSE)) {
  
  spde_lcpbSimStudyCommandArgs = list()
  i = 1
  for(i1 in 1:length(gamma)) {
    thisGamma=gamma[i1]
    
    for(i2 in 1:length(rho)) {
      thisRho = rho[i2]
      
      for(i3 in 1:length(sigmaEpsilon)) {
        thisSigmaEpsilon = sigmaEpsilon[i3]
        
        for(i4 in 1:length(effRange)) {
          thisEffRange = effRange[i4]
          
          for(i5 in 1:length(beta0)) {
            thisBeta0 = beta0[i5]
              
              for(i6 in 1:length(representativeSampling)) {
                thisRepresentativeSampling = representativeSampling[i6]
                
                spde_lcpbSimStudyCommandArgs[[i]] = list(gamma=thisGamma, rho=thisRho, sigmaEpsilon=thisSigmaEpsilon, 
                                                 effRange=thisEffRange, beta0=thisBeta0, 
                                                 representativeSampling=thisRepresentativeSampling)
                i=i+1
              }
          }
        }
      }
    }
  }
  
  save(spde_lcpbSimStudyCommandArgs, file="savedOutput/simStudyResults/spde_lcpbSimStudyCommandArgs.RData")
}



# getSPDE_LCPBCommandArgs()







