# scrip for getting all NMR results for the application

getMortResults = function() {
  
  # first for the risk model
  timeSPDE = system.time(resultsSPDE <- fitSPDEKenyaDat(dat=mort, dataType="mort", 
                                                        significanceCI=.8, 
                                                        nPostSamples=1000, verbose=TRUE, seed=seed, 
                                                        urbanEffect=TRUE, clusterEffect=TRUE, 
                                                        leaveOutRegionName=NULL, doValidation=FALSE))[3]
  
  # now aggregate
  timeAllAgg = system.time(agg <- modLCPB(uDraws=resultsSPDE$uDraws, resultsSPDE$sigmaEpsilonDraws, easpa=thiseaspa, 
                                          includeUrban=TRUE, clusterLevel=FALSE, pixelLevel=TRUE, constituencyLevel=TRUE, countyLevel=TRUE, 
                                          regionLevel=TRUE, nationalLevel=TRUE, doModifiedPixelLevel=FALSE, 
                                          onlyDoModifiedPixelLevel=FALSE, 
                                          doLCPb=TRUE, doLCpb=TRUE, doLcpb=TRUE, urbanEffect=resultsSPDE$fixedEffectSummary[2,1]))[3]
  
  # save results
  save(resultsSPDE, timeSPDE, timeAllAgg, file="savedOutput/application/mortResultsRisk.RData")
  save(agg, timeSPDE, timeAllAgg, file="savedOutput/application/mortResultsAgg.RData")
  
  invisible(NULL)
}

# Make plots for the neonatal mortality application
makeMortPlots = function() {
  # first load the model predictions
  out = load("savedOutput/application/mortResultsAgg.RData")
  
  # calculate the range of predictions and CI widths
  rangePredPixel = c()
  rangePredConstituency = c()
  rangePredCounty = c()
  rangePredProvince = c()
  
  # make the color scales
}



