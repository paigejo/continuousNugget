# script for generating tables of results

generateValidationTables = function(dat=NULL, dataType=c("mort", "ed"), 
                                    significanceCI=.8, urbanEffect=TRUE, clusterEffect=TRUE, kmres=5, 
                                    stratifiedValidation=TRUE) {
  
  # load observations
  dataType = match.arg(dataType)
  if(is.null(dat)) {
    if(dataType == "mort") {
      dat = mort
      dataType2 = "children"
    }
    else {
      dat = ed
      dataType2 = "women"
    }
  }
  if(dataType == "mort")
    fileNameRoot = "Mort"
  else
    fileNameRoot = "Ed"
  
  # first get the full model/in sample results
  fileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationFull", 
                    "_urb", urbanEffect, ".RData")
  out = load(fileName)
  # truthTableFull=truthTableFull, truthTableFullCounty=truthTableFullCounty, truthTableFullProvince=truthTableFullProvince, 
  # timePoint=timePoint, timeAggregationPixel=timeAggregationPixel, timeAggregationCountyProvince=timeAggregationCountyProvince, 
  # fit=fit, fit2=fit2, fit3=fit3, 
  # fullScoresPixel, fullScoresCounty, fullScoresProvince, 
  # file=fileName
  # fullScoresProvince: 
  # fullPooledScoreslcpb=fullPooledScoresProvincelcpb, fullRuralScoreslcpb=fullRuralScoresProvincelcpb, fullUrbanScoreslcpb=fullUrbanScoresProvincelcpb, 
  # fullPooledScoresLcpb=fullPooledScoresProvinceLcpb, fullRuralScoresLcpb=fullRuralScoresProvinceLcpb, fullUrbanScoresLcpb=fullUrbanScoresProvinceLcpb, 
  # fullPooledScoresLCpb=fullPooledScoresProvinceLCpb, fullRuralScoresLCpb=fullRuralScoresProvinceLCpb, fullUrbanScoresLCpb=fullUrbanScoresProvinceLCpb, 
  # fullPooledScoresLCPb=fullPooledScoresProvinceLCPb, fullRuralScoresLCPb=fullRuralScoresProvinceLCPb, fullUrbanScoresLCPb=fullUrbanScoresProvinceLCPb, 
  # fullPooledScoresLCPB=fullPooledScoresProvinceLCPB, fullRuralScoresLCPB=fullRuralScoresProvinceLCPB, fullUrbanScoresLCPB=fullUrbanScoresProvinceLCPB
  
  browser()
  
  ## Pixel level
  fullScoresPixel$fullPooledScoresLCPB$timeAggregationCountyProvince.elapsed = fullScoresPixel$fullPooledScoreslcpb$timeAggregationCountyProvince.elapsed
  allPixelLevelResults = data.frame(rbind(colMeans(fullScoresPixel$fullPooledScoreslcpb), 
                                          colMeans(fullScoresPixel$fullPooledScoresLcpb), 
                                          colMeans(fullScoresPixel$fullPooledScoresLCpb), 
                                          colMeans(fullScoresPixel$fullPooledScoresLCPb), 
                                          colMeans(fullScoresPixel$fullPooledScoresLCPB)))
  rownames(allPixelLevelResults) = c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB")
  allPixelLevelResults$Coverage = 100 * allPixelLevelResults$Coverage
  allPixelLevelResults$Bias = 100 * allPixelLevelResults$Bias
  allPixelLevelResults$RMSE = 100 * allPixelLevelResults$RMSE
  allPixelLevelResults$Width = 100 * allPixelLevelResults$Width
  allPixelLevelResults$CRPS = 100 * allPixelLevelResults$CRPS
  allPixelLevelResults[,c(11:13)] = allPixelLevelResults[,c(11:13)] / 60
  allPixelLevelResults = allPixelLevelResults[,-c(2, 3, 8:10)]
  xtable(allPixelLevelResults)
  
  # County level
  # aggregated
  allCountyLevelResults = data.frame(rbind(colMeans(fullScoresCounty$fullPooledScoreslcpb), 
                                          colMeans(fullScoresCounty$fullPooledScoresLcpb), 
                                          colMeans(fullScoresCounty$fullPooledScoresLCpb), 
                                          colMeans(fullScoresCounty$fullPooledScoresLCPb), 
                                          colMeans(fullScoresCounty$fullPooledScoresLCPB)))
  rownames(allCountyLevelResults) = c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB")
  allCountyLevelResults$Coverage = 100 * allCountyLevelResults$Coverage
  allCountyLevelResults$Bias = 100 * allCountyLevelResults$Bias
  allCountyLevelResults$RMSE = 100 * allCountyLevelResults$RMSE
  allCountyLevelResults$Width = 100 * allCountyLevelResults$Width
  allCountyLevelResults$CRPS = 100 * allCountyLevelResults$CRPS
  allCountyLevelResults[,c(11:13)] = allCountyLevelResults[,c(11:13)] / 60
  allCountyLevelResults = allCountyLevelResults[,-c(2, 3, 8:10)]
  xtable(allCountyLevelResults)
  
  # rural
  allCountyLevelResultsRural = data.frame(rbind(colMeans(fullScoresCounty$fullRuralScoreslcpb), 
                                           colMeans(fullScoresCounty$fullRuralScoresLcpb), 
                                           colMeans(fullScoresCounty$fullRuralScoresLCpb), 
                                           colMeans(fullScoresCounty$fullRuralScoresLCPb), 
                                           colMeans(fullScoresCounty$fullRuralScoresLCPB)))
  rownames(allCountyLevelResultsRural) = c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB")
  allCountyLevelResultsRural$Coverage = 100 * allCountyLevelResultsRural$Coverage
  allCountyLevelResultsRural$Bias = 100 * allCountyLevelResultsRural$Bias
  allCountyLevelResultsRural$RMSE = 100 * allCountyLevelResultsRural$RMSE
  allCountyLevelResultsRural$Width = 100 * allCountyLevelResultsRural$Width
  allCountyLevelResultsRural$CRPS = 100 * allCountyLevelResultsRural$CRPS
  allCountyLevelResultsRural[,c(11:13)] = allCountyLevelResultsRural[,c(11:13)] / 60
  allCountyLevelResultsRural = allCountyLevelResultsRural[,-c(2, 3, 8:10)]
  xtable(allCountyLevelResultsRural)
  
  # urban
  allCountyLevelResultsUrban = data.frame(rbind(colMeans(fullScoresCounty$fullUrbanScoreslcpb), 
                                                colMeans(fullScoresCounty$fullUrbanScoresLcpb), 
                                                colMeans(fullScoresCounty$fullUrbanScoresLCpb), 
                                                colMeans(fullScoresCounty$fullUrbanScoresLCPb), 
                                                colMeans(fullScoresCounty$fullUrbanScoresLCPB)))
  rownames(allCountyLevelResultsUrban) = c("lcpb", "Lcpb", "LCpb", "LCPb", "LCPB")
  allCountyLevelResultsUrban$Coverage = 100 * allCountyLevelResultsUrban$Coverage
  allCountyLevelResultsUrban$Bias = 100 * allCountyLevelResultsUrban$Bias
  allCountyLevelResultsUrban$RMSE = 100 * allCountyLevelResultsUrban$RMSE
  allCountyLevelResultsUrban$Width = 100 * allCountyLevelResultsUrban$Width
  allCountyLevelResultsUrban$CRPS = 100 * allCountyLevelResultsUrban$CRPS
  allCountyLevelResultsUrban[,c(11:13)] = allCountyLevelResultsUrban[,c(11:13)] / 60
  allCountyLevelResultsUrban = allCountyLevelResultsUrban[,-c(2, 3, 8:10)]
  xtable(allCountyLevelResultsUrban)
  
  # now get all of the left out results
  fileName = paste0("savedOutput/validation/resultsSPDE", fileNameRoot, "ValidationAllTemp", 
                    "_urb", urbanEffect, ".RData")
  out2 = load(fileName)
  # completeScoreTable, pooledScoreTable, ruralScoreTable, urbanScoreTable, 
  # binnedScoringRulesuuAll, binnedScoringRulesuUAll, binnedScoringRulesUuAll, binnedScoringRulesUUAll, 
  # binnedScoringRulesAuAll, binnedScoringRulesAUAll, binnedScoringRulesuAAll, binnedScoringRulesUAAll, 
  # binnedScoringRulesAAAll, 
  # singleScores, i
  
  allPixelLevelResultsLCPB = data.frame(rbind(colMeans(pooledScoreTable), colMeans(urbanScoreTable), colMeans(ruralScoreTable))[,-1])
  rownames(allPixelLevelResultsLCPB) = c("All", "Urban", "Rural")
  allPixelLevelResultsLCPB$Coverage = 100 * allPixelLevelResultsLCPB$Coverage
  allPixelLevelResultsLCPB[,1:5] = allPixelLevelResultsLCPB[,1:5] * 10^2
  
  xtable(allPixelLevelResultsLCPB, digits=3)
}

makeAllSimulationStudyTables = function() {
  out = load("savedOutput/simStudyResults/spde_lcpbSimStudyCommandArgs.RData")
  order = c(2, 1, 6, 5, 4, 3, 8, 7)
  
  for(i in 1:length(order)) {
    thisI = order[i]
    theseArguments = spde_lcpbSimStudyCommandArgs[[thisI]]
    print(do.call("makeSimulationStudyTable", theseArguments))
    browser()
  }
  
  invisible(NULL)
}

makeSimulationStudyTable = function(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), 
                                     effRange=150, beta0=-3.9, representativeSampling=FALSE, 
                                     maxDataSets=10, surveyI = 1:10) {
  
  # get file id
  dataID = paste0("Beta", round(beta0, 4), "rho", round(rho, 4), "sigmaEps", 
                  round(sigmaEpsilon, 4), "gamma", round(gamma, 4))
  # out = load(paste0("savedOutput/simDataSets/simDataMulti", dataID, ".RData"))
  fileName = paste0("savedOutput/simStudyResults/scoresLCPB_", dataID, "repSamp", representativeSampling, ".RData")
  # save(allScoreslcpb, allScoresLcpb, allScoresLCpb, allScoresLCPb, allScoresLCPB, file=fileName)
  
  # load the scoring rules
  out = load(file=fileName)
  
  # make the tables
  # require(kableExtra)
  require(xtable)
  
  # first the prevalence scores
  tabPrevalence = rbind(colMeans(allScoreslcpb$prevalenceScores), 
                        colMeans(allScoresLCpb$prevalenceScores), 
                        colMeans(allScoresLCPB$prevalenceScores))[,-c(2,3)]
  
  tabCounts = rbind(colMeans(allScoreslcpb$countScores), 
                    colMeans(allScoresLCpb$countScores), 
                    colMeans(allScoresLCPB$countScores))[,-c(2,3)]
  
  tabRelativePrevalence = rbind(colMeans(allScoreslcpb$relativePrevalenceScores), 
                                colMeans(allScoresLCpb$relativePrevalenceScores), 
                                colMeans(allScoresLCPB$relativePrevalenceScores))[,-c(2,3)]
  
  factorsPrevalence = c(100, 10^3, 10^3, 100, 100)
  factorsCounts = c(1, 10^(-2), 1, 100, 1)
  factorsRelativePrevalence = c(100, 10, 10, 100, 100)
  digitsPrevalence = c()
  
  rownames(tabPrevalence) = c("S", "SC", "SCP")
  rownames(tabCounts) = c("S", "SC", "SCP")
  rownames(tabRelativePrevalence) = c("S", "SC", "SCP")
  
  # xtable(sweep(tabPrevalence, 2, factorsPrevalence, "*"))
  # xtable(sweep(tabCounts, 2, factorsCounts, "*"))
  # xtable(sweep(tabRelativePrevalence, 2, factorsRelativePrevalence, "*"))
  
  tab = rbind(tabPrevalence, tabCounts, tabRelativePrevalence)
  factors = c(1, 1, 1, 100, 1)
  tab = sweep(tab, 2, factors, "*")
  # xtable(tab, digits=4)
  
  # try using kableExtra
  require(kableExtra)
  tab[,4] = round(tab[,4])
  tab[4:6,] = round(tab[4:6,])
  tab = as.data.frame(tab)
  tab = cbind(c(rep(c("S", "SC", "SCP"), 3)), tab)
  names(tab)[1] = ""
  rownames(tab) = NULL
  labelRoot = paste0("Beta", round(beta0, 4), "gamma", round(gamma, 4), "repSample", representativeSampling)
  thisLabel = paste0("tab:thirdProject:simulationStudy", labelRoot)
  print(tab %>% kable("latex", escape = F, booktabs = TRUE, format.args=list(drop0trailing=TRUE, scientific=FALSE), 
  align=c("l", rep("r", ncol(tab) - 1)), digits=3, caption="", label=thisLabel) %>% pack_rows(
    index = c("Prevalence" = 3, "Total Deaths" = 3, "Relative Prevalence" = 3)
  ))
}









