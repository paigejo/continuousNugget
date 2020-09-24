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









