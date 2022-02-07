# function for making fancy colord table of simulation study results. 
# Refer to 
# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_pdf.pdf
# for information about kableExtra
makeFancyTable = function(meanScoresDF, type=c("PvSR", "RvSR", "PvR")) {
  type = match.arg(type)
  
  # meanScoresDF contains the following variables that we actually care about:
  scoreVars = c("RMSE", "Bias", 
                "IntervalScore80", "IntervalScore90", "IntervalScore95", 
                "Coverage80", "Coverage90", "Coverage95", 
                "Width80", "Width90", "Width95", "Time")
  
  # we want to calculate the percent increase of each score value for the 
  # proposed model versus the comparison model. Start by subsetting scores 
  # by relevant model types, then calculate percent increase in scores
  if(type == "PvSR") {
    tab1 = meanScoresDF[meanScoresDF$model == "Prevalence",]
    tab2 = meanScoresDF[meanScoresDF$model == "SmoothRisk",]
  } else if(type == "RvSR") {
    tab1 = meanScoresDF[meanScoresDF$model == "Risk",]
    tab2 = meanScoresDF[meanScoresDF$model == "SmoothRisk",]
  } else if(type == "PvR") {
    tab1 = meanScoresDF[meanScoresDF$model == "Prevalence",]
    tab2 = meanScoresDF[meanScoresDF$model == "Risk",]
  }
  if(type %in% c("PvSR", "RvSR", "PvR")) {
    tab1$model = NULL
    tab2$model = NULL
    meanScoresDF[scoreVars] = (tab1[scoreVars] - tab2[scoreVars])/tab2[scoreVars] * 100
  } 
  
  require(kableExtra)
  
  # make a fancy table for percent increase of each scoreVar
  for(i in 1:length(scoreVars)) {
    thisScore = scoreVars[i]
    thisTab = meanScoresTab[c("beta", "rho", "nClustFac", "nEAsFac", thisScore)]
    
    browser()
    
    # p. 15: nice color use illustration:
    # vs_dt <- iris[1:10, ]
    # vs_dt[1:4] <- lapply(vs_dt[1:4], function(x) {
    #   cell_spec(x, bold = T,
    #             color = spec_color(x, end = 0.9),
    #             font_size = spec_font_size(x))
    # })
    # vs_dt[5] <- cell_spec(vs_dt[[5]], color = "white", bold = T,
    #                       background = spec_color(1:10, end = 0.9, option = "A", direction = -1))
    # kbl(vs_dt, booktabs = T, escape = F, align = "c") %>%
    #   kable_classic("striped", full_width = F)
    
    #  p. 17: 3 headers with nice hlines
    # kbl(dt, booktabs = T) %>%
    #   kable_styling(latex_options = "striped") %>%
    #   add_header_above(c(" ", "Group 1" = 2, "Group 2" = 2, "Group 3" = 2)) %>%
    #   add_header_above(c(" ", "Group 4" = 4, "Group 5" = 2)) %>%
    #   add_header_above(c(" ", "Group 6" = 6), bold = T, italic = T
    
    #  p. 20: collapse_rows with correct hlines
    # collapse_rows_dt <- data.frame(C1 = c(rep("a", 10), rep("b", 5)),
    #                                C2 = c(rep("c", 7), rep("d", 3), rep("c", 2), rep("d", 3)),
    #                                C3 = 1:15,
    #                                C4 = sample(c(0,1), 15, replace = TRUE))
    # kbl(collapse_rows_dt, booktabs = T, align = "c") %>%
    #   column_spec(1, bold=T) %>%
    #   collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")
    
    
    
  }
}

# function for making input to makeFancyTable. 
# meanScoresDF contains the following variables:
# beta, rho, nClustFac, nEAsFac, model, 
# rmse, bias, 
# intScore80, intScore90, intScore95, 
# cvg80, cvg90, cvg95, 
# width80, width90, with95, time
getFullMeanScoresDF = function(iRange=1:54, maxJ=100, coarse=TRUE, areaLevel=c("subarea", "area"), 
                               response=c("p", "Z")) {
  areaLevel = match.arg(areaLevel)
  response = match.arg(response)
  coarseText = ifelse(coarse, "Coarse", "")
  
  allDat = c()
  for(i in iRange) {
    # get population arguments
    out = load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
    theseArgs = spde_prevRiskSimStudyCommandArgs[[i]]
    beta0 = theseArgs$beta0
    rho = theseArgs$rho
    nEAsFac = theseArgs$nEAsFac
    nClustFac = theseArgs$nClustFac
    popPar = list(beta=beta0, rho=rho, nEAsFac=nEAsFac, nClustFac=nClustFac)
    
    thisFile = paste0("savedOutput/simStudyResults/simScoresAll_i", i, "maxJ", maxJ, coarseText, ".RData")
    out = load(thisFile)
    
    if(areaLevel == "subarea" && response == "p") {
      allDat = rbind(allDat, data.frame(c(popPar, list(Model="Prevalence"), as.list(subareaScoresPprevAvg))))
      allDat = rbind(allDat, c(popPar, Model="Risk", subareaScoresPriskAvg))
      allDat = rbind(allDat, c(popPar, Model="SmoothRisk", subareaScoresPsmoothRiskAvg))
    } else if(areaLevel == "subarea" && response == "Z") {
      allDat = rbind(allDat, c(popPar, Model="Prevalence", subareaScoresZprevAvg))
      allDat = rbind(allDat, c(popPar, Model="Risk", subareaScoresZriskAvg))
      allDat = rbind(allDat, c(popPar, Model="SmoothRisk", subareaScoresZsmoothRiskAvg))
    } else if(areaLevel == "area" && response == "p") {
      allDat = rbind(allDat, c(popPar, Model="Prevalence", areaScoresPprevAvg))
      allDat = rbind(allDat, c(popPar, Model="Risk", areaScoresPriskAvg))
      allDat = rbind(allDat, c(popPar, Model="SmoothRisk", areaScoresPsmoothRiskAvg))
    } else if(areaLevel == "area" && response == "Z") {
      allDat = rbind(allDat, c(popPar, Model="Prevalence", areaScoresZprevAvg))
      allDat = rbind(allDat, c(popPar, Model="Risk", areaScoresZriskAvg))
      allDat = rbind(allDat, c(popPar, Model="SmoothRisk", areaScoresZsmoothRiskAvg))
    }
  }
  
  allDat
}





