# function for making fancy colord table of simulation study results. 
# Refer to 
# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_pdf.pdf
# for information about kableExtra
makeFancyTable = function(meanScoresDF, type=c("PvSR", "RvSR", "PvR")) {
  type = match.arg(type)
  
  # meanScoresDF contains the following variables that we actually care about:
  scoreVars = c("RMSE", "Bias", "CRPS", 
                "IntervalScore80", "IntervalScore90", "IntervalScore95", 
                "Coverage80", "Coverage90", "Coverage95", 
                "Width80", "Width90", "Width95", "Time")
  
  # we want to calculate the percent increase of each score value for the 
  # proposed model versus the comparison model. Start by subsetting scores 
  # by relevant model types, then calculate percent increase in scores
  if(type == "PvSR") {
    tab1 = meanScoresDF[meanScoresDF$Model == "Prevalence",]
    tab2 = meanScoresDF[meanScoresDF$Model == "SmoothRisk",]
  } else if(type == "RvSR") {
    tab1 = meanScoresDF[meanScoresDF$Model == "Risk",]
    tab2 = meanScoresDF[meanScoresDF$Model == "SmoothRisk",]
  } else if(type == "PvR") {
    tab1 = meanScoresDF[meanScoresDF$Model == "Prevalence",]
    tab2 = meanScoresDF[meanScoresDF$Model == "Risk",]
  }
  if(type %in% c("PvSR", "RvSR", "PvR")) {
    tab1$Model = NULL
    tab2$Model = NULL
    meanScoresDF = meanScoresDF[meanScoresDF$Model == "Prevalence",]
    meanScoresDF[scoreVars] = (tab1[scoreVars] - tab2[scoreVars])/tab2[scoreVars] * 100
    digits = rep(1, 9)
  } 
  
  require(kableExtra)
  
  # make a fancy table for percent increase of each scoreVar
  for(i in 1:length(scoreVars)) {
    thisScore = scoreVars[i]
    thisTab = meanScoresDF[c("beta", "rho", "nClustFac", "nEAsFac", thisScore)]
    thisTab = thisTab[order(thisTab$beta, decreasing=TRUE),]
    thisTab = thisTab[order(thisTab$nClustFac, decreasing=TRUE),]
    thisTab = thisTab[order(thisTab$rho),]
    thisTab = thisTab[order(thisTab$nEAsFac),]
    thisTab = matrix(unlist(thisTab[thisScore]), nrow=6, ncol=9)
    
    browser()
    
    options(knitr.table.format="latex")
    
    # p. 15: nice color use illustration:
    # vs_dt <- iris[1:10, ]
    # vs_dt[1:4] <- lapply(vs_dt[1:4], function(x) {
    #   cell_spec(x, bold = T,
    #             color = spec_color(x, end = 0.9),
    #             font_size = spec_font_size(x))
    # })
    # vs_dt[5] <- cell_spec(vs_dt[[5]], color = "white", bold = T,
    #                       background = spec_color(1:10, end = 0.9, option = "A", direction = -1))
    # kbl(vs_dt, booktabs = T, escape = F, align = "c", format="latex") %>%
    #   kable_styling(latex_options="striped", full_width = F)
    
    #  p. 17: 3 headers with nice hlines
    # dt <- mtcars[1:5, 1:6]
    # kbl(dt, booktabs = T) %>%
    #   kable_styling(latex_options = "striped") %>%
    #   add_header_above(c(" ", "Group 1" = 2, "Group 2" = 2, "Group 3" = 2)) %>%
    #   add_header_above(c(" ", "Group 4" = 4, "Group 5" = 2)) %>%
    #   add_header_above(c(" ", "Group 6" = 6), bold = T, italic = T)
    
    #  p. 20: collapse_rows with correct hlines
    # collapse_rows_dt <- data.frame(C1 = c(rep("a", 10), rep("b", 5)),
    #                                C2 = c(rep("c", 7), rep("d", 3), rep("c", 2), rep("d", 3)),
    #                                C3 = 1:15,
    #                                C4 = sample(c(0,1), 15, replace = TRUE))
    # kbl(collapse_rows_dt, booktabs = T, align = "c", format="latex") %>%
    #   collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle") %>%
    #   column_spec(1, bold=T)
    
    tempTab = matrix(as.numeric(format(thisTab, digits=rep(1, 9))), nrow=6)
    valRange = range(tempTab)
    tempTab = data.frame(tempTab)
    tempTab <- lapply(tempTab, function(x) {
      cell_spec(x, color = "white", bold = T, format="latex", 
                background = spec_color(x, end = 0.9, option = "A", scale_from=valRange))
    })
    
    # add in buffer columns with small width
    buffCol = list(rep("", 6))
    tempTab = c(list(rep("$r_{mbox{clust}}$", 6), rep(c("1/3", "1", "3"), each=2), 
                     paste(rep("$beta", 6), rep(c("1$", "2$", "3$"), each=2), sep=""), as.character(rep(c(0, -4), 3))), 
                tempTab[1:3], buffCol, tempTab[4:6], buffCol, tempTab[7:9])
    
    rhoVals = c("$1/16$"=1, "$1/4$"=1, "$1/2$"=1)
    kbl(data.frame(tempTab), booktabs = T, escape = F, align = "c", format="latex", 
        linesep=c("", "\\addlinespace"), digits=2, caption="test", col.names=rep("a", 15), 
        toprule=F) %>% 
      column_spec(column=c(8, 12), width="0em") %>%  
      collapse_rows(columns = 2, target=2, headers_to_remove=2, latex_hline = "none", valign = "middle") %>%
      add_header_above(c(" "=4, rep(c(rhoVals, " "=1), 2), rhoVals), escape=F) %>%
      add_header_above(c(" "=4, "$rho$" = 3, " "=1, "$rho$" = 3, " "=1, "$rho$" = 3), escape=F, line=F) %>%
      add_header_above(c(" "=4, "$1/5$" = 3, " "=1, "$1$" = 3, " "=1, "$5$" = 3), escape=F) %>%
      add_header_above(c(" "=4, "$r_{EAs}$" = 11), escape=F, line=F) %>%
      kable_styling(latex_options="basic", full_width = F)
    
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





