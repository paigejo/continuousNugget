# function for making fancy colord table of simulation study results. 
# Refer to 
# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_pdf.pdf
# for information about kableExtra
# valRanges: matrix with 2 rows and ncols length equal to the number of scoring rules with 
#            first row being the low end of the score range and second being the high end.
makeFancyTable = function(meanScoresDF, type=c("PvSR", "RvSR", "PvR", "P", "R", "SR"), valRanges=NULL) {
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
    
    if(type == "PvSR")
      relMods = c("prevalence", "smooth risk")
    else if(type == "PvR")
      relMods = c("prevalence", "risk")
    else if(type == "RvSR")
      relMods = c("risk", "smooth risk")
    captionRoot1 = "Mean percent increase in "
    captionRoot2 = paste0(" of the ", relMods[1], " aggregation model relative to the ", 
                          relMods[2], " aggregation model.")
  } else if(type %in% c("P", "R", "S")) {
    if(type == "P") {
      thisMod = "prevalence"
      meanScoresDF = meanScoresDF[meanScoresDF$Model == "Prevalence",]
    } else if(type == "R") {
      thisMod = "risk"
      meanScoresDF = meanScoresDF[meanScoresDF$Model == "Risk",]
    } else if(type == "SR") {
      thisMod = "smooth risk"
      meanScoresDF = meanScoresDF[meanScoresDF$Model == "SmoothRisk",]
    }
    
    captionRoot1 = "Mean "
    captionRoot2 = paste0(" of the ", thisMod, " aggregation model.")
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
    # kbl(dt, booktabs = T, format="latex") %>%
    #   kable_styling(latex_options = "basic") %>%
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
    
    if(type %in% c("PvSR", "RvSR", "PvR")) {
      # we only care about X.X% increases, no more than one decimal place
      tempTab = matrix(as.numeric(formatC(thisTab, digits=1, format="f")), nrow=6)
      formattedTempTab = tempTab
    } else if(type %in% c("P", "R", "S")) {
      # correct formatting is trickier here. Must make sure same for all models, 
      # but can be different for difference scores
      if(thisScore == "RMSE") {
        tempTab = signif(thisTab, 1)
        formattedTempTab = tempTab
      } else if(thisScore == "Bias") {
        tempTab = matrix(round(thisTab, digits=5), nrow=6)
        formattedTempTab = format(tempTab, format="f", digits=1)
        formattedTempTab = num(tempTab, digits=1, fixed_exponent=Inf, notation="eng")
      } else if(thisScore == "CRPS") {
        
      } else if(thisScore == "IntervalScore80") {
        
      } else if(thisScore == "IntervalScore90") {
        
      } else if(thisScore == "IntervalScore95") {
        
      } else if(thisScore == "Coverage80") {
        
      } else if(thisScore == "Coverage90") {
        
      } else if(thisScore == "Coverage95") {
        
      } else if(thisScore == "Width80") {
        
      } else if(thisScore == "Width90") {
        
      } else if(thisScore == "Width95") {
        
      } else if(thisScore == "Time") {
        
      }
    }
    
    if(!is.null(valRanges)) {
      valRange = valRanges[,i]
    } else {
      valRange = range(tempTab)
    }
    # tempTab = matrix(format(tempTab, digits=1, scientific=FALSE), nrow=6)
    tempTab = data.frame(tempTab)
    if(thisScore %in% c("RMSE", "CRPS", 
                        "IntervalScore80", "IntervalScore90", "IntervalScore95", 
                        "Coverage80", "Coverage90", "Coverage95", 
                        "Width80", "Width90", "Width95", "Time")) {
      # lower is better
      customCols = centerColorScale(256, valRange=c(0,1), center=.46, colScale=makeBlueGreenYellowSequentialColors, rev=TRUE)
      tempTab <- lapply(tempTab, function(x) {
        cell_spec(x, color = "white", bold = T, format="latex", 
                  background = my_spec_color(x, end = 0.9, option = "D", scale_from=valRange, customScale=customCols))
      })
      
    } else if(thisScore %in% c("")) {
      # higher is better (there are no scores like this)
    } else if(thisScore %in% c("Bias", "Coverage80", "Coverage90", "Coverage95")) {
      # closer a particular value is better
      
      if(type %in% c("PvSR", "RvSR", "PvR")) {
        # skip these tables, they're useless: what does pct increase cvg even mean??
        next
      }
      
      # we are making a table of the absolute score, coloring by distance from best value
      if(thisScore == "Bias") {
        closeVal = 0
      } else if(thisScores == "Coverage80") {
        closeVal = 80
      } else if(thisScores == "Coverage90") {
        closeVal = 90
      } else if(thisScores == "Coverage95") {
        closeVal = 95
      }
      valRange=valRange-closeVal
      customCols = centerColorScale(256, colScale=makeBlueGoldDivergingColors, 
                                    valRange=valRange, center=0)
      
      tempTab <- lapply(tempTab, function(x) {
        cell_spec(x, color = "black", bold = T, format="latex", 
                  background = my_spec_color(x-closeVal, scale_from=valRange, customScale=customCols))
      })
    }
    # tempTab <- lapply(tempTab, function(x) {
    #   cell_spec(x, color = "white", bold = T, format="latex", 
    #             background = spec_color(x, end = 0.9, option = "A", scale_from=valRange))
    # })
    
    # add in buffer columns with small width
    buffCol = list(rep("", 6))
    tempTab = c(list(c("\\multirow{6}{*}{$r_{\\tiny \\mbox{clust}}$}", rep(" ", 5)), c("\\multirow{2}{*}{$3$}", "", "\\multirow{2}{*}{$1$}", "", "\\multirow{2}{*}{$1/3$}", ""), 
                     rep(c("\\multirow{2}{*}{$\\beta$}", " "), 3), as.character(rep(c(0, -4), 3))), 
                tempTab[1:3], buffCol, tempTab[4:6], buffCol, tempTab[7:9])
    
    rhoVals = c("$1/16$"=1, "$1/4$"=1, "$1/2$"=1)
    kbl(data.frame(tempTab), booktabs = T, escape = F, align = "c", format="latex", 
        linesep=c("", "\\addlinespace"), digits=2, caption=paste0(captionRoot1, thisScore, captionRoot2), 
        col.names=NULL, bottomrule=FALSE, label=paste0(type, "_", thisScore)) %>% 
      column_spec(column=c(8, 12), width="0em") %>%
      column_spec(column=2, border_right=TRUE) %>%
      add_header_above(c(" "=4, rep(c(rhoVals, " "=1), 2), rhoVals), escape=F) %>%
      add_header_above(c(" "=4, "$\\\\rho$" = 3, " "=1, "$\\\\rho$" = 3, " "=1, "$\\\\rho$" = 3), escape=F, line=F) %>%
      add_header_above(c(" "=4, "$1/5$" = 3, " "=1, "$1$" = 3, " "=1, "$5$" = 3), escape=F) %>%
      add_header_above(c(" "=4, "$r_{\\\\tiny \\\\mbox{EAs}}$" = 11), escape=F, line=F) %>% 
      kable_styling(latex_options="basic", full_width = F) # %>%
    # collapse_rows(columns = 2, latex_hline = "linespace", valign = "middle") # this doesn't work. It's a problem with kableExtra
    
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

# same as kableExtra::spec_color, expect 'customScale' argument can be a custom 
# vector of hex colors of length 256
my_spec_color = function(x, alpha = 1, begin = 0, end = 1, direction = 1, option = "D", 
          na_color = "#BBBBBB", scale_from = NULL, customScale=NULL) 
{
  require(scales)
  if (is.null(scale_from)) {
    x <- round(rescale(x, c(1, 256)))
  }
  else {
    x <- round(rescale(x, to = c(1, 256), from = scale_from))
  }
  if(is.null(customScale)) {
    color_code <- viridisLite::viridis(256, alpha, begin, end, 
                                       direction, option)[x]
  } else {
    color_code <- customScale[x]
  }
  
  color_code[is.na(color_code)] <- na_color
  return(color_code)
}

# display the numbers in engineering notation (x 10^{scale}). For 
# numbers larger than 10^{scale}, scale will continued to be used 
# as the exponent.
# scale: if -Inf, uses the smallest scale, if Inf, uses the largest, 
#        otherwise it is the user specified power of the scale to 
#        display the numbers with
formatEngineering = function(df, scale=-Inf, digits=1) {
  if(scale == -Inf) {
    scale = floor(min(sapply(df, function(x) {log10(abs(x))})))
  } else if(scale == Inf) {
    scale = floor(max(sapply(df, function(x) {log10(abs(x))})))
  }
  
  fac = 10^scale
  
  displayFun = function(x) {
    firstNum = round(x / fac, digits=digits)
    paste0(firstNum, " \\times 10^{", scale, "}")
  }
  
  data.frame(lapply(df, displayFun))
}





