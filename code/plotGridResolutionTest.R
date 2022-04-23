# script for generating results for the grid resolution test

makeGridResTable = function() {
  out = load("savedOutput/simpleExample/gridResolutionTestNairobi_1_100.RData")
  out2 = load("savedOutput/simpleExample/gridResTest_datsAndTruths.RData")
  
  # browser()
  resolutions = c(.2, 1, 5, 25)
  nSamples = c(c(500, 1000, 5000), rep(10000, length(resolutions)-3))
  ##### Plot results ----
  
  # truePrevalenceConstituencyKenya
  # truePrevalenceCountyKenya
  
  # for each truth, calculate coverage and other metrics.
  # In aggregate, generate:
  #   boxplot of individual truth coverages versus res.
  #   scatterplot of average (over all truths) coverage vs res.
  #   pair plot (over all truths and resolutions) of central predictions
  #   boxplot (over all truths) of 80% and 95% CI widths
  
  maxSamples = max(nSamples)
  maxSamples = min(nSamples)
  maxSamples = 1000
  allPredsSmoothRisk = list()
  allPredsRisk = list()
  allPredsPrevalence = list()
  allPredsGriddedRisk = list()
  allCIWidthsSmoothRisk = list()
  allCIWidthsRisk = list()
  allCIWidthsPrevalence = list()
  allCIWidthsGriddedRisk = list()
  allCoveragesSmoothRisk = list()
  allCoveragesRisk = list()
  allCoveragesPrevalence = list()
  allCoveragesGriddedRisk = list()
  for(i in 1:length(truths)) {
    # get truth
    truePrevalenceConstituencyKenya = truths[[i]]$truePrevalenceConstituencyKenya
    
    # Calculate RMSE, 80% and 95% Coverage
    predsSmoothRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
    predsRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
    predsPrevalence = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
    predsGriddedRisk = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
    residsConstituencySmoothRisk = list()
    residsConstituencyRisk = list()
    residsConstituencyPrevalence = list()
    residsConstituencyGriddedRisk = list()
    for(j in 1:length(resolutions)) {
      thisMaxSamples = min(c(nSamples[j], maxSamples))
      predsSmoothRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pSmoothRisk[,1:thisMaxSamples])
      predsRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScaleRisk[,1:thisMaxSamples])
      predsPrevalence[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScalePrevalence[,1:thisMaxSamples])
      predsGriddedRisk[,j] = rowMeans(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pGriddedRisk[,1:thisMaxSamples])
      theseResidsSmoothRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pSmoothRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
      theseResidsRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScaleRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
      theseResidsPrevalence = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pFineScalePrevalence[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
      theseResidsGriddedRisk = sweep(allAggResultsN[[i]][[j]]$subareaPop$aggregationResults$pGriddedRisk[,1:thisMaxSamples], 1, truePrevalenceConstituencyKenya, "-")
      residsConstituencySmoothRisk = c(residsConstituencySmoothRisk, list(theseResidsSmoothRisk))
      residsConstituencyRisk = c(residsConstituencyRisk, list(theseResidsRisk))
      residsConstituencyPrevalence = c(residsConstituencyPrevalence, list(theseResidsPrevalence))
      residsConstituencyGriddedRisk = c(residsConstituencyGriddedRisk, list(theseResidsGriddedRisk))
    }
    allPredsSmoothRisk = c(allPredsSmoothRisk, list(predsSmoothRisk))
    allPredsRisk = c(allPredsRisk, list(predsRisk))
    allPredsPrevalence = c(allPredsPrevalence, list(predsPrevalence))
    allPredsGriddedRisk = c(allPredsGriddedRisk, list(predsGriddedRisk))
    
    lowConstituencySmoothRisk = lapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
    lowConstituencyRisk = lapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
    lowConstituencyPrevalence = lapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
    lowConstituencyGriddedRisk = lapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.025, .05, .1))})})
    highConstituencySmoothRisk = lapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
    highConstituencyRisk = lapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
    highConstituencyPrevalence = lapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
    highConstituencyGriddedRisk = lapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=c(.975, .95, .9))})})
    
    # CIWidthSmoothRisk = highConstituencySmoothRisk - lowConstituencySmoothRisk
    # CIWidthRisk = highConstituencyRisk - lowConstituencyRisk
    # CIWidthPrevalence = highConstituencyPrevalence - lowConstituencyPrevalence
    # CIWidthGriddedRisk = highConstituencyGriddedRisk - lowConstituencyGriddedRisk
    CIWidthSmoothRisk = lapply(residsConstituencySmoothRisk, 
                               function(mat) {
                                 apply(mat, 1, function(x) {
                                   quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                               })
    CIWidthRisk = lapply(residsConstituencyRisk, 
                         function(mat) {
                           apply(mat, 1, function(x) {
                             quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                         })
    CIWidthPrevalence = lapply(residsConstituencyPrevalence, 
                               function(mat) {
                                 apply(mat, 1, function(x) {
                                   quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                               })
    CIWidthGriddedRisk = lapply(residsConstituencyGriddedRisk, 
                                function(mat) {
                                  apply(mat, 1, function(x) {
                                    quantile(x, probs=c(.975, .95, .9)) - quantile(x, probs=c(.025, .05, .1))})
                                })
    allCIWidthsSmoothRisk = c(allCIWidthsSmoothRisk, list(CIWidthSmoothRisk))
    allCIWidthsRisk = c(allCIWidthsRisk, list(CIWidthRisk))
    allCIWidthsPrevalence = c(allCIWidthsPrevalence, list(CIWidthPrevalence))
    allCIWidthsGriddedRisk = c(allCIWidthsGriddedRisk, list(CIWidthGriddedRisk))
    
    # this is never used anyway, so commented out:
    # meanCIWidthSmoothRisk = colMeans(CIWidthSmoothRisk)
    # meanCIWidthRisk = colMeans(CIWidthRisk)
    # meanCIWidthPrevalence = colMeans(CIWidthPrevalence)
    # meanCIWidthGriddedRisk = colMeans(CIWidthGriddedRisk)
    
    # inCISmoothRisk = (0 <= highConstituencySmoothRisk) & (0 >= lowConstituencySmoothRisk)
    # inCIRisk = (0 <= highConstituencyRisk) & (0 >= lowConstituencyRisk)
    # inCIPrevalence = (0 <= highConstituencyPrevalence) & (0 >= lowConstituencyPrevalence)
    # inCIGriddedRisk = (0 <= highConstituencyGriddedRisk) & (0 >= lowConstituencyGriddedRisk)
    inCISmoothRisk = lapply(residsConstituencySmoothRisk, 
                            function(mat) {
                              apply(mat, 1, function(x) {
                                (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                            })
    inCIRisk = lapply(residsConstituencyRisk, 
                      function(mat) {
                        apply(mat, 1, function(x) {
                          (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                      })
    inCIPrevalence = lapply(residsConstituencyPrevalence, 
                            function(mat) {
                              apply(mat, 1, function(x) {
                                (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                            })
    inCIGriddedRisk = lapply(residsConstituencyGriddedRisk, 
                             function(mat) {
                               apply(mat, 1, function(x) {
                                 (0 <= quantile(x, probs=c(.975, .95, .9))) & (0 >= quantile(x, probs=c(.025, .05, .1)))})
                             })
    
    # coverageSmoothRisk = colMeans(inCISmoothRisk)
    # coverageRisk = colMeans(inCIRisk)
    # coveragePrevalence = colMeans(inCIPrevalence)
    # coverageGriddedRisk = colMeans(inCIGriddedRisk)
    coverageSmoothRisk = sapply(inCISmoothRisk, rowMeans)
    coverageRisk = sapply(inCIRisk, rowMeans)
    coveragePrevalence = sapply(inCIPrevalence, rowMeans)
    coverageGriddedRisk = sapply(inCIGriddedRisk, rowMeans)
    allCoveragesSmoothRisk = c(allCoveragesSmoothRisk, list(coverageSmoothRisk))
    allCoveragesRisk = c(allCoveragesRisk, list(coverageRisk))
    allCoveragesPrevalence = c(allCoveragesPrevalence, list(coveragePrevalence))
    allCoveragesGriddedRisk = c(allCoveragesGriddedRisk, list(coverageGriddedRisk))
  }
  
  # compile relevant results
  
  ## central predictions
  allPredsSmoothRiskMat = do.call("rbind", allPredsSmoothRisk)
  allPredsRiskMat = do.call("rbind", allPredsRisk)
  allPredsPrevalenceMat = do.call("rbind", allPredsPrevalence)
  allPredsGriddedRiskMat = do.call("rbind", allPredsGriddedRisk)
  
  ## CI widths
  allCIWidthsSmoothRisk80 = lapply(allCIWidthsSmoothRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[3,]})
  })
  allCIWidthsSmoothRisk80 = do.call("rbind", allCIWidthsSmoothRisk80)
  allCIWidthsSmoothRisk90 = lapply(allCIWidthsSmoothRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[2,]})
  })
  allCIWidthsSmoothRisk90 = do.call("rbind", allCIWidthsSmoothRisk90)
  allCIWidthsSmoothRisk95 = lapply(allCIWidthsSmoothRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[1,]})
  })
  allCIWidthsSmoothRisk95 = do.call("rbind", allCIWidthsSmoothRisk95)
  
  allCIWidthsRisk80 = lapply(allCIWidthsRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[3,]})
  })
  allCIWidthsRisk80 = do.call("rbind", allCIWidthsRisk80)
  allCIWidthsRisk90 = lapply(allCIWidthsRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[2,]})
  })
  allCIWidthsRisk90 = do.call("rbind", allCIWidthsRisk90)
  allCIWidthsRisk95 = lapply(allCIWidthsRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[1,]})
  })
  allCIWidthsRisk95 = do.call("rbind", allCIWidthsRisk95)
  
  allCIWidthsPrevalence80 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[3,]})
  })
  allCIWidthsPrevalence80 = do.call("rbind", allCIWidthsPrevalence80)
  allCIWidthsPrevalence90 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[2,]})
  })
  allCIWidthsPrevalence90 = do.call("rbind", allCIWidthsPrevalence90)
  allCIWidthsPrevalence95 = lapply(allCIWidthsPrevalence, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[1,]})
  })
  allCIWidthsPrevalence95 = do.call("rbind", allCIWidthsPrevalence95)
  
  allCIWidthsGriddedRisk80 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[3,]})
  })
  allCIWidthsGriddedRisk80 = do.call("rbind", allCIWidthsGriddedRisk80)
  allCIWidthsGriddedRisk90 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[2,]})
  })
  allCIWidthsGriddedRisk90 = do.call("rbind", allCIWidthsGriddedRisk90)
  allCIWidthsGriddedRisk95 = lapply(allCIWidthsGriddedRisk, function(listOfResolutions) {
    sapply(listOfResolutions, function(x) {x[1,]})
  })
  allCIWidthsGriddedRisk95 = do.call("rbind", allCIWidthsGriddedRisk95)
  
  meanCIWidthsPrevalence95 = colMeans(allCIWidthsPrevalence95)
  meanCIWidthsRisk95 = colMeans(allCIWidthsRisk95)
  meanCIWidthsSmoothRisk95 = colMeans(allCIWidthsSmoothRisk95)
  meanCIWidthsGriddedRisk95 = colMeans(allCIWidthsGriddedRisk95)
  
  ## coverages
  allCoveragesSmoothRisk80 = sapply(allCoveragesSmoothRisk, function(x) {x[3,]})
  allCoveragesSmoothRisk90 = sapply(allCoveragesSmoothRisk, function(x) {x[2,]})
  allCoveragesSmoothRisk95 = sapply(allCoveragesSmoothRisk, function(x) {x[1,]})
  meanCoveragesSmoothRisk80 = rowMeans(allCoveragesSmoothRisk80)
  meanCoveragesSmoothRisk90 = rowMeans(allCoveragesSmoothRisk90)
  meanCoveragesSmoothRisk95 = rowMeans(allCoveragesSmoothRisk95)
  
  allCoveragesRisk80 = sapply(allCoveragesRisk, function(x) {x[3,]})
  allCoveragesRisk90 = sapply(allCoveragesRisk, function(x) {x[2,]})
  allCoveragesRisk95 = sapply(allCoveragesRisk, function(x) {x[1,]})
  meanCoveragesRisk80 = rowMeans(allCoveragesRisk80)
  meanCoveragesRisk90 = rowMeans(allCoveragesRisk90)
  meanCoveragesRisk95 = rowMeans(allCoveragesRisk95)
  
  allCoveragesPrevalence80 = sapply(allCoveragesPrevalence, function(x) {x[3,]})
  allCoveragesPrevalence90 = sapply(allCoveragesPrevalence, function(x) {x[2,]})
  allCoveragesPrevalence95 = sapply(allCoveragesPrevalence, function(x) {x[1,]})
  meanCoveragesPrevalence80 = rowMeans(allCoveragesPrevalence80)
  meanCoveragesPrevalence90 = rowMeans(allCoveragesPrevalence90)
  meanCoveragesPrevalence95 = rowMeans(allCoveragesPrevalence95)
  
  allCoveragesGriddedRisk80 = sapply(allCoveragesGriddedRisk, function(x) {x[3,]})
  allCoveragesGriddedRisk90 = sapply(allCoveragesGriddedRisk, function(x) {x[2,]})
  allCoveragesGriddedRisk95 = sapply(allCoveragesGriddedRisk, function(x) {x[1,]})
  meanCoveragesGriddedRisk80 = rowMeans(allCoveragesGriddedRisk80)
  meanCoveragesGriddedRisk90 = rowMeans(allCoveragesGriddedRisk90)
  meanCoveragesGriddedRisk95 = rowMeans(allCoveragesGriddedRisk95)
  
  # Make table of mean coverages and CI widths ----
  browser()
  tab = rbind(round(1000*meanCIWidthsPrevalence95, 1), 
              round(1000*meanCIWidthsRisk95, 1), 
              round(1000*meanCIWidthsSmoothRisk95, 1), 
              round(1000*meanCIWidthsGriddedRisk95, 1), 
              round(100* meanCoveragesPrevalence95, 0), 
              round(100* meanCoveragesRisk95, 0), 
              round(100* meanCoveragesSmoothRisk95, 0), 
              round(100* meanCoveragesGriddedRisk95, 0))
  row.names(tab) = rep(c("Empirical", "Latent", "Smooth latent", "Gridded"), 2)
  colnames(tab) = c("200m", "1km", "5km", "25km")
  
  xtable(tab, digits=1)
  
  # Plot CI Widths versus resolution and model
  
  ## 95% CIs
  CIWidth = c(c(allCIWidthsPrevalence95), c(allCIWidthsGriddedRisk95))
  tempRes = resolutions[col(allCIWidthsPrevalence95)]
  tempCon = factor(as.character(poppsubSimple$subarea[col(allCIWidthsPrevalence95)]))
  N=length(tempCon)
  CIWidthFrame = data.frame(Constituency=rep(tempCon, 2), 
                            Resolution=rep(tempRes, 2), 
                            Model=factor(c(rep("Empirical", N), rep("Gridded", N)), 
                                         levels=c("Empirical", "Gridded")), 
                            CIWidth=CIWidth)
  
  # 95% CI Width ----
  pdf(paste0("figures/gridResolutionTest/CIWidthVRes95.pdf"), width=5, height=5)
  ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
    geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
    labs(x="Grid resolution (km)", y="95% credible interval width", fill="Model") + 
    theme_classic() + theme(legend.position = "top", axis.title=element_text(size=12), 
                            axis.text=element_text(size=12), plot.margin=unit(c(0,0,0.45,0), "cm"), 
                            legend.title=element_text(size=12), legend.text=element_text(size=12))
  dev.off()
  
  # Mean EAs per area ----
  eas = meanEAsPerCon2(numToPrint=0, easpa = makeDefaultEASPA())
  eaRange = range(eas$meanTotalEAs)
  require(viridis)
  eaCols = magma(64)
  
  thisKenyaLatRange = c(-4.6, 4.8) # -5, 5.5
  thisKenyaLonRange = c(34.2, 41.5) #33.5, 42.0
  
  pdf(paste0(figDirectory, "gridResolutionTest/EAsPerCon.pdf"), width=5, height=5)
  par(mar=c(3, 3.0, 0, 3.1), oma=c(0, 0, 0, 0.7), mgp=c(1.9,.7,0))
  
  eaTicksSubarea = getLogScaleTicks(eaRange, n=10)
  eaTickLabelsSubarea = as.character(eaTicksSubarea)
  plotMapDat(plotVar=eas$meanTotalEAs, mapDat=adm2compressed, new = TRUE, 
             main="", scaleFun=log, scaleFunInverse=exp, 
             cols=eaCols, ticks=eaTicksSubarea, tickLabels=eaTickLabelsSubarea, 
             zlim=log(eaRange), xlim=thisKenyaLonRange, ylim=thisKenyaLatRange, addColorBar = TRUE, 
             legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.88,.91,.2,.9)), legend.width=3, 
             plotArgs=list(cex.main=1, cex.axis=1, cex.lab=1), legend.mar=0, lwd=.1, border=rgb(.4,.4,.4, .7), 
             xlab="Longitude", ylab="Latitude")
  plotMapDat(mapDat=adm1compressed, lwd=.5, border=rgb(.4,.4,.4))
  
  dev.off()
}
