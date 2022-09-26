library(colorspace)
library(mgcv)

# makes plots as well as parameter estimate tables for the example applications in the manuscript
makeAllPlots = function(dataType=c("ed", "mort"), resultFilenames, modelClasses, modelVariations, 
                        areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                        meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                        relWidthRange, relWidthTicks, relWidthTickLabels, 
                        plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                        widthCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                        ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                        plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                        makeModelPredictions=TRUE, makeCovariograms=TRUE, makePairPlots=TRUE, 
                        loadResults=FALSE, saveResults=!loadResults, singleColorBar=TRUE, doModelClassPlots=FALSE, 
                        col=NULL, lty=NULL, pch=NULL) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  
  # print lattice resolutions
  print("SPDE mesh resolution: ")
  test <- getSPDEMeshKenya()
  # plot(test)
  print(meanSegmentLength(getSPDEMeshKenya(), 50))
  print("three layer coarse resolution:")
  test = makeLatGridsKenya(3,NC=14)
  print(test[[1]]$latWidth)
  print("three layer fine resolution:")
  print(test[[3]]$latWidth)
  print("two layer coarse resolution:")
  test = makeLatGridsKenya(2,NC=c(30, 214))
  print(test[[1]]$latWidth)
  print("two layer fine resolution:")
  print(test[[2]]$latWidth)
  
  if(makeModelPredictions) {
    plotModelPredictions(dat, resultFilenames, modelClasses, modelVariations, areaLevels, 
                         meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                         relWidthRange, relWidthTicks, relWidthTickLabels, 
                         varName, plotNameRoot, resultNameRoot, meanCols, widthCols, 
                         kenyaLatRange, kenyaLonRange, singleColorBar, doModelClassPlots)
  }
  
  if(makeCovariograms) {
    plotCovariograms(dat, resultFilenames, modelClasses, modelVariations, 
                     varName, plotNameRoot, resultNameRoot, 
                     cgramList=NULL, loadResults=loadResults, saveResults=saveResults, doModelClassPlots=doModelClassPlots, 
                     col=col, lty=lty, pch=pch)
  }
  
  if(makePairPlots) {
    makeFinalPercentResidualPlot(dat, resultFilenames, modelClasses, modelVariations, 
                                 areaLevels, meanRange, meanTicks, meanTickLabels, 
                                 plotNameRoot=plotNameRoot, resultNameRoot, 
                                 meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange, doModelClassPlots)
    
    makeFinalPairPlot(dat, resultFilenames, modelClasses, modelVariations, 
                      areaLevels, meanRange, meanTicks, meanTickLabels, 
                      plotNameRoot=plotNameRoot, resultNameRoot, 
                      meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange, doModelClassPlots)
    makePairPlots(dat, resultFilenames, modelClasses, modelVariations, 
                  areaLevels, meanRange, meanTicks, meanTickLabels, 
                  plotNameRoot=plotNameRoot, resultNameRoot, 
                  meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange, doModelClassPlots)
  }
  
  invisible(NULL)
}

plotDataVisualizations = function(dataType=c("mort", "ed"), dat=NULL, 
                                  varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                  sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                                  ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                                  plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    # out = load("../U5MR/kenyaData.RData")
    dat = mort
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
  }
  
  kenyaMap = adm0
  countyMap = adm1compressed
  conMap = adm2compressed
  
  print("generating data visualizations...")
  
  # plot a map of urbanicity if requested (this can take ~10 minutes)
  if(plotUrbanMap) {
    # inside this if statement since it takes around ten minutes to run
    makeUrbanMap(kmres=5, savePlot=TRUE, lonLim=kenyaLonRange, latLim=kenyaLatRange, main="")
    makeDensityMap(kmres=1)
  }
  browser()
  
  thisKenyaLatRange = c(-4.6, 4.8) # -5, 5.5
  thisKenyaLonRange = c(34.2, 41.5) #33.5, 42.0
  
  pdf(file=paste0("figures/application/mort.pdf"), width=5, height=5)
  # par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 3.1, 6))
  par(mar=c(3, 3.0, 0, 3.1), oma=c(0, 0, 0, 0.9), mgp=c(1.9,.7,0))
  cols = makeBlueGreenYellowSequentialColors(64, rev=TRUE)
  prevs = mort$y/mort$n
  # varRange = range(prevs, na.rm=TRUE)
  # par( oma=c( 0,0,0,5)) # save some room for the legend
  plot(mort$lon, mort$lat, type="n", xlim=thisKenyaLonRange, ylim=thisKenyaLatRange, 
       xlab="Longitude", ylab="Latitude", asp=1)
  
  plotMapDat(mapDat=adm2compressed, lwd=.2, border=rgb(.5, .5, .5, .5))
  plotWithColor(mort$lon, mort$lat, prevs, colScale=cols, #scaleFun=shrunkLogit, 
                #scaleFunInverse=shrunkExpit, forceColorsInRange=TRUE, 
                pch=21, cex=.4, ordering="none", new=FALSE, 
                legendArgs=list(axis.args=list(cex.axis=1, tck=-.7, hadj=.1), legend.cex=2, smallplot=c(.9,.93,.2,.9)), 
                xlab="Longitude", ylab="Latitude", asp=1, colorName="bg",
                col=rgb(0, 0, 0, .2))
  plotMapDat(mapDat=adm1compressed, lwd=.5, col=rgb(1, 1, 1, .5))
  
  # world(add=TRUE)
  # par( oma=c(0,0,0,2))
  dev.off()
}

makeUrbanMap = function(popGrid=NULL, kmres=1, savePlot=FALSE, fileName=ifelse(whiteRural, "Figures/UrbanMapWhiteRural.png", "Figures/urbanMap.png"), 
                        nx=850, ny=1050, width=800, height=800, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                        lonLim=kenyaLonRange, latLim=kenyaLatRange, whiteRural=TRUE, main="") {
  # get prediction locations from population grid
  if(is.null(popGrid)) {
    if(kmres == 5)
      load("savedOutput/global/popGrid.RData")
    else {
      load("savedOutput/global/pop.rda")
      popGrid = makePopIntegrationTab(kmRes=kmres, pop=pop, domainPoly=kenyaPoly, 
                                           eastLim=eastLim, northLim=northLim, 
                                           mapProjection=projKenya, poppa=poppaKenya, 
                                           poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                           areaMapDat=adm1, subareaMapDat=adm2)
    }
      
  }
  
  out = load("savedOutput/global/adminMapData.RData")
  countyMap = adm1compressed
  
  # determine which points in Kenya are urban
  urban = popGrid$urban
  
  if(savePlot) {
    png(file=fileName, width=width, height=height)
    par(oma=c( 0,0,0,0), mar=c(5.5, 6.1, 3.5, 6))
  }
  plot(popGrid$lon, popGrid$lat, xlab="Longitude", ylab="Latitude", main=main, xlim=lonLim, ylim=latLim, asp=1, type="n", 
       cex.axis=2, cex.lab=2, tck=-.03, mgp=c(4, 2, 0))
  # quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, 
  #            xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1)
  if(whiteRural)
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c(rgb(0, 0, 0, 0), "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE, FUN=max)
  else {
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
    legend("bottomleft", c("urban", "rural"), col=c("blue", "green"), pch=19)
  }
  # world(add=TRUE)
  plotMapDat(mapDat=adm1compressed, lwd=.5)
  if(savePlot) {
    dev.off()
  }
  browser()
}

makeDensityMap = function(popGrid=NULL, kmres=1, savePlot=TRUE, fileName="figures/application/densityMap.pdf", 
                          nx=850, ny=1050, width=800, height=800, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                          lonLim=kenyaLonRange, latLim=kenyaLatRange, main="", delta=NULL, mean.neighbor=NULL) {
  # get prediction locations from population grid
  if(is.null(popGrid)) {
    if(kmres == 5)
      load("savedOutput/global/popGrid.RData")
    else {
      load("savedOutput/global/pop.rda")
      if(kmres == 1 && is.null(delta) && is.null(mean.neighbor)) {
        delta=0.2
        mean.neighbor=30000
        # delta=.9
        # mean.neighbor=100000
      }
      popGrid = makePopIntegrationTab(kmRes=kmres, pop=pop, domainPoly=kenyaPoly, 
                                      eastLim=eastLim, northLim=northLim, 
                                      mapProjection=SUMMER:::projKenya, poppa=poppaKenya, 
                                      poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                      areaMapDat=adm1, subareaMapDat=adm2, delta=delta, mean.neighbor=mean.neighbor)
    }
    
  }
  
  out = load("savedOutput/global/adminMapData.RData")
  countyMap = adm1compressed
  conMap = adm2compressed
  
  # determine which points in Kenya are urban
  urban = popGrid$urban
  thisKenyaLatRange = c(-4.6, 4.8) # -5, 5.5
  thisKenyaLonRange = c(34.2, 41.5) #33.5, 42.0
  popCols = makeBlueSequentialColors(64)
  
  if(savePlot) {
    # png(file=fileName, width=width, height=height)
    pdf(file=fileName, width=5, height=5)
    # par(oma=c( 0,0,0,3), mar=c(5.5, 6.1, 3.5, 6))
    par(mar=c(3, 3.0, 0, 3.1), oma=c(0, 0, 0, 0.9), mgp=c(1.9,.7,0))
  }
  plot(popGrid$lon, popGrid$lat, xlab="Longitude", ylab="Latitude", main=main, 
       xlim=thisKenyaLonRange, ylim=thisKenyaLatRange, asp=1, type="n")
  # quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, 
  #            xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1)
  
  tempPopGrid = popGrid[popGrid$pop > 10,]
  popTicks = getLogScaleTicks(tempPopGrid$pop, nint=5)
  popTickLabels = as.character(popTicks)
  quilt.plot(tempPopGrid$lon, tempPopGrid$lat, log(tempPopGrid$pop), nx=450, ny=550, add.legend = FALSE, add=TRUE, col=popCols)
  image.plot(zlim=range(log(tempPopGrid$pop)), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=.1), 
             legend.mar = 0, legend.cex=2, legend.width=3,  smallplot=c(.898,.928,.2,.9))
  
  # image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
  #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
  #            legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
  
  plotMapDat(mapDat=adm1compressed, lwd=.5)
  plotMapDat(mapDat=adm2compressed, lwd=.2, border=rgb(.5, .5, .5, .5))
  if(savePlot) {
    dev.off()
  }
  browser()
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
plotModelPredictions = function(dat, resultFilenames, modelClasses, modelVariations, 
                                areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                                meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                                relWidthRange, relWidthTicks, relWidthTickLabels, 
                                varName="education", plotNameRoot="Education", resultNameRoot="Ed", 
                                meanCols=makeRedBlueDivergingColors(64), 
                                widthCols=makeBlueYellowSequentialColors(64), 
                                kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                                singleColorBar=TRUE, doModelClassPlots=FALSE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1 && doModelClassPlots) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
      if(sum(thisI) == 1)
        next
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      plotModelPredictions(dat, thisResultFilenames, thisModelClasses, thisModelVariations, areaLevels, 
                           meanRange, meanTicks, meanTickLabels, widthRange, widthTicks, widthTickLabels, 
                           relWidthRange, relWidthTicks, relWidthTickLabels, 
                           varName, plotNameRoot, resultNameRoot, meanCols, widthCols, 
                           kenyaLatRange, kenyaLonRange, singleColorBar)
    }
  }
  
  # load the fine grid used to approximate continuous prediction
  print("Loading prediction grid and shapefiles...")
  load("../U5MR/popGrid.RData")
  
  # load shape files for plotting
  require(maptools)
  regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1compressed
  
  ##### central estimates and credible interval widths
  
  if(length(uniqueModelClasses) == 1)
    extraPlotNameRoot = uniqueModelClasses[1]
  else
    extraPlotNameRoot = ""
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    # get map for this areal aggregation level
    if(thisArea == "Cluster") {
      next
    } else if(thisArea == "Region") {
      thisMap = regionMap
    } else {
      thisMap = countyMap
    }
    
    print("Plotting central estimates and credible interval widths...")
    width = 400 * numberModels
    
    png(file=paste0("Figures/", resultNameRoot, "/preds", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      
      # load this model and plot results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      predictionList = c(predictionList, list(theseResults$preds))
      widthList = c(widthList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    # plot the predictions
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea %in% c("Region", "County")) {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main="", scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "Pixel"){
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
        plotMapDat(mapDat=thisMap, lwd=.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                     col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      mtext(side = 3, as.expression(modelName), line = 1, cex=2)
      if(j == 1)
        mtext(side = 2, "Estimates", line = 8, cex=2)
    }
    
    # plot the credible interval widths
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea %in% c("Region", "County")) {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main="", scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(widthRange), ticks=widthTicks, tickLabels=widthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "Pixel") {
        
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=range(log(widthRange)))
        plotMapDat(mapDat=thisMap, lwd=.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(log(widthRange)), nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                     col=widthCols, add = TRUE, axis.args=list(at=log(widthTicks), labels=widthTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      if(j == 1)
        mtext(side = 2, "80% CI Width", line = 8, cex=2)
    }
    dev.off()
    
    ##### Do the same thing but with relative credible interval widths
    png(file=paste0("Figures/", resultNameRoot, "/predsRelativeWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,6,0,7), mar=c(5.1, 5.1, 4.1, 2.5))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      
      # load this model and plot results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      predictionList = c(predictionList, list(theseResults$preds))
      widthList = c(widthList, list((theseResults$Q90 - theseResults$Q10)/theseResults$preds))
    }
    
    # plot the predictions
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea %in% c("Region", "County")) {
        plotMapDat(plotVar=predictionList[[j]], new = TRUE, 
                   main="", scaleFun=logit, scaleFunInverse=expit, 
                   cols=meanCols, zlim=logit(meanRange), ticks=meanTicks, tickLabels=meanTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "Pixel"){
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
        plotMapDat(mapDat=thisMap, lwd=.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                     col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      mtext(side = 3, as.expression(modelName), line = 1, cex=2)
      if(j == 1)
        mtext(side = 2, "Estimates", line = 8, cex=2)
    }
    
    # plot the relative credible interval widths
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      thisAddColorBar = (j == numberModels) | !singleColorBar
      
      # load this model and plot results in the given column
      # out = load(resultFilenames[j])
      # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      if(thisArea %in% c("Region", "County")) {
        plotMapDat(plotVar=widthList[[j]], new = TRUE, 
                   main="", scaleFun=log, scaleFunInverse=exp, 
                   cols=widthCols, zlim=log(relWidthRange), ticks=relWidthTicks, tickLabels=relWidthTickLabels, 
                   xlim=kenyaLonRange, ylim=kenyaLatRange, addColorBar = thisAddColorBar, 
                   legendArgs=list(axis.args=list(cex.axis=2, tck=-.7, hadj=-.1), legend.cex=2, smallplot= c(.97,1,.1,.9)), legend.width=3, 
                   plotArgs=list(cex.main=3, cex.axis=2, cex.lab=2), legend.mar=0)
      } else if(thisArea == "Pixel") {
        
        plot(cbind(popGrid$lon, popGrid$lat), type="n", main="", ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1, cex.main=3, cex.axis=2, cex.lab=2)
        quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[j]]), 
                   nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=range(log(relWidthRange)))
        plotMapDat(mapDat=thisMap, lwd=.5)
        points(dat$lon, dat$lat, pch=".")
        
        if(thisAddColorBar) {
          image.plot(zlim=range(log(relWidthRange)), nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                     col=widthCols, add = TRUE, axis.args=list(at=log(relWidthTicks), labels=relWidthTickLabels, cex.axis=2, tck=-.7, hadj=-.1), 
                     legend.mar = 0, legend.cex=2, legend.width=3, smallplot= c(.97,1,.1,.9))
          
          # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
          #            legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
          #            legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
        }
      }
      
      if(j == 1)
        mtext(side = 2, "Relative 80% CI Width", line = 8, cex=2)
    }
    dev.off()
  }
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
makePairPlots = function(dat, resultFilenames, modelClasses, modelVariations, 
                         areaLevels=c("Region", "County", "Pixel"), 
                         meanRange, meanTicks, meanTickLabels, 
                         plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                         ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                         kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), doModelClassPlots=FALSE, 
                         areaMat=NULL, radiiMat=NULL) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1 && doModelClassPlots) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      
      if(sum(thisI) == 1)
        next
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      makePairPlots(dat, thisResultFilenames, thisModelClasses, thisModelVariations, 
                    areaLevels, meanRange, meanTicks, meanTickLabels, 
                    plotNameRoot=paste0(plotNameRoot, thisModelClass), 
                    resultNameRoot, meanCols, ncols, urbCols, kenyaLatRange, kenyaLonRange)
    }
  }
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    if(length(uniqueModelClasses) == 1)
      extraPlotNameRoot = uniqueModelClasses[1]
    else
      extraPlotNameRoot = ""
    
    if(thisArea %in% c("Region", "County")) {
      
    }
    else {
      
    }
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    widthsList = list()
    # modelNames = c()
    modelNames = list()
    for(j in 1:numberModels) {
      # modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      modelNames = c(modelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
      widthsList = c(widthsList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    valMat = do.call("cbind", predsList)
    zlim = range(valMat)
    
    if(thisArea %in% c("Region", "County")) {
      
      # valMat = rbind(1:5, valMat)
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        # else
        #   ylim = zlim2
        # if(diff(range(c(x, y))) > 0.04)
        #   par(usr = c(zlim, zlim))
        # else
        #   par(usr = c(zlim2, zlim2))
        # par(usr = c(xlim, ylim))
        # points(x,y,..., col="blue")
        abline(a = 0,b = 1,...)
        points(x,y,..., col=theseCols)
      }
      
      width = 400 * numberModels
      png(file=paste0("Figures/", resultNameRoot, "/pairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=width)
      
      # pairs(valMat, 
      #       modelNames, 
      #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
      #       main=paste0(thisArea, " estimate comparisons"))
      # lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=2, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " estimate comparisons"), 
              lims=lims, oma=c(3,3,6,10), cex.main=2)
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
                 legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      width = 400 * numberModels
      png(file=paste0("Figures/", resultNameRoot, "/pairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=width)
      
      # pairs(valMat, 
      #       modelNames, 
      #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
      #       main=paste0(thisArea, " estimate comparisons"))
      # lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=2, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " posterior 80% CI width comparisons"), 
              lims=lims, oma=c(3,3,6,10), cex.main=2)
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=7, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=3.0, legend.width=1, legend.shrink=.9, 
                 legend.cex=2, axis.args=list(cex.axis=1, tck=-1, hadj=-.1))
      dev.off()
    } else {
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        # else
        #   ylim = zlim2
        # if(diff(range(c(x, y))) > 0.04)
        #   par(usr = c(zlim, zlim))
        # else
        #   par(usr = c(zlim2, zlim2))
        # par(usr = c(xlim, ylim))
        # points(x,y,..., col="blue")
        abline(a = 0,b = 1,...)
        points(x,y,..., col=theseCols)
      }
      
      width = 2 * numberModels
      pdf(file=paste0("Figures/", resultNameRoot, "/pairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=width)
      
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " estimate comparisons"), 
              lims=lims, oma=c(3,3,6,7))
      # legend("topleft", c("Urban", "Rural"), col=c(urbCols[ncols], urbCols[1]), pch=19)
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      pdf(file=paste0("Figures/", resultNameRoot, "/pairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=width)
      
      lims = rep(list(zlim), numberModels)
      myPairs(valMat, 
              as.expression(unlist(modelNames)), 
              pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
              main=paste0(thisArea, " posterior 80% CI width comparisons"), 
              lims=lims, oma=c(3,3,6,7))
      legend("topleft", c("Urban", "Rural"), col=c(urbCols[ncols], urbCols[1]), pch=19)
      dev.off()
    }
  }
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
makeFinalPairPlot = function(dat, resultFilenames, modelClasses, modelVariations, 
                         areaLevels=c("Region", "County", "Pixel"), 
                         meanRange, meanTicks, meanTickLabels, 
                         plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                         ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                         kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), doModelClassPlots=FALSE, 
                         areaMat=NULL, radiiMat=NULL, areaColors=makePurpleYellowSequentialColors(ncols)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    if(length(uniqueModelClasses) == 1)
      extraPlotNameRoot = uniqueModelClasses[1]
    else
      extraPlotNameRoot = ""
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    widthsList = list()
    # modelNames = c()
    modelNames = list()
    for(j in 1:numberModels) {
      # modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      modelNames = c(modelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
      widthsList = c(widthsList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    valMat = do.call("cbind", predsList)
    zlim = range(valMat)
    
    if(thisArea %in% c("Region", "County")) {
      # determine colorings by area and radius
      theseAreas = getArea(thisArea)
      theseRadii = getRadius(thisArea)
      colI = cut(theseAreas, breaks=seq(min(theseAreas) - .0001, max(theseAreas) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseAreaColors = areaColors[colI]
      colI = cut(theseRadii, breaks=seq(min(theseRadii) - .0001, max(theseRadii) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseRadiusColors = areaColors[colI]
      theseUrbanColors = theseCols
      
      # valMat = rbind(1:5, valMat)
      if(thisArea == "Region")
        cex = 2.5
      else
        cex = 1
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        points(x,y,..., col=theseCols,cex=cex)
        abline(a = 0,b = 1,...)
      }
      
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseAreaColors
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseRadiusColors
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseUrbanColors
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      zlim = range(valMat)
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      width = 400 * (numberModels - 1)
      height = 400
      theseCols = theseUrbanColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
    } else {
      thesePointTypes = sapply(theseResults$urban, function(x) {ifelse(x, 15, 17)})
      
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        # abline(a = 0,b = 1, col=rgb(.4, .4, .4),...)
        points(x,y,..., col=theseCols, pch=thesePointTypes)
        abline(a = 0,b = 1, ...)
      }
      
      width = 4 * (numberModels - 1)
      height = 4
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPairPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(valMat)
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=.8)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPairPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(valMat)
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=.8)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(valMat)
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPairPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        ylab = ""
        
        plot(x, y, main="", type="n", xlab="", ylab="", 
             ylim=zlim, xlim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1)
        if(i == 1)
          mtext(side = 2, as.expression(modelNames[[numberModels]]), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=.8)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
    }
  }
}

# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
makeFinalPercentResidualPlot = function(dat, resultFilenames, modelClasses, modelVariations, 
                                        areaLevels=c("Region", "County", "Pixel"), 
                                        meanRange, meanTicks, meanTickLabels, 
                                        plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                        ncols=29, urbCols=makeGreenBlueSequentialColors(ncols), 
                                        kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), doModelClassPlots=FALSE, 
                                        areaMat=NULL, radiiMat=NULL, areaColors=makePurpleYellowSequentialColors(ncols)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    if(length(uniqueModelClasses) == 1)
      extraPlotNameRoot = uniqueModelClasses[1]
    else
      extraPlotNameRoot = ""
    
    # collect predictions and proportion urban per area/point
    predsList = list()
    widthsList = list()
    # modelNames = c()
    modelNames = list()
    for(j in 1:numberModels) {
      # modelNames = c(modelNames, paste(modelClasses[j], modelVariations[j]))
      modelNames = c(modelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
      
      # load this model and put results in the given column
      out = load(resultFilenames[j])
      theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
      
      # get proportion urban (will be the same for each model)
      colI = cut(as.numeric(theseResults$urban), breaks=seq(0 - .0001, 1 + .0001, l=ncols+1), labels=FALSE)
      theseCols = urbCols[colI]
      
      predsList = c(predsList, list(theseResults$preds))
      widthsList = c(widthsList, list(theseResults$Q90 - theseResults$Q10))
    }
    
    valMat = do.call("cbind", predsList)
    zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
    
    if(thisArea %in% c("Region", "County")) {
      # determine colorings by area, radius, and sparsity
      theseAreas = getArea(thisArea)
      theseRadii = getRadius(thisArea)
      theseSparsities = getAreaPerObservation(thisArea)
      
      colI = cut(theseAreas, breaks=seq(min(theseAreas) - .0001, max(theseAreas) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseAreaColors = areaColors[colI]
      colI = cut(theseRadii, breaks=seq(min(theseRadii) - .0001, max(theseRadii) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseRadiusColors = areaColors[colI]
      colI = cut(theseSparsities, breaks=seq(min(theseSparsities) - .0001, max(theseSparsities) + .0001, l=length(areaColors)+1), labels=FALSE)
      theseSparsityColors = areaColors[colI]
      theseUrbanColors = theseCols
      sparsityTicks = getAreaPerObservationTicks(thisArea)
      sparsityTickLabels = as.character(sparsityTicks)
      sparsityTicks = log(sparsityTicks)
      
      # valMat = rbind(1:5, valMat)
      if(thisArea == "Region")
        cex = 2.5
      else
        cex = 1
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        # else
        #   xlim = zlim2
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        points(x,y,..., col=theseCols,cex=cex)
        abline(h=0, lty=2,...)
      }
      
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseAreaColors
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseRadiusColors
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotSparsityCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      theseCols = theseSparsityColors
      
      for(i in 1:(numberModels-1)) {
        y = valMat[,i]
        x = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseSparsities), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = expression("Sparsity (km"^2 ~" Per Cluster)"), legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(at=sparsityTicks, labels=sparsityTickLabels, cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseUrbanColors
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      width = 400 * (numberModels - 1)
      height = 400
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseSparsityColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsSparsityCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseSparsities), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = expression("Sparsity (km"^2 ~" Per Cluster)"), legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(at=sparsityTicks, labels=sparsityTickLabels, cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      valMat = do.call("cbind", widthsList)
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      width = 400 * (numberModels - 1)
      height = 400
      theseCols = theseUrbanColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=ncols, legend.mar=4, col=urbCols, add=TRUE, 
                 legend.lab = "Urbanicity", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseAreaColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsAreaCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseAreas), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Area (sq. km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseRadiusColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsRadiusCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseRadii), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = "Radius (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
      
      theseCols = theseSparsityColors
      png(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsSparsityCol", plotNameRoot, extraPlotNameRoot, thisArea, ".png"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, pch=19)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      image.plot(legend.only = TRUE, zlim=range(theseSparsities), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                 legend.lab = expression("Sparsity (km"^2 ~" Per Cluster)"), legend.line=5.0, legend.width=3, legend.shrink=.9, 
                 legend.cex=1.5, axis.args=list(at=sparsityTicks, labels=sparsityTickLabels, cex.axis=1.5, tck=-1, hadj=-.1))
      dev.off()
    } else {
      urbanPointTypes = sapply(theseResults$urban, function(x) {ifelse(x, 15, 17)})
      thesePointTypes = urbanPointTypes
      theseUrbanColors = theseCols
      
      # determine colorings by distance to observation
      if(thisArea == "Pixel") {
        theseDistances = getPredictionDistance(doLog=FALSE)
        
        colI = cut(theseDistances, breaks=seq(min(theseDistances) - .0001, max(theseDistances) + .0001, l=length(areaColors)+1), labels=FALSE)
        theseDistanceColors = areaColors[colI]
        theseUrbanColors = theseCols
        distanceTicks = getPredictionDistanceTicks()
        distanceTickLabels = as.character(distanceTicks)
        distanceTicks = distanceTicks
        distancePointTypes = 19
        cex = 0.1
      } else if(thisArea == "Cluster") {
        cex = 0.8
      }
      
      my_line <- function(x,y,...){
        if(diff(range(x)) >= .04)
          xlim = zlim
        
        if(diff(range(y)) >= .04)
          ylim = zlim
        
        # abline(a = 0,b = 1, col=rgb(.4, .4, .4),...)
        points(x,y,..., col=theseCols, pch=thesePointTypes)
        abline(h=0, lty=2, ...)
      }
      
      width = 4 * (numberModels - 1)
      height = 4
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlot", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=cex)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      if(thisArea == "Pixel") {
        thesePointTypes = distancePointTypes
        theseCols = theseDistanceColors
        pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotDistCol", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
        par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
        
        zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
        for(i in 1:(numberModels-1)) {
          x = valMat[,i]
          y = valMat[,numberModels]
          y = 100 * (x-y)/y
          ylab = ""
          
          plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Estimates"), ylab="", 
               ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
          if(i == 1)
            mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
          my_line(x, y, lwd=1, cex=cex)
          mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
        }
        image.plot(legend.only = TRUE, zlim=range(theseDistances), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                   legend.lab = "Distance to Observation (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                   legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
        dev.off()
      }
      
      thesePointTypes = urbanPointTypes
      theseCols = theseUrbanColors
      valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=cex)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      if(thisArea == "Pixel") {
        thesePointTypes = distancePointTypes
        theseCols = theseDistanceColors
        valMat = do.call("cbind", lapply(1:length(widthsList), function(x) {widthsList[[x]] / predsList[[x]]}))
        pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotRelWidthsDistCol", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
        par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
        
        zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
        for(i in 1:(numberModels-1)) {
          x = valMat[,i]
          y = valMat[,numberModels]
          y = 100 * (x-y)/y
          ylab = ""
          
          plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Rel. Width"), ylab="", 
               ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
          if(i == 1)
            mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
          my_line(x, y, lwd=1, cex=cex)
          mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
        }
        image.plot(legend.only = TRUE, zlim=range(theseDistances), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                   legend.lab = "Distance to Observation (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                   legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
        dev.off()
      }
      
      thesePointTypes = urbanPointTypes
      theseCols = theseUrbanColors
      valMat = do.call("cbind", widthsList)
      zlim = range(sweep(valMat[,-ncol(valMat)], 1, valMat[,ncol(valMat)], function(x,y) {100*(x-y)/y}))
      pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidths", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
      par(mfrow=c(1, 3), oma=c(0,4,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
      
      for(i in 1:(numberModels-1)) {
        x = valMat[,i]
        y = valMat[,numberModels]
        y = 100 * (x-y)/y
        ylab = ""
        
        plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
             ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
        if(i == 1)
          mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
        my_line(x, y, lwd=1, cex=cex)
        mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
      }
      dev.off()
      
      if(thisArea == "Pixel") {
        thesePointTypes = distancePointTypes
        theseCols = theseDistanceColors
        pdf(file=paste0("Figures/", resultNameRoot, "/finalPercentResidualPlotWidthsDistCol", plotNameRoot, extraPlotNameRoot, thisArea, ".pdf"), width=width, height=height)
        par(mfrow=c(1, 3), oma=c(0,4,0,5), mar=c(5.1, 4.1, 4.1, 6))
        
        for(i in 1:(numberModels-1)) {
          x = valMat[,i]
          y = valMat[,numberModels]
          y = 100 * (x-y)/y
          ylab = ""
          
          plot(x, y, main="", type="n", xlab=bquote(.(modelNames[[numberModels]]) ~ " Width"), ylab="", 
               ylim=zlim, cex.lab=2, cex.main=2, cex.axis=2, asp=1, log="x")
          if(i == 1)
            mtext(side = 2, bquote("Pct. Diff. from " ~ .(modelNames[[numberModels]])), line = 4, cex=2)
          my_line(x, y, lwd=1, cex=cex)
          mtext(side = 3, as.expression(modelNames[[i]]), line = 1, cex=2)
        }
        image.plot(legend.only = TRUE, zlim=range(theseDistances), nlevel=ncols, legend.mar=4, col=areaColors, add=TRUE, 
                   legend.lab = "Distance to Observation (km)", legend.line=5.0, legend.width=3, legend.shrink=.9, 
                   legend.cex=1.5, axis.args=list(cex.axis=1.5, tck=-1, hadj=-.1))
        dev.off()
      }
    }
  }
}

plotCovariograms = function(dat, resultFilenames, modelClasses, modelVariations, 
                            varName="education", plotNameRoot="Education", resultNameRoot="Ed", 
                            cgramList=NULL, loadResults=FALSE, saveResults=!loadResults, col=NULL, 
                            doModelClassPlots=FALSE, lty=1, pch=19, hNuggetShift=-3) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  numberModels = length(resultFilenames)
  uniqueModelClasses = unique(modelClasses)
  
  if(is.null(col))
    col = rainbow(length(modelClasses))
  
  binomialFamily = grepl("LgtN", resultFilenames[1])
  
  # first get the covariograms if necessary
  if(is.null(cgramList)) {
    cgramList = list()
    for(j in 1:numberModels) {
      modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
      
      # load this model and get the covariogram if necessary
      print(paste0("Loading ", paste(modelClasses[j], modelVariations[j])))
      out = load(resultFilenames[j])
      if(!loadResults || !("cgram" %in% names(results))) {
        hyperDraws = results$fit$hyperMat
        browser()
        
        # determine if this model has a cluster effect
        thisModelClass = modelClasses[j]
        clusterEffect = grepl("includeClusterTRUE", resultFilenames[j])
        
        # hyperparameters will be drawn differently depending on the type of model
        if(thisModelClass == "SPDE") {
          # get hyperparameter draws
          effectiveRangeI = grepl("Range for field", rownames(hyperDraws))
          sdI = grepl("Stdev for field", rownames(hyperDraws))
          effectiveRangeVals = hyperDraws[effectiveRangeI,]
          varVals = hyperDraws[sdI,]^2
          if(clusterEffect)
            nuggetVarVals = 1/hyperDraws[3,]
          else
            nuggetVarVals = rep(0, ncol(hyperDraws))
          
          # get range of the data and the SPDE basis function mesh for which to compute the covariograms
          out = load(paste0("dataPointsKenya.RData"))
          xRangeDat = dataPointsKenya$xRange
          yRangeDat = dataPointsKenya$yRange
          mesh = results$fit$mesh
          
          # compute the covariance function for the different hyperparameter samples
          cgram = covarianceDistributionSPDE(effectiveRangeVals, varVals, nuggetVarVals, mesh, xRangeDat=xRangeDat, yRangeDat=yRangeDat)
        } else if(thisModelClass == "ELK") {
          # get lattice information object, determine whether it's a separateRange model
          latInfo = results$fit$latInfo
          nLayer = length(latInfo)
          separateRanges = grepl("separateRangesTRUE", resultFilenames[j])
          if(separateRanges)
            alphaI = (1 + nLayer+1 + 1):(1 + nLayer+1 + nLayer-1)
          else
            alphaI = 4:(3+nLayer-1)
          zSamples = matrix(hyperDraws[alphaI,], ncol=length(alphaI))
          xSamples = t(matrix(apply(zSamples, 1, multivariateExpit), ncol=length(alphaI)))
          xSamples = rbind(xSamples, 1-colSums(xSamples))
          
          # get hyperparameter draws
          nuggetVarVals = rep(0, ncol(hyperDraws))
          if(separateRanges) {
            kappaVals = sweep(sqrt(8)/exp(hyperDraws[2:(nLayer+1),]), 1, sapply(latInfo, function(x) {x$latWidth}), "*")
            rhoVals = exp(hyperDraws[nLayer+2,])
          } else {
            kappaVals = sqrt(8)/exp(hyperDraws[2,]) * latInfo[[1]]$latWidth
            rhoVals = exp(hyperDraws[3,])
          }
          alphaMat = xSamples
          
          # compute the covariance function for many different hyperparameter samples
          out = covarianceDistributionLKINLA(latInfo, kappaVals, rhoVals, nuggetVarVals, alphaMat)
        } else {
          stop(paste0("Unrecognized model class: ", thisModelClass))
        }
        
        # save results if necessary
        if(saveResults) {
          results$cgram = cgram
          save(results, file=resultFilenames[j])
        }
      } else {
        cgram = results$cgram
      }
      
      # append to our list of covariograms
      cgramList = c(cgramList, list(cgram))
    }
  }
  
  ##### if there are multiple unique model classes, call this function on each one individually 
  ##### as well as on the combination
  if(length(uniqueModelClasses) != 1 && doModelClassPlots) {
    for(i in 1:length(uniqueModelClasses)) {
      thisModelClass = uniqueModelClasses[i]
      thisI = modelClasses == thisModelClass
      if(sum(thisI) == 1)
        next
      
      thisResultFilenames = resultFilenames[thisI]
      thisModelClasses = modelClasses[thisI]
      thisModelVariations = modelVariations[thisI]
      thiscgramList = cgramList[thisI]
      thisColors = col[thisI]
      thispch = pch[thisI]
      thislty = lty[thisI]
      plotCovariograms(dat, thisResultFilenames, thisModelClasses, thisModelVariations, 
                       varName, plotNameRoot, resultNameRoot, thiscgramList, 
                       loadResults=TRUE, saveResults=FALSE, col=thisColors, lty=thislty, 
                       pch=thispch, hNuggetShift=hNuggetShift)
    }
  }
  
  # reorder the models if necessary so that we have 2 blocks of 4 (note that these permutations are 
  # the inverse permutations of themselves)
  # if(numberModels == 8) {
  #   reordering = c(t(rbind(c(1, 2, 5, 6), 
  #                          c(3, 4, 7, 8))))
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # } else if(numberModels == 4) {
  #   reordering = c(t(rbind(c(1, 3), 
  #                          c(2, 4))))
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # }
  # if(numberModels == 6) {
  #   reordering = c(t(rbind(c(1, 3, 4), 
  #                          c(2, 5, 6)))) # urban effects on bottom, first column is SPDE
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # } else {
  #   reordering = 1:numberModels # no need to reorder
  #   resultFilenames = resultFilenames[reordering]
  #   modelClasses = modelClasses[reordering]
  #   modelVariations = modelVariations[reordering]
  #   cgramList = cgramList[reordering]
  # }
  
  # get range of covariogram values
  yRange = c()
  yRangeNoCIs = c()
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    
    thiscgram = cgramList[[j]]
    yRange = range(c(yRange, thiscgram$cov, thiscgram$lowerCov, thiscgram$upperCov))
    yRangeNoCIs = range(c(yRangeNoCIs, thiscgram$cov))
  }
  
  # modify plot file names if necessary
  if(length(uniqueModelClasses) == 1)
    extraPlotNameRoot = uniqueModelClasses[1]
  else
    extraPlotNameRoot = ""
  
  allModelNames = list()
  for(j in 1:numberModels)
    allModelNames = c(allModelNames, bquote(.(modelClasses[j])[.(modelVariations[j])]))
  
  ##### Plot everything
  print("Plotting covariograms...")
  # col = rainbow(numberModels)
  # width = 4 * ceiling(numberModels/2)
  
  pdf(file=paste0("Figures/", resultNameRoot, "/covariogramsAll", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the covariograms together
  # allModelNames = paste(modelClasses, modelVariations)
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      covMean = thiscgram$cov[sortI]
      upperCov=thiscgram$upperCov[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      lowerCov=thiscgram$lowerCov[1,sortI]
    } else {
      covMean = thiscgram$cov[sortI]
      upperCov=thiscgram$upperCov[sortI]
      lowerCov=thiscgram$lowerCov[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.4, main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRange, xlim=range(d), col=col[j])
        lines(d[!d0], covMean[!d0], col=col[j], lwd=2)
      } else {
        points(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.4, col=col[j])
        lines(d[!d0], covMean[!d0], col=col[j], lwd=2)
      }
      
      points(hNuggetShift, mean(lowerCov[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], lowerCov[!d0], col=col[j], lty=2, lwd=2)
      points(hNuggetShift, mean(upperCov[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], upperCov[!d0], col=col[j], lty=2, lwd=2)
    } else {
      if(j == 1) {
        plot(d, covMean, type="l", main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRange, col=col[j], lwd=2)
      } else {
        lines(d, covMean, col=col[j], lwd=2)
      }
      
      lines(d, lowerCov, lty=2, col=col[j], lwd=2)
      lines(d, upperCov, lty=2, col=col[j], lwd=2)
    }
    
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=1, col=col, cex=ifelse(numberModels >= 5, .7, 1))
  dev.off()
  
  pdf(file=paste0("Figures/", resultNameRoot, "/covariogramsAllNoCIs", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the covariograms together
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      covMean = thiscgram$cov[sortI]
      # upperCov=thiscgram$upperCov[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      # lowerCov=thiscgram$lowerCov[1,sortI]
    } else {
      covMean = thiscgram$cov[sortI]
      # upperCov=thiscgram$upperCov[sortI]
      # lowerCov=thiscgram$lowerCov[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.4, main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRangeNoCIs, xlim=range(d), col=col[j], lty=lty[j])
        lines(d[!d0], covMean[!d0], col=col[j], lty=lty[j], lwd=2)
      } else {
        points(hNuggetShift, mean(covMean[d0]), pch=pch[j], cex=.4, col=col[j], lty=lty[j])
        lines(d[!d0], covMean[!d0], col=col[j], lty=lty[j], lwd=2)
        }
    } else {
      if(j == 1) {
        plot(d, covMean, type="l", main="", xlab="Distance (km)", ylab="Covariance", 
             ylim=yRangeNoCIs, col=col[j], lty=lty[j], lwd=2)
      } else {
        lines(d, covMean, col=col[j], lty=lty[j], lwd=2)
      }
    }
    
    # lines(d, lowerCov, lty=2, col=col[j])
    # lines(d, upperCov, lty=2, col=col[j])
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=lty, col=col, cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
  print("Plotting correlograms...")
  
  pdf(file=paste0("Figures/", resultNameRoot, "/correlogramsAll", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the correlograms together
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      corMean = thiscgram$cor[sortI]
      upperCor=thiscgram$upperCor[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      lowerCor=thiscgram$lowerCor[1,sortI]
    } else {
      corMean = thiscgram$cor[sortI]
      upperCor=thiscgram$upperCor[sortI]
      lowerCor=thiscgram$lowerCor[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.4, main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), xlim=range(d), col=col[j])
        lines(d[!d0], corMean[!d0], col=col[j], lwd=2)
      } else {
        points(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.4, col=col[j])
        lines(d[!d0], corMean[!d0], col=col[j], lwd=2)
      }
      
      points(hNuggetShift, mean(lowerCor[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], lowerCor[!d0], col=col[j], lty=2, lwd=2)
      points(hNuggetShift, mean(upperCor[d0]), pch=pch[j], cex=.2, col=col[j], lty=2)
      lines(d[!d0], upperCor[!d0], col=col[j], lty=2, lwd=2)
    } else {
      if(j == 1) {
        plot(d, corMean, type="l", main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), col=col[j], lwd=2)
      } else {
        lines(d, corMean, col=col[j], lwd=2)
      }
      
      lines(d, lowerCor, lty=2, col=col[j], lwd=2)
      lines(d, upperCor, lty=2, col=col[j], lwd=2)
    }
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=1, col=col, cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
  pdf(file=paste0("Figures/", resultNameRoot, "/correlogramsAllNoCIs", plotNameRoot, extraPlotNameRoot, ".pdf"), width=5, height=5)
  
  # plot the correlograms together
  for(j in 1:numberModels) {
    modelName = bquote(.(modelClasses[j])[.(modelVariations[j])])
    thiscgram = cgramList[[j]]
    d = thiscgram$d
    sortI = sort(d, index.return=TRUE)$ix
    d = d[sortI]
    
    if(modelClasses[j] == "SPDE") {
      corMean = thiscgram$cor[sortI]
      # upperCor=thiscgram$upperCor[1,sortI] # second row is the 95% CI, while the first is the 80% CI
      # lowerCor=thiscgram$lowerCor[1,sortI]
    } else {
      corMean = thiscgram$cor[sortI]
      # upperCor=thiscgram$upperCor[sortI]
      # lowerCor=thiscgram$lowerCor[sortI]
    }
    
    if(binomialFamily) {
      d0 = d == 0
      if(j == 1) {
        plot(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.4, main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), xlim=range(d), col=col[j], lty=lty[j])
        lines(d[!d0], corMean[!d0], col=col[j], lty=lty[j], lwd=2)
      } else {
        points(hNuggetShift, mean(corMean[d0]), pch=pch[j], cex=.4, col=col[j], lty=lty[j])
        lines(d[!d0], corMean[!d0], col=col[j], lty=lty[j], lwd=2)
      }
    } else {
      if(j == 1) {
        plot(d, corMean, type="l", main="", xlab="Distance (km)", ylab="Correlation", 
             ylim=c(0,1), col=col[j], lty=lty[j], lwd=2)
      } else {
        lines(d, corMean, col=col[j], lty=lty[j], lwd=2)
      }
    }
    
    # lines(d, lowerCor, lty=2, col=col[j])
    # lines(d, upperCor, lty=2, col=col[j])
    # lines(d, mixtureCovFun(d), col="green")
    # legend("topright", c("Truth", "Estimate", "80% CI"), lty=c(1, 1, 2), col=c("green", "black", "black"))
    
  }
  legend("topright", as.expression(allModelNames), lty=lty, col=col, cex=ifelse(numberModels >= 5, .5, 1))
  dev.off()
  
}

printModelPredictionTables = function(dataType=c("ed", "mort"), resultFilenames=NULL, modelClasses=c("SPDE", "ELK"), modelVariations=c("U", "U"), 
                                      nDigitsPredictions=3, nDigitsParameters=3, byRow=FALSE, areaLevels=c("Region", "County")) {
  
  dataType = match.arg(dataType)
  if(dataType == "mort") {
    out = load("../U5MR/kenyaData.RData")
    dat = mort
    resultNameRoot = "Mort"
  }
  else {
    out = load("../U5MR/kenyaDataEd.RData")
    dat = ed
    resultNameRoot = "Ed"
  }
  resultNameRootLower = tolower(resultNameRoot)
  
  if(is.null(resultFilenames)) {
    resultFilenames = c(
      paste0("savedOutput/", resultNameRoot, "/resultsSPDE", resultNameRootLower, "_urbanEffectTRUE_LgtN.RData"), 
      paste0("savedOutput/", resultNameRoot, "/resultsLKINLA", resultNameRootLower, "_separateRangesTRUE_urbanEffectTRUE_LgtN_noUrbanPrior.RData")
    )
  }
  
  ##### load estimates
  # construct a list of all estimates. First list index is model class, second index is model variation, 
  # third is aggregation level and parameters
  print("getting SPDE estimates...")
  out = load(resultFilenames[1])
  spdeResultsRegion = results$aggregatedResults$predictions$regionPredictions
  spdeResultsCounty = results$aggregatedResults$predictions$countyPredictions
  spdeParameterTable = results$aggregatedResults$parameterSummary
  
  ##### ELK estimates
  print("getting ELK estimates...")
  out = load(resultFilenames[2])
  lkinlaResultsRegion = results$aggregatedResults$predictions$regionPredictions
  lkinlaResultsCounty = results$aggregatedResults$predictions$countyPredictions
  lkinlaParameterTable = results$aggregatedResults$parameterSummary
  
  ##### now we construct the latex table for each areal level
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    # get the model results for this level of areal aggregation
    if(thisArea == "Region") {
      lkinlaResults = lkinlaResultsRegion
      spdeResults = spdeResultsRegion
    } else if(thisArea == "County") {
      lkinlaResults = lkinlaResultsCounty
      spdeResults = spdeResultsCounty
    } else
      stop(paste0("Unrecognized area name: ", thisArea))
    
    # construct result table
    tab = cbind(Estimates=spdeResults$preds, Q10=spdeResults$Q10, Q90=spdeResults$Q90)
    tab = cbind(tab, Estimates=lkinlaResults$preds, Q10=lkinlaResults$Q10, Q90=lkinlaResults$Q90)
    tab = format(tab, digits=nDigitsPredictions) # round to the nearest 1 child per 1000
    tab = cbind(as.character(spdeResults$areaName), tab)
    colnames(tab) = c(thisArea, rep(c("Est", "Q10", "Q90"), 2))
    rownames(tab) = NULL
    
    # print out relevant summary statistics about the predictions
    print(paste0("SPDE U range of ", thisArea, " predictions: ", diff(range(as.numeric(tab[,2])))))
    print(paste0("SPDE U median ", thisArea, " 80% CI width: ", median(as.numeric(tab[,4])-as.numeric(tab[,3]))))
    print(paste0("ELK U range of ", thisArea, " predictions: ", diff(range(as.numeric(tab[,5])))))
    print(paste0("ELK U median 80% ", thisArea, " CI width: ", median(as.numeric(tab[,7])-as.numeric(tab[,6]))))
    
    # print out the reformatted table
    require(kableExtra)
    print(paste0("Predictions at ", tolower(thisArea), " level:"))
    fullTab = tab %>%
      kable("latex", escape = F, booktabs = TRUE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
            align=c("l", rep("r", ncol(tab) - 1)), longtable=thisArea == "County", caption = thisArea) %>% 
      kable_styling(latex_options =c("repeat_header", "scale_down"))
    numberColumns = c(" "=1, "SPDE U"=3, "ELK U"=3)
    print(add_header_above(fullTab, numberColumns, italic=TRUE, bold=TRUE, escape=FALSE, line=TRUE))
  }
  
  
  ##### now construct the parameter tables
  ##### now print the parameter estimates:
  print("printing parameter estimates...")
  
  # SPDE
  
  parameters = round(spdeParameterTable, digits=nDigitsParameters)
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  browser()
  byRow = FALSE
  if(byRow == FALSE)
    nonRangePar = format(parameters[-nrow(parameters),], digits=nDigitsParameters, scientific=FALSE)
  else
    nonRangePar = t(apply(parameters[-nrow(parameters),], 1, format, digits=2, scientific=FALSE))
  formattedParameters = rbind(nonRangePar, format(parameters[nrow(parameters),], digits=0))
  rownames(formattedParameters) = c("Intercept", "Urban", "Total Var", "Spatial Var", "Cluster Var", "Total SD", "Spatial SD", "Cluster SD", "Range")
  
  # rename the columns
  colnames(formattedParameters) = c("Est", "SD", "Q10", "Q50", "Q90")
  
  # print(paste0("Parameter summary table for SPDE ", typeText, " model:"))
  # print(xtable(parameters, digit=3))
  parTable = formattedParameters
  parTable = cbind(rownames(parTable), parTable)
  rownames(parTable) = NULL
  
  # ELK
  
  parameters = round(lkinlaParameterTable, digits=nDigitsParameters)
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  
  if(byRow == FALSE)
    nonRangePar = format(parameters[-(9:10),], digits=nDigitsParameters, scientific=FALSE)
  else
    nonRangePar = t(apply(parameters[-(9:10),], 1, format, digits=2, scientific=FALSE))
  formattedParameters = rbind(nonRangePar[1:8,], format(parameters[9:10,], digits=0), nonRangePar[9:nrow(nonRangePar),])
  rownames(formattedParameters) = c("Intercept", "Urban", "Total Var", "Spatial Var", "Cluster Var", "Total SD", "Spatial SD", "Cluster SD", "Range1", "Range2", 
                                    "alpha1", "alpha2")
  
  # rename the columns
  colnames(formattedParameters) = c("Est", "SD", "Q10", "Q50", "Q90")
  
  # print(paste0("Parameter summary table for SPDE ", typeText, " model:"))
  # print(xtable(parameters, digit=3))
  
  formattedParameters = cbind(rownames(formattedParameters), formattedParameters)
  rownames(formattedParameters) = NULL
  colnames(parTable)[1] = "Parameter"
  colnames(formattedParameters)[1] = "Parameter"
  parTable = rbind(parTable, formattedParameters)
  
  fullTab = kable(parTable, "latex", booktabs = T, escape=FALSE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
                  align=c("l", rep("r", ncol(parTable) - 1)), longtable=FALSE, caption = "Parameter") %>% 
    kable_styling(latex_options="repeat_header")
  
  fullTab = fullTab %>% 
    pack_rows("SPDE U", 1, 9, escape=FALSE, bold=TRUE, italic=TRUE) %>% 
    pack_rows("ELK U", 10, nrow(parTable), escape=FALSE, bold=TRUE, italic=TRUE)
  print(fullTab)
}

# this is mostly a test function for plotting the predictions of a single model
# this function plots central estimates and credible interval widths over a map of Kenya 
# 4 plots are made by default: one for each area level
plotSingleModelPredictions = function(dat=NULL, results, modelName="", targetPop=c("women", "children"), 
                                      areaLevels=c("Region", "County", "Pixel", "Cluster"), 
                                      plotNameRoot="", inExampleFolder=FALSE, 
                                      meanRange=NULL, meanTicks=NULL, meanTickLabels=NULL, widthRange=NULL, widthTicks=NULL, widthTickLabels=NULL, 
                                      meanCols=makeRedBlueDivergingColors(64), 
                                      widthCols=makeBlueYellowSequentialColors(64), 
                                      kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  targetPop = match.arg(targetPop)
  if(targetPop == "women") {
    varName="education"
    plotNameRoot=paste0("Education", plotNameRoot)
    resultNameRoot="Ed"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaDataEd.RData")
      dat = mort
    }
  } else if(targetPop == "children") {
    varName="mort"
    plotNameRoot=paste0("Mort", plotNameRoot)
    resultNameRoot="Mort"
    if(is.null(dat)) {
      out = load("../U5MR/kenyaData.RData")
      dat = mort
    }
  }
  if(!inExampleFolder)
    resultNameRoot = ""
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  # load the fine grid used to approximate continuous prediction (the actual population density values don't matter)
  print("Loading prediction grid and shapefiles...")
  load("../U5MR/popGrid.RData")
  
  # load shape files for plotting
  require(maptools)
  regionMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  out = load("../U5MR/adminMapData.RData")
  kenyaMap = adm0
  countyMap = adm1compressed
  
  ##### central estimates and credible interval widths
  
  for(i in 1:length(areaLevels)) {
    thisArea = areaLevels[i]
    
    # get map for this areal aggregation level
    if(thisArea == "Cluster") {
      next
    } else if(thisArea == "Pixel") {
      thisMap = countyMap
    } else if(thisArea == "Region") {
      thisMap = regionMap
    } else if(thisArea == "County") {
      thisMap = countyMap
    }
    
    print("Plotting central estimates and credible interval widths...")
    numberModels=1
    width = 400 * numberModels
    png(file=paste0("Figures/", resultNameRoot, "/", plotNameRoot, thisArea, ".png"), width=width, height=1000)
    
    if(thisArea %in% c("Region", "County"))
      par(mfrow=c(2,numberModels))
    else
      par(mfrow=c(2,numberModels), oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
    
    # first load in the predictions
    predictionList = list()
    widthList = list()
    
    # load this model and plot results in the given column
    theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
    predictionList = c(predictionList, list(theseResults$preds))
    widthList = c(widthList, list(theseResults$Q90 - theseResults$Q10))
    
    # plot the predictions
    
    # load this model and plot results in the given column
    # out = load(resultFilenames[j])
    # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
    
    if(is.null(meanRange))
      thisMeanRange = logit(range(predictionList[[1]]))
    else
      thisMeanRange = logit(meanRange)
    if(is.null(meanTicks))
      thisMeanTicks = NULL
    else
      thisMeanTicks = logit(meanTicks)
    
    if(thisArea %in% c("Region", "County")) {
      
      plotMapDat(plotVar=predictionList[[1]], new = TRUE, 
                 main=bquote(.(modelName) ~ " estimates"), scaleFun=logit, scaleFunInverse=expit, 
                 cols=meanCols, zlim=thisMeanRange, ticks=meanTicks, tickLabels=meanTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange)
    } else if(thisArea == "Pixel"){
      plot(cbind(popGrid$lon, popGrid$lat), type="n", main=bquote(.(modelName) ~ " estimates"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
      quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(predictionList[[1]]), 
                 nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=thisMeanRange)
      plotMapDat(mapDat=thisMap, lwd=.5)
      points(dat$lon, dat$lat, pch=".")
      if(is.null(thisMeanTicks) || is.null(meanTickLabels)) {
        thisMeanTicks = logit(pretty(expit(thisMeanRange), n=5))
        meanTickLabels = as.character(expit(thisMeanTicks))
        image.plot(zlim=thisMeanRange, nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                   col=meanCols, add = TRUE, axis.args=list(at=thisMeanTicks, labels=meanTickLabels), legend.mar = 0)
      } else {
        image.plot(zlim=thisMeanRange, nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
                   col=meanCols, add = TRUE, axis.args=list(at=thisMeanTicks, labels=meanTickLabels), legend.mar = 0)
      }
      # image.plot(zlim=range(logit(meanRange)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
      #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
    }
    
    # plot the credible interval widths
    
    # plot results in the given column
    # theseResults = results$aggregatedResults$predictions[[paste0(tolower(thisArea), "Predictions")]]
    
    if(is.null(widthRange))
      thisWidthRange = log(range(widthList[[1]]))
    else
      thisWidthRange = log(widthRange)
    if(is.null(widthTicks))
      thisWidthTicks = NULL
    else
      thisWidthTicks = log(widthTicks)
    
    if(thisArea %in% c("Region", "County")) {
      plotMapDat(plotVar=widthList[[1]], new = TRUE, 
                 main=bquote(.(modelName) ~ " 80% CI width"), scaleFun=log, scaleFunInverse=exp, 
                 cols=widthCols, zlim=thisWidthRange, ticks=widthTicks, tickLabels=widthTickLabels, 
                 xlim=kenyaLonRange, ylim=kenyaLatRange)
    } else if(thisArea == "Pixel") {
      
      plot(cbind(popGrid$lon, popGrid$lat), type="n", main=bquote(.(modelName) ~ " 80% CI width"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
      quilt.plot(cbind(popGrid$lon, popGrid$lat), log(widthList[[1]]), 
                 nx=150, ny=150, add.legend=FALSE, add=TRUE, col=widthCols, zlim=thisWidthRange)
      plotMapDat(mapDat=thisMap, lwd=.5)
      points(dat$lon, dat$lat, pch=".")
      if(is.null(thisWidthTicks) || is.null(widthTickLabels)) {
        thisWidthTicks = log(pretty(exp(thisWidthRange), n=5))
        widthTickLabels = as.character(exp(thisWidthTicks))
        image.plot(zlim=thisWidthRange, nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                   col=widthCols, add = TRUE, axis.args=list(at=thisWidthTicks, labels=widthTickLabels), legend.mar = 0)
      } else {
        image.plot(zlim=thisWidthRange, nlevel=length(widthCols), legend.only=TRUE, horizontal=FALSE,
                   col=widthCols, add = TRUE, axis.args=list(at=thisWidthTicks, labels=widthTickLabels), legend.mar = 0)
      }
    }
    dev.off()
  }
}

makeRankPlots = function(postSampleMat, admin=c("admin1","admin2"), highestLowestN=5, 
                         savePlots=FALSE, plotNameSuffix="", CI=.8) {
  admin = match.arg(admin)
  
  if(admin == "admin1") {
    map_shp = adm1compressed
    admin_name_dt = as.data.table(sort(map_shp$NAME_1))
    admin_name_dt$toPlot <- sort(map_shp$NAME_1)
  } else if (admin == "admin2") {
    map_shp = adm2compressed
    sortI = sort(map_shp$CONSTITUEN, index.return=TRUE)$ix
    admin_name_dt = as.data.table(sort(map_shp$CONSTITUEN))
    admin_name_dt$toPlot <- paste0(map_shp$CONSTITUEN[sortI], ", ", map_shp$COUNTY_NAM[sortI])
  }
  
  lowerProb = (1 - CI)/2
  upperProb = 1 - lowerProb
  
  #### rank postsamps ####
  
  rank_mt <- apply(postSampleMat, 2, rank)
  
  pred_dt <- admin_name_dt
  pred_dt[, "ID"] <- 1:nrow(pred_dt)
  pred_dt[, "avg_rank"] <- apply(rank_mt, 1, mean)
  pred_dt[, "low_rank"] <- apply(rank_mt, 1, quantile, probs=lowerProb)
  pred_dt[, "high_rank"] <- apply(rank_mt, 1, quantile, probs=upperProb)
  
  # top 5 state
  if(savePlots) {
    pdf(paste0(figDirectory, "/exploratoryAnalysis/rank_",
               admin, "_rankhigh", highestLowestN, plotNameSuffix, ".pdf"),
        width = 2.5, height = 5)
  }
  
  par(mar = c(2.5, 1, 2, 1), mfrow = c(highestLowestN, 1))
  
  pred_dt_order <- pred_dt[order(avg_rank)]
  
  for (i in 1:highestLowestN){
    
    id <- pred_dt_order[i, ID]
    name <- pred_dt_order[i, toPlot]
    
    rank_vt <- rank_mt[id, ]
    rank_vt <- ifelse(rank_vt <= 10, rank_vt, 10)
    
    avg_rank <- pred_dt_order[i, avg_rank]
    low_rank <- pred_dt_order[i, low_rank]
    high_rank <- pred_dt_order[i, high_rank]
    
    ranktable <- as.data.table(table(rank_vt))
    ranktable <- merge(data.table(rank = as.character(1:10)), ranktable, 
                       by.x = "rank", by.y = "rank_vt", all.x = T)
    ranktable[, "rank" := as.integer(rank)]
    ranktable <- ranktable[order(rank)]
    ranktable[is.na(N), "N"] <- 0
    
    barplot(ranktable$N, width = 0.825, 
            xlim = c(10, 0), xlab = "", ylab = "",
            main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1), 
                          " (", format(round(low_rank, 1), nsmall = 1), ", ", 
                          format(round(high_rank, 1), nsmall = 1), ")"),
            xaxt = "n", yaxt = "n", col = "#31a354", border = F,
            cex.main = 0.75)
    axis(1, at = 10:1-0.5, labels = c("10+", as.character(9:1)), tick = F)
  }
  if(savePlots) {
    dev.off()
  }
  
  
  # bottom 5 state
  if(savePlots) {
    pdf(paste0(figDirectory, "/exploratoryAnalysis/rank_",
               admin, "_ranklow", highestLowestN, plotNameSuffix, ".pdf"),
        width = 2.5, height = 5)
  }
  
  par(mar = c(2.5, 1, 2, 1), mfrow = c(highestLowestN, 1))
  
  pred_dt_order <- pred_dt[order(-avg_rank)]
  
  for (i in 1:highestLowestN){
    
    id <- pred_dt_order[i, ID]
    name <- pred_dt_order[i, toPlot]
    
    avg_rank <- pred_dt_order[i, avg_rank]
    low_rank <- pred_dt_order[i, low_rank]
    high_rank <- pred_dt_order[i, high_rank]
    
    rank_vt <- rank_mt[id, ]
    rank_vt <- ifelse(rank_vt >= (nrow(pred_dt_order)-9), rank_vt, (nrow(pred_dt_order)-9))
    
    ranktable <- as.data.table(table(rank_vt))
    ranktable <- merge(data.table(rank = as.character(nrow(pred_dt_order):(nrow(pred_dt_order)-9))), ranktable, 
                       by.x = "rank", by.y = "rank_vt", all.x = T)
    ranktable[, "rank" := as.integer(rank)]
    ranktable <- ranktable[order(rank)]
    ranktable[is.na(N), "N"] <- 0
    
    barplot(ranktable$N, width = 0.825, 
            xlim = c(10, 0), xlab = "", ylab = "",
            main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1), 
                          " (", format(round(low_rank, 1), nsmall = 1), ", ", 
                          format(round(high_rank, 1), nsmall = 1), ")"),
            xaxt = "n", yaxt = "n", col = "#fc4e2a", border = F,
            cex.main = 0.75)
    axis(1, at = 10:1-0.5, labels = c(as.character(nrow(pred_dt_order):(nrow(pred_dt_order)-8)), paste0(nrow(pred_dt_order)-9, "-")), tick = F)
  }
  if(savePlots) {
    dev.off()
  }
  
  
  #### all states hist ####
  rowcount <- ceiling(nrow(pred_dt_order)/3)
  pdf(paste0(figDirectory, "/exploratoryAnalysis/rank_",
             admin, "_rankall", plotNameSuffix, ".pdf"),
      width = 15, height = rowcount*2)
  
  par(mar = c(2.5, 1, 2, 1), mfcol = c(rowcount, 3))
  
  pred_dt_order <- pred_dt[order(avg_rank)]
  
  for (i in 1:nrow(pred_dt_order)){
    # i <- 1
    
    id <- pred_dt_order[i, ID]
    name <- pred_dt_order[i, toPlot]
    
    rank_vt <- rank_mt[id, ]
    
    avg_rank <- pred_dt_order[i, avg_rank]
    low_rank <- pred_dt_order[i, low_rank]
    high_rank <- pred_dt_order[i, high_rank]
    
    ranktable <- as.data.table(table(rank_vt))
    ranktable <- merge(data.table(rank = as.character(1:nrow(pred_dt_order))), ranktable, 
                       by.x = "rank", by.y = "rank_vt", all.x = T)
    ranktable[, "rank" := as.integer(rank)]
    ranktable <- ranktable[order(rank)]
    ranktable[is.na(N), "N"] <- 0
    
    barplot(ranktable$N, width = 0.825, 
            xlim = c(nrow(pred_dt_order), 0), xlab = "", ylab = "",
            main = paste0(name, "\nER = ", format(round(avg_rank, 1), nsmall = 1), 
                          " (", format(round(low_rank, 1), nsmall = 1), ", ", 
                          format(round(high_rank, 1), nsmall = 1), ")"),
            xaxt = "n", yaxt = "n", col = "#08519c", border = F,
            cex.main = 0.75)
    axis(1, at = nrow(pred_dt_order):1-0.5, labels = as.character(nrow(pred_dt_order):1), tick = F)
  }
  dev.off()
  
  
  invisible(NULL)
}













##### END OF MAIN PLOTTING FUNCTIONS #####

##### The rest of the functions in this script are utility functions for plotting

makeRedBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
  else
    scale_colour_continuous_sequential(h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3, n_interp=n)
}

makePurpleRedSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6)
  else
    scale_colour_continuous_sequential(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, n_interp=n)
}

makePurpleRedDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # if(!ggplot)
  #   sequential_hcl(n, h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6)
  # else
  #   scale_colour_continuous_sequential(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, n_interp=n)
  
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot) {
      divergingx_hcl(n, h1=-87, h3=10, c1=65, c3=79, l1=13, l3=60, p1=3, p2=1.6, h2=0, l2=100, c2=100, p3=3, p4=1.6, rev=rev)
      # diverging_hcl(n, h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, rev=rev)
    }
    else {
      stop("ggplot doesn't work with makePurpleRedDivergingColors")
      scale_colour_continuous_diverging(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, rev=rev, n_interp=n)
    }
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makePurpleRedDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      stop("ggplot doesn't work with makePurpleRedDivergingColors")
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=-87, h2=10, c1=65, c2=79, l1=13, l2=60, p1=3, p2=1.6, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeGreenBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  else
    scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
}

makeGreenBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE, p1=1) {
  # library("colorspace")
  # pal <-choose_palette()
  # if(!ggplot)
  #   sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  # else
  #   scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
  
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev)
    else
      scale_colour_continuous_diverging(h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeGreenBlueDivergingColors(totalColors, rev=rev, p1=p1)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=128, h2=250, c1=100, l1=71, l2=95, p1=p1, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeBlueGoldDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE, p1=1) {
  # library("colorspace")
  # pal <-choose_palette()
  # if(!ggplot)
  #   sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
  # else
  #   scale_colour_continuous_sequential(h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2, n_interp=n)
  
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=265, h2=74, c1=80, l1=66, l2=94, p1=p1, rev=rev)
    else
      scale_colour_continuous_diverging(h1=265, h2=74, c1=80, l1=66, l2=94, p1=p1, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeBlueGoldDivergingColors(totalColors, rev=rev, p1=p1)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=265, h2=74, c1=80, l1=66, l2=94, p1=p1, rev=rev, n_interp=n, mid=center)
    }
  }
}

makePurpleYellowSequentialColors = function(n, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev)
  else
    scale_colour_continuous_sequential(h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev, n_interp=n)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(n_interp=n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown && totalMissingColors > 0)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(n.interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev, mid=center)
    }
  }
}

makeBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE, n_interp=n)
}

makeGreenSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=rev)
  else
    scale_colour_continuous_sequential(h1=128, c1=100, l1=72, l2=95, p1=1.0, rev=rev, n_interp=n)
}

makeYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=86, c1=100, l1=70, l2=95, p1=1.2, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=86, c1=100, l1=70, l2=95, p1=1.2, rev=TRUE, n_interp=n)
}

makeBlueGreenYellowSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, rev=rev)
  else
    scale_colour_continuous_sequential(h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, n_interp=n, rev=rev)
}

makeRedYellowBlueColors = function(n, ggplot=FALSE) {
  if(!ggplot)
    divergingx_hcl(n, palette="RdYlBu")
  else
    scale_colour_continuous_sequential(palette="RdYlBu", n_interp=n)
}

makeRedYellowBlueColorsMod = function(n, ggplot=FALSE) {
  # Original parameters:
  # h1 h2  h3  c1 c2 c3 l1 l2 l3  p1
  # 10 85 260 105 45 70 35 98 35 1.5
  if(!ggplot) {
    divergingx_hcl(n,
                   h1=10, h2=85, h3=260,
                   c1=105, c2=45, c3=70,
                   l1=35, l2=98, l3=35, 
                   p1=1.5, p2=1.2, p3=.6, p4=1.2, 
                   cmax1=150, cmax2=10)
  }
  else
    scale_colour_continuous_sequential(palette="RdYlBu", n_interp=n)
}

makeRedYellowBlueColorsMod2 = function(n, ggplot=FALSE) {
  # Original parameters:
  # h1=10; h2=85; h3=260
  # c1=105; c2=45; c3=70
  # l1=35; l2=98; l3=35
  # p1=1.5; p2=1.2; p3=.6; p4=1.2
  # cmax1=150; cmax2=10
  
  h1=15; h2=79; h3=250
  c1=100; c2=100; c3=40
  l1=60; l2=70; l3=33
  p1=1.2; p3=.5; p4=1
  # 2: p1=.5; p2=1
  
  if(!ggplot) {
    # divergingx_hcl(n, 
    #                h1=10, h2=85, h3=260, 
    #                c1=105, c2=45, c3=70, 
    #                l1=80, l2=55, l3=30)
    if((n %% 2) == 0) {
      lower = sequential_hcl(h1=h1, h2=h2, c1=c1, c2=c2, 
                             l1=l1, l2=l2, p1=p1, n)[seq(1, n, by=2)]
      upper = sequential_hcl(h1=h2, h2=h3, c1=c2, c2=c3, 
                             l1=l2, l2=l3, p1=p3, p2=p4, n)[seq(2, n, by=2)]
    } else {
      lower = sequential_hcl(h1=h1, h2=h2, c1=c1, c2=c2, 
                             l1=l1, l2=l2, p1=p1, ceiling(n/2))
      upper = sequential_hcl(h1=h2, h2=h3, c1=c2, c2=c3, 
                             l1=l2, l2=l3, p1=p3, p2=p4, ceiling(n/2))[-1]
    }
    c(lower, upper)
  }
  else
    stop("ggplot version not supported")
}

makeRedYellowSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=15, h2=79, c1=100, c2=72, l1=40, l2=90, p1=1.2)
  else
    scale_colour_continuous_sequential(h1=15, h2=79, c1=100, c2=52, l1=55, l2=95, p1=1.2, n_interp=n)
}

makeRedGreenDivergingColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
  else
    scale_colour_continuous_sequential(h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5, n_interp=n)
}

# NOTE: this returns only hex colors
makePurpleGreenDivergingColors = function(n, rev=FALSE) {
  require(inlmisc)
  if(!rev) {
    as.character(inlmisc::GetColors(n, scheme = "PRGn"))
  } else {
    rev(as.character(inlmisc::GetColors(n, scheme = "PRGn")))
  }
}

combineTwoScales = function(n, scale1, scale2, args1, args2) {
  if(n %% 2 == 0) {
    n1 = n2 = n/2
  } else {
    n1 = ceiling(n/2)
    n2 = floor(n/2)
  }
  
  c(do.call(scale1, c(args1, list(n=n1))), 
    do.call(scale2, c(args2, list(n=n2))))
}

makeDivergingScale = function(n, scale, ...) {
  do.call("combineTwoScales", list(n=n, scale1=scale, scale2=scale, args1=list(...), args2=list(...)))
}

# centers a color scale at its midpoint. Returns vector of the centered color scale
# colScale a function taking 'n' as input and return a color scale centered in the middle
# ...: other arguments to colScale, such as 'rev'
centerColorScale = function(n, vals=NULL, valRange=NULL, center, colScale, scaleFun=function(x) {x}, 
                            ...) {
  
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  valRange = scaleFun(valRange)
  center = scaleFun(center)
  
  propUp = (valRange[2] - center) / diff(valRange)
  propDown = 1 - propUp
  totalColors = ceiling(2 * max(propUp, propDown) * n)
  tempColors = do.call(colScale, c(list(totalColors), list(...)))
  totalMissingColors = totalColors - n
  
  if(propUp >= propDown && totalMissingColors > 0)
    tempColors[-(1:totalMissingColors)]
  else
    tempColors[1:n]
}

# given continuous color scale and range, chooses colors based on a set of values
getColorsFromScale = function(vals, valRange=NULL, center=NULL, cols, scaleFun=function(x) {x}, 
                              forceValuesInRange=FALSE) {
  
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  if(forceValuesInRange) {
    vals[vals < valRange[1]] = valRange[1]
    vals[vals > valRange[2]] = valRange[2]
  }
  
  valRange = scaleFun(valRange)
  vals = scaleFun(vals)
  vals = vals - valRange[1]
  vals = vals/(valRange[2] - valRange[1])
  
  if(!is.null(center)) {
    center = scaleFun(center)
    n = length(cols)
    
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = cols
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
    
    cols = tempColors
  }
  
  col = cols[round(vals*(length(cols)-1))+1]
  
  col
}

plotWithColor = function(x, y, z, zlim=NULL, colScale=tim.colors(), 
                         legend.mar=7, new=TRUE, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                         n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, legend.width=1.2, addColorBar=TRUE, 
                         legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, orderI=NULL, 
                         ordering=c("none", "increasing", "decreasing"), colorName = c("col", "bg"), ...) {
  ordering = match.arg(ordering)
  colorName = match.arg(colorName)
  
  # remove NA points
  nas = is.na(x) | is.na(y) | is.na(z)
  if(any(nas)) {
    warning("Removing NAs")
    x = x[!nas]
    y = y[!nas]
    z = z[!nas]
  }
  
  # do setup for ploting data if necessary
  if(is.null(zlim)) {
    nas = !is.finite(scaleFun(z))
    zlim = range(z[!nas])
  }
  
  # order the plotting of the points
  if(is.null(orderI)) {
    if(ordering == "increasing") {
      orderI = sort(z, index.return=TRUE)$ix
    } else if(ordering == "decreasing") {
      orderI = sort(z, decreasing=TRUE, index.return=TRUE)$ix
    } else {
      orderI = 1:length(z)
    }
  }
  x = x[orderI]
  y = y[orderI]
  z = z[orderI]
  
  # if(forceColorsInRange) {
  #   z[z > zlim[2]] = zlim[2]
  #   z[z < zlim[1]] = zlim[1]
  # }
  
  # get colors of points
  cols = getColorsFromScale(z, zlim, cols=colScale, scaleFun=scaleFun, 
                            forceValuesInRange=forceColorsInRange)
  
  # generate new plot if necessary
  # browser()
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    if(colorName == "col") {
      do.call("plot", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("plot", c(list(x=x, y=y, bg=cols), list(...)))
    }
  } else {
    if(colorName == "col") {
      do.call("points", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("points", c(list(x=x, y=y, bg=cols), list(...)))
    }
  }
  
  if(addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(tickLabels))
      setTickLabels = TRUE
    
    if(is.null(ticks)) {
      if(setTickLabels)
        tickLabels = pretty(zlim, n=n.ticks, min.n=min.n)
      ticks = scaleFun(tickLabels)
    }
    else {
      if(setTickLabels)
        tickLabels = ticks
      ticks = scaleFun(ticks)
    }
    if(setTickLabels)
      tickLabels = tickLabels[is.finite(ticks)]
    ticks = ticks[is.finite(ticks)]
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=scaleFun(zlim)
    legendArgs$nlevel=length(colScale)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=colScale
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

myPairs = function(x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc, 
                   lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
                   text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3, 
                   cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1, 
                   log = "", lims=NULL) 
{
  if (doText <- missing(text.panel) || is.function(text.panel)) 
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, i, j, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if(!is.null(lims)) {
      x = lims[[i]]
      y = lims[[j]]
    }
    xpd <- NA
    if (side%%2L == 1L && xl[j]) 
      xpd <- FALSE
    if (side%%2L == 0L && yl[i]) 
      xpd <- FALSE
    if (side%%2L == 1L) 
      Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2L) 
    stop("only one column in the argument to 'pairs'")
  if (!all(horInd >= 1L && horInd <= nc)) 
    stop("invalid argument 'horInd'")
  if (!all(verInd >= 1L && verInd <= nc)) 
    stop("invalid argument 'verInd'")
  if (doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      doText <- FALSE
  }
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  main <- if ("main" %in% nmdots) 
    dots$main
  if (is.null(oma)) 
    oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
  opar <- par(mfcol = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
                                                                       4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  xl <- yl <- logical(nc)
  if (is.numeric(log)) 
    xl[log] <- yl[log] <- TRUE
  else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }
  ni <- length(iSet <- if (row1attop) horInd else rev(horInd))
  nj <- length(jSet <- verInd)
  for (j in jSet) for (i in iSet) {
    l <- paste0(if (xl[j]) 
      "x"
      else "", if (yl[i]) 
        "y"
      else "")
    if(is.null(lims)) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ..., log = l)
    } else {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", xlim=lims[[j]], ylim=lims[[i]], ..., log = l)
    }
    if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
      box()
      j.odd <- (match(j, jSet) + !row1attop)%%2L
      i.odd <- (match(i, iSet) + !row1attop)%%2L
      if (i == iSet[1L] && (!j.odd || !has.upper || !has.lower)) 
        localAxis(3L, x[, j], x[, i], j, i, ...)
      if (i == iSet[ni] && (j.odd || !has.upper || !has.lower)) 
        localAxis(1L, x[, j], x[, i], j, i, ...)
      if (j == jSet[1L] && (!i.odd || !has.upper || !has.lower)) 
        localAxis(2L, x[, j], x[, i], j, i, ...)
      if (j == jSet[nj] && (i.odd || !has.upper || !has.lower)) 
        localAxis(4L, x[, j], x[, i], j, i, ...)
      mfg <- par("mfg")
      if (i == j) {
        if (has.diag) 
          localDiagPanel(as.vector(x[, i]), ...)
        if (doText) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          xlp <- if (xl[i]) 
            10^0.5
          else 0.5
          ylp <- if (yl[j]) 
            10^label.pos
          else label.pos
          text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else if (i < j) 
        localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                       i]), ...)
      else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                          i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("the 'panel' function made a new plot")
    }
    else par(new = FALSE)
  }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
          font = font.main)
  }
  invisible(NULL)
}

# plotExampleGaussianProcess(extraPlotName="Nugget", phi=.05, sigma2=0.5^2)
# plotExampleGaussianProcess(sigma2=0)
plotExampleGaussianProcess = function(resGP=512, mu=0, marginalVariance=1^2, phi=.25, kappa=1, sigma2=.1^2, seed=1, extraPlotName="") {
  require(fields)
  require(RandomFields)
  require(spatstat)
  
  set.seed(seed)
  
  # genMaternGP generates a Matern covariance GP on the unit square
  # with the parameters given the in text
  genMaternGP = function(nsim=1, nx=resGP, ny=resGP, asList=TRUE, coords=NULL, method="circulant") {
    #mu = 4, sigma^2=1.5, phi=0.15, kappa=1, beta=2, tau^2 = 0
    
    # use RFmatern and RFsimulate
    obj = RMmatern(nu=kappa, var=marginalVariance, scale=phi)
    # obj = RMwhittle(nu=kappa, var=sigmasq, scale=phi)
    
    if(is.null(coords)) {
      coordsSet=TRUE
      xs = seq(-1, 1, length=nx)
      ys = seq(-1, 1, length=ny)
      coords = make.surface.grid(list(x=xs, y=ys))
    }
    else
      coordsSet = FALSE
    
    if(method == "instrinsic")
      sims = as.matrix(RFsimulate(RPintrinsic(obj), x=coords[,1], y=coords[,2], n=nsim)) + rnorm(nrow(coords)*nsim, sd=sqrt(sigma2)) + mu
    else if(method == "circulant")
      sims = as.matrix(RFsimulate(RPcirculant(obj), x=coords[,1], y=coords[,2], n=nsim)) + rnorm(nrow(coords)*nsim, sd=sqrt(sigma2)) + mu
    else if(method == "cutoff")
      sims = as.matrix(RFsimulate(RPcutoff(obj), x=coords[,1], y=coords[,2], n=nsim)) + rnorm(nrow(coords)*nsim, sd=sqrt(sigma2)) + mu
    
    list(coords=coords, sims=sims)
  }
    
  GP = genMaternGP()
  GPCoords = GP$coords
  # par(mfrow=c(1,1), family="serif")
  png(paste0("Figures/Illustrations/exampleGP", extraPlotName, ".png"), width=800, height=800)
  quilt.plot(GPCoords, GP$sims, nx=resGP, ny=resGP)
  # axis(1, at=seq(-1, 1, l=3))
  # axis(2, at=seq(-1, 1, l=3))
  dev.off()
}

plotExampleMaternCorrelation = function(effectiveScales = c(.1, .5, 1), sigma2=.1^2) {
  cols = rainbow(length(effectiveScales))
  
  pdf("Figures/Illustrations/matern.pdf", width=5, height=5)
  xs = seq(0, 1, l=200)
  plot(xs, (1/(1 + sqrt(sigma2))) * Matern(xs, effectiveScales[1] / sqrt(8), smoothness = 1), type="l", col=cols[1], 
       main="Matern correlation functions", xlab="Distance", ylab="Correlation", ylim=c(0,1))
  if(length(effectiveScales) >= 1) {
    for(i in 2:length(effectiveScales)) {
      lines(xs, (1/(1 + sqrt(sigma2))) * Matern(xs, effectiveScales[i] / sqrt(8), smoothness = 1), col=cols[i])
    }
    legend("topright", paste0("Effective range=", effectiveScales), lty=1, col=cols)
  }
  points(0, 1, pch=19, cex=.7)
  dev.off()
  
  # pdf("Figures/Illustrations/maternINLA.pdf", width=5, height=5)
  # xs = seq(0, 1, l=200)
  # plot(xs, (1 - sqrt(sigma2)) * inla.matern.cov(x=xs, kappa=1 / (effectiveScales[1] / sqrt(8)), nu = 1, corr=TRUE, d=2), type="l", col=cols[1], 
  #      main="Matern correlation functions", xlab="Distance", ylab="Correlation", ylim=c(0,1))
  # if(length(effectiveScales) >= 1) {
  #   for(i in 2:length(effectiveScales)) {
  #     lines(xs, (1 - sqrt(sigma2)) * inla.matern.cov(x=xs, kappa=1 / (effectiveScales[i] / sqrt(8)), nu = 1, corr=TRUE, d=2), col=cols[i])
  #   }
  #   legend("topright", paste0("Effective range=", effectiveScales), lty=1, col=cols)
  # }
  # points(0, 1, pch=19, cex=.7)
  # dev.off()
}

plotBasisKnots = function(NC=14, nBuffer=5, main="Basis Knot Locations") {
  # get lattice points, prediction points
  nLayer=3
  xRangeDat = c(-1, 1)
  yRangeDat = c(-1, 1)
  latInfo = makeLatGrids(xRangeDat, yRangeDat, NC, nBuffer, nLayer=nLayer)
  
  # get first, coarsest layer
  pdf("Figures/illustrations/basisKnots.pdf", width=5, height=5)
  pts = latInfo[[1]]$latCoords
  border = rbind(c(xRangeDat[1], yRangeDat[1]), 
                 c(xRangeDat[1], yRangeDat[2]), 
                 c(xRangeDat[2], yRangeDat[2]), 
                 c(xRangeDat[2], yRangeDat[1]))
  plot(pts[,1], pts[,2], type="n", xlab="x", ylab="y", main=main)
  polygon(border[,1], border[,2], col=rgb(.5, .5, .5, .5), border=rgb(0, 0, 0, 0))
  points(pts[,1], pts[,2], pch="+", cex=1)
  
  # plot second layer
  pts = latInfo[[2]]$latCoords
  points(pts[,1], pts[,2], pch="+", cex=.6, col="darkGreen")
  
  # plot third layer
  pts = latInfo[[3]]$latCoords
  points(pts[,1], pts[,2], pch="+", cex=.25, col="red")
  dev.off()
}

plotDiskDistanceDistribution = function() {
  r=1
  dDiskDist = function(d, r = 1) {
    out = 4 * d / (pi*r^2) * (acos(d / (2*r)) - d / (2*r) * sqrt(1 - (d / (2*r))^2))
    out[d<0] = 0
    out[d>2*r] = 0
    out
  }
  
  ds = seq(0, 2 * r, l=200)
  pdf("Figures/Illustrations/diskDistanceDistribution.pdf", width=5, height=5)
  plot(ds, dDiskDist(ds, r), type="l", main="Distances between points in disk (radius R)", ylab="Probability density", xlab="Distance", xaxt="n", yaxt="n")
  xtick<-seq(0, 2, by=1)
  axis(side=1, at=xtick, labels = FALSE)
  text(x=xtick,  par("usr")[3]-.02, 
       labels = c("0", "R", "2R"), pos = 1, xpd = TRUE)
  
  ytick<-seq(0, max(dDiskDist(ds, r)), l=5)
  axis(side=2, at=ytick, labels = FALSE)
  text(par("usr")[1]-.01, ytick, 
       labels = c("0", TeX("$\\frac{1}{5R}$"), TeX("$\\frac{2}{5R}$"), TeX("$\\frac{3}{5R}$"), TeX("$\\frac{4}{5R}$")), pos = 2, xpd = TRUE)
  dev.off()
  
  print(ds[which.max(dDiskDist(ds, r))])
}

generatePlottingSymbols = function() {
  urbCols = makeGreenBlueSequentialColors(29)
  pdf('Figures/Illustrations/greenTriangle.pdf', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=17, cex=1, col=urbCols[1])
  dev.off()
  
  pdf('Figures/Illustrations/blueSquare.pdf', width=.2, height=.2)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1), xaxs="i", yaxs="i")
  points(0, -0.5, pch=15, cex=1, col=urbCols[29])
  dev.off()
}

getCustomScaleTicks = function(usr, scaleFun=sqrt, nint=5, log=FALSE) {
  # axisTicks(range(scaleFun(x)), log=log, nint=nint)
  axp <- unlist(.axisPars(range(scaleFun(usr)), log = log, nintLog = nint), 
                use.names = FALSE)
  .Call(C_R_CreateAtVector, axp, usr, nint, log)
}


getCustomScaleTicks = function(usr, scaleFun=sqrt, nint=5, log=FALSE) {
  # axisTicks(range(scaleFun(x)), log=log, nint=nint)
  axp <- unlist(.axisPars(range(scaleFun(usr)), log = log, nintLog = nint), 
                use.names = FALSE)
  .Call(C_R_CreateAtVector, axp, usr, nint, log)
}


#' Determine square root axis tick mark positions
#'
#' Determine square root axis tick mark positions
#'
#' This function calculates positions for tick marks for data
#' that has been transformed with `sqrt()`, specifically a directional
#' transformation like `sqrt(abs(x)) * sign(x)`.
#'
#' The main goal of this function is to provide reasonably placed
#' tick marks using integer values.
#'
#' @return
#' Invisibly returns a numeric vector of axis tick positions,
#' named by the display label.
#' The axis values are in square root space while the labels represent
#' the normal space values.
#'
#' @family jam plot functions
#'
#' @param side integer value indicating the axis position, as used
#'    by `axis()`, 1=bottom, 2=left, 3=top, 4=right.
#' @param x optional numeric vector representing the numeric range
#'    to be labeled.
#' @param pretty.n numeric value indicating the number of desired
#'    tick marks, passed to `pretty()`.
#' @param u5.bias numeric value passed to `pretty()` to influence the
#'    frequency of intermediate tick marks.
#' @param big.mark character value passed to `format()` which helps
#'    visually distinguish numbers larger than 1000.
#' @param plot logical indicating whether to plot the axis tick
#'    marks and labels.
#' @param las,cex.axis numeric values passed to `axis()` when drawing
#'    the axis, by default `las=2` plots labels rotated
#'    perpendicular to the axis.
#' @param ... additional parameters are passed to `pretty()`.
#'
#' @export
sqrtAxis <- function
(side=1,
 x=NULL,
 pretty.n=10,
 u5.bias=0,
 big.mark=",",
 plot=TRUE,
 las=2,
 cex.axis=0.6,
 ...)
{
  ## Purpose is to generate a set of tick marks for sqrt
  ## transformed data axes.  It assumes data is already sqrt-transformed,
  ## and that negative values have been treated like:
  ## sqrt(abs(x))*sign(x)
  if (length(side) > 2) {
    x <- side;
    side <- 0;
  }
  if (length(side) == 0) {
    side <- 0;
  }
  if (1 %in% side) {
    xRange <- par("usr")[1:2];
  } else if (2 %in% side) {
    xRange <- par("usr")[3:4];
  } else if (length(x) > 0) {
    xRange <- range(x, na.rm=TRUE);
  }
  
  subdivideSqrt <- function(atPretty1, n=pretty.n, ...) {
    ## Purpose is to take x in form of 0,x1,
    ## and subdivide using pretty()
    atPretty1a <- unique(sort(abs(atPretty1)));
    atPretty1b <- tail(atPretty1a, -2);
    atPretty2a <- pretty(head(atPretty1a,2), n=n, ...);
    return(unique(sort(c(atPretty2a, atPretty1b))));
  }
  
  ## Determine tick positions
  nSubFactor <- 2.44;
  
  atPretty1 <- pretty(xRange^2*sign(xRange),
                      u5.bias=u5.bias,
                      n=(pretty.n)^(1/nSubFactor));
  
  atPretty1old <- atPretty1;
  while (length(atPretty1) <= pretty.n) {
    atPretty1new <- subdivideSqrt(atPretty1,
                                  n=noiseFloor(minimum=2, (pretty.n)^(1/nSubFactor)));
    atPretty1 <- atPretty1new[atPretty1new <= max(abs(xRange^2))];
    atPretty1old <- atPretty1;
  }
  atPretty3 <- unique(sort(
    rep(atPretty1,
        each=length(unique(sign(xRange)))) * sign(xRange)));
  atPretty <- atPretty3[
    (atPretty3 >= head(xRange,1)^2*sign(head(xRange,1)) &
       atPretty3 <= tail(xRange, 1)^2*sign(tail(xRange, 1)))];
  
  xLabel <- sapply(atPretty, function(i){
    format(i,
           trim=TRUE,
           digits=2,
           big.mark=big.mark);
  });
  ## Transform to square root space
  atSqrt <- sqrt(abs(atPretty))*sign(atPretty);
  if (plot) {
    axis(side=side,
         at=atSqrt,
         labels=xLabel,
         las=las,
         cex.axis=cex.axis,
         ...);
  }
  
  invisible(nameVector(atSqrt, xLabel));
}

#' Apply noise floor and ceiling to numeric vector
#'
#' Apply noise floor and ceiling to numeric vector
#'
#' A noise floor is useful when detected numeric values are sometimes below
#' a clear noise threshold, and where some downstream ratio may be calculated
#' using these values. Applying a noise floor ensures the ratios and not
#' artificially higher, especially in cases where the values involved are
#' least reliable. This procedure is expected to produce more conservative
#' and appropriate ratios in that scenario.
#'
#' A ceiling is similar, values above the ceiling are set to the ceiling,
#' which is practical when values above a certain threshold are conceptually
#' similar to those at the threshold. One clear example is plotting
#' `-log10(Pvalue)` when the range of P-values might approach 1e-1000.
#' In this case, setting a ceiling of 50 conceptually equates P-values
#' below 1e-50, while also restricting the axis range of a plot.
#'
#' The ability to set values at the floor to a different value, using
#' `newValue` different from `minimum`, is intended to allow separation
#' of numeric values from the floor for illustrative purposes.
#'
#' @return
#' A numeric vector or matrix, matching the input type `x` where numeric
#' values are fixed to the `minimum` and `ceiling` values as defined
#' by `newValue` and `newCeiling`, respectively.
#'
#' @family jam numeric functions
#'
#' @param x numeric vector or matrix
#' @param minimum numeric floor value
#' @param newValue numeric, by default the same as the floor value. Sometimes
#'    it can be useful to define a different value, one example is to define
#'    values as NA, or another distinct number away from the floor.
#' @param adjustNA logical whether to change NA values to the newValue.
#' @param ceiling numeric value, optionally a ceiling. If defined, then values
#'    above the ceiling value are set to newCeiling.
#' @param newCeiling numeric value when ceiling is defined, values above the
#'    ceiling are set to this numeric value.
#' @param ... additional parameters are ignored.
#'
#' @examples
#' # start with some random data
#' n <- 2000;
#' x1 <- rnorm(n);
#' y1 <- rnorm(n);
#'
#' # apply noise floor and ceiling
#' x2 <- noiseFloor(x1, minimum=-2, ceiling=2);
#' y2 <- noiseFloor(y1, minimum=-2, ceiling=2);
#'
#' # apply noise floor and ceiling with custom replacement values
#' xm <- cbind(x=x1, y=y1);
#' xm3 <- noiseFloor(xm,
#'    minimum=-2, newValue=-3,
#'    ceiling=2, newCeiling=3);
#'
#' parMfrow <- par("mfrow");
#' par("mfrow"=c(2,2));
#' plotSmoothScatter(x1, y1);
#' plotSmoothScatter(x2, y2);
#' plotSmoothScatter(xm3);
#' par("mfrow"=parMfrow);
#'
#' @export
noiseFloor <- function
(x,
 minimum=0,
 newValue=minimum,
 adjustNA=FALSE,
 ceiling=NULL,
 newCeiling=ceiling,
 ...)
{
  ## Purpose is to apply a noise floor, that is, to set all values
  ## to be at least 'minimum' amount.
  ## This function performs no scaling or normalization.
  if (length(x) == 0) {
    return(x);
  }
  if (length(minimum) > 0) {
    if (adjustNA) {
      x[is.na(x) | (!is.na(x) & x < minimum)] <- newValue;
    } else {
      x[!is.na(x) & x < minimum] <- newValue;
    }
  }
  if (length(ceiling) > 0) {
    x[!is.na(x) & x > ceiling] <- newCeiling;
  }
  return(x);
}

#' assign unique names for a vector
#'
#' assign unique names for a vector
#'
#' This function assigns unique names to a vector, if necessary it runs
#' \code{\link{makeNames}} to create unique names. It differs from
#' \code{\link[stats]{setNames}} in that it ensures names are unique,
#' and when no names are supplied, it uses the vector itself to define
#' names. It is helpful to run this function inside an \code{\link[base]{lapply}}
#' function call, which by default maintains names, but does not assign
#' names if the input data did not already have them.
#'
#' When used with a data.frame, it is particularly convenient to pull out
#' a named vector of values. For example, log2 fold changes by gene, where
#' the gene symbols are the name of the vector.
#'
#' \code{nameVector(genedata[,c("Gene","log2FC")])}
#'
#' @return vector with names defined
#'
#' @family jam string functions
#'
#' @param x vector input, or data.frame, matrix, or tibble with two columns,
#'    the second column is used to name values in the first column.
#' @param y NULL or character vector of names. If NULL then x is used.
#'    Note that y is recycled to the length of x, prior to being sent
#'    to the makeNamesFunc.
#'    In fringe cases, y can be a matrix, data.frame, or tibble, in which
#'    case \code{\link{pasteByRow}} will be used to create a character string
#'    to be used for vector names. Note this case is activated only when x
#'    is not a two column matrix, data.frame, or tibble.
#' @param makeNamesFunc function to make names unique, by default
#'    \code{\link{makeNames}} which ensures names are unique.
#' @param ... passed to \code{\link{makeNamesFunc}}, or to
#'    \code{\link{pasteByRow}} if y is a two column data.frame, matrix, or
#'    tibble. Thus, \code{sep} can be defined here as a delimiter between
#'    column values.
#'
#' @examples
#' # it generally just creates names from the vector values
#' nameVector(LETTERS[1:5]);
#'
#' # if values are replicated, the makeNames() function makes them unique
#' V <- rep(LETTERS[1:5], each=3);
#' nameVector(V);
#'
#' # for a two-column data.frame, it creates a named vector using
#' # the values in the first column, and names in the second column.
#' df <- data.frame(seq_along(V), V);
#' df;
#' nameVector(df);
#'
#' # Lastly, admittedly a fringe case, it can take a multi-column data.frame
#' # to generate labels:
#' nameVector(V, df);
#'
#' @export
nameVector <- function
(x,
 y=NULL,
 makeNamesFunc=makeNames,
 ...)
{
  ## Purpose is to name a vector with its own values,
  ## useful for lapply which only names output if the input
  ## vector has names.
  ##
  ## A neat trick is to use the _v# naming scheme in makeNames to
  ## create unique names based upon a single label, e.g.
  ## set1colors <- nameVector(brewer.pal(15, "Set1"), "Set1");
  ##   Set1_v1   Set1_v2   Set1_v3   Set1_v4   Set1_v5   Set1_v6   Set1_v7   Set1_v8   Set1_v9
  ## "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"
  ##
  ## Added bonus, if given a 2-column table, it'll use them as x and y
  if (igrepHas("dataframe", class(x))) {
    x <- as.data.frame(x);
  }
  if (igrepHas("data.frame", class(x)) && ncol(x) == 2) {
    y <- x[[2]];
    x <- x[[1]];
  } else if (igrepHas("matrix", class(x)) && ncol(x) == 2) {
    y <- x[,2];
    x <- x[,1];
  }
  if (length(y) > 0) {
    if (igrepHas("data.frame|matrix", class(y))) {
      ## If given a data.frame use pasteByRow() to create a string
      y <- pasteByRow(y, ...);
    }
    names(x) <- makeNamesFunc(rep(y, length.out=length(x)), ...);
  } else {
    names(x) <- makeNamesFunc(x, ...);
  }
  return(x);
}

#' vector contains any case-insensitive grep match
#'
#' vector contains any case-insensitive grep match
#'
#' This function checks the input vector for any elements matching the
#' grep pattern. The grep is performed case-insensitive (igrep). This function
#' is particularly useful when checking function arguments or object class,
#' where the class(a) might return multiple values, or where the name of
#' the class might be slightly different than expected, e.g. data.frame,
#' data_frame, DataFrame.
#'
#' @param pattern the grep pattern to use with `base::grep()`
#' @param x vector to use in the grep
#' @param ignore.case logical default TRUE, meaning the grep will be performed
#'    in case-insensitive mode.
#' @param minCount integer minimum number of matches required to return TRUE.
#' @param naToBlank logical whether to convert NA to blank, instead of
#'    allowing grep to handle NA values as-is.
#'
#' @return logical indicating whether the grep match criteria were met,
#'    TRUE indicates the grep pattern was present in minCount or more
#'    number of entries.
#'
#' @seealso `base::grep()`
#'
#' @examples
#' a <- c("data.frame","data_frame","tibble","tbl");
#' igrepHas("Data.*Frame", a);
#' igrepHas("matrix", a);
#'
#' @family jam grep functions
#'
#' @export
igrepHas <- function
(pattern, x=NULL, ignore.case=TRUE,
 minCount=1, naToBlank=FALSE,
 ...)
{
  ## Purpose is a quick check for greppable substring, for if() statements
  ##
  ## naToBlank=TRUE will convert NA values to "" prior to running grep
  ##
  ## The special case where minCount is negative (minCount == -1) or larger
  ## than length(x), it will be set to length(x) and therefore
  ## requires all elements of x to meet the grep criteria
  if (minCount < 0 || minCount > length(x)) {
    minCount <- length(x);
  }
  if (length(x) == 0) {
    return(FALSE);
  } else {
    if (naToBlank && any(is.na(x))) {
      x[is.na(x)] <- "";
    }
    length(grep(pattern=pattern,
                x=x,
                ignore.case=ignore.case,
                ...)) >= as.integer(minCount);
  }
}

#' make unique vector names
#'
#' make unique vector names
#'
#' This function extends the basic goal from \code{\link[base]{make.names}}
#' which is intended to make syntactically valid names from a character vector.
#' This makeNames function makes names unique, and offers configurable methods
#' to handle duplicate names. By default, any duplicated entries receive a
#' suffix _v# where # is s running count of entries observed, starting at 1.
#' The \code{\link[base]{make.names}} function, by contrast, renames the
#' second observed entry starting at .1, leaving the original entry
#' unchanged. Optionally, makeNames can rename all entries with a numeric
#' suffix, for consistency.
#'
#' For example:
#' \code{A, A, A, B, B, C}
#' becomes:
#' \code{A_v1, A_v2, A_v3, B_v1, B_v2, C}
#'
#' Also, makeNames always allows "_".
#'
#' This makeNames function is similar to \code{\link[base]{make.unique}}
#' which also converts a vector into a unique vector by adding suffix values,
#' however the \code{\link[base]{make.unique}} function intends to allow
#' repeated operations which recognize duplicated entries and continually
#' increment the suffix number. This makeNames function currently does not
#' handle repeat operations. The recommended approach to workaround having
#' pre-existing versioned names would be to remove suffix values prior to
#' running this function. One small distinction from
#' \code{\link[base]{make.unique}} is that makeNames does version the first
#' entry in a set.
#'
#' @return character vector of unique names
#'
#' @family jam string functions
#'
#' @param x character vector to be used when defining names. All other
#'    vector types will be coerced to character prior to use.
#' @param unique argument which is ignored, included only for
#'    compatibility with `base::make.names`. All results from
#'    `makeNames()` are unique.
#' @param suffix character separator between the original entry and the
#'    version, if necessary.
#' @param renameOnes logical whether to rename single, unduplicated, entries.
#' @param doPadInteger logical whether to pad integer values to a consistent
#'    number of digits, based upon all suffix values needed. This output
#'    allows for more consistent sorting of names. To define a fixed number
#'    of digits, use the useNchar parameter.
#' @param useNchar integer or NULL, number of digits to use when padding
#'    integer values with leading zero, only relevant when usePadInteger=TRUE.
#' @param startN integer number used when numberStyle is "number", this integer
#'    is used for the first entry to be renamed. You can use this value to
#'    make zero-based suffix values, for example.
#' @param numberStyle character style for version numbering
#'    \describe{
#'       \item{"number"}{Use integer numbers to represent each duplicated
#'          entry.}
#'       \item{"letters"}{Use lowercase letters to represent each duplicated
#'          entry. The 27th entry uses the pattern "aa" to represent two
#'          26-base digits. When doPadInteger=TRUE, a zero is still used
#'          to pad the resulting version numbers, again to allow easy sorting
#'          of text values, but also because there is no letter equivalent
#'          for the number zero.
#'          It is usually best to change the suffix to "_" or "" when using
#'          "letters".}
#'       \item{"LETTERS"}{Use uppercase letters to represent each duplicated
#'          entry, with the same rules as applied to "letters".}
#'    }
#' @param renameFirst logical whether to rename the first entry in a set of
#'    duplicated entries. If FALSE then the first entry in a set will not
#'    be versioned, even when renameOnes=TRUE.
#' @param keepNA logical whether to retain NA values using the string "NA".
#'    If keepNA is FALSE, then NA values will remain NA, thus causing some
#'    names to become `<NA>`, which can cause problems with some downstream
#'    functions which assume all names are either NULL or non-NA.
#'
#' @examples
#' V <- rep(LETTERS[1:3], c(2,3,1));
#' makeNames(V);
#' makeNames(V, renameOnes=TRUE);
#' makeNames(V, renameFirst=FALSE);
#' exons <- makeNames(rep("exon", 3), suffix="");
#' makeNames(rep(exons, c(2,3,1)), numberStyle="letters", suffix="");
#'
#' @export
makeNames <- function
(x,
 unique=TRUE,
 suffix="_v",
 renameOnes=FALSE,
 doPadInteger=FALSE,
 startN=1,
 numberStyle=c("number","letters","LETTERS"),
 useNchar=NULL,
 renameFirst=TRUE,
 keepNA=TRUE,
 ...)
{
  ## Purpose is to make unique names without the R mangling that comes
  ## with make.names().
  ## By default, unique entries are not renamed, and entries with two or
  ## more replicates are renamed to NAME_v1, NAME_v2, NAME_v3, etc.
  ##
  ## if renameOnes=TRUE, it will rename singlets to NAME_v1 even if there
  ## is only one entry.
  ##
  ## renameFirst=TRUE will rename each duplicated entry NAME_v1, NAME_v2,
  ## NAME_v3, etc.
  ## renameFirst=FALSE will not rename the first in a set of duplicated
  ## entries, e.g. NAME, NAME_v1, NAME_v2, etc.
  ##
  ## The distinction between renameOnes and renameFirst:
  ## renameOnes=TRUE will rename all singlets and duplicated entries,
  ## starting with the first entry.
  ## renameOnes=FALSE will not rename singlet entries.
  ## renameFirst=TRUE will only rename duplicated entries, starting with
  ## the first entry.
  ## renameFirst=FALSE will not rename the first entry in a set of
  ## duplicated entries.
  ##
  ## the suffix can be changed, e.g. "_r" will name names NAME_r1,
  ## NAME_r2, NAME_r3, etc.
  ##
  ## numberStyle="number" uses integers as the suffix
  ## numberStyle="letters" uses lowercase letters as digits, similar to Excel column names
  ## numberStyle="LETTERS" uses uppercase letters as digits, similar to Excel column names
  ## Be aware that letters can only go to roughly 18,000 entries, given the current implementation
  ## of colNum2excelName
  ##
  ## When useNchar is numeric, it sets doPadInteger=TRUE, and will use at least
  ## that many digits in padding the integer.
  ##
  ##
  ## TODO:
  ## Update logic to be analogous to using make.unique(), which intends
  ## to maintain previous versioning of names without appending deeper
  ## suffices as appropriate.
  ## E.g.     c("",    "",    "",    "_v1", "_v2", "_v3")
  ## becomes  c("_v4", "_v5", "_v6", "_v1", "_v2", "_v3")
  ## instead of
  ##          c("_v1", "_v2", "_v3", "_v1", "_v2", "_v3")
  ## or
  ##          c("_v1_v1", "_v2_v1", "_v3_v1", "_v1_v2", "_v2_v2", "_v3_v2")
  ##
  if (length(x) == 0) {
    return(x);
  }
  numberStyle <- match.arg(numberStyle);
  if (!is.null(useNchar)) {
    useNchar <- as.integer(useNchar);
    doPadInteger=TRUE;
  }
  if (any(c("factor", "ordered") %in% class(x))) {
    x <- as.character(x);
  }
  if (keepNA && any(is.na(x))) {
    x <- rmNA(x,
              naValue="NA");
  }
  
  ## First check for duplicates using anyDuplicated()
  ## version 0.0.35.900, this change speeds assignment
  ## in large vectors when most entries are not duplicated.
  dupes <- duplicated(x);
  
  ## Convert entries to a named count of occurences of each entry
  if (any(dupes)) {
    xSubDupes <- table(x[dupes]) + 1;
    maxCt <- max(c(1,xSubDupes));
    xSubOnes <- setNames(rep(1, sum(!dupes)), x[!dupes]);
    xSub <- c(xSubOnes, xSubDupes);
  } else {
    xSub <- setNames(rep(1, sum(!dupes)), x[!dupes]);
    maxCt <- 1;
  }
  ## version 0.0.34.900 and previous used the method below
  #xSub <- table(as.character(x));
  
  ## Vector of counts to be used
  versionsV <- as.integer(renameFirst):maxCt + startN - 1;
  
  ## If using letters, define the set of letter upfront to save processing
  if (igrepHas("letters", numberStyle)) {
    if (numberStyle %in% "letters") {
      useLetters <- letters[1:26];
      zeroVal <- "A";
    } else {
      useLetters <- LETTERS[1:26];
      zeroVal <- "a";
    }
    num2letters <- colNum2excelName(versionsV,
                                    useLetters=useLetters,
                                    zeroVal=zeroVal,
                                    ...);
    versionsV <- num2letters;
  }
  if (doPadInteger) {
    versionsV <- padInteger(versionsV,
                            useNchar=useNchar,
                            ...);
  }
  
  ## If no duplicated entries
  if (max(xSub) %in% c(-Inf,1)) {
    ## If not renaming the singlet entries, send the same list back
    if (!renameOnes) {
      return(x);
    } else {
      ## If renaming singlets, simply paste the suffix and first entry
      return(paste0(x, suffix, head(versionsV, 1)));
    }
  }
  
  if (renameOnes) {
    xUse <- 1:length(x);
  } else {
    xUse <- (x %in% names(xSub)[xSub > 1]);
  }
  xSub1 <- x[xUse];
  
  ## Preserve the original order
  names(xSub1) <- padInteger(seq_along(xSub1));
  ## Split the vector into a list of vectors, by name
  xSub2 <- split(xSub1, xSub1);
  names(xSub2) <- NULL;
  
  ## Optionally pad the integer to facilitate sorting
  ## Note: This padding only pads integers within each name,
  ## not across all names.
  xSub3 <- lapply(xSub2, function(i){
    versionsV[seq_along(i)];
  });
  
  ## Now simply paste the value, the suffix, and the new version
  xSub1v <- paste0(unlist(xSub2), suffix, unlist(xSub3));
  
  ## Re-order the vector using the original order
  names(xSub1v) <- names(unlist(xSub2));
  xSub1v2 <- xSub1v[names(xSub1)];
  
  ## Assign only the entries we versioned
  x[xUse] <- xSub1v2;
  
  ## Last check for renameFirst=FALSE, in which case we remove the first
  ## versioned entry
  if (!renameFirst) {
    escapeRegex <- function(string){
      gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", string);
    }
    firstVer <- paste0(escapeRegex(paste0(suffix, head(versionsV, 1))), "$");
    x[xUse] <- gsub(firstVer, "", x[xUse]);
  }
  
  return(x);
}






