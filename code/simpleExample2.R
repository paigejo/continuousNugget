# In this script, we look at a simple example of a real area, attempting 
# to estimate NMR accounting for EA level random variation
# NOTE: same as simpleExample.R, but uses recently written code added to SUMMER

# download data ----
# download 2014 Kenya population density and associated TIF file
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/pop.rda?raw=true")
popFilename = paste0(tempDirectory, "/pop.rda")
if(!file.exists(popFilename)) {
  download.file(githubURL, popFilename)
}

githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/worldpop_total_1y_2014_00_00.tif?raw=true")
popTIFFilename = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")
if(!file.exists(popTIFFilename)) {
  download.file(githubURL,popTIFFilename)
}

# load it in
require(raster)
out = load(popFilename)
out

# make sure this is correct for re-projections
pop@file@name = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")

## setup ----
longRangeWajir = c(39, 41)
latRangeWajir = c(0.25, 3.6)
longRangeNorthWajir = c(39, 40.3)
latRangeNorthWajir = c(2.45, 3.7)
constituenciesW = poppsubKenya$subarea[poppsubKenya$area=="Wajir"]
constituencies = poppsubKenya$subarea
offsets = matrix(0, nrow=6, ncol=2)
offsets[1,2] = .1 # shift label for Eldas slightly higher
offsets[6,1] = .15 # shift label for Wajir West slightly farther east
easpsub = meanEAsPerCon2()
easpsub = easpsub[easpsub$area=="Wajir",]
popGridWajir = popGrid[popGrid$area=="Wajir",]
popGridWajirNorth = popGrid[popGrid$subarea=="Wajir North",]
popGridWajirAdjusted = popGridAdjusted[popGridAdjusted$area=="Wajir",]
poppaSimpleAdjusted = poppaKenya[order(poppaKenya$area),]
poppaSimpleAdjusted$popUrb = easpa$EAUrb * 25
poppaSimpleAdjusted$popRur = easpa$EARur * 25
poppaSimpleAdjusted$popTotal = easpa$EATotal * 25
poppaSimpleAdjustedW = poppaSimpleAdjusted[poppaSimpleAdjusted$area == "Wajir",]
popMatSimpleAdjusted = adjustPopGrid(popGrid, poppaSimpleAdjusted)
popMatSimpleAdjusted$subarea = popMatSimpleAdjusted$subarea
popMatSimpleAdjusted$area = popMatSimpleAdjusted$area
popMatSimpleAdjustedNW = popMatSimpleAdjusted[popGrid$subarea=="Wajir North",]
# plotMapDat(mapDat=adm0, lwd=.5, new=TRUE)
# plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], border="green")

# restrict simulation to N Wajir
popMatSimpleAdjustedW = popMatSimpleAdjusted[popMatSimpleAdjusted$area == "Wajir", ]
easpaW = makeDefaultEASPA()
easpaW = easpaW[easpaW$area == "Wajir", ]
popMatW = makeDefaultPopMat()
popMatW = popMatW[popMatW$area == "Wajir", ]
poppsubW = poppsubKenya
poppsubW = poppsubW[poppsubW$area == "Wajir", ]

# set up color scales
redBlueScale = makeRedBlueDivergingColors(64, rev=TRUE)
purpleRedScale = makePurpleRedSequentialColors(64)
yellowScale = makeYellowSequentialColors(64)
yellowRedScale = rev(makeRedYellowSequentialColors(64))
greenScale = makeGreenSequentialColors(64)
summerScale = rev(makeBlueGreenYellowSequentialColors(64))
blueYellowRedScale = rev(makeRedYellowBlueColors(64))
purpleYellowScale = makePurpleYellowSequentialColors(64)
popCols=makeBlueSequentialColors(64)
ruralCols=makeGreenSequentialColors(64)
urbanCols = popCols
riskCols = purpleRedScale
divergingRiskCols = redBlueScale
sdCols = purpleYellowScale

# plot Wajir data ----
pdf(file="figures/simpleExample/wajirData.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
# ticks = c(1, 10, 100, 1000)
# tickLabels = as.character(ticks)
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(ticks), labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
quilt.plot(mort$lon[mort$admin1=="Wajir"], mort$lat[mort$admin1=="Wajir"], mort$y[mort$admin1=="Wajir"]/mort$n[mort$admin1=="Wajir"], 
           nx=50, ny=60, add=TRUE)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# make 1km res pop density grid ----
if(FALSE) {
  popGridFine = makePopIntegrationTab(kmRes=1, pop=pop, domainPoly=kenyaPoly, 
                                     eastLim=eastLim, northLim=northLim, 
                                     mapProjection=projKenya, poppa=poppaKenya, 
                                     poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                     areaMapDat=adm1, subareaMapDat=adm2, 
                                     areaPolygonSubsetI=46)
  
  popGridFineAdjusted = adjustPopMat(popGridFine, poppaTarget=poppaSimpleAdjustedW, adjustBy="area")
  
  # popGridFine = makeInterpPopGrid(kmRes=1, mean.neighbor=500, delta=.05, poppcon=poppcon)
  # popGridFineAdjusted = adjustPopGrid(popGridFine, poppaSimpleAdjusted)
  
  popGridFine = popGridFine[popGridFine$area == "Wajir",]
  popGridFineAdjusted = popGridFineAdjusted[popGridFineAdjusted$area == "Wajir",]
  
  # normalize to have the correct population within the county
  # popGridFine$pop = popGridFine$pop * (poppc$popTotal[poppsubKenya$area=="Wajir"] / sum(popGridFine$pop))
  popRange = range(popGridFine$pop)
  popRangeTarget = range(popGridFineAdjusted$pop)
  
  ## done with setup. Make plots
  
  # plot general population density on fine grid
  pdf(file="figures/simpleExample/wajirPopDensity.pdf", width=4, height=5)
  par(mar=c(4.1, 4.1, 1.1, 4.5))
  plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE,
             kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir,
             leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
  popTicks = c(1, 10, 100, 1000)
  popTickLabels = as.character(popTicks)
  image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1),
             legend.cex=.5, legend.width=1)
  quilt.plot(popGridFine$lon, popGridFine$lat, log(popGridFine$pop),
             col=popCols, nx=230, ny=380, add.legend = FALSE, add=TRUE,
             zlim=log(popRange))
  plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
  addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
  dev.off()
  
  # plot target population density on fine grid
  pdf(file="figures/simpleExample/wajirTargetPopDensity.pdf", width=4, height=5)
  par(mar=c(4.1, 4.1, 1.1, 4.5))
  plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE,
             kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir,
             leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
  popTicks = c(1, 10, 100, 1000)
  popTickLabels = as.character(popTicks)
  image.plot(zlim=log(popRangeTarget), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1),
             legend.cex=.5, legend.width=1)
  quilt.plot(popGridFineAdjusted$lon, popGridFineAdjusted$lat, log(popGridFineAdjusted$pop),
             col=popCols, nx=230, ny=380, add.legend = FALSE, add=TRUE,
             zlim=log(popRangeTarget))
  plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
  addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
  dev.off()
}

# Plot urbanicity ----
# now plot the urban pixels
png(file="figures/simpleExample/wajirUrban.png", width=400, height=500)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=FALSE, legend.mar=0, addColorBar=FALSE)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
popInWajir = popGrid$area == "Wajir"
quilt.plot(popGrid$lon[popInWajir], popGrid$lat[popInWajir], popGrid$urban[popInWajir], col=c(rgb(0, 0, 0, 0), "blue"), nx=60, ny=100, add.legend = FALSE, add=TRUE)
offsets = matrix(0, nrow=6, ncol=2)
offsets[1,2] = .1 # shift label for Eldas slightly higher
offsets[6,1] = .15 # shift label for Wajir West slightly farther east
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.7)
dev.off()

# pdf version:
pdf(file="figures/simpleExample/wajirUrban.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=FALSE, legend.mar=0, addColorBar=FALSE)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
popInWajir = popGrid$area == "Wajir"
quilt.plot(popGrid$lon[popInWajir], popGrid$lat[popInWajir], popGrid$urban[popInWajir], col=c(rgb(0, 0, 0, 0), "blue"), nx=60, ny=100, add.legend = FALSE, add=TRUE)
offsets = matrix(0, nrow=6, ncol=2)
offsets[1,2] = .1 # shift label for Eldas slightly higher
offsets[6,1] = .15 # shift label for Wajir West slightly farther east
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.7)
dev.off()

# now plot urban fraction as a function of constituency
pdf(file="figures/simpleExample/wajirUrbanFrac.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
urbanFraction = (poppsubW$popUrb/poppsubW$popTotal)
ruralUrbanCols=makeGreenBlueSequentialColors(32)
# theseCols = getColorsFromScale(urbanFraction, c(.02, .7), ruralUrbanCols, forceValuesInRange=TRUE)
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(urbanFraction), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=ruralUrbanCols, zlim=log(c(.0045, .7)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(.02, .07, .2, .7)), tickLabels=c(.02, .07, .2, .7))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()


# now plot urban fraction as a function of constituency
pdf(file="figures/simpleExample/wajirUrbanFrac.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
urbanFraction = (poppsubKenya$popUrb/poppsubKenya$popTotal)[poppsubKenya$area=="Wajir"]
ruralUrbanCols=makeGreenBlueSequentialColors(32)
# theseCols = getColorsFromScale(urbanFraction, c(.02, .7), ruralUrbanCols, forceValuesInRange=TRUE)
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(urbanFraction), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=ruralUrbanCols, zlim=log(c(.0045, .7)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(.02, .07, .2, .7)), tickLabels=c(.02, .07, .2, .7))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot expected EAs ----
# plot expected number of urban EAs per constituency

pdf(file="figures/simpleExample/wajirUrbanEAs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
meanUrbanEAs = easpsub$meanUrbanEAs
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(meanUrbanEAs), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=urbanCols, zlim=log(c(2, 200)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(2, 5, 50, 150)), tickLabels=c(2, 5, 50, 150))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot expected number of rural EAs per constituency

pdf(file="figures/simpleExample/wajirRuralEAs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
meanRuralEAs = easpsub$meanRuralEAs
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(meanRuralEAs), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=ruralCols, zlim=log(c(2, 200)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(2, 5, 50, 150)), tickLabels=c(2, 5, 50, 150))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot expected number of rural EAs per constituency

pdf(file="figures/simpleExample/wajirRuralEAsSelfScaled.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
meanRuralEAs = easpsub$meanRuralEAs
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(meanRuralEAs), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=ruralCols, zlim=log(c(40, 200)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(40, 80, 160)), tickLabels=c(40, 80, 160))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

##### Simulate EA locations for example ----
simDat = generateSimDataSetsLCPB2(nsim=1, targetPopMat=popMatSimpleAdjustedW, 
                                 fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                 logisticApproximation=FALSE, gridLevel=TRUE, 
                                 doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
                                 dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                 seed=1, inla.seed=1L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                 easpa=easpaW, popMat=popMatW, poppsub=poppsubW)
eaSamples = simDat$simulatedEAs$eaSamples
eaDat = simDat$simulatedEAs$eaDat
# eaDat = eaDat[eaDat$area == "Wajir",]
eaDatNW = eaDat[eaDat$subarea == "Wajir North", ]
eaSamplesNW = eaSamples[popGridWajir$subarea=="Wajir North",1]
# eaSamples = eaSamples[popGrid$area=="Wajir",1]
NSamples = simDat$simulatedEAs$aggregatedPop$pixelPop$NFineScalePrevalence[,1]
NSamplesNW = NSamples[popGridWajir$subarea=="Wajir North"]
prevalenceSamples = simDat$simulatedEAs$aggregatedPop$pixelPop$pFineScalePrevalence[,1]

# plot EA locations ----
pdf(file="figures/simpleExample/wajirSimEALocs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col="blue")
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot EA locations with color for fine scale prevalence
pRangeMod = range(c(eaDat$pFineScalePrevalence[eaDat$pFineScalePrevalence != 0], eaDat$pFineScaleRisk, eaDat$pSmoothRisk))
ticks = c(.01, .05, .1, .2, .3, .4, .5)
pdf(file="figures/simpleExample/wajirSimEALocsPrevalence.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pFineScalePrevalence, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot EA locations with color for fine scale risk
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRisk.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pFineScaleRisk, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot EA locations with color for fine scale risk on its own scale
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskSelfScaled.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pFineScaleRisk, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot EA locations with color for smooth risk
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRisk.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pSmoothRisk, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot EA locations with color for smooth risk on its own color scale
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRiskSelfScaled.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pSmoothRisk, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot North Wajir EA locations ----
pdf(file="figures/simpleExample/wajirSimEALocsNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           xlim = longRangeNorthWajir, ylim = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.3, col="blue")
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.5)
dev.off()

# plot EA locations with color for fine scale prevalence
pRangeMod = range(c(eaDatNW$pFineScalePrevalence[eaDatNW$pFineScalePrevalence != 0], eaDatNW$pFineScaleRisk, eaDatNW$pSmoothRisk))
# ticks = c(0, .02, .04, .06, .08, .1, .2, .3, .4)
# tickLabels = as.character(ticks)
# tickLabels[c(5)] = ""
# ticks = c(0, .02, .05, .1, .2, .3, .4)
ticks = c(0, .05, .1, .15, .2, .25, .3, .35, .4, .45)
tickLabels = as.character(ticks)
tickLabels[seq(4, 10, by=2)] = ""
pdf(file="figures/simpleExample/wajirSimEALocsPrevalenceNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pFineScalePrevalence, pch=21, cex=1, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, colorName="bg", 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# same but residuals
pdf(file="figures/simpleExample/wajirSimEALocsPrevalenceNWajirResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = eaDatNW$pFineScalePrevalence - eaDatNW$pSmoothRisk
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, prevalenceResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# same but percent residuals
pdf(file="figures/simpleExample/wajirSimEALocsPrevalenceNWajirPctResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = 100 * (eaDatNW$pFineScalePrevalence - eaDatNW$pSmoothRisk) / eaDatNW$pSmoothRisk
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, prevalenceResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# plot EA locations with color for fine scale risk
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pFineScaleRisk, pch=21, cex=1, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# same but different color scale
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajir2.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(eaDatNW$pFineScalePrevalence), center = median(eaDatNW$pFineScalePrevalence))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pFineScaleRisk, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=pRangeMod, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# same but residuals
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajirResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = eaDatNW$pFineScalePrevalence - eaDatNW$pSmoothRisk
riskResiduals = eaDatNW$pFineScaleRisk - eaDatNW$pSmoothRisk
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, riskResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# same but percent residuals
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajirPctResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = 100 * (eaDatNW$pFineScalePrevalence - eaDatNW$pSmoothRisk) / eaDatNW$pSmoothRisk
riskResiduals = 100 * (eaDatNW$pFineScaleRisk - eaDatNW$pSmoothRisk) / eaDatNW$pSmoothRisk
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, riskResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# plot EA locations with color for fine scale risk on its own scale
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskSelfScaledNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pFineScaleRisk, pch=21, cex=.8, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# plot EA locations with color for smooth risk
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRiskNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pSmoothRisk, pch=21, cex=.8, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# plot EA locations with color for smooth risk on its own color scale
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRiskSelfScaledNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pSmoothRisk, pch=21, cex=.8, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# plot pixel numerators ----
# plot pixel prevalence
pdf(file="figures/simpleExample/wajirSimPrevalencePixel.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
pRange = range(prevalenceSamples, na.rm=TRUE)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           prevalenceSamples, 
           col=riskCols, nx=45, ny=60, add.legend = TRUE, add=TRUE, 
           zlim=pRangeMod)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot smooth risk
expectedRisk = simDat$simulatedEAs$aggregatedPop$pixelPop$pSmoothRisk[,1]
pRangeMod = range(c(eaDat$pFineScalePrevalence[eaDat$pFineScalePrevalence != 0], eaDat$pFineScaleRisk, eaDat$pSmoothRisk))

pdf(file="figures/simpleExample/wajirSimSmoothRisk.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           logit(expectedRisk), 
           col=riskCols, nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=logit(pRangeMod))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
ticks = c(.01, .05, .1,.2, .3, .4)
tickLabels = as.character(ticks)
image.plot(zlim=logit(pRangeMod), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=logit(ticks), labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

pdf(file="figures/simpleExample/wajirSimSmoothRiskSelfScaled.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           logit(expectedRisk), 
           col=riskCols, nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=logit(range(expectedRisk)))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
ticks = c(.01, .02, .03, .04, .05, .06, .07, .08, .09, .1)
tickLabels = as.character(ticks)
image.plot(zlim=logit(range(expectedRisk)), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=logit(ticks), labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot smooth risk for Wajir North
expectedRisk = simDat$simulatedEAs$aggregatedPop$pixelPop$pSmoothRisk[,1]
expectedRisk = expectedRisk[popGridWajir$subarea == "Wajir North"]
# pRangeMod = range(c(eaDat$pLCPB[eaDat$pLCPB != 0], eaDat$pLCPb, eaDat$plcpb))
pRangeMod = range(c(eaDatNW$pFineScalePrevalence[eaDatNW$pFineScalePrevalence != 0], eaDatNW$pFineScaleRisk, eaDatNW$pSmoothRisk))

pdf(file="figures/simpleExample/wajirSimSmoothRiskNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           logit(expectedRisk), 
           col=riskCols, nx=29, ny=24, add.legend = FALSE, add=TRUE, 
           zlim=logit(pRangeMod))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
ticks = c(0, .05, .1, .15, .2, .25, .3, .35, .4, .45)
tickLabels = as.character(ticks)
tickLabels[seq(4, 10, by=2)] = ""
image.plot(zlim=logit(pRangeMod), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=logit(ticks), labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

pdf(file="figures/simpleExample/wajirSimSmoothRiskNWajirSelfScaled.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           logit(expectedRisk), 
           col=riskCols, nx=29, ny=24, add.legend = FALSE, add=TRUE, 
           zlim=logit(range(expectedRisk)))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
ticks = c(.02, .03, .04, .05, .06, .07)
tickLabels = as.character(ticks)
image.plot(zlim=logit(range(expectedRisk)), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=logit(ticks), labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

pdf(file="figures/simpleExample/wajirSimSmoothRiskNWajirSelfScaledLinear.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           expectedRisk, # nx=29, ny=24, 
           col=riskCols, add.legend = FALSE, add=TRUE, 
           zlim=range(expectedRisk), grid=gridListNW)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
ticks = c(.02, .03, .04, .05, .06, .07, .08)
tickLabels = as.character(ticks)
image.plot(zlim=range(expectedRisk), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=ticks, labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# plot denominators ----
# counts per pixel
pdf(file="figures/simpleExample/wajirSimEAsPerPixel.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
eaRange = c(1, 60)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           log(eaSamples), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(eaRange))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
eaTicks = c(1, 5, 10, 50)
eaTickLabels = as.character(eaTicks)
image.plot(zlim=log(eaRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(eaTicks), labels=eaTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# same but with number of assigned people per pixel
pdf(file="figures/simpleExample/wajirSimNPerPixel.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
NRange = c(35, 2400)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           log(NSamples), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(NRange))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
NTicks = c(35, 100, 300, 800, 2400)
NTickLabels = as.character(NTicks)
image.plot(zlim=log(NRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(NTicks), labels=NTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# treating denominators as weights
eaSamplesNorm = eaSamples/sum(eaSamples)
NSamples = NSamples/sum(NSamples)
popWeights = popGridWajirAdjusted$pop / sum(popGridWajirAdjusted$pop)

nonZeroRange = range(c(eaSamplesNorm[eaSamplesNorm != 0], 
                       NSamples[NSamples != 0], 
                       popWeights[popWeights != 0]))
weightRange = nonZeroRange
weightTicks = c(1e-04, 1e-03, 1e-02)
weightTickLabels = as.character(weightTicks)

weightRange2 = range(c(eaSamplesNorm[eaSamplesNorm!= 0], NSamples[NSamples != 0]))
weightTicks2 = c(1e-04, 1e-03, 1e-02)
weightTickLabels2 = as.character(weightTicks2)

# plot pop density with same scale
if(FALSE) {
  
  pdf(file="figures/simpleExample/wajirPopDensityNorm.pdf", width=4, height=5)
  par(mar=c(4.1, 4.1, 1.1, 4.5))
  plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
             kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
             leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
  image.plot(zlim=log(weightRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add=TRUE, axis.args=list(at=log(weightTicks), labels=weightTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
             legend.cex=.5, legend.width=1)
  quilt.plot(popGridWajirAdjusted$lon, popGridWajirAdjusted$lat, log(popWeights), 
             col=popCols, nx=45, ny=75, add.legend = FALSE, add=TRUE, 
             zlim=log(weightRange))
  plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
  addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
  dev.off()
}

# same but with counts per pixel
pdf(file="figures/simpleExample/wajirSimEAsPerPixelNorm.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGridWajir$lon, popGridWajir$lat, log(eaSamplesNorm), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
image.plot(zlim=log(weightRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks), labels=weightTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# same but with number of assigned people per pixel
pdf(file="figures/simpleExample/wajirSimNPerPixelNorm.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           log(NSamples), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# NTicks = c(35, 100, 300, 800, 2400)
# NTickLabels = as.character(NTicks)
image.plot(zlim=log(weightRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks), labels=weightTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# same but with number of assigned people per EA
weightsEA = rep(1/815, 815)
pdf(file="figures/simpleExample/wajirSimNPerEANorm.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotWithColor(eaDat$lon, eaDat$lat, weightsEA, 
              colScale=popCols[-(1:5)], add.legend = FALSE, new=FALSE, 
              zlim=weightRange, scaleFun=log, scaleFunInverse=exp, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              pch=19, cex=.2, ticks=weightTicks)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# NTicks = c(35, 100, 300, 800, 2400)
# NTickLabels = as.character(NTicks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# same but with counts per pixel
pdf(file="figures/simpleExample/wajirSimEAsPerPixelNorm2.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           log(eaSamples), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange2))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# same but with number of assigned people per pixel
pdf(file="figures/simpleExample/wajirSimNPerPixelNorm2.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           log(NSamples), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange2))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
NTicks = c(35, 100, 300, 800, 2400)
NTickLabels = as.character(NTicks)
image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# no point in hits plot, since every EA has 25 neonatals
# live births per EA
# pdf(file="figures/simpleExample/wajirSimNPerEA.pdf", width=4, height=5)
# par(mar=c(4.1, 4.1, 1.1, 4.5))
# plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
#            kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
#            leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
# plotWithColor(eaDat$lon, eaDat$lat, eaDat$n, 
#               colScale=popCols[-c(1:5)], new=FALSE, 
#               legend.cex=.5, legend.width=1, legend.mar=5, 
#               pch=19, cex=.2)
# # quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
# #            popGrid$lat[popGrid$area=="Wajir"], 
# #            log(NSamples), 
# #            col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
# #            zlim=log(weightRange2))
# plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# # NTicks = c(35, 100, 300, 800, 2400)
# # NTickLabels = as.character(NTicks)
# # image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
# #            col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
# #            legend.cex=.5, legend.width=1)
# addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
# dev.off()

# plot agg weights North Wajir ----

# treating denominators as weights
eaSamplesNorm = eaSamplesNW/sum(eaSamplesNW)
NSamplesNW = NSamplesNW/sum(NSamplesNW)
popWeights = popMatSimpleAdjustedNW$pop / sum(popMatSimpleAdjustedNW$pop)

nonZeroRange = range(c(eaSamplesNorm[eaSamplesNorm != 0], 
                       NSamplesNW[NSamplesNW != 0], 
                       popWeights[popWeights != 0]))
weightRange = nonZeroRange
weightTicks = c(1e-04, 1e-03, 1e-02)
weightTickLabels = as.character(weightTicks)

weightRange2 = range(c(eaSamplesNorm[eaSamplesNorm!= 0], NSamplesNW[NSamplesNW != 0]))
weightTicks2 = c(1e-04, 1e-03, 1e-02)
weightTickLabels2 = as.character(weightTicks2)

# first generate the North Wajir Grid
gridNWEast = sort(unique(popMatSimpleAdjustedNW$east))
gridNWNorth = sort(unique(popMatSimpleAdjustedNW$north))
gridNWEastMid = median(gridNWEast)
gridNWNorthMid = median(gridNWNorth)
gridNWHorizontal = cbind(gridNWEast, gridNWNorthMid)
gridNWVertical = cbind(gridNWEastMid, gridNWNorth)
gridNWLon = projKenya(gridNWHorizontal, inverse=TRUE)[,1]
gridNWLat = projKenya(gridNWVertical, inverse=TRUE)[,2]
gridListNW = list(x=gridNWLon, y=gridNWLat)

# plot pop density with same scale
if(FALSE) {
  
  pdf(file="figures/simpleExample/wajirPopDensityNormNW.pdf", width=5, height=3.9)
  par(mar=c(4.1, 4.1, 1.1, 4.5))
  plotMapDat(mapDat=adm2[adm2@data$NAME_2=="Wajir North",], new=TRUE, 
             kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
             leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
  image.plot(zlim=log(weightRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add=TRUE, axis.args=list(at=log(weightTicks), labels=weightTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
             legend.cex=.5, legend.width=1)
  quilt.plot(popMatSimpleAdjustedNW$lon, popMatSimpleAdjustedNW$lat, log(popWeights), 
             col=popCols, add.legend = FALSE, add=TRUE, # nx=29, ny=24, 
             zlim=log(weightRange), grid=gridListNW)
  plotMapDat(mapDat=adm2[adm2@data$NAME_2=="Wajir North",], lwd=.5)
  addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
  dev.off()
}

# same but with number of assigned people per EA

nEAsNW = sum(eaSamples[popGridWajir$subarea == "Wajir North"])
weightsEA = rep(1/nEAsNW, nEAsNW)
pdf(file="figures/simpleExample/wajirSimNPerEANormNW.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm2[adm2@data$NAME_2=="Wajir North",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotWithColor(eaDatNW$lon, eaDatNW$lat, weightsEA, 
              colScale=popCols, add.legend = FALSE, new=FALSE, 
              zlim=weightRange, scaleFun=log, scaleFunInverse=exp, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              pch=19, cex=1, ticks=weightTicks)
# plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# NTicks = c(35, 100, 300, 800, 2400)
# NTickLabels = as.character(NTicks)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.8)
dev.off()

# same but with counts per pixel
pdf(file="figures/simpleExample/wajirSimEAsPerPixelNormNW2.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm2[adm2@data$NAME_2=="Wajir North",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           log(eaSamplesNorm), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange2))
# plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# same but with number of assigned people per pixel
pdf(file="figures/simpleExample/wajirSimNPerPixelNormNW.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm2[adm2@data$NAME_2=="Wajir North",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           log(NSamplesNW), grid=gridListNW, # nx=29, ny=24, 
           col=popCols, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange))
# plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
NTicks = c(35, 100, 300, 800, 2400)
NTickLabels = as.character(NTicks)
image.plot(zlim=log(weightRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols, add=TRUE, axis.args=list(at=log(weightTicks), labels=weightTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_2=="Wajir North",], offsets=offsets, cex=.4)
dev.off()

# same but with number of assigned people per pixel and weight range 2
pdf(file="figures/simpleExample/wajirSimNPerPixelNormNW2.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm2[adm2@data$NAME_2=="Wajir North",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           log(NSamplesNW), 
           col=popCols, nx=29, ny=24, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange2))
# plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
NTicks = c(35, 100, 300, 800, 2400)
NTickLabels = as.character(NTicks)
image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_2=="Wajir North",], offsets=offsets, cex=.4)
dev.off()

# no point in hits plot, since every EA has 25 neonatals
# live births per EA
# pdf(file="figures/simpleExample/wajirSimNPerEA.pdf", width=4, height=5)
# par(mar=c(4.1, 4.1, 1.1, 4.5))
# plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
#            kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
#            leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
# plotWithColor(eaDat$lon, eaDat$lat, eaDat$n, 
#               colScale=popCols[-c(1:5)], new=FALSE, 
#               legend.cex=.5, legend.width=1, legend.mar=5, 
#               pch=19, cex=.2)
# # quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
# #            popGrid$lat[popGrid$area=="Wajir"], 
# #            log(NSamplesNW), 
# #            col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
# #            zlim=log(weightRange2))
# plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# # NTicks = c(35, 100, 300, 800, 2400)
# # NTickLabels = as.character(NTicks)
# # image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
# #            col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
# #            legend.cex=.5, legend.width=1)
# addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
# dev.off()

# check pixel constituencies ----
pdf(file="figures/simpleExample/wajirPixelConstituencyCheck.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
tempVals = factor(popGrid$subarea[popGrid$area=="Wajir"])
tempVals = as.numeric(tempVals)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           tempVals, 
           col=rainbow(6), nx=45, ny=60, add.legend = TRUE, add=TRUE, 
           FUN=max)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot constituency denom ----
urbanPixels = popGridWajir$urban
eaSamplesUrban = eaSamples
eaSamplesRural = eaSamples
eaSamplesUrban[!urbanPixels] = 0
eaSamplesRural[urbanPixels] = 0

out = aggregate(c(eaSamplesUrban), by=list(constituency=popGridWajir$subarea), FUN=sum)
nEAs = out$x
nEAsUrban = nEAs

out = aggregate(c(eaSamplesRural), by=list(constituency=popGridWajir$subarea), FUN=sum)
nEAs = out$x
nEAsRural = nEAs

out = aggregate(c(eaSamples), by=list(constituency=popGridWajir$subarea), FUN=sum)
nEAs = out$x
nEAsTotal = nEAs

# plot the number of sampled EAs per constituency
pdf(file="figures/simpleExample/wajirSimEAs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(nEAs), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=urbanCols, zlim=log(c(70, 220)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(80, 120, 160, 200)), tickLabels=c(80, 120, 160, 200))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot the number of sampled Urban EAs per constituency
pdf(file="figures/simpleExample/wajirSimEAsUrban.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(nEAs), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=urbanCols, zlim=log(c(2, 200)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(2, 5, 50, 150)), tickLabels=c(2, 5, 50, 150))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot the number of sampled Urban EAs per constituency
pdf(file="figures/simpleExample/wajirSimEAsRural.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(plotVar=log(nEAs), varCounties=constituenciesW, mapDat=adm2, lwd=.5, 
           cols=ruralCols, zlim=log(c(40, 200)), forceColorsInRange=TRUE, 
           legend.width=1, legend.cex=.5, legend.mar=5, 
           axis.args=list(cex.axis=1, tck=-.7, hadj=-.1), 
           ticks=log(c(40, 80, 160)), tickLabels=c(40, 80, 160))
# image.plot(zlim=log(popRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
#            col=popCols, add=TRUE, axis.args=list(at=log(popTicks), labels=popTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
#            legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2[adm2@data$NAME_1=="Wajir",], offsets=offsets, areaVarName="NAME_2", cex=.4)
dev.off()

# plot prevalence and risk ----
easInWajir = eaDat$area == "Wajir"
pdf(file="figures/simpleExample/wajirSimSmoothRiskVFineScalePrevalence.pdf", width=5, height=5)
plot(eaDat$pSmoothRisk[easInWajir], eaDat$pFineScalePrevalence[easInWajir], pch=19, cex=.1, col="blue", type="n", 
     xlab="Smooth Risk", ylab="Fine Scale Prevalence")
abline(a=0, b=1)
points(eaDat$pSmoothRisk[easInWajir], eaDat$pFineScalePrevalence[easInWajir], pch=19, cex=.1, col="blue")
dev.off()

# plot constituency level prevalence versus risk:
sortI = sort(as.character(poppsubKenya$subarea), index.return=TRUE)$ix
theseCounties = poppsubKenya$area[sortI]
# constituenciesInWajir = theseCounties == "Wajir"
pdf(file="figures/simpleExample/wajirSimSmoothRiskVFineScalePrevalenceConstituency.pdf", width=5, height=5)
plot(simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pSmoothRisk[,1], 
     simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence[,1], 
     pch=19, cex=1, col="blue", type="n", 
     xlab="Smooth Risk", ylab="Fine Scale Prevalence", xlim=c(.03, .1), ylim=c(.03, .1))
abline(a=0, b=1)
points(simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pSmoothRisk[,1], 
       simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence[,1], 
       pch=19, cex=1, col="blue")
text(simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pSmoothRisk[,1], 
     simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence[,1], 
     constituenciesW, cex=.55, pos=4)
dev.off()

# Tabulate prevalence versus risk ----
# calulcate percent different between risk and prevalence
SmoothRisk = simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pSmoothRisk[,1]
Risk = simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScaleRisk[,1]
Prevalence = simDat$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence[,1]
constituencyResults = data.frame(Area=constituenciesW, SmoothRisk=SmoothRisk, Risk=Risk, Prevalence=Prevalence, TotalEAs=nEAsTotal)
wajirCountyI = sort(unique(poppsubKenya$area)) == "Wajir"
SmoothRiskCounty = simDat$simulatedEAs$aggregatedPop$areaPop$aggregationResults$pSmoothRisk[,1]
FineScaleRiskCounty = simDat$simulatedEAs$aggregatedPop$areaPop$aggregationResults$pFineScaleRisk[,1]
FineScalePrevalenceCounty = simDat$simulatedEAs$aggregatedPop$areaPop$aggregationResults$pFineScalePrevalence[,1]
countyResults = data.frame(Area="Wajir", SmoothRisk=SmoothRiskCounty, Risk=FineScaleRiskCounty, Prevalence=FineScalePrevalenceCounty, TotalEAs=sum(nEAsTotal))
combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
#        Area SmoothRisk       Risk Prevalence TotalEAs
#       Eldas 0.05990297 0.05730967 0.06254545       55
#      Tarbaj 0.05111533 0.05100725 0.04200000      120
#  Wajir East 0.02777313 0.02942914 0.03300546      183
# Wajir North 0.05683298 0.05581496 0.06013072      153
# Wajir South 0.05562260 0.05899186 0.05500000      184
#  Wajir West 0.05388420 0.04814104 0.04533333      120
#       Wajir 0.04910678 0.04887061 0.04819632      815

percentDifference = 100*(Prevalence - SmoothRisk) / SmoothRisk
print(data.frame(Constituency=constituenciesW, pctDiff=percentDifference, urbanEAs=nEAsUrban, ruralEAs=nEAsRural, totalEAs=nEAsTotal))
# Constituency    pctDiff urbanEAs ruralEAs totalEAs
#        Eldas   4.411283        0       55       55
#       Tarbaj -17.832874        7      113      120
#   Wajir East  18.839540      135       48      183
#  Wajir North   5.802506        3      150      153
#  Wajir South  -1.119322        5      179      184
#   Wajir West -15.868966       13      107      120

#### Plot Constituency Results ####
ylim=range(c(constituencyResults$SmoothRisk, constituencyResults$Risk, constituencyResults$Prevalence))
smoothRange = range(constituencyResults$SmoothRisk)
pdf(file="figures/simpleExample/pairPlotConstituency.pdf", width=5, height=5)
plot(constituencyResults$SmoothRisk, constituencyResults$Risk, 
     pch=19, col="blue", 
     xlim=smoothRange, ylim=ylim, 
     main="", xlab="Smooth Risk", ylab="Risk or Prevalence")
points(constituencyResults$SmoothRisk, constituencyResults$Prevalence, pch=17, col="purple")
abline(0, 1)
legend("top", c("Risk", "Prevalence"), pch=c(19, 17), col=c("blue", "purple"))
dev.off()

##### Simulate many populations, illustrate uncertainty ----

# run simulation
# nsim = 1000
# simDatMultiOriginal = generateSimDataSetsLCPB(nsim=nsim, adjustedPopMat=popMatSimpleAdjusted, fixPopPerEA=25, logisticApproximation=FALSE, 
#                                               dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
#                                               simPopOnly=TRUE, returnEAinfo=FALSE, seed=1, inla.seed=1L, verbose=TRUE, 
#                                               easpa=easpa, popMat=popMat, constituencyPop=poppsubKenya)
nsim = 10000
simDatMulti = generateSimDataSetsLCPB2(nsim=nsim, targetPopMat=popMatSimpleAdjustedW, 
                                       fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      simPopOnly=TRUE, returnEAinfo=FALSE, seed=1, inla.seed=1L, verbose=TRUE, 
                                      easpa=easpaW, popMat=popMatW, poppsub=poppsubW, gridLevel=TRUE, 
                                      doSmoothRisk=TRUE, doFineScaleRisk=TRUE, logisticApproximation=FALSE)
testnsim = 10000
prevalenceSamples = simDatMulti$subareaPop$aggregationResults$ZFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = simDatMulti$subareaPop$aggregationResults$ZFineScaleRisk[,1:testnsim]
smoothRiskSamples = simDatMulti$subareaPop$aggregationResults$ZSmoothRisk[,1:testnsim]

prevalenceSamples = simDatMulti$subareaPop$aggregationResults$pFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = simDatMulti$subareaPop$aggregationResults$pFineScaleRisk[,1:testnsim]
smoothRiskSamples = simDatMulti$subareaPop$aggregationResults$pSmoothRisk[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
# constituencyResults = rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths)
constituencyResults = data.frame(Area=constituenciesW, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, TotalEAs=nEAsTotal)

smoothRiskCounty = simDatMulti$areaPop$aggregationResults$ZSmoothRisk
fineScaleRiskCounty = simDatMulti$areaPop$aggregationResults$ZFineScaleRisk
fineScalePrevalenceCounty = simDatMulti$areaPop$aggregationResults$ZFineScalePrevalence

smoothRiskCounty = simDatMulti$areaPop$aggregationResults$pSmoothRisk
fineScaleRiskCounty = simDatMulti$areaPop$aggregationResults$pFineScaleRisk
fineScalePrevalenceCounty = simDatMulti$areaPop$aggregationResults$pFineScalePrevalence

prevalenceCIWidths = apply(fineScalePrevalenceCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
countyResults = data.frame(Area="Wajir", SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, TotalEAs=sum(nEAsTotal))

combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
#        Area SmoothRisk       Risk Prevalence TotalEAs
#       Eldas 0.04629711 0.04800643 0.05072920       55
#      Tarbaj 0.04583046 0.04671403 0.04811827      120
#  Wajir East 0.02576206 0.02633897 0.02747312      183
# Wajir North 0.04622396 0.04630701 0.04750958      153
# Wajir South 0.04192406 0.04236640 0.04348542      184
#  Wajir West 0.04240138 0.04367990 0.04497796      120
#       Wajir 0.03533783 0.03554017 0.03577914      815
combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
# pct larger:
100 * (constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk
# 9.573133 4.991903 6.641795 2.781280 3.724255 6.076628
mean(100 * (constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk)
# 5.631499

# 95% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths)
#                            Eldas     Tarbaj Wajir East Wajir North Wajir South  Wajir West
# prevalenceCIWidths    0.07889456 0.07577971 0.04331512  0.07446680  0.06694146  0.07016640
# fineScaleRiskCIWidths 0.07349561 0.07416060 0.04242293  0.07294329  0.06555293  0.06852538
# smoothRiskCIWidths    0.07158451 0.07203241 0.04137768  0.07182707  0.06456175  0.06657485


mean(100* (prevalenceCIWidths - smoothRiskCIWidths) / smoothRiskCIWidths)
# 5.475363

# means
prevalenceMeans = rowMeans(prevalenceSamples)
fineScaleRiskMeans = rowMeans(fineScaleRiskSamples)
smoothRiskMeans = rowMeans(smoothRiskSamples)
rbind(prevalenceMeans, fineScaleRiskMeans, smoothRiskMeans)
#                         Eldas     Tarbaj Wajir East Wajir North Wajir South Wajir West
# prevalenceMeans    0.06437398 0.06314970 0.03486189  0.06359500  0.06244569 0.05976299
# fineScaleRiskMeans 0.06425381 0.06315893 0.03487984  0.06357393  0.06246475 0.05969088
# smoothRiskMeans    0.06424268 0.06314229 0.03488853  0.06364964  0.06243501 0.05969217

# SDs
means = smoothRiskMeans
prevalenceErr = sweep(prevalenceSamples, 1, means, "-")
fineScaleRiskErr = sweep(fineScaleRiskSamples, 1, means, "-")
smoothRiskErr = sweep(smoothRiskSamples, 1, means, "-")
prevalenceSDs = rowMeans(prevalenceErr^2)
fineScaleRiskSDs = rowMeans(fineScaleRiskErr^2)
smoothRiskSDs = rowMeans(smoothRiskErr^2)

rbind(prevalenceSDs, fineScaleRiskSDs, smoothRiskSDs)
#                         Eldas       Tarbaj   Wajir East  Wajir North  Wajir South   Wajir West
# prevalenceSDs    0.0004158730 0.0003728951 0.0001247979 0.0003671224 0.0003012047 0.0003255839
# fineScaleRiskSDs 0.0003737302 0.0003546713 0.0001167018 0.0003524427 0.0002874843 0.0003060737
# smoothRiskSDs    0.0003448184 0.0003413493 0.0001115633 0.0003432380 0.0002788463 0.0002911358

# pct difference
100 * (prevalenceSDs - smoothRiskSDs) / smoothRiskSDs
mean(100 * (prevalenceSDs - smoothRiskSDs) / smoothRiskSDs)
# 11.41997

##### Fit model to SRS cluster data ----
# get the data
dat = simDat$SRSDat$clustDat[[1]]
dat$y = dat$Z
dat$n = dat$N

# fit SPDE cluster level risk model
nSamples = 10000 # 10000 takes about 5 minutes
spdeFit = fitSPDEKenyaDat(dat, nPostSamples=nSamples, popMat=popMatW)

# obtain model output
uDraws = spdeFit$uDraws
sigmaEpsilonDraws = spdeFit$sigmaEpsilonDraws

# apply aggregation models
aggResults = simPopCustom(uDraws, sigmaEpsilonDraws, easpa=easpaW, 
                          popMat=popMatW, targetPopMat=popMatSimpleAdjustedW, 
                          stratifyByUrban=TRUE, 
                          doFineScaleRisk=TRUE, doSmoothRisk=TRUE, doGriddedRisk=TRUE, 
                          doSmoothRiskLogisticApprox=FALSE, 
                          poppsub=poppsubW, subareaLevel=TRUE, gridLevel=TRUE, 
                          fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                          returnEAinfo=FALSE, verbose=TRUE)

# get average number of EAs per area
out = meanEAsPerCon2()
out = out[out$area == "Wajir",]
nEAsTotal = out$meanTotalEAs

# get average number of neonatals per area
nPerCon = aggResults$subareaPop$aggregationResults$NSmoothRisk[,1]

# calculate number of pixels per area (for understanding IHME model variance)
out = aggregate(popMatSimpleAdjustedW$pop, by=list(pixels=popMatSimpleAdjustedW$subarea), FUN=length)
nPixels = out$x

# calculate number of EAs sampled per area (for understanding spatial variance)
out = aggregate(dat$N[dat$area == "Wajir"], by=list(dat$subarea[dat$area == "Wajir"]), FUN=length)
nEAsSampled = out$x

# gather aggregation model output
prevalenceSamples = aggResults$subareaPop$aggregationResults$ZFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = aggResults$subareaPop$aggregationResults$ZFineScaleRisk[,1:testnsim]
smoothRiskSamples = aggResults$subareaPop$aggregationResults$ZSmoothRisk[,1:testnsim]
ihmeSamples = aggResults$subareaPop$aggregationResults$ZGriddedRisk[,1:testnsim]

prevalenceSamples = aggResults$subareaPop$aggregationResults$pFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = aggResults$subareaPop$aggregationResults$pFineScaleRisk[,1:testnsim]
smoothRiskSamples = aggResults$subareaPop$aggregationResults$pSmoothRisk[,1:testnsim]
ihmeSamples = aggResults$subareaPop$aggregationResults$pGriddedRisk[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
constituencyResults = data.frame(Area=constituenciesW, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=nEAsTotal, EAsSampled=nEAsSampled, Pixels=nPixels, Neonatals=nPerCon)

smoothRiskCounty = aggResults$areaPop$aggregationResults$ZSmoothRisk
fineScaleRiskCounty = aggResults$areaPop$aggregationResults$ZFineScaleRisk
fineScalePrevalenceCounty = aggResults$areaPop$aggregationResults$ZFineScalePrevalence
ihmeCounty = aggResults$areaPop$aggregationResults$ZGriddedRisk

smoothRiskCounty = aggResults$areaPop$aggregationResults$pSmoothRisk
fineScaleRiskCounty = aggResults$areaPop$aggregationResults$pFineScaleRisk
fineScalePrevalenceCounty = aggResults$areaPop$aggregationResults$pFineScalePrevalence
ihmeCounty = aggResults$areaPop$aggregationResults$pGriddedRisk

prevalenceCIWidths = apply(fineScalePrevalenceCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
countyResults = data.frame(Area="Wajir", SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=sum(nEAsTotal), EAsSampled=sum(nEAsSampled), Pixels=sum(nPixels), Neonatals=sum(nPerCon))

combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
xtable(combinedResults, 
       digits= c(0, 0, 3, 3, 3, 3, 0, 0, 0, 0), 
       display=c("d", "s", "e", "e", "e", "e", "d", "d", "d", "d"))

# pct larger:
100 * (constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk
mean(100 * (constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk)

# 95% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths, ihmeCIWidths)

mean(100 * (prevalenceCIWidths - smoothRiskCIWidths) / smoothRiskCIWidths)

# means
prevalenceMeans = rowMeans(prevalenceSamples)
fineScaleRiskMeans = rowMeans(fineScaleRiskSamples)
smoothRiskMeans = rowMeans(smoothRiskSamples)
ihmeMeans = rowMeans(ihmeSamples)
rbind(prevalenceMeans, fineScaleRiskMeans, smoothRiskMeans, ihmeMeans)

# SDs
means = smoothRiskMeans
prevalenceErr = sweep(prevalenceSamples, 1, means, "-")
fineScaleRiskErr = sweep(fineScaleRiskSamples, 1, means, "-")
smoothRiskErr = sweep(smoothRiskSamples, 1, means, "-")
ihmeErr = sweep(ihmeSamples, 1, means, "-")
prevalenceSDs = rowMeans(prevalenceErr^2)
fineScaleRiskSDs = rowMeans(fineScaleRiskErr^2)
smoothRiskSDs = rowMeans(smoothRiskErr^2)
ihmeSDs = rowMeans(ihmeErr^2)

rbind(prevalenceSDs, fineScaleRiskSDs, smoothRiskSDs, ihmeSDs)

# pct difference
100 * (prevalenceSDs - smoothRiskSDs) / smoothRiskSDs
mean(100 * (prevalenceSDs - smoothRiskSDs) / smoothRiskSDs)

# Do the same, but in Nairobi ----
longRangeNairobi = c(39, 41)
latRangeNairobi = c(0.25, 3.6)
constituenciesN = poppsubKenya$subarea[poppsubKenya$area=="Nairobi"]
offsets = matrix(0, nrow=4, ncol=2)
# offsets[1,2] = .1 # shift label for Eldas slightly higher
# offsets[6,1] = .15 # shift label for Wajir West slightly farther east
easpsub = meanEAsPerCon2()
easpsub = easpsub[easpsub$area=="Nairobi",]
popGridN = popGrid[popGrid$area=="Nairobi",]
# popGridWajirAdjusted = popGridAdjusted[popGridAdjusted$area=="Wajir",]
# plotMapDat(mapDat=adm0, lwd=.5, new=TRUE)
# plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], border="green")

# restrict simulation to Nairobi
popMatSimpleAdjustedN = popMatSimpleAdjusted[popMatSimpleAdjusted$area == "Nairobi", ]
easpaN = makeDefaultEASPA()
easpaN = easpaN[easpaN$area == "Nairobi", ]
popMatN = makeDefaultPopMat()
popMatN = popMatN[popMatN$area == "Nairobi", ]
poppsubN = poppsubKenya
poppsubN = poppsubN[poppsubN$area == "Nairobi", ]

simDatKenya = generateSimDataSetsLCPB2(nsim=1, targetPopMat=popMatSimpleAdjusted, 
                                      fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                      logisticApproximation=FALSE, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      seed=1, inla.seed=1L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                      easpa=NULL, popMat=NULL, poppsub=poppsubKenya)

# get the data and the true population
dat = simDatKenya$SRSDat$clustDat[[1]]
truePrevalenceConstituencyKenya = simDatKenya$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence[,1]
truePrevalenceCountyKenya = simDatKenya$simulatedEAs$aggregatedPop$areaPop$aggregationResults$pFineScalePrevalence[,1]

# fit SPDE cluster level risk model
nSamples = 10000
spdeFitN = fitSPDEKenyaDat(dat, nPostSamples=nSamples, popMat=popGridN)

# obtain model output
uDraws = spdeFitN$uDraws
sigmaEpsilonDraws = spdeFitN$sigmaEpsilonDraws

# apply aggregation models
aggResultsN = modLCPB(uDraws, sigmaEpsilonDraws, easpaN, popMatN, 
                      popMatSimpleAdjustedN, doLCPb=TRUE, doIHMEModel=TRUE, 
                      constituencyPop=poppsubN, ensureAtLeast1PerConstituency=TRUE, 
                      logisticApproximation=FALSE, verbose=TRUE, 
                      fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1)

# get average number of EAs per area
out = meanEAsPerCon()
out = out[outarea == "Nairobi",]
nEAsTotal = out$meanTotalEAs

# get average number of neonatals per area
nPerCon = aggResultsN$subareaPop$aggregationResults$NSmoothRisk[,1]

# calculate number of pixels per area (for understanding IHME model variance)
out = aggregate(popMatSimpleAdjustedN$pop, by=list(pixels=popMatSimpleAdjustedN$subarea), FUN=length)
nPixels = out$x

# calculate number of EAs sampled per area (for understanding spatial variance)
out = aggregate(dat$N[dat$area == "Nairobi"], by=list(dat$subarea[dat$area == "Nairobi"]), FUN=length)
nEAsSampled = out$x

# gather aggregation model output
prevalenceSamples = aggResultsN$aggregationResults$ZFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = aggResultsN$aggregationResults$ZFineScaleRisk[,1:testnsim]
smoothRiskSamples = aggResultsN$aggregationResults$ZSmoothRisk[,1:testnsim]
ihmeSamples = aggResultsN$aggregationResultsIHME$constituencyMatrices$Z[,1:testnsim]

prevalenceSamples = aggResultsN$subareaPop$aggregationResults$pFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = aggResultsN$subareaPop$aggregationResults$pFineScaleRisk[,1:testnsim]
smoothRiskSamples = aggResultsN$subareaPop$aggregationResults$pSmoothRisk[,1:testnsim]
ihmeSamples = aggResultsN$aggregationResultsIHME$constituencyMatrices$p[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
constituencyResults = data.frame(Area=constituenciesN, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=nEAsTotal, EAsSampled=nEAsSampled, Pixels=nPixels, Neonatals=nPerCon)

smoothRiskCounty = aggResultsN$aggregationResultslcpb$countyMatrices$Z
fineScaleRiskCounty = aggResultsN$aggregationResultsLCPb$countyMatrices$Z
fineScalePrevalenceCounty = aggResultsN$aggregationResultsLCPB$countyMatrices$Z
ihmeCounty = aggResultsN$aggregationResultsIHME$countyMatrices$Z

smoothRiskCounty = aggResultsN$areaPop$aggregationResults$pSmoothRisk
fineScaleRiskCounty = aggResultsN$areaPop$aggregationResults$pFineScaleRisk
fineScalePrevalenceCounty = aggResultsN$areaPop$aggregationResults$pFineScalePrevalence
ihmeCounty = aggResultsN$aggregationResultsIHME$countyMatrices$p

prevalenceCIWidths = apply(fineScalePrevalenceCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
countyResults = data.frame(Area="Nairobi", SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=sum(nEAsTotal), EAsSampled=sum(nEAsSampled), Pixels=sum(nPixels), Neonatals=sum(nPerCon))

combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
combinedResults$Area = as.character(combinedResults$Area)
# combinedResults$Area[1:4] = as.character(1:4)

xtable(combinedResults, 
       digits= c(0, 0, 2, 2, 2, 2, 0, 0, 0, 0), 
       display=c("d", "s", "e", "e", "e", "e", "d", "d", "d", "d"))

# pct larger:
(constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk
mean((constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk)

# 95% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths, ihmeCIWidths)

mean((prevalenceCIWidths - smoothRiskCIWidths) / smoothRiskCIWidths)

# means
prevalenceMeans = rowMeans(prevalenceSamples)
fineScaleRiskMeans = rowMeans(fineScaleRiskSamples)
smoothRiskMeans = rowMeans(smoothRiskSamples)
ihmeMeans = rowMeans(ihmeSamples)
rbind(prevalenceMeans, fineScaleRiskMeans, smoothRiskMeans, ihmeMeans)

# SDs
means = smoothRiskMeans
prevalenceErr = sweep(prevalenceSamples, 1, means, "-")
fineScaleRiskErr = sweep(fineScaleRiskSamples, 1, means, "-")
smoothRiskErr = sweep(smoothRiskSamples, 1, means, "-")
ihmeErr = sweep(ihmeSamples, 1, means, "-")
prevalenceSDs = rowMeans(prevalenceErr^2)
fineScaleRiskSDs = rowMeans(fineScaleRiskErr^2)
smoothRiskSDs = rowMeans(smoothRiskErr^2)
ihmeSDs = rowMeans(ihmeErr^2)

rbind(prevalenceSDs, fineScaleRiskSDs, smoothRiskSDs, ihmeSDs)

# pct difference
(prevalenceSDs - smoothRiskSDs) / smoothRiskSDs
mean((prevalenceSDs - smoothRiskSDs) / smoothRiskSDs)

##### Do the same, but in all of Kenya ----
if(FALSE) {
  # obtain simulated data over all of Kenya
  dat = simDatKenya$SRSDat$clustDat[[1]]
  
  # simulate results over all of Kenya
  # (ONLY DO THIS ON THE CLUSTER)
  nSamples = 10000
  spdeFitKenya = fitSPDEKenyaDat(dat, nPostSamples=nSamples, popMat=popGrid)
  
  # obtain model output
  uDraws = spdeFitKenya$uDraws
  sigmaEpsilonDraws = spdeFitKenya$sigmaEpsilonDraws
  
  # apply aggregation models
  aggResultsKenya = modLCPB(uDraws, sigmaEpsilonDraws, NULL, NULL, 
                            NULL, doLCPb=TRUE, doIHMEModel=TRUE, 
                            constituencyPop=poppsubKenya, ensureAtLeast1PerConstituency=TRUE, 
                            logisticApproximation=FALSE, verbose=TRUE, 
                            fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1)
  
  save(aggResultsKenya, file="savedOutput/simpleExample/aggResultsKenya.RData")
}

# get results
out = load("savedOutput/simpleExample/aggResultsKenya.RData")

# get number of EAs per constituency
out = meanEAsPerCon()
nEAsTotal = out$meanTotalEAs

# get average number of neonatals per area
nPerCon = rowMeans(aggResultsKenya$subareaPop$aggregationResults$NSmoothRisk)

# calculate number of pixels per area (for understanding IHME model variance)
out = aggregate(popMatSimpleAdjustedN$pop, by=list(pixels=popMatSimpleAdjustedN$subarea), FUN=length)
nPixels = out$x

# calculate number of EAs sampled per area (for understanding spatial variance)
out = aggregate(dat$N[dat$area == "Nairobi"], by=list(dat$subarea[dat$area == "Nairobi"]), FUN=length)
nEAsSampled = out$x

# gather aggregation model output
prevalenceSamples = aggResultsKenya$aggregationResults$ZFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = aggResultsKenya$aggregationResults$ZFineScaleRisk[,1:testnsim]
smoothRiskSamples = aggResultsKenya$aggregationResults$ZSmoothRisk[,1:testnsim]
ihmeSamples = aggResultsKenya$aggregationResultsIHME$constituencyMatrices$Z[,1:testnsim]

prevalenceSamples = aggResultsKenya$subareaPop$aggregationResults$pFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = aggResultsKenya$subareaPop$aggregationResults$pFineScaleRisk[,1:testnsim]
smoothRiskSamples = aggResultsKenya$subareaPop$aggregationResults$pSmoothRisk[,1:testnsim]
ihmeSamples = aggResultsKenya$aggregationResultsIHME$constituencyMatrices$p[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
constituencyResults = data.frame(Area=constituenciesN, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=nEAsTotal, EAsSampled=nEAsSampled, Pixels=nPixels, Neonatals=nPerCon)

smoothRiskCounty = aggResultsKenya$aggregationResultslcpb$countyMatrices$Z
fineScaleRiskCounty = aggResultsKenya$aggregationResultsLCPb$countyMatrices$Z
fineScalePrevalenceCounty = aggResultsKenya$aggregationResultsLCPB$countyMatrices$Z
ihmeCounty = aggResultsKenya$aggregationResultsIHME$countyMatrices$Z

smoothRiskCounty = aggResultsKenya$areaPop$aggregationResults$pSmoothRisk
fineScaleRiskCounty = aggResultsKenya$areaPop$aggregationResults$pFineScaleRisk
fineScalePrevalenceCounty = aggResultsKenya$areaPop$aggregationResults$pFineScalePrevalence
ihmeCounty = aggResultsKenya$aggregationResultsIHME$countyMatrices$p

prevalenceCIWidths = apply(fineScalePrevalenceCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
countyResults = data.frame(Area="Nairobi", SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=sum(nEAsTotal), EAsSampled=sum(nEAsSampled), Pixels=sum(nPixels), Neonatals=sum(nPerCon))

combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
combinedResults$Area = as.character(combinedResults$Area)
combinedResults$Area[1:4] = as.character(1:4)

xtable(combinedResults, 
       digits= c(0, 0, 2, 2, 2, 2, 0, 0, 0, 0), 
       display=c("d", "s", "e", "e", "e", "e", "d", "d", "d", "d"))

# pct larger:
(constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk
mean((constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk)

# Tests ----
popMatSimpleAdjustedW2 = popGridAdjusted[popGridAdjusted$area == "Wajir",]
popMatSimpleAdjustedW2$area = popMatSimpleAdjustedW2$area
easpaW2 = easpaW # easpaW is only adjusted to have 25 neonatals/EA in the simulation code itself
simDatMulti2 = generateSimDataSetsLCPB(nsim=nsim, adjustedPopMat=popMatSimpleAdjustedW2, logisticApproximation=FALSE, 
                                       dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                       simPopOnly=TRUE, returnEAinfo=FALSE, seed=123, inla.seed=123L, verbose=TRUE, 
                                       easpa=easpaW2, popMat=popMatW, constituencyPop=poppsubW)

nTest = 100
test = modLCPB(matrix(uDraws[,1:nTest], ncol=nTest), sigmaEpsilonDraws[1:nTest], easpaW, popMatW, 
               popMatSimpleAdjustedW, doLCPb=TRUE, doIHMEModel=TRUE, 
               constituencyPop=poppsubW, ensureAtLeast1PerConstituency=TRUE, 
               logisticApproximation=FALSE, verbose=TRUE, 
               fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1)

simDatMulti3 = generateSimDataSetsLCPB(nsim=nsim, adjustedPopMat=popMatSimpleAdjustedW, fixPopPerEA=57, logisticApproximation=FALSE, 
                                       dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                       simPopOnly=TRUE, returnEAinfo=FALSE, seed=1, inla.seed=1L, verbose=TRUE, 
                                       easpa=easpaW, popMat=popMatW, constituencyPop=poppsubW, stopOnFrameMismatch=FALSE)

testnsim = 10000

prevalenceSamples = simDatMulti3$subareaPop$aggregationResults$pFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = simDatMulti3$subareaPop$aggregationResults$pFineScaleRisk[,1:testnsim]
smoothRiskSamples = simDatMulti3$subareaPop$aggregationResults$pSmoothRisk[,1:testnsim]

prevalenceSamples = simDatMulti2$subareaPop$aggregationResults$pFineScalePrevalence[,1:testnsim]
fineScaleRiskSamples = simDatMulti2$subareaPop$aggregationResults$pFineScaleRisk[,1:testnsim]
smoothRiskSamples = simDatMulti2$subareaPop$aggregationResults$pSmoothRisk[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths)
#                             [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
# prevalenceCIWidths    0.06632373 0.06385187 0.03664516 0.07031925 0.06450509 0.06157220
# fineScaleRiskCIWidths 0.06437442 0.06254643 0.03585227 0.07003291 0.06347952 0.06055693
# smoothRiskCIWidths    0.06326411 0.06199913 0.03558563 0.06956095 0.06322627 0.05932780

# 95% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths)
#                            [,1]       [,2]       [,3]      [,4]       [,5]       [,6]
# prevalenceCIWidths    0.1038467 0.10096725 0.05808698 0.1134532 0.10217768 0.09642105
# fineScaleRiskCIWidths 0.1018004 0.09832631 0.05633098 0.1129878 0.10123629 0.09357837
# smoothRiskCIWidths    0.1002967 0.09726498 0.05627309 0.1126969 0.09972347 0.09276377

# pct larger:
(prevalenceCIWidths - smoothRiskCIWidths) / smoothRiskCIWidths
mean((prevalenceCIWidths - smoothRiskCIWidths) / smoothRiskCIWidths)
# 0.001829163

# means
prevalenceMeans = rowMeans(prevalenceSamples)
fineScaleRiskMeans = rowMeans(fineScaleRiskSamples)
smoothRiskMeans = rowMeans(smoothRiskSamples)
rbind(prevalenceMeans, fineScaleRiskMeans, smoothRiskMeans)
#                          [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
# prevalenceMeans    0.06615164 0.06444411 0.03602050 0.06635717 0.06439644 0.06137775
# fineScaleRiskMeans 0.06610529 0.06436843 0.03599965 0.06635841 0.06436133 0.06140232
# smoothRiskMeans    0.06609431 0.06444141 0.03602750 0.06643923 0.06431593 0.06142097

# SDs
means = smoothRiskMeans
prevalenceErr = sweep(prevalenceSamples, 1, means, "-")
fineScaleRiskErr = sweep(fineScaleRiskSamples, 1, means, "-")
smoothRiskErr = sweep(smoothRiskSamples, 1, means, "-")
prevalenceSDs = rowMeans(prevalenceErr^2)
fineScaleRiskSDs = rowMeans(fineScaleRiskErr^2)
smoothRiskSDs = rowMeans(smoothRiskErr^2)

rbind(prevalenceSDs, fineScaleRiskSDs, smoothRiskSDs)
#                          [,1]         [,2]         [,3]         [,4]         [,5]         [,6]
# prevalenceSDs    0.0007368619 0.0006832285 0.0002270314 0.0008593671 0.0007034295 0.0006356289
# fineScaleRiskSDs 0.0007050809 0.0006552234 0.0002188096 0.0008433078 0.0006901962 0.0006129011
# smoothRiskSDs    0.0006789629 0.0006376858 0.0002157206 0.0008341588 0.0006807164 0.0005927873

testSamples = simDatMulti$aggregationResultsLcpb$constituencyMatrices$p[poppsubKenya$area=="Wajir",]
testCIWidths = apply(testSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})

sum(popMatSimpleAdjusted$pop[popMatSimpleAdjusted$constituency == "Wajir South" & popMatSimpleAdjusted$urban])

temp = popMatSimpleAdjusted
temp$subarea = as.character(temp$subarea)
temp = temp[temp$area == "Wajir",]
aggregate(temp$pop, by=list(admin2=temp$subarea, urban=temp$urban), FUN=sum, drop=FALSE)
#         admin2 urban         x
# 1        Eldas FALSE 1900.89072
# 2       Tarbaj FALSE 2809.62441
# 3   Wajir East FALSE 1134.57550
# 4  Wajir North FALSE 3575.88109
# 5  Wajir South FALSE 4663.99237
# 6   Wajir West FALSE 2215.03592
# 7        Eldas  TRUE         NA
# 8       Tarbaj  TRUE  104.86047
# 9   Wajir East  TRUE 3365.35398
# 10 Wajir North  TRUE   49.45797
# 11 Wajir South  TRUE  267.17647
# 12  Wajir West  TRUE  288.15111

test = easpsub
test$popUrb = 25 * test$meanUrbanEAs
test$popRur = 25 * test$meanRuralEAs
test$popTotal = 25 * test$meanTotalEAs
test
#     Constituency County     popUrb   popRur popTotal meanUrbanEAs meanRuralEAs meanTotalEAs
# 35         Eldas  Wajir    0.00000 1900.891 1900.891     0.000000     76.03563     76.03563
# 241       Tarbaj  Wajir  104.86047 2809.624 2914.485     4.194419    112.38498    116.57940
# 264   Wajir East  Wajir 3365.35398 1134.575 4499.929   134.614159     45.38302    179.99718
# 265  Wajir North  Wajir   49.45797 3575.881 3625.339     1.978319    143.03524    145.01356
# 266  Wajir South  Wajir  267.17647 4663.992 4931.169    10.687059    186.55969    197.24675
# 267   Wajir West  Wajir  288.15111 2215.036 2503.187    11.526044     88.60144    100.12748


# now make the actual table of SDs
prevalenceMeans = rowMeans(prevalenceSamples)
fineScaleRiskMeans = rowMeans(fineScaleRiskSamples)
smoothRiskMeans = rowMeans(smoothRiskSamples)
rbind(prevalenceMeans, fineScaleRiskMeans, smoothRiskMeans)

means = smoothRiskMeans
prevalenceErr = sweep(prevalenceSamples, 1, means, "-")
fineScaleRiskErr = sweep(fineScaleRiskSamples, 1, means, "-")
smoothRiskErr = sweep(smoothRiskSamples, 1, means, "-")
prevalenceSDs = sqrt(rowMeans(prevalenceErr^2))
fineScaleRiskSDs = sqrt(rowMeans(fineScaleRiskErr^2))
smoothRiskSDs = sqrt(rowMeans(smoothRiskErr^2))
countyMean = mean(simDatMulti$areaPop$aggregationResults$pSmoothRisk[wajirCountyI,])
smoothRiskCountySamples = simDatMulti$areaPop$aggregationResults$pSmoothRisk[wajirCountyI,]
riskCountySamples = simDatMulti$areaPop$aggregationResults$pFineScaleRisk[wajirCountyI,]
prevalenceCountySamples = simDatMulti$areaPop$aggregationResults$pFineScalePrevalence[wajirCountyI,]
countySDs = data.frame(Area="Wajir", 
                       SmoothRisk=sqrt(mean((smoothRiskCountySamples-countyMean)^2)), 
                       Risk=sqrt(mean((riskCountySamples-countyMean)^2)), 
                       Prevalence=sqrt(mean((prevalenceCountySamples-countyMean)^2)), 
                       TotalEAs=sum(nEAsTotal))
constituencySDs = data.frame(Area=constituenciesW, 
                             SmoothRisk=smoothRiskSDs, 
                             Risk=fineScaleRiskSDs, 
                             Prevalence=prevalenceSDs, 
                             TotalEAs=nEAsTotal)

combinedResults = rbind(constituencySDs, countySDs)
combinedResults
#          Area SmoothRisk       Risk Prevalence TotalEAs
# 1       Eldas 0.01823582 0.01861409 0.01956040       82
# 2      Tarbaj 0.01813746 0.01826255 0.01882949      113
# 3  Wajir East 0.01029339 0.01072146 0.01090806      168
# 4 Wajir North 0.01853076 0.01917493 0.02003795      148
# 5 Wajir South 0.01649641 0.01672384 0.01687323      194
# 6  Wajir West 0.01702729 0.01778895 0.01796117      110
# 7       Wajir 0.01417645 0.01419317 0.01429907      815

xtable(combinedResults, digits=3)
# \begin{table}[ht]
# \centering
# \begin{tabular}{rlrrrr}
# \hline
# & Area & SmoothRisk & Risk & Prevalence & TotalEAs \\ 
# \hline
# 1 & Eldas & 0.018 & 0.019 & 0.020 &   82 \\ 
# 2 & Tarbaj & 0.018 & 0.018 & 0.019 &  113 \\ 
# 3 & Wajir East & 0.010 & 0.011 & 0.011 &  168 \\ 
# 4 & Wajir North & 0.019 & 0.019 & 0.020 &  148 \\ 
# 5 & Wajir South & 0.016 & 0.017 & 0.017 &  194 \\ 
# 6 & Wajir West & 0.017 & 0.018 & 0.018 &  110 \\ 
# 7 & Wajir & 0.014 & 0.014 & 0.014 &  815 \\ 
# \hline
# \end{tabular}
# \end{table}
xtable(combinedResults * 100, digits=2)

