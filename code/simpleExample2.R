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
summerScale = makeBlueGreenYellowSequentialColors(64)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
  # popGridFine$pop = popGridFine$pop * (poppc$popTotal[poppcarea=="Wajir"] / sum(popGridFine$pop))
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
  addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
  addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.7)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
prevalenceSamples = simDat$simulatedEAs$aggregatedPop$pixelMatricesLCPB$p[,1]

# plot EA locations ----
pdf(file="figures/simpleExample/wajirSimEALocs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col="blue")
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot EA locations with color for fine scale prevalence
pRangeMod = range(c(eaDat$pLCPB[eaDat$pLCPB != 0], eaDat$pLCPb, eaDat$plcpb))
ticks = c(.01, .05, .1, .2, .3, .4, .5)
pdf(file="figures/simpleExample/wajirSimEALocsPrevalence.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pLCPB, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot EA locations with color for fine scale risk
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRisk.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pLCPb, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot EA locations with color for fine scale risk on its own scale
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskSelfScaled.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pLCPb, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot EA locations with color for smooth risk
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRisk.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pLcpb, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot EA locations with color for smooth risk on its own color scale
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRiskSelfScaled.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDat$lon, eaDat$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDat$lon, eaDat$lat, eaDat$pLcpb, pch=19, cex=.2, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="none", ticks=ticks)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.5)
dev.off()

# plot EA locations with color for fine scale prevalence
pRangeMod = range(c(eaDatNW$pLCPB[eaDatNW$pLCPB != 0], eaDatNW$pLCPb, eaDatNW$plcpb))
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
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pLCPB, pch=21, cex=1, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, colorName="bg", 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# same but residuals
pdf(file="figures/simpleExample/wajirSimEALocsPrevalenceNWajirResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = eaDatNW$pLCPB - eaDatNW$plcpb
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, prevalenceResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# same but percent residuals
pdf(file="figures/simpleExample/wajirSimEALocsPrevalenceNWajirPctResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = 100 * (eaDatNW$pLCPB - eaDatNW$plcpb) / eaDatNW$plcpb
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, prevalenceResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# plot EA locations with color for fine scale risk
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pLCPb, pch=21, cex=1, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# same but different color scale
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajir2.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(eaDatNW$pLCPB), center = median(eaDatNW$pLCPB))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pLCPb, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=pRangeMod, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# same but residuals
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajirResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = eaDatNW$pLCPB - eaDatNW$plcpb
riskResiduals = eaDatNW$pLCPb - eaDatNW$plcpb
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, riskResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# same but percent residuals
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskNWajirPctResiduals.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))

plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
prevalenceResiduals = 100 * (eaDatNW$pLCPB - eaDatNW$plcpb) / eaDatNW$plcpb
riskResiduals = 100 * (eaDatNW$pLCPb - eaDatNW$plcpb) / eaDatNW$plcpb
# tempRiskCols = makePurpleRedDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
tempRiskCols = makeRedBlueDivergingColors(64, rev=TRUE, valRange = range(prevalenceResiduals), center = 0)
plotWithColor(eaDatNW$lon, eaDatNW$lat, riskResiduals, pch=21, cex=1, 
              colScale=tempRiskCols, new=FALSE, zlim=range(prevalenceResiduals), 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# plot EA locations with color for fine scale risk on its own scale
pdf(file="figures/simpleExample/wajirSimEALocsFineScaleRiskSelfScaledNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pLCPb, pch=21, cex=.8, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# plot EA locations with color for smooth risk
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRiskNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pLcpb, pch=21, cex=.8, 
              colScale=riskCols, new=FALSE, zlim=pRangeMod, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

# plot EA locations with color for smooth risk on its own color scale
pdf(file="figures/simpleExample/wajirSimEALocsSmoothRiskSelfScaledNWajir.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# points(eaDatNW$lon, eaDatNW$lat, pch=19, cex=.1, col=rgb(1, 0, 0, .2))
plotWithColor(eaDatNW$lon, eaDatNW$lat, eaDatNW$pLcpb, pch=21, cex=.8, 
              colScale=riskCols, new=FALSE, 
              scaleFun=logit, scaleFunInverse=expit, forceColorsInRange=TRUE, 
              legend.cex=.5, legend.width=1, legend.mar=5, colorName="bg", 
              ordering="increasing", ticks=ticks, tickLabels=tickLabels, 
              col=rgb(.5, .5, .5))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot smooth risk
expectedRisk = simDat$simulatedEAs$aggregatedPop$pixelMatriceslcpb$p[,1]
pRangeMod = range(c(eaDat$pLCPB[eaDat$pLCPB != 0], eaDat$pLCPb, eaDat$plcpb))

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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot smooth risk for Wajir North
expectedRisk = simDat$simulatedEAs$aggregatedPop$pixelMatriceslcpb$p[,1]
expectedRisk = expectedRisk[popGridWajir$subarea == "Wajir North"]
# pRangeMod = range(c(eaDat$pLCPB[eaDat$pLCPB != 0], eaDat$pLCPb, eaDat$plcpb))
pRangeMod = range(c(eaDatNW$pLCPB[eaDatNW$pLCPB != 0], eaDatNW$pLCPb, eaDatNW$plcpb))

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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
dev.off()

pdf(file="figures/simpleExample/wajirSimSmoothRiskNWajirSelfScaledLinear.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$subarea=="Wajir North"], 
           popGrid$lat[popGrid$subarea=="Wajir North"], 
           expectedRisk, 
           col=riskCols, nx=29, ny=24, add.legend = FALSE, add=TRUE, 
           zlim=range(expectedRisk))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
ticks = c(.02, .03, .04, .05, .06, .07)
tickLabels = as.character(ticks)
image.plot(zlim=range(expectedRisk), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=ticks, labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
  addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
# addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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

# plot pop density with same scale
if(FALSE) {
  
  pdf(file="figures/simpleExample/wajirPopDensityNormNW.pdf", width=5, height=3.9)
  par(mar=c(4.1, 4.1, 1.1, 4.5))
  plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
             kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
             leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
  image.plot(zlim=log(weightRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add=TRUE, axis.args=list(at=log(weightTicks), labels=weightTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
             legend.cex=.5, legend.width=1)
  quilt.plot(popMatSimpleAdjustedNW$lon, popMatSimpleAdjustedNW$lat, log(popWeights), 
             col=popCols, nx=29, ny=24, add.legend = FALSE, add=TRUE, 
             zlim=log(weightRange))
  plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
  addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
  dev.off()
}

# same but with number of assigned people per EA

nEAsNW = sum(eaSamples[popGridWajir$subarea == "Wajir North"])
weightsEA = rep(1/nEAsNW, nEAsNW)
pdf(file="figures/simpleExample/wajirSimNPerEANormNW.pdf", width=5, height=3.9)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeNorthWajir, kenyaLatRange = latRangeNorthWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotWithColor(eaDatNW$lon, eaDatNW$lat, weightsEA, 
              colScale=popCols, add.legend = FALSE, new=FALSE, 
              zlim=weightRange, scaleFun=log, scaleFunInverse=exp, 
              legend.cex=.5, legend.width=1, legend.mar=5, 
              pch=19, cex=1, ticks=weightTicks)
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
# NTicks = c(35, 100, 300, 800, 2400)
# NTickLabels = as.character(NTicks)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.8)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# same but with number of assigned people per pixel
pdf(file="figures/simpleExample/wajirSimNPerPixelNorm2.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$area=="Wajir"], 
           popGrid$lat[popGrid$area=="Wajir"], 
           log(NSamplesNW), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(weightRange2))
plotMapDat(mapDat=adm2[adm2@data$NAME_1=="Wajir",], lwd=.5)
NTicks = c(35, 100, 300, 800, 2400)
NTickLabels = as.character(NTicks)
image.plot(zlim=log(weightRange2), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(weightTicks2), labels=weightTickLabels2, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
# addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
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
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot prevalence and risk ----
easInWajir = eaDat$area == "Wajir"
pdf(file="figures/simpleExample/wajirSimlcpbVLCPB.pdf", width=5, height=5)
plot(eaDat$plcpb[easInWajir], eaDat$pLCPB[easInWajir], pch=19, cex=.1, col="blue", type="n", 
     xlab="Smooth Risk", ylab="Fine Scale Prevalence")
abline(a=0, b=1)
points(eaDat$plcpb[easInWajir], eaDat$pLCPB[easInWajir], pch=19, cex=.1, col="blue")
dev.off()

# plot constituency level prevalence versus risk:
sortI = sort(as.character(poppsubKenya$Constituency), index.return=TRUE)$ix
theseCounties = poppsubKenyaarea[sortI]
constituenciesInWajir = theseCounties == "Wajir"
pdf(file="figures/simpleExample/wajirSimlcpbVLCPBConstituency.pdf", width=5, height=5)
plot(simDat$simulatedEAs$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1], 
     simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1], 
     pch=19, cex=1, col="blue", type="n", 
     xlab="Smooth Risk", ylab="Fine Scale Prevalence", xlim=c(.03, .1), ylim=c(.03, .1))
abline(a=0, b=1)
points(simDat$simulatedEAs$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1], 
       simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1], 
       pch=19, cex=1, col="blue")
text(simDat$simulatedEAs$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1], 
     simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1], 
     constituenciesW, cex=.55, pos=4)
dev.off()

# Tabulate prevalence versus risk ----
# calulcate percent different between risk and prevalence
lcpb = simDat$simulatedEAs$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[,1]
LCPb = simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPb$constituencyMatrices$p[,1]
LCPB = simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[,1]
constituencyResults = data.frame(Area=constituenciesW, SmoothRisk=lcpb, Risk=LCPb, Prevalence=LCPB, TotalEAs=nEAsTotal)
wajirCountyI = sort(unique(poppcarea)) == "Wajir"
lcpbCounty = simDat$simulatedEAs$aggregatedPop$aggregatedResultslcpb$countyMatrices$p[,1]
LCPbCounty = simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPb$countyMatrices$p[,1]
LCPBCounty = simDat$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$countyMatrices$p[,1]
countyResults = data.frame(Area="Wajir", SmoothRisk=lcpbCounty, Risk=LCPbCounty, Prevalence=LCPBCounty, TotalEAs=sum(nEAsTotal))
combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
#          Area SmoothRisk       Risk Prevalence  TotalEAs
# 1       Eldas 0.06557801 0.06354172 0.06342857       70
# 2      Tarbaj 0.06531854 0.06507623 0.06754098      122
# 3  Wajir East 0.04102762 0.04360798 0.04219780      182
# 4 Wajir North 0.06117243 0.05823436 0.06537931      145
# 5 Wajir South 0.07255666 0.07674612 0.08194595      185
# 6  Wajir West 0.07283498 0.06545748 0.06414414      111
# 7       Wajir 0.06191546 0.06163395 0.06395092      815

percentDifference = (LCPB - lcpb) / lcpb
print(data.frame(Constituency=constituenciesW, pctDiff=percentDifference, urbanEAs=nEAsUrban, ruralEAs=nEAsRural, totalEAs=nEAsTotal))
#   Constituency     pctDiff urbanEAs ruralEAs totalEAs
# 1        Eldas -0.05107582        0       79       79
# 2       Tarbaj  0.05118402        4      121      125
# 3   Wajir East  0.01746166      155       43      198
# 4  Wajir North -0.04886100        4      123      127
# 5  Wajir South -0.03683734        0      198      198
# 6   Wajir West  0.04094348        0       88       88

#  Constituency pctDiff urbanEAs ruralEAs totalEAs
#         Eldas    -5.1        0       79       79
#        Tarbaj     5.1        4      121      125
#    Wajir East     1.7      155       43      198
#   Wajir North    -4.9        4      123      127
#   Wajir South    -3.7        0      198      198
#    Wajir West     4.1        0       88       88

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
simDatMulti = generateSimDataSetsLCPB(nsim=nsim, adjustedPopMat=popMatSimpleAdjustedW, fixPopPerEA=25, logisticApproximation=FALSE, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      simPopOnly=TRUE, returnEAinfo=FALSE, seed=1, inla.seed=1L, verbose=TRUE, 
                                      easpa=easpaW, popMat=popMatW, constituencyPop=poppsubW)
testnsim = 10000
prevalenceSamples = simDatMulti$aggregatedResultsLCPB$constituencyMatrices$Z[,1:testnsim]
fineScaleRiskSamples = simDatMulti$aggregatedResultsLCPb$constituencyMatrices$Z[,1:testnsim]
smoothRiskSamples = simDatMulti$aggregatedResultslcpb$constituencyMatrices$Z[,1:testnsim]

prevalenceSamples = simDatMulti$aggregatedResultsLCPB$constituencyMatrices$p[,1:testnsim]
fineScaleRiskSamples = simDatMulti$aggregatedResultsLCPb$constituencyMatrices$p[,1:testnsim]
smoothRiskSamples = simDatMulti$aggregatedResultslcpb$constituencyMatrices$p[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
# constituencyResults = rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths)
constituencyResults = data.frame(Area=constituenciesW, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, TotalEAs=nEAsTotal)

lcpbCounty = simDatMulti$aggregatedResultslcpb$countyMatrices$Z
LCPbCounty = simDatMulti$aggregatedResultsLCPb$countyMatrices$Z
LCPBCounty = simDatMulti$aggregatedResultsLCPB$countyMatrices$Z

lcpbCounty = simDatMulti$aggregatedResultslcpb$countyMatrices$p
LCPbCounty = simDatMulti$aggregatedResultsLCPb$countyMatrices$p
LCPBCounty = simDatMulti$aggregatedResultsLCPB$countyMatrices$p

prevalenceCIWidths = apply(LCPBCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(LCPbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(lcpbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
countyResults = data.frame(Area="Wajir", SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, TotalEAs=sum(nEAsTotal))

combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
#         Area SmoothRisk       Risk Prevalence TotalEAs
# 1       Eldas 0.06326411 0.06437442 0.06632373       62
# 2      Tarbaj 0.06199913 0.06254643 0.06385187      122
# 3  Wajir East 0.03558563 0.03585227 0.03664516      173
# 4 Wajir North 0.06956095 0.07003291 0.07031925      138
# 5 Wajir South 0.06322627 0.06347952 0.06450509      215
# 6  Wajir West 0.05932780 0.06055693 0.06157220      105
# 7       Wajir 0.05270074 0.05252068 0.05281472      815
combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
# pct larger:
(constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk
# 0.04836255 0.02988332 0.02977417 0.01090124 0.02022610 0.03783050
mean((constituencyResults$Prevalence - constituencyResults$SmoothRisk) / constituencyResults$SmoothRisk)
# [1] 0.02949631

# 95% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))})
rbind(prevalenceCIWidths, fineScaleRiskCIWidths, smoothRiskCIWidths)
# [,1]       [,2]       [,3]      [,4]       [,5]       [,6]
# prevalenceCIWidths    0.1038467 0.10096725 0.05808698 0.1134532 0.10217768 0.09642105
# fineScaleRiskCIWidths 0.1018004 0.09832631 0.05633098 0.1129878 0.10123629 0.09357837
# smoothRiskCIWidths    0.1002967 0.09726498 0.05627309 0.1126969 0.09972347 0.09276377

mean((prevalenceCIWidths - smoothRiskCIWidths) / smoothRiskCIWidths)

# means
prevalenceMeans = rowMeans(prevalenceSamples)
fineScaleRiskMeans = rowMeans(fineScaleRiskSamples)
smoothRiskMeans = rowMeans(smoothRiskSamples)
rbind(prevalenceMeans, fineScaleRiskMeans, smoothRiskMeans)
# [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
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
# [,1]         [,2]         [,3]         [,4]         [,5]         [,6]
# prevalenceSDs    0.0007368619 0.0006832285 0.0002270314 0.0008593671 0.0007034295 0.0006356289
# fineScaleRiskSDs 0.0007050809 0.0006552234 0.0002188096 0.0008433078 0.0006901962 0.0006129011
# smoothRiskSDs    0.0006789629 0.0006376858 0.0002157206 0.0008341588 0.0006807164 0.0005927873

# pct difference
(prevalenceSDs - smoothRiskSDs) / smoothRiskSDs
mean((prevalenceSDs - smoothRiskSDs) / smoothRiskSDs)

##### Fit model to SRS cluster data ----
# get the data
dat = simDat$SRSDat$clustDat[[1]]

# fit SPDE cluster level risk model
nSamples = 10000 # 10000 takes about 5 minutes
spdeFit = fitSPDEKenyaDat(dat, nPostSamples=nSamples, popMat=popMatW)

# obtain model output
uDraws = spdeFit$uDraws
sigmaEpsilonDraws = spdeFit$sigmaEpsilonDraws

# apply aggregation models
aggResults = modLCPB(uDraws, sigmaEpsilonDraws, easpaW, popMatW, 
                     popMatSimpleAdjustedW, doLCPb=TRUE, doIHMEModel=TRUE, 
                     constituencyPop=poppsubW, ensureAtLeast1PerConstituency=TRUE, 
                     logisticApproximation=FALSE, verbose=TRUE, 
                     fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1)

# get average number of EAs per area
out = meanEAsPerCon()
out = out[outarea == "Wajir",]
nEAsTotal = out$meanTotalEAs

# get average number of neonatals per area
nPerCon = aggResults$aggregatedResultslcpb$constituencyMatrices$N[,1]

# calculate number of pixels per area (for understanding IHME model variance)
out = aggregate(popMatSimpleAdjustedW$pop, by=list(pixels=popMatSimpleAdjustedW$subarea), FUN=length)
nPixels = out$x

# calculate number of EAs sampled per area (for understanding spatial variance)
out = aggregate(dat$N[dat$area == "Wajir"], by=list(dat$subarea[dat$area == "Wajir"]), FUN=length)
nEAsSampled = out$x

# gather aggregation model output
prevalenceSamples = aggResults$aggregatedResultsLCPB$constituencyMatrices$Z[,1:testnsim]
fineScaleRiskSamples = aggResults$aggregatedResultsLCPb$constituencyMatrices$Z[,1:testnsim]
smoothRiskSamples = aggResults$aggregatedResultslcpb$constituencyMatrices$Z[,1:testnsim]
ihmeSamples = aggResults$aggregatedResultsIHME$constituencyMatrices$Z[,1:testnsim]

prevalenceSamples = aggResults$aggregatedResultsLCPB$constituencyMatrices$p[,1:testnsim]
fineScaleRiskSamples = aggResults$aggregatedResultsLCPb$constituencyMatrices$p[,1:testnsim]
smoothRiskSamples = aggResults$aggregatedResultslcpb$constituencyMatrices$p[,1:testnsim]
ihmeSamples = aggResults$aggregatedResultsIHME$constituencyMatrices$p[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
constituencyResults = data.frame(Area=constituenciesW, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=nEAsTotal, EAsSampled=nEAsSampled, Pixels=nPixels, Neonatals=nPerCon)

lcpbCounty = aggResults$aggregatedResultslcpb$countyMatrices$Z
LCPbCounty = aggResults$aggregatedResultsLCPb$countyMatrices$Z
LCPBCounty = aggResults$aggregatedResultsLCPB$countyMatrices$Z
ihmeCounty = aggResults$aggregatedResultsIHME$countyMatrices$Z

lcpbCounty = aggResults$aggregatedResultslcpb$countyMatrices$p
LCPbCounty = aggResults$aggregatedResultsLCPb$countyMatrices$p
LCPBCounty = aggResults$aggregatedResultsLCPB$countyMatrices$p
ihmeCounty = aggResults$aggregatedResultsIHME$countyMatrices$p

prevalenceCIWidths = apply(LCPBCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(LCPbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(lcpbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
countyResults = data.frame(Area="Wajir", SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=sum(nEAsTotal), EAsSampled=sum(nEAsSampled), Pixels=sum(nPixels), Neonatals=sum(nPerCon))

combinedResults = rbind(constituencyResults, countyResults)
print(combinedResults)
xtable(combinedResults, digits=3)
xtable(combinedResults, 
       digits= c(0, 0, 3, 3, 3, 3, 0, 0, 0, 0), 
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

# Do the same, but in Nairobi ----
longRangeNairobi = c(39, 41)
latRangeNairobi = c(0.25, 3.6)
constituenciesN = poppsubKenya$subarea[poppsubKenya$area=="Nairobi"]
offsets = matrix(0, nrow=4, ncol=2)
# offsets[1,2] = .1 # shift label for Eldas slightly higher
# offsets[6,1] = .15 # shift label for Wajir West slightly farther east
easpsub = meanEAsPerCon()
easpsub = easpsub[easpsubarea=="Nairobi",]
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

simDatKenya = generateSimDataSetsLCPB(nsim=1, adjustedPopMat=popMatSimpleAdjusted, 
                                      fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                      logisticApproximation=FALSE, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      seed=1, inla.seed=1L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                      easpa=NULL, popMat=NULL, constituencyPop=poppsubKenya)

# get the data and the true population
dat = simDatKenya$SRSDat$clustDat[[1]]
truePrevalenceConstituencyKenya = simDatKenya$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[,1]
truePrevalenceCountyKenya = simDatKenya$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$countyMatrices$p[,1]

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
nPerCon = aggResultsN$aggregatedResultslcpb$constituencyMatrices$N[,1]

# calculate number of pixels per area (for understanding IHME model variance)
out = aggregate(popMatSimpleAdjustedN$pop, by=list(pixels=popMatSimpleAdjustedN$subarea), FUN=length)
nPixels = out$x

# calculate number of EAs sampled per area (for understanding spatial variance)
out = aggregate(dat$N[dat$area == "Nairobi"], by=list(dat$subarea[dat$area == "Nairobi"]), FUN=length)
nEAsSampled = out$x

# gather aggregation model output
prevalenceSamples = aggResultsN$aggregatedResultsLCPB$constituencyMatrices$Z[,1:testnsim]
fineScaleRiskSamples = aggResultsN$aggregatedResultsLCPb$constituencyMatrices$Z[,1:testnsim]
smoothRiskSamples = aggResultsN$aggregatedResultslcpb$constituencyMatrices$Z[,1:testnsim]
ihmeSamples = aggResultsN$aggregatedResultsIHME$constituencyMatrices$Z[,1:testnsim]

prevalenceSamples = aggResultsN$aggregatedResultsLCPB$constituencyMatrices$p[,1:testnsim]
fineScaleRiskSamples = aggResultsN$aggregatedResultsLCPb$constituencyMatrices$p[,1:testnsim]
smoothRiskSamples = aggResultsN$aggregatedResultslcpb$constituencyMatrices$p[,1:testnsim]
ihmeSamples = aggResultsN$aggregatedResultsIHME$constituencyMatrices$p[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
constituencyResults = data.frame(Area=constituenciesN, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=nEAsTotal, EAsSampled=nEAsSampled, Pixels=nPixels, Neonatals=nPerCon)

lcpbCounty = aggResultsN$aggregatedResultslcpb$countyMatrices$Z
LCPbCounty = aggResultsN$aggregatedResultsLCPb$countyMatrices$Z
LCPBCounty = aggResultsN$aggregatedResultsLCPB$countyMatrices$Z
ihmeCounty = aggResultsN$aggregatedResultsIHME$countyMatrices$Z

lcpbCounty = aggResultsN$aggregatedResultslcpb$countyMatrices$p
LCPbCounty = aggResultsN$aggregatedResultsLCPb$countyMatrices$p
LCPBCounty = aggResultsN$aggregatedResultsLCPB$countyMatrices$p
ihmeCounty = aggResultsN$aggregatedResultsIHME$countyMatrices$p

prevalenceCIWidths = apply(LCPBCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(LCPbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(lcpbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
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
nPerCon = rowMeans(aggResultsKenya$aggregatedResultslcpb$constituencyMatrices$N)

# calculate number of pixels per area (for understanding IHME model variance)
out = aggregate(popMatSimpleAdjustedN$pop, by=list(pixels=popMatSimpleAdjustedN$subarea), FUN=length)
nPixels = out$x

# calculate number of EAs sampled per area (for understanding spatial variance)
out = aggregate(dat$N[dat$area == "Nairobi"], by=list(dat$subarea[dat$area == "Nairobi"]), FUN=length)
nEAsSampled = out$x

# gather aggregation model output
prevalenceSamples = aggResultsKenya$aggregatedResultsLCPB$constituencyMatrices$Z[,1:testnsim]
fineScaleRiskSamples = aggResultsKenya$aggregatedResultsLCPb$constituencyMatrices$Z[,1:testnsim]
smoothRiskSamples = aggResultsKenya$aggregatedResultslcpb$constituencyMatrices$Z[,1:testnsim]
ihmeSamples = aggResultsKenya$aggregatedResultsIHME$constituencyMatrices$Z[,1:testnsim]

prevalenceSamples = aggResultsKenya$aggregatedResultsLCPB$constituencyMatrices$p[,1:testnsim]
fineScaleRiskSamples = aggResultsKenya$aggregatedResultsLCPb$constituencyMatrices$p[,1:testnsim]
smoothRiskSamples = aggResultsKenya$aggregatedResultslcpb$constituencyMatrices$p[,1:testnsim]
ihmeSamples = aggResultsKenya$aggregatedResultsIHME$constituencyMatrices$p[,1:testnsim]

# 80% CIs
prevalenceCIWidths = apply(prevalenceSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(fineScaleRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(smoothRiskSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
ihmeCIWidths = apply(ihmeSamples, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
constituencyResults = data.frame(Area=constituenciesN, SmoothRisk=smoothRiskCIWidths, Risk=fineScaleRiskCIWidths, Prevalence=prevalenceCIWidths, IHMERisk=ihmeCIWidths, EAs=nEAsTotal, EAsSampled=nEAsSampled, Pixels=nPixels, Neonatals=nPerCon)

lcpbCounty = aggResultsKenya$aggregatedResultslcpb$countyMatrices$Z
LCPbCounty = aggResultsKenya$aggregatedResultsLCPb$countyMatrices$Z
LCPBCounty = aggResultsKenya$aggregatedResultsLCPB$countyMatrices$Z
ihmeCounty = aggResultsKenya$aggregatedResultsIHME$countyMatrices$Z

lcpbCounty = aggResultsKenya$aggregatedResultslcpb$countyMatrices$p
LCPbCounty = aggResultsKenya$aggregatedResultsLCPb$countyMatrices$p
LCPBCounty = aggResultsKenya$aggregatedResultsLCPB$countyMatrices$p
ihmeCounty = aggResultsKenya$aggregatedResultsIHME$countyMatrices$p

prevalenceCIWidths = apply(LCPBCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
fineScaleRiskCIWidths = apply(LCPbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
smoothRiskCIWidths = apply(lcpbCounty, 1, function(x){diff(quantile(x, probs=c(0.1, 0.9), na.rm=TRUE))})
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

prevalenceSamples = simDatMulti3$aggregatedResultsLCPB$constituencyMatrices$p[,1:testnsim]
fineScaleRiskSamples = simDatMulti3$aggregatedResultsLCPb$constituencyMatrices$p[,1:testnsim]
smoothRiskSamples = simDatMulti3$aggregatedResultslcpb$constituencyMatrices$p[,1:testnsim]

prevalenceSamples = simDatMulti2$aggregatedResultsLCPB$constituencyMatrices$p[,1:testnsim]
fineScaleRiskSamples = simDatMulti2$aggregatedResultsLCPb$constituencyMatrices$p[,1:testnsim]
smoothRiskSamples = simDatMulti2$aggregatedResultslcpb$constituencyMatrices$p[,1:testnsim]

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

testSamples = simDatMulti$aggregatedResultsLcpb$constituencyMatrices$p[poppsubKenya$area=="Wajir",]
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
countyMean = mean(simDatMulti$aggregatedResultslcpb$countyMatrices$p[wajirCountyI,])
smoothRiskCountySamples = simDatMulti$aggregatedResultslcpb$countyMatrices$p[wajirCountyI,]
riskCountySamples = simDatMulti$aggregatedResultsLCPb$countyMatrices$p[wajirCountyI,]
prevalenceCountySamples = simDatMulti$aggregatedResultsLCPB$countyMatrices$p[wajirCountyI,]
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

