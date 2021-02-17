# In this script, we look at a simple example of a real area, attempting 
# to estimate NMR accounting for EA level random variation

## do the setup
longRangeWajir = c(39, 41)
latRangeWajir = c(0.25, 3.6)
constituenciesW = poppcon$Constituency[poppcon$County=="Wajir"]
offsets = matrix(0, nrow=6, ncol=2)
offsets[1,2] = .1 # shift label for Eldas slightly higher
offsets[6,1] = .15 # shift label for Wajir West slightly farther east
easpcon = meanEAsPerCon()
easpcon = easpcon[easpcon$County=="Wajir",]
popGridWajir = popGrid[popGrid$admin1=="Wajir",]

# plotMapDat(mapDat=adm0, lwd=.5, new=TRUE)
# plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], border="green")

# plot population density (1km x 1km)
popCols=makeBlueSequentialColors(64)
ruralCols=makeGreenSequentialColors(64)
urbanCols = popCols
riskCols = makeRedBlueDivergingColors(64, rev=TRUE)

# generate the fine scale pop density grid
popGridFine = makeInterpPopGrid(kmRes=1, mean.neighbor=500, delta=.05)
popGridFine = popGridFine[popGridFine$admin1 == "Wajir",]

# normalize to have the correct population within the county
popGridFine$popOrig = popGridFine$popOrig * (poppc$popTotal[poppc$County=="Wajir"] / sum(popGridFine$popOrig))
popRange = range(popGridFine$popOrig)

## done with setup. Make plots

# plot population density on fine grid
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
quilt.plot(popGridFine$lon, popGridFine$lat, log(popGridFine$popOrig), 
           col=popCols, nx=230, ny=380, add.legend = FALSE, add=TRUE, 
           zlim=log(popRange))
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# now plot the urban pixels
png(file="figures/simpleExample/wajirUrban.png", width=400, height=500)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=FALSE, legend.mar=0, addColorBar=FALSE)
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
popInWajir = popGrid$admin1 == "Wajir"
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
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
popInWajir = popGrid$admin1 == "Wajir"
quilt.plot(popGrid$lon[popInWajir], popGrid$lat[popInWajir], popGrid$urban[popInWajir], col=c(rgb(0, 0, 0, 0), "blue"), nx=60, ny=100, add.legend = FALSE, add=TRUE)
offsets = matrix(0, nrow=6, ncol=2)
offsets[1,2] = .1 # shift label for Eldas slightly higher
offsets[6,1] = .15 # shift label for Wajir West slightly farther east
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# now plot urban fraction as a function of constituency
pdf(file="figures/simpleExample/wajirUrbanFrac.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
urbanFraction = (poppcon$popUrb/poppcon$popTotal)[poppcon$County=="Wajir"]
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
urbanFraction = (poppcon$popUrb/poppcon$popTotal)[poppcon$County=="Wajir"]
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

# plot expected number of urban EAs per constituency

pdf(file="figures/simpleExample/wajirUrbanEAs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
meanUrbanEAs = easpcon$meanUrbanEAs
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
meanRuralEAs = easpcon$meanRuralEAs
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
meanRuralEAs = easpcon$meanRuralEAs
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

##### simulate EA locations for example
simDat = generateSimDataSetsLCPB(nsim=1, seed=123, 
                                 dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/")
simulatedEAs = simDat$simulatedEAs
eaSamples = simDat$simulatedEAs$eaSamples
SRSDat = simDat$SRSDat # or overSampDat?
eaDat = SRSDat$eaDat
eaDat = eaDat[eaDat$admin1 == "Wajir",]
eaSamples = eaSamples[popGrid$admin1=="Wajir",1]

# plot EA locations
pdf(file="figures/simpleExample/wajirSimEALocs.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
points(eaDat$lon, eaDat$lat, pch=19, cex=.2, col=rgb(1, 0, 0, .2))
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# same but with counts per pixel
pdf(file="figures/simpleExample/wajirSimEAsPerPixel.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
eaRange = c(1, 60)
quilt.plot(popGrid$lon[popGrid$admin1=="Wajir"], 
           popGrid$lat[popGrid$admin1=="Wajir"], 
           log(eaSamples), 
           col=popCols[-(1:5)], nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=log(eaRange))
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
eaTicks = c(1, 5, 10, 50)
eaTickLabels = as.character(eaTicks)
image.plot(zlim=log(eaRange), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=popCols[-(1:5)], add=TRUE, axis.args=list(at=log(eaTicks), labels=eaTickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# check pixel constituencies
pdf(file="figures/simpleExample/wajirPixelConstituencyCheck.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
tempVals = factor(popGrid$admin2[popGrid$admin1=="Wajir"])
tempVals = as.numeric(tempVals)
quilt.plot(popGrid$lon[popGrid$admin1=="Wajir"], 
           popGrid$lat[popGrid$admin1=="Wajir"], 
           tempVals, 
           col=rainbow(6), nx=45, ny=60, add.legend = TRUE, add=TRUE, 
           FUN=max)
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot the number of sampled EAs per constituency
out = aggregate(eaSamples, by=list(constituency=popGridWajir$admin2), FUN=sum)
nEAs = out$x
nEAsTotal = nEAs
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

urbanPixels = popGridWajir$urban
eaSamplesUrban = eaSamples
eaSamplesRural = eaSamples
eaSamplesUrban[!urbanPixels] = 0
eaSamplesRural[urbanPixels] = 0

# plot the number of sampled Urban EAs per constituency
out = aggregate(eaSamplesUrban, by=list(constituency=popGridWajir$admin2), FUN=sum)
nEAs = out$x
nEAsUrban = nEAs
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
out = aggregate(eaSamplesRural, by=list(constituency=popGridWajir$admin2), FUN=sum)
nEAs = out$x
nEAsRural = nEAs
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

# plot expected risk surface for each pixel
expectedRisk = SRSDat$aggregatedPop$pixelMatriceslcpb$p[,1]
expectedRisk = expectedRisk[popGrid$admin1 == "Wajir"]

pdf(file="figures/simpleExample/wajirSimExpectedRisk.pdf", width=4, height=5)
par(mar=c(4.1, 4.1, 1.1, 4.5))
plotMapDat(mapDat=adm1[adm1@data$NAME_1=="Wajir",], new=TRUE, 
           kenyaLonRange = longRangeWajir, kenyaLatRange = latRangeWajir, 
           leaveRoomForLegend=TRUE, addColorBar=FALSE, legend.mar=5)
quilt.plot(popGrid$lon[popGrid$admin1=="Wajir"], 
           popGrid$lat[popGrid$admin1=="Wajir"], 
           logit(expectedRisk), 
           col=riskCols, nx=45, ny=60, add.legend = FALSE, add=TRUE, 
           zlim=logit(c(.03, .2)))
plotMapDat(mapDat=adm2[adm2@data$COUNTY_NAM=="Wajir",], lwd=.5)
ticks = c(.03, .08, .13, .18)
tickLabels = as.character(ticks)
image.plot(zlim=logit(c(.03, .2)), nlevel=length(popCols), legend.only=TRUE, horizontal=FALSE,
           col=riskCols, add=TRUE, axis.args=list(at=logit(ticks), labels=tickLabels, cex.axis=1, tck=-.7, hadj=-.1), 
           legend.cex=.5, legend.width=1)
addMapLabels(constituenciesW, mapDat=adm2, offsets=offsets, cex=.4)
dev.off()

# plot EA level prevalence versus risk:
easInWajir = eaDat$admin1 == "Wajir"
pdf(file="figures/simpleExample/wajirSimlcpbVLCPB.pdf", width=5, height=5)
plot(eaDat$plcpb[easInWajir], eaDat$pLCPB[easInWajir], pch=19, cex=.1, col="blue", type="n", 
     xlab="Risk (S)", ylab="Prevalence (SCP)")
abline(a=0, b=1)
points(eaDat$plcpb[easInWajir], eaDat$pLCPB[easInWajir], pch=19, cex=.1, col="blue")
dev.off()

# plot constituency level prevalence versus risk:
sortI = sort(as.character(poppcon$Constituency), index.return=TRUE)$ix
theseCounties = poppcon$County[sortI]
constituenciesInWajir = theseCounties == "Wajir"
pdf(file="figures/simpleExample/wajirSimlcpbVLCPBConstituency.pdf", width=5, height=5)
plot(SRSDat$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1], 
     SRSDat$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1], 
     pch=19, cex=1, col="blue", type="n", 
     xlab="Risk (S)", ylab="Prevalence (SCP)", xlim=c(.065, .15), ylim=c(.065, .15))
abline(a=0, b=1)
points(SRSDat$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1], 
       SRSDat$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1], 
       pch=19, cex=1, col="blue")
text(SRSDat$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1], 
     SRSDat$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1], 
     constituenciesW, cex=.3, pos=4)
dev.off()

# calulcate percent different between risk and prevalence
lcpb = SRSDat$aggregatedPop$aggregatedResultslcpb$constituencyMatrices$p[constituenciesInWajir,1]
LCPB = SRSDat$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesInWajir,1]
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

##### 

