# this script is for comparing performance between models
# download 2014 Kenya population density and associated TIF file
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/pop.rda?raw=true")
popFilename = paste0(tempDirectory, "/pop.rda")
if(!file.exists(popFilename)) {
  download.file(githubURL,popFilename)
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

popMatSimple = makePopIntegrationTab(kmRes=5, pop=pop, domainPoly=kenyaPoly, 
                                     eastLim=eastLim, northLim=northLim, 
                                     mapProjection=projKenya, poppa=poppaKenya, 
                                     poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                     areaMapDat=adm1, subareaMapDat=adm2, 
                                     areaPolygonSubsetI=30)

popMatSimpleNeonatal = adjustPopMat(popMatSimple, poppaTarget=poppsubKenyaNeonatal, adjustBy="subarea")
easpaSimple = makeDefaultEASPA()
easpaSimple = easpaSimple[easpaSimple$area == "Nairobi",]
poppsubSimple = poppsubKenya
poppsubSimple = poppsubSimple[poppsubSimple$area == "Nairobi",]
simDatKenya = generateSimDataSetsLCPB2(nsim=1, targetPopMat=popMatSimpleNeonatal, 
                                       popMat=popMatSimple, 
                                      fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                      logisticApproximation=FALSE, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      seed=1, inla.seed=1L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                      easpa=easpaSimple, poppsub=poppsubSimple)

# get the data and the true population
dat = simDatKenya$SRSDat$clustDat[[1]]
dat$y = dat$Z
dat$n = dat$N
dat$admin1 = dat$area
dat$admin2 = dat$subarea
constituenciesN = poppsubSimple$area == "Nairobi"
countyN = sort(unique(poppc$County)) == "Nairobi"
truePrevalenceConstituencyKenya = simDatKenya$simulatedEAs$aggregatedPop$subareaPop$aggregationResults$pFineScalePrevalence
truePrevalenceCountyKenya = simDatKenya$simulatedEAs$aggregatedPop$areaPop$aggregationResults$pFineScalePrevalence

# construct integration grids at different resolutions
# resolutions = c(1, 10, 20, 40, 60, 80)
# deltas = c(.05, .2, .4, .8, 1.2, 1.4)
# meanNeighbors = c(500, rep(50, 5))
resolutions = c(.1, .2, .5, 1, 2, 5, 10, 20)
deltas = c(.005, .01, .025, .05, .1, .1, .2, .4)
meanNeighbors = c(rep(500, 5), rep(50, 3))
popGrids = list()
popGridsAdjusted = list()
for(i in 1:length(resolutions)) {
  print(paste0("Creating integration grid at ", resolutions[i], " km resolution"))
  # thisPopGrid = makePopIntegrationTab(kmRes=resolutions[i], mapProjection=SUMMER::projKenya, 
  #                                     domainPoly=kenyaPoly, 
  #                                     subareaMapDat=adm2, poppsub=poppsub, 
  #                                     eastLim=eastLim, northLim=northLim, 
  #                                     mean.neighbor=meanNeighbors[i], 
  #                                     delta=deltas[i], 
  #                                     areaMapDat=adm1, areaPolygonSubsetI=30) # Nairobi is the 30th one
  # thisPopGrid$area = thisPopGrid$area
  # thisPopGrid$constituency = thisPopGrid$subarea
  # thisPopGridAdjusted = adjustPopGrid(thisPopGrid, poppconAdjusted, "Constituency")
  
  thisPopGrid = makePopIntegrationTab(kmRes=resolutions[i], pop=pop, domainPoly=kenyaPoly, 
                                      eastLim=eastLim, northLim=northLim, 
                                      mapProjection=SUMMER::projKenya, poppa=poppaKenya, 
                                      poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                      areaMapDat=adm1, subareaMapDat=adm2, 
                                      mean.neighbor=meanNeighbors[i], 
                                      delta=deltas[i], 
                                      areaPolygonSubsetI=30)
  
  thisPopGridAdjusted = adjustPopMat(thisPopGrid, poppaTarget=poppsubKenyaNeonatal, adjustBy="subarea")
  
  popGrids = c(popGrids, list(thisPopGrid))
  popGridsAdjusted = c(popGridsAdjusted, list(thisPopGridAdjusted))
}
sortI = sort(resolutions, index.return=TRUE)$ix
allResolutions = resolutions[sortI]
popGrids = popGrids[sortI]
popGridsAdjusted = popGridsAdjusted[sortI]
names(popGrids) = paste0("popGrid", allResolutions)
names(popGridsAdjusted) = paste0("popGridAdjusted", allResolutions)
ns = sapply(popGrids, nrow)
sum(sapply(popGrids, nrow))
endIs = cumsum(ns)
startIs = c(1, endIs[-length(endIs)]+1)

save(popGrids, popGridsAdjusted, ns, endIs, startIs, file="savedOutput/simpleExample/popGridResTest.RData")
out = load("savedOutput/simpleExample/popGridResTest.RData")

# combine the grids into 1 matrix for evaluating the SPDE model all at once
popMatCombined = do.call("rbind", popGrids)

# fit SPDE cluster level risk model over all integration points
spdeFitN = fitSPDEKenyaDat(dat, nPostSamples=10000, popMat=popMatCombined)

# remove unnecessary parts of the spdeFitN object to save space
spdeFitN$mod$.args = NULL
spdeFitN$mod$all.hyper = NULL

# apply aggregation models at each resolution
nSamples = c(c(500, 1000), rep(10000, length(resolutions)-1))
separateUDraws = list()
for(i in 1:length(popGrids)) {
  # obtain the grids at this resolution
  thisNSamples = nSamples[i]
  
  # obtain model output at this resolution
  thisResolutionI = startIs[i]:endIs[i]
  separateUDraws = c(separateUDraws, spdeFitN$uDraws[thisResolutionI,1:thisNSamples])
}
spdeFitN$uDraws = NULL

# save spdeFitN and its relevant uDraws
save(spdeFitN, separateUDraws, file="savedOutput/simpleExample/spdeFitNandUDraws.RData")
out = load("savedOutput/simpleExample/spdeFitNandUDraws.RData")

# run the actual analysis
aggResultsN = list()
for(i in 1:length(popGrids)) {
  # obtain the grids at this resolution
  thisNSamples = nSamples[i]
  thisPopMat = popGrids[[i]]
  thisPopMatAdjusted = popGridsAdjusted[[i]]
  
  # obtain model output at this resolution
  thisResolutionI = startIs[i]:endIs[i]
  # thisUDraws = spdeFitN$uDraws[thisResolutionI,1:thisNSamples]
  thisUDraws = separateUDraws[[i]]
  sigmaEpsilonDraws = spdeFitN$sigmaEpsilonDraws[1:thisNSamples]
  
  # thisAggResultsN = modLCPB(thisUDraws, sigmaEpsilonDraws, easpaN, thisPopMat, 
  #                           thisPopMatAdjusted, doLCPb=TRUE, doIHMEModel=TRUE, 
  #                           constituencyPop=poppconN, ensureAtLeast1PerConstituency=TRUE, 
  #                           logisticApproximation=FALSE, verbose=TRUE, 
  #                           fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
  #                           stopOnFrameMismatch=FALSE)
  
  thisAggResultsN = simPopCustom(thisUDraws, sigmaEpsilonDraws, easpaSimple, thisPopMat, 
                            thisPopMatAdjusted, doFineScaleRisk=TRUE, doIHMERisk=TRUE, 
                            doSmoothRisk=TRUE, 
                            poppsub=poppsubSimple, min1PerSubarea=TRUE, 
                            doSmoothRiskLogisticApprox=FALSE, 
                            fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1)
  
  aggResultsN = c(aggResultsN, list(thisAggResultsN))
}
names(aggResultsN) = paste0("aggResultsN", resolutions)
save(popGrids, popGridsAdjusted, aggResultsN, file="savedOutput/simpleExample/gridResolutionTestNairobi.RData")

out = load("savedOutput/simpleExample/gridResolutionTestNairobi.RData")

##### Plot results ----

# truePrevalenceConstituencyKenya
# truePrevalenceCountyKenya

# Calculate RMSE, 80% Coverage
predsConstituency = matrix(nrow=length(truePrevalenceConstituencyKenya), ncol=length(resolutions))
residsConstituencySmoothRisk = list()
residsConstituencyRisk = list()
residsConstituencyPrevalence = list()
residsConstituencyGriddedRisk = list()
for(i in 1:length(resolutions)) {
  thesePreds = rowMeans(aggResultsN[[i]]$aggregatedResultslcpb$constituencyMatrices$p)
  theseResidsSmoothRisk = sweep(aggResultsN[[i]]$aggregatedResultslcpb$constituencyMatrices$p, 1, truePrevalenceConstituencyKenya, "-")
  theseResidsRisk = sweep(aggResultsN[[i]]$aggregatedResultsLCPb$constituencyMatrices$p, 1, truePrevalenceConstituencyKenya, "-")
  theseResidsPrevalence = sweep(aggResultsN[[i]]$aggregatedResultsLCPB$constituencyMatrices$p, 1, truePrevalenceConstituencyKenya, "-")
  theseResidsGriddedRisk = sweep(aggResultsN[[i]]$aggregatedResultsIHME$constituencyMatrices$p, 1, truePrevalenceConstituencyKenya, "-")
  predsConstituency[,i] = thesePreds
  residsConstituencySmoothRisk = c(residsConstituencySmoothRisk, list(theseResidsSmoothRisk))
  residsConstituencyRisk = c(residsConstituencyRisk, list(theseResidsRisk))
  residsConstituencyPrevalence = c(residsConstituencyPrevalence, list(theseResidsPrevalence))
  residsConstituencyGriddedRisk = c(residsConstituencyGriddedRisk, list(theseResidsGriddedRisk))
}
residsConstituency = sweep(predsConstituency, 1, truePrevalenceConstituencyKenya, "-")

lowConstituencySmoothRisk = sapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
lowConstituencyRisk = sapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
lowConstituencyPrevalence = sapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
lowConstituencyGriddedRisk = sapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.1)})})
highConstituencySmoothRisk = sapply(residsConstituencySmoothRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
highConstituencyRisk = sapply(residsConstituencyRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
highConstituencyPrevalence = sapply(residsConstituencyPrevalence, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})
highConstituencyGriddedRisk = sapply(residsConstituencyGriddedRisk, function(mat) {apply(mat, 1, function(x) {quantile(x, probs=.9)})})

CIWidthSmoothRisk = highConstituencySmoothRisk - lowConstituencySmoothRisk
CIWidthRisk = highConstituencyRisk - lowConstituencyRisk
CIWidthPrevalence = highConstituencyPrevalence - lowConstituencyPrevalence
CIWidthGriddedRisk = highConstituencyGriddedRisk - lowConstituencyGriddedRisk

meanCIWidthSmoothRisk = colMeans(CIWidthSmoothRisk)
meanCIWidthRisk = colMeans(CIWidthRisk)
meanCIWidthPrevalence = colMeans(CIWidthPrevalence)
meanCIWidthGriddedRisk = colMeans(CIWidthGriddedRisk)

inCISmoothRisk = (0 <= highConstituencySmoothRisk) & (0 >= lowConstituencySmoothRisk)
inCIRisk = (0 <= highConstituencyRisk) & (0 >= lowConstituencyRisk)
inCIPrevalence = (0 <= highConstituencyPrevalence) & (0 >= lowConstituencyPrevalence)
inCIGriddedRisk = (0 <= highConstituencyGriddedRisk) & (0 >= lowConstituencyGriddedRisk)

coverageSmoothRisk = colMeans(inCISmoothRisk)
coverageRisk = colMeans(inCIRisk)
coveragePrevalence = colMeans(inCIPrevalence)
coverageGriddedRisk = colMeans(inCIGriddedRisk)

# Plot central predictions versus resolution
ylim = range(c(predsConstituency))
pdf("figures/gridResolutionTest/predictionVRes.pdf", width=5, height=5)
cols = rainbow(4)
thisFrame = data.frame(predsConstituency)
boxplot(predsConstituency, names=resolutions, col="skyblue", 
        main="", xlab="Resolution", ylab="Prediction")
dev.off()

# Plot CI Widths versus resolution and model
CIWidth = c(c(CIWidthSmoothRisk), c(CIWidthRisk), c(CIWidthPrevalence), c(CIWidthGriddedRisk))
tempRes = resolutions[col(CIWidthSmoothRisk)]
tempCon = factor(as.character(poppcon$Constituency[constituenciesN][col(CIWidthSmoothRisk)]))
N=length(tempCon)
CIWidthFrame = data.frame(Constituency=rep(tempCon, 4), 
                          Resolution=rep(tempRes, 4), 
                          Model=factor(c(rep("Smooth risk", N), rep("Risk", N), 
                                         rep("Prevalence", N), rep("Gridded risk", N)), 
                                       levels=c("Smooth risk", "Risk", "Prevalence", "Gridded risk")), 
                          CIWidth=CIWidth)

pdf("figures/gridResolutionTest/CIWidthVRes.pdf", width=7, height=5)
ggplot(CIWidthFrame, aes(factor(Resolution), CIWidth, fill=factor(Model))) + 
  geom_boxplot(position="dodge2") + scale_y_continuous(trans="log10") +
  labs(x="Grid resolution (km)", y="80% credible interval width", fill="Model") + 
  theme_classic()
dev.off()


pdf("figures/gridResolutionTest/CoverageVRes.pdf", width=5, height=5)
pchs = 15:18
cols = rainbow(4)
ylim = range(c(coverageSmoothRisk), c(coverageRisk), c(coveragePrevalence), c(coverageGriddedRisk))
ylim = c(.2, 1)
plot(resolutions*.97, coverageSmoothRisk, pch=pchs[1], col=cols[1], 
     ylim=ylim, ylab="80% coverage", xlab="Grid resolution (km)", 
     log="x")
abline(a=.8, b=0, lty=2)
points(resolutions*.99, coverageRisk, pch=pchs[2], col=cols[2])
points(resolutions*1.01, coveragePrevalence, pch=pchs[3], col=cols[3])
points(resolutions*1.03, coverageGriddedRisk, pch=pchs[4], col=cols[4])
legend("right", c("Smooth risk", "Risk", "Prevalence", "Gridded risk"), 
       pch=pchs, col=cols)
dev.off()

