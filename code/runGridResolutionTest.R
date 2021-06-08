# this script is for comparing performance between models

simDatKenya = generateSimDataSetsLCPB(nsim=1, adjustedPopMat=popMatSimpleAdjusted, 
                                      fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                                      logisticApproximation=FALSE, 
                                      dataSaveDirectory="~/git/continuousNugget/savedOutput/simpleExample/", 
                                      seed=1, inla.seed=1L, simPopOnly=FALSE, returnEAinfo=TRUE, 
                                      easpa=NULL, popMat=NULL, constituencyPop=poppcon)

# get the data and the true population
dat = simDatKenya$SRSDat$clustDat[[1]]
constituenciesN = poppcon$County == "Nairobi"
countyN = sort(unique(poppc$County)) == "Nairobi"
truePrevalenceConstituencyKenya = simDatKenya$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$constituencyMatrices$p[constituenciesN,1]
truePrevalenceCountyKenya = simDatKenya$simulatedEAs$aggregatedPop$aggregatedResultsLCPB$countyMatrices$p[countyN,1]

# construct integration grids at different resolutions
# resolutions = c(1, 10, 20, 40, 60, 80)
# deltas = c(.05, .2, .4, .8, 1.2, 1.4)
# meanNeighbors = c(500, rep(50, 5))
resolutions = c(.5, 1, 2, 5, 10, 20)
deltas = c(.025, .05, .1, .1, .2, .4)
meanNeighbors = c(rep(500, 3), rep(50, 3))
popGrids = list()
popGridsAdjusted = list()
for(i in 1:length(resolutions)) {
  print(paste0("Creating integration grid at ", resolutions[i], " km resolution"))
  thisPopGrid = makeInterpPopGrid(kmRes=resolutions[i], 
                                  mean.neighbor=meanNeighbors[i], 
                                  delta=deltas[i], conMap=adm2, poppcon=poppcon, 
                                  mapDat=adm1, polygonSubsetI=30) # Nairobi is the 30th one
  thisPopGrid$area = thisPopGrid$admin1
  thisPopGrid$constituency = thisPopGrid$admin2
  thisPopGridAdjusted = adjustPopGrid(thisPopGrid, poppcon, "Constituency")
  
  # subset grids to be over only Nairobi
  thisPopGrid = thisPopGrid[thisPopGrid$admin1 == "Nairobi",]
  thisPopGridAdjusted = thisPopGridAdjusted[thisPopGridAdjusted$admin1 == "Nairobi",]
  
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

# combine the grids
popMatCombined = do.call("rbind", popGrids)

# fit SPDE cluster level risk model over all integration points
nSamples = 10000
spdeFitN = fitSPDEKenyaDat(dat, nPostSamples=nSamples, popMat=popMatCombined)
sigmaEpsilonDraws = spdeFitN$sigmaEpsilonDraws

# apply aggregation models at each resolution
aggResultsN = list()
for(i in 1:length(popGrids)) {
  # obtain the grids at this resolution
  thisPopMat = popGrids[[i]]
  thisPopMatAdjusted = popGridsAdjusted[[i]]
  
  # obtain model output at this resolution
  thisResolutionI = startIs[i]:endIs[i]
  thisUDraws = spdeFitN$uDraws[thisResolutionI,]
  
  thisAggResultsN = modLCPB(thisUDraws, sigmaEpsilonDraws, easpaN, thisPopMat, 
                            thisPopMatAdjusted, doLCPb=TRUE, doIHMEModel=TRUE, 
                            constituencyPop=poppconN, ensureAtLeast1PerConstituency=TRUE, 
                            logisticApproximation=FALSE, verbose=TRUE, 
                            fixPopPerEA=25, fixHHPerEA=25, fixPopPerHH=1, 
                            stopOnFrameMismatch=FALSE)
  
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

