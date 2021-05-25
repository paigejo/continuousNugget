
# install_github("https://github.com/richardli/SUMMER/tree/dev")
library(SUMMER)

data(KenyaPopulationData)
# out = load("~/Downloads/SUMMER-master/data/kenyaPopulationData.rda")
# save(easpaKenya, easpaKenyaNeonatal, popMatKenya, popMatKenyaNeonatal,
#      poppaKenya, poppsubKenya,
#      file="~/Downloads/SUMMER-master/data/kenyaPopulationData.rda")

##### first, we need the shapefiles for the areas (Admin1 areas) and 
##### subareas (Admin2 areas) along with the country, Kenya.
data(kenyaMaps)
# out = load("~/Downloads/SUMMER-master/data/kenyaMaps.rda")
# save(kenyaMesh, adm2Kenya, adm1Kenya, kenyaPoly, provinceMapKenya,
#      file="~/Downloads/SUMMER-master/data/kenyaMaps.rda")
# save(adm2Kenya, adm1Kenya, file="~/Downloads/SUMMER-master/R/sysdata.rda")

# set east/north limits in km (based on SUMMER::projKenya)
eastLim = c(-110.6405, 832.4544)
northLim = c(-555.1739, 608.7130)

# using census frame and map information, create 1 km resolution population 
# density map, normalized within urban/rural x Admin1 areas to get correct 
# total populations. Finer resolution grids can be used to generate poppsub 
# for more accurate population simulations, especially at the Admin2 level
if(FALSE) {
  # WARNING: this code takes a while to run. Can optionally just load in poppsubKenya
  pixelGridFine = generatePixelGrid(kmRes=1, kenyaPoly, eastLim, northLim, mapProjection=projKenya)
  subareas = getAdmin2Kenya(cbind(pixelGridFine$lon, pixelGridFine$lat), mean.neighbor=500, delta=.05)
  subareas = subareas$constituencyNames
  areas = admin2ToAdmin1Kenya(subareas)
  
  # popKenya is a population density raster that isn't in SUMMER
  popMatFine = makeInterpPopMat(popKenya, pixelGridFine, areas=areas, subareas=subareas, 
                                 poppa=poppaKenya, poppsub=NULL, stratifyByUrban=TRUE)
  poppsub = poppRegionFromPopMat(popMatFine, popMatFine$subareas)
  poppsubKenya = cbind(subarea=poppsub$region, area=admin2ToAdmin1Kenya(poppsub$region), poppsub[,-1])
  
  # make the 5 km resolution pixellated grid and get associated population densities, urbanicitites, 
  # areas, and subareas using precomputed populations of subareas from extra fine scale grid
  pixelGrid = generatePixelGrid(kmRes=5, kenyaPoly, eastLim, northLim, mapProjection=projKenya)
  popMat = makeInterpPopMat(popKenya, pixelGrid, areas=popMatKenya$area, subareas=popMatKenya$subarea, 
                            poppa=poppaKenya, poppsub=poppsubKenya, stratifyByUrban=TRUE)
  # check to make sure population grid is normalized correctly
  # test = poppRegionFromPopMat(popMat, popMat$subareas)
  # head(test)
  # head(poppsubKenya)
  # any(abs(test$popUrb - poppsubKenya$popUrb) > .001)
  popMatKenya = popMat
  
  # now adjust the census frame using estimated number of neonatals per household in 
  # urban and rural areas. These numbers are based on microcensus data
  # householdsPerEAUrb = 92.81475
  mothersPerHouseholdUrb = 0.3497151
  childrenPerMotherUrb = 1.295917
  # householdsPerEARur = 87.7717
  mothersPerHouseholdRur = 0.4787696
  childrenPerMotherRur = 1.455222
  # targetPopPerStratumUrban = easpaKenya$EAUrb * householdsPerEAUrb * mothersPerHouseholdUrb *
  #   childrenPerMotherUrb
  # targetPopPerStratumRural = easpaKenya$EARur * householdsPerEARur * mothersPerHouseholdRur *
  #   childrenPerMotherRur
  targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * childrenPerMotherUrb
  targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * childrenPerMotherRur
  easpaKenyaNeonatal = easpaKenya
  easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
  easpaKenyaNeonatal$popRur = targetPopPerStratumRural
  easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + easpaKenyaNeonatal$popRur
  easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / easpaKenyaNeonatal$popTotal
  easpaKenyaNeonatal$pctTotal = 100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)
  
  # generate the target population density by adjusting the current population density grid 
  # at the Admin1 level
  popMatKenyaNeonatal = adjustPopMat(popMatKenya, easpaKenyaNeonatal)
}

##### We are done making the relevant census frame and population density info. 
##### Now we make a model for the risk. We will use an SPDE model with these 
##### parameters for the linear predictor on the logist scale, which are chosen 
##### to be of practical interest:
beta0=-2.9 # intercept
gamma=-1 # urban effect
rho=(1/3)^2 # spatial variance
effRange = 400 # effective spatial range in km
sigmaEpsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation

# simulate the population! Note that this produces multiple dense nEA x nsim and nPixel x nsim 
# matrices. In the future sparse matrices will and chunk by chunk computations may be incorporated.
simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
                    poppsub=poppsubKenya, spdeMesh=kenyaMesh, margVar=rho, tausq=sigmaEpsilon^2, 
                    gamma=gamma, effRange=effRange, beta0=beta0, 
                    seed=123, inla.seed=12, nHHSampled=25, 
                    stratifyByUrban=TRUE, subareaLevel=TRUE, 
                    doFineScaleRisk=TRUE, 
                    min1PerSubarea=TRUE)

# plot EA level (aggregated to pixels by quilt.plot) simulated risk and prevalence
require(fields)
eaDat = simPop$pixelPop$eaDat[[1]]
quilt.plot(eaDat$lon, eaDat$lat, eaDat$Z)
quilt.plot(eaDat$lon, eaDat$lat, eaDat$N)
quilt.plot(eaDat$lon, eaDat$lat, eaDat$pFineScalePrevalence)
quilt.plot(eaDat$lon, eaDat$lat, eaDat$pFineScaleRisk)
quilt.plot(eaDat$lon, eaDat$lat, eaDat$pFineScaleRisk, zlim=c(0, .09))

# plot pixel level risk (without nugget effect)
uDraws = simPop$uDraws
quilt.plot(popMat$lon, popMat$lat, expit(uDraws), zlim=c(0, .1))
quilt.plot(popMat$lon, popMat$lat, popMatKenyaNeonatal$pop)

# plot Admin2 level simulated risk and prevalence
subareaPop = simPop$subareaPop
datAdmin2 = data.frame(subarea=sort(unique(adm2Kenya@data$CONSTITUEN)), 
                 prevalence=subareaPop$fineScalePrevalence$p, 
                 risk=subareaPop$fineScaleRisk$p, 
                 Z=subareaPop$fineScalePrevalence$Z, 
                 N=subareaPop$fineScalePrevalence$N)
prevalenceRange = range(datAdmin2$prevalence)
pdf("figures/bookdown/prevalenceAdm2.pdf", width=5, height=5)
mapPlot(data=datAdmin2, variables=names(dat)[2], by.data=names(dat)[1], 
        geo=adm2Kenya, by.geo="CONSTITUEN", is.long=FALSE, ylim=prevalenceRange)
dev.off()
pdf("figures/bookdown/riskAdm2.pdf", width=5, height=5)
mapPlot(data=datAdmin2, variables=names(dat)[3], by.data=names(dat)[1], 
        geo=adm2Kenya, by.geo="CONSTITUEN", is.long=FALSE, ylim=prevalenceRange)
dev.off()
pdf("figures/bookdown/ZAdm2.pdf", width=5, height=5)
mapPlot(data=datAdmin2, variables=names(dat)[4], by.data=names(dat)[1], 
        geo=adm2Kenya, by.geo="CONSTITUEN", is.long=FALSE)
dev.off()
pdf("figures/bookdown/NAdm2.pdf", width=5, height=5)
mapPlot(data=datAdmin2, variables=names(dat)[5], by.data=names(dat)[1], 
        geo=adm2Kenya, by.geo="CONSTITUEN", is.long=FALSE)
dev.off()

pdf("figures/bookdown/prevalenceVsRiskAdm2.pdf", width=5, height=5)
plot(datAdmin2$risk, datAdmin2$prevalence, 
     xlim=prevalenceRange, ylim=prevalenceRange, 
     main="Admin2 Risk Versus Prevalence", xlab="Risk", ylab="Prevalence", 
     cex=.5, pch=19, col="blue")
abline(0, 1)
dev.off()

# plot Admin1 level simulated risk and prevalence
areaPop = simPop$areaPop
datAdmin1 = data.frame(area=sort(unique(adm1Kenya@data$NAME_1)), 
                 prevalence=areaPop$aggregatedResultsPrevalence$p, 
                 risk=areaPop$aggregatedResultsRisk$p, 
                 Z=areaPop$aggregatedResultsPrevalence$Z, 
                 N=areaPop$aggregatedResultsPrevalence$N)
pdf("figures/bookdown/prevalenceAdm1.pdf", width=5, height=5)
mapPlot(data=datAdmin1, variables=names(datAdmin1)[2], by.data=names(datAdmin1)[1], 
        geo=adm1Kenya, by.geo="NAME_1", is.long=FALSE, ylim=prevalenceRange)
dev.off()
pdf("figures/bookdown/riskAdm1.pdf", width=5, height=5)
mapPlot(data=datAdmin1, variables=names(datAdmin1)[3], by.data=names(datAdmin1)[1], 
        geo=adm1Kenya, by.geo="NAME_1", is.long=FALSE, ylim=prevalenceRange)
dev.off()
pdf("figures/bookdown/ZAdm1.pdf", width=5, height=5)
mapPlot(data=datAdmin1, variables=names(datAdmin1)[4], by.data=names(datAdmin1)[1], 
        geo=adm1Kenya, by.geo="NAME_1", is.long=FALSE)
dev.off()
pdf("figures/bookdown/NAdm1.pdf", width=5, height=5)
mapPlot(data=datAdmin1, variables=names(datAdmin1)[5], by.data=names(datAdmin1)[1], 
        geo=adm1Kenya, by.geo="NAME_1", is.long=FALSE)
dev.off()

pdf("figures/bookdown/prevalenceVsRiskAdm1.pdf", width=5, height=5)
plot(datAdmin1$risk, datAdmin1$prevalence, 
     xlim=prevalenceRange, ylim=prevalenceRange, 
     main="Admin1 Risk Versus Prevalence", xlab="Risk", ylab="Prevalence", 
     cex=.8, pch=19, col="blue")
abline(0, 1)
dev.off()


head(adm2Kenya@data)



