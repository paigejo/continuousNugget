# this script is for setting up basic global variables from preprocessing

# map data
adm1 = readRDS(paste0(dataDirectory, "mapData/KEN_adm1.rds"))
adm0 = readRDS(paste0(dataDirectory, "mapData/KEN_adm0.rds"))
save(adm1, adm0, file=paste0(globalDirectory, "adminMapData.RData"))
load(paste0(globalDirectory, "adminMapData.RData"))

# county to region mapping
ctp = read.csv(paste0(dataDirectory, "mapData/kenya-prov-county-map.csv"))
save(ctp, file=paste0(globalDirectory, "kenya-prov-county-map.RData"))

# generate 5km population density grid over Kenya
popGrid = makeInterpPopGrid(kmRes=5)
save(popGrid, file=paste0(globalDirectory, "popGrid.RData"))
popGrid = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE)
save(popGrid, file=paste0(globalDirectory, "popGridAdjusted.RData"))
popGrid = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE, "women")
save(popGrid, file=paste0(globalDirectory, "popGridAdjustedWomen.RData"))

# empirical distributions
empiricalDistributions = getSurveyEmpiricalDistributions2(maxAge=4) # maxAge is 4 years old because we want number of children in the 5 year period
empiricalDistributions = c(empiricalDistributions, getSurveyEmpiricalDistributionsWomen())
save(empiricalDistributions, file=paste0(globalDirectory, "empiricalDistributions.RData"))

# datasets (these were already created in readKenyaData)
out = load(paste0(globalDirectory, "kenyaData.RData"))
# out = load(paste0(globalDirectory, "kenyaDataEd.RData")) # (education dataset not needed)

# set enumeration areas based on overall population density
load(paste0(globalDirectory, "popGrid.RData"))
numEAs = 96251
totalKenyaPop = 43*10^6 # from DHS 2014 survey final report page 2
kenyaEAs = simEAs(popGrid, numEAs, totalKenyaPop, fixNumUrbanAtTruth=TRUE, seed=123)
save(kenyaEAs, file=paste0(globalDirectory, "kenyaEAs.RData"))
load(paste0(globalDirectory, "kenyaEAs.RData"))

# load 5km population density grid over Kenya
load(paste0(globalDirectory, "popGrid.RData"))
load(paste0(globalDirectory, "popGridAdjusted.RData"))

# get limits of easting/northing for plotting
kenyaLonRange = c(33.5, 42)
kenyaLatRange = c(-5, 5.5)
tmp = projKenya(kenyaLonRange, kenyaLatRange)
eastLim = tmp[,1]
northLim = tmp[,2]
save(eastLim, northLim, file=paste0(globalDirectory, "lims.RData"))
load(paste0(globalDirectory, "lims.RData"))

# simulate 1000 draws of N at each EA:
# out = simNsFull(1, includeUrban=TRUE)
out = simNsFull(1000, includeUrban=TRUE)
NcsUrban = out$NcsUrban
NcsRural = out$NcsRural
save(NcsUrban, NcsRural, file=paste0(globalDirectory, "NcsUrbanRural.RData"))

# Ncs = simNsFull(1, includeUrban=FALSE)
Ncs = simNsFull(1000, includeUrban=FALSE)
save(Ncs, file=paste0(globalDirectory, "Ncs.RData"))

# poppc modified to contain sampled population totals
out = aggregate(mort$n, by=list(admin1=mort$admin1, urban=mort$urban), FUN=sum, drop=FALSE)
out[is.na(out[,3]), 3] = 0
# urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
out2 = cbind(out, rural=0)[48:94,]
out2[, 4] = out$x[1:47]
unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
poppcMort = poppc
poppcMort[unsortedToSorted, 2:3] = out2[,3:4]
poppcMort$popTotal = poppcMort$popUrb + poppcMort$popRur
poppcMort$pctTotal = 100 * poppcMort$popTotal * (1 / sum(poppcMort$popTotal))
poppcMort$pctUrb = 100 * poppcMort$popUrb / poppcMort$popTotal

## clustpc modified to contain sampled cluster and household totals
# first the number of sampled clusters
out = aggregate(mort$n, by=list(admin1=mort$admin1, urban=mort$urban), FUN=length, drop=FALSE)
out[is.na(out[,3]), 3] = 0
# urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
out2 = cbind(out, rural=0)[48:94,]
out2[, 4] = out$x[1:47]
unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
easpcMort = clustpc
names(easpcMort) = names(easpc)
easpcMort[unsortedToSorted, 2:3] = out2[,3:4]
easpcMort$EATotal = easpcMort$EAUrb + easpcMort$EARur

save(poppcMort, file=paste0(globalDirectory, "poppcMort.RData"))
save(easpcMort, file=paste0(globalDirectory, "easpcMort.RData"))


