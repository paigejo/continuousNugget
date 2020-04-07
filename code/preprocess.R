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
empiricalDistributions = getSurveyEmpiricalDistributions2(maxAge=0)
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









