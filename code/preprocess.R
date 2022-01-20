# this script is for setting up basic global variables from preprocessing

# # map data
# adm0 = readRDS(paste0(dataDirectory, "mapData/KEN_adm0.rds"))
# adm1 = readRDS(paste0(dataDirectory, "mapData/KEN_adm1.rds"))
# # adm2 = readRDS(paste0(dataDirectory, "mapData/KEN_adm2.rds"))
# require(rgdal)
# adm2 = readOGR(dsn = "data/mapData/constituencies", layer = "constituencies")
# # MARAKWET WEST was mislabeled as being in WEST POKOT instead of Elgeyo-Marakwet
# # we will also relabel Lugari and Likuyani as being in Bungoma for consistency's sake
# adm2@data[which(grepl("MARAKWET WEST", adm2@data$CONSTITUEN)),]$COUNTY_NAM = "ELGEYO-MARAKWET"
# adm2@data[which(grepl("LUGARI", adm2@data$CONSTITUEN)),]$COUNTY_NAM = "BUNGOMA"
# adm2@data[which(grepl("LIKUYANI", adm2@data$CONSTITUEN)),]$COUNTY_NAM = "BUNGOMA"
# require(spatialEco)
# # fix discrepancies due to border disputes, combine small (<25 km^2) constituencies, remove NAs
# adm2 = sp.na.omit(adm2, col.name="CONSTITUEN")
# require(tools)
# # fix constituency and county strings. County strings must match those of other data frames
# adm2@data$CONSTITUEN = toTitleCase(tolower(adm2@data$CONSTITUEN))
# adm2@data$COUNTY_NAM = toTitleCase(tolower(adm2@data$COUNTY_NAM))
# adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Elegeyo-Marakwet"] = "Elgeyo Marakwet"
# adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Elgeyo-Marakwet" ] = "Elgeyo Marakwet"
# adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Tharaka - Nithi" ] = "Tharaka-Nithi"
# adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Trans Nzoia" ] = "Trans-Nzoia"
# 
# # adm2 = combineConstituencies(adm2, threshold=50)
# adm2 = makeInBorder(adm2)
# adm2@data$Shape_Area = getArea(thisMap=adm2, areaLevel="Constituency")
# 
# # test = removeConstituencyGaps(adm2)
# save(adm2, adm1, adm0, file=paste0(globalDirectory, "adminMapData.RData"))
# load(paste0(globalDirectory, "adminMapData.RData"))
# 
# if(FALSE) {
#   # test mapping of shapefiles
#   plotMapDat(mapDat=adm2, new=TRUE)
#   plotMapDat(mapDat=adm1, new=TRUE)
#   plotMapDat(mapDat=adm0, new=TRUE)
# }

githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "kenyaMaps.rda?raw=true")
tempDirectory = "~/"
popMatFilename = paste0(tempDirectory, "/kenyaMaps.RData")
if(!file.exists(popMatFilename)) {
  download.file(githubURL,popMatFilename)
}

# load it in
out = load(popMatFilename)
out

# save it
save(adm1, adm2, kenyaPoly, "savedOoutput/global/kenyaMapData.RData")

# download 2014 Kenya population density and associated TIF file
tempDirectory = "~/git/continuousNugget/data/popData"
saveDirectory = "~/git/continuousNugget/savedOutput/global"
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "Kenya2014Pop/pop.rda?raw=true")
popFilename = paste0(saveDirectory, "/pop.rda")
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

# make sure this is correct for re-projections
pop@file@name = "~/git/continuousNugget/data/popData/worldpop_total_1y_2014_00_00.tif"


# save it
save(pop, file=popFilename)
save(pop, file="savedOutput/global/pop.RData")

# county to region mapping
ctp = read.csv(paste0(dataDirectory, "mapData/kenya-prov-county-map.csv"))
save(ctp, file=paste0(globalDirectory, "kenya-prov-county-map.RData"))

# for precomputating populations of constituencies possibly crossed with urban/rural
popGridFine = makeInterpPopGrid(kmRes=1, mean.neighbor=500, delta=.05)
constituencies = sort(unique(popGridFine$admin2))
newRows = popGridFine[1:(2*length(constituencies)),]
newRows$popOrig = 0
newRows[1:length(constituencies),]$urban = FALSE
newRows[(length(constituencies) + 1):(2*length(constituencies)),]$urban = TRUE
newRows$admin2 = rep(constituencies, 2)
popGridFine = rbind(popGridFine, 
                    newRows)
out = aggregate(popGridFine$popOrig, by=list(constituency=as.character(popGridFine$admin2), urban=popGridFine$urban), FUN=sum, drop=FALSE)
poppcon = data.frame(Constituency=constituencies, County=constituencyToCounty(constituencies), popUrb=out[(length(constituencies) + 1):(2*length(constituencies)), 3], 
                     popRur=out[1:length(constituencies), 3])

# 
# # make sure constituencies without any urban or rural population in the coarse grid have none in the fine grid as well
# for(i in 1:nrow(poppcon)) {
#   thisRow = poppcon[i,]
#   thisIUrb = (popGrid$admin2 == thisRow$Constituency) & (popGrid$urban == TRUE)
#   thisIRur = (popGrid$admin2 == thisRow$Constituency) & (popGrid$urban == FALSE)
#   if(sum(popGrid$popOrig[thisIUrb]) == 0) {
#     poppcon$popUrb[i] = 0
#   }
#   if(sum(popGrid$popOrig[thisIRur]) == 0) {
#     poppcon$popRur[i] = 0
#   }
# }

# normalize population within each stratum to sum to the stratum population
countyI = match(poppcon$County, easpc$County) # county -> constituency
popInUrbanStratum = poppc$popUrb[countyI]
popInRuralStratum = poppc$popRur[countyI]

out = aggregate(poppcon$popUrb, by=list(County=poppcon$County), FUN=sum, drop=FALSE)
sortI = match(poppc$County, out$County) # sorted county -> county
out = out[sortI,]
normFactorUrban = popInUrbanStratum / out$x[countyI]
normFactorUrban[is.na(normFactorUrban)] = 0
poppcon$popUrb = poppcon$popUrb * normFactorUrban

# # test to make sure we did it right:
# out = aggregate(poppcon$popUrb, by=list(County=poppcon$County), FUN=sum)[sortI,]
# cbind(out, poppc)

out = aggregate(poppcon$popRur, by=list(County=poppcon$County), FUN=sum, drop=FALSE)[sortI,]
normFactorRural = popInRuralStratum / out$x[countyI]
normFactorRural[is.na(normFactorRural)] = 0
poppcon$popRur = poppcon$popRur * normFactorRural

# # test to make sure we did it right:
# out = aggregate(poppcon$popRur, by=list(County=poppcon$County), FUN=sum)[sortI,]
# cbind(out, poppc)

# normalize population so it sums to the total population of Kenya in urban/rural areas. Also calculate total population
poppcon$popUrb = poppcon$popUrb * (sum(poppc$popUrb) / sum(poppcon$popUrb))
poppcon$popUrb[is.na(poppcon$popUrb)] = 0
poppcon$popRur = poppcon$popRur * (sum(poppc$popRur) / sum(poppcon$popRur))
poppcon$popRur[is.na(poppcon$popRur)] = 0
poppcon$popTotal = poppcon$popUrb + poppcon$popRur

save(poppcon, file=paste0(globalDirectory, "poppcon.RData"))

# check that constituency crossed with urbanicity of the dataset make sense
# out = aggregate(mort$admin2, by=list(mort$admin2, mort$urban), FUN=length, drop=FALSE)
# out = cbind(out[1:273,1:3], out[274:546,3])
# out[is.na(out)] = 0
# names(out) = c("Constituency", "urban", "nClustRur", "nClustUrb")
# out = out[,c(1:2, 4, 3)]
# out$nClustTotal = out$nClustRur + out$nClustUrb

cbind(poppcon[,1:2], out[,3:5]/poppcon[,3:5])[,5]

# empirical distributions
empiricalDistributions = getSurveyEmpiricalDistributions2(maxAge=4) # maxAge is 4 years old because we want number of children in the 5 year period
empiricalDistributions = c(empiricalDistributions, getSurveyEmpiricalDistributionsWomen())
save(empiricalDistributions, file=paste0(globalDirectory, "empiricalDistributions.RData"))

# make adjusted poppcon (target population per constituency)
# easpcon = meanEAsPerCon()
# poppconAdjusted = easpcon
hhspcon = meanHHPerCon()
poppconAdjusted = hhspcon
# poppconAdjusted$popUrb = poppconAdjusted$meanUrbanEAs * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) *
#   ecdfExpectation(empiricalDistributions$childrenUrban)
# poppconAdjusted$popRur = poppconAdjusted$meanRuralEAs * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) *
#   ecdfExpectation(empiricalDistributions$childrenRural)
poppconAdjusted$popUrb = poppconAdjusted$meanUrbanHHs * ecdfExpectation(empiricalDistributions$mothersUrban) * 
  ecdfExpectation(empiricalDistributions$childrenUrban)
poppconAdjusted$popRur = poppconAdjusted$meanRuralHHs * ecdfExpectation(empiricalDistributions$mothersRural) * 
  ecdfExpectation(empiricalDistributions$childrenRural)
poppconAdjusted$popTotal = poppconAdjusted$popUrb + poppconAdjusted$popRur
save(poppconAdjusted, file="savedOutput/global/poppconAdjusted.RData")

# generate 5km population density grid over Kenya

# make the 5 km resolution pixellated grid and get associated population densities, urbanicitites, 
# counties, and constituencies using precomputed populations of constituencies from extra fine scale grid
popGrid = makeInterpPopGrid(kmRes=5, poppcon=poppcon, conMap=adm2)
save(popGrid, file=paste0(globalDirectory, "popGrid.RData"))
popGridAdjusted = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE, poppcon=poppcon, conMap=adm2)
save(popGridAdjusted, file=paste0(globalDirectory, "popGridAdjusted.RData"))
popGridAdjustedWomen = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE, "women", poppcon=poppcon, conMap=adm2)
save(popGridAdjustedWomen, file=paste0(globalDirectory, "popGridAdjustedWomen.RData"))

# make the 25 km resolution pixellated grid and get associated population densities, urbanicitites, 
# counties, and constituencies using precomputed populations of constituencies from extra fine scale grid
out = load("savedOutput/global/pop.RData")
popGridCoarse = makePopIntegrationTab(kmRes=25, pop=pop, domainPoly=kenyaPoly, 
                                    eastLim=eastLim, northLim=northLim, 
                                    mapProjection=SUMMER::projKenya, poppa=poppaKenya, 
                                    poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                    areaMapDat=adm1, subareaMapDat=adm2)

popGridCoarseAdjusted = adjustPopMat(popGridCoarse, poppaTarget=poppsubKenyaNeonatal, adjustBy="subarea")
save(popGridCoarse, file=paste0(globalDirectory, "popGridCoarse.RData"))
save(popGridCoarseAdjusted, file=paste0(globalDirectory, "popGridCoarseAdjusted.RData"))
# popGridCoarseAdjustedWomen = makeInterpPopGrid(kmRes=25, adjustPopSurface=TRUE, "women", poppcon=poppcon, conMap=adm2)
# save(popGridCoarseAdjustedWomen, file=paste0(globalDirectory, "popGridCoarseAdjustedWomen.RData"))

# datasets (these were already created in readDat3)
source("code/readDat3.R")
out = load(paste0(globalDirectory, "kenyaData.RData"))
# out = load(paste0(globalDirectory, "kenyaDataEd.RData")) # (education dataset not needed)

# update datasets' admin1 and admin2 values based on polygons
constituencies = getConstituency(cbind(mort$lon, mort$lat))$constituencyNames
counties = constituencyToCounty(constituencies)
mort$admin1 = counties
mort$admin2 = constituencies
save(mort, file=paste0(globalDirectory, "kenyaData.RData"))

# constituencies = getConstituency(cbind(ed$lon, ed$lat))$constituencyNames
# counties = constituencyToCounty(constituencies)
# ed$admin1 = counties
# ed$admin2 = constituencies
# save(ed, file=paste0(globalDirectory, "kenyaDataEd.RData"))

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

# # simulate 1000 draws of N at each EA:
# # out = simNsFull(1, includeUrban=TRUE)
# out = simNsFull(1000, includeUrban=TRUE)
# NcsUrban = out$NcsUrban
# NcsRural = out$NcsRural
# save(NcsUrban, NcsRural, file=paste0(globalDirectory, "NcsUrbanRural.RData"))
# 
# # Ncs = simNsFull(1, includeUrban=FALSE)
# Ncs = simNsFull(1000, includeUrban=FALSE)
# save(Ncs, file=paste0(globalDirectory, "Ncs.RData"))

# change urbanicity to take the pixel level values rather than the cluster level values:
test = getPixelIndex(cbind(mort$east, mort$north), popGrid, mort$admin1)
temp = mort
temp$urban = popGrid$urban[test]

# poppc modified to contain sampled population totals
out = aggregate(temp$n, by=list(admin1=temp$admin1, urban=temp$urban), FUN=sum, drop=FALSE)
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
out = aggregate(temp$n, by=list(admin1=temp$admin1, urban=temp$urban), FUN=length, drop=FALSE)
out[is.na(out[,3]), 3] = 0
# urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
out2 = cbind(out, rural=0)[48:94,]
out2[, 4] = out$x[1:47]
unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
easpcMort = clustpc

names(easpcMort) = names(easpc)
easpcMort[unsortedToSorted, 2:3] = out2[,3:4]
easpcMort$EATotal = easpcMort$EAUrb + easpcMort$EARur

# set the number of households (there are 25 per enumeration area)
easpcMort$HHUrb = easpcMort$EAUrb*25
easpcMort$HHRur = easpcMort$EARur*25
easpcMort$HHTotal = easpcMort$EATotal*25

save(poppcMort, file=paste0(globalDirectory, "poppcMort.RData"))
save(easpcMort, file=paste0(globalDirectory, "easpcMort.RData"))

# adjust poppc to be the target population
load(paste0(globalDirectory, "poppc.RData"))

# calculate the number of children per stratum using true total eas and empirical children per ea from census data
load("../U5MR/empiricalDistributions.RData")

# targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) *
#   ecdfExpectation(empiricalDistributions$childrenUrban)
# targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) *
#   ecdfExpectation(empiricalDistributions$childrenRural)
targetPopPerStratumUrban = easpc$HHUrb * ecdfExpectation(empiricalDistributions$mothersUrban) *
  ecdfExpectation(empiricalDistributions$childrenUrban)
targetPopPerStratumRural = easpc$HHRur * ecdfExpectation(empiricalDistributions$mothersRural) *
  ecdfExpectation(empiricalDistributions$childrenRural)
poppcAdjusted = poppc
poppcAdjusted$popUrb = targetPopPerStratumUrban
poppcAdjusted$popRur = targetPopPerStratumRural
poppcAdjusted$popTotal = targetPopPerStratumUrban + targetPopPerStratumRural
poppcAdjusted$pctTotal = round(100 * poppcAdjusted$popTotal / sum(poppcAdjusted$popTotal), digits=1)
poppcAdjusted$pctUrb = round(100 * poppcAdjusted$popUrb / poppcAdjusted$popTotal, digits=1)
save(poppcAdjusted, file=paste0(globalDirectory, "poppcAdjusted.RData"))






