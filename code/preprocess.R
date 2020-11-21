# this script is for setting up basic global variables from preprocessing

# map data
adm0 = readRDS(paste0(dataDirectory, "mapData/KEN_adm0.rds"))
adm1 = readRDS(paste0(dataDirectory, "mapData/KEN_adm1.rds"))
# adm2 = readRDS(paste0(dataDirectory, "mapData/KEN_adm2.rds"))
require(rgdal)
adm2 = readOGR(dsn = "data/mapData/constituencies", layer = "constituencies")
# MARAKWET WEST was mislabeled as being in WEST POKOT instead of Elgeyo-Marakwet
# we will also relabel Lugari and Likuyani as being in Bungoma for consistency's sake
adm2@data[which(grepl("MARAKWET WEST", adm2@data$CONSTITUEN)),]$COUNTY_NAM = "ELGEYO-MARAKWET"
adm2@data[which(grepl("LUGARI", adm2@data$CONSTITUEN)),]$COUNTY_NAM = "BUNGOMA"
adm2@data[which(grepl("LIKUYANI", adm2@data$CONSTITUEN)),]$COUNTY_NAM = "BUNGOMA"
require(spatialEco)
# fix discrepancies due to border disputes, combine small (<25 km^2) constituencies, remove NAs
adm2 = sp.na.omit(adm2, col.name="CONSTITUEN")
require(tools)
# fix constituency and county strings. County strings must match those of other data frames
adm2@data$CONSTITUEN = toTitleCase(tolower(adm2@data$CONSTITUEN))
adm2@data$COUNTY_NAM = toTitleCase(tolower(adm2@data$COUNTY_NAM))
adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Elegeyo-Marakwet"] = "Elgeyo Marakwet"
adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Elgeyo-Marakwet" ] = "Elgeyo Marakwet"
adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Tharaka - Nithi" ] = "Tharaka-Nithi"
adm2@data$COUNTY_NAM[adm2@data$COUNTY_NAM == "Trans Nzoia" ] = "Trans-Nzoia"

adm2 = combineConstituencies(adm2, threshold=50)
adm2 = makeInBorder(adm2)
adm2@data$Shape_Area = getArea(thisMap=adm2, nameVar="CONSTITUEN")

regionMap = constructRegions()

# test = removeConstituencyGaps(adm2)
save(adm2, adm1, adm0, regionMap, file=paste0(globalDirectory, "adminMapData.RData"))
load(paste0(globalDirectory, "adminMapData.RData"))

# county to region mapping
ctp = read.csv(paste0(dataDirectory, "mapData/kenya-prov-county-map.csv"))
save(ctp, file=paste0(globalDirectory, "kenya-prov-county-map.RData"))

# generate 5km population density grid over Kenya
popGrid = makeInterpPopGrid(kmRes=5)
save(popGrid, file=paste0(globalDirectory, "popGrid.RData"))
popGridAdjusted = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE)
save(popGridAdjusted, file=paste0(globalDirectory, "popGridAdjusted.RData"))
popGridAdjustedWomen = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE, "women")
save(popGridAdjustedWomen, file=paste0(globalDirectory, "popGridAdjustedWomen.RData"))

# empirical distributions
empiricalDistributions = getSurveyEmpiricalDistributions2(maxAge=4) # maxAge is 4 years old because we want number of children in the 5 year period
empiricalDistributions = c(empiricalDistributions, getSurveyEmpiricalDistributionsWomen())
save(empiricalDistributions, file=paste0(globalDirectory, "empiricalDistributions.RData"))

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

# simulate 1000 draws of N at each EA:
# out = simNsFull(1, includeUrban=TRUE)
out = simNsFull(1000, includeUrban=TRUE)
NcsUrban = out$NcsUrban
NcsRural = out$NcsRural
save(NcsUrban, NcsRural, file=paste0(globalDirectory, "NcsUrbanRural.RData"))

# Ncs = simNsFull(1, includeUrban=FALSE)
Ncs = simNsFull(1000, includeUrban=FALSE)
save(Ncs, file=paste0(globalDirectory, "Ncs.RData"))

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
out = aggregate(popGridFine$popOrig, by=list(constituency=as.character(popGridFine$admin2), urban=popGridFine$urban), FUN=sum)
poppcon = data.frame(Constituency=constituencies, County=constituencyToCounty(constituencies), popUrb=out[(length(constituencies) + 1):(2*length(constituencies)), 3], 
                     popRur=out[1:length(constituencies), 3])

# make sure constituencies without any urban or rural population in the coarse grid have none in the fine grid as well
for(i in 1:nrow(poppcon)) {
  thisRow = poppcon[i,]
  thisIUrb = (popGrid$admin2 == thisRow$Constituency) & (popGrid$urban == TRUE)
  thisIRur = (popGrid$admin2 == thisRow$Constituency) & (popGrid$urban == FALSE)
  if(sum(popGrid$popOrig[thisIUrb]) == 0) {
    poppcon$popUrb[i] = 0
  }
  if(sum(popGrid$popOrig[thisIRur]) == 0) {
    poppcon$popRur[i] = 0
  }
}

# normalize population so it sums to the total population of Kenya in urban/rural areas. Also calculate total population
poppcon$popUrb = poppcon$popUrb * (sum(poppc$popUrb) / sum(poppcon$popUrb))
poppcon$popRur = poppcon$popRur * (sum(poppc$popRur) / sum(poppcon$popRur))
poppcon$popTotal = poppcon$popUrb + poppcon$popRur

save(poppcon, file=paste0(globalDirectory, "poppcon.RData"))

# check that constituency crossed with urbanicity of the dataset make sense
out = aggregate(mort$admin2, by=list(mort$admin2, mort$urban), FUN=length, drop=FALSE)
out = cbind(out[1:273,1:3], out[274:546,3])
out[is.na(out)] = 0
names(out) = c("Constituency", "urban", "nClustRur", "nClustUrb")
out = out[,c(1:2, 4, 3)]
out$nClustTotal = out$nClustRur + out$nClustUrb

cbind(poppcon[,1:2], out[,3:5]/poppcon[,3:5])[,5]

