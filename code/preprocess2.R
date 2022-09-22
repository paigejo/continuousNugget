##### This script is based on the SUMMER package example in ?simPopSPDE

## In this script we will create 5km resolution pixellated grid over Kenya, 
## and generate tables of estimated (both target and general) population 
## totals at the area (e.g. Admin-1) and subarea (e.g. Admin-2) levels.

# download Kenya GADM shapefiles from SUMMERdata github repository
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "kenyaMaps.rda?raw=true")
tempDirectory = "~/"
mapsFilename = paste0(tempDirectory, "/kenyaMaps.rda")
if(!file.exists(mapsFilename)) {
  download.file(githubURL,mapsFilename)
}

# load it in
out = load(mapsFilename)
out
adm1@data$NAME_1 = as.character(adm1@data$NAME_1)
adm1@data$NAME_1[adm1@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
adm1@data$NAME_1[adm1@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
adm2@data$NAME_1 = as.character(adm2@data$NAME_1)
adm2@data$NAME_1[adm2@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
adm2@data$NAME_1[adm2@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"

# some Admin-2 areas have the same name
adm2@data$NAME_2 = as.character(adm2@data$NAME_2)
adm2@data$NAME_2[(adm2@data$NAME_1 == "Bungoma") & 
                   (adm2@data$NAME_2 == "Lugari")] = "Lugari, Bungoma"
adm2@data$NAME_2[(adm2@data$NAME_1 == "Kakamega") & 
                   (adm2@data$NAME_2 == "Lugari")] = "Lugari, Kakamega"
adm2@data$NAME_2[(adm2@data$NAME_1 == "Meru") & 
                   (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Meru"
adm2@data$NAME_2[(adm2@data$NAME_1 == "Tharaka-Nithi") & 
                   (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Tharaka-Nithi"

# The spatial area of unknown 8 is so small, it causes problems unless its removed or 
# unioned with another subarea. Union it with neighboring Kakeguria:
newadm2 = adm2
unknown8I = which(newadm2$NAME_2 == "unknown 8")
newadm2$NAME_2[newadm2$NAME_2 %in% c("unknown 8", "Kapenguria")] <- 
  "Kapenguria + unknown 8"
admin2.IDs <- newadm2$NAME_2

library(maptools)
temp <- unionSpatialPolygons(newadm2, admin2.IDs)
tempData = newadm2@data[-unknown8I,]
tempData = tempData[order(tempData$NAME_2),]
newadm2 <- SpatialPolygonsDataFrame(temp, tempData, match.ID = F)
adm2 = newadm2

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

eastLim = c(-110.6405, 832.4544)
northLim = c(-555.1739, 608.7130)

# Instead of the code in the following if statement, we can just use:
library(SUMMER)
data(kenyaPopulationData)
if(FALSE) {
  ## Construct poppsubKenya, a table of urban/rural general population totals 
  ## in each subarea. Technically, this is not necessary since we can load in 
  ## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
  ## the areas in km^2 of the areas and subareas
  
  library(rgdal)
  library(sp)
  
  # use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
  midLon = mean(adm1@bbox[1,])
  midLat = mean(adm1@bbox[2,])
  p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
               " +lat_0=", midLat, " +units=km")
  
  library(rgdal)
  
  adm1proj <- spTransform(adm1, CRS(p4s))
  adm2proj <- spTransform(adm2, CRS(p4s))
  
  # now calculate spatial area in km^2
  library(rgeos)
  admin1Areas = gArea(adm1proj, TRUE)
  admin2Areas = gArea(adm2proj, TRUE)
  areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
  areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
                             spatialArea=admin2Areas)
  
  # Calculate general population totals at the subarea (Admin-2) x urban/rural 
  # level and using 1km resolution population grid. Assign urbanicity by 
  # thresholding population density based on estimated proportion population 
  # urban/rural, making sure total area (Admin-1) urban/rural populations in 
  # each area matches poppaKenya.
  require(fields)
  # NOTE: the following function will typically take ~20 minutes. Can speed up 
  #       by setting kmRes to be higher, but we recommend fine resolution for 
  #       this step, since it only needs to be done once. Instead of running this, 
  #       you can simply run data(kenyaPopulationData)
  system.time(poppsubKenya <- getPoppsub(
    kmRes=1, pop=pop, domainPoly=kenyaPoly,
    eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
    poppa = poppaKenya, areapa=areapaKenya, areapsub=areapsubKenya, 
    areaMapDat=adm1, subareaMapDat=adm2, 
    areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
  
  # threshold poppsub so that subareas with too small an urban or rural population 
  # don't have any (for mainly computational purposes)
  poppsubKenyaThresh = thresholdPoppsub(poppsubKenya, threshProp = .005)
}

# Now generate a general population integration table at 5km resolution, 
# based on subarea (Admin-2) x urban/rural population totals. This takes 
# ~1 minute
system.time(popMatKenya <- makePopIntegrationTab(
  kmRes=5, pop=pop, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = poppaKenya, poppsub=poppsubKenya, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

system.time(popMatKenyaThresh <- makePopIntegrationTab(
  kmRes=5, pop=pop, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = poppaKenya, poppsub=poppsubKenyaThresh, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

## Adjust popMat to be target (neonatal) rather than general population 
## density. First create the target population frame
## (these numbers are based on IPUMS microcensus data)
mothersPerHouseholdUrb = 0.3497151
childrenPerMotherUrb = 1.295917
mothersPerHouseholdRur = 0.4787696
childrenPerMotherRur = 1.455222
targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * 
  childrenPerMotherUrb
targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * 
  childrenPerMotherRur
easpaKenyaNeonatal = easpaKenya
easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
easpaKenyaNeonatal$popRur = targetPopPerStratumRural
easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + 
  easpaKenyaNeonatal$popRur
easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / 
  easpaKenyaNeonatal$popTotal
easpaKenyaNeonatal$pctTotal = 
  100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)

# Generate the target population density by scaling the current 
# population density grid at the Admin1 x urban/rural level
areaSubareaTab = data.frame(area=adm2@data$NAME_1, subarea = adm2@data$NAME_2)

popMatKenyaNeonatal = adjustPopMat(popMatKenya, easpaKenyaNeonatal)
popMatKenyaNeonatalThresh = adjustPopMat(popMatKenyaThresh, easpaKenyaNeonatal)

# Generate neonatal population table from the neonatal population integration 
# matrix. This is technically not necessary for population simulation purposes, 
# but is here for illustrative purposes
poppsubKenyaNeonatal = poppRegionFromPopMat(popMatKenyaNeonatal, 
                                            popMatKenyaNeonatal$subarea)
poppsubKenyaNeonatal = cbind(subarea=poppsubKenyaNeonatal$region, area=areaSubareaTab$area[match(poppsubKenyaNeonatal$region, areaSubareaTab$subarea)], poppsubKenyaNeonatal[,-1])
poppsubKenyaNeonatalThresh = poppRegionFromPopMat(popMatKenyaNeonatalThresh, 
                                            popMatKenyaNeonatalThresh$subarea)
poppsubKenyaNeonatalThresh = cbind(subarea=poppsubKenyaNeonatalThresh$region, area=areaSubareaTab$area[match(poppsubKenyaNeonatalThresh$region, areaSubareaTab$subarea)], poppsubKenyaNeonatalThresh[,-1])

# save files
save(kenyaPoly, adm1, adm2, kenyaMesh, file="savedOutput/global/kenyaMapData.RData")
save(popMatKenya, popMatKenyaNeonatal, popMatKenyaThresh, popMatKenyaNeonatalThresh, file="savedOutput/global/kenyaPopulationMats.RData")

# make facebook popMat/poppsub files ----
popFace = raster("data/popData/facebook_2020.tif")

if(FALSE) {
  ## Construct poppsubKenya, a table of urban/rural general population totals 
  ## in each subarea. Technically, this is not necessary since we can load in 
  ## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
  ## the areas in km^2 of the areas and subareas
  
  library(rgdal)
  library(sp)
  
  # use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
  midLon = mean(adm1@bbox[1,])
  midLat = mean(adm1@bbox[2,])
  p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
               " +lat_0=", midLat, " +units=km")
  
  library(rgdal)
  
  adm1proj <- spTransform(adm1, CRS(p4s))
  adm2proj <- spTransform(adm2, CRS(p4s))
  
  # now calculate spatial area in km^2
  library(rgeos)
  admin1Areas = gArea(adm1proj, TRUE)
  admin2Areas = gArea(adm2proj, TRUE)
  areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
  areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
                             spatialArea=admin2Areas)
  
  # Calculate general population totals at the subarea (Admin-2) x urban/rural 
  # level and using 1km resolution population grid. Assign urbanicity by 
  # thresholding population density based on estimated proportion population 
  # urban/rural, making sure total area (Admin-1) urban/rural populations in 
  # each area matches poppaKenya.
  require(fields)
  # NOTE: the following function will typically take ~20 minutes. Can speed up 
  #       by setting kmRes to be higher, but we recommend fine resolution for 
  #       this step, since it only needs to be done once. Instead of running this, 
  #       you can simply run data(kenyaPopulationData)
  system.time(out <- makePopIntegrationTab(
    kmRes=1, pop=popFace, domainPoly=kenyaPoly,
    eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
    poppa = poppaKenya, areapa=areapaKenya, areapsub=areapsubKenya, 
    areaMapDat=adm1, subareaMapDat=adm2, returnPoppTables=TRUE, 
    areaNameVar = "NAME_1", subareaNameVar="NAME_2", 
    fixZeroPopDensitySubareas=TRUE))
  
  # make sure we don't get any zero population subareas. If so, substitute in 
  # population totals and density from the WorldPop population estimates 
  poppsubKenyaFace = out$poppsub
  popMatKenyaFace = out$popMat
  
  poppsubKenyaFace = modifyPoppsub(poppsubKenyaFace, poppsubKenya[is.nan(poppsubKenyaFace$pctUrb),])
  
  # system.time(out <- makePopIntegrationTab(
  #   kmRes=1, pop=popFace, poppsub=poppsubKenyaFace, domainPoly=kenyaPoly,
  #   eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  #   poppa = poppaKenya, areapa=areapaKenya, areapsub=areapsubKenya, 
  #   areaMapDat=adm1, subareaMapDat=adm2, returnPoppTables=TRUE, 
  #   areaNameVar = "NAME_1", subareaNameVar="NAME_2", extractChunkSize=1000, 
  #   fixZeroPopDensitySubareas=TRUE))
  # poppsubKenyaFace = out$poppsub
  # popMatKenyaFace = out$popMat
  # 
  # popMatKenyaFace$pop[popMatKenyaFace$subarea %in% poppsubKenyaFace$subarea[is.nan(poppsubKenyaFace$pctUrb)]] =
  #   popMatKenya$pop[popMatKenyaFace$subarea %in% poppsubKenyaFace$subarea[is.nan(poppsubKenyaFace$pctUrb)]]
  # 
  
  # threshold poppsub so that subareas with too small an urban or rural population 
  # don't have any (for mainly computational purposes)
  poppsubKenyaFaceThresh = thresholdPoppsub(poppsubKenyaFace, threshProp = .005)
}

# Now generate a general population integration table at 5km resolution, 
# based on subarea (Admin-2) x urban/rural population totals. This takes 
# ~1 minute

# apparently raster and terra don't use chunked operations for extract function...
# instead we must use ad hoc solution of aggregating facebook population density 
# raster to coarser resolution before making the integration table to avoid memory 
# issues. :(
# defaultRasterOptions = rasterOptions()
# rasterOptions(maxmemory=1e8, chunksize=1e5, todisk=TRUE, memfrac=.1)
# rasterOptions(, todisk=FALSE)
# terraOptions(memfrac=.2, todisk=TRUE, memmax=1)

popFaceAgg = raster::aggregate(popFace, fact=100, fun=sum)
system.time(popMatKenyaFace <- makePopIntegrationTab(
  kmRes=5, pop=popFaceAgg, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = poppaKenya, poppsub=poppsubKenya, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2", 
  fixZeroPopDensitySubareas=TRUE))

system.time(popMatKenyaFaceThresh <- makePopIntegrationTab(
  kmRes=5, pop=popFaceAgg, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = poppaKenya, poppsub=poppsubKenyaFaceThresh, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2", 
  fixZeroPopDensitySubareas=TRUE))

## Adjust popMat to be target (neonatal) rather than general population 
## density. First create the target population frame
## (these numbers are based on IPUMS microcensus data)
mothersPerHouseholdUrb = 0.3497151
childrenPerMotherUrb = 1.295917
mothersPerHouseholdRur = 0.4787696
childrenPerMotherRur = 1.455222
targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * 
  childrenPerMotherUrb
targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * 
  childrenPerMotherRur
easpaKenyaNeonatal = easpaKenya
easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
easpaKenyaNeonatal$popRur = targetPopPerStratumRural
easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + 
  easpaKenyaNeonatal$popRur
easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / 
  easpaKenyaNeonatal$popTotal
easpaKenyaNeonatal$pctTotal = 
  100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)

# Generate the target population density by scaling the current 
# population density grid at the Admin1 x urban/rural level
areaSubareaTab = data.frame(area=adm2@data$NAME_1, subarea = adm2@data$NAME_2)

popMatKenyaFaceNeonatal = adjustPopMat(popMatKenyaFace, easpaKenyaNeonatal)
popMatKenyaFaceNeonatalThresh = adjustPopMat(popMatKenyaFaceThresh, easpaKenyaNeonatal)

# Generate neonatal population table from the neonatal population integration 
# matrix. This is technically not necessary for population simulation purposes, 
# but is here for illustrative purposes
poppsubKenyaFaceNeonatal = poppRegionFromPopMat(popMatKenyaNeonatal, 
                                            popMatKenyaNeonatal$subarea)
poppsubKenyaFaceNeonatal = cbind(subarea=poppsubKenyaFaceNeonatal$region, area=areaSubareaTab$area[match(poppsubKenyaFaceNeonatal$region, areaSubareaTab$subarea)], poppsubKenyaFaceNeonatal[,-1])
poppsubKenyaFaceNeonatalThresh = poppRegionFromPopMat(popMatKenyaFaceNeonatalThresh, 
                                                  popMatKenyaFaceNeonatalThresh$subarea)
poppsubKenyaFaceNeonatalThresh = cbind(subarea=poppsubKenyaFaceNeonatalThresh$region, area=areaSubareaTab$area[match(poppsubKenyaFaceNeonatalThresh$region, areaSubareaTab$subarea)], poppsubKenyaFaceNeonatalThresh[,-1])

# save files
save(popMatKenyaFace, popMatKenyaFaceNeonatal, popMatKenyaFaceThresh, popMatKenyaFaceNeonatalThresh, file="savedOutput/global/kenyaFacePopulationMats.RData")
save(poppsubKenyaFace, poppsubKenyaFaceThresh, poppsubKenyaFaceNeonatal, poppsubKenyaFaceNeonatalThresh, file="savedOutput/global/poppsubFace.RData")

# Make 2019 census popMat/poppsub files ----
source("code/modAgg2.R")
easpa2019 = makeEASPA2019(neonatal=FALSE)
easpa2019$pctUrb = 100 * easpa2019$popUrb / easpa2019$popTotal
easpa2019$pctTotal = 100 * easpa2019$popTotal / sum(easpa2019$popTotal)
easpa2019Neonatal = makeEASPA2019(neonatal=TRUE)
easpa2019Neonatal$pctUrb = 100 * easpa2019Neonatal$popUrb / easpa2019Neonatal$popTotal
easpa2019Neonatal$pctTotal = 100 * easpa2019Neonatal$popTotal / sum(easpa2019Neonatal$popTotal)

if(FALSE) {
  ## Construct poppsubKenya, a table of urban/rural general population totals 
  ## in each subarea. Technically, this is not necessary since we can load in 
  ## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
  ## the areas in km^2 of the areas and subareas
  
  library(rgdal)
  library(sp)
  
  # use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
  midLon = mean(adm1@bbox[1,])
  midLat = mean(adm1@bbox[2,])
  p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
               " +lat_0=", midLat, " +units=km")
  
  library(rgdal)
  
  adm1proj <- spTransform(adm1, CRS(p4s))
  adm2proj <- spTransform(adm2, CRS(p4s))
  
  # now calculate spatial area in km^2
  library(rgeos)
  admin1Areas = gArea(adm1proj, TRUE)
  admin2Areas = gArea(adm2proj, TRUE)
  areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
  areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
                             spatialArea=admin2Areas)
  
  # Calculate general population totals at the subarea (Admin-2) x urban/rural 
  # level and using 1km resolution population grid. Assign urbanicity by 
  # thresholding population density based on estimated proportion population 
  # urban/rural, making sure total area (Admin-1) urban/rural populations in 
  # each area matches poppaKenya.
  require(fields)
  # NOTE: the following function will typically take ~20 minutes. Can speed up 
  #       by setting kmRes to be higher, but we recommend fine resolution for 
  #       this step, since it only needs to be done once. Instead of running this, 
  #       you can simply run data(kenyaPopulationData)
  system.time(poppsubKenya2019 <- getPoppsub(
    kmRes=1, pop=pop, domainPoly=kenyaPoly,
    eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
    poppa = easpa2019, areapa=areapaKenya, areapsub=areapsubKenya, 
    areaMapDat=adm1, subareaMapDat=adm2, 
    areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
  
  # threshold poppsub so that subareas with too small an urban or rural population 
  # don't have any (for mainly computational purposes)
  poppsubKenya2019Thresh = thresholdPoppsub(poppsubKenya2019, threshProp = .005)
}

# Now generate a general population integration table at 5km resolution, 
# based on subarea (Admin-2) x urban/rural population totals. This takes 
# ~1 minute
system.time(popMatKenya2019 <- makePopIntegrationTab(
  kmRes=5, pop=pop, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = easpa2019, poppsub=poppsubKenya2019, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

system.time(popMatKenya2019Thresh <- makePopIntegrationTab(
  kmRes=5, pop=pop, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = easpa2019, poppsub=poppsubKenya2019Thresh, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

# Generate the target population density by scaling the current 
# population density grid at the Admin1 x urban/rural level
areaSubareaTab = data.frame(area=adm2@data$NAME_1, subarea = adm2@data$NAME_2)

popMatKenya2019Neonatal = adjustPopMat(popMatKenya2019, easpa2019Neonatal)
popMatKenya2019NeonatalThresh = adjustPopMat(popMatKenya2019Thresh, easpa2019Neonatal)

# Generate neonatal population table from the neonatal population integration 
# matrix. This is technically not necessary for population simulation purposes, 
# but is here for illustrative purposes
poppsubKenya2019Neonatal = poppRegionFromPopMat(popMatKenya2019Neonatal, 
                                            popMatKenya2019Neonatal$subarea)
poppsubKenya2019Neonatal = cbind(subarea=poppsubKenya2019Neonatal$region, area=areaSubareaTab$area[match(poppsubKenya2019Neonatal$region, areaSubareaTab$subarea)], poppsubKenya2019Neonatal[,-1])
poppsubKenya2019NeonatalThresh = poppRegionFromPopMat(popMatKenya2019NeonatalThresh, 
                                                  popMatKenya2019NeonatalThresh$subarea)
poppsubKenya2019NeonatalThresh = cbind(subarea=poppsubKenya2019NeonatalThresh$region, area=areaSubareaTab$area[match(poppsubKenya2019NeonatalThresh$region, areaSubareaTab$subarea)], poppsubKenya2019NeonatalThresh[,-1])

# save files
save(popMatKenya2019, popMatKenya2019Neonatal, popMatKenya2019Thresh, popMatKenya2019NeonatalThresh, file="savedOutput/global/kenya2019PopulationMats.RData")
save(poppsubKenya2019, poppsubKenya2019Thresh, poppsubKenya2019Neonatal, poppsubKenya2019NeonatalThresh, file="savedOutput/global/poppsub2019.RData")


# Make jittered census popMat/poppsub files ----
source("code/modAgg2.R")
easpaJittered = makeEASPAJittered(neonatal=FALSE)
easpaJittered$pctUrb = 100 * easpaJittered$popUrb / easpaJittered$popTotal
easpaJittered$pctTotal = 100 * easpaJittered$popTotal / sum(easpaJittered$popTotal)
easpaJitteredNeonatal = makeEASPAJittered(neonatal=TRUE)
easpaJitteredNeonatal$pctUrb = 100 * easpaJitteredNeonatal$popUrb / easpaJitteredNeonatal$popTotal
easpaJitteredNeonatal$pctTotal = 100 * easpaJitteredNeonatal$popTotal / sum(easpaJitteredNeonatal$popTotal)

if(FALSE) {
  ## Construct poppsubKenya, a table of urban/rural general population totals 
  ## in each subarea. Technically, this is not necessary since we can load in 
  ## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
  ## the areas in km^2 of the areas and subareas
  
  library(rgdal)
  library(sp)
  
  # use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
  midLon = mean(adm1@bbox[1,])
  midLat = mean(adm1@bbox[2,])
  p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
               " +lat_0=", midLat, " +units=km")
  
  library(rgdal)
  
  adm1proj <- spTransform(adm1, CRS(p4s))
  adm2proj <- spTransform(adm2, CRS(p4s))
  
  # now calculate spatial area in km^2
  library(rgeos)
  admin1Areas = gArea(adm1proj, TRUE)
  admin2Areas = gArea(adm2proj, TRUE)
  areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
  areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
                             spatialArea=admin2Areas)
  
  # Calculate general population totals at the subarea (Admin-2) x urban/rural 
  # level and using 1km resolution population grid. Assign urbanicity by 
  # thresholding population density based on estimated proportion population 
  # urban/rural, making sure total area (Admin-1) urban/rural populations in 
  # each area matches poppaKenya.
  require(fields)
  # NOTE: the following function will typically take ~20 minutes. Can speed up 
  #       by setting kmRes to be higher, but we recommend fine resolution for 
  #       this step, since it only needs to be done once. Instead of running this, 
  #       you can simply run data(kenyaPopulationData)
  system.time(poppsubKenyaJittered <- getPoppsub(
    kmRes=1, pop=pop, domainPoly=kenyaPoly,
    eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
    poppa = easpaJittered, areapa=areapaKenya, areapsub=areapsubKenya, 
    areaMapDat=adm1, subareaMapDat=adm2, 
    areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
  
  # threshold poppsub so that subareas with too small an urban or rural population 
  # don't have any (for mainly computational purposes)
  poppsubKenyaJitteredThresh = thresholdPoppsub(poppsubKenyaJittered, threshProp = .005)
}

# Now generate a general population integration table at 5km resolution, 
# based on subarea (Admin-2) x urban/rural population totals. This takes 
# ~1 minute
system.time(popMatKenyaJittered <- makePopIntegrationTab(
  kmRes=5, pop=pop, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = easpaJittered, poppsub=poppsubKenyaJittered, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

system.time(popMatKenyaJitteredThresh <- makePopIntegrationTab(
  kmRes=5, pop=pop, domainPoly=kenyaPoly,
  eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
  poppa = easpaJittered, poppsub=poppsubKenyaJitteredThresh, 
  areaMapDat = adm1, subareaMapDat = adm2,
  areaNameVar = "NAME_1", subareaNameVar="NAME_2"))

# Generate the target population density by scaling the current 
# population density grid at the Admin1 x urban/rural level
areaSubareaTab = data.frame(area=adm2@data$NAME_1, subarea = adm2@data$NAME_2)

popMatKenyaJitteredNeonatal = adjustPopMat(popMatKenyaJittered, easpaJitteredNeonatal)
popMatKenyaJitteredNeonatalThresh = adjustPopMat(popMatKenyaJitteredThresh, easpaJitteredNeonatal)

# Generate neonatal population table from the neonatal population integration 
# matrix. This is technically not necessary for population simulation purposes, 
# but is here for illustrative purposes
poppsubKenyaJitteredNeonatal = poppRegionFromPopMat(popMatKenyaJitteredNeonatal, 
                                                popMatKenyaJitteredNeonatal$subarea)
poppsubKenyaJitteredNeonatal = cbind(subarea=poppsubKenyaJitteredNeonatal$region, area=areaSubareaTab$area[match(poppsubKenyaJitteredNeonatal$region, areaSubareaTab$subarea)], poppsubKenyaJitteredNeonatal[,-1])
poppsubKenyaJitteredNeonatalThresh = poppRegionFromPopMat(popMatKenyaJitteredNeonatalThresh, 
                                                      popMatKenyaJitteredNeonatalThresh$subarea)
poppsubKenyaJitteredNeonatalThresh = cbind(subarea=poppsubKenyaJitteredNeonatalThresh$region, area=areaSubareaTab$area[match(poppsubKenyaJitteredNeonatalThresh$region, areaSubareaTab$subarea)], poppsubKenyaJitteredNeonatalThresh[,-1])

# save files
save(popMatKenyaJittered, popMatKenyaJitteredNeonatal, popMatKenyaJitteredThresh, popMatKenyaJitteredNeonatalThresh, file="savedOutput/global/kenyaJitteredPopulationMats.RData")
save(poppsubKenyaJittered, poppsubKenyaJitteredThresh, poppsubKenyaJitteredNeonatal, poppsubKenyaJitteredNeonatalThresh, file="savedOutput/global/poppsubJittered.RData")

