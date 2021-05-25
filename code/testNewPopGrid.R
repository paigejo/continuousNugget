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

##### generate poppcon with the new constituencies
# for precomputating populations of constituencies possibly crossed with urban/rural
popGridFine = makeInterpPopGrid(kmRes=1, mean.neighbor=500, delta=.05, conMap=adm2)
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

# make coarse popGrid with custom points for constituencies without enough pixels
popGridNew = makeInterpPopGrid(mean.neighbor=50, delta=.1, conMap=adm2, poppcon=poppcon)

dim(popGrid)
dim(popGridNew)

##### tests ----
# check how many pixels there are per area
lonLatGrid = cbind(popGrid$lon, popGrid$lat)
constituencies = getConstituency(lonLatGrid, mean.neighbor=50, delta=.1, conMap=adm2)$constituencyNames
constituencies = factor(constituencies, levels=sort(adm2@data$CONSTITUEN))
constituencies2 = popGridNew$admin2
out = aggregate(constituencies2, by=list(constituency=constituencies2), FUN=length, drop=FALSE)
out

lonLatGrid = cbind(popGridNew$lon, popGridNew$lat)
constituencies = getConstituency(lonLatGrid, mean.neighbor=50, delta=.1, conMap=adm2)$constituencyNames
constituencies = factor(constituencies, levels=sort(adm2@data$CONSTITUEN))
constituencies2 = popGridNew$admin2
out2 = aggregate(constituencies2, by=list(constituency=constituencies2), FUN=length, drop=FALSE)
out2

cbind(out, out2$x)
nas = is.na(out$x)
onePixel = out$x == 1
cbind(out, out2$x)[nas | onePixel,]
#         constituency  x out2$x
# 35   Dagoretti North  1      1
# 36   Dagoretti South  1      1
# 39  Embakasi Central NA      1
# 41    Embakasi North NA      1
# 42    Embakasi South  1      1
# 43     Embakasi West  1      1
# 70             Jomvu  1      1
# 83         Kamukunji NA      1
# 100            Kibra  1      1
# 158         Makadara NA      1
# 172          Mathare NA      1
# 194            Mvita  1      1
# 219            Nyali  1      1
# 235          Ruaraka NA      1

# make sure every constituency has at least one urban and rural pixel if necessary. 
# Not necessary: all of the bad constituencies were entirely urban, and each is 
# given 1 pixel as shown in previous test
# hasUrbanPop = poppcon$popUrb > 0
# hasRuralPop = poppcon$popRur > 0

# determine which constituencies need to be fixed
nas = is.na(out$x)
naConstituencies = out$constituency[nas]
dataNAConstituencies = which(adm2@data$CONSTITUEN %in% as.character(naConstituencies))

# get pixels per constituency for the fine grid
library(SUMMER)
lonLatGridFine = generatePixelGrid(1, kenyaPoly, eastLim, northLim, projKenya)
lonLatGridFine = lonLatGridFine[,1:2]
constituenciesFine = getConstituency(lonLatGridFine, mean.neighbor=500, delta=.05)$constituencyNames
constituenciesFine = factor(constituenciesFine, levels=sort(adm2@data$CONSTITUEN))
outFine = aggregate(constituenciesFine, by=list(constituency=constituenciesFine), FUN=length, drop=FALSE)
outFine
min(outFine$x) # 2

# calculate constituency area (this step is not necessary)
areas = getArea("Constituency", adm2, sortAreas=TRUE)
cbind(outFine, area=areas)

# calculate constituency area again, but using number of fine grid cells, 
# since each fine scale grid cell is 1km^2
fineAreas = outFine$x

# set coarse pixel areas
coarsePixelAreas = rep(25, nrow(popGrid))

# calculate constituency centroids for the problem constituencies
thisSpatialPolyList = as.SpatialPolygons.PolygonsList(adm2@polygons)
centroids = matrix(ncol=2, nrow=nrow(adm2@data))
for(i in 1:nrow(adm2@data)) {
  thisCentroid = coordinates(thisSpatialPolyList[i])
  centroids[i,] = thisCentroid
}
naCentroids = centroids[dataNAConstituencies,]

# 




