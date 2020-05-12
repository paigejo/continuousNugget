# an update of the readDat.R script with cleaner code for first year mortality


# read in the STATA data:
# DHS recode manual:
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# ?:
# https://dhsprogram.com/pubs/pdf/SAR8/SAR8.pdf
# final report:
# https://dhsprogram.com/pubs/pdf/FR308/FR308.pdf
library(haven)
library(fields)
library(zoo)
library(latex2exp)
library(maptools)
library(data.table)

wd = getwd()
setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")

# lines represent births
data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

# extract the columns of interest
# b1 - month of birth of child,
# b2 - year of birth of child, 
# b4 - gender of child, # TODO: important?
# b5 - child alive at time of interview, 
# b7 - age of child at death in months completed, 
# v024 - region of residence
# v025 - type of place of residence urban rural (32.3% of pop. is urban in 2009 and 2014, see FR p.30/2)
# v001 - cluster number (in theory there are 1,612 total.  In actuality there are 1,593)
# v002 - household number
# v005 - sample weight out of 1000000. normalized weights sum to number of households
# v006 - month of interview
# V007 - year of interview
# V136 - total number of household members
# V138 - number of women in the household aged 15 to 49
# V137 - number of children in the household aged 5 or under
subdata <- data.frame(data[,c('b2', 'b1', 'b5', 'b7', 'v024', 'v025', 'v001', 'v002', 'v005', 'v006', 'v007')])

# get first month mortality rate for urban and rural areas
urbanData = subdata[subdata$v025 == 1,]
ruralData = subdata[subdata$v025 != 1,]
diedUrban = urbanData[!is.na(urbanData$b7),]
diedRural = ruralData[!is.na(ruralData$b7),]
# 
sum(diedUrban$b7 == 0) / nrow(urbanData)
sum(diedRural$b7 == 0) / nrow(ruralData)
sum(diedUrban$v005[diedUrban$b7 == 0]) / sum(urbanData$v005)
sum(diedRural$v005[diedRural$b7 == 0]) / sum(ruralData$v005)

# get first month mortality rate for urban and rural areas (try the other version of this variable)
urbanData = data[data$v140 == 1,]
ruralData = data[data$v140 != 1,]
diedUrban = urbanData[!is.na(urbanData$b7),]
diedRural = ruralData[!is.na(ruralData$b7),]
# naïve empirical rate
sum(diedUrban$b7 == 0) / nrow(urbanData)
sum(diedRural$b7 == 0) / nrow(ruralData)
# weighted empirical rate
sum(diedUrban$v005[diedUrban$b7 == 0]) / sum(urbanData$v005)
sum(diedRural$v005[diedRural$b7 == 0]) / sum(ruralData$v005)

# check to make sure the right proportion is urban (which is about 30%)
mean(data$v140 == 1)
mean(data$v025 == 1)

# extract births in the range 2005 to 2010 (most recent five years)
lowYear <- 2005
highYear <- 2009
subdata <- subdata[(subdata[,'b2'] >= lowYear & subdata[,'b2'] <= highYear),]

# add a column for the stratification variable as an interaction between
# the urban/rural indicator 'v025' (1: urban, 2:rural) and the region indicator 'v024'
subdata$regionUral <- with(subdata, interaction(v024, v025), drop=TRUE)

# add a column for the unique households with interaction between
# the household indicator 'v002' and the cluster indicator 'v001'
subdata$hhold <- with(subdata, interaction(v001, v002), drop=TRUE)

# find for each cluster the regionUral indicator
clStrat = subdata[,c("v001", "regionUral", "v024", "v025", "v005")]
clStrat = clStrat[!duplicated(clStrat), ]
colnames(clStrat) = c("clusterID", "regionRural", "region", "urban", "samplingWeight")
clStrat$urban =  clStrat$urban == 1

# determine whether each child survived for at least one month
lived = subdata$b7 > 0
lived[is.na(lived)] = TRUE
died = !lived

# get the number of birth by cluster
n <- table(subdata[,'v001'])
clusterid <- dimnames(n)[[1]]
n.data = data.frame(clusterID=clusterid, n=as.vector(n))

# get the number of deaths by cluster
# y <- table(subdata[,'v001'])
y = aggregate(died, list(subdata[,'v001']), sum)
n.data$y = y$x

# add in strata
mort <- merge(n.data, clStrat, by='clusterID', all=TRUE, sort=TRUE)

# Read geographical information
library(rgdal)
spObj = readOGR(dsn = "Kenya2014gps/", layer = "KEGE71FL")

# Extract (lon, lat) coordinates of all clusters
geoObj = data.frame(cId = spObj$DHSCLUST, lon = spObj$LONGNUM, lat = spObj$LATNUM)

# Extract coordinates of clusters with data
idx = match(mort$clusterID, geoObj$cId)
mort$lon = geoObj$lon[idx]
mort$lat = geoObj$lat[idx]

# Missing geographical information is assigned value (0,0)
# Remove these
missIdx = which(mort$lon == 0)
mort = mort[-missIdx,]

library(SUMMER)
library(foreign)

gpsDat = readShapePoints("Kenya2014gps/KEGE71FL.shp")
coords = attr(gpsDat, "coords")
plot(coords)
world(add=TRUE)
test = coords[,1] < 20 #  these are the observations  whose source is missing.  remove these
sum(test)

#  remove observations with unknown locations and set unknown data to NA
gpsDat=gpsDat[!test,]
names(gpsDat)=c("ID","countryID", "year", "clustID", "countryIDFIPS", "countryIDAdminID", 
                "admin1FIPS","admin1IDSALB", 
                "admin1SALB", "admin1ID", "admin1", "regionID", "region", 
                "source", "urban", "lat", "lon", "altGPS", "altRadar", "coordRef")
gpsDat$altGPS[gpsDat$altGPS == 9999] = NA
gpsDat$altGPS[gpsDat$altRadar == 9999] = NA
gpsDat$urban = gpsDat$urban == "U"
gpsDat$countryIDFIPS[gpsDat$countryIDFIPS == "NULL"] = NA
gpsDat$admin1FIPS[gpsDat$admin1FIPS == "NULL"] = NA
gpsDat$admin1IDSALB[gpsDat$admin1IDSALB == "NULL"] = NA
gpsDat$admin1SALB[gpsDat$admin1SALB == "NULL"] = NA
gpsDat$countryIDAdminID[gpsDat$countryIDAdminID == "NULL"] = NA

setwd(wd)
# save(gpsDat, file=paste0(dataDirectory, "gpsDat.RData"))

# get region and admin data from gps data, add to clusters in mort dataset
gpsI = match(data.frame(rbind(mort$lon, mort$lat)), data.frame(rbind(gpsDat$lon, gpsDat$lat)))
mort$admin1 = gpsDat$admin1[gpsI]
mort$region = gpsDat$region[gpsI]

# get easting and northing using projection
tmp = projKenya(mort$lon, mort$lat)
mort$east = tmp[,1]
mort$north = tmp[,2]
save(mort, file=paste0(globalDirectory, "kenyaData.RData"))

# find average urban and rural neonatal mortality rates
sum(mort$y[mort$urban])/sum(mort$n[mort$urban])
sum(mort$y[!mort$urban])/sum(mort$n[!mort$urban])

# make an extremely rudimentary linear regression model, examining the urban effect
perc = mort$y/mort$n
mod = lm(perc ~ urban + admin1, data=mort)
summary(mod)
