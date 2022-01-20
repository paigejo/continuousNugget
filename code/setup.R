# install requied packages if necessary
if(FALSE) {
  install.packages("Matrix")
  install.packages("spam")
  install.packages("fields")
  install.packages("LatticeKrig")
  install.packages("invgamma")
  install.packages("latex2exp")
  install.packages("xtable")
  install.packages("profvis")
  install.packages("geosphere")
  install.packages("viridis")
  install.packages("sp")
  install.packages("raster")
  install.packages("MCMCpack")
  install.packages("numDeriv")
  install.packages("INLA")
  install.packages("edfun") # for drawing from empirical distributions quickly
  install.packages("data.table")
  install.packages("sampling")
  install.packages("haven")
  install.packages("survey")
  install.packages("abind")
  install.packages("devtools")
  install.packages("ggplot2")
  
  library(devtools)
  install_github("https://github.com/richardli/SUMMER/tree/dev")
}

# load required packages and R scripts
library(Matrix)
library(spam)
library(fields)
library(LatticeKrig)
library(invgamma)
library(latex2exp)
library(xtable)
library(profvis)
library(geosphere)
library(viridis)
library(sp)
library(raster)
library(MCMCpack)
library(numDeriv)
library(INLA)
library(edfun) # for drawing from empirical distributions quickly
library(data.table)
library(sampling)
library(haven)
library(survey)
library(abind)
# install_github("https://github.com/richardli/SUMMER/tree/dev")
library(SUMMER)
# library(Rcpp)
library(ggplot2)

codeDirectory <<- "~/git/continuousNugget/code/"
figDirectory <<- "~/git/continuousNugget/figures/"
dataDirectory <<- "~/git/continuousNugget/data/"
outputDirectory <<- "~/git/continuousNugget/savedOutput/"
globalDirectory <<- "~/git/continuousNugget/savedOutput/global/"
elkDirectory <<- "~/git/LK-INLA/"

inf = sessionInfo()
if(inf$platform == "x86_64-apple-darwin17.0 (64-bit)") {
  setwd("~/git/continuousNugget/")
  options(error=recover)
} else if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)" && inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
  setwd("~/git/continuousNugget/")
} else if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/continuousNugget/")
  options(error=recover)
} else {
  setwd("~/git/continuousNugget/")
  inla.setOption(num.threads=1) # consider raising
  options(error=recover)
}

setwd("~/git/continuousNugget")
# source(paste0(elkDirectory, 'LKinla.R'))
# source(paste0(elkDirectory, 'LKinla_rgeneric.R'))
source('code/modSPDE.R')
source('code/modELK.R')
# source('code/modAgg.R')
source('code/modAgg2.R')
source("code/scores.R")
source('code/compareModels.R')
source("code/utilityFuns.R")
source('code/getSimulationDataSets.R')
# source('code/generateSimDataSets.R')
source('code/generateSimDataSets2.R')
source('code/simStudy.R')
source('code/test.R')
source('code/plotGenerator.R')
source('code/modDirect.R')
source('code/nmr.R')
# Rcpp::sourceCpp("code/Rcpp/rmultinomProbMat_rcpp.cpp")

# set some basic parameters to hold constant throughout analysis
kenyaLonRange = c(33.5, 42)
kenyaLatRange = c(-5, 5.5)
kenyaLonLength = kenyaLonRange[2] - kenyaLonRange[1]
kenyaLatLength = kenyaLatRange[2] - kenyaLatRange[1]
totalKenyaPop = 43*10^6 # from DHS 2014 survey final report page 2
totalRows = 4320 # latitude
totalCols = 8640 # longitude
increaseFac=1
resPerDeg = 24
extentCols = round(kenyaLonLength*resPerDeg)
extentRows = round(kenyaLatLength*resPerDeg)
lonsInterp = seq(kenyaLonRange[1], kenyaLonRange[2], l=round(extentCols*increaseFac))
latsInterp = seq(kenyaLatRange[1], kenyaLatRange[2], l=round(extentRows*increaseFac))
lonRes = lonsInterp[2] - lonsInterp[1]
latRes = latsInterp[2] - latsInterp[1]

## load in global variables made from the following script: 
if(FALSE) {
  source('code/preprocess.R')
}

# load(paste0(globalDirectory, "adminMapData.RData"))
load(paste0(globalDirectory, "kenyaMapData.RData"))

# county to region mapping
ctp = load(paste0(globalDirectory, "kenya-prov-county-map.RData"))

# number of EAs and clusters per county/region
numEAs = 96251
load(paste0(globalDirectory, "easpc.RData"))
load(paste0(globalDirectory, "easpr.RData"))
load(paste0(globalDirectory, "clustpc.RData"))
clustpc$area = clustpc$County
load(paste0(globalDirectory, "clustpr.RData"))
load(paste0(globalDirectory, "poppc.RData"))
load(paste0(globalDirectory, "poppr.RData"))
# load(paste0(globalDirectory, "poppcon.RData"))
load(paste0(globalDirectory, "poppcMort.RData"))
load(paste0(globalDirectory, "easpcMort.RData"))
load(paste0(globalDirectory, "poppcAdjusted.RData"))
# load(paste0(globalDirectory, "poppconAdjusted.RData"))
load(paste0(globalDirectory, "empiricalDistributions.RData"))

data(kenyaPopulationData)
githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
                    "popMatKenya5km.RData?raw=true")
tempDirectory = "~/"
popMatFilename = paste0(tempDirectory, "/popMatKenya5km.RData")
if(!file.exists(popMatFilename)) {
  download.file(githubURL,popMatFilename)
}

# load it in
out = load(popMatFilename)
out

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
spdeMesh = kenyaMesh

# cleanup adm2 and adm1 maps
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

poppsubKenyaNeonatal = SUMMER:::poppRegionFromPopMat(popMatKenyaNeonatal, 
                                            popMatKenyaNeonatal$subarea)
poppsubKenyaNeonatal = 
  cbind(subarea=poppsubKenyaNeonatal$region, 
        area=adm2@data$NAME_1[match(poppsubKenyaNeonatal$region, adm2@data$NAME_2)], 
        poppsubKenyaNeonatal[,-1])

poppcon = poppsubKenya
poppconAdjusted = poppsubKenyaNeonatal

# load enumeration areas and neonatal mortality data
load(paste0(globalDirectory, "kenyaEAs.RData"))
load(paste0(globalDirectory, "kenyaData.RData"))

# load 5km population density grid (of neonatals) over Kenya
# load(paste0(globalDirectory, "popGrid.RData"))
# load(paste0(globalDirectory, "popGridAdjusted.RData"))
load(paste0(globalDirectory, "kenyaPopulationMats.RData"))
popGrid = popMatKenya
popGridAdjusted = popMatKenyaNeonatal
load(paste0(globalDirectory, "popGridCoarse.RData"))
load(paste0(globalDirectory, "popGridCoarseAdjusted.RData"))

# load limits
load(paste0(globalDirectory, "lims.RData"))

# parallelization
if(!exists("doParallel") || (exists("doParallel") && doParallel == FALSE)) {
  assign("doParallel", FALSE, envir=.GlobalEnv)
  assign("cores", NULL, envir=.GlobalEnv)
  assign("cl", NULL, envir=.GlobalEnv)
}

# determine version of PROJ.4
ver = rgdal::rgdal_extSoftVersion()
theseNames = names(ver)
thisI = which(grepl("PROJ", theseNames))
PROJ6 <- as.numeric(substr(ver[thisI], 1, 1)) >= 6




