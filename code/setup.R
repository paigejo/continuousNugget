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

codeDirectory <<- "~/git/continuousNugget/code/"
figDirectory <<- "~/git/continuousNugget/figures/"
dataDirectory <<- "~/git/continuousNugget/data/"
outputDirectory <<- "~/git/continuousNugget/savedOutput/"
globalDirectory <<- "~/git/continuousNugget/savedOutput/global/"
elkDirectory <<- "~/git/LK-INLA/"

inf = sessionInfo()
if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)" && inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
  setwd("~/git/continuousNugget/")
} else if(inf$platform != "x86_64-w64-mingw32/x64 (64-bit)" && inf$platform != "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/continuousNugget/")
  options(error=recover)
} else if(inf$platform == "x86_64-pc-linux-gnu (64-bit)") {
  setwd("~/git/continuousNugget/")
  inla.setOption(num.threads=1) # consider raising
  options(error=recover)
} else {
  setwd("U:/git/continuousNugget/")
  inla.setOption(num.threads=1)
  options(error=recover)
}

setwd("~/git/continuousNugget")
source(paste0(elkDirectory, 'LKinla.R'))
source(paste0(elkDirectory, 'LKinla_rgeneric.R'))
source('code/modSPDE.R')
source('code/modLKinla.R')
source("code/scores.R")
source('code/compareModels.R')
source("code/utilityFuns.R")
source('code/getSimulationDataSets.R')
source('code/simStudy.R')
# source('code/plotGenerator.R')

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

load(paste0(globalDirectory, "adminMapData.RData"))

# county to region mapping
ctp = load(paste0(globalDirectory, "kenya-prov-county-map.RData"))

# number of EAs and clusters per county/region
numEAs = 96251
load(paste0(globalDirectory, "easpc.RData"))
load(paste0(globalDirectory, "easpr.RData"))
load(paste0(globalDirectory, "clustpc.RData"))
load(paste0(globalDirectory, "clustpr.RData"))
load(paste0(globalDirectory, "poppc.RData"))
load(paste0(globalDirectory, "poppr.RData"))

# set enumeration areas
load(paste0(globalDirectory, "kenyaEAs.RData"))

# load 5km population density grid (of neonatals) over Kenya
load(paste0(globalDirectory, "popGridAdjusted.RData"))

# parallelization
if(!exists("doParallel") || (exists("doParallel") && doParallel == FALSE)) {
  assign("doParallel", FALSE, envir=.GlobalEnv)
  assign("cores", NULL, envir=.GlobalEnv)
  assign("cl", NULL, envir=.GlobalEnv)
}