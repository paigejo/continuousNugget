library("rgdal", lib.loc="/sw/ebpkgs/software/rgdal/1.4-8-foss-2019b-R-3.6.2")
source("setup.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("savedOutput/simDataSets/spde_lcpbCommandArgs.RData")
argList = spde_lcpbCommandArgs[[index]]
invisible(do.call("resultsSPDE_LCPB", argList))
