source("setup.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("savedOutput/simDataSets/spde_lcpbCommandArgs.RData")
argList = spde_lcpbCommandArgs[[index]]
invisible(do.call("resultsSPDE_LCPB", argList))
