source("setup.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))

jobInds = getJobIndices(index, rev=TRUE)
i = jobInds[1]
j = jobInds[2]

load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
argList = spde_prevRiskSimStudyCommandArgs[[i]]
argList$j = j

do.call("runSimStudy", argList)
stopCluster(cl)