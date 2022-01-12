# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=20000 --pty bash

source("setup.R")
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == 700

jobInds = getJobIndices(index, rev=TRUE)
i = jobInds[1]
j = jobInds[2]

# load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
# argList = spde_prevRiskSimStudyCommandArgs[[i]]
# argList$j = j

system.time(invisible(runSimStudyij(i, j)))