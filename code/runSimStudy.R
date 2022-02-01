# srun --partition=CPUQ --time=02:00:00 --mem-per-cpu=20000 --pty bash

source("setup.R")
index = as.numeric(commandArgs(trailingOnly = TRUE)) # test with index == 700
coarse=TRUE

jobInds = getJobIndices(index, rev=TRUE)
i = jobInds[1]
j = jobInds[2]

# load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
# argList = spde_prevRiskSimStudyCommandArgs[[i]]
# argList$j = j

# Rprof("savedOutput/simStudyResults/tempFiles/data.Rprof", interval = 0.01, line.profiling = TRUE,
#       gc.profiling = TRUE, memory.profiling = TRUE)

source("code/simPopTest.R")
# p = profvis({
system.time(out <- runSimStudyij(i, j, coarse=coarse, doGC=TRUE))
# })
# save(p, file="savedOutput/simStudyResults/tempFiles/profFile.RData")

# Rprof(NULL)
# profvis(prof_input = "data.Rprof")