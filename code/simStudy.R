# simulation study scripts for UM5R study

# for kmres=2.5: nx=350, ny=430
makeUrbanMap = function(popGrid=NULL, kmres=2.5, savePlot=FALSE, fileName=ifelse(whiteRural, "figures/UrbanMapWhiteRural.png", "figures/urbanMap.png"), 
                        nx=850, ny=1050, width=500, height=500, lonLim=kenyaLonRange, latLim=kenyaLatRange, whiteRural=TRUE) {
  # get prediction locations from population grid
  if(is.null(popGrid)) {
    if(kmres == 5)
      load("popGrid.RData")
    else
      popGrid = makeInterpPopGrid(kmres)
  }
  
  # determine which points in Kenya are urban
  threshes = setThresholds2()
  popThreshes = sapply(1:nrow(popGrid), function(i) {threshes$threshes[threshes$counties == popGrid$admin1[i]]})
  urban = popGrid$popOrig > popThreshes
  
  if(savePlot) {
    png(file=fileName, width=width, height=height)
    par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  }
  plot(popGrid$lon, popGrid$lat, xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1, type="n")
  # quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, 
  #            xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1)
  if(whiteRural)
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c(rgb(0, 0, 0, 0), "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
  else {
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
    legend("bottomleft", c("urban", "rural"), col=c("blue", "green"), pch=19)
  }
  # world(add=TRUE)
  plotMapDat(adm1)
  if(savePlot) {
    dev.off()
  }
}

# simulate enumeration areas from population data.  stratified by urban 
# and rural and county
simEAs = function(kenyaPop, numEAs=96251, totalKenyaPop=43.0 * 10^6, seed=123, 
                  sampleByPop=TRUE, fixNumUrbanAtTruth=TRUE, thisEasPC=easpc) {
  set.seed(seed)
  
  # determine which points in Kenya are urban
  threshes = setThresholds()
  popThreshes = sapply(1:nrow(kenyaPop), function(i) {threshes$threshes[threshes$counties == kenyaPop$admin1[i]]})
  urban = kenyaPop$popOrig > popThreshes
  
  counties = unique(kenyaPop$admin1)
  allSamples = c()
  allUrban = c()
  for(i in 1:length(counties)) {
    countyName = counties[i]
    countyI = kenyaPop$admin1 == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # urban samples
    thisCountyI = countyI & urban
    if(sampleByPop)
      probs = thisCountyI * kenyaPop$popOrig
    else
      probs = thisCountyI
    probs = probs/sum(probs)
    
    # number of rural and urban clusters to sample: either fix it or choose it 
    # randomly based on sample probabilities
    if(fixNumUrbanAtTruth) {
      numUrban = easpc$EAUrb[easpc$County == countyName]
      numRural = easpc$EARur[easpc$County == countyName]
    }
    else {
      numSamplesTotal = easpc$EATotal[easpc$County == countyName]
      countyProbs = countyI
      if(sampleByPop)
        countyProbs = countyProbs * kenyaPop$popOrig
      countyProbs = countyProbs/sum(countyProbs)
      sampleIs = sample(1:nrow(kenyaPop), numSamplesTotal, TRUE, prob=countyProbs)
      numUrban = sum(urban[sampleIs])
      numRural = numSamplesTotal - numUrban
    }
    allSamples = c(allSamples, sample(1:nrow(kenyaPop), numUrban, TRUE, prob=probs))
    allUrban = c(allUrban, rep(TRUE, numUrban))
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      if(sampleByPop)
        probs = thisCountyI * kenyaPop$popOrig
      else
        probs = thisCountyI
      probs = probs/sum(probs)
      allSamples = c(allSamples, sample(1:nrow(kenyaPop), numRural, TRUE, prob=probs))
      allUrban = c(allUrban, rep(FALSE, numRural))
    }
  }
  
  EAIs = allSamples
  kenyaEAs = kenyaPop[EAIs,]
  coordsX = runif(numEAs, min=kenyaEAs$lon - lonRes/2, max=kenyaEAs$lon + lonRes/2)
  coordsY = runif(numEAs, min=kenyaEAs$lat - latRes/2, max=kenyaEAs$lat + latRes/2)
  kenyaEAs$lon = coordsX
  kenyaEAs$lat = coordsY
  decreaseFac = totalKenyaPop/sum(kenyaEAs$pop)
  kenyaEAs$pop = kenyaEAs$pop*decreaseFac
  kenyaEAs$region = countyToRegion(kenyaEAs$admin1)
  kenyaEAs$urban = allUrban
  
  # project coordinates to northing/easting
  projCoords = projKenya(coordsX, coordsY)
  kenyaEAs$east = projCoords[,1]
  kenyaEAs$north = projCoords[,2]
  kenyaEAs$pixelIs = EAIs
  
  kenyaEAs
}

# get empirical mothers and children distribution per cluster
getSurveyEmpiricalDistributions = function(data=NULL, yearOfBirth = 2010) {
  # read in DHS data
  if(is.null(data))
    data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))
  
  # initialize lists of numbers for the empirical distributions
  clusters = unique(data$v001)
  nHouseholds = rep(0, length(clusters))
  nMothers = list()
  nChildren = list()
  
  for(i in 1:length(clusters)) {
    if(i %% 25 == 0)
      print(paste0("household ", i, "/", length(nHouseholds)))
    
    thisCluster = clusters[i]
    clusterData = data[data$v001 == thisCluster, ]
    households = unique(data$v002[data$v001 == thisCluster])
    nHouseholds[i] = max(households)
    
    for(j in 1:length(households)) {
      thisHousehold = households[j]
      thisData = clusterData[clusterData$v002 == thisHousehold, ]
      nMothers = c(nMothers, list(thisData$v138[1]))
      nChildren = c(nChildren, list(sum(thisData$b2 == yearOfBirth)))
    }
  }
  
  nMothers = unlist(nMothers)
  nChildren = unlist(nChildren)
  householdDistribution = ecdf(nHouseholds)
  motherDistribution = ecdf(nMothers)
  childrenDistribution = ecdf(nChildren)
  list(households=householdDistribution, mothers=motherDistribution, children=childrenDistribution)
}

# get empirical mothers and children distribution per cluster
# maxAge: the maximum age of the children being included in the distribution rounded down
# doNeonatal: if TRUE and maxAge == 0, randomly remove 11/12 of births
getSurveyEmpiricalDistributions2 = function(data=NULL, dataDHS=NULL, maxAge=0, doNeonatal=FALSE) {
  
  # check if inputs make sense
  if(doNeonatal && maxAge != 0)
    stop("If doNeonatal is set to TRUE, maxAge must be 0")
  
  # read in DHS data for the estimate of the number of households within each cluster
  print("reading in DHS data")
  if(is.null(dataDHS))
    dataDHS <- data.frame(read_dta("../U5MR/Kenya2014BirthRecode/KEBR70FL.DTA"))

  # initialize lists of numbers for the empirical distributions
  print("calculating households per cluster empirical distribution")
  urbanDat = dataDHS[dataDHS$v025 == 1, ]
  ruralDat = dataDHS[dataDHS$v025 == 2, ]
  nHouseholds = aggregate(dataDHS$v002, data.frame(v001=dataDHS$v001), function(x) {max(x)})$x
  nHouseholdsUrban = aggregate(urbanDat$v002, data.frame(v001=urbanDat$v001), function(x) {max(x)})$x
  nHouseholdsRural = aggregate(ruralDat$v002, data.frame(v001=ruralDat$v001), function(x) {max(x)})$x
  
  nHouseholds=nHouseholds[(nHouseholds >= 50) & (nHouseholds <= 150)]
  nHouseholdsUrban=nHouseholdsUrban[(nHouseholdsUrban >= 50) & (nHouseholdsUrban <= 150)]
  nHouseholdsRural=nHouseholdsRural[(nHouseholdsRural >= 50) & (nHouseholdsRural <= 150)]

  # read in census data for the other empirical distributions
  print("reading and census data if necessary")
  if(is.null(data))
    data <- data.frame(read_dta("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/popSurvey2009/Population_2009KPHC_10PCT_STATA.dta"))
  data = data[, c("PROVINCE", "DISTRICT", "DIVISION", "LOCATION", "SUBLOC", "COUNTY", 
                  "RecreatedEANO", "HHNO", "EATYPE", "BARCODE", "TIF", "P12", "P13")]
  # data2 = data[, c(1:16, 19, 20)]
  
  print("beginning data table operations")
  dat = as.data.table(data)
  # dat2 = as.data.table(data2)
  print(system.time(outMothers <- dat[, {
    .("nMothers"=uniqueN(P13[P13 != 0 & P12 <= maxAge]), "urban"=any(EATYPE != 1))
  }
  , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  # print(system.time(outMothers2 <- dat2[, {
  #   .("nMothers"=uniqueN(P13[P13 != 0 & P12 <= maxAge]), "urban"=any(EATYPE != 1))
  # }
  # , by = .(REC_TYPE, PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, 
  #          BARCODE, HHTYPE, EASTATUS, EATYPE, CONST, TIF)]))
  
  dat0 = dat[P13 != 0 & P12 <= maxAge,]
  # dat02 = dat2[P13 != 0 & P12 <= maxAge,]
  print(system.time(outChildren <- dat0[, {
    if(any(unlist(tapply(P13, P13, length, simplify = FALSE))) >= 15)
      print(.SD)
    .("nChildren"=unlist(tapply(P13, P13, length, simplify = FALSE)), "urban"=rep(any(EATYPE != 1), uniqueN(P13)))
  }
  , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  # print(system.time(outChildren2 <- dat02[, {
  #   if(any(unlist(tapply(P13, P13, length, simplify = FALSE))) >= 15)
  #     print(.SD)
  #   .("nChildren"=unlist(tapply(P13, P13, length, simplify = FALSE)), "urban"=rep(any(EATYPE != 1), uniqueN(P13)))
  # }
  # , by = .(REC_TYPE, PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, 
  #          BARCODE, HHTYPE, EASTATUS, EATYPE, CONST, TIF)]))
  
  nMothers = outMothers[["nMothers"]]
  nMothersUrban = outMothers[urban == TRUE, nMothers]
  nMothersRural = outMothers[urban == FALSE, nMothers]
  nChildren = outChildren[["nChildren"]]
  if(doNeonatal) {
    nChildren = rbinom(length(nChildren), size=nChildren, prob=1/12)
  }
  nChildrenUrban = outChildren[urban == TRUE, nChildren]
  nChildrenRural = outChildren[urban == FALSE, nChildren]
  
  householdDistribution = ecdf(nHouseholds)
  motherDistribution = ecdf(nMothers)
  childrenDistribution = ecdf(nChildren)
  householdDistributionUrban = ecdf(nHouseholdsUrban)
  motherDistributionUrban = ecdf(nMothersUrban)
  childrenDistributionUrban = ecdf(nChildrenUrban)
  householdDistributionRural = ecdf(nHouseholdsRural)
  motherDistributionRural = ecdf(nMothersRural)
  childrenDistributionRural = ecdf(nChildrenRural)
  list(households=householdDistribution, mothers=motherDistribution, children=childrenDistribution,
       householdsUrban=householdDistributionUrban, mothersUrban=motherDistributionUrban, childrenUrban=childrenDistributionUrban,
       householdsRural=householdDistributionRural, mothersRural=motherDistributionRural, childrenRural=childrenDistributionRural)
}

# get empirical women distribution per household
# minAge: the minimum age of the women being included in the distribution rounded down
# maxAge: the maximum age of the women being included in the distribution rounded down
getSurveyEmpiricalDistributionsWomen = function(data=NULL, dataDHS=NULL, minAge=20, maxAge=29, 
                                                datawd="~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/") {
  require(haven)
  thiswd = getwd()
  setwd(datawd)
  
  # read in census data for the other empirical distributions
  print("reading in census data if necessary")
  if(is.null(data))
    data <- data.frame(read_dta("popSurvey2009/Population_2009KPHC_10PCT_STATA.dta"))
  # data = data[, c("PROVINCE", "DISTRICT", "DIVISION", "LOCATION", "SUBLOC", "COUNTY", 
  #                 "RecreatedEANO", "HHNO", "EATYPE", "BARCODE", "P12", "P13")]
  data = data[, c("PROVINCE", "DISTRICT", "DIVISION", "LOCATION", "SUBLOC", "COUNTY", 
                  "RecreatedEANO", "HHNO", "EATYPE", "BARCODE", "P11", "P12")]
  
  # P11: sex (1: male, 2: female)
  # P12: age in years completed
  
  print("beginning data table operations")
  dat = as.data.table(data)
  # print(system.time(outMothers <- dat[, {
  #   .("nMothers"=uniqueN(P13[P13 != 0 & P12 <= maxAge]), "urban"=any(EATYPE != 1))
  # }
  # , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  print(system.time(outWomen <- dat[, {
    .("nWomen"=length(P11[P11 != 1 & P12 <= maxAge & P12 >= minAge]), "urban"=any(EATYPE != 1))
  }
  , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  
  nWomen = outWomen[["nWomen"]]
  nWomenUrban = outWomen[urban == TRUE, nWomen]
  nWomenRural = outWomen[urban == FALSE, nWomen]
  
  womenDistribution = ecdf(nWomen)
  womenDistributionUrban = ecdf(nWomenUrban)
  womenDistributionRural = ecdf(nWomenRural)
  
  setwd(thiswd)
  list(women=womenDistribution, womenUrban=womenDistributionUrban, womenRural=womenDistributionRural)
}

# set thresholds within each county based on percent population urban
setThresholds = function() {
  require(raster)
  
  load("../U5MR/kenyaPopProj.RData")
  load("../U5MR/adminMapData.RData")
  load("../U5MR/poppc.RData")
  
  getCountyThresh = function(countyName) {
    # if Nairobi or Mombasa, always urban
    if((countyName == "Nairobi") || (countyName == "Mombasa"))
      return(-Inf)
    
    # do the setup
    thisCounty = as.character(kenyaPop$admin1) == countyName
    thisPop = kenyaPop$popOrig[thisCounty]
    thisTot = sum(thisPop)
    pctUrb = poppc$pctUrb[poppc$County == countyName]/100
    pctRural = 1 - pctUrb
    
    # objective function to minimize
    # objFun = function(thresh) {
    #   curPctUrb = sum(thisPop[thisPop > thresh])/thisTot
    #   (curPctUrb - pctUrb)^2
    # }
    
    # do optimization
    # out = optim(10, objFun)
    # thresh = out$par
    # out = optimize(objFun, c(.01, 50000))
    # thresh = out$par
    
    # calculate threshold by integrating ecdf via sorted value cumulative sum
    sortedPop = sort(thisPop)
    cumsumPop = cumsum(sortedPop)
    threshI = match(1, cumsumPop >= thisTot*pctRural)
    thresh = sortedPop[threshI]
    
    # print(paste0("pctUrb: ", pctUrb, "; resPctUrb: ", sum(thisPop[thisPop > thresh])/thisTot, "; thresh: ", thresh, "; obj: ", out$objective))
    thresh
  }
  
  # compute threshold for each county
  counties = poppc$County
  threshes = sapply(counties, getCountyThresh)
  
  list(counties=counties, threshes=threshes)
}

# simulate EAs household, mother, and child values in the clusters based on empirical distributions 
# from getSurveyEmpiricalDistributions. First simulate the "long form", where each row is 
# a household, then simulate the "short form", were each row is a EA
simEAEmpirical = function(empiricalDistributions, nEAs=423) {
  # simulate the long form
  householdValues = knots(empiricalDistributions$households)
  householdProbs = empiricalDistributions$households(householdValues) - empiricalDistributions$households(householdValues - 1)
  motherValues = knots(empiricalDistributions$mothers)
  motherProbs = empiricalDistributions$mothers(motherValues) - empiricalDistributions$mothers(motherValues - 1)
  childrenValues = knots(empiricalDistributions$children)
  childrenProbs = empiricalDistributions$children(childrenValues) - empiricalDistributions$children(childrenValues - 1)
  
  # simulate the number of households per numeration area, and mothers and children per household
  nHouseholds = sample(householdValues, nEAs, TRUE, householdProbs)
  nMothers = sample(motherValues, sum(nHouseholds), TRUE, motherProbs)
  nChildren = sample(childrenValues, sum(nHouseholds), TRUE, childrenProbs)
  
  # generate full data sets
  longform = data.frame(eaI=rep(1:nEAs, nHouseholds), mothers=nMothers, children=nChildren)
  totalMothers = aggregate(longform$mothers, data.frame(list(house=longform$eaI)), sum)$x
  totalChildren = aggregate(longform$children, data.frame(list(house=longform$eaI)), sum)$x
  shortform = data.frame(eaI=1:nEAs, mothers=totalMothers, children=totalChildren)
  
  list(longform=longform, shortform=shortform)
}

# same as simEAEmpirical, except use empirical distributions from getSurveyEmpiricalDistributions2, 
# which stratify by urban/rural
simEAEmpirical2 = function(empiricalDistributions, eaDat, eaDatLong) {
  urban = eaDat$urban
  nEAs = length(urban)
  nUrban = sum(urban)
  nRural = nEAs - nUrban
  
  # simulate the number of households per numeration area
  nHouseholdsUrban = recdf(nUrban, empiricalDistributions$householdsUrban)
  nHouseholdsRural = recdf(nRural, empiricalDistributions$householdsRural)
  nHouseholds = urban
  nHouseholds[urban] = nHouseholdsUrban
  nHouseholds[!urban] = nHouseholdsRural
  nMothersUrban = recdf(nUrban, empiricalDistributions$mothersUrban)
  nMothersRural = recdf(nRural, empiricalDistributions$mothersRural)
  nMothers = urban
  nMothers[urban] = nMothersUrban
  nMothers[!urban] = nMothersRural
  nChildrenUrban = recdf(nUrban, empiricalDistributions$childrenUrban)
  nChildrenRural = recdf(nRural, empiricalDistributions$childrenRural)
  nChildren = urban
  nChildren[urban] = nChildrenUrban
  nChildren[!urban] = nChildrenRural
  
  # generate full data sets
  longform = data.frame(eaI=rep(1:nEAs, nHouseholds), mothers=nMothers, children=nChildren)
  totalMothers = aggregate(longform$mothers, data.frame(list(house=longform$eaI)), sum)$x
  totalChildren = aggregate(longform$children, data.frame(list(house=longform$eaI)), sum)$x
  shortform = data.frame(eaI=1:nEAs, mothers=totalMothers, children=totalChildren)
  
  list(longform=longform, shortform=shortform)
}

saveSimDat = function(seed=123, HHoldVar=0) {
  # out = simDat(kenyaEAs, beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5, effRange=300, seed=seed)
  out = simDat(kenyaEAs, beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5, effRange=300, seed=seed, HHoldVar=HHoldVar)
  eaDat = out$eaDat
  clustDat = out$clustDat
  eaID = out$eaI
  save(eaDat, clustDat, eaID, file=paste0("simDatHH", round(HHoldVar, 3), ".RData"))
}

# given a vector of populations with associated counties, and a threshes 
# object such as an output to a setThresholds function, this computes 
# a vector of length equal to the pops vector whether the point is urban.
setUrbanByThreshes = function(pops, counties, threshes) {
  allThreshes = sapply(1:length(pops), function(i) {threshes$threshes[threshes$counties == counties[i]]})
  urban = pops > allThreshes
  urban
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample
simClusters = function(eaDat, numClusters = 423, seed=123) {
  set.seed(seed)
  
  # set number of clusters empirically if not otherwise set, scale by a factor depending on how many clusters we sample
  thisclustpc = clustpc
  thisclustpc[,2:ncol(thisclustpc)] = round(thisclustpc[,2:ncol(thisclustpc)] * (numClusters / sum(thisclustpc[,4])))
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  counties = as.character(thisclustpc$County)
  urban = eaDat$urban
  allSamples = c()
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = as.character(eaDat$admin1) == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # number of rural and urban clusters to sample
    numUrban = thisclustpc$clustUrb[thisclustpc$County == countyName]
    numRural = thisclustpc$clustRur[thisclustpc$County == countyName]
    
    # urban samples
    thisCountyI = countyI & urban
    # probs = thisCountyI * eaDat$popOrig
    probs = thisCountyI
    probs = probs/sum(probs)
    allSamples = c(allSamples, sample(1:nrow(eaDat), numUrban, TRUE, prob=probs))
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      # probs = thisCountyI * eaDat$popOrig
      probs = thisCountyI
      probs = probs/sum(probs)
      allSamples = c(allSamples, sample(1:nrow(eaDat), numRural, TRUE, prob=probs))
    }
  }
  
  clustDat = eaDat[allSamples,]
  list(sampleI = allSamples, clustDat = clustDat)
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample.  Simulate same number of clusters per county (9)
# urbanProps: either a single proportion or a sequence of proportions of length 47 in order corresponding to counties
simClusters2 = function(eaDat, numClusters=423, urbanProps=NULL, counties=as.character(clustpc$County), seed=NULL, doJitter=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  # compute number of clusters to sample from each county
  numPerCounty = round(numClusters/47)
  if(round(numClusters/47) != numClusters/47)
    warning("number of clusters doens't divide evenly into number of counties")
  
  # set the sampling frame
  if(is.null(urbanProps)) {
    # set proportion of urban clusters empirically if not otherwise set, scale by a factor depending on how many clusters we sample
    thisclustpc = clustpc
    thisclustpc[,2:ncol(thisclustpc)] = round(thisclustpc[,2:ncol(thisclustpc)] * (numClusters / sum(thisclustpc$clustTotal)))
  }
  else {
    # set sampling frame based on user proportions that are urban
    
    if(length(urbanProps) == 1)
      urbanProps = rep(urbanProps, 47)
    
    thisclustpc = data.frame(list(County = counties, clustUrb=numPerCounty*urbanProps, clustRur = numPerCounty*(1-urbanProps), 
                                  clustTotal = numPerCounty))
  }
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  urban = eaDat$urban
  allSamples = c()
  sampleWeights = c()
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = as.character(eaDat$admin1) == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # number of rural and urban clusters to sample (sample numPerCounty total per county, so must rescale each)
    numUrban = thisclustpc$clustUrb[thisclustpc$County == countyName] * (numPerCounty/thisclustpc$clustTotal[thisclustpc$County == countyName])
    numRural = thisclustpc$clustRur[thisclustpc$County == countyName] * (numPerCounty/thisclustpc$clustTotal[thisclustpc$County == countyName])
    
    # make sure there's a whole number of rural and urban summing to numPerCounty (9 by default)
    if(round(numUrban) + round(numRural) != numPerCounty) {
      if(runif(1) > .5) {
        numUrban = ceiling(numUrban)
        numRural = floor(numRural)
      }
      else {
        numUrban = floor(numUrban)
        numRural = ceiling(numRural)
      }
    }
    else {
      numRural = round(numRural)
      numUrban = round(numUrban)
    }
    
    # urban samples
    thisCountyI = countyI & urban
    # probs = thisCountyI * eaDat$popOrig
    probs = thisCountyI
    probs = probs/sum(probs)
    newSamples = sample(1:nrow(eaDat), numUrban, FALSE, prob=probs)
    newWeights = rep(1/(numUrban/sum(eaDat$urban[thisCountyI])), numUrban)
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      # probs = thisCountyI * eaDat$popOrig
      probs = thisCountyI
      probs = probs/sum(probs)
      newSamples = c(newSamples, sample(1:nrow(eaDat), numRural, FALSE, prob=probs))
      newWeights = c(newWeights, rep(1/(numRural/sum(!eaDat$urban[thisCountyI])), numRural))
    }
    
    # update clusters and sampling weights for this county
    allSamples = c(allSamples, newSamples)
    sampleWeights = c(sampleWeights, newWeights)
  }
  
  clustDat = eaDat[allSamples,]
  clustDat$samplingWeight = sampleWeights
  
  # computer cluster locations if necessary
  if(doJitter) {
    stop("how do we want to jitter the data?")
  }
  
  list(sampleI = allSamples, clustDat = clustDat)
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample.  Simulate same number of clusters per county (9)
# urbanOverSample: oversample in urban areas by this ratio.  Any urban EA is urbanProp times more likely 
#            to be sampled than a rural EA. If urbanOverSample is 1 than no oversampling occurs
simClusters3 = function(eaDat, numClusters=423, urbanOverSample=1, nsim=1, seed=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  # compute number of clusters to sample from each county
  numPerCounty = round(numClusters/47)
  if(round(numClusters/47) != numClusters/47)
    warning("number of clusters doens't divide evenly into number of counties")
  
  # set some basic variables
  thisclustpc = clustpc
  counties = clustpc$County
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  urban = eaDat$urban
  eaIs = matrix(nrow=numClusters, ncol=nsim)
  sampleWeights = matrix(nrow=numClusters, ncol=nsim)
  curRow = 1
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = as.character(eaDat$admin1) == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # sampling probabilities depend on urban and rural strata
    thisUrban = countyI & urban
    thisRural = countyI & !urban
    probs = thisRural + thisUrban * urbanOverSample
    probs = probs/sum(probs)
    endRow = curRow - 1 + numPerCounty
    getSamples = function(i) {
      sample(1:nrow(eaDat), numPerCounty, FALSE, prob=probs)
    }
    thisEAIs = sapply(1:nsim, getSamples)
    
    eaIs[curRow:endRow, ] = thisEAIs
    sampleWeights[curRow:endRow, ] <- matrix(1/probs[thisEAIs], nrow=numPerCounty, ncol=nsim)
    curRow = endRow + 1
  }
  
  list(eaIs=eaIs, sampleWeights=sampleWeights)
}

# convert from output of simClusters3 to the result of generateSimDataSets.R
genAndreaFormatFromEAIs = function(eaDat, eaIs, sampleWeights) {
  lapply(1:ncol(eaIs), genClustDatFromEAIs, eaDat=eaDat, eaIs=eaIs, sampleWeights=sampleWeights)
}

# given the above function (simClusters2Mult), generates a dataframe following the format 
# of simClusters2 given a simulation index, i
genClustDatFromEAIs = function(eaDat, eaIs, sampleWeights, i) {
  clustDat = eaDat[eaIs[,i],]
  clustDat$samplingWeight = sampleWeights[,i]
  clustDat$eaIs = eaIs[,i]
  clustDat
}

# convert from output of simClusters3 to the result of generateSimDataSets.R
genAndreaFormatFromEAIsLong = function(eaDat, eaIs, eaDatLong, HHIs, sampleWeights) {
  lapply(1:ncol(eaIs), genClustDatFromEAIsLong, eaDat=eaDat, eaIs=eaIs, eaDatLong=eaDatLong, 
         HHIs=HHIs, sampleWeights=sampleWeights)
}

# convert from output of simClusters3 to the result of generateSimDataSets.R
genAndreaFormatFromEAIsLong2 = function(eaDat, eaIs, eaDatLong, HHIs, sampleWeights, 
                                        doSmoothRisk=TRUE, doFineScaleRisk=TRUE, doGriddedRisk=TRUE) {
  lapply(1:ncol(eaIs), genClustDatFromEAIsLong2, eaDat=eaDat, eaIs=eaIs, eaDatLong=eaDatLong, 
         HHIs=HHIs, sampleWeights=sampleWeights, 
         doSmoothRisk=doSmoothRisk, doFineScaleRisk=doFineScaleRisk, 
         doGriddedRisk=doGriddedRisk)
}

# given the above function (simClusters2Mult), generates a dataframe following the format 
# of simClusters2 given a simulation index, i
genClustDatFromEAIsLong = function(eaDat, eaIs, eaDatLong, HHIs, sampleWeights, i) {
  clustDatLong = as.data.table(eaDatLong[HHIs[,i],])
  
  # aggregate the long format to the short format
  # TODO: figure out what columns are in the cluster data and how to aggregate
  # clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], pop=pop[1], popOrig=popOrig[1], 
  #                             admin1=admin1[1], region=region[1], 
  #                             east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
  #                             numWomen=sum(numWomen), numChildren=sum(numChildren), 
  #                             died=sum(died), 
  #                             probsHH=probsHH[1], trueProbDeathNoNug=trueProbDeathNoNug[1], 
  #                             trueProbDeath=trueProbDeath[1]), 
  #                         by=eaIs]
  clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                              popOverall=popOverall[1], popTarget=popTarget[1], 
                              admin1=admin1[1], admin2=admin2[1], region=region[1], 
                              east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                              n=sum(n), y=sum(y), 
                              plcpb=plcpb[1], pLcpb=pLcpb[1], pLCpb=pLCpb[1], pLCPb=pLCPb[1], pLCPB=pLCPB[1]), 
                          by=eaIs]
  
  clustDat$samplingWeight = sampleWeights[,i]
  clustDat$eaIs = eaIs[,i]
  clustDat
}

# given the above function (simClusters2Mult), generates a dataframe following the format 
# of simClusters2 given a simulation index, i
genClustDatFromEAIsLong2 = function(eaDat, eaIs, eaDatLong, HHIs, sampleWeights, i, 
                                    doSmoothRisk=TRUE, doFineScaleRisk=TRUE, doGriddedRisk=TRUE) {
  clustDatLong = as.data.table(eaDatLong[HHIs[,i],])
  
  # aggregate the long format to the short format
  # TODO: figure out what columns are in the cluster data and how to aggregate
  # clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], pop=pop[1], popOrig=popOrig[1], 
  #                             admin1=admin1[1], region=region[1], 
  #                             east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
  #                             numWomen=sum(numWomen), numChildren=sum(numChildren), 
  #                             died=sum(died), 
  #                             probsHH=probsHH[1], trueProbDeathNoNug=trueProbDeathNoNug[1], 
  #                             trueProbDeath=trueProbDeath[1]), 
  #                         by=eaIs]
  if(doSmoothRisk && doFineScaleRisk && doGriddedRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pSmoothRisk=pSmoothRisk[1], pFineScaleRisk=pFineScaleRisk[1], 
                                pFineScalePrevalence=sum(y)/sum(n), 
                                pGriddedRisk=pGriddedRisk[1]), 
                            by=eaIs]
  } else if(doSmoothRisk && doFineScaleRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pSmoothRisk=pSmoothRisk[1], pFineScaleRisk=pFineScaleRisk[1], 
                                pFineScalePrevalence=sum(y)/sum(n)), 
                            by=eaIs]
  } else if(doSmoothRisk && doGriddedRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pSmoothRisk=pSmoothRisk[1], 
                                pFineScalePrevalence=sum(y)/sum(n), 
                                pGriddedRisk=pGriddedRisk[1]), 
                            by=eaIs]
  } else if(doFineScaleRisk && doGriddedRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pFineScaleRisk=pFineScaleRisk[1], 
                                pFineScalePrevalence=sum(y)/sum(n), 
                                pGriddedRisk=pGriddedRisk[1]), 
                            by=eaIs]
  } else if(doSmoothRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pSmoothRisk=pSmoothRisk[1], 
                                pFineScalePrevalence=sum(y)/sum(n)), 
                            by=eaIs]
  } else if(doFineScaleRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pFineScaleRisk=pFineScaleRisk[1], 
                                pFineScalePrevalence=sum(y)/sum(n)), 
                            by=eaIs]
  } else if(doGriddedRisk) {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pFineScalePrevalence=sum(y)/sum(n), 
                                pGriddedRisk=pGriddedRisk[1]), 
                            by=eaIs]
  } else {
    clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], 
                                popDensity=popDensity[1], popDensityTarget=popDensityTarget[1], 
                                area=area[1], subarea=subarea[1], 
                                east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                                n=sum(n), y=sum(y), 
                                pFineScalePrevalence=sum(y)/sum(n)), 
                            by=eaIs]
  }
  
  clustDat$samplingWeight = sampleWeights[,i]
  clustDat$eaIs = eaIs[,i]
  clustDat$pFineScalePrevalence[clustDat$n == 0] = NA
  clustDat
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample.  Simulate same number of clusters per county (9)
# eaDatLong: long format enumeration area table, where each row is a household
# nPerStrata: if fixedPerStrata == TRUE, this is the number of clusters sampled per strata
# SRS: either SRS or PPS sampling
# pctUrb: the percentage for each county's general population that is urban
simClustersEmpirical = function(eaDat, eaDatLong, nsim=1, seed=NULL, urbanOverSamplefrac=0, 
                                nHHSampled=25, fixedPerStrata=FALSE, nPerStrata=3, representativeSampling=FALSE, verbose=TRUE, 
                                thisclustpc=clustpc, pctUrb=poppc$pctUrb) {
  if(!is.null(seed))
    set.seed(seed)
  
  # get the number of clusters for each strata
  nUrbanClusters = thisclustpc$clustUrb
  nRuralClusters = thisclustpc$clustRur
  
  # make sure the ea tables are a data.tables
  if(!("data.table" %in% class(eaDatLong)))
    eaDatLong = as.data.table(eaDatLong)
  if(!("data.table" %in% class(eaDat)))
    eaDat = as.data.table(eaDat)
  eaDat = cbind(eaIs=1:nrow(eaDat), eaDat)
  # admin1 = eaDat$admin1
  admin1 = eaDat$area
  
  if(!fixedPerStrata) {
    # set number of clusters empirically if not otherwise set, scale by a factor depending on how much to oversample urban areas
    # thisclustpc = clustpc
    # thisclustpc$clustUrb = round(thisclustpc$clustUrb * (1 + urbanOverSamplefrac))
    # thisclustpc$clustUrb[thisclustpc$clustUrb >= thisclustpc$clustTotal] = thisclustpc$clustTotal[thisclustpc$clustUrb >= thisclustpc$clustTotal]
    # thisclustpc$clustRur = thisclustpc$clustTotal - thisclustpc$clustUrb
    tempclustpc = thisclustpc
    tempclustpc$clustRur = round(tempclustpc$clustRur)
    tempclustpc$clustRur[tempclustpc$clustRur < 0] = 0
    tempclustpc$clustUrb = tempclustpc$clustTotal - tempclustpc$clustRur
    
    if(representativeSampling) {
      # in this case, we want the naive model to be unbiased so we need the proportion of urban clusters and rural clusters to 
      # be representative of the population within each county
      tempclustpc$clustUrb = (pctUrb / 100) * tempclustpc$clustTotal
      tempclustpc$clustRur = (1 - pctUrb / 100) * tempclustpc$clustTotal
      
      # NOTE: the number of urban and rural clusters will be fractions, so we need to pick the number randomly for each county 
      #       so that on average the estimate within each county will be unbiased
    }
  } else {
    tempclustpc = thisclustpc
    tempclustpc$clustUrb = nPerStrata
    tempclustpc$clustRur = nPerStrata
    tempclustpc$clustRur[tempclustpc$County == "Mombasa"] = 0
    tempclustpc$clustRur[tempclustpc$County == "Nairobi"] = 0
  }
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  counties = as.character(tempclustpc$County)
  urban = eaDat[["urban"]]
  
  # check if we need to randomly allocate clusters
  roundedUrban = round(tempclustpc$clustUrb, 0)
  roundedRural = round(tempclustpc$clustRur, 0)
  randomClusters = any(roundedUrban != tempclustpc$clustUrb) | any(roundedRural != tempclustpc$clustRur)
  
  # if we have random clusters, precalculate certain probabilities
  if(randomClusters) {
    extraUrban = tempclustpc$clustUrb - floor(tempclustpc$clustUrb)
    extraRural = tempclustpc$clustRur - floor(tempclustpc$clustRur)
    probUrbanTotal = sum(extraUrban) / sum(extraUrban + extraRural)
    probUrban = extraUrban / sum(extraUrban)
    probRural = extraRural / sum(extraRural)
  }
  
  generateSample = function() {
    tempclustpcI = tempclustpc
    
    if(randomClusters) {
      # we need to randomly allot clusters to strata so the naive case remains unbiased
      
      # calculate the number of clusters we need to randomly allot
      nExtra = round(sum(extraUrban + extraRural), 0)
      
      # calculate the number of extra urban and rural clusters
      nExtraUrbanTotal = rbinom(1, nExtra, probUrbanTotal)
      nExtraRuralTotal = nExtra - nExtraUrbanTotal
      
      # calculate the number of extra clusters in each strata conditional on the number of urban and rural
      nExtraUrban = rmultinom(1, nExtraUrbanTotal, probUrban)
      nExtraRural = rmultinom(1, nExtraRuralTotal, probRural)
      
      # add extra clusters to their corresponding strata
      tempclustpcI$clustUrb = floor(tempclustpcI$clustUrb) + nExtraUrban
      tempclustpcI$clustRur = floor(tempclustpcI$clustRur) + nExtraRural
      tempclustpcI$clustTotal = tempclustpcI$clustUrb + tempclustpcI$clustRur
    }
    
    allSamples = c()
    # allSamples = list()
    for(i in 1:nrow(tempclustpcI)) {
      countyName = counties[i]
      
      # if(as.character(countyName) == "Isiolo") {
      #   # Lamu
      #   # Elgeyo Marakwet
      #   # Marsabit
      #   # Taita Taveta
      #   # Isiolo
      #   browser()
      # }
      
      countyI = eaDat[,as.character(admin1) == countyName]
      if(verbose)
        print(paste0("County ", i, "/", length(counties), ": ", countyName))
      
      # number of rural and urban clusters to sample
      numUrban = tempclustpcI$clustUrb[tempclustpcI$County == countyName]
      numRural = tempclustpcI$clustRur[tempclustpcI$County == countyName]
      
      # urban samples
      thisCountyI = countyI & urban
      
      # probs = thisCountyI * eaDat$popOrig
      if(representativeSampling)
        probs = thisCountyI
      else
        probs = thisCountyI * eaDat$nHH
      probs = probs/sum(probs)
      # allSamples = c(allSamples, sample(1:nrow(eaDat), numUrban, FALSE, prob=probs))
      zs = probs[thisCountyI] * numUrban # the inclusion probabilities
      if(any(zs >= 1-1e-6) && sum(thisCountyI) != numUrban) {
        warning("inclusion probabilities greater than one, adjusting... Sampling will not quite be PPS")
        while(any(zs >= 1-1e-6)) {
          zs[zs >= 1-1e-6] = max(zs[zs < 1-1e-6])
          zs = zs * (numUrban / sum(zs))
        }
      }
      
      # first check if we are trying to sample all of the stratum EAs. Otherwise 
      # use the Midzuno method for PPS sampling
      if(sum(thisCountyI) == numUrban) {
        newSamples = which(thisCountyI)
      } else {
        # newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPtille(zs))]
        newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPmidzuno(zs))]
      }
      
      allSamples = c(allSamples, newSamples)
      
      # allSamples = c(allSamples, pps.sampling(eaDat$nHH[thisCountyI], numUrban, (1:nrow(eaDat))[thisCountyI], FALSE, prob=probs))
      
      # rural samples (if the county has any)
      if((countyName != "Mombasa") && (countyName != "Nairobi")) {
        thisCountyI = countyI & !urban
        # probs = thisCountyI * eaDat$popOrig
        if(representativeSampling)
          probs = thisCountyI
        else
          probs = thisCountyI * eaDat$nHH
        probs = probs/sum(probs)
        
        zs = probs[thisCountyI] * numRural # the inclusion probabilities
        if(any(zs >= 1-1e-6) && sum(thisCountyI) != numRural) {
          warning("inclusion probabilities greater than one, adjusting... Sampling will not quite be PPS")
          while(any(zs >= 1-1e-6)) {
            zs[zs >= 1-1e-6] = max(zs[zs < 1-1e-6])
            zs = zs * (numRural / sum(zs))
          }
        }
        
        # first check if we are trying to sample all of the stratum EAs. Otherwise 
        # use the Midzuno method for PPS sampling
        if(sum(thisCountyI) == numRural) {
          newSamples = which(thisCountyI)
        } else {
          # newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPtille(zs))]
          newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPmidzuno(zs))]
        }
        
        
        # allSamples = c(allSamples, sample(1:nrow(eaDat), numRural, FALSE, prob=probs))
        allSamples = c(allSamples, newSamples)
      }
    }
    sort(allSamples)
  }
  # the comment out line below works, accept it does not enable you to print out progress
  # eaIs = replicate(nsim, generateSample())
  eaIs = matrix(ncol=nsim, nrow=sum(tempclustpc$clustUrb) + sum(tempclustpc$clustRur))
  for(i in 1:nsim) {
    if(verbose && i %% 10 == 0) {
      print(paste0("simulating data set ", i, "/", nsim))
    }
    eaIs[, i] = generateSample()
  }
  
  ## decide which households within each cluster to sample
  eaDatLong = data.table(HHI=1:nrow(eaDatLong), eaDatLong)
  sampleHouseholds = function(sim) {
    assign("thisEAIs", eaIs[,sim], globalenv())
    clusterDatLong = eaDatLong[eaIs %in% get("thisEAIs", globalenv()), ]
    assign("nHHSampled", nHHSampled, globalenv())
    HHIs = clusterDatLong[, {
      .("HHIs"=sort(sample(HHI, get("nHHSampled", globalenv()), FALSE)))
    }, by=eaIs]
    # cumulativeHouseholds = c(0, cumsum(eaDat[["nHH"]]))
    # HHIs[["HHIs"]] = HHIs[["HHIs"]] + cumulativeHouseholds[HHIs[["eaIs"]] - 1]
    HHIs[["HHIs"]]
  }
  HHIs = sapply(1:nsim, sampleHouseholds)
  
  # compute the probability of each enumeration area and household being in our sample
  probsEA = rep(0, nrow(eaDat))
  for(i in 1:nrow(tempclustpc)) {
    countyName = counties[i]
    countyI = eaDat[,as.character(admin1) == countyName]
    
    # number of rural and urban clusters to sample
    numUrban = tempclustpc$clustUrb[tempclustpc$County == countyName]
    numRural = tempclustpc$clustRur[tempclustpc$County == countyName]
    
    # urban samples
    thisCountyI = countyI & urban
    if(representativeSampling)
      probs = thisCountyI
    else
      probs = thisCountyI * eaDat[["nHH"]]
    probs = probs * (1/sum(probs))
    zs = probs[thisCountyI] * numUrban # the inclusion probabilities
    if(any(zs >= 1-1e-6)) {
      warning("inclusion probabilities greater than one, modifying to be .999")
      zs[zs >= 1-1e-6] = .999
    }
    probsEA[thisCountyI] = zs
    
    if(representativeSampling) {
      if(any(abs(unique(zs) - tempclustpc$clustUrb[i]/sum(eaDat$admin1 == countyName & eaDat$urban == TRUE)) > .0001)) {
        print("woah there")
      }
    }
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      # probs = thisCountyI * eaDat$popOrig
      if(representativeSampling)
        probs = thisCountyI
      else
        probs = thisCountyI * eaDat[["nHH"]]
      probs = probs * (1/sum(probs))
      zs = probs[thisCountyI] * numRural # the inclusion probabilities
      if(any(zs >= 1-1e-6)) {
        warning("inclusion probabilities greater than one, modifying to be .999")
        zs[zs >= 1-1e-6] = .999
      }
      probsEA[thisCountyI] = zs
      
      if(representativeSampling) {
        if(any(abs(unique(zs) - tempclustpc$clustRur[i]/sum(eaDat$admin1 == countyName & eaDat$urban == FALSE)) > .0001)) {
          print("woah there")
        }
      }
    }
  }
  probsHH = probsEA * (nHHSampled / eaDat[["nHH"]])
  
  # calculate the sampling weights. Each child has probability of being sampled same as its household, 
  # calculate the child's weight then multiply by number in HH or EA
  sampleWeightsLong = eaDatLong[["n"]] / probsHH[eaDatLong[["eaIs"]]]
  # sampleWeights = eaDat[["numChildren"]] / probsHH
  # sampleWeightsEA = sampleWeights
  # sampleWeights = matrix(sampleWeights[eaIs], ncol=nsim)
  getSampleWeightsCluster = function(i) {
    sampleWeightsClusterLong = sampleWeightsLong[HHIs[,i]]
    thisEAIs = eaDatLong[HHIs[,i], eaIs]
    aggregate(sampleWeightsClusterLong, by=list(thisEAIs), FUN="sum")$x
  }
  sampleWeights = sapply(1:nsim, getSampleWeightsCluster)
  list(eaIs=eaIs, HHIs=HHIs, sampleWeights=sampleWeights)
}

simClustersSRS = function(nsim=1, eaDat, eaDatLong, nEASampled, nHHSampled=25, seed=NULL) {
  
  # first sample the EAs
  eaIs = replicate(nsim, sample(1:nrow(eaDat), nEASampled, replace=FALSE))
  
  # now sample households
  getOneHHIs = function(i) {
    thisEAIs = eaIs[,i]
    
    thisHHIs = (1:nrow(eaDatLong))[eaDatLong$eaIs %in% thisEAIs]
    thisEADatLong = eaDatLong[thisHHIs,]
    temp = aggregate(thisHHIs, by=list(eaI=thisEADatLong$eaI), 
                     FUN=function(hhIs) {sample(hhIs, nHHSampled, replace=FALSE)})
    unlist(temp[[2]])
  }
  HHIs = sapply(1:nsim, getOneHHIs)
  
  # calculate sampling weights
  eaProb = nEASampled/nrow(eaDat)
  hhProb = nHHSampled/eaDat$nHH
  sampleProbs = matrix(eaProb * hhProb[c(eaIs)], ncol=nsim)
  sampleWeights = 1/sampleProbs
  
  list(eaIs=eaIs, HHIs=HHIs, sampleWeights=sampleWeights)
}

expit = function(x) { exp(x)/(1+exp(x)) }

# function for simulating data given enumeration areas and their info.  Simulate using 
# LatticeKrig model on logit scale
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
simDat = function(eaDat, clustDat=NULL, urbanProps=NULL, counties=as.character(clustpc$County), 
                  nsim=1, margVar=4, nu=1.5, NC=5, effRange=300, nLayer=3, beta0=0, gamma=-1, 
                  tausq=1, normalize=TRUE, nBuffer=5, womenFrac=1/3, 
                  numClusters=423, seed=NULL, fixedWomenPerClust=TRUE, HHoldVar=0) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(nsim != 1)
    stop("multiple simulations not yet implemented")
  
  ### first generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  LKArgs = list(coords=eaCoords, nsim=nsim, NC=NC, margVar=margVar, effRange=effRange, nu=nu, 
                     nLayer=nLayer, normalize=normalize, nBuffer=nBuffer)
  simVals = do.call("LKSimulator2", LKArgs)
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=nsim)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  if(is.null(clustDat))
    clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  
  # collect indices corresponding to EA and dataframe
  clustI = clustDat$sampleI
  clustDat = clustDat$clustDat
  
  # if fewer than 25 women in the EA, sample exactly all of them (only matters if fixedWomenPerClust == FALSE)
  numWomenClust = apply(cbind(numWomen, 25), 1, min)
  numChildrenClust = numWomenClust
  
  # helper function: for any EA index, get the number of mortalities sampled
  getClustMort = function(ind) {
    # convert from binomial data to multiple Bernouli trials
    thisDied = c(rep(1, eaDied[ind]), rep(0, numChildren[ind] - eaDied[ind]))
    
    # sum of sampled mortalities
    sum(sample(thisDied, numChildrenClust[ind], FALSE))
  }
  
  # get number of mortalities within each cluster out of 25 sampled households in cluster
  clustDied = sapply(clustI, getClustMort)
  
  # return simulated data
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  clustDat$died = clustDied
  clustDat$numWomen = numWomenClust[clustI]
  clustDat$numChildren = numChildrenClust[clustI]
  clustDat$trueProbDeath = probs[clustI]
  clustDat$trueProbDeathNoNug = probsNoNug[clustI]
  
  list(eaDat=eaDat, clustDat=clustDat, eaI=clustI)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# SPDE or LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
# useSPDE: simulate the datasets with SPDE model or LatticeKrig model
simDat2 = function(eaDat, clustDat=NULL, nsim=1, margVar=4, nu=1, NC=5, effRange=300, 
                   nLayer=3, beta0=0, gamma=-1, tausq=1, normalize=TRUE, nBuffer=5, 
                   womenFrac=1/3, urbanOverSample=1, numClusters=423, seed=NULL, 
                   HHoldVar=0, fixedWomenPerClust=TRUE, useSPDE=TRUE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### first generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  print("Simulating nationwide mortality rates and data")
  if(!useSPDE) {
  LKArgs = list(coords=eaCoords, nsim=1, NC=NC, margVar=margVar, effRange=effRange, nu=nu, 
                nLayer=nLayer, normalize=normalize, nBuffer=nBuffer)
  simVals = do.call("LKSimulator2", LKArgs)
  }
  else {
    if(nu != 1)
      stop("SPDE model only supports nu=1")
    
    SPDEArgs = list(coords=eaCoords, nsim=1, margVar=margVar, effRange=effRange)
    simVals = do.call("simSPDE", SPDEArgs)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
  }
  
  # return simulated data
  print("finishing up...")
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIs(eaDat, clustDat$eaIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
simDatLK = function(eaDat, clustDat=NULL, nsim=1, margVar=4, nu=1.5, NC=5, effRange=300, 
                      nLayer=3, beta0=0, gamma=-1, tausq=1, normalize=TRUE, nBuffer=5, 
                      womenFrac=1/3, urbanOverSample=1, numClusters=423, seed=NULL, fullEADat=NULL, 
                      HHoldVar=0, fixedWomenPerClust=TRUE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### first generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  print("Simulating nationwide mortality rates and data")
  LKArgs = list(coords=eaCoords, nsim=1, NC=NC, margVar=margVar, effRange=effRange, nu=nu, 
                nLayer=nLayer, normalize=normalize, nBuffer=nBuffer)
  simVals = do.call("LKSimulator2", LKArgs)
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
  }
  
  # return simulated data
  print("finishing up...")
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIs(eaDat, clustDat$eaIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# SPDE or LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# empiricalDistributions: the set of empirical distributions for households per cluster, 
#                         mothers per household, and children per mother stratified by urban
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDatShort: enumeration area data set in "short form" were each row is a cluster. 
#             "long form" is where each row is a household.
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the spatial process, excluding household end cluster effects. 
#          If 0, no spatial component is included
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
simDatEmpirical = function(empiricalDistributions, eaDat, clustDat=NULL, nsim=1, 
                           margVar=4, effRange=150, beta0=0, gamma=-1, tausq=1, 
                           urbanOverSamplefrac=0, seed=NULL, 
                           fullEADat=NULL, HHoldVar=0, nHHSampled=25) {
  if(!is.null(seed))
    set.seed(seed)
  
  # make edfun versions of the ecdfs for speed considerations
  fastDistributions = list()
  for(i in 1:length(empiricalDistributions)) {
    fastDistributions = c(fastDistributions, list(ecdf2edfun(empiricalDistributions[[i]])))
  }
  names(fastDistributions) = names(empiricalDistributions)
  
  ### simulate binomial data for each enumeration area
  # first generate the number of households
  numHouseholds = rep(0, nrow(eaDat))
  numHouseholds[eaDat$urban] = recdf(sum(eaDat$urban), fastDistributions$householdsUrban)
  numHouseholds[!eaDat$urban] = recdf(sum(!eaDat$urban), fastDistributions$householdsRural)
  eaDat$nHH = numHouseholds
  
  # now expand the eaDat table to be in the long format, were each row is a house
  rowsLong = rep(1:nrow(eaDat), numHouseholds)
  eaDatLong = eaDat[rowsLong, ]
  eaDatLong$eaIs = rowsLong
  eaUrbanLong = eaDatLong$urban
  
  # generate how many women and childen are in each cluster
  numWomenLong = rep(0, length(eaUrbanLong))
  numWomenLong[eaUrbanLong] = recdf(sum(eaUrbanLong), fastDistributions$mothersUrban)
  numWomenLong[!eaUrbanLong] = recdf(sum(!eaUrbanLong), fastDistributions$mothersRural)
  eaDatLong$numWomen = numWomenLong
  
  # now expand the eaDat table to be in the extra long format, were each row is a mother
  rowsExtraLong = rep(1:length(eaUrbanLong), numWomenLong)
  eaUrbanExtraLong = eaUrbanLong[rowsExtraLong]
  
  # generate the number of children per mother
  numChildrenExtraLong = rep(0, length(eaUrbanExtraLong))
  numChildrenExtraLong[eaUrbanExtraLong] = recdf(sum(eaUrbanExtraLong), fastDistributions$childrenUrban)
  numChildrenExtraLong[!eaUrbanExtraLong] = recdf(sum(!eaUrbanExtraLong), fastDistributions$childrenRural)
  
  # aggregate up from mother to household level
  numChildrenLong = rep(0, length(eaUrbanLong))
  numChildrenLong[numWomenLong!=0] = aggregate(numChildrenExtraLong, list(rowsExtraLong), sum)$x
  eaDatLong$numChildren = numChildrenLong
  
  ### generate Binomial probabilities from transformed logit scale GP
  # generate SPDE simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  print("Simulating nationwide mortality rates and data")
  if(margVar != 0) {
    SPDEArgs = list(coords=eaCoords, nsim=1, margVar=margVar, effRange=effRange, kenya=TRUE)
    simVals = do.call("simSPDE", SPDEArgs)
  } else {
    simVals = matrix(rep(0, nrow(eaCoords)), ncol=1)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  simValsNugLong = simValsNug[rowsLong]
  
  if(HHoldVar != 0) {
    simValsFinalLong = simValsNugLong + matrix(rnorm(length(simValsNugLong), sd=sqrt(HHoldVar)), ncol=1)
    probsFinalLong = expit(simValsFinalLong)
    eaDiedLong = matrix(rbinom(nrow(eaDatLong)*nsim, eaDatLong$numChildren, probsFinalLong), ncol=nsim)
  }
  else {
    probsFinalLong = expit(simValsNugLong)
    eaDiedLong = matrix(rbinom(nrow(eaDatLong)*nsim, eaDatLong$numChildren, rep(probsFinalLong, nsim)), ncol=nsim)
  }
  # transform back to original scale for the cluster level probabilities
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  probsNoNugLong = probsNoNug[rowsLong]
  probsNugLong = probs[rowsLong]
  
  eaDatLong$died = eaDiedLong
  eaDatLong$probsHH = probsFinalLong
  eaDatLong$trueProbDeathNoNug = probsNoNugLong
  eaDatLong$trueProbDeath = probsNugLong # technically, probsHH is the true probability at the household level
  
  # aggregate up from household to cluster level
  numChildren = aggregate(numChildrenLong, list(rowsLong), sum)$x
  numWomen = aggregate(numWomenLong, list(rowsLong), sum)$x
  eaDied = aggregate(c(eaDiedLong), list(rowsLong), sum)$x
  
  eaDat = as.data.frame(eaDat)
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numChildren
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    # clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
    clustDat = simClustersEmpirical(eaDat, eaDatLong, nsim, NULL, urbanOverSamplefrac, nHHSampled)
  }
  
  # return simulated data
  print("finishing up...")
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIsLong(eaDat, clustDat$eaIs, eaDatLong, clustDat$HHIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# SPDE or LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# empiricalDistributions: the set of empirical distributions for households per cluster, 
#                         mothers per household, and children per mother stratified by urban
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDatShort: enumeration area data set in "short form" were each row is a cluster. 
#             "long form" is where each row is a household.
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the spatial process, excluding household end cluster effects. 
#          If 0, no spatial component is included
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# simPopOnly: Do we just care about the simulated population, or also the cluster surveys?
simDatLCPB = function(nsim=1, margVar=0.243, tausq=0.463, 
                      gamma=0.009, effRange=406.51, beta0=-3.922, 
                      urbanOverSamplefrac=0, seed=NULL, inla.seed=0L, 
                      nHHSampled=25, fixPopPerEA=NULL, fixHHPerEA=NULL, fixPopPerHH=NULL, 
                      easpa=NULL, popMat=NULL, targetPopMat=NULL, 
                      includeUrban=TRUE, pixelLevel=TRUE, 
                      constituencyLevel=TRUE, countyLevel=TRUE, 
                      regionLevel=TRUE, nationalLevel=TRUE, 
                      doLcpb=TRUE, doLCpb=TRUE, doLCPb=TRUE, doIHME=TRUE, poppsub=poppcon, 
                      min1PerSubarea=TRUE, clustDat=NULL, 
                      spreadEAsInPixels=FALSE, logisticApproximation=TRUE, 
                      simPopOnly=FALSE, returnEAinfo=!simPopOnly, verbose=TRUE, 
                      stopOnFrameMismatch=TRUE, thisclustpc=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(!simPopOnly && nsim != 1) {
    stop("simulating multiple populations and surveys (nsim != 1) not yet supported if simPopOnly is FALSE")
  }
  
  ## set default inputs
  # construct default household and population count per county table
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA()
    # area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
    # HHUrb: the number of households in the urban part of the area
    # HHRur: the number of households in the rural part of the area
    # HHTotal: the number of households in the the area
    # popUrb: the number of people in the urban part of the area
    # popRur: the number of people in the rural part of the area
    # popTotal: the number of people in the the area
  }
  if(!is.null(fixPopPerEA)) {
    easpa$popUrb = easpa$EAUrb * fixPopPerEA
    easpa$popRur = easpa$EARur * fixPopPerEA
    easpa$popTotal = easpa$EATotal * fixPopPerEA
  }
  if(!is.null(fixHHPerEA)) {
    easpa$HHUrb = easpa$EAUrb * fixHHPerEA
    easpa$HHRur = easpa$EARur * fixHHPerEA
    easpa$HHTotal = easpa$EATotal * fixHHPerEA
  }
  if(!is.null(fixPopPerHH)) {
    if(fixPopPerHH * fixHHPerEA != fixPopPerEA) {
      stop("fixPopPerHH * fixHHPerEA != fixPopPerEA")
    }
  }
  totalEAs = sum(easpa$EATotal)
  totalHouseholds = sum(easpa$HHTotal)
  
  # construct overall population density surface
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  
  # construct target population density surface, adjusted based on known or estimated stratum populations (or "faux" stratum populations)
  if(is.null(targetPopMat)) {
    targetPopMat = adjustPopGrid(popMat, poppcAdjusted=easpa)
  }
  
  # set the number of sampled clusters based on DHS and available constituencies in the sampling frame
  if(is.null(thisclustpc)) {
    thisclustpc = clustpc
    usedCounties = unique(easpa$area)
    thisclustpc = thisclustpc[thisclustpc$County %in% usedCounties,]
  }
  
  ### generate Binomial probabilities from transformed logit scale GP
  # generate SPDE simulations
  pixelCoords = cbind(popMat$east, popMat$north)
  print("Simulating spatial risk")
  if(margVar != 0) {
    SPDEArgs = list(coords=pixelCoords, nsim=nsim, margVar=margVar, effRange=effRange, kenya=TRUE, inla.seed=inla.seed)
    simVals = do.call("simSPDE", SPDEArgs)
  } else {
    simVals = matrix(rep(0, nrow(eaCoords)), ncol=nsim)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*popMat$urban, "+")
  
  # simulate nugget/cluster effect
  epsc = matrix(rnorm(totalEAs * nsim, sd=sqrt(tausq)), ncol=nsim)
  
  # transform back to original scale for the pixel level probabilities
  probsNoNug = expit(simVals)
  
  # simulate the enumeration areas
  uDraws = simVals
  sigmaEpsilonDraws = rep(sqrt(tausq), nsim)
  
  print("Simulating population: EA locations, aggregate risk and prevalence")
  outLCPB = modLCPB(uDraws=uDraws, sigmaEpsilonDraws=sigmaEpsilonDraws, 
                    easpa=easpa, popMat=popMat, targetPopMat=targetPopMat, empiricalDistributions=NULL, 
                    includeUrban=includeUrban, pixelLevel=pixelLevel, 
                    constituencyLevel=constituencyLevel, countyLevel=countyLevel, 
                    regionLevel=regionLevel, nationalLevel=nationalLevel, doModifiedPixelLevel=FALSE, 
                    doLcpb=doLcpb, doLCpb=doLCpb, doLCPb=doLCPb, doIHME=doIHME, poppsub=poppsub, 
                    min1PerSubarea=min1PerSubarea, urbanEffect=gamma, 
                    returnEAinfo=returnEAinfo, epsc=epsc, fixPopPerEA=fixPopPerEA, 
                    fixHHPerEA=fixHHPerEA, fixPopPerHH=fixPopPerHH, 
                    logisticApproximation=logisticApproximation, verbose=verbose, 
                    stopOnFrameMismatch=stopOnFrameMismatch)
  
  if(returnEAinfo) {
    eaDatList = outLCPB$eaDatList
    eaSamples = outLCPB$eaSamples
  }
  
  # spread EAs uniformly within each 5km x 5km pixel
  if(spreadEAsInPixels && returnEAinfo) {
    for(i in 1:nsim) {
      eaDat = eaDatList[[i]]
      eaDat$east = jitter(eaDat$east, amount=2.5)
      eaDat$north = jitter(eaDat$north, amount=2.5)
      lonlat = projKenya(eaDat$east, eaDat$north, inverse=TRUE)
      eaDat$lon = lonlat[,1]
      eaDat$lat = lonlat[,2]
      
      eaDatList[[i]] = eaDat
    }
  }
  
  if(!simPopOnly && returnEAinfo) {
    if(nsim != 1) {
      stop("simulating multiple populations and surveys (nsim != 1) not supported if simPopOnly is FALSE")
    }
    
    print("Simulating surveys")
    eaDat = eaDatList[[1]]
    
    ### Simulate household survey if necessary
    # first generate the number of households
    numHouseholds = eaDat$nHH
    
    # now expand the eaDat table to be in the long format, were each row is a house
    rowsLong = rep(1:nrow(eaDat), numHouseholds)
    eaDatLong = eaDat[rowsLong, ]
    eaDatLong$eaIs = rowsLong
    eaUrbanLong = eaDatLong$urban
    
    # function for randomly spreading people among households in long form data:
    extendThisDat = function(xs, nHH) {
      revenMultinom = function(sizeK) {
        size = sizeK[1]
        k = sizeK[2]
        prob = rep(1/k, k)
        rmultinom(1, size, prob)
      }
      unlist(apply(cbind(xs, nHH), 1, revenMultinom))
    }
    
    # generate how many of the target population are in each cluster
    if(is.null(fixPopPerHH)) {
      lived = extendThisDat(eaDat$N - eaDat$Z, numHouseholds)
      died = extendThisDat(eaDat$Z, numHouseholds)
    } else if(fixPopPerHH == 1) {
      extendDatEven = function(Ns, Zs) {
        
        spreadAmongHHs = function(thisRow) {
          thisN = thisRow[1]
          thisZ = thisRow[2]
          
          # spread population evenly among households
          # nHH = thisRow$nHH
          # hhI = sample(1:nHH, nHH, replace=FALSE)
          c(rep(1, thisZ), rep(0, thisN - thisZ))
        }
        c(unlist(apply(cbind(Ns, Zs), 1, spreadAmongHHs)))
      }
      died = extendDatEven(eaDat$N, eaDat$Z)
      lived = 1 - died
    } else {
      stop("If fixPopPerHH is not NULL it must be 1")
    }
    
    eaDatLong$n = died + lived
    eaDatLong$y = died
    eaDatLong$nHH = 1
    eaDatLong$pLCPB = eaDatLong$y / eaDatLong$n
    eaDatLong$pLCPB[eaDatLong$n == 0] = eaDatLong$pLCPb[eaDatLong$n == 0]
    
    ### sample clusters and households within EAs
    # first generate clusters
    # if(is.null(clustDat))
    #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
    if(is.null(clustDat)) {
      print("simulating cluster locations:")
      # clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
      clustDat = simClustersEmpirical(eaDat, eaDatLong, nsim, NULL, urbanOverSamplefrac, nHHSampled, thisclustpc=thisclustpc)
    }
    
    # return simulated data
    print("finishing up...")
    
    # return cluster data in Andrea's format:
    clustList = genAndreaFormatFromEAIsLong(eaDat, clustDat$eaIs, eaDatLong, clustDat$HHIs, clustDat$sampleWeights)
    
    list(eaDat=eaDat, eaSamples=eaSamples, clustDat=clustList, aggregatedPop=outLCPB, thisclustpc=thisclustpc)
  } else if(returnEAinfo) {
    list(eaDatList=eaDatList, eaSamples=eaSamples, aggregatedPop=outLCPB, thisclustpc=thisclustpc)
  } else {
    outLCPB
  }
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# SPDE or LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# empiricalDistributions: the set of empirical distributions for households per cluster, 
#                         mothers per household, and children per mother stratified by urban
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDatShort: enumeration area data set in "short form" were each row is a cluster. 
#             "long form" is where each row is a household.
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the spatial process, excluding household end cluster effects. 
#          If 0, no spatial component is included
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# simPopOnly: Do we just care about the simulated population, or also the cluster surveys?
simDatLCPB2 = function(nsim=1, margVar=0.243, sigmaEpsilon=sqrt(0.463), 
                      gamma=0.009, effRange=406.51, beta0=-3.922, 
                      urbanOverSamplefrac=0, seed=NULL, inla.seed=0L, 
                      nHHSampled=25, fixPopPerEA=NULL, fixHHPerEA=NULL, fixPopPerHH=NULL, 
                      easpa=NULL, popMat=NULL, targetPopMat=NULL, 
                      stratifyByUrban=TRUE, gridLevel=TRUE, subareaLevel=TRUE, 
                      doSmoothRisk=TRUE, doFineScaleRisk=TRUE, 
                      poppsub=poppsubKenya, min1PerSubarea=TRUE, 
                      spreadEAsInPixels=FALSE, logisticApproximation=TRUE, 
                      simPopOnly=FALSE, returnEAinfo=!simPopOnly, verbose=TRUE, 
                      stopOnFrameMismatch=TRUE, thisclustpc=NULL, 
                      nEAsFac=1, nClustFac=1, representativeSampling=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ## set default inputs
  # construct default household and population count per county table
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA()
    # area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
    # HHUrb: the number of households in the urban part of the area
    # HHRur: the number of households in the rural part of the area
    # HHTotal: the number of households in the the area
    # popUrb: the number of people in the urban part of the area
    # popRur: the number of people in the rural part of the area
    # popTotal: the number of people in the the area
  }
  
  # scale easpa and poppsub using nEAsFac
  easpa[,c("EAUrb", "EARur", "EATotal", 
           "HHUrb", "HHRur", "HHTotal", 
           "popUrb", "popRur", "popTotal")] = 
    round(nEAsFac * easpa[,c("EAUrb", "EARur", "EATotal", 
                       "HHUrb", "HHRur", "HHTotal", 
                       "popUrb", "popRur", "popTotal")])
  easpa$EATotal = easpa$EAUrb + easpa$EARur
  easpa$HHTotal = easpa$HHUrb + easpa$HHRur
  easpa$popTotal = easpa$popUrb + easpa$popRur
  
  if(!is.null(poppsub)) {
    poppsub[,c("popUrb", "popRur", "popTotal")] = 
      nEAsFac * poppsub[,c("popUrb", "popRur", "popTotal")]
    poppsub$popTotal = poppsub$popUrb + poppsub$popRur
  } else if(nEAsFac != 1) {
    stop("nEAsFac is not 1, but poppsub is NULL")
  }
    # adjust easpa based on fixed number of people/household or households/EA if specified
  if(!is.null(fixPopPerEA)) {
    easpa$popUrb = easpa$EAUrb * fixPopPerEA
    easpa$popRur = easpa$EARur * fixPopPerEA
    easpa$popTotal = easpa$EATotal * fixPopPerEA
  }
  if(!is.null(fixHHPerEA)) {
    easpa$HHUrb = easpa$EAUrb * fixHHPerEA
    easpa$HHRur = easpa$EARur * fixHHPerEA
    easpa$HHTotal = easpa$EATotal * fixHHPerEA
  }
  if(!is.null(fixPopPerHH)) {
    if(fixPopPerHH * fixHHPerEA != fixPopPerEA) {
      stop("fixPopPerHH * fixHHPerEA != fixPopPerEA")
    }
  }
  totalEAs = sum(easpa$EATotal)
  totalHouseholds = sum(easpa$HHTotal)
  
  # construct overall population density surface up to a scalar factor
  if(is.null(popMat)) {
    popMat = makeDefaultPopMat()
    # lon: longitude
    # lat: latitude
    # east: easting (km)
    # north: northing (km)
    # pop: proportional to population density for each grid cell
    # area: an id or area name in which the grid cell corresponding to each row resides
    # urban: whether the grid cell is urban or rural
  }
  kmres = getPopMatResolution(popMat)
  
  # construct target population density surface, adjusted based on known or estimated stratum populations (or "faux" stratum populations)
  if(is.null(targetPopMat)) {
    targetPopMat = adjustPopGrid(popMat, poppcAdjusted=easpa)
  }
  
  # set the number of sampled clusters based on DHS and available constituencies in the sampling frame
  if(is.null(thisclustpc)) {
    thisclustpc = clustpc
    usedCounties = unique(easpa$area)
    thisclustpc = thisclustpc[thisclustpc$area %in% usedCounties,]
    
    thisclustpc[,c("clustUrb", "clustRur", "clustTotal")] = 
      round(nClustFac * thisclustpc[,c("clustUrb", "clustRur", "clustTotal")])
    thisclustpc$clustTotal = thisclustpc$clustUrb + thisclustpc$clustRur
    thisclustpc$HHUrb = thisclustpc$clustUrb * nHHSampled
    thisclustpc$HHRur = thisclustpc$clustRur * nHHSampled
    thisclustpc$HHTotal = thisclustpc$HHUrb + thisclustpc$HHRur
  }
  
  # ensure we're not sampling more clusters than EAs in any stratum
  thisclustpc = thisclustpc[order(thisclustpc$area),]
  easpa = easpa[order(easpa$area),]
  if(any(thisclustpc$clustUrb > easpa$EAUrb) | 
     any(thisclustpc$clustRur > easpa$EARur) | 
     any(thisclustpc$clustTotal > easpa$EATotal)) {
    warning("In some areas the number of requested clusters to sample is larger than the number of EAs. Setting equal to the number of EAs")
    # first make sure number of sampled EAs in individual strata are within limits
    thisclustpc$clustUrb[thisclustpc$clustUrb > easpa$EAUrb] = easpa$EAUrb[thisclustpc$clustUrb > easpa$EAUrb]
    thisclustpc$clustRur[thisclustpc$clustRur > easpa$EARur] = easpa$EARur[thisclustpc$clustRur > easpa$EARur]
    thisclustpc$clustTotal = thisclustpc$clustUrb + thisclustpc$clustRur
    
    # now make sure county totals are within limits. If not, scale strata in county appropriately
    if(any(thisclustpc$clustTotal > easpa$EATotal)) {
      badCounties = thisclustpc$clustTotal > easpa$EATotal
      badTotals = thisclustpc$clustTotal[badCounties]
      badTargets = easpa$EATotal[badCounties]
      thisclustpc$clustUrb[badCounties] = floor(thisclustpc$clustUrb[badCounties] * badTargets / badTotals)
      thisclustpc$clustRur[badCounties] = floor(thisclustpc$clustUrb[badCounties] * badTargets / badTotals)
      thisclustpc$clustTotal = thisclustpc$clustUrb + thisclustpc$clustRur
    }
  }
  
  ### generate Binomial probabilities from transformed logit scale GP
  # generate SPDE simulations
  pixelCoords = cbind(popMat$east, popMat$north)
  
  outLCPB = simPopSPDE(nsim=nsim, easpa=easpa, popMat=popMat, targetPopMat=targetPopMat, 
                       poppsub=poppsub, spdeMesh=spdeMesh, 
                       margVar=margVar, sigmaEpsilon=sigmaEpsilon, 
                       gamma=gamma, effRange=effRange, beta0=beta0, 
                       seed=seed, inla.seed=inla.seed, nHHSampled=nHHSampled, 
                       stratifyByUrban=stratifyByUrban, 
                       subareaLevel=subareaLevel, gridLevel=gridLevel, 
                       doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                       doSmoothRiskLogisticApprox=logisticApproximation, 
                       min1PerSubarea=min1PerSubarea, 
                       fixPopPerEA=fixPopPerEA, fixHHPerEA=fixHHPerEA, fixPopPerHH=fixPopPerHH)
  
  if(returnEAinfo) {
    eaDatList = outLCPB$eaPop$eaDatList
    eaSamples = outLCPB$eaPop$eaSamples
  }
  
  # spread EAs uniformly within each 5km x 5km pixel
  if(spreadEAsInPixels && returnEAinfo) {
    for(i in 1:nsim) {
      eaDat = eaDatList[[i]]
      eaDat$east = jitter(eaDat$east, amount=kmres/2)
      eaDat$north = jitter(eaDat$north, amount=kmres/2)
      lonlat = projKenya(eaDat$east, eaDat$north, inverse=TRUE)
      eaDat$lon = lonlat[,1]
      eaDat$lat = lonlat[,2]
      
      eaDatList[[i]] = eaDat
    }
  }
  
  if(!simPopOnly && returnEAinfo) {
    print("Simulating surveys")
    
    surveys = list()
    for(i in 1:nsim) {
      if(i %% 10 == 0) {
        print(paste0("simulating survey ", i, "/", nsim))
      }
      eaDat = eaDatList[[i]]
      
      ### Simulate household survey if necessary
      # first generate the number of households
      numHouseholds = eaDat$nHH
      
      # now expand the eaDat table to be in the long format, were each row is a house
      rowsLong = rep(1:nrow(eaDat), numHouseholds)
      eaDatLong = eaDat[rowsLong, ]
      eaDatLong$eaIs = rowsLong
      eaUrbanLong = eaDatLong$urban
      
      # function for randomly spreading people among households in long form data:
      extendThisDat = function(xs, nHH) {
        revenMultinom = function(sizeK) {
          size = sizeK[1]
          k = sizeK[2]
          prob = rep(1/k, k)
          rmultinom(1, size, prob)
        }
        unlist(apply(cbind(xs, nHH), 1, revenMultinom))
      }
      
      # generate how many of the target population are in each cluster
      if(is.null(fixPopPerHH)) {
        lived = extendThisDat(eaDat$N - eaDat$Z, numHouseholds)
        died = extendThisDat(eaDat$Z, numHouseholds)
      } else if(fixPopPerHH == 1) {
        extendDatEven = function(Ns, Zs) {
          
          spreadAmongHHs = function(thisRow) {
            thisN = thisRow[1]
            thisZ = thisRow[2]
            
            # spread population evenly among households
            # nHH = thisRow$nHH
            # hhI = sample(1:nHH, nHH, replace=FALSE)
            c(rep(1, thisZ), rep(0, thisN - thisZ))
          }
          c(unlist(apply(cbind(Ns, Zs), 1, spreadAmongHHs)))
        }
        died = extendDatEven(eaDat$N, eaDat$Z)
        lived = 1 - died
      } else {
        stop("If fixPopPerHH is not NULL it must be 1")
      }
      
      eaDatLong$n = died + lived
      eaDatLong$y = died
      eaDatLong$nHH = 1
      eaDatLong$pFineScalePrevalence = eaDatLong$y / eaDatLong$n
      eaDatLong$pFineScalePrevalence[eaDatLong$n == 0] = 0
      
      ### sample clusters and households within EAs
      # first generate clusters
      # clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
      # clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
      clustDat = simClustersEmpirical(eaDat, eaDatLong, nsim=1, NULL, urbanOverSamplefrac, nHHSampled, 
                                      thisclustpc=thisclustpc, representativeSampling=representativeSampling, 
                                      verbose=FALSE)
      
      # return cluster data in Andrea's format:
      clustList = genAndreaFormatFromEAIsLong2(eaDat, clustDat$eaIs, eaDatLong, clustDat$HHIs, 
                                               clustDat$sampleWeights, doFineScaleRisk=doFineScaleRisk, 
                                               doSmoothRisk=doSmoothRisk, doGriddedRisk=FALSE)[[1]]
      clustList$Z = clustList$y
      clustList$N = clustList$n
      
      surveys = c(surveys, list(clustList))
    }
    
    list(eaDatList=eaDatList, eaSamples=eaSamples, clustDat=surveys, aggregatedPop=outLCPB, thisclustpc=thisclustpc)
  } else if(returnEAinfo) {
    list(eaDatList=eaDatList, eaSamples=eaSamples, aggregatedPop=outLCPB, thisclustpc=thisclustpc)
  } else {
    outLCPB
  }
}


simNsFull = function(n=1, includeUrban=TRUE, easpa=makeDefaultEASPA(), empiricalDistributions=NULL) {
  if(is.null(empiricalDistributions)) {
    out = load(paste0(globalDirectory, "empiricalDistributions.RData"))
    # list(households=householdDistribution, mothers=motherDistribution, children=childrenDistribution,
    #      householdsUrban=householdDistributionUrban, mothersUrban=motherDistributionUrban, childrenUrban=childrenDistributionUrban,
    #      householdsRural=householdDistributionRural, mothersRural=motherDistributionRural, childrenRural=childrenDistributionRural)
    # May also have popUrban, popRural, and popTotal, in which case they are used
    
    # make edfun versions of the ecdfs
    fastDistributions = list()
    for(i in 1:length(empiricalDistributions)) {
      fastDistributions = c(fastDistributions, list(ecdf2edfun(empiricalDistributions[[i]])))
    }
    names(fastDistributions) = names(empiricalDistributions)
  }
  
  # NOTE1: we assume the values for each EA within urban/rural boundaries are iid
  if(includeUrban) {
    if(!is.null(fastDistributions$popUrban)) {
      NcsUrban = matrix(recdf(sum(easpa$EAUrb)*n, fastDistributions$popUrban), ncol=n)
      NcsRural = matrix(recdf(sum(easpa$EARur)*n, fastDistributions$popRural), ncol=n)
    } else {
      NcsUrban = matrix(recdfComposed(sum(easpa$EAUrb)*n, list(fastDistributions$householdsUrban, 
                                                               fastDistributions$mothersUrban, 
                                                               fastDistributions$childrenUrban)), 
                        ncol=n)
      NcsRural = matrix(recdfComposed(sum(easpa$EARur)*n, list(fastDistributions$householdsRural, 
                                                               fastDistributions$mothersRural, 
                                                               fastDistributions$childrenRural)), 
                        ncol=n)
    }
    list(NcsUrban=NcsUrban, NcsRural=NcsRural)
  } else {
    if(!is.null(fastDistributions$popTotal)) {
      Ncs = matrix(recdf(sum(easpa$EATotal)*n, fastDistributions$popTotal), ncol=n)
    } else {
      Ncs = matrix(recdfComposed(sum(easpa$EATotal)*n, list(fastDistributions$households, 
                                                            fastDistributions$mothers, 
                                                            fastDistributions$children)), 
                   ncol=n)
    }
    Ncs
  }
}

# Same as simDatLK, but uses SPDE model to simulate data set instead of LatticeKrig.
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
simDatSPDE = function(eaDat, clustDat=NULL, nsim=1, margVar=1, effRange=300, 
                      beta0=0, gamma=-1, tausq=1, 
                      womenFrac=1/3, urbanOverSample=1, numClusters=423, seed=NULL, fullEADat=NULL, 
                      HHoldVar=0, fixedWomenPerClust=TRUE, plotProbs=FALSE, 
                      savePlots=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### first generate Binomial probabilities from transformed logit scale GP
  ## generate SPDE simulations
  # generate mesh grid
  eaCoords = cbind(eaDat$east, eaDat$north)
  mesh = getSPDEMeshGrid(eaCoords, doPlot = FALSE)
  
  print("Simulating nationwide mortality rates and data")
  SPDEArgs = list(coords=eaCoords, nsim=1, margVar=margVar, effRange=effRange)
  simVals = do.call("simSPDE", SPDEArgs)
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  # plot first simulation if requested by user on both logit and linear scales
  if(plotProbs) {
    if(savePlots)
      pdf("figures/spdeSimTest.pdf", width=8, height=5)
    
    par(mfrow=c(1, 2))
    quilt.plot(eaCoords, simValsNug[, 1], main ="Simulated Logit Mortality Rates", 
               xlab="Easting", ylab="Northing")
    plotMapDat(project=TRUE)
    
    quilt.plot(eaCoords, probs[, 1], main ="Simulated Mortality Rates", 
               xlab="Easting", ylab="Northing")
    plotMapDat(project=TRUE)
    
    if(savePlots)
      dev.off()
  }
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
  }
  
  # return simulated data
  print("finishing up...")
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIs(eaDat, clustDat$eaIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

## TODO: add in multiple seeds, one population for each
runSimStudy = function(gamma=0, rho=(1/3)^2, sigmaEpsilon=sqrt(1/2.5), 
                       effRange=400, beta0=-3.9, nEAsFac=1, nClustFac=1, 
                       seed=123, i=NULL, jRange=1:100, 
                       clust=NULL, nClustProc=16) {
  set.seed(seed)
  
  if(!is.null(i)) {
    out = load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
    
    theseArgs = spde_prevRiskSimStudyCommandArgs[[i]]
    beta0 = theseArgs$beta0
    gamma = theseArgs$gamma
    rho = theseArgs$rho
    sigmaEpsilon = theseArgs$sigmaEpsilon
    effRange = theseArgs$effRange
    nEAsFac = theseArgs$nEasFac
    nClustFac = theseArgs$nClustFac
    
    allSeeds = matrix(sample(1:100000000, length(spde_prevRiskSimStudyCommandArgs)*100, replace=FALSE), ncol=100)
  }
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Bet", signif(beta0, 3), "rho", signif(rho, 3), "sigEps", 
                  signif(sigmaEpsilon, 3), "gam", signif(gamma, 3), 
                  "nEAsFac", signif(nEAsFac, 3), "nClustFac", signif(nClustFac, 3))
  
  # load data (overSampDat, SRSDat)
  out = load(paste0(dataSaveDirectory, "simData", dataID, ".RData"))
  clustDat = stratDat$clustDatList[j]
  
  # TODO:
  # make easpa based on nEAsFac
  # make other necessary inputs to simPopSPDE
  
  # make easpa based on nEAsFac
  easpa = makeDefaultEASPA()
  
  # scale easpa and poppsub using nEAsFac
  easpa[,c("EAUrb", "EARur", "EATotal", 
           "HHUrb", "HHRur", "HHTotal", 
           "popUrb", "popRur", "popTotal")] = 
    round(nEAsFac * easpa[,c("EAUrb", "EARur", "EATotal", 
                             "HHUrb", "HHRur", "HHTotal", 
                             "popUrb", "popRur", "popTotal")])
  easpa$EATotal = easpa$EAUrb + easpa$EARur
  easpa$HHTotal = easpa$HHUrb + easpa$HHRur
  easpa$popTotal = easpa$popUrb + easpa$popRur
  
  poppsub = poppsubKenya
  poppsub[,c("popUrb", "popRur", "popTotal")] = 
    nEAsFac * poppsub[,c("popUrb", "popRur", "popTotal")]
  poppsub$popTotal = poppsub$popUrb + poppsub$popRur
  
  # make parallel cluster, and export relevant values
  if(is.null(clust)) {
    library(parallel)
    cl <- makeCluster(nClustProc)
    clusterEvalQ(cl, {
      setwd("~/git/continuousNugget/")
      source("code/setup.R")
    })
    clusterExport(cl, c("allSeeds", "beta0", "gamma", "rho", "sigmaEpsilon", 
                        "effRange", "nEAsFac", "nClustFac", "easpa", "poppsub"))
  }
  
  simOne = function(j) {
    
    thisSeed = allSeeds[i,j]
    set.seed(thisSeed)
    
    # run the models
    out = simPopSPDE(nsim=1000, easpa=easpa, popMat=popMatKenya, 
                     targetPopMat=popMatKenyaNeonatal, poppsub=poppsub, 
                     spdeMesh=NULL, margVar=rho, sigmaEpsilon=sigmaEpsilon, 
                     gamma=gamma, effRange=effRange, beta0=beta0, 
                     stratifyByUrban=TRUE, subareaLevel=TRUE, gridLevel=FALSE, 
                     doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
                     doGriddedRisk=FALSE, doSmoothRiskLogisticApprox=FALSE, 
                     min1PerSubarea=TRUE)
    
    subareaPop = out$subareaPop$aggregationResults
    areaPop = out$areaPop$aggregationResults
    timings = out$timings
    
    # subset results we care about (save burden to scratch directory)
    subareaPopP = subareaPop[c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk")]
    areaPopP = areaPop[c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk")]
    subareaPopZ = subareaPop[c("ZFineScalePrevalence", "ZFineScaleRisk", "ZSmoothRisk")]
    areaPopZ = areaPop[c("ZFineScalePrevalence", "ZFineScaleRisk", "ZSmoothRisk")]
    
    # save file
    dataIDout = paste0(dataID, "_i", i, "j", j)
    save(subareaPopP, areaPopP, file=paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_p.RData"))
    save(subareaPopZ, areaPopZ, file=paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_Z.RData"))
    
    list(subareaPopAggRes=subareaPop, areaPopAggRes=areaPop)
  }
  
}

runSimStudyij = function(i, j, seed=123, doSmoothRiskLogisticApprox=FALSE, coarse=TRUE, doGC=FALSE) {
  time1 = proc.time()[3]
  set.seed(seed)
  
  out = load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
  theseArgs = spde_prevRiskSimStudyCommandArgs[[i]]
  beta0 = theseArgs$beta0
  gamma = theseArgs$gamma
  rho = theseArgs$rho
  sigmaEpsilon = theseArgs$sigmaEpsilon
  effRange = theseArgs$effRange
  nEAsFac = theseArgs$nEAsFac
  nClustFac = theseArgs$nClustFac
  
  allSeeds = matrix(sample(1:100000000, length(spde_prevRiskSimStudyCommandArgs)*100, replace=FALSE), ncol=100)
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Bet", signif(beta0, 3), "rho", signif(rho, 3), "sigEps", 
                  signif(sigmaEpsilon, 3), "gam", signif(gamma, 3), 
                  "nEAsFac", signif(nEAsFac, 3), "nClustFac", signif(nClustFac, 3))
  
  # load data (overSampDat, SRSDat)
  dataSaveDirectory = "savedOutput/simDataSets/"
  out = load(paste0(dataSaveDirectory, "simData", dataID, ".RData"))
  clustDat = stratDat$clustDatList[[j]]
  
  # TODO:
  # make easpa based on nEAsFac
  # make other necessary inputs to simPopSPDE
  
  # make easpa based on nEAsFac
  easpa = makeDefaultEASPA()
  
  # scale easpa and poppsub using nEAsFac
  easpa[,c("EAUrb", "EARur", "EATotal", 
           "HHUrb", "HHRur", "HHTotal", 
           "popUrb", "popRur", "popTotal")] = 
    round(nEAsFac * easpa[,c("EAUrb", "EARur", "EATotal", 
                             "HHUrb", "HHRur", "HHTotal", 
                             "popUrb", "popRur", "popTotal")])
  easpa$EATotal = easpa$EAUrb + easpa$EARur
  easpa$HHTotal = easpa$HHUrb + easpa$HHRur
  easpa$popTotal = easpa$popUrb + easpa$popRur
  
  poppsub = poppsubKenya
  poppsub[,c("popUrb", "popRur", "popTotal")] = 
    nEAsFac * poppsub[,c("popUrb", "popRur", "popTotal")]
  poppsub$popTotal = poppsub$popUrb + poppsub$popRur
  
  thisSeed = allSeeds[i,j]
  set.seed(thisSeed)
  
  if(coarse) {
    popMat = popGridCoarse
    popMatAdjusted = popGridCoarseAdjusted
  } else {
    popMat = popGrid
    popMatAdjusted = popGridAdjusted
  }
  
  # Fit the risk model:
  time2 = proc.time()[3]
  riskOut = fitSPDEKenyaDat(dat=clustDat, nPostSamples=1000, popMat=popMat)
  logitDraws = riskOut$uDraws # uDraws includes fixed effects and is on logit scale
  sigmaEpsilonDraws = riskOut$sigmaEpsilonDraws
  if(doGC) {
    rm(riskOut)
    gc(TRUE)
  }
  time3 = proc.time()[3]
  
  # run the models
  out = simPopCustom(logitRiskDraws=logitDraws, sigmaEpsilonDraws=sigmaEpsilonDraws, 
                     easpa=easpa, popMat=popMat, 
                     targetPopMat=popMatAdjusted, poppsub=poppsub, 
                     stratifyByUrban=TRUE, subareaLevel=TRUE, gridLevel=FALSE, 
                     doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
                     doGriddedRisk=FALSE, doSmoothRiskLogisticApprox=doSmoothRiskLogisticApprox, 
                     min1PerSubarea=TRUE, doGC=doGC)
  
  subareaPop = out$subareaPop$aggregationResults
  areaPop = out$areaPop$aggregationResults
  allTimings = out$allTimings
  processedTimings = out$processedTimings
  aggregationTimings=list(allTimings=allTimings, processedTimings=processedTimings)
  time4 = proc.time()[3]
  
  # get timings
  rawTimes = c(time1, time2, time3, time4)
  totalTimes = diff(rawTimes)
  totalTimes = c(totalTimes, sum(totalTimes))
  names(totalTimes) = c("setup", "SPDEmodel", "aggregationModel", "totalTime")
  
  # subset results we care about (save burden to scratch directory)
  subareaPopP = subareaPop[c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk", 
                             "pUrbanFineScalePrevalence", "pUrbanFineScaleRisk", "pUrbanSmoothRisk", 
                             "pRuralFineScalePrevalence", "pRuralFineScaleRisk", "pRuralSmoothRisk")]
  areaPopP = areaPop[c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk", 
                       "pUrbanFineScalePrevalence", "pUrbanFineScaleRisk", "pUrbanSmoothRisk", 
                       "pRuralFineScalePrevalence", "pRuralFineScaleRisk", "pRuralSmoothRisk")]
  subareaPopZ = subareaPop[c("ZFineScalePrevalence", "ZFineScaleRisk", "ZSmoothRisk", 
                             "ZUrbanFineScalePrevalence", "ZUrbanFineScaleRisk", "ZUrbanSmoothRisk", 
                             "ZRuralFineScalePrevalence", "ZRuralFineScaleRisk", "ZRuralSmoothRisk")]
  areaPopZ = areaPop[c("ZFineScalePrevalence", "ZFineScaleRisk", "ZSmoothRisk", 
                       "ZUrbanFineScalePrevalence", "ZUrbanFineScaleRisk", "ZUrbanSmoothRisk", 
                       "ZRuralFineScalePrevalence", "ZRuralFineScaleRisk", "ZRuralSmoothRisk")]
  
  # save files
  coarseText = ifelse(coarse, "Coarse", "")
  dataIDout = paste0("simOut_i", i, "j", j, coarseText)
  save(subareaPopP, areaPopP, aggregationTimings, rawTimes, totalTimes, file=paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_p.RData"))
  save(subareaPopZ, areaPopZ, aggregationTimings, rawTimes, totalTimes, file=paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_Z.RData"))
  
  list(subareaPopAggRes=subareaPop, areaPopAggRes=areaPop, 
       aggregationTimings=aggregationTimings, 
       rawTimes=rawTimes, totalTimes=totalTimes)
}

processSimStudyResultsij = function(i, j, coarse=TRUE) {
  
  out = load("savedOutput/simStudyResults/spde_prevRiskSimStudyCommandArgs.RData")
  theseArgs = spde_prevRiskSimStudyCommandArgs[[i]]
  beta0 = theseArgs$beta0
  gamma = theseArgs$gamma
  rho = theseArgs$rho
  sigmaEpsilon = theseArgs$sigmaEpsilon
  effRange = theseArgs$effRange
  nEAsFac = theseArgs$nEAsFac
  nClustFac = theseArgs$nClustFac
  
  allSeeds = matrix(sample(1:100000000, length(spde_prevRiskSimStudyCommandArgs)*100, replace=FALSE), ncol=100)
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Bet", signif(beta0, 3), "rho", signif(rho, 3), "sigEps", 
                  signif(sigmaEpsilon, 3), "gam", signif(gamma, 3), 
                  "nEAsFac", signif(nEAsFac, 3), "nClustFac", signif(nClustFac, 3))
  
  # load data (overSampDat, SRSDat)
  dataSaveDirectory = "savedOutput/simDataSets/"
  out = load(paste0(dataSaveDirectory, "simData", dataID, ".RData"))
  # clustDat = stratDat$clustDatList[[j]]
  thisPop = stratDat$aggregatedPop
  
  # load results
  coarseText = ifelse(coarse, "Coarse", "")
  dataIDout = paste0("simOut_i", i, "j", j, coarseText)
  # load(subareaPopP, areaPopP, aggregationTimings, rawTimes, totalTimes, file=paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_p.RData"))
  # load(subareaPopZ, areaPopZ, aggregationTimings, rawTimes, totalTimes, file=paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_Z.RData"))
  load(paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_p.RData"))
  load(paste0("savedOutput/simStudyResults/tempFiles/", dataIDout, "_Z.RData"))
  
  ## compare results to population via scores
  
  # first obtain central estimates, which are the same for all models
  subareaUrbanEstsP = rowMeans(subareaPopP$pUrbanSmoothRisk)
  subareaUrbanEstsZ = rowMeans(subareaPopZ$ZUrbanSmoothRisk)
  subareaRuralEstsP = rowMeans(subareaPopP$pRuralSmoothRisk)
  subareaRuralEstsZ = rowMeans(subareaPopZ$ZRuralSmoothRisk)
  subareaEstsP = rowMeans(subareaPopP$pSmoothRisk)
  subareaEstsZ = rowMeans(subareaPopZ$ZSmoothRisk)
  areaEstsP = rowMeans(areaPopP$pSmoothRisk)
  areaEstsZ = rowMeans(areaPopZ$ZSmoothRisk)
  
  # fine scale prevalence scores
  subareaStratScoresPprev = getScores(truth=c(thisPop$subareaPop$aggregationResults$pUrbanFineScalePrevalence[,j], 
                                              thisPop$subareaPop$aggregationResults$pRuralFineScalePrevalence[,j]), 
                                      est=c(subareaUrbanEstsP, subareaRuralEstsP), estMat=rbind(subareaPopP$pUrbanFineScalePrevalence, subareaPopP$pRuralFineScalePrevalence), 
                                      significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  subareaStratScoresZprev = getScores(truth=c(thisPop$subareaPop$aggregationResults$ZUrbanFineScalePrevalence[,j], 
                                              thisPop$subareaPop$aggregationResults$ZRuralFineScalePrevalence[,j]), 
                                      est=c(subareaUrbanEstsZ, subareaRuralEstsZ), estMat=rbind(subareaPopZ$ZUrbanFineScalePrevalence, subareaPopZ$ZRuralFineScalePrevalence), 
                                      significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  
  subareaScoresPprev = getScores(truth=thisPop$subareaPop$aggregationResults$pFineScalePrevalence[,j], 
                                 est=subareaEstsP, estMat=subareaPopP$pFineScalePrevalence, 
                             significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  subareaScoresZprev = getScores(truth=thisPop$subareaPop$aggregationResults$ZFineScalePrevalence[,j], 
                                 est=subareaEstsZ, estMat=subareaPopZ$ZFineScalePrevalence, 
                             significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  # every once in a while a subarea x stratum will have nonzero population but zero deaths. Hence, throwOutAllNAs must be set to TRUE for sensible scores
  subareaScoresPrelprev = getScores(truth=thisPop$subareaPop$aggregationResults$pUrbanFineScalePrevalence[,j]/
                                   thisPop$subareaPop$aggregationResults$pRuralFineScalePrevalence[,j], 
                                 estMat=subareaPopP$pUrbanFineScalePrevalence/subareaPopP$pRuralFineScalePrevalence, 
                                 significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  
  areaScoresPprev = getScores(truth=thisPop$areaPop$aggregationResults$pFineScalePrevalence[,j], 
                              est=areaEstsP, estMat=areaPopP$pFineScalePrevalence, 
                             significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  areaScoresZprev = getScores(truth=thisPop$areaPop$aggregationResults$ZFineScalePrevalence[,j], 
                              est=areaEstsZ, estMat=areaPopZ$ZFineScalePrevalence, 
                             significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  areaScoresPrelprev = getScores(truth=thisPop$areaPop$aggregationResults$pUrbanFineScalePrevalence[,j]/thisPop$areaPop$aggregationResults$pRuralFineScalePrevalence[,j], 
                              estMat=areaPopP$pUrbanFineScalePrevalence/areaPopP$pRuralFineScalePrevalence, 
                              significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE)
  
  # fine scale risk scores
  subareaStratScoresPrisk = getScores(truth=c(thisPop$subareaPop$aggregationResults$pUrbanFineScalePrevalence[,j], 
                                                    thisPop$subareaPop$aggregationResults$pRuralFineScalePrevalence[,j]), 
                                            est=c(subareaUrbanEstsP, subareaRuralEstsP), estMat=rbind(subareaPopP$pUrbanFineScaleRisk, subareaPopP$pRuralFineScaleRisk), 
                                            significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  subareaStratScoresZrisk = getScores(truth=c(thisPop$subareaPop$aggregationResults$ZUrbanFineScalePrevalence[,j], 
                                                    thisPop$subareaPop$aggregationResults$ZRuralFineScalePrevalence[,j]), 
                                            est=c(subareaUrbanEstsZ, subareaRuralEstsZ), estMat=rbind(subareaPopZ$ZUrbanFineScaleRisk, subareaPopZ$ZRuralFineScaleRisk), 
                                            significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  
  subareaScoresPrisk = getScores(truth=thisPop$subareaPop$aggregationResults$pFineScalePrevalence[,j], 
                                 est=subareaEstsP, estMat=subareaPopP$pFineScaleRisk, 
                                 significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  subareaScoresZrisk = getScores(truth=thisPop$subareaPop$aggregationResults$ZFineScalePrevalence[,j], 
                                 est=subareaEstsZ, estMat=subareaPopZ$ZFineScaleRisk, 
                                 significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  subareaScoresPrelrisk = getScores(truth=thisPop$subareaPop$aggregationResults$pUrbanFineScalePrevalence[,j]/thisPop$subareaPop$aggregationResults$pRuralFineScalePrevalence[,j], 
                                 estMat=subareaPopP$pUrbanFineScaleRisk/subareaPopP$pRuralFineScaleRisk, 
                                 significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE)
  
  areaScoresPrisk = getScores(truth=thisPop$areaPop$aggregationResults$pFineScalePrevalence[,j], 
                              est=areaEstsP, estMat=areaPopP$pFineScaleRisk, 
                              significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  areaScoresZrisk = getScores(truth=thisPop$areaPop$aggregationResults$ZFineScalePrevalence[,j], 
                              est=areaEstsZ, estMat=areaPopZ$ZFineScaleRisk, 
                              significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  areaScoresPrelrisk = getScores(truth=thisPop$areaPop$aggregationResults$pUrbanFineScalePrevalence[,j]/thisPop$areaPop$aggregationResults$pRuralFineScalePrevalence[,j], 
                              estMat=areaPopP$pUrbanFineScaleRisk/areaPopP$pRuralFineScaleRisk, 
                              significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE)
  
  # smooth risk scores
  subareaStratScoresPsmoothRisk = getScores(truth=c(thisPop$subareaPop$aggregationResults$pUrbanFineScalePrevalence[,j], 
                                                    thisPop$subareaPop$aggregationResults$pRuralFineScalePrevalence[,j]), 
                                      est=c(subareaUrbanEstsP, subareaRuralEstsP), estMat=rbind(subareaPopP$pUrbanSmoothRisk, subareaPopP$pRuralSmoothRisk), 
                                      significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  subareaStratScoresZsmoothRisk = getScores(truth=c(thisPop$subareaPop$aggregationResults$ZUrbanFineScalePrevalence[,j], 
                                                    thisPop$subareaPop$aggregationResults$ZRuralFineScalePrevalence[,j]), 
                                            est=c(subareaUrbanEstsZ, subareaRuralEstsZ), estMat=rbind(subareaPopZ$ZUrbanSmoothRisk, subareaPopZ$ZRuralSmoothRisk), 
                                            significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE, throwOutAllNAs=TRUE)
  
  subareaScoresPsmoothRisk = getScores(truth=thisPop$subareaPop$aggregationResults$pFineScalePrevalence[,j], 
                                       est=subareaEstsP, estMat=subareaPopP$pSmoothRisk, 
                                 significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  subareaScoresZsmoothRisk = getScores(truth=thisPop$subareaPop$aggregationResults$ZFineScalePrevalence[,j], 
                                       est=subareaEstsZ, estMat=subareaPopZ$ZSmoothRisk, 
                                 significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  subareaScoresPrelsmoothRisk = getScores(truth=thisPop$subareaPop$aggregationResults$pUrbanFineScalePrevalence[,j]/thisPop$subareaPop$aggregationResults$pRuralFineScalePrevalence[,j], 
                                       estMat=subareaPopP$pUrbanSmoothRisk/subareaPopP$pRuralSmoothRisk, 
                                       significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE)
  
  areaScoresPsmoothRisk = getScores(truth=thisPop$areaPop$aggregationResults$pFineScalePrevalence[,j], 
                                    est=areaEstsP, estMat=areaPopP$pSmoothRisk, 
                              significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  areaScoresZsmoothRisk = getScores(truth=thisPop$areaPop$aggregationResults$ZFineScalePrevalence[,j], 
                                    est=areaEstsZ, estMat=areaPopZ$ZSmoothRisk, 
                              significance=c(.8, .9, .95), doFuzzyReject=TRUE)
  areaScoresPrelsmoothRisk = getScores(truth=thisPop$areaPop$aggregationResults$pUrbanFineScalePrevalence[,j]/thisPop$areaPop$aggregationResults$pRuralFineScalePrevalence[,j], 
                                    estMat=areaPopP$pUrbanSmoothRisk/areaPopP$pRuralSmoothRisk, 
                                    significance=c(.8, .9, .95), doFuzzyReject=TRUE, na.rm=TRUE, setInfToNA=TRUE)
  
  # add computation times to the scores
  subareaStratScoresPprev = cbind(subareaStratScoresPprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  subareaStratScoresZprev = cbind(subareaStratScoresZprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  subareaScoresPprev = cbind(subareaScoresPprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  subareaScoresZprev = cbind(subareaScoresZprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  subareaScoresPrelprev = cbind(subareaScoresPrelprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  areaScoresPprev = cbind(areaScoresPprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  areaScoresZprev = cbind(areaScoresZprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  areaScoresPrelprev = cbind(areaScoresPrelprev, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[1])/60)
  
  subareaStratScoresPrisk = cbind(subareaStratScoresPrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  subareaStratScoresZrisk = cbind(subareaStratScoresZrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  subareaScoresPrisk = cbind(subareaScoresPrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  subareaScoresZrisk = cbind(subareaScoresZrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  subareaScoresPrelrisk = cbind(subareaScoresPrelrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  areaScoresPrisk = cbind(areaScoresPrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  areaScoresZrisk = cbind(areaScoresZrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  areaScoresPrelrisk = cbind(areaScoresPrelrisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[2])/60)
  
  subareaStratScoresPsmoothRisk = cbind(subareaStratScoresPsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  subareaStratScoresZsmoothRisk = cbind(subareaStratScoresZsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  subareaScoresPsmoothRisk = cbind(subareaScoresPsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  subareaScoresZsmoothRisk = cbind(subareaScoresZsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  subareaScoresPrelsmoothRisk = cbind(subareaScoresPrelsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  areaScoresPsmoothRisk = cbind(areaScoresPsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  areaScoresZsmoothRisk = cbind(areaScoresZsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  areaScoresPrelsmoothRisk = cbind(areaScoresPrelsmoothRisk, Time=(totalTimes[1] + totalTimes[2] + aggregationTimings$processedTimings[3])/60)
  
  save(subareaStratScoresPprev, subareaStratScoresZprev, 
       subareaScoresPprev, subareaScoresZprev, subareaScoresPrelprev, 
       areaScoresPprev, areaScoresZprev, areaScoresPrelprev, 
       subareaStratScoresPrisk, subareaStratScoresZrisk, 
       subareaScoresPrisk, subareaScoresZrisk, subareaScoresPrelrisk, 
       areaScoresPrisk, areaScoresZrisk, areaScoresPrelrisk, 
       subareaStratScoresPsmoothRisk, subareaStratScoresZsmoothRisk, 
       subareaScoresPsmoothRisk, subareaScoresZsmoothRisk, subareaScoresPrelsmoothRisk, 
       areaScoresPsmoothRisk, areaScoresZsmoothRisk, areaScoresPrelsmoothRisk, 
       aggregationTimings, totalTimes, 
       file=paste0("savedOutput/simStudyResults/tempFiles/simScores_i", i, "j", j, coarseText, ".RData"))
  
  res = list(subareaStratScoresPprev, subareaStratScoresZprev, 
             subareaScoresPprev, subareaScoresZprev, subareaScoresPrelprev, 
             areaScoresPprev, areaScoresZprev, areaScoresPrelprev, 
             subareaStratScoresPrisk, subareaStratScoresZrisk, 
             subareaScoresPrisk, subareaScoresZrisk, subareaScoresPrelrisk, 
             areaScoresPrisk, areaScoresZrisk, areaScoresPrelrisk, 
             subareaStratScoresPsmoothRisk, subareaStratScoresZsmoothRisk, 
             subareaScoresPsmoothRisk, subareaScoresZsmoothRisk, subareaScoresPrelsmoothRisk, 
             areaScoresPsmoothRisk, areaScoresZsmoothRisk, areaScoresPrelsmoothRisk, 
             aggregationTimings, totalTimes)
  names(res) = c("subareaStratScoresPprev", "subareaStratScoresZprev", 
                 "subareaScoresPprev", "subareaScoresZprev", "subareaScoresPrelprev", 
                 "areaScoresPprev", "areaScoresZprev", "areaScoresPrelprev", 
                 "subareaStratScoresPrisk", "subareaStratScoresZrisk", 
                 "subareaScoresPrisk", "subareaScoresZrisk", "subareaScoresPrelrisk", 
                 "areaScoresPrisk", "areaScoresZrisk", "areaScoresPrelrisk", 
                 "subareaStratScoresPsmoothRisk", "subareaStratScoresZsmoothRisk", 
                 "subareaScoresPsmoothRisk", "subareaScoresZsmoothRisk", "subareaScoresPrelsmoothRisk", 
                 "areaScoresPsmoothRisk", "areaScoresZsmoothRisk", "areaScoresPrelsmoothRisk", 
                 "aggregationTimings", "totalTimes")
  invisible(res)
}

combineProcessedResults = function(is=1:54, maxJ=100, initialProcess=TRUE, combineScores=TRUE, coarse=TRUE) {
  
  if(initialProcess) {
    print("Doing initial processing/scoring...")
    for(i in is) {
      for(j in 1:maxJ) {
        if((j %% 10) == 1) {
          print(paste0("processing for i=", i, ", j=", j))
        }
        processSimStudyResultsij(i, j, coarse=coarse)
      }
    }
  }
  
  if(combineScores) {
    print("Combining scores...")
    for(i in is) {
      
      subareaStratScoresPprevAll = c()
      subareaStratScoresZprevAll = c()
      subareaScoresPprevAll = c()
      subareaScoresZprevAll = c()
      subareaScoresPrelprevAll = c()
      areaScoresPprevAll = c()
      areaScoresZprevAll = c()
      areaScoresPrelprevAll = c()
      
      subareaStratScoresPriskAll = c()
      subareaStratScoresPrelriskAll = c()
      subareaScoresPriskAll = c()
      subareaScoresPrelriskAll = c()
      subareaScoresZriskAll = c()
      areaScoresPriskAll = c()
      areaScoresZriskAll = c()
      areaScoresPrelriskAll = c()
      
      subareaStratScoresPsmoothRiskAll = c()
      subareaStratScoresZsmoothRiskAll = c()
      subareaScoresPsmoothRiskAll = c()
      subareaScoresZsmoothRiskAll = c()
      subareaScoresPrelsmoothRiskAll = c()
      areaScoresPsmoothRiskAll = c()
      areaScoresZsmoothRiskAll = c()
      areaScoresPrelsmoothRiskAll = c()
      
      aggregationTimingsDetailedAll = c()
      aggregationTimingsProcessedAll = c()
      totalTimesAll = c()
      
      for(j in 1:maxJ) {
        if((j %% 10) == 1) {
          print(paste0("combining scores for i=", i, ", j=", j))
        }
        coarseText = ifelse(coarse, "Coarse", "")
        out = load(paste0("savedOutput/simStudyResults/tempFiles/simScores_i", i, "j", j, coarseText, ".RData"))
        
        subareaStratScoresPprevAll = rbind(subareaStratScoresPprevAll, subareaStratScoresPprev)
        subareaStratScoresZprevAll = rbind(subareaStratScoresZprevAll, subareaStratScoresZprev)
        subareaScoresPprevAll = rbind(subareaScoresPprevAll, subareaScoresPprev)
        subareaScoresZprevAll = rbind(subareaScoresZprevAll, subareaScoresZprev)
        subareaScoresPrelprevAll = rbind(subareaScoresPrelprevAll, subareaScoresPrelprev)
        areaScoresPprevAll = rbind(areaScoresPprevAll, areaScoresPprev)
        areaScoresZprevAll = rbind(areaScoresZprevAll, areaScoresZprev)
        areaScoresPrelprevAll = rbind(areaScoresPrelprevAll, areaScoresPrelprev)
        
        subareaStratScoresPriskAll = rbind(subareaStratScoresPriskAll, subareaStratScoresPrisk)
        subareaStratScoresZriskAll = rbind(subareaStratScoresZriskAll, subareaStratScoresZrisk)
        subareaScoresPriskAll = rbind(subareaScoresPriskAll, subareaScoresPrisk)
        subareaScoresZriskAll = rbind(subareaScoresZriskAll, subareaScoresZrisk)
        subareaScoresPrelriskAll = rbind(subareaScoresPrelriskAll, subareaScoresPrelrisk)
        areaScoresPriskAll = rbind(areaScoresPriskAll, areaScoresPrisk)
        areaScoresZriskAll = rbind(areaScoresZriskAll, areaScoresZrisk)
        areaScoresPrelriskAll = rbind(areaScoresPrelriskAll, areaScoresPrelrisk)
        
        subareaStratScoresPsmoothRiskAll = rbind(subareaStratScoresPsmoothRiskAll, subareaStratScoresPsmoothRisk)
        subareaStratScoresZsmoothRiskAll = rbind(subareaStratScoresZsmoothRiskAll, subareaStratScoresZsmoothRisk)
        subareaScoresPsmoothRiskAll = rbind(subareaScoresPsmoothRiskAll, subareaScoresPsmoothRisk)
        subareaScoresZsmoothRiskAll = rbind(subareaScoresZsmoothRiskAll, subareaScoresZsmoothRisk)
        subareaScoresPrelsmoothRiskAll = rbind(subareaScoresPrelsmoothRiskAll, subareaScoresPrelsmoothRisk)
        areaScoresPsmoothRiskAll = rbind(areaScoresPsmoothRiskAll, areaScoresPsmoothRisk)
        areaScoresZsmoothRiskAll = rbind(areaScoresZsmoothRiskAll, areaScoresZsmoothRisk)
        areaScoresPrelsmoothRiskAll = rbind(areaScoresPrelsmoothRiskAll, areaScoresPrelsmoothRisk)
        
        aggregationTimingsDetailedAll = rbind(aggregationTimingsDetailedAll, aggregationTimings$allTimings[,1])
        aggregationTimingsProcessedAll = rbind(aggregationTimingsProcessedAll, aggregationTimings$processedTimings)
        totalTimesAll = rbind(totalTimesAll, totalTimes)
      }
      colnames(aggregationTimingsDetailedAll) = row.names(aggregationTimings$allTimings)
      
      # calculate average scores
      subareaStratScoresPprevAvg = colMeans(as.matrix(subareaStratScoresPprevAll))
      subareaStratScoresZprevAvg = colMeans(as.matrix(subareaStratScoresZprevAll))
      subareaScoresPprevAvg = colMeans(as.matrix(subareaScoresPprevAll))
      subareaScoresZprevAvg = colMeans(as.matrix(subareaScoresZprevAll))
      subareaScoresPrelprevAvg = colMeans(as.matrix(subareaScoresPrelprevAll), na.rm=TRUE)
      areaScoresPprevAvg = colMeans(as.matrix(areaScoresPprevAll))
      areaScoresZprevAvg = colMeans(as.matrix(areaScoresZprevAll))
      areaScoresPrelprevAvg = colMeans(as.matrix(areaScoresPrelprevAll), na.rm=TRUE)
      
      subareaStratScoresPriskAvg = colMeans(as.matrix(subareaStratScoresPriskAll))
      subareaStratScoresZriskAvg = colMeans(as.matrix(subareaStratScoresZriskAll))
      subareaScoresPriskAvg = colMeans(as.matrix(subareaScoresPriskAll))
      subareaScoresZriskAvg = colMeans(as.matrix(subareaScoresZriskAll))
      subareaScoresPrelriskAvg = colMeans(as.matrix(subareaScoresPrelriskAll), na.rm=TRUE)
      areaScoresPriskAvg = colMeans(as.matrix(areaScoresPriskAll))
      areaScoresZriskAvg = colMeans(as.matrix(areaScoresZriskAll))
      areaScoresPrelriskAvg = colMeans(as.matrix(areaScoresPrelriskAll), na.rm=TRUE)
      
      subareaStratScoresPsmoothRiskAvg = colMeans(as.matrix(subareaStratScoresPsmoothRiskAll))
      subareaStratScoresZsmoothRiskAvg = colMeans(as.matrix(subareaStratScoresZsmoothRiskAll))
      subareaScoresPsmoothRiskAvg = colMeans(as.matrix(subareaScoresPsmoothRiskAll))
      subareaScoresZsmoothRiskAvg = colMeans(as.matrix(subareaScoresZsmoothRiskAll))
      subareaScoresPrelsmoothRiskAvg = colMeans(as.matrix(subareaScoresPrelsmoothRiskAll), na.rm=TRUE)
      areaScoresPsmoothRiskAvg = colMeans(as.matrix(areaScoresPsmoothRiskAll))
      areaScoresZsmoothRiskAvg = colMeans(as.matrix(areaScoresZsmoothRiskAll))
      areaScoresPrelsmoothRiskAvg = colMeans(as.matrix(areaScoresPrelsmoothRiskAll), na.rm=TRUE)
      
      # calculate average timing of individual parts of the code
      aggregationTimingsDetailedAvg = colMeans(aggregationTimingsDetailedAll)
      aggregationTimingsProcessedAvg = colMeans(aggregationTimingsProcessedAll)
      totalTimesAvg = colMeans(totalTimesAll)
      
      save(subareaStratScoresPprevAll, subareaStratScoresZprevAll, 
           subareaScoresPprevAll, subareaScoresZprevAll, subareaScoresPrelprevAll, 
           areaScoresPprevAll, areaScoresZprevAll, areaScoresPrelprevAll, 
           subareaStratScoresPriskAll, subareaStratScoresZriskAll, 
           subareaScoresPriskAll, subareaScoresZriskAll, subareaScoresPrelriskAll, 
           areaScoresPriskAll, areaScoresZriskAll, areaScoresPrelriskAll, 
           subareaStratScoresPsmoothRiskAll, subareaStratScoresZsmoothRiskAll, 
           subareaScoresPsmoothRiskAll, subareaScoresZsmoothRiskAll, subareaScoresPrelsmoothRiskAll, 
           areaScoresPsmoothRiskAll, areaScoresZsmoothRiskAll, areaScoresPrelsmoothRiskAll, 
           aggregationTimingsDetailedAll, aggregationTimingsProcessedAll, totalTimesAll, 
           subareaStratScoresPprevAvg, subareaStratScoresZprevAvg, 
           subareaScoresPprevAvg, subareaScoresZprevAvg, subareaScoresPrelprevAvg, 
           areaScoresPprevAvg, areaScoresZprevAvg, areaScoresPrelprevAvg, 
           subareaStratScoresPriskAvg, subareaStratScoresZriskAvg, 
           subareaScoresPriskAvg, subareaScoresZriskAvg, subareaScoresPrelriskAvg, 
           areaScoresPriskAvg, areaScoresZriskAvg, areaScoresPrelriskAvg, 
           subareaStratScoresPsmoothRiskAvg, subareaStratScoresZsmoothRiskAvg, 
           subareaScoresPsmoothRiskAvg, subareaScoresZsmoothRiskAvg, subareaScoresPrelsmoothRiskAvg, 
           areaScoresPsmoothRiskAvg, areaScoresZsmoothRiskAvg, areaScoresPrelsmoothRiskAvg, 
           aggregationTimingsDetailedAvg, aggregationTimingsProcessedAvg, totalTimesAvg, 
           file=paste0("savedOutput/simStudyResults/simScoresAll_i", i, "maxJ", maxJ, coarseText, ".RData"))
    }
  }
}

printSimStudyTablei = function(i=1, maxJ=100, coarse=TRUE) {
  # save(subareaScoresPprevAll, subareaScoresZprevAll, areaScoresPprevAll, areaScoresZprevAll, 
  #      subareaScoresPriskAll, subareaScoresZriskAll, areaScoresPriskAll, areaScoresZriskAll, 
  #      subareaScoresPsmoothRiskAll, subareaScoresZsmoothRiskAll, areaScoresPsmoothRiskAll, areaScoresZsmoothRiskAll, 
  #      aggregationTimingsDetailedAll, aggregationTimingsProcessedAll, totalTimesAll, 
  #      subareaScoresPprevAvg, subareaScoresZprevAvg, areaScoresPprevAvg, areaScoresZprevAvg, 
  #      subareaScoresPriskAvg, subareaScoresZriskAvg, areaScoresPriskAvg, areaScoresZriskAvg, 
  #      subareaScoresPsmoothRiskAvg, subareaScoresZsmoothRiskAvg, areaScoresPsmoothRiskAvg, areaScoresZsmoothRiskAvg, 
  #      aggregationTimingsDetailedAvg, aggregationTimingsProcessedAvg, totalTimesAvg, 
  #      file=paste0("savedOutput/simStudyResults/simScoresAll_i", i, "maxJ", maxJ, ".RData"))
  coarseText = ifelse(coarse, "Coarse", "")
  out = load(paste0("savedOutput/simStudyResults/simScoresAll_i", i, "maxJ", maxJ, coarseText, ".RData"))
  
  # make 4 tables: subareaPScores, subareaZScores, areaPScores, and areaZScores
  
  subareaPScores = rbind(subareaScoresPprevAvg, 
                         subareaScoresPriskAvg, 
                         subareaScoresPsmoothRiskAvg)
  row.names(subareaPScores) = c("Prevalence", "Risk", "SmoothRisk")
  
  subareaZScores = rbind(subareaScoresZprevAvg, 
                         subareaScoresZriskAvg, 
                         subareaScoresZsmoothRiskAvg)
  row.names(subareaZScores) = c("Prevalence", "Risk", "SmoothRisk")
  
  subareaPrelScores = rbind(subareaScoresPrelprevAvg, 
                         subareaScoresPrelriskAvg, 
                         subareaScoresPrelsmoothRiskAvg)
  row.names(subareaPrelScores) = c("Prevalence", "Risk", "SmoothRisk")
  
  areaPScores = rbind(areaScoresPprevAvg, 
                         areaScoresPriskAvg, 
                         areaScoresPsmoothRiskAvg)
  row.names(areaPScores) = c("Prevalence", "Risk", "SmoothRisk")
  
  areaZScores = rbind(areaScoresZprevAvg, 
                         areaScoresZriskAvg, 
                         areaScoresZsmoothRiskAvg)
  row.names(areaZScores) = c("Prevalence", "Risk", "SmoothRisk")
  
  areaPrelScores = rbind(areaScoresPrelprevAvg, 
                      areaScoresPrelriskAvg, 
                      areaScoresPrelsmoothRiskAvg)
  row.names(areaPrelScores) = c("Prevalence", "Risk", "SmoothRisk")
  
  # multiply coverages by 100 to be in percent units, remove variance and MSE
  coverageCols = grep("Coverage", colnames(subareaPScores))
  varCols = grep("Var", colnames(subareaPScores))
  mseCols = grep("MSE", colnames(subareaPScores[,-varCols]))
  rmseCols = grep("RMSE", colnames(subareaPScores[,-varCols]))
  mseCols = mseCols[-which(mseCols == rmseCols)]
  
  subareaPScores[,coverageCols] = 100 * subareaPScores[,coverageCols]
  subareaPScores = subareaPScores[,-varCols]
  subareaPScores = subareaPScores[,-mseCols]
  
  subareaZScores[,coverageCols] = 100 * subareaZScores[,coverageCols]
  subareaZScores = subareaZScores[,-varCols]
  subareaZScores = subareaZScores[,-mseCols]
  
  subareaPrelScores[,coverageCols] = 100 * subareaPrelScores[,coverageCols]
  subareaPrelScores = subareaPrelScores[,-varCols]
  subareaPrelScores = subareaPrelScores[,-mseCols]
  
  areaPScores[,coverageCols] = 100 * areaPScores[,coverageCols]
  areaPScores = areaPScores[,-varCols]
  areaPScores = areaPScores[,-mseCols]
  
  areaZScores[,coverageCols] = 100 * areaZScores[,coverageCols]
  areaZScores = areaZScores[,-varCols]
  areaZScores = areaZScores[,-mseCols]
  
  areaPrelScores[,coverageCols] = 100 * areaPrelScores[,coverageCols]
  areaPrelScores = areaPrelScores[,-varCols]
  areaPrelScores = areaPrelScores[,-mseCols]
  
  # print tables
  displayP = c("s", rep("e", 3), rep("e", 3), rep("d", 3), rep("f", 4))
  digitsP = c(0, rep(-2, 3), rep(2, 3), rep(0, 3), rep(3, 3), 1)
  displayZ = c("s", rep("f", 3), rep("f", 3), rep("d", 3), rep("f", 4))
  digitsZ = c(0, rep(1, 2), 2, rep(1, 3), rep(0, 3), rep(1, 3), 1)
  displayPrel = c("s", rep("f", 3), rep("f", 3), rep("d", 3), rep("f", 4))
  digitsPrel = c(0, rep(2, 2), 2, rep(2, 3), rep(0, 3), rep(3, 3), 1)
  print(xtable(subareaPScores, display=displayP, digits=digitsP), math.style.exponents=TRUE)
  print(xtable(subareaZScores, display=displayZ, digits=digitsZ), math.style.exponents=TRUE)
  print(xtable(subareaPrelScores, display=displayPrel, digits=digitsPrel), math.style.exponents=TRUE)
  
  print(xtable(areaPScores, display=displayP, digits=digitsP), math.style.exponents=TRUE)
  print(xtable(areaZScores, display=displayZ, digits=digitsZ), math.style.exponents=TRUE)
  print(xtable(areaPrelScores, display=displayPrel, digits=digitsPrel), math.style.exponents=TRUE)
}

generateJobList = function(workDir="savedOutput/simStudyResults/tempFiles/", iRange=1:54, jRange=1:100, extensiveCheck=FALSE, extensiveNACheck=FALSE) {
  # save current directory to return to later. Set directory to job file locations
  thisDir = getwd()
  setwd(workDir)
  
  # get indices of relevant characters in the file names
  out = system("ls simOut_*Coarse_p.RData", intern=TRUE)
  iIs = rep(8, length(out))
  jIs = sapply(out, function(x) {gregexpr("j", x)[[1]][1]})
  CIs = sapply(out, function(x) {gregexpr("C", x)[[1]][1]})
  
  # determine whether the relevant file exists or not
  fileExists = matrix(FALSE, nrow=length(iRange), ncol=length(jRange))
  for(i in 1:length(out)) {
    thisFilename = out[i]
    iI = iIs[i]
    jI = jIs[i]
    CI = CIs[i]
    
    thisI = as.numeric(substr(thisFilename, iI+1, jI-1))
    thisJ = as.numeric(substr(thisFilename, jI+1, CI-1))
    if(!extensiveCheck) {
      fileExists[thisI, thisJ] = TRUE
    } else {
      # check to make sure we can actually load the relevant files
      thisFilenameP = thisFilename
      thisFilenameZ = gsub("_p", "_Z", thisFilenameP)
      
      canLoad <<- TRUE
      tmp = tryCatch(load(thisFilenameP), error= function(e) {canLoad <<- FALSE})
      tmp = tryCatch(load(thisFilenameZ), error= function(e) {canLoad <<- FALSE})
      
      if(extensiveNACheck) {
        if(all(is.na(subareaPopP$pFineScalePrevalence))) 
          canLoad <<- FALSE
      }
      
      fileExists[thisI, thisJ] = canLoad
    }
  }
  missingJobInds = (1:length(fileExists))[!c(t(fileExists))]
  
  if(length(missingJobInds) == 0) {
    setwd(thisDir)
    print("No missing jobs")
    return(invsible(NULL))
  }else if(length(missingJobInds) == 1) {
    setwd(thisDir)
    return(missingJobInds)
  } else {
    print(paste0("Total missing job results: ", length(missingJobInds)))
  }
  
  # shorten the string to make it readable:
  
  # add a fake missing job index at the end to make sure we get the final missing jobs
  missingJobInds = c(missingJobInds, length(fileExists) + 100)
  
  # look for groups of consecutive missing job indices
  lastInd = missingJobInds[1]
  firstJobInd = lastInd
  jobList = ""
  for(i in 2:length(missingJobInds)) {
    thisInd = missingJobInds[i]
    thisDiff = thisInd - lastInd
    commaStr = ifelse(firstJobInd == missingJobInds[1], "", ",")
    
    if(thisDiff != 1) {
      # this index and last index aren't in the same group, must add 
      # last index (or its associated group) to the job list string
      
      if(firstJobInd == lastInd) {
        # last index was a singleton
        jobList = paste0(jobList, commaStr, lastInd)
      } else {
        # last index was the last consecutive job in the group
        jobList = paste0(jobList, commaStr, "[", firstJobInd, "-", lastInd, "]")
      }
      
      firstJobInd = thisInd
    } else {
      # this index and last index are in the same group. Don't need to 
      # do anything here
    }
    lastInd = thisInd
  }
  
  setwd(thisDir)
  
  return(jobList)
}





