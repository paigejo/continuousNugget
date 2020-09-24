# this script is for calculating direct estimates

modDirect = function(dat=NULL, dataType=c("mort", "ed"), significanceCI=.8, level=c("county", "stratum", "constituency", "constituencyStratum"), 
                     getContrast=FALSE, logisticApproximation=FALSE) {
  
  level = match.arg(level)
  dataType = match.arg(dataType)
  if(is.null(dat)) {
    if(dataType == "mort") {
      dat = mort
    }
    else {
      dat = ed
    }
  }
  
  # extend the dataset to binary form if necessary
  if(any(dat$n > 1)) {
    extendedDat = extendData(dat, dat$clusterID)
    dat = extendedDat
  }
  
  if(getContrast && level == "stratum") {
    stop("getContrast cannot be set to TRUE if level == 'stratum'")
  }
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  
  if(level %in% c("constituency", "constituencyStratum")) {
    constituencies = dat$admin2
    uniqueConstituencies = sort(unique(constituencies))
    regions = constituencyToCounty(uniqueConstituencies)
  } else {
    regions <- sort(unique(dat$admin1))
  }
  
  if(level == "county") {
    results<-data.frame(admin1=rep(regions,each=1))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  } else if(level == "stratum"){
    urbanResults<-data.frame(admin1=rep(regions,each=1))
    urbanResults$var.est<-urbanResults$logit.est<-urbanResults$upper<-urbanResults$lower<-urbanResults$est<-NA
    urbanResults$converge <- NA
    ruralResults<-data.frame(admin1=rep(regions,each=1))
    ruralResults$var.est<-ruralResults$logit.est<-ruralResults$upper<-ruralResults$lower<-ruralResults$est<-NA
    ruralResults$converge <- NA
  } else if(level == "constituency") {
    results<-data.frame(admin1=regions, admin2=uniqueConstituencies)
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  } else if(level == "constituencyStratum") {
    urbanResults<-data.frame(admin1=regions, admin2=uniqueConstituencies)
    urbanResults$var.est<-urbanResults$logit.est<-urbanResults$upper<-urbanResults$lower<-urbanResults$est<-NA
    urbanResults$converge <- NA
    ruralResults<-data.frame(admin1=regions, admin2=uniqueConstituencies)
    ruralResults$var.est<-ruralResults$logit.est<-ruralResults$upper<-ruralResults$lower<-ruralResults$est<-NA
    ruralResults$converge <- NA
  }
  
  stratVar = with(dat, interaction(admin1, urban), drop=TRUE)
  
  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the clusterID denote
    ##        one stage cluster design (clusterID is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~clusterID,
                              strata =NULL,
                              weights=NULL, data=dat)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    dat$strat <- stratVar
    
    # --- setting up the design object --- #
    ## NOTE: -the clusterID denote
    ##        one stage cluster design (clusterID is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~clusterID,
                              strata=~strat, nest=T, 
                              weights=~samplingWeight, data=dat)
  }
  
  nAreas = length(regions)
  for(i in 1:nAreas){
    if(level == "county") {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat, svydesign=my.svydesign, 
                                           area=results$admin1[i], 
                                           logisticApproximation=logisticApproximation)
    } else if(level == "stratum"){
      urbanResults[i, 2:7] <- region.time.HTStratumDat(dataobj=dat, svydesign=my.svydesign, 
                                                       area=urbanResults$admin1[i], urb = TRUE, 
                                                       logisticApproximation=logisticApproximation)
      ruralResults[i, 2:7] <- region.time.HTStratumDat(dataobj=dat, svydesign=my.svydesign, 
                                                       area=ruralResults$admin1[i], urb = FALSE, 
                                                       logisticApproximation=logisticApproximation)
    } else if(level == "constituency") {
      results[i, 3:8] <- region.time.HTDat(dataobj=dat, svydesign=my.svydesign, 
                                           area=results$admin2[i], areaVarName="admin2", 
                                           logisticApproximation=logisticApproximation)
    } else if(level == "constituencyStratum") {
      urbanResults[i, 3:8] <- region.time.HTStratumDat(dataobj=dat, svydesign=my.svydesign, 
                                                       area=urbanResults$admin2[i], urb = TRUE, areaVarName="admin2", 
                                                       logisticApproximation=logisticApproximation)
      ruralResults[i, 3:8] <- region.time.HTStratumDat(dataobj=dat, svydesign=my.svydesign, 
                                                       area=ruralResults$admin2[i], urb = FALSE, areaVarName="admin2", 
                                                       logisticApproximation=logisticApproximation)
    }
  }
  
  if(level %in% c("stratum", "constituencyStratum")) {
    return(list(urbanResults=urbanResults, ruralResults=ruralResults))
  } else if(level %in% c("county", "constituency")) {
    return(results)
  }
}

# function to extend dataset to binary form one row at a time
extendData <- function(clustDat, clusterID=clustDat$clusterID, divideWeight=TRUE, useNumChildrenDied=FALSE){
  clustDat$clusterID = clusterID
  clustDatRows = lapply(1:nrow(clustDat), function(i) {clustDat[i,]})
  do.call("rbind", lapply(clustDatRows, "extendDataRow", 
                                 divideWeight=divideWeight, useNumChildrenDied=useNumChildrenDied))
}

# function to extend dataset to binary form one row at a time
extendDataRow <- function(clustDatRow, clusterID=clustDatRow$clusterID, divideWeight=TRUE, useNumChildrenDied=FALSE){
  
  # add extra columns for ageMonth, ageGrpD, clusterID, v002
  if(useNumChildrenDied) {
    clustDatRow$n = clustDatRow$numChildren
    clustDatRow$y = clustDatRow$died
  }
  n = clustDatRow$n
  # tmp = data.frame(clustDatRow[c(1, 6:16)])
  
  # get everything that can be extended without modification
  unmodifiedI = match(c("clusterID", "regionRural", "region", "urban", "lon", "lat", "admin1", "admin2", "east", "north", "strat"), names(clustDatRow))
  unmodifiedI = unmodifiedI[!is.na(unmodifiedI)]
  # tmp = data.frame(clustDatRow[c(1, c(4, 6:ncol(clustDatRow)))])
  if(length(dim(clustDatRow)) >= 2) {
    tmp = data.frame(clustDatRow[1, unmodifiedI])
  } else {
    tmp = data.frame(clustDatRow[unmodifiedI])
  }
  
  # tmp$clusterID = clusterID
  
  ageMonth = rep(0, n)
  ageGrpD = rep("[0,1)", n)
  clusterID = rep(clusterID, n)
  # there is only one child and one mother per household.
  # All 25 households are sampled
  householdID = 1:n
  
  y = c(rep(0,n-clustDatRow$y), rep(1, clustDatRow$y))
  if(clustDatRow["urban"][1,1]){
    urbanRural = rep("urban", n)
  } else {
    urbanRural = rep("rural", n)
  }
  # admin1 = rep(clustDatRow$admin1, n)
  
  res = merge(data.frame(y, ageMonth, ageGrpD, clusterID, householdID, urbanRural), tmp, by="clusterID")
  
  # the below line was commented out since each cluster only has one type of admin and urban level. 
  # The equivalent line has been added into the parent function
  # res$regionRural <- with(res, interaction(admin1, urbanRural), drop=TRUE)
  
  if(divideWeight)
    res$samplingWeight = clustDatRow$samplingWeight / n
  return(res)
}

extendDataDat <- function(clustDatRow, v001, divideWeight=TRUE){
  
  # add extra columns for ageMonth, ageGrpD, v001, v002
  n = clustDatRow$n
  # the only things we need are admin1 and sampling weight, but we must get rid of 
  # urban, y, and the number of women since those will be recalculated
  # tmp = data.frame(clustDatRow[c(1, 6:16)])
  # tmp = data.frame(clustDatRow[c(1, c(4, 6:ncol(clustDatRow)))])
  tmp = data.frame(clustDatRow[c(1, c(4, 6:(ncol(clustDatRow) - 2)))])
  tmp$v001 = v001
  
  ageMonth = rep(0, n)
  # ageGrpD = rep("[0,1)", n)
  v001 = rep(v001, n)
  # there is only one child and one mother per household.
  # All 25 households are sampled
  v002 = 1:n
  
  y = c(rep(0,n-clustDatRow$y), rep(1, clustDatRow$y))
  if(clustDatRow["urban"][1,1]){
    urbanRural = rep("urban", n)
  } else {
    urbanRural = rep("rural", n)
  }
  # admin1 = rep(clustDatRow$admin1, n)
  
  # res = merge(data.frame(y, ageMonth, ageGrpD, v001, v002, urbanRural), tmp, by="v001")
  res = merge(data.frame(y, ageMonth, v001, v002, urbanRural), tmp, by="v001")
  
  # the below line was commented out since each cluster only has one type of admin and urban level. 
  # The equivalent line has been added into the parent function
  # res$regionRural <- with(res, interaction(admin1, urbanRural), drop=TRUE)
  
  if(divideWeight)
    res$samplingWeight = res$samplingWeight / n
  return(res)
}


# - a function that reads in a glm or svyglm - #
# - object and returns the estimate and SE - #
# - specifics in the supplementary materials - #
## This function takes care of the delta method
## to calculate the variance of u5m as a function
## of the age specific hazards, \beta_a .

get.est<-function(glm.ob, logisticApproximation=FALSE){
  
  beta<-summary(glm.ob)$coef[,1]
  
  var.est <- vcov(glm.ob)[1,1]
  if(var.est < 1e-15) {
    est <-expit(beta)
  } else {
    if(var.est > 1e15 || !is.finite(var.est)) {
      browser()
    }
    est <-logitNormMean(matrix(c(beta, sqrt(var.est)), nrow=1), logisticApproximation=logisticApproximation)
  }
  # est <-expit(beta)
  
  # compute 80% CI intervals
  lower <- logit(est)+qnorm(c(0.1))*sqrt(var.est)
  upper <- logit(est)+qnorm(c(0.9))*sqrt(var.est)
  return(c(est,lower, upper,logit(est),var.est))
}

# -- a function to subset the design based on a region and time period -- #
# -- and then run the svyglm function and return the get.est() results -- #

## First line in function allows you to subset your data and ALSO the specified
## svydesign object into area (usually v024 variable in DHS) 
## and time (per5 is a variable we construct for the 5-year periods in the Stata step)
## Second line fits the survey-weighted glm

region.time.HT<-function(dataobj, svydesign, area, areaVarName="admin1", logisticApproximation=FALSE){
  # first check if observations are all 0 or 1
  if(sum(dataobj$y[dataobj[[areaVarName]] == area]) == 0) {
    return(c(0, -Inf, Inf, -Inf, Inf, 0))
  } else if(sum(dataobj$y[dataobj[[areaVarName]] == area]) == sum(dataobj$n[dataobj[[areaVarName]] == area])) {
    return(c(1, -Inf, Inf, Inf, Inf, 0))
  }
  
  if(areaVarName == "admin1") {
    tmp<-subset(svydesign, (admin1==area))
  } else if(areaVarName == "admin2") {
    tmp<-subset(svydesign, (admin2==area))
  } else {
    stop("unrecognized areaVarName")
  }
  
  tt2 <- tryCatch(glmob<-svyglm(y.x~1,
                                design=tmp,family=quasibinomial, maxit=50), 
                  error=function(e) e, warning=function(w) w)
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob, logisticApproximation=logisticApproximation)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob, logisticApproximation=logisticApproximation)
    res = c(res, 0)
    return(res)
  }
}

region.time.HTDat<-function(dataobj, svydesign, area, nationalEstimate=FALSE, areaVarName="admin1", logisticApproximation=FALSE){
  # first check if observations are all 0 or 1
  if(sum(dataobj$y[dataobj[[areaVarName]] == area]) == 0) {
    return(c(0, -Inf, Inf, -Inf, Inf, 0))
  } else if(sum(dataobj$y[dataobj[[areaVarName]] == area]) == sum(dataobj$n[dataobj[[areaVarName]] == area])) {
    return(c(1, -Inf, Inf, Inf, Inf, 0))
  }
  
  if(!nationalEstimate) {
    
    if(areaVarName == "admin1") {
      tmp<-subset(svydesign, (admin1==area))
    } else if(areaVarName == "admin2") {
      tmp<-subset(svydesign, (admin2==area))
    } else {
      stop("unrecognized areaVarName")
    }
    
    tt2 <- tryCatch(glmob<-svyglm(y~1,
                                  design=tmp,family=quasibinomial, maxit=50), 
                    error=function(e) e, warning=function(w) w)
  } else {
    thisUrban = area == 1
    tmp<-subset(svydesign, (urban==thisUrban))
    tt2 <- tryCatch(glmob<-svyglm(y~1,
                                  design=tmp,family=quasibinomial, maxit=50), 
                    error=function(e) e, warning=function(w) w)
  }
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob, logisticApproximation=logisticApproximation)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob, logisticApproximation=logisticApproximation)
    res = c(res, 0)
    return(res)
  }
}

region.time.HTStratumDat<-function(dataobj, svydesign, area, urb, areaVarName="admin1", logisticApproximation=FALSE){
  # first check if observations are all 0 or 1
  if(sum(dataobj$y[dataobj[[areaVarName]] == area & dataobj$urban == urb]) == 0) {
    return(c(0, -Inf, Inf, -Inf, Inf, 0))
  } else if(sum(dataobj$y[dataobj[[areaVarName]] == area & dataobj$urban == urb]) == sum(dataobj$n[dataobj[[areaVarName]] == area & dataobj$urban == urb])) {
    return(c(1, -Inf, Inf, Inf, Inf, 0))
  }
  
  # make sure we don't have empirical proportion of zero or one
  thisDat = dataobj[dataobj[[areaVarName]] == area & dataobj$urban == urb,]
  if(sum(thisDat$y) == 0 || sum(thisDat$y) == sum(thisDat$n)) {
    return(data.frame(est=NA, lower=NA, upper=NA, logit.est=NA, var.est=NA, converge=1))
  }
  
  if(areaVarName == "admin1") {
    tmp<-subset(svydesign, (admin1==area) & (urban == urb))
  } else if(areaVarName == "admin2") {
    tmp<-subset(svydesign, (admin2==area) & (urban == urb))
  } else {
    stop("unrecognized areaVarName")
  }
  
  tt2 <- tryCatch(glmob<-svyglm(y~1,
                                design=tmp,family=quasibinomial, maxit=50), 
                  error=function(e) e, warning=function(w) w)
  
  if(is(tt2, "warning")){
    if(grepl("agegroups", tt2)){
      res <- get.est(glmob, logisticApproximation=logisticApproximation)
      res = c(res, 2)
    } else {
      res = c(rep(NA, 5), 3)
    }
    return(res)
  }
  if(is(tt2,"error")){
    res = c(rep(NA, 5), 1)
    return(res)
  } else {
    res <- get.est(glmob, logisticApproximation=logisticApproximation)
    res = c(res, 0)
    return(res)
  }
}


# # get direct estimates
# modDirect <- function(dat_obj, stratVar, useSamplingWeights=TRUE){
#   
#   options(survey.lonely.psu="adjust")
#   
#   # --- setting up a place to store results --- #
#   regions <- sort(unique(dat_obj$admin1))
#   regions_num  <- 1:length(regions)
#   
#   results<-data.frame(admin1=rep(regions,each=1))
#   results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
#   results$converge <- NA
#   
#   if(useSamplingWeights){
#     dat_obj$wt <- dat_obj$samplingWeight
#   } else {
#     dat_obj$wt <- NULL
#   }
#   
#   if(is.null(stratVar)){
#     # --- setting up the design object --- #
#     ## NOTE: -the v001 denote
#     ##        one stage cluster design (v001 is cluster)
#     ##       -This call below specifies our survey design
#     ##        nest = T argument nests clusters within strata
#     my.svydesign <- svydesign(id= ~v001,
#                               strata =NULL,
#                               weights=NULL, data=dat_obj)
#   } else {
#     ## not in all surveys does v022 contain the correct sampling strata
#     ## Thus, the correct vector has to be provided externally
#     dat_obj$strat <- stratVar
#     
#     # --- setting up the design object --- #
#     ## NOTE: -the v001 denote
#     ##        one stage cluster design (v001 is cluster)
#     ##       -This call below specifies our survey design
#     ##        nest = T argument nests clusters within strata
#     my.svydesign <- svydesign(id= ~v001,
#                               strata=~strat, nest=T, 
#                               weights=~wt, data=dat_obj)
#   }
#   
#   for(i in 1:nrow(results)){
#     results[i, 2:7] <- region.time.HT(dataobj=dat_obj, svydesign=my.svydesign, 
#                                       area=results$admin1[i])
#   }
#   return(results)
# }

# get direct estimates 
modDirectDat <- function(dat_obj, stratVar, useSamplingWeights=TRUE, nationalEstimate=FALSE, 
                            getContrast=nationalEstimate){
  
  options(survey.lonely.psu="adjust")
  
  # --- setting up a place to store results --- #
  regions <- sort(unique(dat_obj$admin1))
  regions_num  <- 1:length(regions)
  
  if(!nationalEstimate) {
    results<-data.frame(admin1=rep(regions,each=1))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  }
  else {
    results<-data.frame(urban=c(TRUE, FALSE))
    results$var.est<-results$logit.est<-results$upper<-results$lower<-results$est<-NA
    results$converge <- NA
  }
  
  if(useSamplingWeights){
    dat_obj$wt <- dat_obj$samplingWeight
  } else {
    dat_obj$wt <- NULL
  }
  
  if(is.null(stratVar)){
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata =NULL,
                              weights=NULL, data=dat_obj)
  } else {
    ## not in all surveys does v022 contain the correct sampling strata
    ## Thus, the correct vector has to be provided externally
    dat_obj$strat <- stratVar
    
    # --- setting up the design object --- #
    ## NOTE: -the v001 denote
    ##        one stage cluster design (v001 is cluster)
    ##       -This call below specifies our survey design
    ##        nest = T argument nests clusters within strata
    my.svydesign <- svydesign(id= ~v001,
                              strata=~strat, nest=T, 
                              weights=~wt, data=dat_obj)
  }
  
  for(i in 1:nrow(results)){
    if(!nationalEstimate) {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                           area=results$admin1[i], nationalEstimate=nationalEstimate)
    }
    else {
      results[i, 2:7] <- region.time.HTDat(dataobj=dat_obj, svydesign=my.svydesign, 
                                           area=i, nationalEstimate=nationalEstimate)
    }
  }
  
  if(getContrast) {
    # out = svyby(~y, by = ~urban, design = svydesign, svymean)
    glmob<-svyglm(y~urban,
                  design=my.svydesign,family=quasibinomial, maxit=50)
    
    # get contrast mean and variance
    est = glmob$coefficients[2]
    urbanVar = vcov(glmob)[2,2]
    
    # get confidence interval
    lower = est + qnorm(0.025, sd=sqrt(urbanVar))
    upper = est + qnorm(0.975, sd=sqrt(urbanVar))
    contrastStats = list(est=est, sd=sqrt(urbanVar), lower95=lower, upper95=upper)
    return(list(results=results, contrastStats=contrastStats))
  } else {
    return(results)
  }
  
}