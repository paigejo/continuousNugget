# Functions for simulating populations

# Use the easpc variable to construct a table suitable for input into the aggregation models. 
# output must have the following values: 
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
# Also, the rows of this table will be ordered by Area
makeDefaultEASPA = function(dataType=c("children", "women"), validationClusterI=NULL, 
                            useClustersAsEAs=!is.null(validationClusterI), usePixelUrban=TRUE) {
  dataType = match.arg(dataType)
  
  if(!useClustersAsEAs) {
    out = merge(easpc, adjustPopulationPerCountyTable(dataType))
  } else {
    if(is.null(validationClusterI)) {
      # This case will likely not be used, and corresponds with leave one county out
      out = merge(easpcMort, poppcMort)
    } else {
      ## construct faux population based on left out observations
      thisMort = mort[validationClusterI,]
      
      # change urbanicity to take the pixel level values rather than the cluster level values if necessary:
      test = getPixelIndex(cbind(thisMort$east, thisMort$north), popGrid, thisMort$admin1)
      temp = thisMort
      if(usePixelUrban) {
        temp$urban = popGrid$urban[test]
      }
      
      # faux poppc
      out = aggregate(temp$n, by=list(admin1=temp$admin1, urban=temp$urban), FUN=sum, drop=FALSE)
      out[is.na(out[,3]), 3] = 0
      # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
      out2 = cbind(out, rural=0)[48:94,]
      out2[, 4] = out$x[1:47]
      unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
      poppcMort = poppc
      poppcMort[unsortedToSorted, 2:3] = out2[,3:4]
      poppcMort$popTotal = poppcMort$popUrb + poppcMort$popRur
      poppcMort$pctTotal = 100 * poppcMort$popTotal * (1 / sum(poppcMort$popTotal))
      poppcMort$pctUrb = 100 * poppcMort$popUrb / poppcMort$popTotal
      
      # faux easpc
      out = aggregate(temp$n, by=list(admin1=temp$admin1, urban=temp$urban), FUN=length, drop=FALSE)
      out[is.na(out[,3]), 3] = 0
      # urbanToRuralI = c(1:27, 29, 31:47) # skip mombasa and nairobi
      out2 = cbind(out, rural=0)[48:94,]
      out2[, 4] = out$x[1:47]
      unsortedToSorted = sort(poppc$County, index.return=TRUE)$ix
      easpcMort = clustpc
      names(easpcMort) = names(easpc)
      easpcMort[unsortedToSorted, 2:3] = out2[,3:4]
      easpcMort$EATotal = easpcMort$EAUrb + easpcMort$EARur
      
      easpcMort$HHUrb = 25 * easpcMort$EAUrb
      easpcMort$HHRur = 25 * easpcMort$EARur
      easpcMort$HHTotal = easpcMort$HHUrb + easpcMort$HHRur
      
      # combine results
      out = merge(easpcMort, poppcMort)
    }
  }
  
  names(out)[1] = "area"
  
  # sort by area
  sortI = sort(out$area, index.return=TRUE)$ix
  out = out[sortI,]
  
  out
}

# Use the popGrid variable to construct a matrix suitable for input into the aggregation models. 
# population matrix must have the following values: 
# lon: longitude
# lat: latitude
# east: easting (km)
# north: northing (km)
# pop: proportional to population density for each grid cell
# area: an id or area name in which the grid cell corresponding to each row resides
# urban: whether the grid cell is urban or rural
# constituency: the sub-area
# province: the super-area
makeDefaultPopMat = function(getAdjusted=FALSE) {
  if(getAdjusted) {
    out = popGridAdjusted
  } else {
    out = popGrid
  }
  # out$province = out$region
  
  out
}

#' Simulate populations and areal prevalences
#' 
#' Given a spatial risk model, simulate populations and population prevalences at the 
#' enumeration area level (represented as points), and aggregate to the pixel and 
#' administrative areal level.
#' 
#' @param nsim Number of simulations
#' @param easpa data.frame of enumeration area, households, and target population per area stratified by urban/rural with variables:
#' \describe{
#'   \item{area}{name of area}
#'   \item{EAUrb}{number of urban enumeration areas in the area}
#'   \item{EARur}{number of rural enumeration areas in the area}
#'   \item{EATotal}{total number of enumeration areas in the area}
#'   \item{HHUrb}{number of urban households in the area}
#'   \item{HHRur}{number of rural households in the area}
#'   \item{HHTotal}{total number of households in the area}
#'   \item{popUrb}{total urban (target) population of area}
#'   \item{popRur}{total rural (target) population of area}
#'   \item{popTotal}{total (general) population of area}
#' }
#' @param popMat Pixellated grid data frame with variables `lon`, `lat`, `pop`, `area`, `subareas` (if subareaLevel is TRUE), `urban` (if stratifyByUrban is TRUE), `east`, and `north`
#' @param targetPopMat Same as popMat, but `pop` variable gives target rather than general population
#' @param poppsub data.frame of population per subarea separated by 
#' urban/rural using for population density grid normalization or urbanicity 
#' classification. Often based on extra fine scale population density grid. Has variables:
#' @param spdeMesh Triangular mesh for the SPDE
#' @param margVar Marginal variance of the spatial process, excluding cluster effects. 
#'          If 0, no spatial component is included
#' @param effRange Effective spatial range for the SPDE model
#' @param beta0 Intercept of logit model for risk
#' @param gamma Effect of urban on logit scale for logit model for risk
#' @param sigmaEpsilon Standard deviation on the logit scale for iid Gaussian EA level random effects in the risk model
#' @param seed Random number generator seed
#' @param inla.seed Seed input to inla.qsample. 0L sets seed intelligently, 
#'            > 0 sets a specific seed, < 0 keeps existing RNG
#' @param nHHSampled Number of households sampled per enumeration area. Default is 25 to match DHS surveys
#' @param stratifyByUrban Whether or not to stratify simulations by urban/rural classification
#' @param subareaLevel Whether or not to aggregate the population by subarea
#' @param doFineScaleRisk Whether or not to calculate the fine scale risk at each aggregation level in addition to the prevalence
#' @param doSmoothRisk Whether or not to calculate the smooth risk at each aggregation level in addition to the prevalence
#' @param doSmoothRiskLogisticApprox Whether to use logistic approximation when calculating smooth risk. See 
#' \code{\link{logitNormMean}} for details.
#' @param min1PerSubarea If TRUE, ensures there is at least 1 EA per subarea. If subareas are particularly unlikely to 
#' have enumeration areas since they have a very low proportion of the population in an area, then setting this to TRUE may be 
#' computationally intensive.
#' @param logitRiskDraws nIntegrationPoints x nsim dimension matrix of draws from the pixel leve risk field on logit scale, leaving out 
#' potential nugget/cluster/EA level effects.
#' @param sigmaEpsilonDraws nsim length vector of draws of cluster effect logit scale SD (joint draws with logitRiskDraws)
#' @param validationPixelI CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of pixels for which we want to simulate populations (used for pixel level validation)
#' @param validationClusterI CURRENTLY FOR TESTING PURPOSES ONLY a set of indices of cluster for which we want to simulate populations (used for cluster level validation)
#' @param clustersPerPixel CURRENTLY FOR TESTING PURPOSES ONLY Used for pixel level validation. Fixes the number of EAs per pixel.
#' @param returnEAinfo If TRUE, returns information on every individual EA (BAU) for each simulated population
#' @param epsc nEAs x nsim matrix of simulated EA (BAU) level iid effects representing fine scale variation in 
#'       risk. If NULL, they are simulated as iid Gaussian on a logit scale with 
#'       SD given by sigmaEpsilonDraws
#' list(pixelPop=outPixelLevel, subareaPop=outSubareaLevel, areaPop=outAreaLevel, logitRiskDraws=logitRiskDraws)
#' @return The simulated population aggregated to the enumeration area, 
#' pixel, subarea (generally Admin2), and area (generally Admin1) levels. Output includes:
#' \item{pixelPop}{A list of pixel level population aggregates}
#' \item{subareaPop}{A list of `subarea` level population aggregates}
#' \item{areaPop}{A list of `area` level population aggregates}
#' Each of these contains population numerator and denominator as well as prevalence and risk 
#' information aggregated to the appropriate level.
#' 
#' @details For population simulation and aggregation, we consider three models: smooth  
#' risk, fine scale risk, and the fine scale prevalence. All will be described in detail 
#' in a paper in preparation. In the smooth risk model, pixel level risks are integrated 
#' with respect to target population density when producing areal estimates on a prespecified 
#' set of integration points. The target population may be, for example, neonatals rather 
#' than the general population. In the fine scale models, enumeration areas (EAs) are simulated as 
#' point locations and iid random effects in the EA level risk are allowed. EAs and populations are dispersed conditional on the (possibly 
#' approximately) known number of EAs, households, and target population at a particular 
#' areal level (these we call `areas`) using multilevel multinomial sampling, first 
#' sampling the EAs, then distributing households among the EAs, then the target population 
#' among the households. Any areal level below the `areas` we call `subareas`. For instance, 
#' the `areas` might be Admin-1 if that is the smallest level at which the number of EAs, 
#' households, and people is known, and the `subareas` might be Admin-2. The multilevel 
#' multinomial sampling may be stratified by urban/rural within the areas if the number of 
#' EAs, households, and people is also approximately known at that level.
#' 
#' Within each EA we assume a fixed probability of an event occurring, which is the fine scale `risk`. 
#' The fine scale `prevalence` is the empirical proportion of events within that EA. We assume EA 
#' level logit scale iid N(0, sigmaEpsilon^2) random effects in the risk model. When averaged 
#' with equal weights over all EAs in an areal unit, this forms the fine scale risk. When 
#' instead the population numerators and denominators are aggregated, and are used to 
#' calculate the empirical proportion of events occurring in an areal unit, the resulting 
#' quantity is the fine scale prevalence in that areal unit.
#' 
#' Note that these functions can be used for either simulating populations for simulation 
#' studies, or for generating predictions accounting for uncertainty in EA locations 
#' and fine scale variation occuring at the EA level due to EA level iid random effects. 
#' Required, however, is a seperately fit EA level spatial risk model 
#' and information on the spatial population density and the population frame.
#' 
#' @author John Paige
#' @references In preparation
#' @seealso \code{\link{simPopCustom}}, \code{\link{makePopIntegrationTab}}, \code{\link{adjustPopMat}}, \code{\link{simSPDE}}.
#' @examples 
#' \dontrun{
#' ## In this script we will create 5km resolution pixellated grid over Kenya, 
#' ## and generate tables of estimated (both target and general) population 
#' ## totals at the area (e.g. Admin-1) and subarea (e.g. Admin-2) levels. Then 
#' ## we will use that to simulate populations of 
#' 
#' # download Kenya GADM shapefiles from SUMMERdata github repository
#' githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
#'                     "kenyaMaps.rda?raw=true")
#' tempDirectory = "~/"
#' mapsFilename = paste0(tempDirectory, "/kenyaMaps.rda")
#' if(!file.exists(mapsFilename)) {
#'   download.file(githubURL,mapsFilename)
#' }
#' 
#' # load it in
#' out = load(mapsFilename)
#' out
#' adm1@data$NAME_1 = as.character(adm1@data$NAME_1)
#' adm1@data$NAME_1[adm1@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
#' adm1@data$NAME_1[adm1@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
#' adm2@data$NAME_1 = as.character(adm2@data$NAME_1)
#' adm2@data$NAME_1[adm2@data$NAME_1 == "Trans Nzoia"] = "Trans-Nzoia"
#' adm2@data$NAME_1[adm2@data$NAME_1 == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
#' 
#' # some Admin-2 areas have the same name
#' adm2@data$NAME_2 = as.character(adm2@data$NAME_2)
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Bungoma") & 
#'                    (adm2@data$NAME_2 == "Lugari")] = "Lugari, Bungoma"
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Kakamega") & 
#'                    (adm2@data$NAME_2 == "Lugari")] = "Lugari, Kakamega"
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Meru") & 
#'                    (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Meru"
#' adm2@data$NAME_2[(adm2@data$NAME_1 == "Tharaka-Nithi") & 
#'                    (adm2@data$NAME_2 == "Igembe South")] = "Igembe South, Tharaka-Nithi"
#' 
#' # The spatial area of unknown 8 is so small, it causes problems unless its removed or 
#' # unioned with another subarea. Union it with neighboring Kakeguria:
#' newadm2 = adm2
#' unknown8I = which(newadm2$NAME_2 == "unknown 8")
#' newadm2$NAME_2[newadm2$NAME_2 %in% c("unknown 8", "Kapenguria")] <- 
#'   "Kapenguria + unknown 8"
#' admin2.IDs <- newadm2$NAME_2
#' 
#' library(maptools)
#' temp <- unionSpatialPolygons(newadm2, admin2.IDs)
#' tempData = newadm2@data[-unknown8I,]
#' tempData = tempData[order(tempData$NAME_2),]
#' newadm2 <- SpatialPolygonsDataFrame(temp, tempData, match.ID = F)
#' adm2 = newadm2
#' 
#' # download 2014 Kenya population density and associated TIF file
#' githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
#'                     "Kenya2014Pop/pop.rda?raw=true")
#' popFilename = paste0(tempDirectory, "/pop.rda")
#' if(!file.exists(popFilename)) {
#'   download.file(githubURL,popFilename)
#' }
#' 
#' githubURL <- paste0("https://github.com/paigejo/SUMMERdata/blob/main/data/", 
#'                     "Kenya2014Pop/worldpop_total_1y_2014_00_00.tif?raw=true")
#' popTIFFilename = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")
#' if(!file.exists(popTIFFilename)) {
#'   download.file(githubURL,popTIFFilename)
#' }
#' 
#' # load it in
#' require(raster)
#' out = load(popFilename)
#' out
#' 
#' # make sure this is correct for re-projections
#' pop@file@name = paste0(tempDirectory, "/worldpop_total_1y_2014_00_00.tif")
#' 
#' eastLim = c(-110.6405, 832.4544)
#' northLim = c(-555.1739, 608.7130)
#' 
#' ## Construct poppsubKenya, a table of urban/rural general population totals 
#' ## in each subarea. Technically, this is not necessary since we can load in 
#' ## poppsubKenya via data(kenyaPopulationData). First, we will need to calculate 
#' ## the areas in km^2 of the areas and subareas
#' 
#' library(rgdal)
#' library(sp)
#' 
#' # use Lambert equal area projection of areas (Admin-1) and subareas (Admin-2)
#' midLon = mean(adm1@bbox[1,])
#' midLat = mean(adm1@bbox[2,])
#' p4s = paste0("+proj=laea +x_0=0 +y_0=0 +lon_0=", midLon, 
#'              " +lat_0=", midLat, " +units=km")
#' 
#' library(rgdal)
#' 
#' adm1proj <- spTransform(adm1, CRS(p4s))
#' adm2proj <- spTransform(adm2, CRS(p4s))
#' 
#' # now calculate spatial area in km^2
#' library(rgeos)
#' admin1Areas = gArea(adm1proj, TRUE)
#' admin2Areas = gArea(adm2proj, TRUE)
#' areapaKenya = data.frame(area=adm1proj@data$NAME_1, spatialArea=admin1Areas)
#' areapsubKenya = data.frame(area=adm2proj@data$NAME_1, subarea=adm2proj@data$NAME_2, 
#'                            spatialArea=admin2Areas)
#' 
#' # Calculate general population totals at the subarea (Admin-2) x urban/rural 
#' # level and using 1km resolution population grid. Assign urbanicity by 
#' # thresholding population density based on estimated proportion population 
#' # urban/rural, making sure total area (Admin-1) urban/rural populations in 
#' # each area matches poppaKenya.
#' require(fields)
#' # NOTE: the following function will typically take ~20 minutes. Can speed up 
#' #       by setting kmRes to be higher, but we recommend fine resolution for 
#' #       this step, since it only needs to be done once.
#' system.time(poppsubKenya <- getPoppsub(
#'   kmRes=1, pop=pop, domainPoly=kenyaPoly,
#'   eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
#'   poppa = poppaKenya, areapa=areapaKenya, areapsub=areapsubKenya, 
#'   areaMapDat=adm1, subareaMapDat=adm2, 
#'   areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
#' 
#' # Now generate a general population integration table at 5km resolution, 
#' # based on subarea (Admin-2) x urban/rural population totals. This takes 
#' # ~1 minute
#' system.time(popMatKenya <- makePopIntegrationTab(
#'   kmRes=5, pop=pop, domainPoly=kenyaPoly,
#'   eastLim=eastLim, northLim=northLim, mapProjection=projKenya,
#'   poppa = poppaKenya, poppsub=poppsubKenya, 
#'   areaMapDat = adm1, subareaMapDat = adm2,
#'   areaNameVar = "NAME_1", subareaNameVar="NAME_2"))
#' 
#' ## Adjust popMat to be target (neonatal) rather than general population 
#' ## density. First create the target population frame
#' ## (these numbers are based on IPUMS microcensus data)
#' mothersPerHouseholdUrb = 0.3497151
#' childrenPerMotherUrb = 1.295917
#' mothersPerHouseholdRur = 0.4787696
#' childrenPerMotherRur = 1.455222
#' targetPopPerStratumUrban = easpaKenya$HHUrb * mothersPerHouseholdUrb * 
#'   childrenPerMotherUrb
#' targetPopPerStratumRural = easpaKenya$HHRur * mothersPerHouseholdRur * 
#'   childrenPerMotherRur
#' easpaKenyaNeonatal = easpaKenya
#' easpaKenyaNeonatal$popUrb = targetPopPerStratumUrban
#' easpaKenyaNeonatal$popRur = targetPopPerStratumRural
#' easpaKenyaNeonatal$popTotal = easpaKenyaNeonatal$popUrb + 
#'   easpaKenyaNeonatal$popRur
#' easpaKenyaNeonatal$pctUrb = 100 * easpaKenyaNeonatal$popUrb / 
#'   easpaKenyaNeonatal$popTotal
#' easpaKenyaNeonatal$pctTotal = 
#'   100 * easpaKenyaNeonatal$popTotal / sum(easpaKenyaNeonatal$popTotal)
#' 
#' # Generate the target population density by scaling the current 
#' # population density grid at the Admin1 x urban/rural level
#' popMatKenyaNeonatal = adjustPopMat(popMatKenya, easpaKenyaNeonatal)
#' 
#' # Generate neonatal population table from the neonatal population integration 
#' # matrix. This is technically not necessary for population simulation purposes, 
#' # but is here for illustrative purposes
#' poppsubKenyaNeonatal = poppRegionFromPopMat(popMatKenyaNeonatal, 
#'                                             popMatKenyaNeonatal$subarea)
#' poppsubKenyaNeonatal = 
#'   cbind(subarea=poppsubKenyaNeonatal$region, 
#'         area=adm2@data$NAME_1[match(poppsubKenyaNeonatal$region, adm2@data$NAME_2)], 
#'         poppsubKenyaNeonatal[,-1])
#' print(head(poppsubKenyaNeonatal))
#' 
#' ## Now we're ready to simulate neonatal populations along with neonatal 
#' ## mortality risks and prevalences
#' 
#' # use the following model to simulate the neonatal population based roughly 
#' # on Paige et al. (2020) neonatal mortality modeling for Kenya.
#' beta0=-2.9 # intercept
#' gamma=-1 # urban effect
#' rho=(1/3)^2 # spatial variance
#' effRange = 400 # effective spatial range in km
#' sigmaEpsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # Run a simulation! This produces multiple dense nEA x nsim and nPixel x nsim 
#' # matrices. In the future sparse matrices and chunk by chunk computations 
#' # may be incorporated.
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spdeMesh=kenyaMesh, 
#'                     margVar=rho, sigmaEpsilon=sigmaEpsilon, 
#'                     gamma=gamma, effRange=effRange, beta0=beta0, 
#'                     seed=12, inla.seed=12, nHHSampled=25, 
#'                     stratifyByUrban=TRUE, subareaLevel=TRUE, 
#'                     doFineScaleRisk=TRUE, doSmoothRisk=TRUE, 
#'                     min1PerSubarea=TRUE)
#' 
#' # get average absolute percent error relative to fine scale prevalence at Admin-2 level
#' tempDat = simPop$subareaPop$aggregationResults[,c("region", "pFineScalePrevalence", 
#'                                                   "pFineScaleRisk", "pSmoothRisk")]
#' 100*mean(abs(tempDat$pFineScalePrevalence - tempDat$pFineScaleRisk) / 
#'            tempDat$pFineScalePrevalence)
#' 100*mean(abs(tempDat$pFineScalePrevalence - tempDat$pSmoothRisk) / 
#'            tempDat$pFineScalePrevalence)
#' 100*mean(abs(tempDat$pFineScaleRisk - tempDat$pSmoothRisk) / 
#'            tempDat$pFineScalePrevalence)
#' 
#' # verify number of EAs per area and subarea
#' cbind(aggregate(simPop$eaPop$eaSamples[,1], by=list(area=popMatKenya$area), FUN=sum), 
#'       trueNumEAs=easpaKenya$EATotal[order(easpaKenya$area)])
#' aggregate(simPop$eaPop$eaSamples[,1], by=list(area=popMatKenya$subarea), FUN=sum)
#' 
#' ## plot simulated population
#' # directory for plotting 
#' # (mapPlot takes longer when it doesn't save to a file)
#' tempDirectory = "~/"
#' 
#' # pixel level
#' zlim = c(0, quantile(probs=.995, c(simPop$pixelPop$pFineScalePrevalence, 
#'                                    simPop$pixelPop$pFineScaleRisk, 
#'                                    simPop$pixelPop$pSmoothRisk), na.rm=TRUE))
#' pdf(file=paste0(tempDirectory, "simPopSPDEPixel.pdf"), width=8, height=8)
#' par(mfrow=c(2,2))
#' plot(adm2, asp=1)
#' points(simPop$eaPop$eaDatList[[1]]$lon, simPop$eaPop$eaDatList[[1]]$lat, pch=".", col="blue")
#' plot(adm2, asp=1)
#' quilt.plot(popMatKenya$lon, popMatKenya$lat, simPop$pixelPop$pFineScalePrevalence, 
#'            zlim=zlim, add=TRUE, FUN=function(x) {mean(x, na.rm=TRUE)})
#' plot(adm2, asp=1)
#' quilt.plot(popMatKenya$lon, popMatKenya$lat, simPop$pixelPop$pFineScaleRisk, 
#'            zlim=zlim, add=TRUE, FUN=function(x) {mean(x, na.rm=TRUE)})
#' quilt.plot(popMatKenya$lon, popMatKenya$lat, simPop$pixelPop$pSmoothRisk, 
#'            zlim=zlim, FUN=function(x) {mean(x, na.rm=TRUE)}, asp=1)
#' plot(adm2, add=TRUE)
#' dev.off()
#' 
#' range(simPop$eaPop$eaDatList[[1]]$N)
#' 
#' # areal (Admin-1) level: these results should look essentially identical
#' 
#' tempDat = simPop$areaPop$aggregationResults[,c("region", "pFineScalePrevalence", 
#'                                                "pFineScaleRisk", "pSmoothRisk")]
#' pdf(file=paste0(tempDirectory, "simPopSPDEAdmin-1.pdf"), width=7, height=7)
#' mapPlot(tempDat, 
#'         variables=c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk"), 
#'         geo=adm1, by.geo="NAME_1", by.data="region", is.long=FALSE)
#' dev.off()
#' 
#' # subareal (Admin-2) level: these results should look subtley different 
#' # depending on the type of prevalence/risk considered
#' tempDat = simPop$subareaPop$aggregationResults[,c("region", "pFineScalePrevalence", 
#'                                                   "pFineScaleRisk", "pSmoothRisk")]
#' pdf(file=paste0(tempDirectory, "simPopSPDEAdmin-2.pdf"), width=7, height=7)
#' mapPlot(tempDat, 
#'         variables=c("pFineScalePrevalence", "pFineScaleRisk", "pSmoothRisk"), 
#'         geo=adm2, by.geo="NAME_2", by.data="region", is.long=FALSE)
#' dev.off()
#' }
#' @name simPop
NULL

#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses SPDE model for generating spatial risk and can include iid cluster 
#' level effect.
#' 
#' @export
simPopSPDE = function(nsim=1, easpa, popMat, targetPopMat, poppsub, spdeMesh, 
                      margVar=0.243, sigmaEpsilon=sqrt(0.463), 
                      gamma=0.009, effRange=406.51, beta0=-3.922, 
                      seed=NULL, inla.seed=-1L, nHHSampled=25, 
                      stratifyByUrban=TRUE, subareaLevel=TRUE, gridLevel=FALSE, 
                      doFineScaleRisk=FALSE, doSmoothRisk=FALSE, 
                      doGriddedRisk=FALSE, 
                      doSmoothRiskLogisticApprox=FALSE, 
                      min1PerSubarea=TRUE, 
                      fixPopPerEA=NULL, fixHHPerEA=NULL, fixPopPerHH=NULL) {
  if(!is.null(seed))  {
    set.seed(seed)
    
    if(inla.seed < 0) {
      stop("seed specified, but not inla.seed. Set inla.seed to a positive integer to ensure reproducibility")
    }
  }
  
  if(nsim > 1) {
    warning("nsim > 1. eaDat will only be generated for first simulation")
  }
  
  totalEAs = sum(easpa$EATotal)
  totalHouseholds = sum(easpa$HHTotal)
  
  ### generate Binomial probabilities from transformed logit scale GP
  # generate SPDE simulations
  pixelCoords = cbind(popMat$east, popMat$north)
  
  print("Using SPDE model to simulate EA and pixel level risk")
  if(margVar != 0) {
    SPDEArgs = list(coords=pixelCoords, nsim=nsim, margVar=margVar, effRange=effRange, 
                    mesh=spdeMesh, inla.seed=inla.seed)
    simVals = do.call("simSPDE", SPDEArgs)
  } else {
    simVals = matrix(rep(0, nrow(pixelCoords)), ncol=1)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*popMat$urban, "+")
  
  # simulate nugget/cluster effect
  epsc = matrix(stats::rnorm(totalEAs*nsim, sd=sigmaEpsilon), ncol=nsim)
  
  # transform back to original scale for the pixel level probabilities
  probsNoNug = expit(simVals)
  
  # simulate the enumeration areas
  logitRiskDraws = simVals
  sigmaEpsilonDraws = rep(sigmaEpsilon, nsim)
  
  out = simPopCustom(logitRiskDraws=logitRiskDraws, sigmaEpsilonDraws=sigmaEpsilonDraws, easpa=easpa, 
                     popMat=popMat, targetPopMat=targetPopMat, 
                     stratifyByUrban=stratifyByUrban, doSmoothRisk=doSmoothRisk, 
                     doGriddedRisk=doGriddedRisk, 
                     doSmoothRiskLogisticApprox=doSmoothRiskLogisticApprox, 
                     doFineScaleRisk=doFineScaleRisk, poppsub=poppsub, 
                     subareaLevel=subareaLevel, gridLevel=gridLevel, 
                     min1PerSubarea=min1PerSubarea, returnEAinfo=TRUE, epsc=epsc, 
                     fixPopPerEA=fixPopPerEA, fixHHPerEA=fixHHPerEA, fixPopPerHH=fixPopPerHH)
  eaPop = list(eaDatList=out$eaDatList, eaSamples=out$eaSamples)
  out$eaDatList = NULL
  out$eaSamples = NULL
  outPixelLevel = out$pixelPop
  outSubareaLevel = out$subareaPop
  outAreaLevel = out$areaPop
  
  list(eaPop=eaPop, pixelPop=outPixelLevel, subareaPop=outSubareaLevel, areaPop=outAreaLevel, logitRiskDraws=logitRiskDraws)
}

#' Aggregate populations to the specified areal level
#' 
#' Takes simulated populations and aggregates 
#' them to the specified areal level. Also calculates the aggregated risk and prevalence.
#' 
#' @param pixelLevelPop pixel level population information that we want aggregate. In the same format as output from \code{\link{simPopCustom}}
#' @param eaSamples nIntegrationPoint x nsim matrix of the number of enumeration areas per pixel sampled in the input pixel level population
#' @param areas character vector of length nIntegrationPoints of area names over which we 
#'        want to aggregate. Can also be subareas
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param targetPopMat pixellated grid data frame with variables `lon`, `lat`, `pop` (target population), `area`, `subareas` (if subareaLevel is TRUE), `urban` (if stratifyByUrban is TRUE), `east`, and `north`
#' @param doFineScaleRisk whether or not to calculate the fine scale risk in addition to the prevalence. See details
#' @param doSmoothRisk Whether or not to calculate the smooth risk in addition to the prevalence. See details
#' @param easpa see \code{\link{simPopSPDE}}
#' @param areaLevelPop output of \code{\link{simPopCustom}} containing pixel level information 
#'               about the population of interest
#' @param areasFrom character vector of length equal to the number of areas from which 
#'            we would like to aggregate containing the unique names of the areas. 
#'            Can also be subareas, but these are smaller than the "to areas", and 
#'            each "from area" must be entirely contained in a single "to area"
#' @param areasTo character vector of length equal to the number of areas from which 
#'          we would like to aggregate containing the names of the areas containing 
#'          with each respective `from' area. Can also be a set of subareas, 
#'          but these are larger than the "from areas".
#' 
#' @author John Paige
#' @references In Preparation
#' @return A list containing elements `fineScalePrevalence` and `fineScaleRisk`. Each 
#' of these are in turn lists with aggregated prevalence and risk for the area of 
#' interest, containg the following elements, were paranethesis indicate the elements 
#' for the fineScaleRisk model rather than fineScalePrevalence:
#' \item{p}{Aggregated prevalence (risk), calculated as aggregate of Z divided by 
#' aggregate of N}
#' \item{Z}{Aggregated (expected) population numerator}
#' \item{N}{Aggregated (expected) population denominator}
#' \item{pUrban}{Aggregated prevalence (risk) in urban part of the area, calculated 
#' as aggregate of Z divided by aggregate of N}
#' \item{ZUrban}{Aggregated (expected) population numerator in urban part of the area}
#' \item{NUrban}{Aggregated (expected) population denominator in urban part of the area}
#' \item{pRural}{Aggregated prevalence (risk) in rural part of the area, calculated 
#' as aggregate of Z divided by aggregate of N}
#' \item{ZRural}{Aggregated (expected) population numerator in rural part of the area}
#' \item{NRural}{Aggregated (expected) population denominator in rural part of the area}
#' \item{A}{Aggregation matrix used to aggregate from pixel level to areal level}
#' \item{AUrban}{Aggregation matrix used to aggregate from pixel level to urban part of the areal level}
#' \item{ARural}{Aggregation matrix used to aggregate from pixel level to rural part of the areal level}
#' @seealso \code{\link{areaPopToArea}}
#' @name aggPop
#' @examples
#' \dontrun{
#' ##### Now we make a model for the risk. We will use an SPDE model with these 
#' ##### parameters for the linear predictor on the logist scale, which are chosen 
#' ##### to be of practical interest:
#' beta0=-2.9 # intercept
#' gamma=-1 # urban effect
#' rho=(1/3)^2 # spatial variance
#' effRange = 400 # effective spatial range in km
#' sigmaEpsilon=sqrt(1/2.5) # cluster (nugget) effect standard deviation
#' 
#' # simulate the population! Note that this produces multiple dense 
#' # nEA x nsim and nIntegrationPoint x nsim matrices. In the future 
#' # sparse matrices will and chunk by chunk computations may be incorporated.
#' simPop = simPopSPDE(nsim=1, easpa=easpaKenyaNeonatal, 
#'                     popMat=popMatKenya, targetPopMat=popMatKenyaNeonatal, 
#'                     poppsub=poppsubKenya, spdeMesh=kenyaMesh, 
#'                     margVar=rho, sigmaEpsilonSq=sigmaEpsilon^2, 
#'                     gamma=gamma, effRange=effRange, beta0=beta0, 
#'                     seed=123, inla.seed=12, nHHSampled=25, 
#'                     stratifyByUrban=TRUE, subareaLevel=TRUE, 
#'                     doFineScaleRisk=TRUE, 
#'                     min1PerSubarea=TRUE)
#' 
#' pixelPop = simPop$pixelPop
#' subareaPop = pixelPopToArea(pixelLevelPop=pixelPop, eaSamples=pixelPop$eaSamples, 
#'   areas=popMatKenya$subarea, stratifyByUrban=TRUE, 
#'   targetPopMat=popMatKenyaNeonatal, doFineScaleRisk=TRUE)
#' 
#' # get areas associated with each subarea for aggregation
#' tempAreasFrom = popMatKenya$subarea
#' tempAreasTo = popMatKenya$area
#' areasFrom = sort(unique(tempAreasFrom))
#' areasToI = match(areasFrom, tempAreasFrom)
#' areasTo = tempAreasTo[areasToI]
#' 
#' # do the aggregation from subareas to areas
#' outAreaLevel = areaPopToArea(areaLevelPop=subareaPop, 
#'   areasFrom=areasFrom, areasTo=areasTo, 
#'   stratifyByUrban=TRUE, doFineScaleRisk=TRUE)
#' }
NULL

#' @describeIn aggPop Aggregate from EA to areal level
#' @export
EAPopToArea = function(EALevelPop, areaMat, areaLevels, urbanMat=NULL, 
                       doFineScalePrevalence=!is.null(EALevelPop$pFineScalePrevalence), 
                       doFineScaleRisk=!is.null(EALevelPop$pFineScaleRisk)) {
  
  aggregationResults = list()
  
  # fine scale prevalence aggregation model
  if(doFineScalePrevalence) {
    nSamples = EALevelPop$NFineScalePrevalence
    zSamples = EALevelPop$ZFineScalePrevalence
    zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works
    out = aggPredsVariablePerArea(popNumerators=zSamples, popDenominators=nSamples, 
                                  areaMat=areaMat, urbanMat=urbanMat)
    aggregationResultsFineScalePrevalence = out
    names(aggregationResultsFineScalePrevalence)[-1] = paste(names(aggregationResultsFineScalePrevalence)[-1], "FineScalePrevalence", sep="")
    
    aggregationResults = c(aggregationResults, aggregationResultsFineScalePrevalence[-1])
  }
  
  # fine scale risk aggregation model
  if(doFineScaleRisk) {
    nSamplesFineScaleRisk = EALevelPop$NFineScaleRisk
    zSamplesFineScaleRisk = EALevelPop$ZFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPredsVariablePerArea(popNumerators=zSamples, popDenominators=nSamples, 
                                  areaMat=areaMat, urbanMat=urbanMat)
    aggregationResultsFineScaleRisk = out
    names(aggregationResultsFineScaleRisk)[-1] = paste(names(aggregationResultsFineScaleRisk)[-1], "FineScaleRisk", sep="")
    
    aggregationResults = c(aggregationResults, aggregationResultsFineScaleRisk[-1])
  }
  
  aggregationResults
}

#' @describeIn aggPop Aggregate from grid to areal level
#' @export
pixelPopToArea = function(pixelLevelPop, eaSamples, areas, stratifyByUrban=TRUE, targetPopMat=NULL, 
                          doFineScalePrevalence=!is.null(pixelLevelPop$pFineScalePrevalence), 
                          doFineScaleRisk=!is.null(pixelLevelPop$pFineScaleRisk), 
                          doSmoothRisk=!is.null(pixelLevelPop$pSmoothRisk), 
                          doGriddedRisk=!is.null(pixelLevelPop$pGriddedRisk)) {
  
  aggregationResults = list(region=sort(unique(areas)))
  aggregationMatrices = list()
  
  # fine scale prevalence aggregation model
  if(doFineScalePrevalence) {
    nSamples = pixelLevelPop$NFineScalePrevalence
    zSamples = pixelLevelPop$ZFineScalePrevalence
    zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works
    out = aggPixelPreds(Zg=zSamples, Ng=nSamples, areas=areas, targetPopMat=targetPopMat, 
                        useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
    aggregationResultsFineScalePrevalence = out$aggregationResults
    aggregationMatricesFineScalePrevalence = out$aggregationMatrices
    names(aggregationResultsFineScalePrevalence)[-1] = paste(names(aggregationResultsFineScalePrevalence)[-1], "FineScalePrevalence", sep="")
    names(aggregationMatricesFineScalePrevalence) = paste(names(aggregationMatricesFineScalePrevalence), "FineScalePrevalence", sep="")
    
    aggregationResults = c(aggregationResults, aggregationResultsFineScalePrevalence[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesFineScalePrevalence)
  }
  
  
  # fine scale risk aggregation model
  if(doFineScaleRisk) {
    nSamplesFineScaleRisk = pixelLevelPop$NFineScaleRisk
    zSamplesFineScaleRisk = pixelLevelPop$ZFineScaleRisk
    zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPreds(Zg=zSamplesFineScaleRisk, Ng=nSamplesFineScaleRisk, areas=areas, targetPopMat=targetPopMat, 
                        useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
    aggregationResultsFineScaleRisk = out$aggregationResults
    aggregationMatricesFineScaleRisk = out$aggregationMatrices
    names(aggregationResultsFineScaleRisk)[-1] = paste(names(aggregationResultsFineScaleRisk)[-1], "FineScaleRisk", sep="")
    names(aggregationMatricesFineScaleRisk) = paste(names(aggregationMatricesFineScaleRisk), "FineScaleRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsFineScaleRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsFineScaleRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesFineScaleRisk)
  }
  
  if(doSmoothRisk) {
    # NOTE: although useDensity is set to FALSE, that's only because the density is already 
    # directly incorporated into nSamplesSmoothRisk
    nSamplesSmoothRisk = pixelLevelPop$NSmoothRisk
    zSamplesSmoothRisk = pixelLevelPop$ZSmoothRisk
    zSamplesSmoothRisk[is.na(zSamplesSmoothRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPreds(Zg=zSamplesSmoothRisk, Ng=nSamplesSmoothRisk, areas=areas, targetPopMat=targetPopMat, 
                        useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
    aggregationResultsSmoothRisk = out$aggregationResults
    aggregationMatricesSmoothRisk = out$aggregationMatrices
    names(aggregationResultsSmoothRisk)[-1] = paste(names(aggregationResultsSmoothRisk)[-1], "SmoothRisk", sep="")
    names(aggregationMatricesSmoothRisk) = paste(names(aggregationMatricesSmoothRisk), "SmoothRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsSmoothRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsSmoothRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesSmoothRisk)
  }
  
  if(doGriddedRisk) {
    # IHME risk model
    nSamplesGriddedRisk = pixelLevelPop$NgriddedRisk
    zSamplesGriddedRisk = pixelLevelPop$ZgriddedRisk
    zSamplesGriddedRisk[is.na(zSamplesGriddedRisk)] = 0 # must set to zero temporarily so matrix multiplication works out
    out = aggPixelPreds(Zg=zSamplesGriddedRisk, Ng=nSamplesGriddedRisk, areas=areas, targetPopMat=targetPopMat, 
                        useDensity=FALSE, stratifyByUrban=stratifyByUrban, normalize=FALSE)
    aggregationResultsGriddedRisk = out$aggregationResults
    aggregationMatricesGriddedRisk = out$aggregationMatrices
    names(aggregationResultsGriddedRisk)[-1] = paste(names(aggregationResultsGriddedRisk)[-1], "GriddedRisk", sep="")
    names(aggregationMatricesGriddedRisk) = paste(names(aggregationMatricesGriddedRisk), "GriddedRisk", sep="")
    
    # aggregationResults = merge(aggregationResults, aggregationResultsGriddedRisk, by="region")
    aggregationResults = c(aggregationResults, aggregationResultsGriddedRisk[-1])
    aggregationMatrices = c(aggregationMatrices, aggregationMatricesGriddedRisk)
  }
  
  list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
}

#' Helper function of \code{\link{pixelPopToArea}}
#' 
#' Aggregates population from the 
#' pixel level to the level of the area of interest.
#' 
#' @param popNumerators nIntegrationPoint x nsim matrix of simulated response (population numerators) for each pixel and sample
#' @param popDenominators nIntegrationPoint x nsim matrix of simulated counts (population denominators) for each pixel and sample
#' @param areas Either a nIntegrationPoint length character vector of areas (or subareas) or a nIntegrationPoint row x 
#' nDraws column matrix whose columns are vectors of areas associated with each value being aggregated.
#' @param urban nIntegrationPoint length vector of indicators specifying whether or not pixels are urban or rural. If areas is a 
#' matrix, this should be a matrix with equal dimension
#' @param targetPopMat same as in \code{\link{simPopCustom}}
#' @param stratifyByUrban whether or not to stratify simulations by urban/rural classification
#' @param normalize if TRUE, pixel level aggregation weights within specified area are normalized to sum to 1. This produces an 
#' average of the values in popNumerators rather than a sum. In general, should only be set to TRUE for smooth integrals of risk, say over 
#' target population density (i.e. if popDenominators is set to the target population density and popNumerators is set to the risk).
aggPreds = function(popNumerators, popDenominators, areas, urban=targetPopMat$urban, targetPopMat=NULL, 
                    stratifyByUrban=TRUE, normalize=FALSE) {
  
  if(is.matrix(areas)) {
    if(normalize) {
      stop("normalize must be set to FALSE if areas is a matrix")
    }
    if(stratifyByUrban && !is.matrix(urban)) {
      stop("Either both area and urban must be matrices, or both must be vectors")
    } else if(!stratifyByUrban) {
      urban = NULL
    }
    return(aggPredsVariablePerArea(popNumerators=popNumerators, popDenominators=popDenominators, 
                                   areaMat=areas, urbanMat=urban))
  }
  predsUrban = urban
  predsArea = areas
  
  # set NAs and pixels without any sample size to 0
  popDenominators[is.na(popDenominators)] = 0
  popNumerators[popDenominators == 0] = 0
  
  # function to aggregate predictions to the given areal level. Use the 
  # following function to get numerical integration matrix for a given 
  # level of areal aggregation. Returned matrices have dimension 
  # length(unique(areaNames)) x length(areaNames)
  # areaNames: length(popNumerators) length vector of area names associated with each value 
  #            being aggregated
  # normalize: whether or not to normalize the rows of the matrices to sum to 1 or to instead 
  #            contain only binary values (or non-binary values based on the binary values if 
  #            urbanProportions is not NULL)
  getIntegrationMatrix = function(areaNames, urbanProportions=NULL, normalize=FALSE) {
    
    equalDensities = rep(1, nrow(popNumerators))
    densities = equalDensities
    
    uniqueNames = sort(unique(areaNames))
    getMatrixHelper = function(i, thisUrban=NULL, thisNormalize=normalize) {
      areaI = areaNames == uniqueNames[i]
      
      theseDensities = equalDensities
      
      # make sure we only include pixels in the given area and, if necessary, with the given urbanicity
      theseDensities[!areaI] = 0
      if(!is.null(thisUrban))
        theseDensities[predsUrban != thisUrban] = 0
      thisSum = sum(theseDensities)
      if(thisSum != 0 && thisNormalize)
        theseDensities * (1/thisSum)
      else if(thisSum == 0)
        rep(0, length(theseDensities))
      else
        theseDensities
    }
    
    if(!stratifyByUrban) {
      integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      
      integrationMatrix
    } else {
      integrationMatrixUrban = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=TRUE), ncol=length(uniqueNames)))
      integrationMatrixRural = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper, thisUrban=FALSE), ncol=length(uniqueNames)))
      if(!is.null(urbanProportions)) {
        integrationMatrix = sweep(integrationMatrixUrban, 1, urbanProportions, "*") + sweep(integrationMatrixRural, 1, 1-urbanProportions, "*")
      } else {
        integrationMatrix = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
      }
      
      rownames(integrationMatrix) = uniqueNames
      rownames(integrationMatrixUrban) = uniqueNames
      rownames(integrationMatrixRural) = uniqueNames
      
      list(integrationMatrix=integrationMatrix, 
           integrationMatrixUrban=integrationMatrixUrban, 
           integrationMatrixRural=integrationMatrixRural)
    }
  }
  
  # Use the following function to perform the
  # aggregations
  getIntegratedPredictions = function(areaNames) {
    # get numerical integration matrix
    A = getIntegrationMatrix(areaNames, normalize=normalize)
    
    # aggregate the prediction and denominator matrices (for whole areas and also urban/rural strata if necessary)
    if(!stratifyByUrban) {
      ZAggregated = A %*% popNumerators
      NAggregated = A %*% popDenominators
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      aggregationResults = list(p=pAggregated, Z=ZAggregated, N=NAggregated)
      aggregationMatrices = list(A=A, AUrban=NULL, ARural=NULL)
    } else {
      AUrban = A$integrationMatrixUrban
      ARural = A$integrationMatrixRural
      A = A$integrationMatrix
      
      # first aggregate the numerator. The denominator will depend on the aggregation method
      ZAggregated = A %*% popNumerators
      ZAggregatedUrban = AUrban %*% popNumerators
      ZAggregatedRural = ARural %*% popNumerators
      
      # we must also aggregate the denominator to calculate 
      # the aggregated empirical proportions
      NAggregated = A %*% popDenominators
      pAggregated = ZAggregated / NAggregated
      pAggregated[NAggregated == 0] = NA
      
      NAggregatedUrban = AUrban %*% popDenominators
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = NA
      
      NAggregatedRural = ARural %*% popDenominators
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = NA
      
      aggregationResults = list(region=sort(unique(areas)), p=pAggregated, Z=ZAggregated, N=NAggregated, 
                                pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                                pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
      aggregationMatrices = list(A=A, AUrban=AUrban, ARural=ARural)
    }
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  areaResults = getIntegratedPredictions(predsArea)
  
  # return results
  areaResults
}

# This is very similar to aggPixelPreds. The main difference is that 
# a given row of popNumerators and popDenominators may be associated 
# with different areas depending on the column/draw. Aggregating is 
# therefore not as simple as a few matrix multiplications. If every 
# column of areaMat is identical, use aggPixelPreds, since it is 
# faster. Typically used for aggregating EA level predictions to 
# pixel or subareal levels. Make sure areaMat is a factor with levels 
# equal to all possible areas even ones in which there are no 
# predictions to aggregate, and that the levels are in the correct 
# order
aggPredsVariablePerArea = function(popNumerators, popDenominators, 
                                   areaMat, areaLevels, urbanMat=NULL) {
  
  # function for aggregating values for each grid cell for draw i
  getDrawColumn = function(i, vals) {
    # calculate levels over which to aggregate
    indices = factor(as.character(areaMat[,i]), levels=areaLevels)
    
    # aggregate
    out = c(tapply(vals[,i], indices, FUN=sum, na.rm=TRUE))
    unlist(out)
  }
  
  # aggregate numerators and denominators, calculate prevalence/risk
  popNumerators[popDenominators == 0] = 0
  if(!is.null(urbanMat)) {
    # do the same for urban/rural strata:
    ## urban:
    popNumeratorsUrban = popNumerators * urbanMat
    popDenominatorsUrban = popDenominators * urbanMat
    NareaUrban <- sapply(1:ncol(popDenominators), getDrawColumn, vals=popDenominatorsUrban)
    NareaUrban[is.na(NareaUrban)] = 0
    ZareaUrban <- sapply(1:ncol(popNumerators), getDrawColumn, vals=popNumeratorsUrban)
    ZareaUrban[is.na(ZareaUrban)] = 0
    pAreaUrban = ZareaUrban / NareaUrban
    pAreaUrban[NareaUrban == 0] = NA
    
    ## rural:
    popNumeratorsRural = popNumerators * (!urbanMat)
    popDenominatorsRural = popDenominators * (!urbanMat)
    NareaRural <- sapply(1:ncol(popDenominators), getDrawColumn, vals=popDenominatorsRural)
    NareaRural[is.na(NareaRural)] = 0
    ZareaRural <- sapply(1:ncol(popNumerators), getDrawColumn, vals=popNumeratorsRural)
    ZareaRural[is.na(ZareaRural)] = 0
    pAreaRural = ZareaRural / NareaRural
    pAreaRural[NareaRural == 0] = NA
    
    ## combined:
    Narea = NareaUrban + NareaRural
    Zarea = ZareaUrban + ZareaRural
    pArea = Zarea / Narea
    pArea[Narea == 0] = NA
    
    list(region=areaLevels, Z=Zarea, N=Narea, p=pArea, 
         ZUrban=ZareaUrban, NUrban=NareaUrban, pUrban=pAreaUrban, 
         ZRural=ZareaRural, NRural=NareaRural, pRural=pAreaRural)
  } else {
    # in this case, we don't care about stratification, only the overall areal results
    Narea <- sapply(1:ncol(popDenominators), getDrawColumn, vals=popDenominators)
    Narea[is.na(Narea)] = 0
    Zarea <- sapply(1:ncol(popNumerators), getDrawColumn, vals=popNumerators)
    Zarea[is.na(Zarea)] = 0
    pArea = Zarea / Narea
    pArea[Narea == 0] = NA
    
    list(region=areaLevels, Z=Zarea, N=Narea, p=pArea)
  }
}

#' @describeIn aggPop Aggregate areal populations to another areal level
#' @export
areaPopToArea = function(areaLevelPop, areasFrom, areasTo, 
                         stratifyByUrban=TRUE, 
                         doFineScaleRisk=!is.null(areaLevelPop$aggregationResults$pFineScaleRisk), 
                         doSmoothRisk=!is.null(areaLevelPop$aggregationResults$pSmoothRisk), 
                         doGriddedRisk=!is.null(areaLevelPop$aggregationResults$pGriddedRisk)) {
  
  if(length(areasFrom) != length(unique(areasFrom))) {
    stop("areasFrom must contain only unique names of areas to which we want to aggregate")
  }
  
  uniqueNames = sort(unique(areasTo))
  
  # construct row of the aggregation matrix given toArea index
  getMatrixHelper = function(i) {
    # get areasTo associated with this fromArea
    thisToArea = uniqueNames[i]
    thisFromAreas = unique(areasFrom[areasTo == thisToArea])
    areaI = areasFrom %in% thisFromAreas
    
    areaI
  }
  
  # construct the aggregation matrix from areasFrom to areasTo
  A = t(matrix(sapply(1:length(uniqueNames), getMatrixHelper), ncol=length(uniqueNames)))
  rownames(A) = uniqueNames
  colnames(A) = areasFrom
  
  ##### aggregate populations
  # fine scale prevalence aggregation model
  getAggregationResults = function(resultNameRoot="FineScalePrevalence") {
    if(! (paste("N", resultNameRoot, sep="") %in% names(areaLevelPop$aggregationResults))) {
      stop(paste0(resultNameRoot, " was not computed in input areaLevelPop"))
    }
    
    capitalResultNameRoot = resultNameRoot
    
    nSamples = areaLevelPop$aggregationResults[[paste("N", resultNameRoot, sep="")]]
    zSamples = areaLevelPop$aggregationResults[[paste("Z", resultNameRoot, sep="")]]
    zSamples[is.na(zSamples)] = 0 # must set to zero temporarily so matrix multiplication works out
    
    ZAggregated =  A %*% zSamples
    NAggregated =  A %*% nSamples
    pAggregated = ZAggregated / NAggregated
    pAggregated[NAggregated == 0] = NA
    if(!is.null(areaLevelPop$aggregationMatrices)) {
      thisA=A %*% areaLevelPop$aggregationMatrices[[paste("A", resultNameRoot, sep="")]]
      rownames(thisA) = uniqueNames
    } else {
      thisA = NULL
    }
    
    if(stratifyByUrban) {
      nSamplesUrban = areaLevelPop$aggregationResults[[paste("NUrban", resultNameRoot, sep="")]]
      zSamplesUrban = areaLevelPop$aggregationResults[[paste("ZUrban", resultNameRoot, sep="")]]
      zSamplesUrban[is.na(zSamplesUrban)] = 0 # must set to zero temporarily so matrix multiplication works out
      
      nSamplesRural = areaLevelPop$aggregationResults[[paste("NRural", resultNameRoot, sep="")]]
      zSamplesRural = areaLevelPop$aggregationResults[[paste("ZRural", resultNameRoot, sep="")]]
      zSamplesRural[is.na(zSamplesRural)] = 0 # must set to zero temporarily so matrix multiplication works out
      
      ZAggregatedUrban =  A %*% zSamplesUrban
      NAggregatedUrban =  A %*% nSamplesUrban
      pAggregatedUrban = ZAggregatedUrban / NAggregatedUrban
      pAggregatedUrban[NAggregatedUrban == 0] = NA
      if(!is.null(areaLevelPop$aggregationMatrices)) {
        thisAUrban=A %*% areaLevelPop$aggregationMatrices[[paste("AUrban", resultNameRoot, sep="")]]
        rownames(thisAUrban) = uniqueNames
      } else {
        thisAUrban = NULL
      }
      
      ZAggregatedRural =  A %*% zSamplesRural
      NAggregatedRural =  A %*% nSamplesRural
      pAggregatedRural = ZAggregatedRural / NAggregatedRural
      pAggregatedRural[NAggregatedRural == 0] = NA
      if(!is.null(areaLevelPop$aggregationMatrices)) {
        thisARural=A %*% areaLevelPop$aggregationMatrices[[paste("ARural", resultNameRoot, sep="")]]
        rownames(thisARural) = uniqueNames
      } else {
        thisARural = NULL
      }
    } else {
      ZAggregatedUrban = NULL
      NAggregatedUrban = NULL
      pAggregatedUrban = NULL
      thisAUrban=NULL
      
      ZAggregatedRural = NULL
      NAggregatedRural = NULL
      pAggregatedRural = NULL
      thisARural=NULL
    }
    
    aggregationResults = list(region=uniqueNames, p=pAggregated, Z=ZAggregated, N=NAggregated, 
                              pUrban=pAggregatedUrban, ZUrban=ZAggregatedUrban, NUrban=NAggregatedUrban, 
                              pRural=pAggregatedRural, ZRural=ZAggregatedRural, NRural=NAggregatedRural)
    if(!is.null(thisA)) {
      aggregationMatrices = list(A=thisA, AUrban=thisAUrban, ARural=thisARural)
      names(aggregationMatrices) = paste(names(aggregationMatrices)[-1], capitalResultNameRoot, sep="")
    } else {
      aggregationMatrices = NULL
    }
    
    # capitalResultNameRoot = paste(toupper(substr(capitalResultNameRoot, 1, 1)), substr(capitalResultNameRoot, 2, nchar(capitalResultNameRoot)), sep="")
    names(aggregationResults)[-1] = paste(names(aggregationResults)[-1], capitalResultNameRoot, sep="")
    
    list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
  }
  
  # fine scale prevalence model 
  out = getAggregationResults("FineScalePrevalence")
  aggregationResults = out$aggregationResults
  aggregationMatrices = out$aggregationMatrices
  
  # fine scale risk aggregation model
  if(doFineScaleRisk) {
    out = getAggregationResults("FineScaleRisk")
    resFineScaleRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resFineScaleRisk, by="region")
    aggregationResults = c(aggregationResults, resFineScaleRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  if(doSmoothRisk) {
    out = getAggregationResults("SmoothRisk")
    resSmoothRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resSmoothRisk, by="region")
    aggregationResults = c(aggregationResults, resSmoothRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  if(doGriddedRisk) {
    out = getAggregationResults("GriddedRisk")
    resGriddedRisk = out$aggregationResults
    resAggregationMatrices = out$aggregationMatrices
    
    # aggregationResults = merge(aggregationResults, resGriddedRisk, by="region")
    aggregationResults = c(aggregationResults, resGriddedRisk)
    aggregationMatrices = c(aggregationMatrices, resAggregationMatrices)
  }
  
  list(aggregationResults=aggregationResults, aggregationMatrices=aggregationMatrices)
}

#' @describeIn simPop
#' Simulate populations and population prevalences given census frame and population density 
#' information. Uses custom spatial logit risk function and can include iid cluster 
#' level effect.
#' 
#' @export
simPopCustom = function(logitRiskDraws, sigmaEpsilonDraws, easpa, popMat, targetPopMat, 
                        stratifyByUrban=TRUE, validationPixelI=NULL, validationClusterI=NULL, 
                        clustersPerPixel=NULL, 
                        doFineScaleRisk=FALSE, doSmoothRisk=FALSE, doGriddedRisk=FALSE, 
                        doSmoothRiskLogisticApprox=FALSE, 
                        poppsub=NULL, subareaLevel=FALSE, gridLevel=FALSE, 
                        min1PerSubarea=TRUE, 
                        fixPopPerEA=NULL, fixHHPerEA=NULL, fixPopPerHH=NULL, 
                        returnEAinfo=FALSE, epsc=NULL, stopOnFrameMismatch=FALSE, tol=1e-3, 
                        verbose=TRUE) {
  
  if(!is.null(validationPixelI) || !is.null(validationClusterI) || !is.null(clustersPerPixel)) {
    stop("validationPixelI, validationClusterI, and clustersPerPixel not yet fully implemented")
  }
  
  if(any(c(!is.null(fixPopPerEA), !is.null(fixHHPerEA), !is.null(fixPopPerHH)))) {
    if(verbose) {
      print("adjusting easpa population and households based on fixPopPerEA, fixHHPerEA, and fixPopPerHH")
    }
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
  
  # make sure the population integration grids match the population frame
  out = checkPopFrameAndIntWeights(popMat=popMat, targetPopMat=targetPopMat, 
                                   easpa=easpa, poppsub=poppsub, 
                                   stopOnFrameMismatch=stopOnFrameMismatch, 
                                   tol=tol)
  popMat = out$popMat
  targetPopMat = out$targetPopMat
  
  nDraws = ncol(logitRiskDraws)
  
  # set default inputs
  totalEAs = sum(easpa$EATotal)
  if(!is.null(clustersPerPixel)) {
    emptyPixels = clustersPerPixel == 0
    if(totalEAs != sum(clustersPerPixel))
      stop("sum(easpa$EATotal) != sum(clustersPerPixel)")
  }
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  # determine if we care about subareas (smallest areas we care about. No info of EAs per subarea)
  sampleBySubarea = !is.null(popMat$subarea) && !is.null(poppsub)
  
  ##### Line 1 (of the algorithm): take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  print("drawing EAs")
  pixelProbs = popMat$pop
  
  # take draws from the stratified binomial process for each posterior sample
  if(is.null(clustersPerPixel)) {
    if(verbose) {
      print("Sampling EAs per pixel...")
    }
    if(sampleBySubarea) {
      eaSamples = SUMMER:::rStratifiedMultnomialBySubarea(nDraws, popMat, easpa, stratifyByUrban, poppsub=poppsub, 
                                                 min1PerSubarea=min1PerSubarea)
    } else {
      eaSamples = SUMMER:::rStratifiedMultnomial(nDraws, popMat, easpa, stratifyByUrban)
    }
  }
  
  if(!is.null(clustersPerPixel) && !exists("eaSamples")) {
    eaSamples = matrix(rep(clustersPerPixel, nDraws), ncol=nDraws)
  }
  
  # make matrix (or list) of pixel indices mapping matrices of EA values to matrices of pixel values
  if(!is.null(clustersPerPixel)) {
    pixelIndices = rep(1:nrow(popMat), times=clustersPerPixel) # this contains repetitions and has length == nEAs
    uniquePixelIndices = sort(unique(pixelIndices))
  } else {
    pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  }
  
  # determine which EAs are urban if necessary
  if(stratifyByUrban) {
    # urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
    if(!is.null(clustersPerPixel)) {
      urbanVals = popMat$urban[pixelIndices]
      uniqueUrbanVals = popMat$urban[uniquePixelIndices]
    }else {
      urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=nDraws)
    }
  } else {
    urbanMat = NULL
  }
  
  # determine which EAs are from which area
  if(!is.null(clustersPerPixel)) {
    areaVals = popMat$area[pixelIndices]
    uniqueAreaVals = popMat$area[uniquePixelIndices]
  } else {
    areaMat = matrix(popMat$area[pixelIndexMat], ncol=nDraws)
  }
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  if(verbose) {
    print("simulating EA level risks, numerators, and denominators")
  }
  if(is.null(epsc)) {
    epsc = matrix(stats::rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
  }
  
  ##### Line 3: draw EA population denominators, N
  
  if(!is.null(clustersPerPixel)) {
    if(is.null(validationPixelI))
      stop("clustersPerPixel must only be set for validation, but validationPixelI is NULL")
    
    # in this case, every left out cluster has exactly 25 households. Simply sample target population 
    # with equal probability from each cluster/faux EA
    if(is.null(fixPopPerEA)) {
      Ncs = SUMMER:::sampleNMultilevelMultinomialFixed(clustersPerPixel, nDraws=nDraws, pixelIndices=pixelIndices, 
                                                       urbanVals=urbanVals, areaVals=areaVals, easpa=easpa, popMat=popMat, stratifyByUrban=stratifyByUrban, 
                                                       verbose=verbose)
    } else {
      if(verbose) {
        print(paste0("drawing Ns for each EA..."))
      }
      Ncs = matrix(rep(fixPopPerEA, totalEAs*nDraws), ncol=nDraws)
    }
  } else {
    if(returnEAinfo) {
      out = SUMMER:::sampleNMultilevelMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                         popMat=popMat, stratifyByUrban=stratifyByUrban, verbose=verbose, returnEAinfo=returnEAinfo)
      householdDraws = out$householdDraws
      Ncs = out$targetPopDraws
    } else {
      Ncs <- SUMMER:::sampleNMultilevelMultinomial(pixelIndexMat=pixelIndexMat, urbanMat=urbanMat, areaMat=areaMat, easpaList=list(easpa), 
                                          popMat=popMat, stratifyByUrban=stratifyByUrban, verbose=verbose, returnEAinfo=returnEAinfo)
      householdDraws = NULL
    }
  }
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  if(verbose) {
    print("Calculating EA level risk and number of deaths...")
  }
  if(!is.null(clustersPerPixel)) {
    uc = logitRiskDraws[pixelIndices,]
    muc = expit(uc + epsc)
  } else {
    uc = matrix(logitRiskDraws[cbind(rep(rep(1:nrow(logitRiskDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
    muc = expit(uc + epsc)
  }
  
  # calculate Z_{ic} for each EA in each pixel
  Zcs = matrix(stats::rbinom(n=totalEAs * nDraws, size=Ncs, prob=as.matrix(muc)), ncol=nDraws)
  
  ##### Line 4: Aggregate appropriate values from EAs to the grid cell level
  
  # function for aggregating values for each grid cell
  getPixelColumnFromEAs = function(i, vals, applyFun=sum, popWeightMatrix=NULL) {
    # calculate levels over which to aggregate
    if(!is.null(clustersPerPixel)) {
      indices = pixelIndices
    } else {
      indices = factor(as.character(pixelIndexMat[,i]))
    }
    
    # in this case (the LCPb model), we calculate weighted means within factor levels using popWeightMatrix
    if(!is.null(popWeightMatrix)) {
      stop("using popWeightMatrix is no longer support, since this is much slower than calculating normalized 
           weights separately and multiplying values by them outside this function")
      # applyFun = function(x) {stats::weighted.mean(x, popWeightMatrix[,i], na.rm=TRUE)}
      
      Data = data.frame(v=vals[,i], w=popWeightMatrix[,i])
      out = sapply(split(Data, indices), function(x) stats::weighted.mean(x$v,x$w))
    } else {
      if(!is.null(clustersPerPixel)) {
        # out = tapply(vals[,i], factor(as.character(pixelIndices)), FUN=applyFun)
        out = tapply(vals[,i], indices, FUN=applyFun)
      } else {
        out = tapply(vals[,i], indices, FUN=applyFun)
      }
    }
    
    if(!is.null(clustersPerPixel)) {
      returnValues = out
    } else {
      indices = as.numeric(names(out))
      
      returnValues = rep(NA, nrow(logitRiskDraws))
      returnValues[indices] = out
    }
    
    returnValues
  }
  
  ##### Line 5: We already did this, resulting in logitRiskDraws input
  
  ##### Line 6: aggregate population denominators for each grid cell to get N_{ig}
  if(verbose) {
    print("Aggregating from EA level to the pixel level...")
  }
  
  # Ng <- sapply(1:ncol(Ncs), getPixelColumnFromEAs, vals=Ncs)
  # Ng[is.na(Ng)] = 0
  # 
  # ##### Line 7: aggregate response for each grid cell to get Z_{ig}
  # Zg <- sapply(1:ncol(Zc), getPixelColumnFromEAs, vals=Zc)
  # 
  # ##### Line 8: Calculate empirical mortality proportions for each grid cell, p_{ig}. 
  # #####         Whenever N_{ig} is 0, set p_{ig} to NA as well
  # pg = Zg / Ng
  # pg[Ng == 0] = NA
  
  # even if we aren't producing grid level results in general, gridded and smooth risk 
  # require it
  griddedRisk = NULL
  zSamplesGriddedRisk = NULL
  nSamplesGriddedRisk = NULL
  smoothRisk = NULL
  zSamplesSmoothRisk = NULL
  nSamplesSmoothRisk = NULL
  if(doGriddedRisk) {
    # NOTE: these cluster effects aren't necessarily in agreement with the 
    #       cluster effects of the fine scale models, since they're drawn 
    #       independently
    pixelEps = matrix(rnorm(length(logitRiskDraws), mean=0, sd=rep(sigmaEpsilonDraws, each=nrow(logitRiskDraws))), ncol=nDraws)
    griddedRisk = expit(logitRiskDraws + pixelEps)
    
    # get expected numerator and denominators
    nSamplesGriddedRisk = matrix(rep(targetPopMat$pop, nDraws), ncol=nDraws)
    zSamplesGriddedRisk = sweep(griddedRisk, 1, targetPopMat$pop, "*")
  }
  
  if(doSmoothRisk) {
    # integrate out the cluster effect in the risk
    smoothRisk = matrix(SUMMER::logitNormMean(cbind(c(as.matrix(logitRiskDraws)), rep(sigmaEpsilonDraws, each=nrow(logitRiskDraws))), logisticApprox=doSmoothRiskLogisticApprox), nrow=nrow(logitRiskDraws))
    
    # get expected numerator and denominators
    nSamplesSmoothRisk = matrix(rep(targetPopMat$pop, nDraws), ncol=nDraws)
    zSamplesSmoothRisk = sweep(smoothRisk, 1, targetPopMat$pop, "*")
  }
  
  ZcsFineScaleRisk = NULL
  NcsSamplesFineScaleRisk = NULL
  if(gridLevel) {
    # calculate pixel/grid level results
    out = aggPredsVariablePerArea(popNumerators=Zcs, popDenominators=Ncs, 
                                  areaMat=pixelIndexMat, areaLevels=1:nrow(popMat))
    Ng = out$N
    Zg = out$Z
    pg = out$p
    
    ##### calculate results for the other models if necessary
    
    if(doFineScaleRisk) {
      ## the following code is for a previous model that is no longer used, and so is commented out:
      # in order to get valid count estimates, we also need the expected denominator per EA in each stratum:
      # fineScaleRisk = sapply(1:ncol(muc), getPixelColumnFromEAs, vals=muc, applyFun=function(x) {mean(x, na.rm=TRUE)})
      # fineScaleRisk[!is.finite(fineScaleRisk)] = NA
      # nPerEA = SUMMER:::getExpectedNperEA(easpa, targetPopMat)
      # nSamplesFineScaleRisk = sweep(eaSamples, 1, nPerEA, "*")
      # zSamplesFineScaleRisk[is.na(zSamplesFineScaleRisk)] = 0
      
      # We aggregate just like for the fine scale prevalence model, but using the expected population 
      # numerators conditional on the true denominators.
      NcsFineScaleRisk = Ncs
      ZcsFineScaleRisk = muc * NcsFineScaleRisk
      out = aggPredsVariablePerArea(popNumerators=ZcsFineScaleRisk, popDenominators=Ncs, 
                                    areaMat=pixelIndexMat, areaLevels=1:nrow(popMat))
      nSamplesFineScaleRisk = out$N
      zSamplesFineScaleRisk = out$Z
      fineScaleRisk = out$p
    } else {
      nSamplesFineScaleRisk = NULL
      zSamplesFineScaleRisk = NULL
    }
    
    ##### Extra steps: collect draws at each level and generate:
    ##### areas, preds, 
    pixelLevelPop = list(pFineScalePrevalence=pg, ZFineScalePrevalence=Zg, NFineScalePrevalence=Ng, 
                         pFineScaleRisk=fineScaleRisk, ZFineScaleRisk=zSamplesFineScaleRisk, NFineScaleRisk=nSamplesFineScaleRisk, 
                         pSmoothRisk=smoothRisk, ZSmoothRisk=zSamplesSmoothRisk, NSmoothRisk=nSamplesSmoothRisk, 
                         pGriddedRisk=griddedRisk, ZGriddedRisk=zSamplesGriddedRisk, NGriddedRisk=nSamplesGriddedRisk)
  } else if(doGriddedRisk || doSmoothRisk) {
    pixelLevelPop = list(pFineScalePrevalence=NULL, ZFineScalePrevalence=NULL, NFineScalePrevalence=NULL, 
                         pFineScaleRisk=NULL, ZFineScaleRisk=NULL, NFineScaleRisk=NULL, 
                         pSmoothRisk=smoothRisk, ZSmoothRisk=zSamplesSmoothRisk, NSmoothRisk=nSamplesSmoothRisk, 
                         pGriddedRisk=griddedRisk, ZGriddedRisk=zSamplesGriddedRisk, NGriddedRisk=nSamplesGriddedRisk)
  } else {
    pixelLevelPop = NULL
  }
  
  if(subareaLevel) {
    # aggregate up from the most aggregated set of results we have: pixel/grid level if we have them, 
    # EA level results otherwise.
    if(verbose) {
      print("Aggregating to the subarea level...")
    }
    
    subareas = sort(unique(popMat$subarea))
    if(gridLevel) {
      # first get results for fine scale prevalence
      subareaLevelPop = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$subarea, 
                                       stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                                       doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                                       doGriddedRisk=doGriddedRisk)
    } else {
      # we aggregate EA level results to subarea level. This takes a bit more leg work
      subareaLevelPop = list()
      
      ## create matrices of relevant variable values at the EA level
      subareaMat = matrix(popMat$subarea[pixelIndexMat], ncol=ncol(pixelIndexMat))
      if(stratifyByUrban) {
        urbanMat = matrix(popMat$urban[pixelIndexMat], ncol=ncol(pixelIndexMat))
      } else {
        urbanMat = NULL
      }
      
      # fine scale prevalence
      fineScalePrevalenceSubarea = aggPredsVariablePerArea(popNumerators=Zcs, popDenominators=Ncs, 
                                                           areaMat=subareaMat, urbanMat=urbanMat, 
                                                           areaLevels=sort(as.character(poppsub$subarea)))
      names(fineScalePrevalenceSubarea)[-1] = paste(names(fineScalePrevalenceSubarea)[-1], "FineScalePrevalence", sep="")
      subareaLevelPop = c(subareaLevelPop, fineScalePrevalenceSubarea)
      
      # fine scale risk
      if(doFineScaleRisk) {
        NcsFineScaleRisk = Ncs
        ZcsFineScaleRisk = NcsFineScaleRisk * muc
        fineScaleRiskSubarea = aggPredsVariablePerArea(popNumerators=ZcsFineScaleRisk, 
                                                       popDenominators=NcsFineScaleRisk, 
                                                       areaMat=subareaMat, urbanMat=urbanMat, 
                                                       areaLevels=sort(as.character(poppsub$subarea)))
        names(fineScaleRiskSubarea)[-1] = paste(names(fineScaleRiskSubarea)[-1], "FineScaleRisk", sep="")
        subareaLevelPop = c(subareaLevelPop, fineScaleRiskSubarea)
      }
      
      # smooth risk
      if(doSmoothRisk || doGriddedRisk) {
        # we have pixel level results, so use them instead of EA level results since it's faster
        
        smoothGriddedRiskSubarea = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$subarea, 
                                           stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                                           doFineScalePrevalence=FALSE, doFineScaleRisk=FALSE, 
                                           doSmoothRisk=doSmoothRisk, doGriddedRisk=doGriddedRisk)
        smoothGriddedRiskSubarea$aggregationResults$region = NULL
        subareaLevelPop = c(subareaLevelPop, smoothGriddedRiskSubarea$aggregationResults)
      }
      
      # aggregation matrices are not applicable for aggregating from EA level to subarea level
      subareaLevelPop = list(aggregationResults=subareaLevelPop, aggregationMatrices=NULL)
    }
  } else {
    subareaLevelPop = NULL
  }
  
  # aggregate to the area level. Aggregate up from largest possible level to save time
  if(verbose) {
    print("Aggregating to the area level...")
  }
  if(subareaLevel) {
    areasFrom = subareaLevelPop$aggregationResults$region
    areasTo = poppsub$area[match(subareaLevelPop$aggregationResults$region, 
                                 poppsub$subarea)]
    areaLevelPop = areaPopToArea(subareaLevelPop, areasFrom=areasFrom, areasTo=areasTo, 
                                 stratifyByUrban=stratifyByUrban, doFineScaleRisk=doFineScaleRisk, 
                                 doSmoothRisk=doSmoothRisk, doGriddedRisk=doGriddedRisk)
  } else if(gridLevel) {
    areaLevelPop = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$areas, 
                                  stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                                  doFineScaleRisk=doFineScaleRisk, doSmoothRisk=doSmoothRisk, 
                                  doGriddedRisk=doGriddedRisk)
  } else {
    # In this case we must aggregate up from EA level for all risk/prevalence types except 
    # those only defined at the pixel level, and als aggregate those from the pixel level
    
    # First aggregate from EA level
    EALevelPop = list(pixelIndexMat=pixelIndexMat, 
                      ZFineScalePrevalence=Zcs, NFineScalePrevalence=Ncs, 
                      ZFineScaleRisk=ZcsFineScaleRisk, NFineScaleRisk=Ncs)
    
    areaLevelPop = EAPopToArea(EALevelPop, areaMat=areaMat, areaLevels=sort(as.character(poppsub$subarea)), 
                               urbanMat=urbanMat, doFineScalePrevalence=TRUE, 
                               doFineScaleRisk=doFineScaleRisk)
    
    # now aggregate grid level models
    if(doSmoothRisk || doGriddedRisk) {
      out = pixelPopToArea(pixelLevelPop, eaSamples=eaSamples, areas=popMat$area, 
                           stratifyByUrban=stratifyByUrban, targetPopMat=targetPopMat, 
                           doFineScalePrevalence=FALSE, doFineScaleRisk=FALSE, 
                           doSmoothRisk=doSmoothRisk, doGriddedRisk=doGriddedRisk)
      out = out$aggregationResults
      out$region = NULL
    }
    
    areaLevelPop = c(areaLevelPop, out)
  }
  
  if(!returnEAinfo) {
    list(pixelPop=pixelLevelPop, subareaPop=subareaLevelPop, areaPop=areaLevelPop)
  } else {
    
    # return list of eaDat objects
    getEADat = function(i) {
      theseI = pixelIndexMat[,i]
      
      eaDat = data.frame(lon=popMat$lon[theseI], lat=popMat$lat[theseI], 
                         area=popMat$area[theseI], subarea=rep("temp", length(theseI)), 
                         urban=popMat$urban[theseI], east=popMat$east[theseI], north=popMat$north[theseI], 
                         popDensity=popMat$pop[theseI], popDensityTarget=targetPopMat$pop[theseI], pixelIs=theseI, 
                         nHH=householdDraws[,i], N=Ncs[,i], Z=Zcs[,i], 
                         pFineScalePrevalence=Zcs[,i]/Ncs[,i])
      if(doFineScaleRisk) {
        eaDat$pFineScaleRisk = muc[,i]
        eaDat$ZFineScaleRisk = muc[,i] * Ncs[,i]
        eaDat$NFineScaleRisk = Ncs[,i]
      }
      if(doSmoothRisk) {
        eaDat$pSmoothRisk=smoothRisk[theseI,i]
        # N and Z not defined at the EA level
      }
      if(doGriddedRisk) {
        eaDat$pGriddedRisk=griddedRisk[theseI,i]
        # N and Z not defined at the EA level
      }
      
      if("subarea" %in% names(popMat)) {
        eaDat$subarea = popMat$subarea[theseI]
      } else {
        eaDat$subarea = NULL
      }
      
      eaDat$pFineScalePrevalence[eaDat$N == 0] = NA
      
      eaDat
    }
    
    if(verbose) {
      print("Constructing list of simulated EA data.frames...")
    }
    
    eaDatList = lapply(1:nDraws, getEADat)
    
    list(pixelPop=pixelLevelPop, subareaPop=subareaLevelPop, areaPop=areaLevelPop, 
         eaDatList=eaDatList, eaSamples=eaSamples)
  }
}

checkPopFrameAndIntWeights = function(popMat, targetPopMat, easpa, poppsub, stopOnFrameMismatch=TRUE, 
                                      tol=1e-3) {
  # first check the adjusted (target) population density
  temp = adjustPopGrid(targetPopMat, easpa)
  maxDiff = max(abs(temp$pop - targetPopMat$pop))
  targetPopMat$pop = temp$pop
  
  if(maxDiff > tol) {
    changedTargetPopMat = TRUE
    print(paste0("input targetPopMat did not match with ", 
                 "input easpa, with maximum pixel population ", 
                 "difference of ", maxDiff, ". Modifying ", 
                 "targetPopMat accordingly."))
    if(stopOnFrameMismatch) {
      stop(paste0("input targetPopMat did not match with ", 
                  "input easpa, with maximum pixel population ", 
                  "difference of ", maxDiff, ". Modifying ", 
                  "targetPopMat accordingly."))
    }
    warning(paste0("input targetPopMat did not match with ", 
                   "input easpa, with maximum pixel population ", 
                   "difference of ", maxDiff, ". Modifying ", 
                   "targetPopMat accordingly."))
  } else {
    changedTargetPopMat = FALSE
  }
  
  # aggregate general population density
  stratifyByUrban = "urban" %in% names(popMat)
  temp = adjustPopMat(popMat, poppsub, adjustBy="subarea", stratifyByUrban=stratifyByUrban)
  maxDiff = max(abs(temp$pop - popMat$pop))
  popMat$pop = temp$pop
  if(maxDiff > tol) {
    changedPopMat = TRUE
    print(paste0("input popMat did not match with ", 
                 "input easpa, which maximum pixel population", 
                 "difference of ", maxDiff, ". Modifying ", 
                 "popMat accordingly."))
    if(stopOnFrameMismatch) {
      stop(print(paste0("input popMat did not match with ", 
                        "input easpa, which maximum pixel population", 
                        "difference of ", maxDiff, ". Modifying ", 
                        "popMat accordingly.")))
    }
    warning(print(paste0("input popMat did not match with ", 
                         "input easpa, which maximum pixel population", 
                         "difference of ", maxDiff, ". Modifying ", 
                         "popMat accordingly.")))
  } else {
    changedPopMat = FALSE
  }
  
  list(popMat=popMat, targetPopMat=targetPopMat, changed=list(changedPopMat=changedPopMat, changedTargetPopMat=changedTargetPopMat))
}






