### Aggregation Models
## Each aggregation model must at least use the following inputs: 
# uDraws: posterior draws from the spatial field. nLocs x nsim dimension matrix.
# sigmaEpsilonDraws: posterior draws of cluster effect SD (joint draws with uDraws). Of length nsim
# easpa: table of number of EAs per area similar to easpc
# popMat: nLoc row matrix of population densities/EA densities (not necessarily normalized) and strata 
#         for pixel grid such as from the popGrid variable. Eventually, this may include 'subarea' ID 
#         for aggregation to
# empiricalDistributions: A list of empirical distributions for the EA denominators for each urban/rural 
#                         value 
# includeUrban: if TRUE, uses EAsUrb and EAsRur in easpa.  Otherwise uses EAsTotal
## In addition, here are a few other possible inputs:
# maxEmptyFraction: the maximum fraction of posterior samples for which a given subarea can have 
#                   no EAs before separate marginal approximations are performed. If 1, never 
#                   perform separate marginal approximations
## NOTE1: For aggregation consistent models: 'areas' should be the finest spatial areas for which number 
##        of urban and rural EAs is known. For coarser areas, further aggregate the draws of each 
##        aggregation function as necessary. The same goes for the pixel grid. For larger 'subareas' 
##        containing an unknown number of EAs, in which case aggregate pixel predictions as need be. 
## NOTE2: Aggregation consistent models:
#           pixels/subareas (unknown EAs): CPBL, cpbl (any others?)
#           areas/regions (known EAs): all of them?

# Use the popGrid variable to construct a matrix suitable for input into the aggregation models. 
# population matrix must have the following values: 
# lon: longitude
# lat: latitude
# east: easting (km)
# north: northing (km)
# pop: proportional to population density for each grid cell
# area: an id or area name in which the grid cell corresponding to each row resides
# urban: whether the grid cell is urban or rural
makeDefaultPopMat = function() {
  out = popGrid
  out$pop = out$popOrig
  out$area = out$admin1
  out$popOrig = NULL
  out$admin1 = NULL
  
  out
}

# Use the easpc variable to construct a table suitable for input into the aggregation models. 
# output must have the following values: 
# area: the name or id of the area
# EAUrb: the number of EAs in the urban part of the area
# EARur: the number of EAs in the rural part of the area
# EATotal: the number of EAs in the the area
# Also, the rows of this table will be ordered by Area
makeDefaultEASPA = function() {
  out = easpc
  out$pop = out$popOrig
  out$popOrig = NULL
  names(out) = c("area", "EAUrb", "EARur", "EATotal")
  
  # sort by area
  sortI = sort(out$area, index.return=TRUE)$ix
  out = out[sortI,]
  
  out
}

# in this model, we assume a binomial process for the EA locations. Follows algorithm 2 from the 
# outline
modCPBL = function(uDraws, sigmaEpsilonDraws, easpa=NULL, popMat=NULL, empiricalDistributions=NULL, 
                   includeUrban=TRUE, maxEmptyFraction=1) {
  nDraws = ncol(uDraws)
  
  # set default inputs
  if(is.null(easpa)) {
    easpa = makeDefaultEASPA()
    # Area: the name or id of the area
    # EAUrb: the number of EAs in the urban part of the area
    # EARur: the number of EAs in the rural part of the area
    # EATotal: the number of EAs in the the area
  }
  totalEAs = sum(easpa$EATotal)
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
  
  # get area names
  areas = sort(unique(popMat$area))
  if(any(areas != easpa$area))
    stop("area names and easpa do not match popMat or are not in the correct order")
  
  ##### Line 1: take draws from the binomial process for each stratum (each row of easpa)
  # get probabilities for each pixel (or at least something proportional within each stratum)
  pixelProbs = popMat$pop
  rMyMultinomial = function(n, i, urban=TRUE) {
    # determine which pixels and how many EAs are in this stratum
    if(includeUrban) {
      includeI = popMat$area == areas[i] & popMat$urban == urban
      nEA = ifelse(urban, easpa$EAUrb[i], easpa$EARur[i])
    }
    else {
      includeI = popMat$area == areas[i]
      nEA = ifelse(urban, easpa$EAUrb[i], easpa$EATotal[i])
    }
    
    # sample from the pixels if this stratum exists
    if(sum(includeI) == 0)
      return(matrix(nrow=0, ncol=n))
    thesePixelProbs = popMat$pop[includeI]
    rmultinom(n, nEA, prob=thesePixelProbs)
  }
  
  # function for determining how to recombine separate multinomials into the draws over all pixels
  getSortIndices = function(i, urban=TRUE) {
    # determine which pixels and how many EAs are in this stratum
    if(includeUrban) {
      includeI = popMat$area == areas[i] & popMat$urban == urban
    }
    else {
      includeI = popMat$area == areas[i]
    }
    
    which(includeI)
  }
  
  # gives nPixels x n matrix of draws from the stratified multinomial with values 
  # corresponding to the value of |C^g| for each pixel, g (the number of EAs/pixel)
  rStratifiedMultnomial = function(n) {
    # we will need to draw separate multinomial for each stratum. Start by 
    # creating matrix of all draws of |C^g|
    eaSamples = matrix(NA, nrow=nrow(popMat), ncol=n)
    
    # now draw multinomials
    if(includeUrban) {
      # draw for each area crossed with urban/rural
      urbanSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=TRUE))
      ruralSamples = do.call("rbind", lapply(1:length(areas), rMyMultinomial, n=n, urban=FALSE))
      
      # get the indices used to recombine into the full set of draws
      urbanIndices = unlist(sapply(1:length(areas), getSortIndices, urban=TRUE))
      ruralIndices = unlist(sapply(1:length(areas), getSortIndices, urban=FALSE))
      
      # recombine into eaSamples
      eaSamples[urbanIndices,] = urbanSamples
      eaSamples[ruralIndices,] = ruralSamples
    } else {
      # draw for each area
      stratumSamples = rbind(sapply(1:length(areas), n=n, rMyMultinomial))
      
      # get the indices used to recombine into the full set of draws
      stratumIndices = c(sapply(1:length(areas), getSortIndices))
      
      # recombine into eaSamples
      eaSamples[stratumIndices,] = stratumSamples
    }
    
    # return results
    eaSamples
  }
  
  # take draws from the stratified binomial process for each posterior sample
  eaSamples = rStratifiedMultnomial(nDraws)
  
  # determine which EAs are urban if necessary
  if(includeUrban) {
    urbanMat = matrix(rep(rep(popMat$urban, nDraws), times=c(eaSamples)), ncol=nDraws)
  }
  
  # make matrix of pixel indices mapping matrices of EA values to matrices of pixel values
  pixelIndexMat = matrix(rep(rep(1:nrow(popMat), nDraws), times=eaSamples), ncol=nDraws)
  
  ##### Line 2: draw cluster effects, epsilon
  # NOTE1: we assume there are many more EAs then sampled clusters, so that 
  #       the cluster effects for each EA, including those sampled, are iid
  epsc = matrix(rnorm(totalEAs*nDraws, sd=rep(sigmaEpsilonDraws, each=totalEAs)), ncol=nDraws)
  
  ##### Line 3: draw EA population denominators, N
  # NOTE1: we assume the values for each EA within urban/rural boundaries are iid
  if(includeUrban) {
    if(!is.null(fastDistributions$popUrban)) {
      NcsUrban = matrix(recdf(sum(easpa$EAUrb), fastDistributions$popUrban), ncol=nDraws)
      NcsRural = matrix(recdf(sum(easpa$EARur), fastDistributions$popRural), ncol=nDraws)
    } else {
      NcsUrban = matrix(recdfComposed(sum(easpa$EAUrb)*nDraws, list(fastDistributions$householdsUrban, 
                                                                    fastDistributions$mothersUrban, 
                                                                    fastDistributions$childrenUrban)), 
                        ncol=nDraws)
      NcsRural = matrix(recdfComposed(sum(easpa$EARur)*nDraws, list(fastDistributions$householdsRural, 
                                                                    fastDistributions$mothersRural, 
                                                                    fastDistributions$childrenRural)), 
                        ncol=nDraws)
    }
    
    Ncs[urbanMat] = NcsUrban
    Ncs[!urbanMat] = NcsRural
  } else {
    if(!is.null(fastDistributions$popTotal)) {
      Ncs = matrix(recdf(sum(easpa$EAsTotal), fastDistributions$popTotal), ncol=nDraws)
    } else {
      Ncs = matrix(recdfComposed(sum(easpa$EAsTotal), list(fastDistributions$households, 
                                                         fastDistributions$mothers, 
                                                         fastDistributions$children)), 
                   ncol=nDraws)
    }
  }
  
  ##### do part of Line 7 in advance
  # calculate mu_{ic} for each EA in each pixel
  browser()
  uc = matrix(uDraws[cbind(rep(rep(1:nrow(uDraws), nDraws), times=c(eaSamples)), rep(1:nDraws, each=totalEAs))], ncol=nDraws)
  muc = expit(uc + epsc)
  
  # calculate Z_{ic} for each EA in each pixel
  Zc = matrix(rbinom(totalEAs * nDraws, n=Ncs, prob=muc), ncol=nDraws)
  
  ##### Line 4: calculate the empirical mortality proportion, p_{ig} for each grid cell and draw
  # calculate population total for each grid cell, N_{ig}
  Ng = 
  
  getPig = function(c) {
    # get the number of enumeration areas
    nea = eaSamples[,c]
    
    ##### Line 5: get the value of the spatial term, u. This is already done via uDraws
    
    ##### Line 6: draw population denominators for each EA and pixel depending on pixel stratum
    if(includeUrban) {
      if(!is.null(fastDistributions$popUrban)) {
        NcsUrban = lapply(nea[popMat$urban], recdf, fastDistributions$popUrban)
        NcsRural = lapply(nea[!popMat$urban], rcedf, fastDistributions$popRural)
      } else {
        NcsUrban = lapply(nea[popMat$urban], recdfComposed, list(fastDistributions$householdsUrban, 
                                                                    fastDistributions$mothersUrban, 
                                                                    fastDistributions$childrenUrban))
        NcsRural = lapply(nea[!popMat$urban], recdfComposed, list(fastDistributions$householdsRural, 
                                                                     fastDistributions$mothersRural, 
                                                                     fastDistributions$childrenRural))
      }
      Ngs = numeric(nrow(popMat))
      Ngs[popMat$urban] = c(sapply(NcsUrban, sum))
      Ngs[!popMat$urban] = c(sapply(NcsRural, sum))
    } else {
      if(!is.null(fastDistributions$popTotal)) {
        Ncs = lapply(nea, rcedf, fastDistributions$popTotal)
      } else {
        Ncs = lapply(nea, recdfComposed, list(fastDistributions$households, 
                                              fastDistributions$mothers, 
                                              fastDistributions$children))
      }
      Ngs = c(sapply(Ncs, sum))
    }
    
    ##### Line 7: calculate pixel numerators, Zig
    # we already calculated  mu_{ic} for each EA in each pixel. Now draw binomials
    lapply(1:length(Ngs), function(i) {
      rbinom(Ncs[i], length(Ncs[i]), muc[i])
    })
  }
}









