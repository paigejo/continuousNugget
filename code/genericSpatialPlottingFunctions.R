# helpful spatial plotting functions using base R and fields:

# continuous plotting ----
# x, y: horizontal and vertical spatial coordinates
# z: response
# zlim: range of the response
# cols: color vector representing the color scale
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# legend.mar, n.ticks, legend.width, legendArgs (as legend.args): see ?image.plot
# min.n: approximate number of ticks in color scale. See ?pretty
# orderI: a specific ordering to plot the points in
# ordering: in what order to plot the points, where the ordering is based on the response z
plotWithColor = function(x, y, z, zlim=NULL, colScale=tim.colors(), 
                         legend.mar=7, new=TRUE, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                         n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, legend.width=1.2, addColorBar=TRUE, 
                         legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, orderI=NULL, 
                         ordering=c("none", "increasing", "decreasing"), colorName = c("col", "bg"), ...) {
  ordering = match.arg(ordering)
  colorName = match.arg(colorName)
  
  # remove NA points
  nas = is.na(x) | is.na(y) | is.na(z)
  if(any(nas)) {
    warning("Removing NAs")
    x = x[!nas]
    y = y[!nas]
    z = z[!nas]
  }
  
  # do setup for ploting data if necessary
  if(is.null(zlim)) {
    nas = !is.finite(scaleFun(z))
    zlim = range(z[!nas])
  }
  
  # order the plotting of the points
  if(is.null(orderI)) {
    if(ordering == "increasing") {
      orderI = sort(z, index.return=TRUE)$ix
    } else if(ordering == "decreasing") {
      orderI = sort(z, decreasing=TRUE, index.return=TRUE)$ix
    } else {
      orderI = 1:length(z)
    }
  }
  x = x[orderI]
  y = y[orderI]
  z = z[orderI]
  
  # if(forceColorsInRange) {
  #   z[z > zlim[2]] = zlim[2]
  #   z[z < zlim[1]] = zlim[1]
  # }
  
  # get colors of points
  cols = getColorsFromScale(z, zlim, cols=colScale, scaleFun=scaleFun, 
                            forceValuesInRange=forceColorsInRange)
  
  # generate new plot if necessary
  # browser()
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    if(colorName == "col") {
      do.call("plot", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("plot", c(list(x=x, y=y, bg=cols), list(...)))
    }
  } else {
    if(colorName == "col") {
      do.call("points", c(list(x=x, y=y, col=cols), list(...)))
    } else {
      do.call("points", c(list(x=x, y=y, bg=cols), list(...)))
    }
  }
  
  if(addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(tickLabels))
      setTickLabels = TRUE
    
    if(is.null(ticks)) {
      if(setTickLabels)
        tickLabels = pretty(zlim, n=n.ticks, min.n=min.n)
      ticks = scaleFun(tickLabels)
    }
    else {
      if(setTickLabels)
        tickLabels = ticks
      ticks = scaleFun(ticks)
    }
    if(setTickLabels)
      tickLabels = tickLabels[is.finite(ticks)]
    ticks = ticks[is.finite(ticks)]
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=scaleFun(zlim)
    legendArgs$nlevel=length(colScale)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=colScale
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

# TODO: make nice version of quilt.plot including scales for color scale
# base it on this example. Make sure color scale range is centered correctly for
# diverging color scales
# quilt.plot(popMat$lon, popMat$lat, pixelMean, FUN=function(x){log(mean(x, na.rm=TRUE))}, 
#            zlim=log(meanRangePixel), nx=160, ny=160, main="", cex.main=3, col=meanCols, 
#            add.legend=FALSE, cex.axis=2, xlab="", ylab="Latitude", 
#            xlim=kenyaLonRange, ylim=c(-5.5, 5.8), asp=1, cex.lab=3)
# plotMapDat(mapDat=countyMap, lwd=.5, border=rgb(.4,.4,.4))
# points(mort$lon, mort$lat, pch=19, cex=.1)
# 
# meanTicksPixel = getLogScaleTicks(meanRangePixel)
# meanTickLabelsPixel = as.character(meanTicksPixel)
# meanTicks = pretty(meanRange, n=5)
# meanTickLabels = as.character(meanTicks)
# image.plot(zlim=range(log(meanRangePixel)), nlevel=length(meanCols), legend.only=TRUE, horizontal=FALSE,
#            col=meanCols, add = TRUE, axis.args=list(at=log(meanTicksPixel), labels=meanTickLabelsPixel, cex.axis=2, tck=-.7, hadj=-.1),
#            legend.mar = 0, legend.cex=2, legend.width=3, smallplot=c(.88,.91,.1,.9))

# discrete plotting ----

# for plotting areal values assuming plotVar is in alphabetical order of the area names. 
# mapDat: a SpatialPolygonsDataFrame object, where each polygon may have a value you wish to plot
# plotVar: if null, plot the areal boundaries only. Else plot plotVar values for each area
# varAreas: area names associated with plotVar
# zlim: range of the response
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using myProjection function.  This can be used when plotting the projected `easting' 
#          and `northing' variables for instance.
# cols: color vector representing the color scale
# legend.mar, legend.args, n.ticks: see ?image.plot
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# crosshatchNADensity: Adds crosshatching for areas with NA values. See ?polygon density argument.
# # min.n: approximate number of ticks in color scale. See ?pretty
# myProjection: a map projection function taking a 2 column matrix of coordinates 
#   and projects them.
# ...: arguments to polygon function
plotMapDat = function(mapDat, plotVar=NULL, varAreas, regionNames=sort(unique(varAreas)), zlim=NULL, project=FALSE, cols=tim.colors(), 
                      legend.mar=7, new=TRUE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                      ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, addColorBar=TRUE, 
                      legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, 
                      crosshatchNADensity=10, myProjection=NULL, ...) {
  require(fields)
  
  # do setup for plotting data by area if necessary
  if(!is.null(plotVar)) {
    if(is.null(zlim)) {
      zlim = range(plotVar)
    }
    
    if(forceColorsInRange) {
      plotVar[plotVar > zlim[2]] = zlim[2]
      plotVar[plotVar < zlim[1]] = zlim[1]
    }
  }
  
  # generate new plot if necessary
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(project) {
      if(is.null(xlab))
        xlab = "East (km)"
      if(is.null(xlim))
        xlim = eastLim
      if(is.null(ylab))
        ylab = "North (km)"
      if(is.null(ylim))
        ylim = northLim
    }
    else {
      if(is.null(xlab))
        xlab = "Longitude"
      if(is.null(ylab))
        ylab = "Latitude"
      if(is.null(xlim)) {
        xlim = mapDat@bbox[1,]
      }
      if(is.null(ylim)) {
        ylim = mapDat@bbox[2,]
      }
    }
    if(is.null(main))
      main = ""
    
    if(is.null(plotArgs)) {
      plotArgs = list(main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=asp)
    } else {
      plotArgs$main = main
      plotArgs$xlab = xlab
      plotArgs$ylab = ylab
      plotArgs$xlim = xlim
      plotArgs$ylim = ylim
      plotArgs$asp = asp
    }
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    do.call("plot", c(list(1, 2, type="n"), plotArgs))
  }
  
  # add polygons to plot
  polys = mapDat@polygons
  plotArea = function(i) {
    areaPolys = polys[[i]]@Polygons
    
    if(is.null(plotVar)) {
      if(!project)
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(areaPolys[[x]]@coords), list(...)))})
      else
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(myProjection(areaPolys[[x]]@coords)), list(...)))})
    }
    else {
      # get index of plotVar corresponding to this area
      thisI = which(varAreas == regionNames[i])
      
      # if there is no matching region name, do nothing
      if(length(thisI) == 0) {
        return(NULL)
      }
      
      # get color to plot
      vals = c(zlim, scaleFun(plotVar[thisI]))
      vals = vals-zlim[1]
      vals = vals/(zlim[2] - zlim[1])
      col = cols[round(vals[3]*(length(cols)-1))+1]
      if(is.na(vals[3])) {
        thisDensity = crosshatchNADensity
      } else {
        thisDensity = NULL
      }
      
      if(!project)
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(areaPolys[[x]]@coords, col=col, density=thisDensity), list(...)))})
      else
        sapply(1:length(areaPolys), function(x) {do.call("polygon", c(list(myProjection(areaPolys[[x]]@coords), col=col, density=thisDensity), list(...)))})
    }
    
  }
  
  sapply(1:length(polys), plotArea)
  
  if(!is.null(plotVar) && addColorBar) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(scaleFunInverse(zlim), n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    
    # par( oma=c( 0,0,0,3))
    
    # set list of arguments to image.plot
    legendArgs$zlim=zlim
    legendArgs$nlevel=length(cols)
    legendArgs$legend.only=TRUE
    legendArgs$horizontal=FALSE
    legendArgs$col=cols
    legendArgs$add = TRUE
    if(is.null(legendArgs$axis.args))
      legendArgs$axis.args=list(at=ticks, labels=tickLabels)
    else {
      legendArgs$axis.args$at=ticks
      legendArgs$axis.args$labels=tickLabels
    }
    legendArgs$legend.mar=legend.mar
    legendArgs$legend.width=legend.width
    do.call("image.plot", legendArgs)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

# for plotting administration area names (for use with plotMapDat)
# varAreas: the names of the areas we wish to plot
# mapDat: spatial polygon file containing information on the areas
# ...: additional arguments to the text function
addMapLabels = function(varAreas, mapDat, offsets=NULL, areaVarName, ...) {
  
  # determine the names of the mapDat areas
  regionNames = as.character(mapDat@data[[areaVarName]])
  
  # plot map labels
  xs = coordinates(mapDat)
  includeI = match(varAreas, regionNames)
  xs = xs[includeI,]
  
  if(!is.null(offsets)) {
    xs = xs + offsets
  }
  
  text(xs, as.character(varAreas), ...)
  
  invisible(NULL)
}

# color scales ----
makePurpleYellowSequentialColors = function(n, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev)
  else
    scale_colour_continuous_sequential(h1=-100, h2=100, c1=60, cmax=74, c2=100, l1=15, l2=95, p1=2, p2=0.9, rev=rev, n_interp=n)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev, n_interp=n, mid=center)
    }
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange) && is.null(center)) {
    if(!ggplot)
      diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    else
      scale_colour_continuous_diverging(n_interp=n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange if center is NULL
    if(!ggplot) {
      propUp = (valRange[2] - center) / diff(valRange)
      propDown = 1 - propUp
      totalColors = ceiling(2 * max(propUp, propDown) * n)
      tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
      totalMissingColors = totalColors - n
      
      if(propUp >= propDown && totalMissingColors > 0)
        tempColors[-(1:totalMissingColors)]
      else
        tempColors[1:n]
    } else {
      if(is.null(center))
        center = min(valRange) + abs(diff(valRange))/2
      scale_colour_continuous_diverging(n.interp, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev, mid=center)
    }
  }
}

makeBlueSequentialColors = function(n, ggplot=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  if(!ggplot)
    sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
  else
    scale_colour_continuous_sequential(h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE, n_interp=n)
}

makeBlueGreenYellowSequentialColors = function(n, ggplot=FALSE, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(!ggplot)
    sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, rev=rev)
  else
    scale_colour_continuous_sequential(h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1, n_interp=n, rev=rev)
}

makeRedYellowBlueColors = function(n, ggplot=FALSE) {
  if(!ggplot)
    divergingx_hcl(n, palette="RdYlBu")
  else
    scale_colour_continuous_sequential(palette="RdYlBu", n_interp=n)
}

# color scale utility functions ----
# combine two color scale functions (that return vector of colors given number of colors), 
# given the number of colors in the scale desired
combineTwoScales = function(n, scale1, scale2, args1, args2) {
  if(n %% 2 == 0) {
    n1 = n2 = n/2
  } else {
    n1 = ceiling(n/2)
    n2 = floor(n/2)
  }
  
  c(do.call(scale1, c(args1, list(n=n1))), 
    do.call(scale2, c(args2, list(n=n2))))
}

# convert a single color sequential scale into a diverging scale
makeDivergingScale = function(n, scale, ...) {
  do.call("combineTwoScales", list(n=n, scale1=scale, scale2=scale, args1=list(...), args2=list(...)))
}

# centers a color scale at its midpoint. Returns vector of the centered color scale. 
# Useful when using diverging scales centered at 0 for data with asymmetric range
# colScale a function taking 'n' as input and return a color scale centered in the middle
# n: number of colors
# vals: values to plot (unscaled)
# valRange: the range of vals you want to plot
# center: the center of the range of values (this is the center of the diverging color scale)
# colScale: function return a vector of colors taking inputs 'n' and '...', e.g. makeRedBlueDivergingColors
# scaleFun: scales the color scale to be on, e.g., log scale
# ...: other arguments to colScale, such as 'rev'
# Example:
# # construct data:
# testX = exp(rnorm(100))
# # construct centered color scale on log scale
# test = centerColorScale(64, testX, center=1, colScale=makeRedBlueDivergingColors, scaleFun=log)
# # get the colors associated with each value
# testCols = getColorsFromScale(testX, center=1, cols=test, scaleFun=log)
# # plot the result
# plot(testX, col=testCols, pch=19)
centerColorScale = function(n, vals=NULL, valRange=NULL, center, colScale, scaleFun=function(x) {x}, 
                            ...) {
  require("colorspace")
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  valRange = scaleFun(valRange)
  center = scaleFun(center)
  
  propUp = (valRange[2] - center) / diff(valRange)
  propDown = 1 - propUp
  totalColors = ceiling(2 * max(propUp, propDown) * n)
  tempColors = do.call(colScale, c(list(totalColors), list(...)))
  totalMissingColors = totalColors - n
  
  if(propUp >= propDown && totalMissingColors > 0)
    tempColors[-(1:totalMissingColors)]
  else
    tempColors[1:n]
}

# Long story short, use this function when using a log, logit, or some other fancy 
# transformation of a color scale.
# given data to be assigned colors from scale (vals), the range of the values, the center 
# of the value range (only for for diverging scale), the color scale as a vector 
# of colors (cols), the scale of the color scale (e.g. identity, log, logit), 
# and whether or not to force the colors into valRange, this function returns the 
# colors associated with each value in vals.
# Example:
# # construct data:
# testX = exp(rnorm(100))
# # construct centered color scale on log scale
# test = centerColorScale(64, testX, center=1, colScale=makeRedBlueDivergingColors, scaleFun=log)
# # get the colors associated with each value
# testCols = getColorsFromScale(testX, center=1, cols=test, scaleFun=log)
# # plot the result
# plot(testX, col=testCols, pch=19)
getColorsFromScale = function(vals, valRange=NULL, center=NULL, cols, scaleFun=function(x) {x}, 
                              forceValuesInRange=FALSE) {
  
  if(is.null(valRange)) {
    nas = !is.finite(scaleFun(vals))
    valRange = range(vals[!nas])
  }
  
  if(forceValuesInRange) {
    vals[vals < valRange[1]] = valRange[1]
    vals[vals > valRange[2]] = valRange[2]
  }
  
  valRange = scaleFun(valRange)
  vals = scaleFun(vals)
  vals = vals - valRange[1]
  vals = vals/(valRange[2] - valRange[1])
  
  if(!is.null(center)) {
    center = scaleFun(center)
    n = length(cols)
    
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = cols
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
    
    cols = tempColors
  }
  
  col = cols[round(vals*(length(cols)-1))+1]
  
  col
}

# axis scale ticks ---- 
# x: numbers on probability scale
getLogitScaleTicks = function(x, nint=3, add.5=FALSE) {
  minX = min(x)
  maxX = max(x)
  
  # check if range contains .5. If so, make sure .5 is in range of data
  if(minX <= .5 && maxX >= .5) {
    x = c(x, .5)
  }
  
  # first generate ticks below .5, then flip about .5
  rx = x
  rx[rx > .5] = 1 - rx[rx > .5]
  
  # now add log scale ticks
  lowerTicks = axisTicks(range(log10(rx)), log=TRUE, nint=nint)
  upperTicks = rev(1 - lowerTicks)
  
  if(add.5) {
    c(lowerTicks, .5, upperTicks)
  } else {
    c(lowerTicks, upperTicks)
  }
}

# x: numbers on positive real scale
getLogScaleTicks = function(x, nint=5) {
  axisTicks(range(log10(x)), log=TRUE, nint=nint)
}

## Functions from jamba package ----
getCustomScaleTicks = function(usr, scaleFun=sqrt, nint=5, log=FALSE) {
  # axisTicks(range(scaleFun(x)), log=log, nint=nint)
  axp <- unlist(.axisPars(range(scaleFun(usr)), log = log, nintLog = nint), 
                use.names = FALSE)
  .Call(C_R_CreateAtVector, axp, usr, nint, log)
}


#' Determine square root axis tick mark positions
#'
#' Determine square root axis tick mark positions
#'
#' This function calculates positions for tick marks for data
#' that has been transformed with `sqrt()`, specifically a directional
#' transformation like `sqrt(abs(x)) * sign(x)`.
#'
#' The main goal of this function is to provide reasonably placed
#' tick marks using integer values.
#'
#' @return
#' Invisibly returns a numeric vector of axis tick positions,
#' named by the display label.
#' The axis values are in square root space while the labels represent
#' the normal space values.
#'
#' @family jam plot functions
#'
#' @param side integer value indicating the axis position, as used
#'    by `axis()`, 1=bottom, 2=left, 3=top, 4=right.
#' @param x optional numeric vector representing the numeric range
#'    to be labeled.
#' @param pretty.n numeric value indicating the number of desired
#'    tick marks, passed to `pretty()`.
#' @param u5.bias numeric value passed to `pretty()` to influence the
#'    frequency of intermediate tick marks.
#' @param big.mark character value passed to `format()` which helps
#'    visually distinguish numbers larger than 1000.
#' @param plot logical indicating whether to plot the axis tick
#'    marks and labels.
#' @param las,cex.axis numeric values passed to `axis()` when drawing
#'    the axis, by default `las=2` plots labels rotated
#'    perpendicular to the axis.
#' @param ... additional parameters are passed to `pretty()`.
#'
#' @export
sqrtAxis <- function
(side=1,
 x=NULL,
 pretty.n=10,
 u5.bias=0,
 big.mark=",",
 plot=TRUE,
 las=2,
 cex.axis=0.6,
 ...)
{
  ## Purpose is to generate a set of tick marks for sqrt
  ## transformed data axes.  It assumes data is already sqrt-transformed,
  ## and that negative values have been treated like:
  ## sqrt(abs(x))*sign(x)
  if (length(side) > 2) {
    x <- side;
    side <- 0;
  }
  if (length(side) == 0) {
    side <- 0;
  }
  if (1 %in% side) {
    xRange <- par("usr")[1:2];
  } else if (2 %in% side) {
    xRange <- par("usr")[3:4];
  } else if (length(x) > 0) {
    xRange <- range(x, na.rm=TRUE);
  }
  
  subdivideSqrt <- function(atPretty1, n=pretty.n, ...) {
    ## Purpose is to take x in form of 0,x1,
    ## and subdivide using pretty()
    atPretty1a <- unique(sort(abs(atPretty1)));
    atPretty1b <- tail(atPretty1a, -2);
    atPretty2a <- pretty(head(atPretty1a,2), n=n, ...);
    return(unique(sort(c(atPretty2a, atPretty1b))));
  }
  
  ## Determine tick positions
  nSubFactor <- 2.44;
  
  atPretty1 <- pretty(xRange^2*sign(xRange),
                      u5.bias=u5.bias,
                      n=(pretty.n)^(1/nSubFactor));
  
  atPretty1old <- atPretty1;
  while (length(atPretty1) <= pretty.n) {
    atPretty1new <- subdivideSqrt(atPretty1,
                                  n=noiseFloor(minimum=2, (pretty.n)^(1/nSubFactor)));
    atPretty1 <- atPretty1new[atPretty1new <= max(abs(xRange^2))];
    atPretty1old <- atPretty1;
  }
  atPretty3 <- unique(sort(
    rep(atPretty1,
        each=length(unique(sign(xRange)))) * sign(xRange)));
  atPretty <- atPretty3[
    (atPretty3 >= head(xRange,1)^2*sign(head(xRange,1)) &
       atPretty3 <= tail(xRange, 1)^2*sign(tail(xRange, 1)))];
  
  xLabel <- sapply(atPretty, function(i){
    format(i,
           trim=TRUE,
           digits=2,
           big.mark=big.mark);
  });
  ## Transform to square root space
  atSqrt <- sqrt(abs(atPretty))*sign(atPretty);
  if (plot) {
    axis(side=side,
         at=atSqrt,
         labels=xLabel,
         las=las,
         cex.axis=cex.axis,
         ...);
  }
  
  invisible(nameVector(atSqrt, xLabel));
}

#' Apply noise floor and ceiling to numeric vector
#'
#' Apply noise floor and ceiling to numeric vector
#'
#' A noise floor is useful when detected numeric values are sometimes below
#' a clear noise threshold, and where some downstream ratio may be calculated
#' using these values. Applying a noise floor ensures the ratios and not
#' artificially higher, especially in cases where the values involved are
#' least reliable. This procedure is expected to produce more conservative
#' and appropriate ratios in that scenario.
#'
#' A ceiling is similar, values above the ceiling are set to the ceiling,
#' which is practical when values above a certain threshold are conceptually
#' similar to those at the threshold. One clear example is plotting
#' `-log10(Pvalue)` when the range of P-values might approach 1e-1000.
#' In this case, setting a ceiling of 50 conceptually equates P-values
#' below 1e-50, while also restricting the axis range of a plot.
#'
#' The ability to set values at the floor to a different value, using
#' `newValue` different from `minimum`, is intended to allow separation
#' of numeric values from the floor for illustrative purposes.
#'
#' @return
#' A numeric vector or matrix, matching the input type `x` where numeric
#' values are fixed to the `minimum` and `ceiling` values as defined
#' by `newValue` and `newCeiling`, respectively.
#'
#' @family jam numeric functions
#'
#' @param x numeric vector or matrix
#' @param minimum numeric floor value
#' @param newValue numeric, by default the same as the floor value. Sometimes
#'    it can be useful to define a different value, one example is to define
#'    values as NA, or another distinct number away from the floor.
#' @param adjustNA logical whether to change NA values to the newValue.
#' @param ceiling numeric value, optionally a ceiling. If defined, then values
#'    above the ceiling value are set to newCeiling.
#' @param newCeiling numeric value when ceiling is defined, values above the
#'    ceiling are set to this numeric value.
#' @param ... additional parameters are ignored.
#'
#' @examples
#' # start with some random data
#' n <- 2000;
#' x1 <- rnorm(n);
#' y1 <- rnorm(n);
#'
#' # apply noise floor and ceiling
#' x2 <- noiseFloor(x1, minimum=-2, ceiling=2);
#' y2 <- noiseFloor(y1, minimum=-2, ceiling=2);
#'
#' # apply noise floor and ceiling with custom replacement values
#' xm <- cbind(x=x1, y=y1);
#' xm3 <- noiseFloor(xm,
#'    minimum=-2, newValue=-3,
#'    ceiling=2, newCeiling=3);
#'
#' parMfrow <- par("mfrow");
#' par("mfrow"=c(2,2));
#' plotSmoothScatter(x1, y1);
#' plotSmoothScatter(x2, y2);
#' plotSmoothScatter(xm3);
#' par("mfrow"=parMfrow);
#'
#' @export
noiseFloor <- function
(x,
 minimum=0,
 newValue=minimum,
 adjustNA=FALSE,
 ceiling=NULL,
 newCeiling=ceiling,
 ...)
{
  ## Purpose is to apply a noise floor, that is, to set all values
  ## to be at least 'minimum' amount.
  ## This function performs no scaling or normalization.
  if (length(x) == 0) {
    return(x);
  }
  if (length(minimum) > 0) {
    if (adjustNA) {
      x[is.na(x) | (!is.na(x) & x < minimum)] <- newValue;
    } else {
      x[!is.na(x) & x < minimum] <- newValue;
    }
  }
  if (length(ceiling) > 0) {
    x[!is.na(x) & x > ceiling] <- newCeiling;
  }
  return(x);
}