# script for the sensitivity analysis of the application 
# with respect to the population density map


# load facebook population file ----
facePop = raster("data/popData/facebook_2020.tif")
save(facePop, file="savedOutput/global/facePop.RData")

# get frame/pop info ----
popMatSimple = makePopIntegrationTab(kmRes=5, pop=pop, domainPoly=kenyaPoly, 
                                     eastLim=eastLim, northLim=northLim, 
                                     mapProjection=projKenya, poppa=poppaKenya, 
                                     poppsub=poppsubKenya, stratifyByUrban=TRUE, 
                                     areaMapDat=adm1, subareaMapDat=adm2, 
                                     areaPolygonSubsetI=30)

popMatSimpleNeonatal = adjustPopMat(popMatSimple, poppaTarget=poppsubKenyaNeonatal, adjustBy="subarea")
easpaSimple = makeDefaultEASPA()
easpaSimple = easpaSimple[easpaSimple$area == "Nairobi",]
poppsubSimple = poppsubKenya
poppsubSimple = poppsubSimple[poppsubSimple$area == "Nairobi",]