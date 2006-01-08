require(maptools)
gSRDF <- readShapePoly(system.file("shapes/georgia.shp", package="spgwr")[1], proj4string=CRS("+proj=longlat +datum=NAD27"))
